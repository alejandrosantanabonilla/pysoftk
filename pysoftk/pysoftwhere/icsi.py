import MDAnalysis as mda
import networkx as nx
import MDAnalysis.analysis.distances
import scipy.stats as stats
import numpy as np
from itertools import product as product


class icsi:
    """
    Calculates the **Intrinsic Surface Density (ISD)** of a defined cluster
    within a molecular dynamics simulation.

    The Intrinsic Surface Density (ISD) provides a quantitative measure of
    the local density fluctuations at the interface of a molecular cluster.
    This class implements the ISD calculation, as described in scientific
    literature, to characterize the roughness or compactness of a cluster's
    surface. It processes molecular positions, converts them to spherical
    coordinates, and performs binning to determine the intrinsic radius,
    spherical radius, and interface values.

    Parameters
    ----------
    u : MDAnalysis.Universe
        An **MDAnalysis Universe object** containing the trajectory and
        topology information of the molecular system. This universe provides
        the context for atom selections and positions.
    cluster_resids : list of int
        A **list of residue IDs** (integers) that collectively define the
        cluster of interest. All atoms belonging to these residues will be
        considered part of the cluster.
    cluster_atoms_positions : numpy.ndarray
        A **3D NumPy array** of shape `(n_atoms, 3)` representing the Cartesian
        coordinates (x, y, z) of all atoms within the `cluster_resids`
        for a specific frame.
    core_sel_atoms_positions : numpy.ndarray
        A **3D NumPy array** of shape `(n_core_atoms, 3)` representing the
        Cartesian coordinates (x, y, z) of atoms defined as the **core region**
        of the cluster. These atoms are used to establish the "inner" boundary
        for the intrinsic radius calculation.
    shell_sel_atoms_positions : numpy.ndarray
        A **3D NumPy array** of shape `(n_shell_atoms, 3)` representing the
        Cartesian coordinates (x, y, z) of atoms defined as the **shell region**
        of the cluster. These atoms represent the "outer" boundary and are
        used to probe the surface.
    frame : int, optional
        The **index of the trajectory frame** to use for the ISD calculation.
        A value of -1 (default) indicates the last frame of the trajectory.
    no_bins : int, optional
        The **number of bins** to use for constructing the 2D histogram in
        spherical coordinate space ($\cos(\theta)$ and $\phi$). A higher
        number of bins provides finer resolution but may lead to sparser data.
        Defaults to 31.
    no_random_points : int, optional
        The **number of random points** to generate for normalization purposes
        when `normalisation_run` is set to `True`. These points are uniformly
        distributed within the simulation box. Defaults to 10000.
    normalisation_run : bool, optional
        A flag indicating whether the calculation should use **randomly
        generated points for normalization** instead of the actual `shell_sel_atoms_positions`.
        Set to `True` for normalization runs, `False` (default) for actual
        ISD calculation.

    Attributes
    ----------
    intrinsic_r : list of float
        A list containing the **calculated intrinsic radii** for each shell atom.
        The intrinsic radius represents the distance from the cluster's effective
        surface to the shell atom.
    spherical_r : list of float
        A list containing the **calculated spherical radii** for each shell atom.
        The spherical radius is the radial distance of the shell atom from the
        center of mass of the cluster.
    interface_vals : numpy.ndarray
        A **2D NumPy array** representing the calculated interface values
        across the binned spherical coordinate space. This array effectively
        maps the "edge" of the core region in different angular directions.

    Raises
    ------
    ValueError
        If an invalid `frame` index is provided (e.g., out of trajectory bounds).

    See Also
    --------
    MDAnalysis.Universe : For details on the MDAnalysis Universe object.
    scipy.stats.binned_statistic_2d : For details on the 2D histogramming function.

    Notes
    -----
    The calculation of the intrinsic surface density involves:
    1.  Defining a cluster, its core, and its shell regions.
    2.  Translating all atomic coordinates to the center of mass of the cluster.
    3.  Converting Cartesian coordinates of core and shell atoms to spherical
        coordinates (r, $\theta$, $\phi$).
    4.  Binning the spherical coordinates to create a 2D histogram.
    5.  Determining the "interface" position, typically by finding the maximum
        radial distance of core atoms within each angular bin.
    6.  Calculating the intrinsic radius for each shell atom as the difference
        between its spherical radius and the interface position in its
        corresponding angular bin.

    The method handles missing data (NaNs) in the interface values by
    averaging from neighboring bins, ensuring a continuous interface definition.
    """
    def __init__(self,
                 u,
                 cluster_resids,
                 cluster_atoms_positions,
                 core_sel_atoms_positions,
                 shell_sel_atoms_positions,
                 frame=-1,
                 no_bins=31,
                 no_random_points=10000,
                 normalisation_run=False):

        self.u = u
        self.frame = frame
        self.no_bins = no_bins
        self.cluster_resids = cluster_resids
        self.cluster_atoms_positions = cluster_atoms_positions
        self.core_sel_atoms_positions = core_sel_atoms_positions
        self.shell_sel_atoms_positions = shell_sel_atoms_positions
        self.no_random_points = no_random_points
        self.normalisation_run = normalisation_run

        self._spherical_intrinsic_density()

    def _spherical_intrinsic_density(self):
        """
        Performs the core calculations for the **spherical intrinsic density (ISD)**.

        This method encapsulates the multi-step process of computing the ISD,
        including coordinate transformations, atom selections, 2D histogramming,
        and the calculation of intrinsic and spherical radii. It also
        handles the imputation of missing interface values and, optionally,
        uses random points for normalization runs.

        The steps involved are:
        1.  Loading the specified trajectory frame from the MDAnalysis Universe.
        2.  Calculating the **center of mass (COM)** for the entire cluster.
        3.  Translating all atomic positions to the COM.
        4.  Converting Cartesian coordinates of both core and shell atoms
            (or random points during a normalization run) to spherical coordinates (r, $\theta$, $\phi$).
        5.  Performing 2D binning of shell atom positions (or random points)
            and core atom positions in ($\cos(\theta)$, $\phi$) space to determine
            `shell_pos` counts and `interface` values (max radial distance of core atoms).
        6.  Iteratively filling any `NaN` (Not-a-Number) values in the `interface_vals`
            array by averaging from valid neighboring bins, ensuring a smooth
            interface definition.
        7.  Calculating the `intrinsic_r` and `spherical_r` for each shell atom
            based on its position relative to the calculated interface.

        This method populates the `self.intrinsic_r`, `self.spherical_r`,
        and `self.interface_vals` attributes of the class.

        Returns
        -------
        None
            This method does not return any value directly. It updates the
            instance attributes: `self.intrinsic_r`, `self.spherical_r`, and
            `self.interface_vals`.
        """

        def cartesian_to_spherical(xyz):
            """
            Converts a set of 3D Cartesian coordinates to spherical coordinates.

            The spherical coordinates are defined as:
            - **r (radius)**: The Euclidean distance from the origin.
            - **theta ($\theta$)**: The polar angle (inclination) from the positive Z-axis
              (range: 0 to $\pi$ radians).
            - **phi ($\phi$)**: The azimuthal angle in the XY-plane from the positive X-axis
              (range: $-\pi$ to $\pi$ radians).

            Parameters
            ----------
            xyz : numpy.ndarray
                A **3D NumPy array** of shape `(N, 3)` where `N` is the number
                of points, and each row contains the `[x, y, z]` Cartesian
                coordinates.

            Returns
            -------
            numpy.ndarray
                A **3D NumPy array** of shape `(N, 6)`. Each row contains
                `[x, y, z, r, theta, phi]` where `r`, `theta`, and `phi` are
                the calculated spherical coordinates.
            """
            ptsnew = np.hstack((xyz, np.zeros(xyz.shape)))
            xy = xyz[:, 0]**2 + xyz[:, 1]**2
            ptsnew[:, 3] = np.sqrt(xy + xyz[:, 2]**2)
            ptsnew[:, 4] = np.arctan2(np.sqrt(xy), xyz[:, 2])  # for elevation angle defined from Z-axis down
            #ptsnew[:,4] = np.arctan2(xyz[:,2], np.sqrt(xy)) # for elevation angle defined from XY-plane up
            ptsnew[:, 5] = np.arctan2(xyz[:, 1], xyz[:, 0])
            return ptsnew

        self.u.trajectory[self.frame]
        box_ = self.u.dimensions

        cluster_sel = self.u.select_atoms('resid ' + " ".join([str(i) for i in self.cluster_resids]))
        com = sum([(1 / (sum(cluster_sel.masses))) * cluster_sel.masses[i] * self.cluster_atoms_positions[i] for i in range(len(cluster_sel))])

        ####calc spherical polars of core
        core_array = cartesian_to_spherical(self.core_sel_atoms_positions - com)
        r_core_tmp = core_array[:, 3]
        theta_core_tmp = core_array[:, 4]
        phi_core_tmp = core_array[:, 5]

        ####calc spherical polars of shell
        if self.normalisation_run == False:
            shell_array = cartesian_to_spherical(self.shell_sel_atoms_positions - com)
        elif self.normalisation_run == True:
            shell_array = cartesian_to_spherical(np.random.uniform(-self.u.dimensions[0] / 2, self.u.dimensions[0] / 2, size=(self.no_random_points, 3)))

        r_shell_sphere_tmp = shell_array[:, 3]
        theta_shell_sphere_tmp = shell_array[:, 4]
        phi_shell_sphere_tmp = shell_array[:, 5]

        ###calc shell bins

        shell_pos = stats.binned_statistic_2d([np.cos(i) for i in theta_shell_sphere_tmp], phi_shell_sphere_tmp, r_shell_sphere_tmp,
                                              bins=[np.linspace(-1, 1, self.no_bins), np.linspace(-np.pi, np.pi, self.no_bins)], expand_binnumbers=True, statistic='count')

        ####calc interface positions
        ###np.max finds edge of core!
        interface = stats.binned_statistic_2d([np.cos(i) for i in theta_core_tmp], phi_core_tmp, r_core_tmp,
                                              bins=[np.linspace(-1, 1, self.no_bins), np.linspace(-np.pi, np.pi, self.no_bins)],
                                              expand_binnumbers=True, statistic=np.max)

        interface_vals = interface.statistic.copy()

        while np.count_nonzero(~np.isnan(interface_vals)) != (len(np.linspace(0, box_[0], num=self.no_bins)) - 1) * (len(np.linspace(0, box_[0], num=self.no_bins)) - 1):
            for i in range(len(interface_vals)):
                for j in range(len(interface_vals)):

                    if np.isnan(interface_vals[i][j]):
                        n_i = [i - 1, i, i + 1]
                        n_j = [j - 1, j, j + 1]
                        interface_vals[i][j] = np.nanmean([interface_vals[divmod(ip, len(interface_vals))[1]][divmod(jp, len(interface_vals))[1]] for ip, jp in product(n_i, n_j)])

        intrinsic_r, spherical_r = [], []

        for i in range(len(r_shell_sphere_tmp)):
            bin_X = shell_pos.binnumber[0][i] - 1
            bin_Y = shell_pos.binnumber[1][i] - 1

            interface_pos = interface_vals[bin_X][bin_Y]
            intrinsic_r.append(-interface_pos + r_shell_sphere_tmp[i])
            spherical_r.append(r_shell_sphere_tmp[i])

        self.intrinsic_r, self.spherical_r, self.interface_vals = intrinsic_r, spherical_r, interface_vals

    def __call__(self):
        """
        Enables the `icsi` object to be called like a function, returning
        the computed intrinsic radius, spherical radius, and interface values.

        This magic method provides a convenient way to access the results
        of the ISD calculation after the `icsi` object has been initialized
        and its internal `_spherical_intrinsic_density` method has executed.

        Returns
        -------
        tuple
            A tuple containing three elements:
            - **intrinsic_r** (list of float): The calculated intrinsic radii.
            - **spherical_r** (list of float): The calculated spherical radii.
            - **interface_vals** (numpy.ndarray): The 2D NumPy array of interface values.
        """

        return self.intrinsic_r, self.spherical_r, self.interface_vals
