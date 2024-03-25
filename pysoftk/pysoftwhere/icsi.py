import MDAnalysis as mda
import networkx as nx
import MDAnalysis.analysis.distances
import scipy.stats as stats
import numpy as np
from itertools import product as product


class icsi:
    """
    Calculates the intrinsic surface density (ISD) of a cluster in a molecular
    simulation.

    This class implements the Intrinsic Surface Density (ISD) calculation as
    described in [Reference]. It takes a universe (`u`), cluster definition
    information, and options for frame and binning as input and returns the
    intrinsic radius, spherical radius, and interface values.

    Args:
        u: An MDAnalysis universe object.
        cluster_resids: A list of residue indices that define the cluster.
        cluster_atoms_positions: A 3D NumPy array of atomic positions in the
            cluster across all frames.
        core_sel_atoms_positions: A 3D NumPy array of atomic positions for the
            core region of the cluster.
        shell_sel_atoms_positions: A 3D NumPy array of atomic positions for the
            shell region of the cluster.
        frame: The frame index to use for the calculation (default: -1 for
            the last frame).
        no_bins: The number of bins for the 2D histogram (default: 31).
        no_random_points: The number of random points to use for
            normalization (default: 10000).
        normalisation_run: Whether to use random points for normalization
            (default: False).

    Returns:
        intrinsic_r: A list of intrinsic radii.
        spherical_r: A list of spherical radii.
        interface_vals: A 2D NumPy array of interface values.

    Raises:
        ValueError: If an invalid frame index is provided.


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
        self.cluster_atoms_positions=cluster_atoms_positions
        self.core_sel_atoms_positions=core_sel_atoms_positions
        self.shell_sel_atoms_positions=shell_sel_atoms_positions
        self.no_random_points=no_random_points
        self.normalisation_run=normalisation_run
        
        self._spherical_intrinsic_density()
        
    def _spherical_intrinsic_density(self):
        """
        Calculates the spherical intrinsic density of the cluster.

        This method performs the core calculations for the ISD, including:

        - Converting Cartesian coordinates to spherical coordinates.
        - Selecting core and shell atoms.
        - Calculating binned statistics for shell and interface positions.
        - Calculating intrinsic radii.

        """
        
        def cartesian_to_spherical(xyz):
            """
            Converts Cartesian coordinates to spherical coordinates.

            Args:
                xyz: A 3D NumPy array of Cartesian coordinates.

            Returns:
                A 3D NumPy array of spherical coordinates.
            """
            ptsnew = np.hstack((xyz, np.zeros(xyz.shape)))
            xy = xyz[:,0]**2 + xyz[:,1]**2
            ptsnew[:,3] = np.sqrt(xy + xyz[:,2]**2)
            ptsnew[:,4] = np.arctan2(np.sqrt(xy), xyz[:,2]) # for elevation angle defined from Z-axis down
            #ptsnew[:,4] = np.arctan2(xyz[:,2], np.sqrt(xy)) # for elevation angle defined from XY-plane up
            ptsnew[:,5] = np.arctan2(xyz[:,1], xyz[:,0])
            return ptsnew
        
        self.u.trajectory[self.frame]
        box_=self.u.dimensions
        
        cluster_sel = self.u.select_atoms('resid '+" ".join([str(i) for i in self.cluster_resids]))
        com = sum([(1/(sum(cluster_sel.masses)))*cluster_sel.masses[i]*self.cluster_atoms_positions[i] for i in range(len(cluster_sel))])
        ####calc spherical polars of core
        core_array = cartesian_to_spherical(self.core_sel_atoms_positions-com)
        r_core_tmp=core_array[:,3]
        theta_core_tmp=core_array[:,4]
        phi_core_tmp=core_array[:,5]
        ####calc spherical polars of shell
        
        if self.normalisation_run == False:
            shell_array = cartesian_to_spherical(self.shell_sel_atoms_positions-com)
        elif self.normalisation_run == True:
            shell_array = cartesian_to_spherical(np.random.uniform(-self.u.dimensions[0]/2,self.u.dimensions[0]/2,size=(self.no_random_points,3)))
        
        r_shell_sphere_tmp=shell_array[:,3]
        theta_shell_sphere_tmp=shell_array[:,4]
        phi_shell_sphere_tmp=shell_array[:,5]

        ###calc shell bins
        
        shell_pos=stats.binned_statistic_2d([np.cos(i) for i in theta_shell_sphere_tmp],phi_shell_sphere_tmp,r_shell_sphere_tmp,
                                            bins=[np.linspace(-1,1,self.no_bins),np.linspace(-np.pi,np.pi,self.no_bins)],expand_binnumbers=True,statistic='count')
    
    
        ####calc interface positions
        ###np.max finds edge of core!
        interface=stats.binned_statistic_2d([np.cos(i) for i in theta_core_tmp], phi_core_tmp,r_core_tmp,
                                            bins=[np.linspace(-1,1,self.no_bins),np.linspace(-np.pi,np.pi,self.no_bins)],
                                            expand_binnumbers=True,statistic=np.max)

        interface_vals=interface.statistic.copy()    
    
        while np.count_nonzero(~np.isnan(interface_vals))!=(len(np.linspace(0,box_[0],num=self.no_bins))-1)*(len(np.linspace(0,box_[0],num=self.no_bins))-1):
            for i in range(len(interface_vals)):
                for j in range(len(interface_vals)):
    
                    if np.isnan(interface_vals[i][j]):
                        n_i=[i-1,i,i+1]
                        n_j=[j-1,j,j+1]
                        interface_vals[i][j]=np.nanmean([interface_vals[divmod(ip,len(interface_vals))[1]][divmod(jp,len(interface_vals))[1]] for ip,jp in product(n_i,n_j)])
    
    
        intrinsic_r,spherical_r = [],[]
    
        for i in range(len(r_shell_sphere_tmp)):
            bin_X=shell_pos.binnumber[0][i]-1
            bin_Y=shell_pos.binnumber[1][i]-1
    
    
            interface_pos=interface_vals[bin_X][bin_Y]
            intrinsic_r.append(-interface_pos+r_shell_sphere_tmp[i])
            spherical_r.append(r_shell_sphere_tmp[i])
            
        
        self.intrinsic_r,self.spherical_r, self.interface_vals = intrinsic_r, spherical_r, interface_vals
       
    def __call__(self):
        """
        .
        """     

        return self.intrinsic_r,self.spherical_r, self.interface_vals


