import os
import time
from functools import wraps

def timeit(func):
    @wraps(func)
    def timeit_wrapper(*args, **kwargs):
        import time
        start_time = time.perf_counter()
        result = func(*args, **kwargs)
        end_time = time.perf_counter()
        total_time = end_time - start_time
        print(f'Function {func.__name__} Took {total_time:.4f} seconds')
        return result
    return timeit_wrapper

class umap_analysis(object):
    """ A class used to compute to do umap cluster identification 
        of polymer conformations within an spatial configuration 
        of moieties. 
    """

    def __init__(self, all_pos):
        """Initiating the class umap_analysis.
         
        Parameters
        -----------

        all_pos : str
          Name of numpy array with the distances 
          to be analysed.
        
        """
        
        self.all_pos = all_pos

    def get_array(self):
        """Function to load the all pos distances.

        Parameters
        -----------

        None
        

        Returns
        --------

        all_pos: np.array
            np.array to load

        """
        import numpy as np

        all_pos=np.load(self.all_pos)

        return all_pos
        

    def umap_run(self, n_neigh, cwd, n_comp=2, min_d=0.0,
                 verb=True, rs=9):
        """Function to run umap on all_pos array.

        Parameters
        -----------

        n_neigh : class.int
            larger numbers forcus more on global properties 
            (and increase computational effort required)

        cwd : class.str
            path to save output file

        n_comp : class.int
            set the number of dimensions to reduce to

        min_d : class.int
            must be 0.0 if points are to be clustered 
            later with HDBSCAN.

        verb : class.boolean
            To show the progress bar

        rs : class.int
            random seed. Must be set to the same number 
            for reproducibility. 
            
        Returns
        --------

        X_embedded: np.array
            Clustered data

        """

        from sklearn.preprocessing import StandardScaler
        import numpy as np
        import matplotlib.pyplot as plt
        import umap

        all_pos=self.get_array()

        X_embedded = umap.UMAP(
        n_neighbors=n_neigh,  # larger numbers forcus more on global properties (and increase computational effort required)
        n_components=n_comp,  # set the number of dimensions to reduce to
        min_dist=min_d,    # must be 0.0 if points are to be clustered later with HDBSCAN
        verbose=verb,
        random_state=rs).fit(all_pos)


        plt.rcParams.update({'font.size': 18})
        fig=plt.figure(figsize=(7,6))
        ax = plt.subplot(111)

        plt.scatter(X_embedded.embedding_[:,0], X_embedded.embedding_[:,1], s=5)
        plt.xlabel("$ UMAP_{1} $")
        plt.ylabel("$UMAP_{2} $")

        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)

        plt.savefig('{cwd}/umpap_output.png'.format(cwd=cwd), bbox_inches='tight', dpi=200)

        return X_embedded


    def hdbscan_cluster(self, X_embedded, min_cs, epsilon, cwd, method="leaf"):
        """Function to hdbscan on umap clustered data.

        Parameters
        -----------

        X_embedded : np.array
            umap clustered data

        min_cs : class.int
            minimum number of points for a group to be considered 
            a cluster.

        epsilon : class.int
            clusters closer than this distance apart will be merged.

        cwd : class.str
            path to save output file.
            
        method : class.str
            cluster_method_selection from hdbscan to determine how it 
            selects flat clusters from the cluster tree hierarchy.

        Returns
        --------

        y_pred : np.array
            cluster label assignment for each data point in the umap 
            cluster data.

        """

        from sklearn.preprocessing import StandardScaler
        import numpy as np
        import matplotlib.pyplot as plt
        import hdbscan

        hdb = hdbscan.HDBSCAN(
        min_cluster_size= min_cs,              # minimum number of points for a group to be considered a cluster
        cluster_selection_epsilon=epsilon,      # clusters closer than this distance apart will be merged
        cluster_selection_method=method)   # changes the way clusters are initialised. "eom" tends to lead to fewer clusters than "leaf")

        hdb.fit(X_embedded.embedding_)


        y_pred = hdb.fit_predict(X_embedded.embedding_)
        print(f"Not in a cluster: {sum(y_pred==-1)} ({100 * sum(y_pred==-1) / y_pred.size:.2f}%)")
  
        plt.figure(figsize=(6,4))
        fig, ax = plt.subplots()

        ax.scatter(X_embedded.embedding_[:,0], X_embedded.embedding_[:,1], c=y_pred, cmap='Paired', s=5)

        for cluster in range(max(y_pred)+1):     
            mean = np.mean(X_embedded.embedding_[y_pred==cluster], axis=0)
            ax.text(
                mean[0], mean[1],
                cluster,
                fontsize=20, 
                fontweight = 'bold'
                )
    
        # Definitions of the plot
        plt.rcParams['font.family'] = 'sans-serif'
        plt.rcParams['axes.spines.top'] = False
        plt.rcParams['axes.spines.left'] = False
        plt.rcParams['axes.spines.right'] = False
        plt.rcParams['axes.spines.bottom'] = False
        plt.rcParams['font.size'] = 10

        # Options for the plot
        ax.get_yaxis().set_visible(False)
        ax.get_xaxis().set_visible(False)
        ax.xaxis.set_visible(False)
        ax.yaxis.set_visible(False)

        plt.savefig('{cwd}/hdbscan_output.png'.format(cwd=cwd), bbox_inches='tight', dpi=200)

        return y_pred


    def cluster_order_appearance(self, X_embedded, y_pred, cwd):
        """Function to ouput hdbscan with the groups labeled ordered 
           with their order of appearance in the trajectory.

        Parameters
        -----------

        X_embedded : np.array
            umap clustered data

        y_pred : np.array
            cluster label assignment for each data point in the 
            umap cluster data.

        cwd : class.str
            path to save output file.
            
        Returns
        --------

        None:
            A plot with the results of the calculation.
            
        """

        import numpy as np
        import matplotlib.pyplot as plt

        #ignore noise
        labels = np.unique(y_pred) 

        unique_clusters, frame_first_seen = np.unique(y_pred, return_index=True)
        sorter = np.argsort(frame_first_seen)
        cluster_order = unique_clusters[sorter]

        cluster_labels = np.full_like(y_pred, fill_value=-1)
        for index, cluster in enumerate(cluster_order, start=1):   
            if cluster == -1:
                cluster_labels[y_pred==cluster] = -1
                index-= 1
                continue
    
        cluster_labels[y_pred==cluster] = index
        mean_pdb = []
        plt.scatter(X_embedded.embedding_[:,0], X_embedded.embedding_[:,1],c=y_pred, cmap='Paired')

        for cluster in np.unique(cluster_labels):
            if cluster == -1:
                continue
    
            mean = np.mean(X_embedded.embedding_[cluster_labels==cluster], axis=0)
            mean_pdb.append(mean)
            plt.text(
            mean[0], mean[1],
            cluster,
            fontsize=13,
            weight="bold"
            )

        plt.savefig('{cwd}/hdbscan_output_ordered.png'.format(cwd=cwd), bbox_inches='tight', dpi=200)


    def find_closest_point(self, array, target_point):
        """Function to find the closest point from an array to a given target point

        Parameters
        -----------

        array : np.array
            array of numbers from which the value will be found

        target_point : class.float
            target number

        
        Returns
        --------

        closest_point: np.array
            numpy array with the closest number from array to the target point.
            
        """

        import numpy as np
        from scipy.spatial.distance import cdist

        # Calculate pairwise distances
        distances = cdist(array, np.array([target_point]))
        # Find the index of the closest row
        closest_index = np.argmin(distances)
        # Get the closest point
        closest_point = array[closest_index]  

        return closest_point


    def get_average_rep_cluster(self, X_embedded, y_pred):
        """Function to get the X_embedded point that represents the 
           mean of each umap cluster.

        Parameters
        -----------

        X_embedded : np.array
            umap clustered data

        y_pred : np.array
            cluster label assignment for each data point in 
            the umap cluster data.

        cwd : class.str
            path to save output file.
            
        Returns
        --------

        cluster_points : np.array
            np.array with the X_embedded values that represent 
            the mean of each cluster.

        """

        import numpy as np

        unique_clusters = np.unique(y_pred)
        cluster_points=[]

        for cluster in unique_clusters:
            if cluster == -1:
                continue

            center = np.mean(X_embedded.embedding_[y_pred == cluster], axis=0)
            point = self.find_closest_point(X_embedded.embedding_, center)
            cluster_points.append(point)

        return cluster_points
 
