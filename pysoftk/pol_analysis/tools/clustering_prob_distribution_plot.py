import MDAnalysis as mda 
import MDAnalysis.analysis.distances
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats
import MDAnalysis.analysis.rms
import networkx as nx
import MDAnalysis.analysis.rdf as rdf
import matplotlib.patches as mpatches
import concurrent.futures
import os
import time
from functools import wraps
from tqdm.auto import tqdm
import itertools
import pandas as pd


def total_number_clusters(dataframe, start, color ='#fde725'):

    '''The dataframe input should be the output of the spatial_clustering'''

    '''It is recommended to do this calculation from the frame when the micelle is equilibrated'''

    df = dataframe.iloc[int(start):]


    df_exploded = df.explode('micelle_size')

    

    counts = df_exploded['micelle_size'].value_counts().reset_index()
    counts.columns = ['micelle_size', 'counts']
    counts_sorted = counts.sort_values(by=['micelle_size'], ascending=True)

    max_cluster_size = counts['micelle_size'].max()
    mean_micelle_size = df_exploded['micelle_size'].mean() 

    print(df_exploded['micelle_size'].iloc[:10])
    print(df_exploded['micelle_size'].iloc[:10].mean())

    total_clusters = counts['counts'].sum()
    bins= list(range(max_cluster_size  +1))
    prob_size = np.zeros(max_cluster_size+1)
 
    for i in range(len(prob_size)):
        for k in range(len(counts_sorted['micelle_size'])):
            if counts_sorted['micelle_size'][k] == i+1:
                prob_size[i+1] = (counts_sorted['counts'][k]/total_clusters)*counts_sorted['micelle_size'][k]
    

    prob_size_cluster_size = prob_size/mean_micelle_size
    plt.rcParams.update({'font.size': 18})
    plt.rcParams['axes.labelsize'] = 18


    plt.rcParams['axes.labelsize'] = 18
    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['axes.spines.top'] = False
    plt.rcParams['axes.spines.right'] = False
    plt.rcParams['xtick.major.pad']='5'
    plt.rcParams['ytick.major.pad']='12'
    plt.rcParams['font.size'] = 18
    plt.rcParams['axes.titlesize'] = 15
    plt.rcParams['axes.labelsize'] = 15
    plt.rcParams['axes.linewidth'] = 2 
    fig, ax = plt.subplots()
    ax.bar(bins, prob_size_cluster_size, edgecolor=color, color=color)
    ax.xaxis.set_tick_params(width=2)
    ax.yaxis.set_tick_params(width=2)

    ax.set_ylabel('Probability', fontsize = 20, labelpad = 10)


    plt.xlabel('$N_{\mathrm{agg}}$', fontsize = 20, labelpad = 10)
    ax.set_xlim(0.05, max_cluster_size +1)
    ax.set_yticks(np.arange(0, 0.5, 0.1))
    ax.set_ylim(0,0.4)
    ax.set_xticks(np.arange(2, (max_cluster_size +1), 2))

    ax.tick_params(axis='x', labelsize=18)

    fig.tight_layout()
    


