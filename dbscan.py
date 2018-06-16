# -------------------------------------------------------------------------------
# Density-based spatial clustering of applications with noise (DBSCAN) applied
# on functional connectivity similarity matrix
# -------------------------------------------------------------------------------

from sklearn.cluster import DBSCAN
import numpy as np 
from matplotlib import pyplot as plt 
from nilearn.plotting import plot_matrix, plot_connectome
%matplotlib inline

def dbscan(similarity, static_fc, sub):
    dbscan = DBSCAN(metric = "euclidean").fit_predict(similarity[sub])
    plt.plot(np.asarray(dbscan) + 2)
    plt.title("State transitions", fontsize = 20)
    plot_matrix(static_fc[sub])
    plt.title("Static connectivity", fontsize = 20)  
    plot_matrix(similarity[sub])
    plt.title("Similarity matrix", fontsize = 20)
        
    n = 1
    
    for cluster in range(-1, dbscan.max()+1):
        label = dbscan == cluster
        percent = np.sum((label)/len(label) * 100)
        max_seq = maximum_sequence(label, cluster)
        mean_mat = np.mean(all_FC_sl[1][label], axis = 0)    #plot_matrix(mean_mat)
        mean_matrices.append(mean_mat)
        plot_matrix(mean_mat, auto_fit= True, vmin = static_fc[sub].min(), vmax = static_fc[sub].max()) #, axes = ax
        plt.title("State %d" %n, fontsize = 20)
        plt.suptitle("Percent of occurence: %d%%" %percent, backgroundcolor = "white")

        n += 1
