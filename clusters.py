# -*- coding: utf-8 -*-
"""
Created on Mon Jun  4 19:49:03 2018

@author: majd
"""

from scipy.spatial.distance import cdist, pdist
from sklearn.cluster import KMeans
from scipy .stats.mstats import zscore
from sklearn.metrics import silhouette_score

def elbow(data,K_range=range(1,20)):
    data=upper_tri(data)
    data=np.asarray(data)
    for k in K:
        kmeanModel = KMeans(n_clusters=k).fit(data)
        clusters=kmeanModel.cluster_centers_
        avgWithinSS.append(sum(np.min(cdist(data, clusters, 'euclidean'), axis=1)) / data.shape[0])
    fig = plt.figure(figsize=(15,15))
    ax = fig.add_subplot(111)
    ax.plot(K, avgWithinSS, 'b*-')
    plt.grid(True)
    plt.xlabel('Number of clusters')
    plt.ylabel('Average within-cluster sum of squares')
    plt.title('Elbow for KMeans clustering')
    plt.savefig('average com')
    plt.show()
    
def silhouette(data,K_range=range(1,20)):
    s = []
    data=upper_tri(data)
    data=np.asarray(data)
    for k in K_range:
        kmeans = KMeans(n_clusters=k)
        kmeans.fit(data)
        labels = kmeans.labels_
        centroids = kmeans.cluster_centers_
        s.append(silhouette_score(data, labels, metric='euclidean'))
    fig=plt.figure(figsize=(15,15))
    plt.plot(K_range,s)
    plt.ylabel("Silhouette")
    plt.xlabel("k")
    plt.title("Silhouette for K-means cell's behaviour")
    plt.savefig('Silhouette Kmeans')
    plt.show()