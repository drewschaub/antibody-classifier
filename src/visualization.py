import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import scipy.cluster.hierarchy as sch

def plot_rmsd_heatmap(rmsd_matrix, structure_names, title='RMSD Heatmap'):
    fig, ax = plt.subplots(figsize=(10, 8))
    sns.heatmap(rmsd_matrix, 
                xticklabels=structure_names,
                yticklabels=structure_names, 
                cmap='Blues_r',
                square=True,
                annot=True,
                fmt='.2f')
    plt.title(title)
    plt.tight_layout()
    plt.show()

def plot_dendrogram(rmsd_matrix, labels=None, method='average'):
    # Linkage
    linked = sch.linkage(rmsd_matrix, method=method)
    # Dendrogram
    plt.figure()
    sch.dendrogram(linked, labels=labels, orientation='right')
    plt.show()
