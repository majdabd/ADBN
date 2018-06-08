# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
import matplotlib.pyplot as plt 
import numpy as np
def plot_dfc(states,dwell,freq,mod='correlation',labels=None,title='DFC'):
    fig=plt.figure(figsize=(20,20))
    for i,s in enumerate(states[mod]):
        np.fill_diagonal(s,0)
        ax=fig.add_subplot(3,4,i+1) 
        if dwell and freq is not None:
            ax.set_title("state %d \n Mean Dwell Time: " %(i+1) + 
                            "%.f TR \n Mean frequency: %d %% "
                            %(dwell[mod][i],freq[mod][i]),fontsize=12)
        im=plt.imshow(s, interpolation='nearest', cmap=plt.cm.RdBu_r,vmin=-abs(s.max()),vmax=abs(s.max()))
        plt.colorbar(im, fraction=0.07, pad=0.04)
        if labels is not None:
            x_ticks = plt.xticks(range(len(labels)),labels, rotation=90)
            y_ticks = plt.yticks(range(len(labels)),labels)
            plt.rc('xtick', labelsize=10) 
            plt.rc('ytick', labelsize=10) 
        fig.suptitle(title,fontsize=30)
        fig.subplots_adjust(hspace=1,wspace=0.4)
        plt.savefig(title,bbox_inches='tight') 
    plt.show()