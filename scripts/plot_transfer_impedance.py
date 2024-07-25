import sys, math
import h5py
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

font = {'family' : 'sans-serif',
        'sans-serif': 'Arial',
        'style': 'normal',
        'weight': 'normal',
        'size'   : 24 }
matplotlib.rc('font', **font)

def plot_lower_tri_heatmap(df, output="transfer_impedance_matrix.png"):
    mask = np.zeros_like(df, dtype=bool)
    mask[np.triu_indices_from(mask)] = True

    # Want diagonal elements as well
    mask[np.diag_indices_from(mask)] = False

    # Set up the matplotlib figure
    f, ax = plt.subplots(figsize=(15, 9))

    # Generate a custom diverging colormap
    cmap = sns.diverging_palette(220, 10, as_cmap=True)

    # Draw the heatmap with the mask and correct aspect ratio
    sns_plot = sns.heatmap(df, mask=mask, cmap=cmap, center=0,
                           square=True, linewidths=.5, cbar_kws={"shrink": .5})
    # save to file
    fig = sns_plot.get_figure()
    fig.savefig(output)

    

fname = sys.argv[1]

f = h5py.File(fname,"r")
z_array = f["Populations"]["GC"]["Cell Clamp Results"]["transfer impedance"]["Attribute Value"][:]
f.close()

N = int(math.sqrt(len(z_array)))
z_matrix = z_array.reshape((N, N))

plot_lower_tri_heatmap(z_matrix)

plt.show()
    
