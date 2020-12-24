from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import pdist

import pandas as pd
import numpy as np
import seaborn as sns


def plot_heatmap(
    dataframe,
    cmap = sns.cubehelix_palette(
        16,
        start = 2.1,
        rot = -0.2,
        gamma = 0.4,
        hue = 1,
        light = 0.8,
        dark = 0.2,
        as_cmap = True
    )
):
    row_linkage = linkage(
        pdist(dataframe.to_numpy()), method = "ward", optimal_ordering = True
    )

    row_cluster = fcluster(
        row_linkage, 8, criterion = 'maxclust'
    )

    fun_palette = dict(
        zip(
            np.unique(row_cluster), sns.color_palette("Paired", len(np.unique(row_cluster)))
        )
    )

    row_color = pd.DataFrame(
        pd.Series(
            row_cluster,
            index = dataframe.index,
            name = 'clusters').map(fun_palette)
    )

    col_linkage = linkage(
        pdist(dataframe.T.to_numpy()), method = "weighted", optimal_ordering = True
    )

    fig_heatmap = sns.clustermap(
        dataframe, row_linkage = row_linkage,
        col_linkage = col_linkage, cmap = cmap,
        cbar_kws = {'orientation': 'horizontal'},
        cbar_pos=(0.05, 1.01, 0.4, 0.02), method = "weighted",
        row_colors = row_color, dendrogram_ratio=0.05,
        colors_ratio=0.02, yticklabels = True
    )

    return fig_heatmap
