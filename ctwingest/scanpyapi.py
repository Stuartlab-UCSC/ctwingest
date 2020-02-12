import scanpy as sc
import anndata
import pandas as pd
import numpy as np
from scipy.sparse.csr import csr_matrix
from geosketch import gs

def n_clusters(adata, cluster_solution_name):
    return len(adata.obs[cluster_solution_name].unique().dropna())


def get_expression(adata, use_raw=True, use_geosketch=False, geosketch_N=10000, geosketch_transform=None):
    """Grab expression and put into pandas dataframe."""
    if use_raw:
        ad = adata.raw
    else:
        ad = adata

    if isinstance(ad.X, csr_matrix):
        df = pd.DataFrame(ad.X.toarray(), index=ad.obs_names, columns=ad.var_names)
    else:
        df = pd.DataFrame(ad.X, index=ad.obs_names, columns=ad.var_names)

    if use_geosketch:
        if geosketch_transform==None:
            sc.tl.umap(adata, n_components=2)
            geosketch_transform = "umap"
        sketch_index = gs(adata.obsm[("X_" + geosketch_transform)], geosketch_N, replace=False)
        df = df.loc[sketch_index]
        
    return df.transpose()

from collections import Counter
from math import ceil


def multi_sketch(dimRed, fractions, clusters):
    """
    Do geometric sketches of the data given in dimRed, one per fraction
    """
    total_cells = dimRed.shape[0]
    percentages = ["pct" + str(int(i*100)) for i in fractions]
    all_counts = Counter(clusters)
    cluster_names = [str(i) for i in sorted([int(float(i)) for i in list(set(clusters))])]
    sketch_index_by_percentage = {}
    sketch_N = []

    sketch_df = pd.DataFrame(columns = percentages + ['full'], index = cluster_names)
    for key, value in all_counts.items():
        sketch_df.loc[key, 'full'] = value
    for i, fraction in enumerate(fractions):
        N = ceil(total_cells * fraction)
        print("total number of cells: ", total_cells, "; fraction: ", fraction, "; fraction # cells: ", N)
        sketch_N = sketch_N + [N]
        this_sketch_index = gs(dimRed, N, replace=False)
        sketch_index_by_percentage[percentages[i]] = this_sketch_index
        subset_counts = Counter(clusters.iloc[this_sketch_index])
        for key, value in subset_counts.items():
            sketch_df.loc[key, percentages[i]] = value
    return sketch_df, sketch_index_by_percentage, sketch_N


def std_gt_0_genes(centroids):
    """returns genes that have std != 0"""
    # remove clusters with std = 0 first
    centroids = centroids.loc[:, centroids.std(axis=0) != 0]
    # now remove genes with std = 0
    centroids = centroids.loc[centroids.std(axis=1) != 0, :]
    return centroids.index


def proportion_expressed_cluster(adata, clustering, use_raw=True):
    """
    outputs a dataframe [genes x cluster] that is the percentage that gene is expressed in a cluster
    :param ad: scanpy.Anndata
    :param cluster_solution_name: string, key accessor for ad.obs cluster_solution (
    :return: pandas.DataFrame
    """
    if use_raw:
        ad = adata.raw
    else:
        ad = adata

    pcent_df = pd.DataFrame(index=ad.var_names)
    for cluster_name in clustering.unique():
        cells_in_cluster = clustering.index[clustering == cluster_name]

        gt0 = (ad[cells_in_cluster].X > 0).sum(axis=0).transpose()
        pcent_df[cluster_name] = gt0 / float(len(cells_in_cluster))

    pcent_df.columns = pcent_df.columns.astype("str")

    return pcent_df


def all_genes_n(ad, use_raw):
    if use_raw:
        n = ad.raw.n_vars
    else:
        n = ad.n_vars

    return n


def get_obs(ad, obs_key):
    """Returns pandas series of observation annotations"""
    return ad.obs[obs_key]


def readh5ad(filename):
    return sc.read(filename)


def celltype_mapping(adata, cluster_name="louvain", mapping_name=None):
    """
    Assumes clusters have a many -> 1 mapping to cell types.
    :param adata:
    :param cluster_name:
    :param mapping_name: the key for the cell type assignment in ad.obs
    :return: pandas.Series
    """
    try:
        cluster_celltype = adata.obs[[cluster_name, mapping_name]]
        mapping = cluster_celltype.groupby(cluster_name).first()[mapping_name]
    except KeyError:
        mapping = None
        pass

    return mapping


def has_raw(ad):
    if ad.raw is None:
        return False
    elif isinstance(ad.raw, anndata.core.anndata.Raw):
        return True
    else:
        raise TypeError("anndata.Raw is of an unauthorized type: %s" % type(ad.raw))


def centroids(ad, cs_name="louvain", use_raw=True):
    cluster_solution = ad.obs[cs_name]
    # Calculate each centroid.

    if use_raw:
        genes = ad.raw.var_names
        adata = ad.raw
    else:
        genes = ad.var_names
        adata = ad

    centers = pd.DataFrame(index=genes)

    for cluster_name in cluster_solution.unique():
        cells_in_cluster = ad.obs.index[ad.obs[cs_name] == cluster_name]
        if isinstance(adata.X, np.ndarray):
            means = adata[cells_in_cluster].X.mean(axis=0).tolist()
        else: # might be sparse array?
            means = adata[cells_in_cluster].X.mean(axis=0).tolist()[0]

        centroid = pd.Series(
            means,
            index=genes
        )
        centers[cluster_name] = centroid

    centers.columns = centers.columns.astype("str")
    return centers


def get_xys(adata, key="X_umap"):
    return pd.DataFrame(adata.obsm[key], index=adata.obs_names)


def get_expression(adata, use_raw):
    if use_raw:
        ad = adata.raw
    else:
        ad = adata

    if isinstance(ad.X, csr_matrix):
        df = pd.DataFrame(ad.X.toarray(), index=ad.obs_names, columns=ad.var_names)
    else:
        df = pd.DataFrame(ad.X, index=ad.obs_names, columns=ad.var_names)

    return df.transpose()


def mito_genes(gene_symbols):
    """
    Filters a list of hugo gene symbols to mitochonrial genes.
    :param gene_symbols: a list of hugo gene symbols
    :return: a list of mitochondrial genes
    """
    mito_genes = prefixed(gene_symbols, "MT-")
    #mito_genes.extend(prefixed(gene_symbols, "MRPS"))
    #mito_genes.extend(prefixed(gene_symbols, "MRPL"))
    return mito_genes


def prefixed(strlist, prefix):
    """
    Filter a list to values starting with the prefix string
    :param strlist: a list of strings
    :param prefix: str
    :return: a subset of the original list to values only beginning with the prefix string
    """
    return [g for g in strlist if str(g).startswith(prefix)]
