from pathlib import Path
import scanpy as sc
import conga
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from conga import util
import scipy
import os
import sys
import anndata
from conga.preprocess import retrieve_tcrs_from_adata
from conga.preprocess import read_adata
from conga.preprocess import normalize_and_log_the_raw_matrix
from conga.tcrdist.tcr_distances import TcrDistCalculator

def run_CoNGA(features_file, tcr_datafile, gex_datafile, organism, outfile_prefix,
         gex_datatype, clones_file, outfile_prefix_for_qc_plots, working_directory, print_versions = True):

    if print_versions:
        print('scanpy version: ' + sc.__version__)
        print('scipy version: ' + scipy.__version__)
        print('pandas version: ' + pd.__version__)
        print('anndata version: ' + anndata.__version__)

    #set working directory (inherited from the R session)
    os.chdir(working_directory)
    tcrdist_calculator = TcrDistCalculator(organism)
    #make directories
    os.makedirs(os.path.dirname(os.path.join(os.getcwd(),outfile_prefix)), exist_ok=True)
    #initializing clonotype file and kPCA
    conga.tcrdist.make_10x_clones_file.make_10x_clones_file( tcr_datafile, organism, clones_file)
    conga.preprocess.make_tcrdist_kernel_pcs_file_from_clones_file( clones_file, organism )

    # TODO: this is for debugging:
    adataTest = sc.read_10x_h5(gex_datafile, gex_only=True )
    print(adataTest)
    print('Index length: ' + str(len(adataTest.obs.index)))
    print('Index unique length: ' + str(len(adataTest.obs.index.unique())))

    adataTest.var_names_make_unique()
    print(adataTest)
    print('Index length after make_unique: ' + str(len(adataTest.obs.index)))
    print('Index unique length after make_unique: ' + str(len(adataTest.obs.index.unique())))

    x = np.repeat([True, False], 1556/2)
    print('After subset')
    adata2 = adataTest[x,:].copy()
    print(adata2.obs.index.is_unique)
    print(adata2)
    print('Index length after make_unique: ' + str(len(adata2.obs.index)))
    print('Index unique length after make_unique: ' + str(len(adata2.obs.index.unique())))

    adata = conga.preprocess.read_dataset(gex_datafile, gex_datatype, clones_file )
    genes_df = pd.read_csv(features_file, header=None)
    os.makedirs(os.path.dirname(os.path.join(os.getcwd(),outfile_prefix_for_qc_plots)), exist_ok=True)

    # store the organism info in adata
    adata.uns['organism'] = organism
    adata.uns['force_variable_genes'] = list(genes_df[0].values)
    #scanpy gene expression preprocessing wrapper
    adata = conga.preprocess.filter_and_scale( 
        adata, 
        min_genes_per_cell=200,
        max_genes_per_cell=3500,
        max_percent_mito=0.1,
        outfile_prefix_for_qc_plots = outfile_prefix_for_qc_plots
    )
    
    # print out number of cells in adata
    print('Total cells in adata object: ' + str(adata.shape[0]))
    
    #create a second adata object that stores one representative clone per clonotype (specifically for gene expression analysis)
    adata2 = conga.preprocess.reduce_to_single_cell_per_clone(adata)
    adata2 = conga.preprocess.cluster_and_tsne_and_umap( adata2 )

    os.makedirs(os.path.dirname(os.path.join(os.getcwd(),outfile_prefix)), exist_ok=True)
    plt.figure(figsize=(12,6))
    plt.subplot(121)
    xy = adata2.obsm['X_gex_2d']
    clusters = np.array(adata2.obs['clusters_gex'])
    cmap = plt.get_cmap('tab20')
    colors = [ cmap.colors[x%20] for x in clusters]
    plt.scatter( xy[:,0], xy[:,1], c=colors)
    plt.title('GEX UMAP colored by GEX clusters')

    plt.subplot(122)
    xy = adata2.obsm['X_tcr_2d']
    clusters = np.array(adata2.obs['clusters_tcr'])
    cmap = plt.get_cmap('tab20')
    colors = [ cmap.colors[x%20] for x in clusters]
    plt.scatter( xy[:,0], xy[:,1], c=colors)
    plt.title('TCR UMAP colored by TCR clusters');

    # these are the nbrhood sizes, as a fraction of the entire dataset:
    #TODO: maybe we should expose these as arguments?
    nbr_fracs = [0.01, 0.1]

    # we use this nbrhood size for computing the nndists
    nbr_frac_for_nndists = 0.01

    all_nbrs, nndists_gex, nndists_tcr = conga.preprocess.calc_nbrs(
        adata2, nbr_fracs, also_calc_nndists=True, nbr_frac_for_nndists=nbr_frac_for_nndists)

    # stash these in obs array, they are used in a few places...
    adata2.obs['nndists_gex'] = nndists_gex
    adata2.obs['nndists_tcr'] = nndists_tcr

    conga.preprocess.setup_tcr_cluster_names(adata2) #stores in adata.uns

    results = conga.correlations.run_graph_vs_graph(
        adata2, all_nbrs, outfile_prefix=outfile_prefix)

    #put the conga hits on top
    conga_scores = adata2.obs['conga_scores']
    colors = np.sqrt(np.maximum(-1*np.log10(conga_scores),0.0))
    reorder = np.argsort(colors)

    plt.figure(figsize=(12,6))
    plt.subplot(121)
    xy = adata2.obsm['X_gex_2d']
    plt.scatter( xy[reorder,0], xy[reorder,1], c=colors[reorder], vmin=0, vmax=np.sqrt(5))
    plt.title('GEX UMAP colored by conga score')

    plt.subplot(122)
    xy = adata2.obsm['X_tcr_2d']
    plt.scatter( xy[reorder,0], xy[reorder,1], c=colors[reorder], vmin=0, vmax=np.sqrt(5))
    plt.title('TCR UMAP colored by conga score');

    conga.correlations.run_graph_vs_features(
        adata2, all_nbrs, outfile_prefix=outfile_prefix)

    conga.plotting.make_graph_vs_features_plots(
        adata2, all_nbrs, outfile_prefix)

    nbrs_gex, nbrs_tcr = all_nbrs[0.1]

    min_cluster_size = 2

    conga.plotting.make_graph_vs_graph_logos(
        adata2,
        outfile_prefix,
        min_cluster_size,
        nbrs_gex,
        nbrs_tcr,
    )

    html_file = outfile_prefix+'_results_summary.html'
    conga.plotting.make_html_summary(adata2, html_file)
    adata2.obs.to_csv(outfile_prefix+"_adata2obs.csv")

