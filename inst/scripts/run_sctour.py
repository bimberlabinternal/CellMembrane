import sctour as sct
import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import json

def run_sctour(GEXfile, metafile, exclusion_json_path, ptime_out_file):
    print('scanpy version: ' + sc.__version__)
    print('pandas version: ' + pd.__version__)
    print('sctour version: ' + sct.__version__)

    adataObj = sc.read_10x_h5(GEXfile)
    info = pd.read_csv(metafile, sep=',', index_col=0)
    with open(exclusion_json_path) as f:
        exclusionList = json.load(f)
    cells = adataObj.obs_names.intersection(info.index)
    if len(cells)==0:
        print("No cells in intersection of adata object and metadata supplied in info. Please ensure these objects are correct.")
        exit(1)
    adataObj = adataObj[cells, :]
    adataObj.obs['ClusterNames_0.2'] = info.loc[cells, 'ClusterNames_0.2'].copy()
    adataObj.obs['SubjectId'] = info.loc[cells, 'SubjectId'].copy()

    adataObj.X = round(adataObj.X).astype(np.float32)
    adataObj.obs['Population'] = info.loc[cells, 'Population'].copy()
    adataObj.obs['TandNK_ActivationCore_UCell'] = info.loc[cells, 'TandNK_ActivationCore_UCell'].copy()
    adataObj.obs['Timepoint'] = info.loc[cells, 'Timepoint'].copy()
    sc.pp.calculate_qc_metrics(adataObj, percent_top=None, log1p=False, inplace=True)
    sc.pp.filter_genes(adataObj, min_cells=20)
    sc.pp.highly_variable_genes(adataObj, flavor='seurat_v3', n_top_genes=2000, subset=True, inplace=False)
    
    adataObj = adataObj[:, list(set(adataObj.var_names) - set(exclusionList))].copy()

    # Added to avoid: 'SparseCSRView' object has no attribute 'A' error. See: https://github.com/LiQian-XC/sctour/issues/10
    if adataObj.is_view:
        print('AnnData object is a view, converting')
        adataObj = adataObj.copy()

    if adataObj.X.is_view:
        print('AnnData.X object is a view, converting')
        adataObj.X = adataObj.X.copy()

    print('Anndata object:')
    print(adataObj)
    print(adataObj.X)

    tnode = sct.train.Trainer(adataObj)
    tnode.train()
    adataObj.obs['ptime'] = tnode.get_time()
    mix_zs, zs, pred_zs = tnode.get_latentsp(alpha_z=0.2, alpha_predz=0.8)
    adataObj.obsm['X_TNODE'] = mix_zs
    adataObj = adataObj[np.argsort(adataObj.obs['ptime'].values), :]
    sc.pp.neighbors(adataObj, use_rep='X_TNODE', n_neighbors=15)
    sc.tl.umap(adataObj, min_dist=0.1)
    adataObj.obs['ptime'] = sct.train.reverse_time(adataObj.obs['ptime'].values)
    adataObj.obsm['X_VF'] = tnode.get_vector_field(adataObj.obs['ptime'].values, adataObj.obsm['X_TNODE'])
    
    fig, axs = plt.subplots(ncols=4, nrows=1, figsize=(15, 3))
    adataObj.obs = adataObj.obs.astype({'ClusterNames_0.2':'category', 'SubjectId':'category'})
    sc.pl.umap(adataObj, color='Timepoint', size=20, ax=axs[0], show=False)
    sc.pl.umap(adataObj, color='ptime', size=20, ax=axs[1], show=False)
    sc.pl.umap(adataObj, color='TandNK_ActivationCore_UCell', size=20, ax=axs[2], show=False)
    sc.pl.umap(adataObj, color='SubjectId', size=20, ax=axs[3], show=False)
    plt.show()

    df = pd.DataFrame(adataObj.obs['ptime'])
    df.to_csv(ptime_out_file)

    return adataObj
