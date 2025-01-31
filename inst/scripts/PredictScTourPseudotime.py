import sctour as sct
import scanpy as sc
import numpy as np
import pandas as pd
import torch
import anndata
from anndata import AnnData

torch.serialization.add_safe_globals([AnnData, anndata._core.file_backing.AnnDataFileManager, np.core.multiarray._reconstruct, np.ndarray, np.dtype, np])

def PredictPseudotime(GEXfile, model_file, ptime_out_file, embedding_out_file):
    #read count data and variable genes
    adataObj = sc.read_10x_h5(GEXfile)
    adataObj.X = round(adataObj.X).astype(np.float32)
    
    # Added to avoid: 'SparseCSRView' object has no attribute 'A' error. See: https://github.com/LiQian-XC/sctour/issues/10
    if adataObj.X.getformat() == 'csr':
        print('AnnData object is a csr matrix, converting to dense because scipy depreciated the .A shorthand')
        adataObj.X = adataObj.X.toarray()

    checkpoint = torch.load(model_file, map_location=torch.device('cpu'), weights_only = False)
    model_adata = checkpoint['adata']

    genes_in_model = model_adata.var.index.values.tolist()

    #basic preprocessing
    sc.pp.calculate_qc_metrics(adataObj, percent_top=None, log1p=False, inplace=True)
    sc.pp.filter_genes(adataObj, min_cells=20)
    sc.pp.highly_variable_genes(adataObj, flavor='seurat_v3', n_top_genes=2000, subset=True, inplace=False)
    
    #subset to genes found in the pretrained model.
    adataObj = adataObj[:, genes_in_model]
    #initalize a trainer and pull a previously saved model from model_file

    tnode = sct.predict.load_model(model_file)
    pred_t = sct.predict.predict_time(adata = adataObj, model = tnode)
    adataObj.obs['ptime'] = pred_t
    mix_zs, zs, pred_zs  = sct.predict.predict_latentsp(adata = adataObj, model = tnode)
    adataObj.obsm['X_TNODE'] = mix_zs
    np.savetxt(embedding_out_file, adataObj.obsm['X_TNODE'] , delimiter=",")

    df = pd.DataFrame(adataObj.obs['ptime'])
    df.to_csv(ptime_out_file)

    return adataObj
