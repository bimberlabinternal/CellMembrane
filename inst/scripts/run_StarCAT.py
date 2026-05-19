import os
import importlib.metadata
import sys

import numpy as np
import scanpy as sc

from starcat import starCAT


def run_StarCAT(gex_datafile, reference, output_dir, name, cachedir,
                print_versions=True):
    if print_versions:
        print('scanpy version: ' + importlib.metadata.version('scanpy'))
        print('starcatpy version: ' + importlib.metadata.version('starcatpy'))
        print('anndata version: ' + importlib.metadata.version('anndata'))
        print('numpy version: ' + importlib.metadata.version('numpy'))
        print('Python version: ' + sys.version)

    os.makedirs(output_dir, exist_ok=True)

    adata = sc.read_10x_h5(gex_datafile)
    adata.var_names_make_unique()

    adata.X = np.round(adata.X).astype(np.float32)

    print('Total cells passed to starCAT: ' + str(adata.shape[0]))
    print('Total genes passed to starCAT: ' + str(adata.shape[1]))

    cat = starCAT(reference=reference, cachedir=cachedir)
    usage_norm, scores = cat.fit_transform(adata)
    cat.save_results(output_dir, name)

    print('starCAT outputs written to ' + output_dir)
