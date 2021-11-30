import scanpy as sc
from os import listdir
from os.path import join
import anndata
import scanpy as sc
import numpy as np
import pandas as pd
from os.path import exists
# Read Morse et all 

def get_budinger():
    print('budinger')
    path_budinger = '../data/budinger/budinger_input_scvi_scanpy_norm.h5ad'
    print(exists(path_budinger), path_budinger)
    assert exists(path_budinger)
    budinger = sc.read_h5ad(path_budinger)

    budinger.obs['cell.type'] = budinger.obs['Cluster'].str[2:]
    budinger.obs['study'] = 'budinger'
    budinger.layers["counts"] = budinger.raw.X.copy()
    budinger.obs['patient.id'] = budinger.obs['Sample Name']
    budinger.obs['disease.status'] = 'covid'

    # MP scoring
    score_marker_genes_ipf(budinger)
    
    return budinger

def get_adams():
    print('adams')

    path_adams = '../data/adams/adams_input_scvi_Macrophage_scanpy_norm.h5ad'
    print(exists(path_adams), path_adams)
    assert exists(path_adams)

    adams = sc.read_h5ad(path_adams)
    adams.layers["counts"] = adams.raw.X.copy()
    adams.obs['study'] = 'adams'
    adams.obs['patient.id'] = adams.obs['Subject_Identity']
    adams.obs['disease.status'] = adams.obs['Disease_Identity']
    # MP scoring
    score_marker_genes_ipf(adams)

    return adams

def get_morse():
    # MORSE
    print('morse')
    rootdir = '../data/morse/by_patient_scanpy_norm_mac'
    ad_all = []
    for f in listdir(rootdir):
        
        next_path_morse = join(rootdir, f)
        print(exists(next_path_morse), next_path_morse)
        assert exists(next_path_morse)
    
        ad_next = sc.read_h5ad(next_path_morse)
        # print(f)
        if ad_next.shape[0] == 0:
            # print(f, 'empty')
            continue        

        ad_next.var_names_make_unique();
        ad_next.raw.var.index = ad_next.var.index
        # del ad_next.raw
        ad_all.append(ad_next)
    morse = anndata.concat(ad_all)
    morse.obs['cell.type'] = morse.obs['cell.type'].str[2:]
    morse.obs['study'] = 'morse'
    morse.var = ad_all[0].var.copy()
    morse.layers["counts"] = morse.raw.X.copy()
    morse.obs['disease.status'] = np.where(morse.obs['sample.id'].str.contains('IPF'), 'IPF', 'control')

    # MP scoring
    score_marker_genes_ipf(morse)

    return morse

def get_bal():
    print('bal')
    # BAL
    path_bal = '../data/bal/bal.h5ad'
    path_bal_feats = '../data/bal/bal_feature_names.tsv.gz'
    
    print(exists(path_bal), path_bal)
    print(exists(path_bal_feats), path_bal)
    
    assert exists(path_bal)
    assert exists(path_bal_feats)
    
    bal = sc.read_h5ad(path_bal)    
    bal_df = pd.read_csv(path_bal_feats, compression='gzip', sep='\t')
    bal.obs['patient.id'] = bal.obs['patient']
    bal.var['ensembl'] = bal.var.index


    bal.var['symbol'] = bal.var.index.map(bal_df.set_index('ensID')['feature'].to_dict())
    # print(sum(pd.isnull(bal.var['symbol'])))
    bal = bal[:,~pd.isnull(bal.var['symbol'])]
    bal.var.index = np.array(bal.var['symbol'])
    bal.var_names_make_unique()
    # print(bal.var.index.value_counts())
    bal.raw = bal.copy()
    bal.layers["counts"] = bal.raw.X.copy()
    bal.obs['study'] = 'BAL'
    bal.obs['disease.status'] = 'BAL (NA)'
    bal.obs['cell.type'] = bal.obs['Celltype_2']
    
    # QC for BAL and only MPs
    bal = bal[bal.obs['percent.mt'] < 10,:]
    bal = bal[bal.obs['nCount_RNA'] > 1000,:]
    bal = bal[bal.obs['nFeature_RNA'] > 1000,:]
    bal = bal[bal.obs['Celltype_2'].str.contains('Macrophages'),:]
    sc.pp.normalize_per_cell(bal, counts_per_cell_after=1e6);
    sc.pp.log1p(bal);
    
    # remove unnecessary keys
    for k in list(bal.obsm):
        if k != 'RNA_MNN_40_UMAP':
            del bal.obsm[k]
    for k in list(bal.obs):
        if 'res' in k and not 'RNA_mnn_40' in k:
            del bal.obs[k]

    print('adding annotation from bal mdm...')
    # add annotation for bal_mdm
    bal_mdm = sc.read_h5ad('../data/bal/bal_mdm.h5ad')
    bal.obs['mdm.type'] = bal.obs.index.map(bal_mdm.obs.Cluster.to_dict())
    
    # MP scoring
    print('scoring genes...')
    score_marker_genes_ipf(bal)
    return bal
    
def get_concatenated_dataset(N=1000):
    # Read Morse et all file by file
    # ADAMS
    adams = get_adams()
    morse = get_morse()
    bal = get_bal()
    budinger = get_budinger()
    
    #HVGs
    sc.pp.highly_variable_genes(adams, n_top_genes=4000)
    sc.pp.highly_variable_genes(morse, n_top_genes=4000)
    sc.pp.highly_variable_genes(bal, n_top_genes=4000)
    sc.pp.highly_variable_genes(budinger, n_top_genes=4000)
    
    a = set(adams.var[adams.var['highly_variable']].index)
    b = set(morse.var[morse.var['highly_variable']].index)
    c = set(bal.var[bal.var['highly_variable']].index)
    d = set(budinger.var[budinger.var['highly_variable']].index)
    
    hvg = set(a).union(b).union(c).union(d)
    
    print(adams.layers['counts'].max())
    print(morse.layers['counts'].max())
    print(bal.layers['counts'].max())
    print(budinger.layers['counts'].max())

    print('morse+adams')
    adata = morse.concatenate(adams)
    print('(morse+adams) + bal')
    adata = adata.concatenate(bal)
    print('(morse+adams+bal) + budinger')
    adata = adata.concatenate(budinger)

    adata.var['highly_variable'] = adata.var.index.isin(hvg)
    print(adata.var['highly_variable'].value_counts())
    
    adata.obs['disease.status'] = adata.obs['disease.status'].str.replace('Control', 'control')
    adata.obs['cell.type'] = np.where(adata.obs['cell.type'].str.contains('nan'),
                                      adata.obs['Subclass_Cell_Identity'], adata.obs['cell.type'])
    
    if not 'cell.type' in adata.obs:
        adata.obs['cell.type'] = adata.obs['Subclass_Cell_Identity']
    
    return adata

def get_marker_genes_ipf():
    marker_genes = {}
    marker_genes['MP.markers'] = {"TREM2", "CD9", "SPP1", "GPNMB", "LGALS3", "LGALS1", "FABP4", "FABP5", "ACP5", "PSAP",
                                     "FTH1", "LIPA", "CTSD", "CTSB", "CSTB", "CTSL", "APOE", "APOC1", "CD63", "LPL"}
    marker_genes['MP.others'] = {"FTL1", "CTSC", "HEXA", "HEXB", "MMP12", "MERTK", "TYROBP", "MFGE8", "DOCK1", "ADGRB1", "IL7R"}
    marker_genes['MP.SPP1'] = {"SPP1"}
    marker_genes['MP.SPP2'] = {"SPP2"}
    marker_genes['MP.all'] = {g for k in marker_genes for g in marker_genes[k]}
    return marker_genes

def score_marker_genes_ipf(adata):
    marker_genes = get_marker_genes_ipf()
    for k in marker_genes:
        marker_genes[k] = marker_genes[k].intersection(set(adata.var.symbol))
    
    for k in ['MP.markers', 'MP.others', 'MP.all']:
        print('scoring', k,  len(marker_genes[k]))
        next_marker_genes = marker_genes[k];
        sc.tl.score_genes(adata, list(next_marker_genes), use_raw=False, score_name=k + ".score");
        
def subset_anndata(adata, N=15000):
    import random
    ad = adata
    if N is not None:
        adata_sel = adata[adata.obs.index.isin(set(random.sample(list(adata.obs.index), N))) | (adata.obs['study'].isin({'BAL', 'morse'})),:]
        # adata_sel = adata[(adata.obs['study'].isin({'BAL', 'morse'})),:]
        ad = adata_sel
    print(ad.shape)
    # print(ad.obs['patient.id'].value_counts())
    print(ad.obs['study'].value_counts())
    return ad

def cartesian(arrays, out=None):
    """
    Generate a cartesian product of input arrays.
    
    Based on discussion from here
    https://stackoverflow.com/questions/1208118/using-numpy-to-build-an-array-of-all-combinations-of-two-arrays
    
    Parameters
    ----------
    arrays : list of array-like
        1-D arrays to form the cartesian product of.
    out : ndarray
        Array to place the cartesian product in.

    Returns
    -------
    out : ndarray
        2-D array of shape (M, len(arrays)) containing cartesian products
        formed of input arrays.

    Examples
    --------
    >>> cartesian(([1, 2, 3], [4, 5], [6, 7]))
    array([[1, 4, 6],
           [1, 4, 7],
           [1, 5, 6],
           [1, 5, 7],
           [2, 4, 6],
           [2, 4, 7],
           [2, 5, 6],
           [2, 5, 7],
           [3, 4, 6],
           [3, 4, 7],
           [3, 5, 6],
           [3, 5, 7]])

    """

    arrays = [np.asarray(x) for x in arrays]
    dtype = arrays[0].dtype

    n = np.prod([x.size for x in arrays])
    if out is None:
        out = np.zeros([n, len(arrays)], dtype=dtype)

    m = int(n / arrays[0].size)
    out[:,0] = np.repeat(arrays[0], m)
    if arrays[1:]:
        cartesian(arrays[1:], out=out[0:m, 1:])
        for j in range(1, arrays[0].size):
            out[j*m:(j+1)*m, 1:] = out[0:m, 1:]
    return out