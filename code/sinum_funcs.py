# -*- coding: utf-8 -*-
"""
Created on Thu Aug 18 21:21:04 2022

@author: TingYi
"""

import scipy.sparse
import math
import numpy as np
import os
import pandas as pd
from sklearn.preprocessing import KBinsDiscretizer
from itertools import combinations
import warnings
warnings.filterwarnings("ignore", category=RuntimeWarning)

def mkdir(addDir):
    if not os.path.exists(addDir):
        os.makedirs(addDir)
def run_normalize(ndm,data):
    a = np.mean(np.sum(np.sign(data), axis=0))
    c = 2000
    sumOfCol = np.sum(ndm, axis=0)
    ndm = ndm / sumOfCol * a ** 2 / c
    ndm = np.round(ndm, 6)
    ndm[np.isnan(ndm)] = 0 # deal with RuntimeWarning: invalid value encountered in true_divide
    return ndm
def get_upper_and_lower_file(data, boxsize):
    geneNum, cellNum = data.shape
    upper = np.zeros((geneNum, cellNum))
    lower = np.zeros((geneNum, cellNum))
    for i in range(geneNum):
        s2 = np.argsort(data[i, :])
        s1 = np.array(sorted(data[i, :]))
        # sign_c = np.sum(np.where(s1 == 0, 0, 1))
        sign_c = np.sum(s1>0)
        n0 = cellNum - sign_c
        h=int(np.round(boxsize / 2 * sign_c))
        k = 1
        while k <= cellNum:
            s = 0
            while k + s + 1 <= cellNum and s1[(k - 1) + s + 1] == s1[(k - 1)]:
                s += 1
            if s >= h:
                upper[i, s2[k - 1:k + s]] = data[i, s2[k - 1]]
                lower[i, s2[k - 1:k + s]] = data[i, s2[k - 1]]
            else:
                upper[i, s2[k - 1:k + s]] = data[i, s2[min(cellNum - 1, (k - 1) + s + h)]]
                lower[i, s2[k - 1:k + s]] = data[i, s2[max(n0 * (n0 > h) + 1 - 1, (k - 1) - h)]]
            k = k + s + 1
    return (upper,lower)
#%%
def convert_zScore(mutual_info, geneNum):
    mu = np.average(mutual_info, axis=1).reshape(geneNum, 1)
    std = np.std(mutual_info, axis=1).reshape(geneNum, 1)
    mutual_info = (mutual_info-mu) / std
    mutual_info[np.isnan(mutual_info)] = 0
    mutual_info = np.round(mutual_info, 2)
    return mutual_info
def get_entropyX_matrix(Xt,geneTot,cellTot,binsNum):
    freq = np.zeros((geneTot,binsNum))
    for r in range(geneTot):
        bin_id = Xt[r,:]
        unique, counts = np.unique(bin_id, return_counts=True)
        freq[r,unique.astype('int32')] = counts
    prob = freq/cellTot
    EntropyX = prob*np.log2(prob, where=prob != 0)
    EntropyX = np.round(EntropyX, 6)
    EntropyX[np.isnan(EntropyX)] = 0
    return EntropyX
def get_entropyXY_matrix(Xt, freq, g1_site, g2_site, cellTot):
    freq.fill(0)
    for i in range(cellTot):
        freq[Xt[g1_site,i],Xt[g2_site,i]] += 1
    prob = freq/cellTot
    entropies = prob*np.log2(prob, where=prob!=0)
    entropies = np.round(entropies,6)
    entropies[np.isnan(entropies)] = 0
    return entropies
def symmetrize(a):
    return a + a.T - np.diag(a.diagonal())
def sinumnet(GEM, outdir, save, boxsize, cutoff):
    """
    :param data: GEM (gene * cell matrix)
    :param save: ouput file name
    :param start: Count from the xth gene
    :param stop: Stop counting until the yth gene
    :param boxsize: box size
    """
    geneLst=GEM.index
    cellLst=GEM.columns
    data=GEM.values
    del(GEM)
    geneTot, cellTot = data.shape
    binsNum = int(np.round(cellTot**0.5))
    est = KBinsDiscretizer(n_bins=binsNum, encode='ordinal', strategy='uniform').fit(data.T)
    Xt = est.transform(data.T).T.astype('int32')
    
    upper, lower = get_upper_and_lower_file(data, boxsize)
    Xupp = est.transform(upper.T).T.astype('int32')
    Xupp = Xupp+1
    Xlow = est.transform(lower.T).T.astype('int32')

    EntropyX = get_entropyX_matrix(Xt,geneTot,cellTot,binsNum)
    
    freq = np.zeros((binsNum,binsNum))
    # signal = np.where(data, 1, 0)
    
    genePairnum=math.comb(geneTot,2)
    for s, genePair in enumerate(combinations(range(geneTot),2)):
        if s%500==0:print(f'{s}/{genePairnum}')
        g1=genePair[0]
        g2=genePair[1]
        mutual_info = np.zeros((genePairnum, cellTot)) # genepair vs cell matrix
        entropies = get_entropyXY_matrix(Xt, freq, g1, g2, cellTot)
        for cell in range(cellTot):
            if data[g1,cell]*data[g2,cell]==0:
                mi=0 # mi = 0 if any of the gene in target cell does not express
                continue
            start1 = Xlow[g1, cell]
            end1 = Xupp[g1, cell]
            start2 = Xlow[g2, cell]
            end2 = Xupp[g2, cell]
            HX = -np.sum(EntropyX[g1, start1:end1])
            HY = -np.sum(EntropyX[g2, start2:end2])
            HXY = -np.sum(entropies[start1:end1, start2:end2])
            mi = np.round(HX+HY-HXY, 6)
            mutual_info[s, cell] = mi
    zscore = convert_zScore(mutual_info, genePairnum)
    
    print('Constructing adjacency matrix and degree matrix')
    dm=np.zeros((geneTot, cellTot)) # initialize degree matrix
    for cell in range(cellTot):
        print(f'start cell {cell+1}...')
        adj_matrix=np.zeros((geneTot, geneTot)) #initialize adjacency matrix
        adj_matrix[np.triu_indices(geneTot, 1)] = np.where(zscore[:,cell]>cutoff,1,0)

        adj_matrix=symmetrize(adj_matrix)
        sparse_matrix=scipy.sparse.csc_matrix(adj_matrix)
        mkdir(f'{outdir}/output_adj_network')
        scipy.sparse.save_npz(f'{outdir}/output_adj_network/{save}_adj_network.npz', sparse_matrix)
        ## save as degree matrix
        dm[:,cell]=np.sum(adj_matrix, axis=1)
    ndm=run_normalize(dm, data)
    ndm_df=pd.DataFrame(ndm,index=geneLst,columns=cellLst)
    mkdir(f'{outdir}/output_degree_matrix')
    ndm_df.to_csv(f'{outdir}/output_degree_matrix/{save}_degree_matrix.csv')
    
    #%%
    # dm=np.zeros((geneTot, cellTot)) # initialize degree matrix
    # for cell in range(cellTot):
    #     print(f'Constructing cell {cell+1}...')
    #     MI_matrix=np.zeros((geneTot, geneTot)) # initialize network matrix
    #     for genePair in combinations(range(geneTot),2):
    #         g1=genePair[0]
    #         g2=genePair[1]
            
            
    #         MI_matrix[g1,g2]=mi
    #     ## save as unweighted adjacency matrix
    #     adj_matrix=np.where(MI_matrix > cutoff, 1, 0)
    #     adj_matrix=symmetrize(adj_matrix)
    #     sparse_matrix=scipy.sparse.csc_matrix(adj_matrix)
    #     mkdir(f'{outdir}/output_adj_network')
    #     scipy.sparse.save_npz(f'{outdir}/output_adj_network/{save}_adj_network.npz', sparse_matrix)
    #     ## save as degree matrix
    #     dm[:,cell]=np.sum(adj_matrix, axis=1)
    # dm_df=pd.DataFrame(dm,index=GEM.index,columns=GEM.columns)
    # mkdir(f'{outdir}/output_degree_matrix')
    # dm_df.to_csv(f'{outdir}/output_degree_matrix/{save}_degree_matrix.csv')
    # #%%
        
    #     geneNum = geneTot-g1_site-1
    #     mutual_info = np.zeros((geneNum, cellTot))
    #     for s, g2_site in enumerate(range(g1_site+1, geneTot)):
    #         entropies = get_entropyXY_matrix(Xt, freq, g1_site, g2_site, cellTot)
    #         for cell_site in range(cellTot):
    #             start1 = Xlow[g1_site, cell_site]
    #             end1 = Xupp[g1_site, cell_site]
    #             start2 = Xlow[g2_site, cell_site]
    #             end2 = Xupp[g2_site, cell_site]
    #             HX = -np.sum(EntropyX[g1_site, start1:end1])
    #             HY = -np.sum(EntropyX[g2_site, start2:end2])
    #             HXY = -np.sum(entropies[start1:end1, start2:end2])
    #             mi = round(HX+HY-HXY, 6)
    #             mutual_info[s, cell_site] = mi

    #     mutual_info = convert_zScore(mutual_info, geneNum)

    #     mutual_info = convert_P(mutual_info, signal, g1_site, geneNum, geneTot)

    #     save_pmin_zScored(save, mutual_info, g1_site)