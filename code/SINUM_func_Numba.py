import numpy as np
import pandas as pd
from sklearn.preprocessing import KBinsDiscretizer
import math, sys, os, time
from scipy import sparse
from numba import njit

@njit
def get_upper_and_lower_file(data, boxsize):
    geneNum, cellNum = data.shape
    upper = np.zeros((geneNum, cellNum))
    lower = np.zeros((geneNum, cellNum))
    for i in range(geneNum):
        s2 = np.argsort(data[i, :])
        s1 = np.array(sorted(data[i, :]))
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

@njit
def get_entropyX_matrix(Xt,geneTot,cellTot,binsNum):
    freq = np.zeros((geneTot,binsNum))
    for r in range(geneTot):
        bin_id = Xt[r,:]
        # unique, counts = np.unique(bin_id, return_counts=True)
        # freq[r,unique.astype('int32')] = counts
        for i in range(cellTot):
            freq[r, bin_id[i]] += 1
    prob = freq/cellTot
    # EntropyX = prob*np.log2(prob, where=prob != 0)
    # EntropyX = np.round(EntropyX, 6)
    # EntropyX[np.isnan(EntropyX)] = 0
    EntropyX = np.where(prob != 0, prob * np.log2(prob), 0)
    return np.round(EntropyX, 6)

@njit
def get_entropyXY_matrix(Xt, freq, g1_site, g2_site, cellTot):
    freq.fill(0)
    for i in range(cellTot):
        freq[Xt[g1_site,i],Xt[g2_site,i]] += 1
    prob = freq/cellTot
    # entropies = prob*np.log2(prob, where=prob!=0)
    # entropies = np.round(entropies,6)
    # entropies[np.isnan(entropies)] = 0
    entropies = np.where(prob != 0, prob * np.log2(prob), 0)
    return np.round(entropies, 6)

@njit
def compute_mutual_info_matrix(cell_k, Xt, EntropyX, freq, data, Xlow, Xupp, cellTot, geneTot):
    MI_matrix = np.zeros((geneTot, geneTot))
    for g1 in range(geneTot):
        for g2 in range(g1 + 1, geneTot):
            if data[g1,cell_k]*data[g2,cell_k]==0:
                continue # mi = 0 if any of the gene in target cell does not express
            entropies = get_entropyXY_matrix(Xt, freq, g1, g2, cellTot)
            start1 = Xlow[g1, cell_k]
            end1 = Xupp[g1, cell_k]
            start2 = Xlow[g2, cell_k]
            end2 = Xupp[g2, cell_k]
            HX = -np.sum(EntropyX[g1, start1:end1])
            HY = -np.sum(EntropyX[g2, start2:end2])
            HXY = -np.sum(entropies[start1:end1, start2:end2])
            mi = np.round(HX+HY-HXY, 6)
            MI_matrix[g1,g2] = mi
    
    return MI_matrix

def sinumnet(data, boxsize, JobStart, JobEnd, savedir):
    """
    :param data: GEM (gene * cell matrix)
    :param save: ouput file name
    :param start: Count from the xth gene
    :param stop: Stop counting until the yth gene
    :param boxsize: box size
    """
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
    
    genePairnum=math.comb(geneTot,2)
    
    for cell_k in range(JobStart-1, JobEnd):
        MI_matrix = compute_mutual_info_matrix(cell_k, Xt, EntropyX, freq, data, Xlow, Xupp, cellTot, geneTot, genePairnum)
        MI_matrix = np.round(MI_matrix, 6)
        sparse_matrix = sparse.csc_matrix(MI_matrix)
        savepath = os.path.join(savedir, f'sinum_c{cell_k+1}.MI.npz')
        sparse.save_npz(savepath, sparse_matrix)

# if __name__ == "__main__":
#     dataSet = sys.argv[1]
#     cellJobListRange = sys.argv[2]#1-100
#     boxsize = 0.2
#     jobstart = int(cellJobListRange.split('-')[0])
#     jobend = int(cellJobListRange.split('-')[1])
#     fmt = "count"
#     GEMpath = f'/data/WenYao/locCSN/scRNAseq_datasets/{dataSet}/{dataSet}_{fmt}_scRNA_rmZero_log2.txt'
#     GEM = pd.read_csv(GEMpath, sep='\t', index_col=0)
#     data=GEM.values
#     del(GEM)
#     savedir = f"./{dataSet}/net/"
#     print(f"Cell range: {jobstart}-{jobend}")
#     print(f"[{time.ctime()}] Start!")
#     start = time.time()
#     sinumnet(data, boxsize, jobstart, jobend, savedir)
#     end = time.time()
#     print(f"[{time.ctime()}] End!")
#     print(f"Time cost: {end-start:.3f} secs.")
