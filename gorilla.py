from scipy import stats
import numpy as np

def getLabelsVector(g, gene_set, ranks):
	return [1 if x in gene_set else 0 for (y, x) in sorted(zip(ranks, list(g.nodes.keys())))]
	
def mHG(labels):
    minHGT = 1
    best_idx = -1
    B = np.sum(labels)
    N = len(labels)
    bn = 0
    for n in range(N):
        bn += labels[n]
        HGT = stats.hypergeom.sf(bn, N, B, n+1) + stats.hypergeom.pmf(bn, N, B, n+1) 
        if minHGT > HGT:
            minHGT = HGT
            best_idx = n
    return minHGT, best_idx