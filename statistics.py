from scipy import stats, misc
from propagation import *
import numpy as np
import bisect

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

def mHG_pval(p, labels):
    B = np.sum(labels)
    N = len(labels)
    W = N - B
    """
    hg_mat = np.zeros(shape=(N + 1, B + 1))
    hg_mat[0, 0] = 1
    for b in range(1, B + 1):
        mat[0, b] = 0
    for n in range(0, N):
    """     
    mat = np.zeros(shape=(W + 1, B + 1))
    mat[0, 0] = 1
    for w in range(0, W + 1):
        if w % 100 == 0:
            print(w)
        for b in range(0, B + 1):
            HGT = stats.hypergeom.sf(b, N, B, w + b) + stats.hypergeom.pmf(b, N, B, w + b)
            if HGT <= p:
                mat[w, b] = 0
            else:
                if b == 0:
                    mat[w, b] = mat[w - 1, b] * ((W - w + 1) / (B + W - b - w + 1))
                elif w == 0:
                    mat[w, b] = ((W - w + 1) / (B + W - b - w + 1)) + mat[w, b - 1] * ((B - b + 1) / (B + W - b - w + 1))
                else:
                    mat[w, b] = mat[w - 1, b] * ((W - w + 1) / (B + W - b - w + 1)) + mat[w, b - 1] * ((B - b + 1) / (B + W - b - w + 1))
    return mat[W, B]

def getDiffValue(g, sub_g, g_ranks, sub_g_ranks, genes):
    dist = 0
    for i in range(len(genes)):
            dist += (g_ranks[g.gene2index[genes[i]]] - sub_g_ranks[sub_g.gene2index[genes[i]]]) ** 2
    return dist

def generateHealthyDist(g, genes, mutation_num):
    gene_dists = {}
    for gene in genes:
        gene_dists[gene] = []
    for i in range(1000):
        scores, iterations = propagation.propagate(g, 'RAND_MT', RANDOM_PRIORS=mutation_num)
        ranks = len(g.nodes) - scores.argsort().argsort()
        for gene in genes:
            gene_dists[gene].append(ranks[g.gene2index[gene]])
    for gene in genes:
        gene_dists[gene] = sorted(gene_dists[gene])
    return gene_dists

def getB2HValue(g, sub_g, g_ranks, sub_g_ranks, gene_dists):
        for gene in gene_dists:
                disease_rank = g_ranks[g.gene2index[gene]]
                new_rank = sub_g_ranks[g.gene2index[gene]]
                disease_q = bisect.bisect(disease_rank, gene_dists[gene]) / len(gene_dists[gene])
                new_q = bisect.bisect(new_rank, gene_dists[gene]) / len(gene_dists[gene])
                print('gene', gene)
                print('old rank', disease_rank, 'new rank', new_rank)
                print('old q', disease_q, 'new q', new_q)
                b2h = abs(new_q - disease_q)
                
def mHG_pval_old(p, labels):
    B = np.sum(labels)
    N = len(labels)
    mat = np.zeros(shape=(N + 1, B + 1))
    mat[0, 0] = 1
    for n in range(0, B):
        mat[n, n + 1] = 0
    for n in range(1, N + 1):
        for b in range(max(B - N + n, 0), min(B, n) + 1):
            HGT = stats.hypergeom.sf(b, N, B, n) + stats.hypergeom.pmf(b, N, B, n)
            if HGT <= p:
                mat[n, b] = 0
            else:
                if b == 0:
                    mat[n, b] = mat[n - 1, b]
                else:
                    mat[n, b] = mat[n - 1, b] + mat[n - 1, b - 1]
    NchooseB = misc.comb(N, B)
    return 1 - mat[N,B]/NchooseB
