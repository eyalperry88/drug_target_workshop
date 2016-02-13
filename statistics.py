""" Functions for statistic calculations, for analysis or validations """

from scipy import stats
import propagation
import numpy as np

# Number of iterations used to generate the "healthy" distribution
HEALTHY_ITERATIONS = 1000

def generateHealthyDist(g, genes, mutation_num):
    """ Generates a random distribution for each gene by propagating from
        mutation_num random nodes. """
    gene_dists = {}
    for gene in genes:
        gene_dists[gene] = []
    for i in range(HEALTHY_ITERATIONS):
        scores, iterations = propagation.propagate(g, 'RAND_MT', RANDOM_PRIORS=mutation_num)
        ranks = len(g.nodes) - stats.rankdata(scores)
        for gene in genes:
            gene_dists[gene].append(ranks[g.gene2index[gene]])
    for gene in genes:
        gene_dists[gene] = sorted(gene_dists[gene])
    return gene_dists

def getB2HValue(g, sub_g, g_ranks, sub_g_ranks, gene_dists):
    """ Calculate the B2H score (see paper) """
    sum_b2h = 0
    for gene in gene_dists:
            disease_rank = g_ranks[g.gene2index[gene]]
            new_rank = sub_g_ranks[sub_g.gene2index[gene]]
            disease_q = np.searchsorted(gene_dists[gene], disease_rank) / HEALTHY_ITERATIONS
            new_q = np.searchsorted(gene_dists[gene], new_rank) / HEALTHY_ITERATIONS
            b2h = abs(new_q - disease_q)
            sum_b2h += b2h
    avg_b2h = sum_b2h / len(gene_dists)
    return avg_b2h
                
