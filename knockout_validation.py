from statistics import *
from parsing_utils import *

f = open('diff_per_gene_sorted.txt', 'r')
all_genes = []
for line in f:
    gene = line.split('\t')[0]
    all_genes.append(gene)
f.close()

g = DTWGraph()
loadPPIData("data/PPI_HIPPIE.txt", g)
causal_genes = loadCausalGenes("data/cancer_cosmic_genes.txt", g)
labels = [1 if x in causal_genes else 0 for x in all_genes]
p, mHG_idx = mHG(labels)
print('mHG', p, mHG_idx)

gene_num = len(all_genes)
k = round(gene_num / 10) # using top 10 percent
causal_gene_hits = 0
above_threshold_genes = 0
i = 0
causal_genes_in_all_genes = 0
for gene in all_genes:
    if gene in causal_genes:
        causal_genes_in_all_genes += 1
        if i < k:
            causal_gene_hits += 1
    i += 1
    
print('hit ' + str(causal_gene_hits) + ' out of ' + str(causal_genes_in_all_genes) + '(above threshold: ' + str(above_threshold_genes) + ')')
print('score: ' + str(float(causal_gene_hits) / len(causal_genes)))
p = stats.hypergeom.sf(causal_gene_hits, len(all_genes), causal_genes_in_all_genes, k) + stats.hypergeom.pmf(causal_gene_hits, len(all_genes), causal_genes_in_all_genes, k)
print('p-value: ' + str(p))
print('-log(p-value): ' + str(-np.log10(p)))
