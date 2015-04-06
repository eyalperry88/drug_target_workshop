from parsing_utils import *
from propagation import *
import numpy
from scipy import stats

"""
patients = ["TCGA-AB-2803-03",
"TCGA-AB-2805-03",
"TCGA-AB-2806-03",
"TCGA-AB-2807-03",
"TCGA-AB-2810-03",
"TCGA-AB-2811-03",
"TCGA-AB-2812-03",
"TCGA-AB-2813-03",
"TCGA-AB-2814-03",
"TCGA-AB-2815-03"]
"""
patients = [
"TCGA-AB-2980-03",
"TCGA-AB-2981-03",
"TCGA-AB-2982-03",
"TCGA-AB-2983-03",
"TCGA-AB-2984-03",
"TCGA-AB-2985-03",
"TCGA-AB-2986-03",
"TCGA-AB-2987-03",
"TCGA-AB-2988-03",
"TCGA-AB-2990-03",
"TCGA-AB-2991-03",
"TCGA-AB-2992-03",
"TCGA-AB-2993-03",
"TCGA-AB-2994-03",
"TCGA-AB-2995-03",
"TCGA-AB-2996-03",
"TCGA-AB-2998-03",
"TCGA-AB-2999-03",
"TCGA-AB-3000-03",
"TCGA-AB-3001-03",
"TCGA-AB-3002-03",
"TCGA-AB-3005-03",
"TCGA-AB-3006-03",
"TCGA-AB-3007-03",
"TCGA-AB-3008-03",
"TCGA-AB-3009-03",
"TCGA-AB-3011-03",
"TCGA-AB-3012-03",
"TCGA-AB-2921-03",
"TCGA-AB-2882-03",
"TCGA-AB-2858-03",
"TCGA-AB-2927-03",
"TCGA-AB-2909-03",
"TCGA-AB-2808-03",
"TCGA-AB-2946-03",
"TCGA-AB-2944-03",
"TCGA-AB-2853-03",
"TCGA-AB-2828-03",
"TCGA-AB-2822-03",
"TCGA-AB-2935-03"]

g = DTWGraph()
print('Loading PPI...')
loadPPIData("data/PPI_HIPPIE.txt", g)
print('Normalizing weights...')
g.normalizeWeights()
gene_num = len(g.nodes)
results_avg = {}
results_min = {}

for patient in patients:
    print("Working on patient " + patient)
    g.initGraph()
    print('Loading expression...')
    ge_loaded = loadExpressionData("data/AML_Expression.txt", g, patient)
    if (ge_loaded == 0):
        print('no expression data')
        continue
    print('Loading mutations...')
    mt_loaded = loadMutationData("data/AML_Mutations.txt", g, patient)
    if (mt_loaded == 0):
        print('no mutation data')        
        continue
    print('Propagating from expression...')
    GEscores = propagate(g, 'GE')
    GE_average = numpy.average(GEscores.values())
    GE_std = numpy.std(GEscores.values())
    print('Before score normalization - mean: ' + str(GE_average) + ', std: ' + str(GE_std))
    GEscores.update((key, (val - GE_average) / GE_std)  for key, val in GEscores.items())
    print('Propagating from mutation...')
    MTscores = propagate(g, 'MT')
    MT_average = numpy.average(MTscores.values())
    MT_std = numpy.std(MTscores.values())
    print('Before score normalization - mean: ' + str(MT_average) + ', std: ' + str(MT_std))
    MTscores.update((key, (val - MT_average) / MT_std)  for key, val in MTscores.items())
    
    # u, prob = stats.mannwhitneyu(GEscores.values(), MTscores.values())
    # print('MWU prob is ' + str(prob))
    
    results_avg[patient] = {}
    results_min[patient] = {}
    for gene in GEscores:
        results_avg[patient][gene] = (GEscores[gene] + MTscores[gene]) / 2
        results_min[patient][gene] = min(GEscores[gene], MTscores[gene])

    """
    import operator
    joint_scores = {}
    for gene in GEscores:
        joint_scores[gene] = (GEscores[gene] + MTscores[gene]) / 2
    sorted_scores = sorted(joint_scores.items(), key=operator.itemgetter(1), reverse=True)
    print("Top 10: ")
    print(sorted_scores[0:10])
    """
actual_patients = results_avg.keys()
for i in range(0, len(actual_patients) - 1):
    for j in range(i + 1, len(actual_patients)):
        print('Comparing patient ' + actual_patients[i] + ' vs. ' + actual_patients[j])
        T, p1 = stats.wilcoxon(results_avg[actual_patients[i]].values(), results_avg[actual_patients[j]].values())
        print('Wilcoxon prob for avg is ' + str(p1))
        
        T, p2 = stats.wilcoxon(results_min[actual_patients[i]].values(), results_min[actual_patients[j]].values())
        print('Wilcoxon prob for min is ' + str(p2))
        

# using top 10 percent
k = gene_num / 10
expected = {}
observed = {}
for gene in g.nodes:
    expected[gene] = float(k * len(actual_patients)) / gene_num
    observed[gene] = 0
for i in range(0, len(actual_patients)):
    sorted_scores = sorted(results_min[actual_patients[i]].items(), key=operator.itemgetter(1), reverse=True)
    for gene, score in sorted_scores[0:k]:
        observed[gene] += 1
chi = 0
threshold = 5
f_causal = open("data/cancer_census_genes.txt", 'r')
causal_genes = []
for line in f_causal:
    causal_genes.append(line.strip())
f_causal.close()
causal_gene_hits = 0
for gene in observed:
    if observed[gene] > threshold:
        if gene in causal_genes:
            print('hit ' + gene)
            causal_gene_hits += 1
        # print(gene + ' - ' + str(observed[gene]))
    chi += (float(observed[gene] - expected[gene]) ** 2) / expected[gene]
print('chi statistic: ' + str(chi))

print('hit ' + str(causal_gene_hits) + ' out of ' + str(len(causal_genes)))
print('score: ' + str(float(causal_gene_hits) / len(causal_genes)))
