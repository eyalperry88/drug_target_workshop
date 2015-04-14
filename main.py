from parsing_utils import *
import propagation
import numpy
from scipy import stats


# top 10:
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
# bottom 40
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
"""

# top 100
# patients = ["TCGA-AB-2803-03", "TCGA-AB-2805-03", "TCGA-AB-2806-03", "TCGA-AB-2807-03", "TCGA-AB-2810-03", "TCGA-AB-2811-03", "TCGA-AB-2812-03", "TCGA-AB-2813-03", "TCGA-AB-2814-03", "TCGA-AB-2815-03", "TCGA-AB-2816-03", "TCGA-AB-2817-03", "TCGA-AB-2818-03", "TCGA-AB-2819-03", "TCGA-AB-2820-03", "TCGA-AB-2821-03", "TCGA-AB-2823-03", "TCGA-AB-2824-03", "TCGA-AB-2825-03", "TCGA-AB-2826-03", "TCGA-AB-2830-03", "TCGA-AB-2832-03", "TCGA-AB-2833-03", "TCGA-AB-2834-03", "TCGA-AB-2835-03", "TCGA-AB-2836-03", "TCGA-AB-2837-03", "TCGA-AB-2838-03", "TCGA-AB-2839-03", "TCGA-AB-2840-03", "TCGA-AB-2841-03", "TCGA-AB-2842-03", "TCGA-AB-2843-03", "TCGA-AB-2844-03", "TCGA-AB-2845-03", "TCGA-AB-2846-03", "TCGA-AB-2847-03", "TCGA-AB-2848-03", "TCGA-AB-2849-03", "TCGA-AB-2851-03", "TCGA-AB-2854-03", "TCGA-AB-2855-03", "TCGA-AB-2856-03", "TCGA-AB-2857-03", "TCGA-AB-2859-03", "TCGA-AB-2860-03", "TCGA-AB-2861-03", "TCGA-AB-2862-03", "TCGA-AB-2863-03", "TCGA-AB-2865-03", "TCGA-AB-2866-03", "TCGA-AB-2867-03", "TCGA-AB-2868-03", "TCGA-AB-2869-03", "TCGA-AB-2870-03", "TCGA-AB-2871-03", "TCGA-AB-2872-03", "TCGA-AB-2873-03", "TCGA-AB-2874-03", "TCGA-AB-2875-03", "TCGA-AB-2877-03", "TCGA-AB-2879-03", "TCGA-AB-2880-03", "TCGA-AB-2881-03", "TCGA-AB-2884-03", "TCGA-AB-2885-03", "TCGA-AB-2886-03", "TCGA-AB-2887-03", "TCGA-AB-2888-03", "TCGA-AB-2889-03", "TCGA-AB-2890-03", "TCGA-AB-2891-03", "TCGA-AB-2895-03", "TCGA-AB-2896-03", "TCGA-AB-2897-03", "TCGA-AB-2898-03", "TCGA-AB-2899-03", "TCGA-AB-2900-03", "TCGA-AB-2901-03", "TCGA-AB-2903-03", "TCGA-AB-2904-03", "TCGA-AB-2908-03", "TCGA-AB-2910-03", "TCGA-AB-2911-03", "TCGA-AB-2912-03", "TCGA-AB-2913-03", "TCGA-AB-2914-03", "TCGA-AB-2915-03", "TCGA-AB-2916-03", "TCGA-AB-2917-03", "TCGA-AB-2918-03", "TCGA-AB-2919-03", "TCGA-AB-2920-03", "TCGA-AB-2924-03", "TCGA-AB-2925-03", "TCGA-AB-2928-03", "TCGA-AB-2929-03", "TCGA-AB-2930-03", "TCGA-AB-2931-03", "TCGA-AB-2932-03"]


g = DTWGraph()
print('Loading PPI...')
loadPPIData("data/PPI_HIPPIE.txt", g)
print('Normalizing weights...')
g.normalizeWeights()
gene_num = len(g.nodes)
results_avg = {}
results_min = {}

for patient in patients:
    print("Working on patient", patient)
    if patient in results_min:
        print("patient already exists - skipping")
        continue;
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
    GEscores = propagation.propagate(g, 'GE')
    GE_average = numpy.average(list(GEscores.values()))
    GE_std = numpy.std(list(GEscores.values()))
    print('Before score normalization - mean:', GE_average, ', std:', GE_std)
    GEscores.update((key, (val - GE_average) / GE_std)  for key, val in GEscores.items())
    print('Propagating from mutation...')
    MTscores = propagation.propagate(g, 'MT')
    MT_average = numpy.average(list(MTscores.values()))
    MT_std = numpy.std(list(MTscores.values()))
    print('Before score normalization - mean:', MT_average, ', std:', MT_std)
    MTscores.update((key, (val - MT_average) / MT_std)  for key, val in MTscores.items())
    
    # u, prob = stats.mannwhitneyu(GEscores.values(), MTscores.values())
    # print('MWU prob is ' + str(prob))
    
    # results_avg[patient] = {}
    results_min[patient] = {}
    for gene in GEscores:
        # results_avg[patient][gene] = (GEscores[gene] + MTscores[gene]) / 2
        results_min[patient][gene] = min(GEscores[gene], MTscores[gene])

actual_patients = list(results_min.keys())
"""

for i in range(0, len(actual_patients) - 1):
    for j in range(i + 1, len(actual_patients)):
        print('Comparing patient ' + actual_patients[i] + ' vs. ' + actual_patients[j])
        T, p1 = stats.wilcoxon(results_avg[actual_patients[i]].values(), results_avg[actual_patients[j]].values())
        print('Wilcoxon prob for avg is ' + str(p1))
        
        T, p2 = stats.wilcoxon(results_min[actual_patients[i]].values(), results_min[actual_patients[j]].values())
        print('Wilcoxon prob for min is ' + str(p2))
"""       

# using top 10 percent
import operator
k = round(gene_num / 10)
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
threshold = len(actual_patients) / 2
f_causal = open("data/cancer_census_genes.txt", 'r')
causal_genes = []
for line in f_causal:
    causal_genes.append(line.strip())
f_causal.close()
causal_gene_hits = 0
for gene in observed:
    if observed[gene] > threshold:
        if gene in causal_genes:
            causal_gene_hits += 1
        # print(gene + ' - ' + str(observed[gene]))
    chi += (float(observed[gene] - expected[gene]) ** 2) / expected[gene]
print('chi statistic: ' + str(chi))

print('hit ' + str(causal_gene_hits) + ' out of ' + str(len(causal_genes)))
print('score: ' + str(float(causal_gene_hits) / len(causal_genes)))
