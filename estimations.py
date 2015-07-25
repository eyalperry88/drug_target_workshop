from parsing_utils import *
import propagation
import numpy as np
import operator
from scipy import stats
from statistics import *

"""
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

# top 50
# patients = ["TCGA-AB-2803-03", "TCGA-AB-2805-03", "TCGA-AB-2806-03", "TCGA-AB-2807-03", "TCGA-AB-2810-03", "TCGA-AB-2811-03", "TCGA-AB-2812-03", "TCGA-AB-2813-03", "TCGA-AB-2814-03", "TCGA-AB-2815-03", "TCGA-AB-2816-03", "TCGA-AB-2817-03", "TCGA-AB-2818-03", "TCGA-AB-2819-03", "TCGA-AB-2820-03", "TCGA-AB-2821-03", "TCGA-AB-2823-03", "TCGA-AB-2824-03", "TCGA-AB-2825-03", "TCGA-AB-2826-03", "TCGA-AB-2830-03", "TCGA-AB-2832-03", "TCGA-AB-2833-03", "TCGA-AB-2834-03", "TCGA-AB-2835-03", "TCGA-AB-2836-03", "TCGA-AB-2837-03", "TCGA-AB-2838-03", "TCGA-AB-2839-03", "TCGA-AB-2840-03", "TCGA-AB-2841-03", "TCGA-AB-2842-03", "TCGA-AB-2843-03", "TCGA-AB-2844-03", "TCGA-AB-2845-03", "TCGA-AB-2846-03", "TCGA-AB-2847-03", "TCGA-AB-2848-03", "TCGA-AB-2849-03", "TCGA-AB-2851-03", "TCGA-AB-2854-03", "TCGA-AB-2855-03", "TCGA-AB-2856-03", "TCGA-AB-2857-03", "TCGA-AB-2859-03", "TCGA-AB-2860-03", "TCGA-AB-2861-03", "TCGA-AB-2862-03", "TCGA-AB-2863-03", "TCGA-AB-2865-03"]
# top 100
# patients = ["TCGA-AB-2803-03", "TCGA-AB-2805-03", "TCGA-AB-2806-03", "TCGA-AB-2807-03", "TCGA-AB-2810-03", "TCGA-AB-2811-03", "TCGA-AB-2812-03", "TCGA-AB-2813-03", "TCGA-AB-2814-03", "TCGA-AB-2815-03", "TCGA-AB-2816-03", "TCGA-AB-2817-03", "TCGA-AB-2818-03", "TCGA-AB-2819-03", "TCGA-AB-2820-03", "TCGA-AB-2821-03", "TCGA-AB-2823-03", "TCGA-AB-2824-03", "TCGA-AB-2825-03", "TCGA-AB-2826-03", "TCGA-AB-2830-03", "TCGA-AB-2832-03", "TCGA-AB-2833-03", "TCGA-AB-2834-03", "TCGA-AB-2835-03", "TCGA-AB-2836-03", "TCGA-AB-2837-03", "TCGA-AB-2838-03", "TCGA-AB-2839-03", "TCGA-AB-2840-03", "TCGA-AB-2841-03", "TCGA-AB-2842-03", "TCGA-AB-2843-03", "TCGA-AB-2844-03", "TCGA-AB-2845-03", "TCGA-AB-2846-03", "TCGA-AB-2847-03", "TCGA-AB-2848-03", "TCGA-AB-2849-03", "TCGA-AB-2851-03", "TCGA-AB-2854-03", "TCGA-AB-2855-03", "TCGA-AB-2856-03", "TCGA-AB-2857-03", "TCGA-AB-2859-03", "TCGA-AB-2860-03", "TCGA-AB-2861-03", "TCGA-AB-2862-03", "TCGA-AB-2863-03", "TCGA-AB-2865-03", "TCGA-AB-2866-03", "TCGA-AB-2867-03", "TCGA-AB-2868-03", "TCGA-AB-2869-03", "TCGA-AB-2870-03", "TCGA-AB-2871-03", "TCGA-AB-2872-03", "TCGA-AB-2873-03", "TCGA-AB-2874-03", "TCGA-AB-2875-03", "TCGA-AB-2877-03", "TCGA-AB-2879-03", "TCGA-AB-2880-03", "TCGA-AB-2881-03", "TCGA-AB-2884-03", "TCGA-AB-2885-03", "TCGA-AB-2886-03", "TCGA-AB-2887-03", "TCGA-AB-2888-03", "TCGA-AB-2889-03", "TCGA-AB-2890-03", "TCGA-AB-2891-03", "TCGA-AB-2895-03", "TCGA-AB-2896-03", "TCGA-AB-2897-03", "TCGA-AB-2898-03", "TCGA-AB-2899-03", "TCGA-AB-2900-03", "TCGA-AB-2901-03", "TCGA-AB-2903-03", "TCGA-AB-2904-03", "TCGA-AB-2908-03", "TCGA-AB-2910-03", "TCGA-AB-2911-03", "TCGA-AB-2912-03", "TCGA-AB-2913-03", "TCGA-AB-2914-03", "TCGA-AB-2915-03", "TCGA-AB-2916-03", "TCGA-AB-2917-03", "TCGA-AB-2918-03", "TCGA-AB-2919-03", "TCGA-AB-2920-03", "TCGA-AB-2924-03", "TCGA-AB-2925-03", "TCGA-AB-2928-03", "TCGA-AB-2929-03", "TCGA-AB-2930-03", "TCGA-AB-2931-03", "TCGA-AB-2932-03"]
# all patients 
patients = ["TCGA-AB-2803-03", "TCGA-AB-2805-03", "TCGA-AB-2806-03", "TCGA-AB-2807-03", "TCGA-AB-2810-03", "TCGA-AB-2811-03", "TCGA-AB-2812-03", "TCGA-AB-2813-03", "TCGA-AB-2814-03", "TCGA-AB-2815-03", "TCGA-AB-2816-03", "TCGA-AB-2817-03", "TCGA-AB-2818-03", "TCGA-AB-2819-03", "TCGA-AB-2820-03", "TCGA-AB-2821-03", "TCGA-AB-2823-03", "TCGA-AB-2824-03", "TCGA-AB-2825-03", "TCGA-AB-2826-03", "TCGA-AB-2830-03", "TCGA-AB-2832-03", "TCGA-AB-2833-03", "TCGA-AB-2834-03", "TCGA-AB-2835-03", "TCGA-AB-2836-03", "TCGA-AB-2837-03", "TCGA-AB-2838-03", "TCGA-AB-2839-03", "TCGA-AB-2840-03", "TCGA-AB-2841-03", "TCGA-AB-2842-03", "TCGA-AB-2843-03", "TCGA-AB-2844-03", "TCGA-AB-2845-03", "TCGA-AB-2846-03", "TCGA-AB-2847-03", "TCGA-AB-2848-03", "TCGA-AB-2849-03", "TCGA-AB-2851-03", "TCGA-AB-2854-03", "TCGA-AB-2855-03", "TCGA-AB-2856-03", "TCGA-AB-2857-03", "TCGA-AB-2859-03", "TCGA-AB-2860-03", "TCGA-AB-2861-03", "TCGA-AB-2862-03", "TCGA-AB-2863-03", "TCGA-AB-2865-03", "TCGA-AB-2866-03", "TCGA-AB-2867-03", "TCGA-AB-2868-03", "TCGA-AB-2869-03", "TCGA-AB-2870-03", "TCGA-AB-2871-03", "TCGA-AB-2872-03", "TCGA-AB-2873-03", "TCGA-AB-2874-03", "TCGA-AB-2875-03", "TCGA-AB-2877-03", "TCGA-AB-2879-03", "TCGA-AB-2880-03", "TCGA-AB-2881-03", "TCGA-AB-2884-03", "TCGA-AB-2885-03", "TCGA-AB-2886-03", "TCGA-AB-2887-03", "TCGA-AB-2888-03", "TCGA-AB-2889-03", "TCGA-AB-2890-03", "TCGA-AB-2891-03", "TCGA-AB-2895-03", "TCGA-AB-2896-03", "TCGA-AB-2897-03", "TCGA-AB-2898-03", "TCGA-AB-2899-03", "TCGA-AB-2900-03", "TCGA-AB-2901-03", "TCGA-AB-2903-03", "TCGA-AB-2904-03", "TCGA-AB-2908-03", "TCGA-AB-2910-03", "TCGA-AB-2911-03", "TCGA-AB-2912-03", "TCGA-AB-2913-03", "TCGA-AB-2914-03", "TCGA-AB-2915-03", "TCGA-AB-2916-03", "TCGA-AB-2917-03", "TCGA-AB-2918-03", "TCGA-AB-2919-03", "TCGA-AB-2920-03", "TCGA-AB-2924-03", "TCGA-AB-2925-03", "TCGA-AB-2928-03", "TCGA-AB-2929-03", "TCGA-AB-2930-03", "TCGA-AB-2931-03", "TCGA-AB-2932-03", "TCGA-AB-2933-03", "TCGA-AB-2934-03", "TCGA-AB-2936-03", "TCGA-AB-2937-03", "TCGA-AB-2938-03", "TCGA-AB-2939-03", "TCGA-AB-2940-03", "TCGA-AB-2941-03", "TCGA-AB-2942-03", "TCGA-AB-2943-03", "TCGA-AB-2948-03", "TCGA-AB-2949-03", "TCGA-AB-2950-03", "TCGA-AB-2952-03", "TCGA-AB-2954-03", "TCGA-AB-2955-03", "TCGA-AB-2956-03", "TCGA-AB-2959-03", "TCGA-AB-2963-03", "TCGA-AB-2964-03", "TCGA-AB-2965-03", "TCGA-AB-2966-03", "TCGA-AB-2967-03", "TCGA-AB-2969-03", "TCGA-AB-2970-03", "TCGA-AB-2971-03", "TCGA-AB-2972-03", "TCGA-AB-2973-03", "TCGA-AB-2975-03", "TCGA-AB-2976-03", "TCGA-AB-2977-03", "TCGA-AB-2978-03", "TCGA-AB-2979-03", "TCGA-AB-2980-03", "TCGA-AB-2981-03", "TCGA-AB-2982-03", "TCGA-AB-2983-03", "TCGA-AB-2984-03", "TCGA-AB-2985-03", "TCGA-AB-2986-03", "TCGA-AB-2987-03", "TCGA-AB-2988-03", "TCGA-AB-2990-03", "TCGA-AB-2991-03", "TCGA-AB-2992-03", "TCGA-AB-2993-03", "TCGA-AB-2994-03", "TCGA-AB-2995-03", "TCGA-AB-2996-03", "TCGA-AB-2998-03", "TCGA-AB-2999-03", "TCGA-AB-3000-03", "TCGA-AB-3001-03", "TCGA-AB-3002-03", "TCGA-AB-3005-03", "TCGA-AB-3006-03", "TCGA-AB-3007-03", "TCGA-AB-3008-03", "TCGA-AB-3009-03", "TCGA-AB-3011-03", "TCGA-AB-3012-03", "TCGA-AB-2921-03", "TCGA-AB-2882-03", "TCGA-AB-2858-03", "TCGA-AB-2927-03", "TCGA-AB-2909-03", "TCGA-AB-2808-03", "TCGA-AB-2946-03", "TCGA-AB-2944-03", "TCGA-AB-2853-03", "TCGA-AB-2828-03", "TCGA-AB-2822-03", "TCGA-AB-2935-03"]


g = DTWGraph()
print('Loading PPI...')
loadPPIData("data/PPI_HIPPIE.txt", g)
print('Normalizing weights...')
g.normalizeWeights()
gene_num = len(g.nodes)
results_avg = {}
results_max = {}
results_mut = {}
results_de = {}
run_stats = {'expression_iterations' : [], 'mutation_iterations' : []}

count_patient = 1
for patient in patients:
    print("Working on patient", patient, '#' + str(count_patient))
    count_patient += 1
    if patient in results_mut:
        print("patient already exists - skipping")
        continue;
        
    g.initGraph()
    """
    print('Loading expression...')
    ge_loaded = loadExpressionData("data/AML_Expression.txt", g, patient)
    if (ge_loaded == 0):
        print('no expression data')
        continue
    """
    print('Loading mutations...')
    mt_loaded = loadMutationData("data/AML_Mutations.txt", g, patient)
    if (mt_loaded == 0):
        print('no mutation data')        
        continue
    """
    print('Propagating from expression...')
    GEscores, GE_iterations = propagation.propagate(g, 'GE')
    run_stats['expression_iterations'].append(GE_iterations)
    GEranks = gene_num - stats.rankdata(GEscores)
    """
    print('Propagating from mutation...')
    MTscores, MT_iterations = propagation.propagate(g, 'MT')
    run_stats['mutation_iterations'].append(MT_iterations)
    MTranks = gene_num - stats.rankdata(MTscores)
    
    results_avg[patient] = {}
    results_max[patient] = {}
    results_mut[patient] = {}
    results_de[patient] = {}
    for gene in g.nodes:
        results_mut[patient][gene] = MTranks[g.gene2index[gene]]

actual_patients = list(results_mut.keys())    



k = round(gene_num / 10) # using top 10 percent
expected = {}
observed = {}
"""
for gene in g.nodes:
    expected[gene] = float(k * len(actual_patients)) / gene_num
    observed[gene] = 0
for i in range(0, len(actual_patients)):
    sorted_scores = sorted(results_mut[actual_patients[i]].items(), key=operator.itemgetter(1), reverse=False)
    for gene, score in sorted_scores[0:k]:
        observed[gene] += 1
threshold = len(actual_patients) / 2
"""
sum_of_ranks = np.zeros(gene_num)
count_below = 0
count_above = 0
for patient in actual_patients:
    for gene in results_mut[patient]:
        sum_of_ranks[g.gene2index[gene]] += results_mut[patient][gene]
        if gene == 'FLT3':
            if results_mut[patient][gene] < k:
                count_below += 1
            else:
                count_above += 1
ranked_sum_of_ranks = sum_of_ranks.argsort().argsort()
for gene in g.nodes:
    if ranked_sum_of_ranks[g.gene2index[gene]] <= k:
        observed[gene] = 1
threshold = 0


igenes = ["NPM","FLT3","DNM3A","RUNX1","IDHP","IDHC","TET2","RARA","PML","RASN",
"P53","CEBPA","WT1","PEBB","MYH11","KMT2A","PTN11","KIT","MTG8","U2AF1",
"RASK","SMC1A","SMC3"]

for igene in igenes:
    if igene in g.nodes:
        if igene in observed and observed[igene] == 1:
            print(igene,'is in the top')
    else:
        print('bad name!', igene)
"""
causal_genes = loadCausalGenes("data/AML_KEGG_genes.txt", g)
causal_gene_hits = 0
above_threshold_genes = 0
for gene in observed:
    if observed[gene] > threshold:
        above_threshold_genes += 1
        if gene in causal_genes:
            causal_gene_hits += 1
print('hit ' + str(causal_gene_hits) + ' out of ' + str(len(causal_genes)) + '(above threshold: ' + str(above_threshold_genes) + ')')
print('score: ' + str(float(causal_gene_hits) / len(causal_genes)))
p = stats.hypergeom.sf(causal_gene_hits, gene_num, len(causal_genes), above_threshold_genes)
print('p-value: ' + str(p))
print('-log(p-value): ' + str(-np.log10(p)))

causal_genes = loadCausalGenes("data/AML_cosmic_genes.txt", g)
causal_gene_hits = 0
above_threshold_genes = 0
for gene in observed:
    if observed[gene] > threshold:
        above_threshold_genes += 1
        if gene in causal_genes:
            causal_gene_hits += 1
print('hit ' + str(causal_gene_hits) + ' out of ' + str(len(causal_genes)) + '(above threshold: ' + str(above_threshold_genes) + ')')
print('score: ' + str(float(causal_gene_hits) / len(causal_genes)))
p = stats.hypergeom.sf(causal_gene_hits, gene_num, len(causal_genes), above_threshold_genes)
print('p-value: ' + str(p))
print('-log(p-value): ' + str(-np.log10(p)))

causal_genes = loadCausalGenes("data/cancer_cosmic_genes.txt", g)
causal_gene_hits = 0
above_threshold_genes = 0
for gene in observed:
    if observed[gene] > threshold:
        above_threshold_genes += 1
        if gene in causal_genes:
            causal_gene_hits += 1
print('hit ' + str(causal_gene_hits) + ' out of ' + str(len(causal_genes)) + '(above threshold: ' + str(above_threshold_genes) + ')')
print('score: ' + str(float(causal_gene_hits) / len(causal_genes)))
p = stats.hypergeom.sf(causal_gene_hits, gene_num, len(causal_genes), above_threshold_genes)
print('p-value: ' + str(p))
print('-log(p-value): ' + str(-np.log10(p)))

causal_genes = loadCausalGenes("data/AML_drug_targets.txt", g)
causal_gene_hits = 0
above_threshold_genes = 0
for gene in observed:
    if observed[gene] > threshold:
        above_threshold_genes += 1
        if gene in causal_genes:
            causal_gene_hits += 1
print('hit ' + str(causal_gene_hits) + ' out of ' + str(len(causal_genes)) + '(above threshold: ' + str(above_threshold_genes) + ')')
print('score: ' + str(float(causal_gene_hits) / len(causal_genes)))
p = stats.hypergeom.sf(causal_gene_hits, gene_num, len(causal_genes), above_threshold_genes)
print('p-value: ' + str(p))
print('-log(p-value): ' + str(-np.log10(p)))


total_ranks = np.zeros(gene_num)
for patient in actual_patients:
    for gene in results_max[patient]:
        total_ranks[g.gene2index[gene]] += results_avg[patient][gene]
labels = getLabelsVector(g, causal_genes, total_ranks)
p, mHG_idx = mHG(labels)
print('mHG for average: ', p, 'index: ', mHG_idx)

total_ranks = np.zeros(gene_num)
for patient in actual_patients:
    for gene in results_max[patient]:
        total_ranks[g.gene2index[gene]] += results_max[patient][gene]
labels = getLabelsVector(g, causal_genes, total_ranks)
p, mHG_idx = mHG(labels)
print('mHG for max: ', p, 'index: ', mHG_idx)

total_ranks = np.zeros(gene_num)
for patient in actual_patients:
    for gene in results_max[patient]:
        total_ranks[g.gene2index[gene]] += results_mut[patient][gene]
labels = getLabelsVector(g, causal_genes, total_ranks)
p, mHG_idx = mHG(labels)
print('mHG for mutation_only: ', p, 'index: ', mHG_idx)

total_ranks = np.zeros(gene_num)
for patient in actual_patients:
    for gene in results_max[patient]:
        total_ranks[g.gene2index[gene]] += results_de[patient][gene]
labels = getLabelsVector(g, causal_genes, total_ranks)
p, mHG_idx = mHG(labels)
print('mHG for de only: ', p, 'index: ', mHG_idx)
"""