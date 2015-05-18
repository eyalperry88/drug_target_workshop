import numpy as np
import operator
from statistics import *
from parsing_utils import *

# all patients 
patients = ["TCGA-AB-2803-03", "TCGA-AB-2805-03", "TCGA-AB-2806-03", "TCGA-AB-2807-03", "TCGA-AB-2810-03", "TCGA-AB-2811-03", "TCGA-AB-2812-03", "TCGA-AB-2813-03", "TCGA-AB-2814-03", "TCGA-AB-2815-03", "TCGA-AB-2816-03", "TCGA-AB-2817-03", "TCGA-AB-2818-03", "TCGA-AB-2819-03", "TCGA-AB-2820-03", "TCGA-AB-2821-03", "TCGA-AB-2823-03", "TCGA-AB-2824-03", "TCGA-AB-2825-03", "TCGA-AB-2826-03", "TCGA-AB-2830-03", "TCGA-AB-2832-03", "TCGA-AB-2833-03", "TCGA-AB-2834-03", "TCGA-AB-2835-03", "TCGA-AB-2836-03", "TCGA-AB-2837-03", "TCGA-AB-2838-03", "TCGA-AB-2839-03", "TCGA-AB-2840-03", "TCGA-AB-2841-03", "TCGA-AB-2842-03", "TCGA-AB-2843-03", "TCGA-AB-2844-03", "TCGA-AB-2845-03", "TCGA-AB-2846-03", "TCGA-AB-2847-03", "TCGA-AB-2848-03", "TCGA-AB-2849-03", "TCGA-AB-2851-03", "TCGA-AB-2854-03", "TCGA-AB-2855-03", "TCGA-AB-2856-03", "TCGA-AB-2857-03", "TCGA-AB-2859-03", "TCGA-AB-2860-03", "TCGA-AB-2861-03", "TCGA-AB-2862-03", "TCGA-AB-2863-03", "TCGA-AB-2865-03", "TCGA-AB-2866-03", "TCGA-AB-2867-03", "TCGA-AB-2868-03", "TCGA-AB-2869-03", "TCGA-AB-2870-03", "TCGA-AB-2871-03", "TCGA-AB-2872-03", "TCGA-AB-2873-03", "TCGA-AB-2874-03", "TCGA-AB-2875-03", "TCGA-AB-2877-03", "TCGA-AB-2879-03", "TCGA-AB-2880-03", "TCGA-AB-2881-03", "TCGA-AB-2884-03", "TCGA-AB-2885-03", "TCGA-AB-2886-03", "TCGA-AB-2887-03", "TCGA-AB-2888-03", "TCGA-AB-2889-03", "TCGA-AB-2890-03", "TCGA-AB-2891-03", "TCGA-AB-2895-03", "TCGA-AB-2896-03", "TCGA-AB-2897-03", "TCGA-AB-2898-03", "TCGA-AB-2899-03", "TCGA-AB-2900-03", "TCGA-AB-2901-03", "TCGA-AB-2903-03", "TCGA-AB-2904-03", "TCGA-AB-2908-03", "TCGA-AB-2910-03", "TCGA-AB-2911-03", "TCGA-AB-2912-03", "TCGA-AB-2913-03", "TCGA-AB-2914-03", "TCGA-AB-2915-03", "TCGA-AB-2916-03", "TCGA-AB-2917-03", "TCGA-AB-2918-03", "TCGA-AB-2919-03", "TCGA-AB-2920-03", "TCGA-AB-2924-03", "TCGA-AB-2925-03", "TCGA-AB-2928-03", "TCGA-AB-2929-03", "TCGA-AB-2930-03", "TCGA-AB-2931-03", "TCGA-AB-2932-03", "TCGA-AB-2933-03", "TCGA-AB-2934-03", "TCGA-AB-2936-03", "TCGA-AB-2937-03", "TCGA-AB-2938-03", "TCGA-AB-2939-03", "TCGA-AB-2940-03", "TCGA-AB-2941-03", "TCGA-AB-2942-03", "TCGA-AB-2943-03", "TCGA-AB-2948-03", "TCGA-AB-2949-03", "TCGA-AB-2950-03", "TCGA-AB-2952-03", "TCGA-AB-2954-03", "TCGA-AB-2955-03", "TCGA-AB-2956-03", "TCGA-AB-2959-03", "TCGA-AB-2963-03", "TCGA-AB-2964-03", "TCGA-AB-2965-03", "TCGA-AB-2966-03", "TCGA-AB-2967-03", "TCGA-AB-2969-03", "TCGA-AB-2970-03", "TCGA-AB-2971-03", "TCGA-AB-2972-03", "TCGA-AB-2973-03", "TCGA-AB-2975-03", "TCGA-AB-2976-03", "TCGA-AB-2977-03", "TCGA-AB-2978-03", "TCGA-AB-2979-03", "TCGA-AB-2980-03", "TCGA-AB-2981-03", "TCGA-AB-2982-03", "TCGA-AB-2983-03", "TCGA-AB-2984-03", "TCGA-AB-2985-03", "TCGA-AB-2986-03", "TCGA-AB-2987-03", "TCGA-AB-2988-03", "TCGA-AB-2990-03", "TCGA-AB-2991-03", "TCGA-AB-2992-03", "TCGA-AB-2993-03", "TCGA-AB-2994-03", "TCGA-AB-2995-03", "TCGA-AB-2996-03", "TCGA-AB-2998-03", "TCGA-AB-2999-03", "TCGA-AB-3000-03", "TCGA-AB-3001-03", "TCGA-AB-3002-03", "TCGA-AB-3005-03", "TCGA-AB-3006-03", "TCGA-AB-3007-03", "TCGA-AB-3008-03", "TCGA-AB-3009-03", "TCGA-AB-3011-03", "TCGA-AB-3012-03", "TCGA-AB-2921-03", "TCGA-AB-2882-03", "TCGA-AB-2858-03", "TCGA-AB-2927-03", "TCGA-AB-2909-03", "TCGA-AB-2808-03", "TCGA-AB-2946-03", "TCGA-AB-2944-03", "TCGA-AB-2853-03", "TCGA-AB-2828-03", "TCGA-AB-2822-03", "TCGA-AB-2935-03"]

genes_total_diff = {}
genes_total_b2h = {}

for patient in patients:
    try:
        f = open('output/' + patient + '_diff.txt', 'r')
        for line in f:
            row_data = line.strip().split('\t')
            gene = row_data[0]
            if gene in genes_total_diff:
                genes_total_diff[gene] += np.sqrt(float(row_data[1]))
            else:
                genes_total_diff[gene] = np.sqrt(float(row_data[1]))
        f.close()
        
        f = open('output/' + patient + '_b2h.txt', 'r')
        for line in f:
            row_data = line.strip().split('\t')
            gene = row_data[0]
            if gene in genes_total_b2h:
                genes_total_b2h[gene] += float(row_data[1])
            else:
                genes_total_b2h[gene] = float(row_data[1])
        f.close()
    except FileNotFoundError:
        pass
    
g = DTWGraph()
loadPPIData("data/PPI_HIPPIE.txt", g)
drug_targets = loadCausalGenes("data/AML_drug_targets.txt", g)

sorted_diff = sorted(genes_total_diff.items(), key=operator.itemgetter(1), reverse=True)
f = open('output/total_diff.txt', 'w')
sorted_genes = []
for gene, diff in sorted_diff:
    f.write(gene + '\t' + str(diff) + '\n')
    sorted_genes.append(gene)
f.close()

print('stats for diff score')

labels = [1 if x in drug_targets else 0 for x in sorted_genes]
p, mHG_idx = mHG(labels)
print('dmHG', p, mHG_idx)

gene_num = len(sorted_genes)
k = round(gene_num / 10) # using top 10 percent
causal_gene_hits = 0
i = 0
causal_genes_in_all_genes = 0
for gene in sorted_genes:
    if gene in drug_targets:
        causal_genes_in_all_genes += 1
        if i < k:
            causal_gene_hits += 1
    i += 1
    
print('hit ' + str(causal_gene_hits) + ' out of ' + str(causal_genes_in_all_genes))
print('score: ' + str(float(causal_gene_hits) / len(drug_targets)))
p = stats.hypergeom.sf(causal_gene_hits, len(sorted_genes), causal_genes_in_all_genes, k) + stats.hypergeom.pmf(causal_gene_hits, len(sorted_genes), causal_genes_in_all_genes, k)
print('p-value: ' + str(p))
print('-log(p-value): ' + str(-np.log10(p)))

sorted_b2h = sorted(genes_total_b2h.items(), key=operator.itemgetter(1), reverse=True)
sorted_genes = []
f = open('output/total_b2h.txt', 'w')
for gene, diff in sorted_b2h:
    f.write(gene + '\t' + str(diff) + '\n')
    sorted_genes.append(gene)
f.close()

print('stats for b2h score')

labels = [1 if x in drug_targets else 0 for x in sorted_genes]
p, mHG_idx = mHG(labels)
print('dmHG', p, mHG_idx)

gene_num = len(sorted_genes)
k = round(gene_num / 10) # using top 10 percent
causal_gene_hits = 0
i = 0
causal_genes_in_all_genes = 0
for gene in sorted_genes:
    if gene in drug_targets:
        causal_genes_in_all_genes += 1
        if i < k:
            causal_gene_hits += 1
    i += 1
    
print('hit ' + str(causal_gene_hits) + ' out of ' + str(causal_genes_in_all_genes))
print('score: ' + str(float(causal_gene_hits) / len(drug_targets)))
p = stats.hypergeom.sf(causal_gene_hits, len(sorted_genes), causal_genes_in_all_genes, k) + stats.hypergeom.pmf(causal_gene_hits, len(sorted_genes), causal_genes_in_all_genes, k)
print('p-value: ' + str(p))
print('-log(p-value): ' + str(-np.log10(p)))

