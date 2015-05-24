from parsing_utils import *
import propagation
import statistics
import numpy
from scipy import stats
import sys
import operator
import os

# all patients 
patients = ["TCGA-AB-2803-03", "TCGA-AB-2805-03", "TCGA-AB-2806-03", "TCGA-AB-2807-03", "TCGA-AB-2810-03", "TCGA-AB-2811-03", "TCGA-AB-2812-03", "TCGA-AB-2813-03", "TCGA-AB-2814-03", "TCGA-AB-2815-03", "TCGA-AB-2816-03", "TCGA-AB-2817-03", "TCGA-AB-2818-03", "TCGA-AB-2819-03", "TCGA-AB-2820-03", "TCGA-AB-2821-03", "TCGA-AB-2823-03", "TCGA-AB-2824-03", "TCGA-AB-2825-03", "TCGA-AB-2826-03", "TCGA-AB-2830-03", "TCGA-AB-2832-03", "TCGA-AB-2833-03", "TCGA-AB-2834-03", "TCGA-AB-2835-03", "TCGA-AB-2836-03", "TCGA-AB-2837-03", "TCGA-AB-2838-03", "TCGA-AB-2839-03", "TCGA-AB-2840-03", "TCGA-AB-2841-03", "TCGA-AB-2842-03", "TCGA-AB-2843-03", "TCGA-AB-2844-03", "TCGA-AB-2845-03", "TCGA-AB-2846-03", "TCGA-AB-2847-03", "TCGA-AB-2848-03", "TCGA-AB-2849-03", "TCGA-AB-2851-03", "TCGA-AB-2854-03", "TCGA-AB-2855-03", "TCGA-AB-2856-03", "TCGA-AB-2857-03", "TCGA-AB-2859-03", "TCGA-AB-2860-03", "TCGA-AB-2861-03", "TCGA-AB-2862-03", "TCGA-AB-2863-03", "TCGA-AB-2865-03", "TCGA-AB-2866-03", "TCGA-AB-2867-03", "TCGA-AB-2868-03", "TCGA-AB-2869-03", "TCGA-AB-2870-03", "TCGA-AB-2871-03", "TCGA-AB-2872-03", "TCGA-AB-2873-03", "TCGA-AB-2874-03", "TCGA-AB-2875-03", "TCGA-AB-2877-03", "TCGA-AB-2879-03", "TCGA-AB-2880-03", "TCGA-AB-2881-03", "TCGA-AB-2884-03", "TCGA-AB-2885-03", "TCGA-AB-2886-03", "TCGA-AB-2887-03", "TCGA-AB-2888-03", "TCGA-AB-2889-03", "TCGA-AB-2890-03", "TCGA-AB-2891-03", "TCGA-AB-2895-03", "TCGA-AB-2896-03", "TCGA-AB-2897-03", "TCGA-AB-2898-03", "TCGA-AB-2899-03", "TCGA-AB-2900-03", "TCGA-AB-2901-03", "TCGA-AB-2903-03", "TCGA-AB-2904-03", "TCGA-AB-2908-03", "TCGA-AB-2910-03", "TCGA-AB-2911-03", "TCGA-AB-2912-03", "TCGA-AB-2913-03", "TCGA-AB-2914-03", "TCGA-AB-2915-03", "TCGA-AB-2916-03", "TCGA-AB-2917-03", "TCGA-AB-2918-03", "TCGA-AB-2919-03", "TCGA-AB-2920-03", "TCGA-AB-2924-03", "TCGA-AB-2925-03", "TCGA-AB-2928-03", "TCGA-AB-2929-03", "TCGA-AB-2930-03", "TCGA-AB-2931-03", "TCGA-AB-2932-03", "TCGA-AB-2933-03", "TCGA-AB-2934-03", "TCGA-AB-2936-03", "TCGA-AB-2937-03", "TCGA-AB-2938-03", "TCGA-AB-2939-03", "TCGA-AB-2940-03", "TCGA-AB-2941-03", "TCGA-AB-2942-03", "TCGA-AB-2943-03", "TCGA-AB-2948-03", "TCGA-AB-2949-03", "TCGA-AB-2950-03", "TCGA-AB-2952-03", "TCGA-AB-2954-03", "TCGA-AB-2955-03", "TCGA-AB-2956-03", "TCGA-AB-2959-03", "TCGA-AB-2963-03", "TCGA-AB-2964-03", "TCGA-AB-2965-03", "TCGA-AB-2966-03", "TCGA-AB-2967-03", "TCGA-AB-2969-03", "TCGA-AB-2970-03", "TCGA-AB-2971-03", "TCGA-AB-2972-03", "TCGA-AB-2973-03", "TCGA-AB-2975-03", "TCGA-AB-2976-03", "TCGA-AB-2977-03", "TCGA-AB-2978-03", "TCGA-AB-2979-03", "TCGA-AB-2980-03", "TCGA-AB-2981-03", "TCGA-AB-2982-03", "TCGA-AB-2983-03", "TCGA-AB-2984-03", "TCGA-AB-2985-03", "TCGA-AB-2986-03", "TCGA-AB-2987-03", "TCGA-AB-2988-03", "TCGA-AB-2990-03", "TCGA-AB-2991-03", "TCGA-AB-2992-03", "TCGA-AB-2993-03", "TCGA-AB-2994-03", "TCGA-AB-2995-03", "TCGA-AB-2996-03", "TCGA-AB-2998-03", "TCGA-AB-2999-03", "TCGA-AB-3000-03", "TCGA-AB-3001-03", "TCGA-AB-3002-03", "TCGA-AB-3005-03", "TCGA-AB-3006-03", "TCGA-AB-3007-03", "TCGA-AB-3008-03", "TCGA-AB-3009-03", "TCGA-AB-3011-03", "TCGA-AB-3012-03", "TCGA-AB-2921-03", "TCGA-AB-2882-03", "TCGA-AB-2858-03", "TCGA-AB-2927-03", "TCGA-AB-2909-03", "TCGA-AB-2808-03", "TCGA-AB-2946-03", "TCGA-AB-2944-03", "TCGA-AB-2853-03", "TCGA-AB-2828-03", "TCGA-AB-2822-03", "TCGA-AB-2935-03"]


g = DTWGraph()
print('Loading PPI...')
loadPPIData("data/PPI_HIPPIE.txt", g)
print('Normalizing weights...')
g.normalizeWeights()
gene_num = len(g.nodes)

causal_genes_raw = loadCausalGenes("data/AML_cosmic_genes.txt", g)
causal_genes = causal_genes_raw[9:]
print(causal_genes)
left_out = causal_genes_raw[:9]
print(left_out)
f_left_out = open('output_cv_1/left_out.txt', 'w')
for i in range(len(left_out)):
    f_left_out.write(left_out[i] + '\n')
f_left_out.close()

for i in range(len(patients)):
    patient = patients[i]
    if os.path.isfile('output_cv_1/' + patient + '_diff.txt'):
        print("Skipping patient", patient, ' #' + str(i))
        continue
    
    print("Working on patient", patient, ' #' + str(i))
    g.initGraph()
    """
    print('Loading expression...')
    ge_loaded = loadExpressionData("data/AML_Expression.txt", g, patient)
    if (ge_loaded == 0):
        print('no expression data')
        # exit(0)
    """
    print('Loading mutations...')
    mt_loaded = loadMutationData("data/AML_Mutations.txt", g, patient)
    if (mt_loaded == 0):
        print('no mutation data')        
        continue
    """
    print('Propagating from expression...')
    GEscores, GE_iterations = propagation.propagate(g, 'GE')
    GEranks = gene_num - stats.rankdata(GEscores)
    """ 
    print('Propagating from mutation...')
    MTscores, MT_iterations = propagation.propagate(g, 'MT')
    MTranks = gene_num - stats.rankdata(MTscores)
    '''
    g_ranks = [0 for i in range(gene_num)]
    print(len(GEranks))
    for gene in g.nodes:
        gene_index = g.gene2index[gene]
        g_ranks[gene_index] = max(GEranks[gene_index], MTranks[gene_index])
    '''

    k = round(gene_num / 10)
    diff_per_gene = {}
    diff_per_gene_b2h = {}
    
    count = 0
    mutation_num = 0
    for gene in g.nodes:
        if g.nodes[gene].mutation_type != None:
            mutation_num += 1
    print('mutations:', mutation_num)
    healthy_dists = statistics.generateHealthyDist(g, causal_genes, mutation_num)

    for gene in g.nodes:
        if gene in causal_genes:
            continue
        if g.nodes[gene].mutation_type == None and MTranks[g.gene2index[gene]] < k:
            print('Knocking out', gene, '(', count, 'out of', k, ')')
            count += 1
        
            sub_MTscores, sub_MT_iterations = propagation.propagate(g, 'MT', KNOCKOUT_IDX = g.gene2index[gene])
            sub_MTranks = gene_num - stats.rankdata(sub_MTscores)
            
            diff_b2h = statistics.getB2HValue(g, g, MTranks, sub_MTranks, healthy_dists)
            diff = statistics.getDiffValue(g, g, MTranks, sub_MTranks, causal_genes)
            diff_per_gene_b2h[gene] = diff_b2h
            diff_per_gene[gene] = diff
            print('diff:', diff, 'diff_b2h:', diff_b2h)


    sorted_diffs = sorted(diff_per_gene.items(), key=operator.itemgetter(1), reverse=True)
    sorted_diffs_b2h = sorted(diff_per_gene_b2h.items(), key=operator.itemgetter(1), reverse=True)
    f = open('output_cv_1/' + patient + '_diff.txt', 'w')
    f_b2h = open('output_cv_1/' + patient + '_b2h.txt', 'w')
    for gene, diff in sorted_diffs:
        f.write(gene + '\t' + str(diff) + '\n')
    for gene, diff in sorted_diffs_b2h:
        f_b2h.write(gene + '\t' + str(diff) + '\n')
    f.close()
    f_b2h.close()
