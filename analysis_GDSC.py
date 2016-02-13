""" Checks B2H scores for different drugs """

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
FLT_patients = ["TCGA-AB-2870-03","TCGA-AB-2998-03","TCGA-AB-2875-03","TCGA-AB-2840-03","TCGA-AB-2980-03","TCGA-AB-2906-03","TCGA-AB-2942-03","TCGA-AB-2994-03","TCGA-AB-3007-03","TCGA-AB-2957-03","TCGA-AB-2811-03","TCGA-AB-2812-03","TCGA-AB-2909-03","TCGA-AB-2963-03","TCGA-AB-2869-03","TCGA-AB-2965-03","TCGA-AB-2871-03","TCGA-AB-2816-03","TCGA-AB-2818-03","TCGA-AB-2988-03","TCGA-AB-2986-03","TCGA-AB-2913-03","TCGA-AB-2918-03","TCGA-AB-2825-03","TCGA-AB-2877-03","TCGA-AB-2879-03","TCGA-AB-2880-03","TCGA-AB-2974-03","TCGA-AB-2830-03","TCGA-AB-2921-03","TCGA-AB-2925-03","TCGA-AB-2926-03","TCGA-AB-2970-03","TCGA-AB-2922-03","TCGA-AB-2836-03","TCGA-AB-2969-03","TCGA-AB-2924-03","TCGA-AB-2930-03","TCGA-AB-2931-03","TCGA-AB-2934-03","TCGA-AB-2851-03","TCGA-AB-2976-03","TCGA-AB-2853-03","TCGA-AB-3006-03","TCGA-AB-2981-03","TCGA-AB-2947-03","TCGA-AB-2993-03","TCGA-AB-2900-03","TCGA-AB-2895-03","TCGA-AB-2910-03","TCGA-AB-2814-03","TCGA-AB-2915-03","TCGA-AB-2834-03","TCGA-AB-2844-03","TCGA-AB-2928-03","TCGA-AB-2945-03"]

g = DTWGraph()
print('Loading PPI...')
loadPPIData("data/PPI_HIPPIE.txt", g)
print('Normalizing weights...')
g.normalizeWeights()
gene_num = len(g.nodes)

#causal_genes = loadCausalGenes("data/AML_cosmic_genes.txt", g)

mut_f = open("data/AML_Mutations.txt")
mut_per_patient = {}
for line in mut_f:
    gene_data = line.split('\t')
    pat = gene_data[1]
    gen = gene_data[0]
    if pat in mut_per_patient:
        mut_per_patient[pat].append(gen)
    else:
        mut_per_patient[pat] = [gen]
mut_f.close()

drugs = ['AP-24534', 'AZD7762', 'Bosutinib', 'CEP-701', 
         'FH535', 'GW441756', 'Sorafenib', 'Sunitinib']
         
f_outputs = {}
targets = {}
for drug in drugs:
    f_outputs[drug] = open('GDSC/results_per_drug/' + drug + '.txt', 'w')
    
    f_tmp = open('GDSC/' + drug + '.txt', 'r')
    tmp_targets = []
    for line in f_tmp:
        tmp_targets.append(line.strip())
    f_tmp.close()
    
    targets[drug] = tmp_targets

for i in range(len(patients)):
    patient = patients[i]
    
    """
    if patient not in mut_per_patient or 'FLT3' not in mut_per_patient[patient]:
        f_skipped.write('Skipping ' + patient + ' (1)\n')+
        continue
    """
    if patient not in FLT_patients:
        continue
    
    print("Working on patient", patient, ' #' + str(i))
    g.initGraph()
    
    print('Loading mutations...')
    mt_loaded = loadMutationData("data/AML_Mutations.txt", g, patient)
    if (mt_loaded == 0):
        print('no mutation data') 
        # f_skipped.write('Skipping ' + patient + ' (2)\n')
        continue
    
    """
    if g.nodes['FLT3'].mutation_type == None:
        print('Skipping (no FLT mutation)')
        f_skipped.write('Skipping ' + patient + ' (3)\n')
        continue
    """
    
    print('Propagating from mutation...')
    MTscores, MT_iterations = propagation.propagate(g, 'MT')
    MTranks = gene_num - stats.rankdata(MTscores)
    
    try:
        causal_genes = loadCausalGenes("data/aliases/exp/" + patient + "_exp_aliases.txt", g)
    except:
        # f_skipped.write('Skipping ' + patient + ' (4)\n')
        continue
    for drug in drugs:
        for x in targets[drug]:
            if x in causal_genes:
                causal_genes.remove(x)
    
    diff_per_gene_b2h = {}
    
    count = 0
    mutation_num = 0
    for gene in g.nodes:
        if g.nodes[gene].mutation_type != None:
            mutation_num += 1
    print('mutations:', mutation_num)
    healthy_dists = statistics.generateHealthyDist(g, causal_genes, mutation_num)

    for drug in drugs:
        print('\n', drug, 'knockout')
        
        gene_indices = []
        for gene in targets[drug]:
            """
            if g.nodes[gene].mutation_type != None:
                print(gene, 'is mutated')
                break
            """
            gene_indices.append(g.gene2index[gene])
        if len(gene_indices) == len(targets[drug]):
            sub_MTscores, sub_MT_iterations = propagation.propagate(g, 'MT', KNOCKOUT_IDX = gene_indices)
            sub_MTranks = gene_num - stats.rankdata(sub_MTscores)
            
            diff_b2h = statistics.getB2HValue(g, g, MTranks, sub_MTranks, healthy_dists)
            
            f_outputs[drug].write(str(diff_b2h) + '\n')
            print('B2H:', diff_b2h)
        print()

for drug in drugs:
    f_outputs[drug].close()
