from parsing_utils import *
from propagation import *

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

for patient in patients:
    print("Working on patient ", patient)
    g = DTWGraph()
    print('Loading PPI...')
    loadPPIData("data/PPI_HIPPIE.txt", g)
    print('Normalizing weights...')
    g.normalizeWeights()
    print('Loading expression...')
    loadExpressionData("data/AML_Expression.txt", g, patient)
    print('Loading mutations...')
    loadMutationData("data/AML_Mutations.txt", g, patient)
    

    print('Propagating from expression...')
    GEscores = propagate(g, 'GE')
    print('Propagating from mutation...')
    MTscores = propagate(g, 'MT')

    import operator
    joint_scores = {}
    for gene in GEscores:
        joint_scores[gene] = (GEscores[gene] + MTscores[gene]) / 2
    sorted_scores = sorted(joint_scores.items(), key=operator.itemgetter(1), reverse=True)
    print("Top 10: ")
    print(sorted_scores[0:10])

