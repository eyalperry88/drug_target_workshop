from parsing_utils import *
from propagation import *

g = DTWGraph()
print('Loading genes...')
loadGenes("data/genes.txt", g)
print('Loading expression...')
loadExpressionData("data/AML_Expression.txt", g, "TCGA-AB-2963-03")
print('Loading mutations...')
loadMutationData("data/AML_Mutations.txt", g, "TCGA-AB-2963-03")
print('Loading PPI...')
loadPPIData("data/PPI_HIPPIE.txt", g)
print('Normalizing weights...')
g.normalizeWeights()

print('Propagating from expression...')
GEscores = propagate(g, 'GE')
print('Propagating from mutation...')
MTscores = propagate(g, 'MT')
