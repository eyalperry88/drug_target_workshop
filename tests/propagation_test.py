from parsing_utils import *
from propagation import *

g = DTWGraph()
loadGenes("test_data/genes_test.txt", g)
loadExpressionData("test_data/exp_test.txt", g)
loadMutationData("test_data/mut_test.txt", g)
loadPPIData("test_data/simple_ppi_test.txt", g)

scores = propagate(g, 1, 0.5)
assert(scores["GENE1"] == 0.5)
assert(scores["GENE2"] == 0.5)
assert(scores["GENE3"] == 0)
scores = propagate(g, 2, 0.5)
assert(scores["GENE1"] == 0.5)
assert(scores["GENE2"] == 0.5)
assert(scores["GENE3"] == 0.5)
scores = propagate(g, 3, 0.5)
assert(scores["GENE1"] == 0.75)
assert(scores["GENE2"] == 0.75)
assert(scores["GENE3"] == 0.5)
print scores
