from parsing_utils import *
from propagation import *

def runTest():
    g = DTWGraph()
    loadPPIData("test_data/simple_ppi_test.txt", g)
    loadExpressionData("test_data/exp_test.txt", g, "TCGA-AB-2803-03")
    loadMutationData("test_data/mut_test.txt", g, "TCGA-AB-2963-03")
    
    # TODO: make this test work...
    scores, it = propagate(g, 'GE', ALPHA = 0.5, MAX_ITERATIONS = 1)
    assert(scores["GENE1"] == 0.5)
    assert(scores["GENE2"] == 0.5)
    assert(scores["GENE3"] == 0)
    scores, it = propagate(g, 'GE', ALPHA = 0.5, MAX_ITERATIONS = 2)
    assert(scores["GENE1"] == 0.5)
    assert(scores["GENE2"] == 0.5)
    assert(scores["GENE3"] == 0.5)
    scores, it = propagate(g, 'GE', ALPHA = 0.5, MAX_ITERATIONS = 3)
    assert(scores["GENE1"] == 0.75)
    assert(scores["GENE2"] == 0.75)
    assert(scores["GENE3"] == 0.5)
