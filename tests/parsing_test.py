from parsing_utils import *

def runTest():
    g = DTWGraph()
    loadPPIData("test_data/ppi_test.txt", g)
    assert("GENE1" in g.nodes)
    assert("GENE2" in g.nodes)
    assert("GENE3" in g.nodes)
    assert(g.getEdgeWeight("GENE1", "GENE2") == 0.75)
    loadExpressionData("test_data/exp_test.txt", g, "TCGA-AB-2803-03")
    assert(g.nodes["GENE1"].expression_level == 'under')
    assert(g.nodes["GENE2"].expression_level == 'over')
    assert(g.nodes["GENE3"].expression_level == None)
    loadMutationData("test_data/mut_test.txt", g, "TCGA-AB-2963-03")
    assert(g.nodes["GENE1"].mutation_type == 'Substitution - Missense')
    assert(g.nodes["GENE2"].mutation_type == 'Substitution - Missense')
    assert(g.nodes["GENE3"].mutation_type == None)
    
