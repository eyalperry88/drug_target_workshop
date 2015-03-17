from parsing_utils import *

# this changes the working dir to the current one
# so we can load the local test data files
import os
abspath = os.path.abspath(__file__)
dname = os.path.dirname(abspath)
os.chdir(dname)

g = DTWGraph()
loadGenes("genes_test.txt", g)
assert("GENE1" in g.nodes)
assert("GENE2" in g.nodes)
loadExpressionData("exp_test.txt", g)
assert(g.nodes["GENE1"].expression_level == 'under')
assert(g.nodes["GENE1"].pr_value == -2.097)
assert(g.nodes["GENE2"].expression_level == 'over')
assert(g.nodes["GENE2"].pr_value == 2.058)
loadMutationData("mut_test.txt", g)
assert(g.nodes["GENE1"].mutation_type == 'Substitution - Missense')
assert(g.nodes["GENE2"].mutation_type == 'Substitution - Missense')
loadPPIData("ppi_test.txt", g)
assert(g.getEdgeWeight("GENE1", "GENE2") == 0.75)
