from parsing_utils import *

g = DTWGraph()
g.addNode("GENE1")
g.addNode("GENE2")
g.addEdge("GENE1", "GENE2", 0.7)
loadExpressionData("exp_test.txt", g)
assert(g.nodes["GENE1"].expression_level == 'under')
assert(g.nodes["GENE1"].pr_value == -2.097)
assert(g.nodes["GENE2"].expression_level == 'over')
assert(g.nodes["GENE2"].pr_value == 2.058)
loadMutationData("mut_test.txt", g)
assert(g.nodes["GENE1"].mutation_type == 'Substitution - Missense')
assert(g.nodes["GENE2"].mutation_type == 'Substitution - Missense')
