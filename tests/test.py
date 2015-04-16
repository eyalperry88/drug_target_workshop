from data_structures import DTWGraph

def runTest():
    g = DTWGraph()
    g.addNode("GENE1")
    g.addNode("GENE2")
    g.addEdge("GENE1", "GENE2", 0.7)
    g.removeNode("GENE1")
    try:
        g.addEdge("GENE1", "GENE2", 0.7)
        print(':(')
    except:
        pass
    g.addNode("GENE3")
    g.addEdge("GENE3", "GENE2", 0.6)
    assert(g.getEdgeWeight("GENE2", "GENE3") == 0.6)
    g.normalizeWeights()
    g.removeEdge("GENE2", "GENE3")
    assert(len(g.nodes["GENE2"].neighbors) == 0)
    g.nodes["GENE2"].mutation_type = "some_type"
    g.nodes["GENE2"].expression_level = "some_exp"
    g.initGraph()
    assert(g.nodes["GENE2"].mutation_type == None)
    assert(g.nodes["GENE2"].expression_level == None)
    
    
