
class DTWNode:
    'Node in the DTW graph'

    def __init__(self, pr_value, mutation_type, expression_level):
        self.pr_value = pr_value
        self.mutation_type = mutation_type
        self.expression_level = expression_level
        self.neighbors = []

    def addNeighbor(self, gene_id):
        self.neighbors.append(gene_id)

    def removeNeighbor(self, gene_id):
        self.neighbors.remove(gene_id)

class DTWGraph:

    def __init__(self):
        self.nodes = {}
        self.weights = {}

    def addNode(self, gene_id, pr_value=0, mutation_type=None, expression_level=None):
        node = DTWNode(pr_value, mutation_type, expression_level)
        self.nodes[gene_id] = node

    def addEdge(self, gene1, gene2, weight):
        node1 = self.nodes[gene1]
        node1.addNeighbor(gene2)
        node2 = self.nodes[gene2]
        node2.addNeighbor(gene1)
        self.weights[tuple(sorted((gene1, gene2)))] = weight

    def removeNode(self, gene_id):
        node = self.nodes[gene_id]
        for neighbor_id in node.neighbors:
            neighbor = self.nodes[neighbor_id]
            neighbor.removeNeighbor(gene_id)
            del self.weights[tuple(sorted((neighbor_id, gene_id)))]
        del self.nodes[gene_id]

    def removeEdge(self, gene1, gene2):
        node1 = self.nodes[gene1]
        node1.removeNeighbor(gene2)
        node2 = self.nodes[gene2]
        node2.removeNeighbor(gene1)
        del self.weights[tuple(sorted((gene1, gene2)))]

    def getEdgeWeight(self, gene1, gene2):
        key = tuple(sorted((gene1, gene2)))
        if key in self.weights:
            return self.weights.get(key)
        raise KeyError('EdgeWeight')

    def normalizeWeights(self):
        pass

    

    
                

    
    

    
