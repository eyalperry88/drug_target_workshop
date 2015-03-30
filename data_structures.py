import numpy as np

class DTWNode:
    'Node in the DTW graph'

    def __init__(self, mutation_type, expression_level):
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

    def addNode(self, gene_id, mutation_type=None, expression_level=None):
        node = DTWNode(mutation_type, expression_level)
        self.nodes[gene_id] = node

    def addEdge(self, gene1, gene2, weight):
        node1 = self.nodes[gene1]
        node1.addNeighbor(gene2)
        node2 = self.nodes[gene2]
        node2.addNeighbor(gene1)
        self.setEdgeWeight(gene1, gene2, weight)

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
        
    def setEdgeWeight(self, gene1, gene2, weight):
        self.weights[tuple(sorted((gene1, gene2)))] = weight

    def normalizeWeights(self):
        n = len(self.nodes)
        node_to_idx = {}
        i = 0
        for node in self.nodes:
            node_to_idx[node] = i
            i += 1
        sum_vec = np.zeros(n)
        for node in self.nodes:
            for neighbor in self.nodes[node].neighbors:
                sum_vec[node_to_idx[node]] += self.getEdgeWeight(node, neighbor)
        sum_vec = np.reciprocal(np.sqrt(sum_vec))
        new_weights = {}
        for node in self.nodes:
            for neighbor in self.nodes[node].neighbors:
                new_w = sum_vec[node_to_idx[node]] * self.getEdgeWeight(node, neighbor) * sum_vec[node_to_idx[neighbor]]
                new_weights[tuple(sorted((node, neighbor)))] = new_w
        for key in new_weights:
            self.weights[key] = new_weights[key]