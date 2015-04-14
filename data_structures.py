import numpy as np
import copy

class DTWNode:
    'Node in the DTW graph'

    def __init__(self, mutation_type, expression_level):
        self.mutation_type = mutation_type
        self.expression_level = expression_level
        self.neighbors = {}

    def addNeighbor(self, gene_id, weight):
        self.neighbors[gene_id] = weight

    def removeNeighbor(self, gene_id):
        del self.neighbors[gene_id]

class DTWGraph:

    def __init__(self):
        self.nodes = {}
        self.weights = {}

    def addNode(self, gene_id, mutation_type=None, expression_level=None):
        node = DTWNode(mutation_type, expression_level)
        self.nodes[gene_id] = node

    def addEdge(self, gene1, gene2, weight):
        node1 = self.nodes[gene1]
        node1.addNeighbor(gene2, weight)
        node2 = self.nodes[gene2]
        node2.addNeighbor(gene1, weight)
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
            self.nodes[key[0]].neighbors[key[1]] = new_weights[key]
            self.nodes[key[1]].neighbors[key[0]] = new_weights[key]


    def initGraph(self):
        for node in self.nodes:
            self.nodes[node].mutation_type = None
            self.nodes[node].expression_level = None

    def createSubGraph(self, pr_values, threshold):
        print('original graph contains ' + len(self.nodes) + ' nodes and ' + len(self.weights) + ' edges')
        newGraph = copy.deepcopy(self)
        for node in newGraph.nodes:
            if pr_values[node] < threshold:
                newGraph.removeNode(node)
        print('sub graph contains ' + len(newGraph.nodes) + ' nodes and ' + len(newGraph.weights) + ' edges')
        return newGraph