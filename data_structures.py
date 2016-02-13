""" Graph data structures """
import numpy as np
import copy
from scipy import sparse
import time

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
        neighbors_copy = list(node.neighbors.keys())
        for neighbor_id in neighbors_copy:
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
        self.mapIndices()
        sum_vec = np.zeros(n)
        for node in self.nodes:
            for neighbor in self.nodes[node].neighbors:
                sum_vec[self.gene2index[node]] += self.getEdgeWeight(node, neighbor)
        sum_vec = np.reciprocal(np.sqrt(sum_vec))
        new_weights = {}
        for node in self.nodes:
            for neighbor in self.nodes[node].neighbors:
                new_w = sum_vec[self.gene2index[node]] * self.getEdgeWeight(node, neighbor) * sum_vec[self.gene2index[neighbor]]
                new_weights[tuple(sorted((node, neighbor)))] = new_w
        self.W = sparse.lil_matrix((n, n))
        for key, w in new_weights.items():
            self.weights[key] = w
            self.nodes[key[0]].neighbors[key[1]] = w
            self.nodes[key[1]].neighbors[key[0]] = w
            self.W[self.gene2index[key[0]], self.gene2index[key[1]]] = w
            self.W[self.gene2index[key[1]], self.gene2index[key[0]]] = w
        self.W = self.W.tocsr()
    
    def mapIndices(self):
        self.index2gene = []
        self.gene2index = {}
        idx = 0
        for gene in self.nodes:
            self.gene2index[gene] = idx
            self.index2gene.append(gene)
            idx += 1

    def initGraph(self):
        for node in self.nodes:
            self.nodes[node].mutation_type = None
            self.nodes[node].expression_level = None

    def createSubGraph(self, node):
        start = time.clock()
        newGraph = copy.deepcopy(self)
        """
        newGraph.mapIndices()
        for gene in newGraph.nodes:
            if newGraph.gene2index[gene] != 
        """ 
        n = len(newGraph.nodes)
        newGraph.removeNode(node)
        indices = [i for i in range(n)]
        indices.remove(newGraph.gene2index[node])
        newGraph.W = self.W[indices, :]
        newGraph.W = newGraph.W.tocsc()[:, indices]
        #newGraph.mapIndices()
        newGraph.index2gene.remove(node)
        newGraph.gene2index = {}
        for i in range(len(newGraph.index2gene)):
            newGraph.gene2index[newGraph.index2gene[i]] = i
        end = time.clock()
        print ('time2:', end-start)
        return newGraph
