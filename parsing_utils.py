from data_structures import *

'reading and parsing data helpers'

def loadGenes(filename, graph):
    file = open(filename, 'r')
    count = 0
    for line in file:
        graph.addNode(line.strip())
        count += 1
    print('Loaded ' + str(count) + ' genes.')
    file.close()
    
def loadExpressionData(filename, graph, patient): #add option to do it by patient?
    file = open(filename, 'r')
    next(file) #first line has no data
    count = 0
    for line in file:
        gene_data = line.split('\t')
        if (gene_data[1] == patient) and (gene_data[0] in graph.nodes):
            node = graph.nodes[gene_data[0]]
            node.expression_level = gene_data[2]
            count += 1
    print('Loaded ' + str(count) + ' differentially expressed genes.')
    file.close()

def loadMutationData(filename, graph, patient): #same...
    file = open(filename, 'r')
    next(file) #first line has no data
    count = 0
    for line in file:
        gene_data = line.split('\t')
        if (gene_data[1] == patient) and (gene_data[0] in graph.nodes):
            #filter by status?
            graph.nodes[gene_data[0]].mutation_type = gene_data[4]
            count += 1
    print('Loaded ' + str(count) + ' mutated genes.')
    file.close()

def loadPPIData(filename, graph):
    file = open(filename, 'r')
    count = 0
    for line in file:
        gene_data = line.strip().split('\t')
        if len(gene_data) != 3 or float(gene_data[2]) == 0:
            continue
        count += 1
        graph.addEdge(gene_data[0].strip(), gene_data[1].strip(), float(gene_data[2]))
    print('Loaded ' + str(count) + ' edges.')
    file.close()
