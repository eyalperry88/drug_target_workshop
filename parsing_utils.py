from data_structures import *

'reading and parsing data helpers'
    
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
    return count

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
    return count

def loadPPIData(filename, graph):
    file = open(filename, 'r')
    count_nodes = 0
    count_edges = 0
    for line in file:
        gene_data = line.strip().split('\t')
        if len(gene_data) != 3 or float(gene_data[2]) == 0:
            continue
        for gene_id in [gene_data[0].strip(), gene_data[1].strip()]:
            if gene_id not in graph.nodes:
                graph.addNode(gene_id)
                count_nodes += 1
        count_edges += 1
        graph.addEdge(gene_data[0].strip(), gene_data[1].strip(), float(gene_data[2]))
    print('Loaded ' + str(count_nodes) + ' nodes and ' + str(count_edges) + ' edges.')
    file.close()
