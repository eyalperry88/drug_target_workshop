from data_structures import *

'reading and parsing data helpers'

def loadExpressionData(filename, graph): #add option to do it by patient?
    file = open(filename, 'r')
    next(file) #first line has no data
    for line in file:
        gene_data = line.split('\t')
        node = graph.nodes[gene_data[0]]
        node.expression_level = gene_data[2]
        node.pr_value = float(gene_data[3].replace('\n', '')) #in the future - normalize z-score 
    file.close()

def loadMutationData(filename, graph): #same...
    file = open(filename, 'r')
    next(file) #first line has no data
    for line in file:
        gene_data = line.split('\t')
        #filter by status?
        graph.nodes[gene_data[0]].mutation_type = gene_data[4] 
    file.close()
