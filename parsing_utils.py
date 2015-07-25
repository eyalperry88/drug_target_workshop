from data_structures import *
import urllib.request
import os

'reading and parsing data helpers'

def checkAliases(gene, all_genes):
    # check aliases - using UNIPROT API
    try:
        with urllib.request.urlopen('http://www.uniprot.org/uniprot/?query='+gene+'+AND+organism:9606&format=tab&columns=genes') as response:
            data = response.read()
            first_row = True
            for row in data.splitlines():
                if first_row:
                    first_row = False
                    continue
                genes = row.decode("utf-8").split()
                if gene in genes:
                    for other_gene in genes:
                        if other_gene in all_genes:
                            print('Found alias', other_gene, 'instead of', gene)
                            return other_gene
    except:
        print("error for: " + gene)
    return None

def loadExpressionData(filename, graph, patient): #add option to do it by patient?
    file = open(filename, 'r')
    next(file) #first line has no data
    count = 0
    for line in file:
        gene_data = line.split('\t')
        if gene_data[1] == patient or patient == 'ALL':
            if gene_data[0] in graph.nodes:
                node = graph.nodes[gene_data[0]]
                node.expression_level = gene_data[2]
                count += 1
            else:
                gene_true_name = checkAliases(gene_data[0], graph.nodes)
                if gene_true_name:
                    node = graph.nodes[gene_true_name]
                    node.expression_level = gene_data[2]
                    count += 1
                else:
                    print('Could not find suitable alias for', gene_data[0])
    print('Loaded ' + str(count) + ' differentially expressed genes.')
    file.close()
    return count
    
def loadExpressionData2(filename, graph):
    file = open(filename, 'r')
    count = 0
    for line in file:
        gene = line.strip()
        if gene in graph.nodes:
            graph.nodes[gene].expression_level = 'Y'
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
        if gene_data[1] == patient or patient == 'ALL':
            if gene_data[0] in graph.nodes:
                #filter by status?
                graph.nodes[gene_data[0]].mutation_type = gene_data[4]
                count += 1
            else:
                gene_true_name = checkAliases(gene_data[0], graph.nodes)
                if gene_true_name:
                    graph.nodes[gene_true_name].mutation_type = gene_data[4]
                    count += 1
                else:
                    pass
                    # print('Could not find suitable alias for', gene_data[0])
                    
    print('Loaded ' + str(count) + ' mutated genes.')
    file.close()
    return count
    
def getMutations(filename, graph, patient): #same...
    file = open(filename, 'r')
    next(file) #first line has no data
    count = 0
    mut = []
    for line in file:
        gene_data = line.split('\t')
        if gene_data[1] == patient or patient == 'ALL':
            if gene_data[0] in graph.nodes:
                #filter by status?
                mut.append(gene_data[0])
            else:
                gene_true_name = checkAliases(gene_data[0], graph.nodes)
                if gene_true_name:
                    mut.append(gene_true_name)
                else:
                    print('Could not find suitable alias for', gene_data[0])
    file.close()
    return mut

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
    

def loadCausalGenes(filename, graph):
    file = open(filename, 'r')
    count_genes = 0
    count_bad = 0
    causal_genes = []
    for line in file:
        gene = line.strip()
        if gene in graph.nodes:
            causal_genes.append(gene)
            count_genes += 1
        else:
            gene_true_name = checkAliases(gene, graph.nodes)
            if gene_true_name:
                causal_genes.append(gene_true_name)
                count_genes += 1
            else:
                print('Could not find suitable alias for', gene)
                count_bad += 1
    print('Loaded ' + str(count_genes) + ' genes (' + str(count_bad) + ' not found on graph)')
    file.close()
    return causal_genes

def getExpressionAliases(filename, graph):
    file = open(filename, 'r')
    next(file) #first line has no data
    count = 0
    skip = 0
    curr_patient = ''
    for line in file:
        gene_data = line.split('\t')
        last_patient = curr_patient
        curr_patient = gene_data[1]
        if curr_patient != last_patient:
            if os.path.isfile('data/aliases/exp/' + curr_patient + '_exp_aliases.txt'):
                print("Skipping patient" + curr_patient)
                skip = 1
                continue
            skip = 0
            print(curr_patient)
            print(count)
            if count != 0:
                f.close()
            f = open('data/aliases/exp/' + curr_patient + '_exp_aliases.txt' ,'w')
        if (skip == 0) :
            if gene_data[0] in graph.nodes:
                print(gene_data[0])
                f.write('\n' + gene_data[0])
                count += 1
            else:
                gene_true_name = checkAliases(gene_data[0], graph.nodes)
                if gene_true_name:
                    f.write('\n' + gene_true_name)
                    count += 1
                #else:
                #   print('Could not find suitable alias for', gene_data[0])
    print('Loaded ' + str(count) + ' differentially expressed genes.')
    file.close()
    f.close()
    return count
    
