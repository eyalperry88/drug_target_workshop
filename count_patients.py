patient_data = {}
file = open("data/AML_Expression.txt", 'r')
next(file) #first line has no data
for line in file:
    gene_data = line.split('\t')
    if gene_data[1] in patient_data:
        patient_data[gene_data[1]]['expression'] += 1
    else:
        patient_data[gene_data[1]] = {'expression' : 1, 'mutation' : 0}
file.close()
file = open("data/AML_Mutations.txt", 'r')
next(file) #first line has no data
for line in file:
    gene_data = line.split('\t')
    if gene_data[1] in patient_data:
        patient_data[gene_data[1]]['mutation'] += 1
    else:
        patient_data[gene_data[1]] = {'expression' : 0, 'mutation' : 1}
file.close()

count_patients = {'both' : 0, 'mutation_only' : 0, 'expression_only' : 0}
for patient in patient_data:
    if patient_data[patient]['expression'] == 0:
        count_patients['mutation_only'] += 1
    elif patient_data[patient]['mutation'] == 0:
        count_patients['expression_only'] += 1
    else:
        count_patients['both'] += 1
print(str(count_patients['both']) + ' patients with sufficient data')
print(str(count_patients['expression_only']) + ' patients with only expression data')
print(str(count_patients['mutation_only']) + ' patients with only mutation data')