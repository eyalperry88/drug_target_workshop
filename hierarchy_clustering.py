import numpy as np
import os
from scipy.spatial.distance import pdist
from scipy.cluster.hierarchy import linkage, dendrogram
import matplotlib.pyplot as plt

DIR = 'output_b2h'

# load categories
f_categories = open('data/AML_categories.txt', 'r')
categories = {}
for line in f_categories:
    row_data = line.strip().split('\t')
    categories[row_data[0]] = row_data[1]
f_categories.close()

# first, get the top k results from each patient
consensus_genes = []
k = 150
filenames = os.listdir(DIR)
for filename in filenames:
    f = open(DIR + '/' + filename, 'r')
    i = 0
    for line in f:
        if i == k:
            break
        row_data = line.strip().split('\t')
        if row_data[0] not in consensus_genes:
            consensus_genes.append(row_data[0])
        i += 1
    f.close()

patient_files = []
subtypes = ['M3', 'M5']
label_colors = {'M3': 'r', 'M5': 'g', 'M6' : 'b', 'M7' : 'm'}
filenames = os.listdir(DIR)
for filename in filenames:
    patient = filename[:12]
    if categories[patient] not in subtypes:
        continue
    patient_files.append(filename)

num_samples = len(patient_files)
print('Patients: ', num_samples)
num_genes = len(consensus_genes)
print('Genes: ', len(consensus_genes))
# now create a matrix whos rows are the patients and columns are the genes
X = np.zeros(shape=(num_samples, num_genes))
labels = []
i = 0
for filename in patient_files:
    patient = filename[:12]
    labels.append(categories[patient])
        
    f = open(DIR + '/' + filename, 'r')
    j = 0
    top_k_genes = []
    for line in f:
        if j == k:
            break
        row_data = line.strip().split('\t')
        top_k_genes.append(row_data[0])
        j += 1
    f.close()
    
    for j in range(num_genes):
        if consensus_genes[j] in top_k_genes:
            X[i, j] = 1
    i += 1

Y = pdist(X, metric='jaccard')

Z = linkage(Y, method='complete', metric='jaccard')

plt.clf()
dendrogram(Z, labels=labels, orientation='left', color_threshold=0)
ax = plt.gca()
xlbls = ax.get_ymajorticklabels()
for lbl in xlbls:
    if lbl.get_text() in label_colors:
        lbl.set_color(label_colors[lbl.get_text()])
plt.savefig("test.png", dpi = 800)