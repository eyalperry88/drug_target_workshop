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
label_colors = {'M0' : 'r', 'M1' : 'g', 'M2' : 'b', 'M3' : 'c', 'M4' : 'm', 'M5' : 'y'}
subtypes = list(label_colors.keys()) + ['M6', 'M7']
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
dendrogram(Z, labels=labels, orientation='left', color_threshold=0, leaf_font_size=6)
ax = plt.gca()
xlbls = ax.get_ymajorticklabels()
for lbl in xlbls:
    if lbl.get_text() in label_colors:
        lbl.set_color(label_colors[lbl.get_text()])
plt.savefig("test.png", dpi = 800)


dists_per_type = {'M0' : [], 'M1' : [], 'M2' : [], 'M3' : [], 'M4' : [], 'M5' : [],
                 'M6' : [], 'M7' : []}
dists_between_M3 = []
dists_M3_and_others = []
dists_others = []
dists_between = []
dists_among = []
idx = 0
for i in range(len(labels) - 1):
    for j in range(i + 1, len(labels)):
        if labels[i] == labels[j]:
            if labels[i] == 'M3':
                dists_between_M3.append(Y[idx])
            else:
                dists_M3_and_others.append(Y[idx])
            dists_per_type[labels[i]].append(Y[idx])
            dists_between.append(Y[idx])
        else:
            if labels[i] == 'M3' or labels[j] == 'M3':
                dists_M3_and_others.append(Y[idx])
            else:
                dists_others.append(Y[idx])
            dists_among.append(Y[idx])
        idx += 1
mean_per_type = {}
for subtype in dists_per_type:
    mean_per_type[subtype] = np.mean(dists_per_type[subtype])
var_per_type = {}
for subtype in dists_per_type:
    var_per_type[subtype] = np.var(dists_per_type[subtype])

plt.clf()
plt.hist(dists_between, bins=18, histtype='step', color='b', normed=True, label='Between Same Type')
plt.hist(dists_among, bins=18, histtype='step', color='r', normed=True, label='Among Different Types')
plt.hist(dists_between_M3, bins=18, histtype='step', color='y', normed=True, label='Between M3 Type')
plt.hist(dists_M3_and_others, bins=18, histtype='step', color='g', normed=True, label='M3 Type vs. Other Types')
plt.title("Distributions of Similiarities")
plt.legend(loc=2)
plt.savefig("hist.png", dpi = 800)

np.random.shuffle(Y)

dists_between_M3 = []
dists_M3_and_others = []
dists_others = []
dists_between = []
dists_among = []
idx = 0
for i in range(len(labels) - 1):
    for j in range(i + 1, len(labels)):
        if labels[i] == labels[j]:
            if labels[i] == 'M3':
                dists_between_M3.append(Y[idx])
            else:
                dists_M3_and_others.append(Y[idx])
            dists_between.append(Y[idx])
        else:
            if labels[i] == 'M3' or labels[j] == 'M3':
                dists_M3_and_others.append(Y[idx])
            else:
                dists_others.append(Y[idx])
            dists_among.append(Y[idx])
        idx += 1

plt.clf()
plt.hist(dists_between, bins=18, histtype='step', color='b', normed=True, label='Between Same Type')
plt.hist(dists_among, bins=18, histtype='step', color='r', normed=True, label='Among Different Types')
plt.hist(dists_between_M3, bins=18, histtype='step', color='y', normed=True, label='Between Pseudo-M3 Type')
plt.hist(dists_M3_and_others, bins=18, histtype='step', color='g', normed=True, label='Pseudo-M3 Type vs. Other Types')
plt.title("Distributions of Random Similiarities")
plt.legend(loc=2)
plt.savefig("hist_rand.png", dpi = 800)
