import matplotlib.pyplot as plt
import numpy as np
import mpld3
from mpl_toolkits.axes_grid1 import make_axes_locatable



### Read in file, store gene pairs in edges, genes in nodes, and weights in weights
edges = []
nodes = []
weights = []
for item in file('../output/gene_pairs_count_chr1_s0.05_EUR.txt') :
	item = item.split('\t')
	if int(item[1]) > 5 :
		genes = item[0].replace('(', '').replace(')', '').replace("'", '').split(', ')
		edges.append((genes[0], genes[1]))
		nodes.append(genes[0])
		nodes.append(genes[1])
		weights.append(float(item[1]))

print len(edges), 'pairs of genes listed'

### Create dictionary of all genes with how often they are co-mutated
all_genes = {}
for ii, gene in enumerate(nodes) :
	if gene in all_genes :
		all_genes[gene]['weight'] += weights[ii/2]
	else : 
		all_genes[gene] = {}
		all_genes[gene]['weight'] = weights[ii/2]

### Rank order genes by mutation frequency and add to dictionary
ordered_genes = []	
for w in sorted(all_genes, key=all_genes.get, reverse=True):
	ordered_genes.append(w)
for ii, gene in enumerate(ordered_genes) :
	all_genes[gene]['rank'] = ii
  
### Create array of gene mutation pairs
n_genes = len(all_genes)
gene_array = np.zeros((n_genes, n_genes))

for ii, edge in enumerate(edges) :
	xx = all_genes[edge[0]]['rank']
	yy = all_genes[edge[1]]['rank']
	if xx < yy :
		gene_array[xx][yy] += weights[ii]
	else : 
		gene_array[yy][xx] += weights[ii]

### Create array of labels for gene mutation pairs
x_points = [n_genes - ii - 1 for ii in range(n_genes) for jj in range(n_genes)]
y_points = [jj for ii in range(n_genes) for jj in range(n_genes)]
point_labels = []
for ii in range(n_genes) :
	for jj in range(n_genes) :
		g1 = ordered_genes[ii]
		g2 = ordered_genes[jj]
		count = str(int(gene_array[jj][ii]))
		point_labels.append('('+ g1 +', '+ g2 +') = ' + count)
print x_points[:10], point_labels[:10]
print x_points[-10:], point_labels[-10:]
point_labels.reverse()
print point_labels[:10]

print len(all_genes), 'unique genes'
max_key = max(all_genes, key=all_genes.get)
print 'gene', max_key, 'has', all_genes[max_key], 'connections'

### Plot array with interactive labels
fig, ax = plt.subplots(figsize=(10,10))
ax.set_xticks([])
ax.set_yticks([])
plt.title("Gene pairs for European Population", fontsize=30)

scatter = ax.scatter(x_points, y_points, s=5,
	alpha=0.1, cmap=plt.cm.jet, marker = 's')

divider = make_axes_locatable(ax)

tooltip = mpld3.plugins.PointLabelTooltip(scatter, labels=point_labels)

A = plt.imshow(np.log(gene_array), interpolation='nearest', origin = 'lower')
mpld3.plugins.connect(fig, tooltip)

mpld3.show()
