import time

tic = time.clock()
import numpy as np
from random import sample
import os.path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cm
from mpl_toolkits.axes_grid1 import make_axes_locatable
import itertools
import argparse

c_help = 'type a chromosome 1-22'
pop_help = 'type a population 0-5; 0:ALL, 1:EUR, 2:EAS, 3:AFR, 4:AMR, 5:SAS'
description = 'Process mutation sets (-c and -POP are required).'
parser = argparse.ArgumentParser(description = description)
parser.add_argument("-c", type=int,
                    help=c_help)
parser.add_argument("-pop", type=int,
                    help=pop_help)
args = parser.parse_args()
c = args.c

dataset = '20130502'
SIFT = 0.05
n_runs = 100
pops = ('ALL', 'EUR', 'EAS', 'AFR', 'AMR', 'SAS')
siftfile = '../data/' + dataset + '/sifted.SIFT.chr' + str(c) + '.txt'
data_dir = '../data/' + dataset + '/'
pop_dir = '../populations/'
outdata_dir = '../output/'
plot_dir = '../plots/'
OutputFormat = '.png'

POP = pops[args.pop]
chrom = 'chr' + str(c)

font = {'family':'serif',
	'size':14	}
plt.rc('font', **font)


print time.clock() - tic, 'modules loaded'

tic = time.clock()

class ReadData :
	def read_names(self, POP) :
		tic = time.clock()
		print 'reading in IDs for population', POP
		namefile = pop_dir + '1kg_' + POP + '.txt'
		f = open(namefile, 'r')
		text = f.read()
		f.close()
		text = text.split()
		all_ids = text[0:]
		file = data_dir + 'columns.txt'
		f = open(file, 'r')
		text = f.read()
		f.close()
		genome_ids = text.split()
		
		ids = list(set(all_ids) & set(genome_ids))
		
		print len(ids), 'ids in population', POP, 'out of', len(all_ids), 'listed'
		print time.clock() - tic
		return ids

	def read_sifted_indices(self, siftfile) :
		print 'reading in rs with sift scores below', SIFT
		## NB This file is in the format of:
		## line number, rs number, ENSG number, SIFT, Phenotype
		tic = time.clock()
		sift_indices = []
		genes = {}
		all_genes = []
		for item in file(siftfile, 'r') :
			item = item.split()
			if len(item) == 6 :
				sift = float(item[3])
				if sift < SIFT :
					sift_indices.append(item[0])
					genes[item[0]] = item[2]
					all_genes.append(item[2])
		
		seen = set()			
		[x for x in all_genes if x not in seen and not seen.add(x)]
		all_genes = seen

		print len(sift_indices), 'variants have sift scores'
		print len(all_genes), 'genes are included'
		print 'time: ', time.clock() - tic

		return sift_indices, genes, all_genes

	def read_individuals(self, ids, sift_indices) :
		print 'reading in individual mutation files'
		tic = time.clock()
		mutation_index_array = []
		for name in ids :
			filename = data_dir + chrom + 'n/' + chrom + '.' + name
			f = open(filename, 'r')
			text = f.read()
			f.close()
			text = text.split()
			sifted_mutations = list(set(sift_indices).intersection(text))
			sifted_genes = ()
			sifted_genes = [genes[x] for x in sifted_mutations]
			seen = set()
			[x for x in sifted_genes if x not in seen and not seen.add(x)]
			mutation_index_array.append(list(seen))
#			print len(sifted_genes), len(seen)
		
		print 'mutation index array for', ids[0], ':', mutation_index_array[0]
		print 'time: ', time.clock() - tic
		return mutation_index_array

	def read_pairs_overlap(self, indpairsfile) :
		print 'reading in individual crossover mutations'
		tic = time.clock()
		pairs_overlap = np.loadtxt(indpairsfile, unpack=True)
		pairs_overlap = np.transpose(pairs_overlap)

		print 'time: ', time.clock() - tic
		return pairs_overlap


class Results :

	def pair_individuals(self, mutation_index_array) :
		print 'cross matching mutations in individuals'
		tic = time.clock()
	
		n_p = len(mutation_index_array)
		n_pairs = n_p/2
		list_p = np.linspace(0, n_p - 1, n_p).astype(int)
		print 'calculating individual crossover for', n_p, 'individuals' 
		pairs_overlap = np.zeros((n_runs, n_pairs))
		for run in range(n_runs) :
			randomized_list = sample(list_p, n_p)
			for pq in range(n_pairs) :
				array1 = mutation_index_array[randomized_list[2*pq]]
				array2 = mutation_index_array[randomized_list[2*pq + 1]]
				pair_array = set(array1) & set(array2)
				pairs_overlap[run][pq] = len(pair_array)

		print 'time: ', time.clock() - tic
		return pairs_overlap

	def gene_pairs(self, mutation_index_array) :
		print 'cross matching pairs of genes'

		tic = time.clock()
		n_p = len(mutation_index_array)
		gene_pair_list = {}
		for pp in range(n_p) :	
			pairs = itertools.combinations(mutation_index_array[pp], 2)
			for pair in pairs :
				key = str(pair)
				if key not in gene_pair_list : gene_pair_list[key] = 1
				else : gene_pair_list[key] += 1

		toc = time.clock()
		print 'time: ', toc - tic, 'seconds to cross correlate pairs'
		
		return gene_pair_list

class PlotData :		

	def individual_overlap(self, POP, pairs_overlap) :
		print 'plotting cross matched individuals'
		tic = time.clock()
		
		pairs_overlap = np.array(pairs_overlap)		
		print pairs_overlap.shape

		min_p = np.min(pairs_overlap)
		max_p = np.max(pairs_overlap)
		print 'min, max = ', min_p, max_p
		nbins = int(max_p) + 1
#		n_runs = len(pairs_overlap)


		print np.min(pairs_overlap), np.max(pairs_overlap)
		nbins = int(np.max(pairs_overlap))
		bin_centres = np.linspace(0, nbins, nbins)
		bin_edges = np.linspace(-0.5, nbins + 0.5, nbins + 1)

		fig = plt.figure(frameon=False, figsize=(10, 9))
		ax = fig.add_subplot(111)
		ax.set_xlim([0,35])
		hists = []
		max_h = 0
		for run in range(n_runs) :
			h, edges = np.histogram(pairs_overlap[run], bins = bin_edges)
			ax.plot(bin_centres, h, alpha = 0.5)
			max_h = max(max_h, max(h))

		plt.xlabel('Number of overlapping gene mutations', fontsize = 24)
		plt.ylabel(r'frequency', fontsize = 28)
		text1 = 'population ' + POP + '\n' +\
			'chromosome ' + str(c) + '\n' + \
			'SIFT < ' + str(SIFT) + '\n' + \
			str(n_runs) + ' runs'
		plt.text(.95, .95, text1, fontsize = 24, 
			verticalalignment='top', horizontalalignment='right',
			transform = ax.transAxes)

		outputFile = plot_dir + 'pair_distribution_c' + str(c) + '_s' + \
			str(SIFT) + '_' + POP + OutputFormat
		plt.savefig(outputFile)  
		print 'Saved file to', outputFile
		plt.close()
		print time.clock() - tic


class WriteData :
	def write_pair_individuals(self, indpairsfile, pairs_overlap) :	
		print 'writing mutation index array to', indpairsfile
		tic = time.clock()
		np.savetxt(indpairsfile, pairs_overlap, fmt = '%i')
		print time.clock() - tic
		
	def write_gene_pairs(self, genepairsfile, gene_pair_list) :
		print 'writing gene pair list to', genepairsfile
		f = open(genepairsfile, 'w')
		for key, count in gene_pair_list.iteritems() :
			f.write(key + '\t' + str(count) + '\n')
		f.close()

############################################################



if __name__ == '__main__':

	rd = ReadData()
	res = Results()
	wr = WriteData()
	pd = PlotData()
	
	indpairsfile = outdata_dir + 'individual_pairs_overlap_chr' + str(c) + '_s' + \
		str(SIFT) + '_' + POP + '.txt'
	genepairsfile = outdata_dir + 'gene_pairs_count_chr' + str(c) + '_s' + \
		str(SIFT) + '_' + POP + '.txt'

	print time.clock() - tic


	ids = rd.read_names(POP)
	n_pairs = len(ids)/2
	print n_pairs, 'unique pairs of individuals'
	
	
	sift_indices, genes, all_genes = rd.read_sifted_indices(siftfile)
	mutation_index_array = rd.read_individuals(ids, sift_indices)

	pairs_overlap = res.pair_individuals(mutation_index_array)
	gene_pair_list = res.gene_pairs(mutation_index_array)
	
	
	wr.write_pair_individuals(indpairsfile, pairs_overlap)
	wr.write_gene_pairs(genepairsfile, gene_pair_list)

#	pairs_overlap = rd.read_pairs_overlap(indpairsfile)
	pd.individual_overlap(POP, pairs_overlap)
	