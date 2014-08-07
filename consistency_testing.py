#!/usr/bin/python

from feast import *
import numpy as np
from sklearn import naive_bayes
from sklearn import cross_validation as xval
from sklearn.svm import SVC
from sklearn.svm import LinearSVC
from sklearn.preprocessing import StandardScaler
from sklearn.grid_search import GridSearchCV


import output_missing_data_matrix as omdm
import data_analysis as da
import utils
import process_results_otus as pro

import json 

import argparse
# from biom.parse import parse_biom_table 
import os
import sys
import pdb, traceback

__author__ = "Calvin Morrison"
__copyright__ = "Copyright 2013, EESI Lab"
__credits__ = ["Calvin Morrison", "Gregory Ditzler", "Jean-Luc Bouchot"]
__license__ = "GPL"
__version__ = "0.1.0-dev"
__maintainer__ = "Calvin Morrison"
__email__ = "mutantturkey@gmail.com"
__status__ = "Development"


def consistency(list_a, list_b, larger_set_size):
	idx1, idx2 = da.find_common_samples(list_a, list_b)
	# Note: We should make sure that both lists have the same size
	k = len(list_a)
	c = float(len(idx1)*larger_set_size-k*k)/float(k*(larger_set_size-k))
	
	return c
	
	
def launch_tests_feature_selection(biom_matrix, labels, site_names, otu_ids, otu_phylo, obj_function='mim', nb_feat_list = range(50,150), nb_rdm_tests=10, output_folder = '.', k_fold=10, obj_fun_parameter = [], compare_classifier=True):
	# go over all feature sizes
	# go over all the random tests
	# Save the lowest, highest, and average consistency 
	
	# Replace nb_feat_list with max_nb_feat
	# We can also increase nb_rdm_tests
	
	max_nb_feat = np.max(nb_feat_list) # Might need to check that we are dealing with numpy arrays, and not lists
	
	otu_fold_nb = 0
	skf = xval.StratifiedShuffleSplit(labels,nb_rdm_tests,1.0/k_fold)
	cur_folder = os.path.join(output_folder, 'tmp')
	if not os.path.isdir(cur_folder):
		os.makedirs(cur_folder);
		
	
	for train, test in skf:
		train_set = biom_matrix[train,:]
		test_set = biom_matrix[test,:]
		train_labels = labels[train]
		relevant_otus = da.find_relevant_otus(train_set, train_labels, obj_function, max_nb_feat, obj_fun_parameter)
		da.save_found_otus(relevant_otus, otu_ids, otu_phylo, cur_folder, 'otu_fold' + str(otu_fold_nb) + '.txt')
		
		with open(os.path.join(cur_folder, 'samples_fold'+str(otu_fold_nb) + '.txt'), 'w') as f:
			for item in test:
				print>>f, item
		
		otu_fold_nb = otu_fold_nb+1
		
	
	if compare_classifier:
	
		test_accuracies(biom_matrix, labels, site_names, otu_ids, cur_folder, nb_feat_list, nb_rdm_tests, k_fold)
	
#		otu_fold_nb = 0
#		
#		for train, test in skf:
#			
#			print otu_fold_nb
#			otu_fname = os.path.join(cur_folder, 'otu_fold' + str(otu_fold_nb) + '.txt')
#			# First read the whole set of OTUs
#			whole_otu_set, otu_ids = pro.read_relevant_otus(otu_fname)
#			for idx, val in enumerate(nb_feat_list):
#				relevant_otus = whole_otu_set[0:val]
#				train_set = biom_matrix[train,:]
#				test_set = biom_matrix[test,:]
#				train_set = train_set[:,relevant_otus]
#				test_set = test_set[:,relevant_otus]
#		
#				cur_nbc = da.learn_nbc(train_set, train_labels)
#				predicted_output = cur_nbc.predict(test_set)
#				
#				cur_error_mnbc = cur_nbc.score(test_set, labels[test])
#				
#				# Gaussian NBC 
#				gnbc = naive_bayes.GaussianNB()
#				gnbc.fit(train_set, train_labels)
#				predicted_output_g = gnbc.predict(test_set)
#				cur_error_mnbc = gnbc.score(test_set, labels[test])
#				
#				# SVM - RBF kernel
#				C_range = 10.0 ** np.arange(-2, 4)
#				gamma_range = 10.0 ** np.arange(-5, 4)
#				param_grid = dict(gamma=gamma_range, C=C_range)
#
#				
#				cv = xval.StratifiedKFold(y=train_labels, n_folds=5)
#				grid = GridSearchCV(SVC(cache_size=800), param_grid=param_grid, cv=cv)
#				grid.fit(train_set, train_labels)
#				
#				
#				predicted_output_svm = grid.predict(test_set)
#				
#				cur_error_rbf = grid.score(test_set, labels[test])
#				# Still have to write the optimal parameters!
#				f = open(os.path.join(cur_folder,"n"+str(val)+"test" + str(otu_fold_nb) + "_optimalrbf.txt"),"w")
#				print >> f, grid.best_estimator_
#				f.close()
#				
#				# Linear classifier
#				linear_classifier = LinearSVC()
#				linear_classifier.fit(train_set, train_labels)
#				predicted_output_linear = linear_classifier.predict(test_set)
#				
#				cur_error_linear = linear_classifier.score(test_set, labels[test])
#				
#				# Print in a file:
#				# SampleID \t True Label \t nbcMN Prediction \t gnbc prediction
#				# We need to make sure we save the whole phylogeny of the OTU
#				f = open(os.path.join(cur_folder, "n"+str(val)+"test" + str(otu_fold_nb) + ".txt"),"w")
#				f.write("SampleID\tLabel\tMultinomial\tGaussianNBC\tLinear\tRBF\n")
#				for idx, l in enumerate(test):
#					f.write(site_names[l] + "\t" + str(labels[l]) + "\t" + str(predicted_output[idx]) + "\t" + str(predicted_output_g[idx]) + "\t" + str(predicted_output_linear[idx]) + "\t" + str(predicted_output_svm[idx]) + "\n")
#				f.close() 
#				
#			otu_fold_nb = otu_fold_nb + 1
#				
				
				
				
			

def test_accuracies(biom_matrix, labels, site_names, otu_ids, cur_folder, nb_feat_list = range(50,150), nb_rdm_tests=10, k_fold=10): # Later we will add a parameter for the choice of the classifier
	
	
	for idx, val in enumerate(nb_feat_list):
		otu_fold_nb = 0
		
		
		for otu_fold_nb in xrange(0,k_fold):
			
			# Read the elements in the test sets
			test_set_fname = os.path.join(cur_folder, 'samples_fold' + str(otu_fold_nb) + '.txt')
			with open(test_set_fname, 'r') as f:
				test = [int(line.rstrip('\n')) for line in f]
			
			#print test
			train = [i for i in xrange(0,biom_matrix.shape[0]) if i not in test]
			
			otu_fname = os.path.join(cur_folder, 'otu_fold' + str(otu_fold_nb) + '.txt')
			whole_otu_set, otu_ids = pro.read_relevant_otus(otu_fname)
			
			relevant_otus = whole_otu_set[0:val]
			train_set = biom_matrix[train,:]
			test_set = biom_matrix[test,:]
			train_set = train_set[:,relevant_otus]
			test_set = test_set[:,relevant_otus]
			train_labels = labels[train]
		
			cur_nbc = da.learn_nbc(train_set, train_labels)
			predicted_output = cur_nbc.predict(test_set)
				
			cur_error_mnbc = cur_nbc.score(test_set, labels[test])
				
			# Gaussian NBC 
			gnbc = naive_bayes.GaussianNB()
			gnbc.fit(train_set, train_labels)
			predicted_output_g = gnbc.predict(test_set)
			cur_error_mnbc = gnbc.score(test_set, labels[test])
				
			# SVM - RBF kernel
			C_range = 10.0 ** np.arange(-2, 4)
			gamma_range = 10.0 ** np.arange(-5, 4)
			param_grid = dict(gamma=gamma_range, C=C_range)

				
			cv = xval.StratifiedKFold(y=train_labels, n_folds=5)
			grid = GridSearchCV(SVC(cache_size=800, class_weight='auto'), param_grid=param_grid, cv=cv)
			grid.fit(train_set, train_labels)
				
				
			predicted_output_svm = grid.predict(test_set)
				
			cur_error_rbf = grid.score(test_set, labels[test])
			# Still have to write the optimal parameters!
			f = open(os.path.join(cur_folder,"n"+str(val)+"test" + str(otu_fold_nb) + "_optimalrbf.txt"),"w")
			print >> f, grid.best_estimator_
			f.close()
				
			# Linear classifier
			linear_classifier = LinearSVC(class_weight='auto')
			linear_classifier.fit(train_set, train_labels)
			predicted_output_linear = linear_classifier.predict(test_set)
				
			cur_error_linear = linear_classifier.score(test_set, labels[test])
				
			# Print in a file:
			# SampleID \t True Label \t nbcMN Prediction \t gnbc prediction
			# We need to make sure we save the whole phylogeny of the OTU
			f = open(os.path.join(cur_folder, "n"+str(val)+"test" + str(otu_fold_nb) + ".txt"),"w")
			f.write("SampleID\tLabel\tMultinomial\tGaussianNBC\tLinear\tRBF\n")
			for idx, l in enumerate(test):
				f.write(site_names[l] + "\t" + str(labels[l]) + "\t" + str(predicted_output[idx]) + "\t" + str(predicted_output_g[idx]) + "\t" + str(predicted_output_linear[idx]) + "\t" + str(predicted_output_svm[idx]) + "\n")
			f.close() 
			
		otu_fold_nb	 = otu_fold_nb + 1
	

def get_consistencies(set_sizes, larger_set_size, nb_tests, output_folder):
	
	avg_consistency = {}
	max_consistency = {}
	min_consistency = {}
	std_consistency = {}
	
	for one_set_size in set_sizes:
	
		consistencies = []
		print 'Consistencies with set size ' + str(one_set_size)
	
		# avg_consistency[one_set_size] = 0
		# max_consistency[one_set_size] = 0
		# min_consistency[one_set_size] = 10000
		
		tests_considered = 0; # The grand total should be nfold*nb_tests*(nfold*nb_tests-1)/2
		# look up in the containing folder 
		# files are in for instance in tmp/n150/rdmtest1/otus_fold2.txt
		for one_fold in range(0,nb_tests-1):
			# open this file and look at the consistency with all the other tests of that size
			path_to_file_a = os.path.join(output_folder, 'tmp', 'n'+str(one_set_size))
			fname_a =  'otu_fold' + str(one_fold) + '.txt'
			
			# print 'list A taken from: ' + os.path.join(path_to_file_a, fname_a)
			
			# /!\ Define list_a here
			list_a = []
			with open(os.path.join(path_to_file_a, fname_a), 'r') as file_a:
				for i,row in enumerate(file_a):
					# -1 to get rid of the newline
					split_row = row[:-1].split("\t")
					# list_a.append(int(split_row[0]))
					list_a.append(split_row[0])
				
				
			
			for another_rdm_test in range(one_fold+1,nb_tests):
				# /!\ Define list_b here
				path_to_file_b = os.path.join(output_folder, 'tmp')
				fname_b =  'n'+str(one_set_size) + 'otu_fold' + str(another_rdm_test) + '.txt'
				
				# print 'list B taken from: ' + os.path.join(path_to_file_b, fname_b)
				
				list_b = []
				with open(os.path.join(path_to_file_b, fname_b), 'r') as file_b:
					for i,row in enumerate(file_b):
						# -1 to get rid of the newline
						split_row = row[:-1].split("\t")
						# list_b.append(int(split_row[0]))
						list_b.append(split_row[0])
				
				
				
				cur_consistency = consistency(list_a, list_b, larger_set_size)
				consistencies.append(cur_consistency)
				# max_consistency[one_set_size] = np.max([max_consistency[one_set_size],cur_consistency])
				# min_consistency[one_set_size] = np.min([min_consistency[one_set_size],cur_consistency])
				# avg_consistency[one_set_size] = avg_consistency[one_set_size]+cur_consistency
				tests_considered = tests_considered+1
		
		avg_consistency[one_set_size] = sum(consistencies)/float(len(consistencies))
		max_consistency[one_set_size] = max(consistencies)
		min_consistency[one_set_size] = min(consistencies)
		print consistencies
		std_consistency[one_set_size] = np.std(np.array(consistencies))
		# max_consistency[one_set_size] = max_consistency[one_set_size]
		# min_consistency[one_set_size] = min_consistency[one_set_size]
		# avg_consistency[one_set_size] = avg_consistency[one_set_size]/tests_considered
		
		
	return avg_consistency, max_consistency, min_consistency, std_consistency
				



def get_accuracies(set_sizes, larger_set_size, nb_tests, output_folder):
	
	avg_accuracy_g = {}
	max_accuracy_g = {}
	min_accuracy_g = {}
	std_accuracy_g = {}
	
	avg_accuracy = {}
	max_accuracy = {}
	min_accuracy = {}
	std_accuracy = {}
	
	
	for one_set_size in set_sizes:
		# avg_accuracy_g[one_set_size] = 0
		# max_accuracy_g[one_set_size] = 0
		# min_accuracy_g[one_set_size] = 10000
		
		# avg_accuracy[one_set_size] = 0
		# max_accuracy[one_set_size] = 0
		# min_accuracy[one_set_size] = 10000
		print 'Accuracies for a set size of ' + str(one_set_size)
		
		acc_g = []
		acc_nb = []
		
		tests_considered = 0; # The grand total should be nfold*nb_tests*(nfold*nb_tests-1)/2
		# look up in the containing folder 
		# files are in for instance in tmp/n150/rdmtest1/otus_fold2.txt
		for one_fold in range(0,nb_tests-1):
			# open this file and look at the consistency with all the other tests of that size
			path_to_file_a = os.path.join(output_folder, 'tmp')
			fname_a =  'n'+str(one_set_size)+'test' + str(one_fold) + '.txt'
						
			sample_names = []
			gaussian_prediction = []
			mn_prediction = []
			true_label = []
			
			with open(os.path.join(path_to_file_a, fname_a), 'r') as fh:
				next(fh)
				for i,row in enumerate(fh):
					# -1 to get rid of the newline
					split_row = row[:-1].split("\t")
					# Make sure we have 3 rows.
					sample_names.append(split_row[0])
					true_label.append(split_row[1])
					gaussian_prediction.append(split_row[2])
					mn_prediction.append(split_row[3])

			
			labels_equal_g = [1 if true_label[i]==gaussian_prediction[i] else 0 for i in xrange(0,len(true_label))]
			#for i in xrange(0,len(true_label)):
			#	print 'True label: {} - Gaussian Prediction: {} - Isequal? {}'.format(str(true_label[i]), str(gaussian_prediction[i]), true_label[i]==gaussian_prediction[i])
			labels_equal = [1 if true_label[i]==mn_prediction[i] else 0 for i in xrange(0,len(true_label))]
			# print labels_equal_g
			# max_accuracy_g[one_set_size] = np.max([max_accuracy_g[one_set_size],sum(labels_equal_g)/float(len(true_label))])
			# min_accuracy_g[one_set_size] = np.min([min_accuracy_g[one_set_size],sum(labels_equal_g)/float(len(true_label))])
			# avg_accuracy_g[one_set_size] = avg_accuracy_g[one_set_size]+sum(labels_equal_g)/float(len(true_label))
			
			# max_accuracy[one_set_size] = np.max([max_accuracy[one_set_size],sum(labels_equal)/float(len(true_label))])
			# min_accuracy[one_set_size] = np.min([min_accuracy[one_set_size],sum(labels_equal)/float(len(true_label))])
			# avg_accuracy[one_set_size] = avg_accuracy[one_set_size]+sum(labels_equal)/float(len(true_label))
			
			acc_g.append(sum(labels_equal_g)/float(len(true_label)))
			acc_nb.append(sum(labels_equal)/float(len(true_label)))
			
			tests_considered = tests_considered+1
				
				
		
		avg_accuracy_g[one_set_size] = sum(acc_g)/float(len(acc_g))
		max_accuracy_g[one_set_size] = max(acc_g)
		min_accuracy_g[one_set_size] = min(acc_g)
		std_accuracy_g[one_set_size] = np.std(np.array(acc_g))
		
		avg_accuracy[one_set_size] = sum(acc_nb)/float(len(acc_nb))
		max_accuracy[one_set_size] = max(acc_nb)
		min_accuracy[one_set_size] = min(acc_nb)
		std_accuracy[one_set_size] = np.std(np.array(acc_nb))
		
		# max_accuracy_g[one_set_size] = max_accuracy_g[one_set_size]
		# min_accuracy_g[one_set_size] = min_accuracy_g[one_set_size]
		# avg_accuracy_g[one_set_size] = avg_accuracy_g[one_set_size]/tests_considered
		
		# max_accuracy[one_set_size] = max_accuracy[one_set_size]
		# min_accuracy[one_set_size] = min_accuracy[one_set_size]
		# avg_accuracy[one_set_size] = avg_accuracy[one_set_size]/tests_considered
		
		
	return avg_accuracy_g, max_accuracy_g, min_accuracy_g, std_accuracy_g, avg_accuracy, max_accuracy, min_accuracy, std_accuracy
				


def save_results(vartype, fname, avg_value, max_val, min_val, std_val):
	# Note: avg, min, max, and std are dictionaries containing the same keys: the number of features
	# dict_keys = avg_value.keys()
	with open(fname,'w') as fh:
		fh.write('NbFeatures\tAvg'+vartype+'\tMin'+vartype+'\tMax'+vartype+'\tStd'+vartype+'\n')
		for a_key in avg_value.keys():
			fh.write(str(a_key)+'\t'+str(avg_value[a_key])+'\t'+str(min_val[a_key])+'\t'+str(max_val[a_key])+'\t'+str(std_val[a_key])+'\n')
	return None
	


def main():
	parser = argparse.ArgumentParser(description = "")
	parser.add_argument("-b", "--biom-file", help="An input biom file", required=True)
	parser.add_argument("-m", "--mapping-file", help="A mapping file", required=True)
	parser.add_argument("-c", "--class-label", help="Which data are we trying to analyze", required=True)
	parser.add_argument("-d", "--subclass", action="append", help="Subselect only some of the data - if specified, this should appear at least twice with the adequate options. ex: -c SEX -d male -d female", required=False)
	parser.add_argument("-o", "--output-folder", help="The folder to output our data to", required=True)
	parser.add_argument("-p", "--min-features", help="Minimum number of features to test", default=50, required=False)
	parser.add_argument("-q", "--max-features", help="Maximum number of features to test", default=150, required=False)
	parser.add_argument("-s", "--step-size", help="Step size within the range of the number of features to be tested", default=1, required=False)
	#parser.add_argument("-p", "--predictor", help="Classifier/Predictor used", default="nbc", required=False) # As of today, contains only nbc
	parser.add_argument("-j", "--objective-function", help="Objective function for the feature selection algorithm", default="mim", required=False)
	parser.add_argument("-t", "--output-type", help="data output format. default: CSV options: csv, matlab, r, numpy", default="csv", required=False)
	parser.add_argument("-f", "--select-field", help="Field to extract a subset of the data. e.g. EN_BIOME, COUNTRY. The default considers the whole dataset", default=None, required=False)
	parser.add_argument("-g", "--value-field", action="append", help="When used with -f specifies the value of the field to filter - THIS IS REQUIRED IF -f if present", default=None, required=False)
	parser.add_argument("-k", "--cluster", action="append", help="Allows to subgroup some of the labels. Ex: -k 'Vegan Vegan+Seafood'. The different values are separated with semi colon. Requires at least two appearances. This cannot be used in conjunction with the -d option", default=None, required=False)
## Need to be continued!!!!
	print "Definition of the arguments done"
	global output_type
	print "Start of the program"

	args = parser.parse_args()
	
	output_type = args.output_type.lower()
		
	# if our folder doesn't exist create it
	if not os.path.isdir(args.output_folder):
		os.mkdir(args.output_folder);

	nb_features = range(int(args.min_features), int(args.max_features)+1, int(args.step_size))
	print "nb_features prepared"
	
	matrix, site_names, otu_ids, otu_phylo = utils.load_biom(args.biom_file)
	metadata = utils.load_map(args.mapping_file)
	class_labels = []
	for sample in site_names:
		class_labels.append(metadata[sample][args.class_label])
	
	print "class_labels loaded"
	
	interesting_samples = range(0,len(site_names))
	
	if args.select_field is not None:
		interesting_fields = [it.lower() for it in args.value_field]
		print interesting_fields
		subsample_habitat = [i for i,sample in enumerate(site_names) if metadata[sample][args.select_field].lower() in interesting_fields]
		interesting_samples = list(set(interesting_samples).intersection(set(subsample_habitat)))
	
	if args.subclass is not None:
		target_labels = [it.lower() for it in args.subclass]
		subsamples = [i for i in xrange(0,len(class_labels)) if class_labels[i].lower() in target_labels]
		interesting_samples = list(set(interesting_samples).intersection(set(subsamples)))
	
	if (args.cluster is not None) and (args.subclass is None):
		print "In da cluster separation"
		clusters = [it for it in args.cluster]
		clusters_dict = {}
		print "Initial Dictionary created"
		for idx, a_cluster in enumerate(clusters):
			print "In da loop"
			#keys = a_cluster.split()
			keys = a_cluster.split(';')
			print keys
			for a_key in keys:
				clusters_dict[a_key.lower()] = idx
		
		subsamples = [i for i in xrange(0,len(class_labels)) if class_labels[i].lower() in clusters_dict]
		interesting_samples = list(set(interesting_samples).intersection(set(subsamples)))		
		for i in subsamples:
			class_labels[i] = 'cluster'+str(clusters_dict[class_labels[i].lower()])
	
	matrix = matrix[interesting_samples,:]
	class_labels = [class_labels[i] for i in interesting_samples]
	
	class_labels,labels_key = utils.discretize(class_labels)
	
	matrix = matrix+1 
	row_sums = matrix.sum(axis=1)
	
	matrix = matrix / row_sums[:, np.newaxis]
	matrix = np.ceil(matrix/matrix.min())
	
	
	
	
	# So far, we have the biom file open and the environment parameters
	# We can now launch our feature selection algorithm
	further_param = [] # This has to be adapted to the case we are using other objective functions
	
	
	nb_tests = 10
	nb_folds = 5
		
	launch_tests_feature_selection(matrix, np.array(map(int,class_labels)), site_names, otu_ids, otu_phylo, args.objective_function, nb_features, nb_tests, args.output_folder, nb_folds)
	avg_consistency, max_consistency, min_consistency, std_consistency = get_consistencies(nb_features, len(otu_ids), nb_tests, args.output_folder)
	save_results('consistency', os.path.join(args.output_folder,'consistencyresults.txt'), avg_consistency,max_consistency,min_consistency,std_consistency)
	
	
	avg_accuracy_g, max_accuracy_g, min_accuracy_g, std_accuracy_g, avg_accuracy, max_accuracy, min_accuracy, std_accuracy = get_accuracies(nb_features, len(otu_ids), nb_tests, args.output_folder)
	save_results('Accuracy', os.path.join(args.output_folder,'accuracyGaussianresults.txt'), avg_accuracy_g,max_accuracy_g,min_accuracy_g,std_accuracy_g)

		


if __name__ == "__main__":
	try:
		sys.exit(main())
	except:
		type, value, tb = sys.exc_info()
		traceback.print_exc()
		pdb.post_mortem(tb)
