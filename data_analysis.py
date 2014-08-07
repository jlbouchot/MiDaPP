#!/usr/bin/python

from feast import *
import numpy as np
from sklearn import naive_bayes
import output_missing_data_matrix as omdm

import argparse
#from biom.parse import parse_biom_table 
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

def find_common_samples(biom_sample_ids, mapping_sample_ids):
	idx_from_biom = [idx for idx, val in enumerate(biom_sample_ids) if val in mapping_sample_ids]
	idx_from_mapping = [idx for idx, val in enumerate(mapping_sample_ids) if val in biom_sample_ids]
	
	
	return idx_from_biom, idx_from_mapping
	
def save_found_otus(otu_idx, otu_ids, otu_phylo, output_folder, fname="otunames.txt"):
	
	# We need to make sure we save the whole phylogeny of the OTU
	f = open(os.path.join(output_folder,fname),"w")
	for l in otu_idx:
		f.write(str(int(l)) + '\t' + str(otu_ids[int(l)])+ "\t" + otu_phylo[otu_ids[int(l)]] + "\n")
	f.close()
	
	# save_found_otus(relevant_features, variable_names, args.output_folder)
	
	return

def find_relevant_otus(biom_matrix, class_labels, obj_function, n_select, further_param=[]):
## Need a parser for the relevant parameters in each case.

	data = biom_matrix
	# relevant_features = MIM(data, class_labels, n_select );
	# print "As an input we got the following option: -", obj_function.lower()
	if obj_function.lower() == 'jmi':
		print "Let's try JMI"
	 	relevant_features = JMI(data, class_labels, n_select)
	elif obj_function.lower() == 'mim':
		print "Let's try MIM"
		relevant_features = MIM(data, class_labels, n_select)
	else:
		# For now, as default, we will use MIM
		print "Let's try MIM as a default option"
		relevant_features = MIM(data, class_labels, n_select)
	
	return [int(relevant_features[i]) for i in range(0,len(relevant_features))]

def learn_nbc(data_matrix, class_labels):
	nbc = naive_bayes.MultinomialNB()
	nbc.fit(data_matrix,class_labels)
	return nbc
	
	
def setup_foi(class_name, mapping_file, output_folder):
	# As of right now, we'll deal only with discrete data
	
	
	# with open(mapping_file, 'r') as fh:
		# header = fh.readline()
		# headers = header[1:-1].split("\t")
		
		# target_idx = headers.index(class_name)
		# samples_with_property = []
		# for row in fh:
			# columns = row[:-1].split("\t")
			# samples_with_property.append(columns[0])
				
	# # We want to keep the samples that are both in the sample_names AND in the samples_with_property
	# samples_in_features = env_features[:,0]
	# # for i in xrange(0,len(samples_in_features)):
	# # 	print('From biom: {} - From mapping: {} - are equal? {}'.format(sample_names[i],samples_in_features[i],sample_names[i]==samples_in_features[i]))
	# biom_idx, blah = find_common_samples(sample_names, samples_with_property)
	# feat_idx, blah = find_common_samples(samples_in_features, samples_with_property)
	
	# # new_sample_names = [sample_names[i] for i in biom_idx]
	# # new_matrix = biom_matrix[biom_idx,:]
	# # new_features = env_features[feat_idx,:]
	# return biom_matrix[biom_idx,:], [sample_names[i] for i in biom_idx], env_features[feat_idx,:]
	
	
	
	features_of_interest = {}
	features_of_interest[class_name] = {}
	features_of_interest[class_name]['discrete'] = 'D'
	features_of_interest[class_name]['type'] = 'T'
	features_dic = omdm.load_data_dic(mapping_file)
	
	# Have to deal with that better: 
	(desired_headers, desired_data) = omdm.load_desired_features(features_of_interest, features_dic, output_folder)
	env_features = np.array(desired_data)
	
	return env_features

	
def setup_biom(biom_file):
	ih = open(biom_file, 'rb')
	biom_table = parse_biom_table(ih)
	ih.close()

	matrix, site_names, variable_names = omdm.biom_table_to_array(biom_table)
	matrix = np.array(matrix).transpose()
	return matrix, site_names, variable_names
	

def select_submatrix(biom_matrix, sample_names, map_name, env_features, select_field, value_field):
	
	######
	## Still have to deal with when we want to have more than one selected value, e.g. L_Palmn AND R_Palm
	## This should actually be done at the same time as the feature reading
	######
	interesting_fields = value_field.split()
	
	with open(map_name, 'r') as fh:
		header = fh.readline()
		headers = header[1:-1].split("\t")
		
		target_idx = headers.index(select_field)
		samples_with_property = []
		for row in fh:
			columns = row[:-1].split("\t")
			if columns[target_idx] in interesting_fields:
			# if columns[target_idx] == value_field:
				samples_with_property.append(columns[0])
				
	# We want to keep the samples that are both in the sample_names AND in the samples_with_property
	samples_in_features = env_features[:,0]
	# for i in xrange(0,len(samples_in_features)):
	# 	print('From biom: {} - From mapping: {} - are equal? {}'.format(sample_names[i],samples_in_features[i],sample_names[i]==samples_in_features[i]))
	biom_idx, blah = find_common_samples(sample_names, samples_with_property)
	feat_idx, blah = find_common_samples(samples_in_features, samples_with_property)
	
	# new_sample_names = [sample_names[i] for i in biom_idx]
	# new_matrix = biom_matrix[biom_idx,:]
	# new_features = env_features[feat_idx,:]
	return biom_matrix[biom_idx,:], [sample_names[i] for i in biom_idx], env_features[feat_idx,:]


# For some reasons, python seems to reorder the list in a manner that is not interesting to us...
# It might come from the way we constructed everything - need to check that. 	
def align_matrices(label_vector, biom_site_names, label_site_names):
	ordered_labels = [0]*len(biom_site_names)
	blah = list(label_site_names)
	for i in xrange(0,len(biom_site_names)): # In theory, both site names vector should have the same size and same content
		# ordered_labels[i] = label_vector[list(label_site_names).index[biom_site_names)i]]]
		ordered_labels[i] = label_vector[blah.index(biom_site_names[i])]
		# ordered_labels[i] = label_vector[label_site_names.index[biom_site_names)i]]]

		
	# use numpy.where
	return ordered_labels
	
	
def main():
	parser = argparse.ArgumentParser(description = "")
	parser.add_argument("-b", "--biom-file", help="An input biom file", required=True)
	parser.add_argument("-m", "--mapping-file", help="A mapping file", required=True)
	parser.add_argument("-c", "--class-label", help="Which data are we trying to analyze", required=True)
	parser.add_argument("-o", "--output-folder", help="The folder to output our data to", required=True)
	parser.add_argument("-n", "--nb-features", help="Number of features to be selected", required=True)
	parser.add_argument("-p", "--predictor", help="Classifier/Predictor used") # As of today, contains only nbc
	parser.add_argument("-j", "--objective-function", help="Objective function for the feature selection algorithm", default="mim", required=False)
	parser.add_argument("-t", "--output-type", help="data output format. default: CSV options: csv, matlab, r, numpy", default="csv", required=False)
	parser.add_argument("-f", "--select-field", help="Field to extract a subset of the data. e.g. EN_BIOME, COUNTRY. The default considers the whole dataset", default="", required=False)
	parser.add_argument("-g", "--value-field", action="append", help="When used with -f specifies the value of the field to filter - THIS IS REQUIRED IF -f if present", default="", required=False)

## Need to be continued!!!!
	global output_type

	args = parser.parse_args()

	output_type = args.output_type.lower()
		
	# if our folder doesn't exist create it
	if not os.path.isdir(args.output_folder):
		os.makedirs(args.output_folder);

	env_features = setup_foi(args.class_label, args.mapping_file, args.output_folder)
	
	matrix, site_names, variable_names = setup_biom(args.biom_file)
	
	
	idx_from_biom, idx_from_mapping = find_common_samples(site_names, env_features[:,0])
	matrix = matrix[idx_from_biom,:]
	env_features = env_features[idx_from_mapping,:]
	if args.select_field is not None:
		matrix, site_names, env_features = select_submatrix(matrix, site_names, args.mapping_file, env_features, args.select_field, args.value_field)
	# So far, we have the biom file open and the environment parameters
	# We can now launch our feature selection algorithm
	further_param = [] # This has to be adapted to the case we are using other objective functions
	relevant_features = find_relevant_otus(matrix, np.array(map(int, env_features[:,1])), args.objective_function, int(args.nb_features), further_param)
	# this_nbc = learn_nbc(matrix[:,relevant_features], np.array(map(int, env_features[:,1])))
	save_found_otus(relevant_features, variable_names, args.output_folder)
	

		


if __name__ == "__main__":
	try:
		sys.exit(main())
	except:
		type, value, tb = sys.exc_info()
		traceback.print_exc()
		pdb.post_mortem(tb)
