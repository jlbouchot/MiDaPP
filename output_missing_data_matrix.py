#!/usr/bin/python

import argparse
from biom.parse import parse_biom_table 
import numpy
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

def save_discrete_mapping(discrete_file, filetype, obs_dic): 
	
	if(filetype == 'csv'):
		with open(discrete_file, 'w') as fh:
			fh.write("#KEY\tMAPPED_VAL\tORIGINAL_VAL\n")
			for m in obs_dic.keys():
				for i in obs_dic[m].items():
					fh.write(m + "\t" + str(i[0]) + "\t" + str(i[1]) + "\n")

def create_discrete_mapping(features_arr, header_arr, features_of_interest):
	obs_dic = {}

	# skip the "sample"
	for f in range(1, len(header_arr)):
		feature = header_arr[f]
		print feature
		if(features_of_interest[feature]['discrete'] == "D"):
			i = 0
			obs_dic[feature] = {}
			for j in range(len(features_arr)):

				obs = features_arr[j][f]
				if(obs not in obs_dic[feature]):
					obs_dic[feature][obs] = i
					i = i + 1
	return obs_dic

def save_biom_matrix_csv(input_file, output_file):

	ih = open(input_file, 'r')
	biom_table = parse_biom_table(ih)

	oh = open(output_file, 'w')
	matrix, site_names, variable_names = biom_table_to_array(biom_table)

	oh.write("Samples\t")
	for x in xrange(len(matrix[0])):
		oh.write(site_names[x] + "\t")
	oh.write("\n")

	for x in xrange(len(matrix)):
		oh.write(variable_names[x] + "\t")
		for y in xrange(len(matrix[x])):
			oh.write(str(matrix[x][y]) + "\t")
		oh.write("\n")

def biom_table_to_array(biom_table):
	"""
		biom_table_to_array(biom_table)
		@biom_table - BIOM dictionary obtained using the biom parser
		@data_matrix (return) - numpy matrix containing the data: per row: a certain feature, per column: a certain site/sample
		@site_names (return) - name or key to the different sites
		@variable_names (return) - name of the feature used (for instance k-mers/OTUs counts)
	"""
	site_names = list(biom_table.SampleIds)
	variable_names = list(biom_table.ObservationIds)

	data_matrix = [o for o in biom_table.iterObservationData()] 

	return data_matrix, site_names, variable_names


def load_features_of_interest(foi_file):

	features = {}
	with open(foi_file, 'r') as fh:
		for i,row in enumerate(fh):
			# -1 to get rid of the newline
			split_row = row[:-1].split("\t")
			# Make sure we have 3 rows.
			if(len(split_row) != 3):
				raise("Three rows are required in your feature of interest file. Error row: '" + row[:-1] + "'")

			# populate our features of interest dict
			feature_name = split_row[0]
			discrete_val = split_row[1]
			type_val = split_row[2]

			features[feature_name] = {}
			features[feature_name]['discrete'] = discrete_val
			features[feature_name]['type'] = type_val 
	return features


def load_data_dic(mapping_file):
	features_dic = {}

	with open(mapping_file, 'r') as fh:
		header_line = fh.readline()	
		# use 1:-1 for the first # and the newline
		headers = header_line[1:-1].split("\t")

		for row in fh:
			columns = row[:-1].split("\t")
			sample_name =	columns[0]

			if sample_name not in features_dic:
				features_dic[sample_name] = {}

			for i, col in enumerate(columns):
				features_dic[sample_name][headers[i]] = col

	return features_dic


def save_desired_features(input_mapping_file, input_foi_file, output_data_matrix_file, output_discrete_file, filetype):
		desired_features_arr = []

		features_of_interest = load_features_of_interest(input_foi_file)
		features_dic = load_data_dic(input_mapping_file)

		# for each sample, add our desired data to an array
		for x, sample in list(enumerate(features_dic.keys())):
			desired_features_arr.append([])
			desired_features_arr[x].append(sample)

			# for each desired feature, add it from the features_dictionary
			for desired_feature in features_of_interest.keys():
				for sample_feature in features_dic[sample].keys():

					# if our feauture matches, parse it in
					if sample_feature == desired_feature:
						value = features_dic[sample][sample_feature]

						if value.lower() == "none" or value.lower() == 'nan':
							value = None # types.NoneType. ?
						elif features_of_interest[desired_feature]['type'].lower() == "f":
							value = float(value)
						elif features_of_interest[desired_feature]['type'].lower() == "i":
							value = int(value)
						elif features_of_interest[desired_feature]['type'].lower() == "t":
							value = str(value)
						else:
							print "Unknown Type: " + features['type'] + " for " + desired_feature[0]
					
						desired_features_arr[x].append(value)

		# create our header
		desired_features_header = []
		desired_features_header.append("Sample")
		for feature in features_of_interest:
			desired_features_header.append(feature)

		# create our discrete mappings 
		obs_dic = create_discrete_mapping(desired_features_arr, desired_features_header, features_of_interest)

		# convert requested variables into their discrete counterparts
		desired_features_arr = map(list, zip(*desired_features_arr))

		for f in range(1, len(desired_features_header)):
			feature = desired_features_header[f]
			if(features_of_interest[feature]['discrete'] == "D"):
				for o in range(len(desired_features_arr[f])):
					desired_features_arr[f][o] = obs_dic[feature][desired_features_arr[f][o]]

		desired_features_arr = map(list, zip(*desired_features_arr))
				
		# write our observation dictionary out
		save_discrete_mapping(output_discrete_file, filetype, obs_dic);

		
		with open(output_data_matrix_file, 'w') as fh:
			for header in desired_features_header[:-1]:
				fh.write(header + "\t")
			fh.write(desired_features_header[-1] + "\n")
				
			for row in desired_features_arr:
				for col in row[:-1]:
					fh.write(str(col) + "\t")
				fh.write(str(row[-1]) + "\n")

def main():
		parser = argparse.ArgumentParser(description = "")
		parser.add_argument("-i", "--biom-file", help="An input biom file", required=True)
		parser.add_argument("-m", "--mapping-file", help="A mapping file", required=True)
		parser.add_argument("-f", "--feature-of-interest-file", help="this file is \
		tab delimited with 3 columns, the first is the feature name, the second is \
		continuous/discrete and the third is the type: f - float, i - integer, t - \
		test", required=True)
		parser.add_argument("-o", "--output-folder", help="The folder to output our data to", required=True)
		parser.add_argument("-t", "--output-type", help="data output format. default: CSV options: csv, matlab, r, numpy", default="csv", required=False)

		args = parser.parse_args()

		# if our folder doesn't exist create it
		if not os.path.isdir(args.output_folder):
			os.mkdir(args.output_folder);

		save_desired_features(args.mapping_file,
													args.feature_of_interest_file,
													args.output_folder + "/desired_features.csv",
													args.output_folder + "/discrete_feature_mapping.csv",
													args.output_type)

		save_biom_matrix_csv(args.biom_file, args.output_folder + "/data_matrix.csv")


		# output our data matrix
		# if args.output_type == "csv":
		# 	numpy.savetxt(args.output_folder + "/biom_data_matrix.csv", data_matrix, delimiter="\t", fmt="%f")
		# elif args.output_type == "numpy":
		# 	numpy.save(args.output_folder + "/biom_data_matrix.npy", data_matrix, delimiter="\t")
		# elif args.output_type == "matlab":
		# 	import scipy.io
		# 	dic = {}
		# 	dic['matrix'] = data_matrix
		# 	scipy.io.savemat(args.output_folder + "/biom_data_matrix.mat", dic)
		# elif args.output_type == "r":
		#  	print "R is unsupported"	
		# else:
		# 	print "unknown output type"


if __name__ == "__main__":
	try:
		sys.exit(main())
	except:
		type, value, tb = sys.exc_info()
		traceback.print_exc()
		pdb.post_mortem(tb)
