#!/usr/bin/python

import argparse
from biom.parse import parse_biom_table 
import numpy
import os
import sys
import pdb

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
	missing_data_dic = {}

	with open(mapping_file, 'r') as fh:
		header_line = fh.readline()	
		# use 1:-1 for the first # and the newline
		headers = header_line[1:-1].split("\t")

		for row in fh:
			columns = row[:-1].split("\t")
			sample_name =	columns[0]

			if sample_name not in missing_data_dic:
				missing_data_dic[sample_name] = {}

			for i, col in enumerate(columns):
				missing_data_dic[sample_name][headers[i]] = col

	return missing_data_dic

def save_desired_features(input_mapping_file, input_foi_file, output_data_matrix_file, output_discrete_file):
		features = {}
		missing_data_dic = {}
		data_array = []
		obs_dic = {}

		features = load_features_of_interest(input_foi_file)
		missing_data_dic = load_data_dic(input_mapping_file)
		pdb.set_trace()


		# for each sample, add our desired data to an array
		data_array.append([])
		data_array[0].append("Sample")
		for feature in features:
			data_array[0].append(feature)

		for x, sample in list(enumerate(missing_data_dic.keys())):
			data_array.append([])
			data_array[x].append(sample)

			# for each desired feature, add it from the missing_data_dictionary
			for desired_feature in features.keys():
				for sample_feature in missing_data_dic[sample].keys():

					# if our feauture matches, parse it in
					if sample_feature == desired_feature:
						value = missing_data_dic[sample][sample_feature]

						if value.lower() == "none":
							value = None # types.NoneType. ?
						elif desired_feature[2].lower() == "f":
							value = float(value)
						elif desired_feature[2].lower() == "i":
							value = int(value)
					
						data_array[x].append(value)

		# transpose our output array
		data_array = map(list, zip(*data_array))

		# replace our value which are supposed to be discete with discrete values
		for f, feature in list(enumerate(data_array[1:])):
			if(features[feature[0]]['discrete'] == "discrete"):

				obs_dic[feature[0]] = {}
				i = 0
				for obs in feature[1:]:
					if(obs not in obs_dic[feature[0]]):
						obs_dic[feature[0]][obs] = i
						i = i + 1
				for x in xrange(1,len(feature)):
					if(feature[x] is not None):
						feature[x] = obs_dic[feature[0]][feature[x]]

		# write our observation dictionary out
		
		save_discrete_mapping(output_discrete_file, 'csv', obs_dic);

		data_array = map(list, zip(*data_array))
		with open(output_data_matrix_file, 'w') as fh:
			for row in data_array:
				for col in row[:-1]:
					fh.write(str(col) + "\t")
				fh.write(str(row[len(row) - 1]) + "\n")

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
													args.output_folder + "/discrete_feature_mapping.csv")

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
    sys.exit(main())
