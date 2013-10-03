from biom.parse import parse_biom_table 

import argparse
import data_predictor
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

def save_desired_features(mapping_file, foi_file, output_file):
		features = {}
		missing_data_dic = {}
		output_array = []
		output_fh = open(output_file, 'w');

		with open(foi_file, 'r') as fh:
			for i,row in enumerate(fh):
				split_row = row[:-1].split("\t")
				features[split_row[0]] = {}
				features[split_row[0]]['discrete'] = split_row[1]
				features[split_row[0]]['type'] = split_row[2]

		with open(mapping_file, 'r') as fh:
			header_line = fh.readline()	
			headers = header_line[1:-1].split("\t")

			for row in fh:
				columns = row.split("\t")
				sample_name =	columns[0]

				if sample_name not in missing_data_dic:
					missing_data_dic[sample_name] = {}

				for i, col in enumerate(columns):
						missing_data_dic[sample_name][headers[i]] = col 

		# for each sample, add our desired data to an array
		output_array.append([])
		output_array[0].append("Sample")
		for feature in features:
			output_array[0].append(feature)

		for x, sample in list(enumerate(missing_data_dic.keys()))[1:]:
			output_array.append([])
			output_array[x].append(sample)

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
					
						output_array[x].append(value)

		# transpose our output array
		output_array = map(list, zip(*output_array))

		# replace our value which are supposed to be discete with descrete values
		for feature in output_array[1:]: 
			if(features[feature[0]]['discrete'] == "discrete"):
				i = 0
				obs_dic = {}
				for obs in feature[1:]:
					if(obs not in obs_dic):
						obs_dic[obs] = i
						i = i + 1

				for x,obs in enumerate(feature[1:]):
					if feature[x] is not None:
						feature[x] = obs_dic[obs]

		pdb.set_trace()


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

		save_desired_features(args.mapping_file, args.feature_of_interest_file, args.output_folder + "/desired_features.csv");

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
