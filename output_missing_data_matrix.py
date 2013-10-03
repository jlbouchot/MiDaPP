
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

		data_matrix, site_names, variable_names, environmental_param, hashtable_env \
		= data_predictor.create_fused_data_matrix(args.biom_file, args.mapping_file, args.feature_of_interest_file)

		# if our folder doesn't exist create it
		if not os.path.isdir(args.output_folder):
			os.mkdir(args.output_folder);

		# output our data matrix
		if args.output_type == "csv":
			numpy.savetxt(args.output_folder + "/biom_data_matrix.csv", data_matrix, delimiter="\t", fmt="%f")
		elif args.output_type == "numpy":
			numpy.save(args.output_folder + "/biom_data_matrix.npy", data_matrix, delimiter="\t")
		elif args.output_type == "matlab":
			import scipy.io
			dic = {}
			dic['matrix'] = data_matrix
			scipy.io.savemat(args.output_folder + "/biom_data_matrix.mat", dic)
		else:
			print "unknown output type"

		pdb.set_trace()


if __name__ == "__main__":

    sys.exit(main())
