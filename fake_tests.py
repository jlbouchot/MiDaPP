# Script to test the data_predictor module
import data_predictor


biom_fname = "biom_test_file.biom"
map_fname = "mapping_test_file.txt"
foi_fname = "foi_test_file.txt"

data_matrix, site_names, variable_names, environmental_param, hashtable_env = data_predictor.create_fused_data_matrix(biom_fname, map_fname, foi_fname)

print "Size of the data matrix after merging both environmental and metagenomics data: ", data_matrix.shape

# for a_param in variable_names:
# 	print a_param

print "List of the environmental parameter considered (should be similar to the one in the file foi_test_file.txt"
for an_env in environmental_param:
	print an_env


