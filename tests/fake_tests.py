# Script to test the data_predictor module
import data_predictor
import numpy


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

# Test the artificial unknown mask generation on a random full matrix
missing_ratio = 0.2
from_rows = list(numpy.arange(len(variable_names)-len(environmental_param),len(variable_names)))
rdm_matrix = numpy.random.uniform(0,1,data_matrix.size).reshape(data_matrix.shape)
print rdm_matrix.shape
row_none, col_none, bin_mask, new_matrix = data_predictor.create_random_unknown_mask(data_matrix, from_rows, missing_ratio)

print "Number of missing data from the row list: ", len(row_none)
print "Number of missing data from the col list: ", len(col_none)

new_row_none, new_col_none, new_is_none = data_predictor.find_missing_data_locations(new_matrix)


if len(new_row_none) == len(row_none):
	print "Both have the same size" # We should actually make sure that have the same components!

i = 0
while i < 10: # could be tested for all using len(new_row_none)
	print bin_mask[row_none[i],col_none[i]], " should be true and ", new_matrix[row_none[i],col_none[i]], " should be none"
	i = i+1

