# Script to test the data_predictor module
import data_predictor



biom_fname ="~/mount/study550/study_550_closed_reference_otu_table.biom"
map_fname = "~/mount/study550/study_550_mapping_file.txt"
foi_fname = "~/mount/study550/feature_of_interest_study_550.txt"

data_matrix, site_names, variable_names, environmental_param, hashtable_env = data_predictor.create_fused_data_matrix(biom_fname, map_fname, foi_fname)

print data_matrix.shape
print environmental_param

row_none, col_none, is_none = data_predictor.find_missing_data_locations(data_matrix)

print "There are ", len(list(row_none)), " none's in this study"

# Test the artificial unknown mask generation on a random full matrix
missing_ratio = 0.2
from_rows = xrange(len(variable_names),len(variable_names)+len(environmental_param))
rdm_matrix = numpy.random.uniform(0,1,data_matrix.size).reshape(data_matrix.shape)
row_none, col_none, bin_mask, new_matrix = data_predictor.create_random_unknown_mask(data_matrix, from_rows, missing_ratio)

print "Number of missing data from the row list: ", len(row_none)
print "Number of missing data from the col list: ", len(col_none)

new_row_none, new_col_none, new_is_none = data_predictor.find_missing_data_locations(new_matrix)

if len(new_row_none) == len(row_none):
	print "Both have the same size" # We should actually make sure that have the same components!
