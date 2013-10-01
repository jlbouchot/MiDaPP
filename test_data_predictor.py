# Script to test the data_predictor module
import data_predictor
import numpy


biom_fname = "../../study_722_closed_reference_otu_table.biom"
map_fname = "study_722_mapping_file.txt"
foi_fname = "feature_of_interest_study_722.txt"

data_matrix, site_names, variable_names, environmental_param, hashtable_env = data_predictor.create_fused_data_matrix(biom_fname, map_fname, foi_fname)

print data_matrix.shape
print environmental_param

row_none, col_none, is_none = data_predictor.find_missing_data_locations(data_matrix)

print "There are ", len(list(row_none)), " none's in this study"
# print row_none
# print col_none

print "None's for the variable:", variable_names[19233], "are : ", data_matrix[19233,:]

# Add a sorting: number of missing data per site then number of missing data 


## Test the artificial unknown mask generation on a random full matrix
#missing_ratio = 0.2
#from_rows = xrange(len(variable_names)-len(environmental_param),len(variable_names))

#rdm_matrix = numpy.random.uniform(0,1,data_matrix.size).reshape(data_matrix.shape)
#row_none, col_none, bin_mask, new_matrix = data_predictor.create_random_unknown_mask(data_matrix, from_rows, missing_ratio)

#print "Number of missing data from the row list: ", len(row_none)
#print "Number of missing data from the col list: ", len(col_none)

#new_row_none, new_col_none, new_is_none = data_predictor.find_missing_data_locations(new_matrix)
#print len(new_row_none)
#if len(new_row_none) == len(row_none):
#	print "Both have the same size" # We should actually make sure that have the same components!
