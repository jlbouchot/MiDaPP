# Script to test the data_predictor module
import data_predictor


biom_fname = "study_550_closed_reference_otu_table.biom"
map_fname = "study_550_mapping_file.txt"
foi_fname = "feature_of_interest_study_550.txt"

data_matrix, site_names, variable_names, environmental_param, hashtable_env = data_predictor.create_fused_data_matrix(biom_fname, map_fname, foi_fname)

print data_matrix.shape
print environmental_param

row_none, col_none, is_none = data_predictor.find_missing_data_locations(data_matrix)

print "There are ", len(list(row_none)), " none's in this study"

