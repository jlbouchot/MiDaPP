# File created on 22 July 2013
import numpy
from biom.parse import parse_biom_table 
from qiime.parse import parse_mapping_file_to_dict
from sets import Set


__author__ = "Jean-Luc Bouchot"
__copyright__ = "Copyright 2011, The QIIME project"
__credits__ = ["Jean-Luc Bouchot"]
__license__ = "GPL"
__version__ = "1.7.0-dev"
__maintainer__ = "Jean-Luc Bouchot"
__email__ = "jean-luc.bouchot@drexel.edu"
__status__ = "Development"


def merge_matrices(data_matrix, env_table, site_names, site_names_env, variable_names, environmental_param):
	"""
		merge_matrices(data_matrix, env_table, site_names, site_names_env, variable_names, environmental_param)
		@data_matrix - the matrix containing the features; sites are ordered column wise, features are ordered row-wise
		@env_table - the matrix containing the features from the mapping file (same ordering as the data_matrix)
		@site_names - a list of names corresponding to each column of the data_matrix
		@site_names_env - the list of the site names from which the env_table was done
		@variable_names - list containing the names of the features for the data_matrix
		@environmental_param - a list of the different environmental parameters kept
		@complete_matrix (return) - 
		@complete_sites (return) - 
		@complete_variables (return) - 
	"""

	# Merge the two site name lists.
	complete_sites = list(Set(site_names_env).intersection(Set(site_names)))
	# Get the index sets useful for both matrices
	idx_env = []
	for a_name in complete_sites: # There must be a faster and neater way to do this, but I have no idea how (and don't want to look for)
		idx_env.append(site_names_env.index(a_name))
	idx_biom = []
	for a_name in complete_sites:
		idx_biom.append(site_names.index(a_name))

	# Extract AND ORDER submatrices with only the columns that are common to both of them 

	# Merge these two submatrices
	complete_matrix = numpy.vstack((data_matrix[:,idx_biom],env_table[:,idx_env])) # make sure the data type coincide
	return complete_matrix, complete_sites, variable_names.extend(environmental_param)


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
	data_matrix = [] # Maybe need to use a numpy.ndarray... Or we could have something like 
	# data_matrix = numpy.ndarray((biom_table['shape'])) but don't know how to use this.
	for one_observation in biom_table.iterObservationData(): # Observations correspond to what we call variables or features (is that true?)
		data_matrix.append(one_observation)

	return numpy.array(data_matrix), site_names, variable_names
		


def read_environment_table(map_fhandler, foi_fhandler):
	"""
		read_environment_table(map_fname, foi_fname)
		@map_fhandler - file handler to the mapping .txt tab delimited file
		@foi_fhandler - file handler to the tab delimited text file containing two columns: the first one contains the keys to the important environmental factors, and the second column contains whether the variables is continuous (e.g. temperature) or discrete (male/female)
		@data_matrix (return) - dense matrix containing both OTU/k-mer data AND environmental data
		@site_names (return) - names of the different samples
		@environmental_param (return) - names of the environmental parameters (as used in the foi_fname)
		@is_continuous (return) - vector containing whether a particular foi is continuous or not
		@hashtable_env (return) - correspondences between discrete variables and class affected (for instance male->1, female->0)
	"""
	obj, comm = parse_mapping_file_to_dict(map_fhandler)
	site_names = obj.rowKeys()
	environmental_param = []
	is_continuous = []
	for a_line in foi_fhandler:
		param_name, param_type = a_line.split("\t") # /!\ Should check that we actually have a single tab (for nice software engineering)
		environmental_param.append(param_name)
		is_continuous.append(param_type == "C") # Here again some better check should be done (tolerate 'continuous', 'c', 'C', undefined?) or at least throw an exception

	extracted_obj = obj.getCols(environmental_param)
	

	# Now for the discrete variables, we need to create the dictionaries/mapping from one value to a class assignment
	hashtable_env = {}
	for idx, value in enumerate(is_continuous):
		nb_discrete = 0
		if not value:
			# Get all the different potential values
			cur_col = []
			
			for a_sample in site_names:
				cur_col.append(extracted_obj[a_sample][environmental_param[idx]])



			# Single out these values
			unique_values = numpy.unique(cur_col)

			# Iterate over these single values to create the dictionary
			dict_hash = {}
			reverse_dict = {} # I'm pretty sure there is much better way to do that...
			for a_class, a_value in enumerate(unique_values):
				# Do some stuff
				dict_hash[a_class] = a_value # Not sure that's the best way to do it
				reverse_dict[a_value] = a_class


			# Convert the actual data in the data matrix to classes
			for a_sample in site_names:
				extracted_obj[a_sample][environmental_param[idx]] = reverse_dict[extracted_obj[a_sample][environmental_param[idx]]]
			
			
			nb_discrete = nb_discrete+1
			# Add this new map in the environment dictionary
			hashtable_env[environmental_param[idx]] = dict_hash

	# convert into a numpy array:
	data_matrix = numpy.array(extracted_obj.toLists())
	return data_matrix.transpose(), site_names, environmental_param, is_continuous, hashtable_env


def create_fused_data_matrix(biom_fname, map_fname, foi_fname):
	"""
		create_fused_data_matrix(biom_fname, map_fname=None, foi_fname=None)
		@biom_fname - path to the BIOM file containing the interesting counts (shall they be k-mers or OTUs, up to you!)
		@map_fname - contains the path to the mapping file with the environmental factors. default = None, in which case SOMETHING HAPPENS
		@foi_fname - path to the tab delimited text file containing two columns: the first one contains the keys to the important environmental factors, and the second column contains whether the variables is continuous (e.g. temperature) or discrete (male/female)
		@data_matrix (return) - dense matrix containing both OTU/k-mer data AND environmental data
		@site_names (return) - names of the different samples
		@variable_names (return) - Names of the different variables used
		@environmental_param (return) - names of the environmental parameters (as used in the foi_fname)
		@hashtable_env (return) - correspondences between discrete variables and class affected (for instance male->1, female->0)
	"""
	
	biom_table = parse_biom_table(open(biom_fname, 'U'))
	# Function to convert a biom dictionary to a numpy array
	data_matrix, site_names, variable_names = biom_table_to_array(biom_table)
	env_table, site_names_env, environmental_param, is_continuous, hashtable_env = read_environment_table(open(map_fname,'U'), open(foi_fname, 'U'))


	# Merge both data matrices making sure to keep only the sites that are common to both mapping and biom files and making sure that each site from one file is well aligned with the associated one from the other file.
	complete_matrix, complete_sites, complete_variables = merge_matrices(data_matrix, env_table, site_names, site_names_env, variable_names, environmental_param)

	return complete_matrix, complete_sites, variable_names, environmental_param, hashtable_env


def find_missing_data_locations(data_matrix):
	"""
		find_missing_data_locations(data_matrix)
		@data_matrix - the matrix organised such that data_matrix[r][c] corresponds to feature r in site c
		@row_coordinates (return) - an array containing the rows of the different None values in the data matrix
		@col_coordinates (return) - an array containing the cols of the different None values in the data matrix
		@is_none (return) - a numpy ndarray containing True's wherever a none has been found
	"""	
	is_none = numpy.equal(data_matrix,None)
	itemidx = numpy.where(is_none == True)
	return itemidx[0], itemidx[1], is_none


def complete_missing_data(data_matrix, location_list, method):
	return []

def get_learning_methods():
	# Return the different methods we want to implement
	return []
