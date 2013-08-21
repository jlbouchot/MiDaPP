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
	idx_env = [site_names_env.index(a_name) for a_name in complete_sites]
	idx_biom = [site_names.index(a_name) for a_name in complete_sites]
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
	data_matrix = [o for o in biom_table.iterObservationData()] 
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
	environmental_param_ = []
	is_continuous = []
	data_types = []
	for a_line in foi_fhandler:
		param_name, param_type, data_type = a_line.split("\t") # /!\ Should check that we actually have two tabs (for nice software engineering)
		environmental_param_.append(param_name)
		is_continuous.append(param_type.lower() == "c") # Here again some better check should be done (tolerate 'continuous', 'c', 'C', undefined?) or at least throw an exception
		data_types.append(data_type)

	extracted_obj = obj.getCols(environmental_param_)
	environmental_param = extracted_obj.colKeys()

	# Convert from Strings to whatever it should be! -> Should be able to do better
	for idx, a_variable in enumerate(environmental_param_):

		for a_site in site_names:
			if extracted_obj[a_site][a_variable].lower() == "none":
				extracted_obj[a_site][a_variable] = None # types.NoneType. ?
			elif data_types[idx].lower() == "f":
				extracted_obj[a_site][a_variable] = float(extracted_obj[a_site][a_variable])
			elif data_types[idx].lower() == "i":
				extracted_obj[a_site][a_variable] = int(extracted_obj[a_site][a_variable])


	# Now for the discrete variables, we need to create the dictionaries/mapping from one value to a class assignment
	hashtable_env = {}
	for idx, value in enumerate(is_continuous):

		if not value:
			# Get all the different potential values
			cur_col = []
			
			for a_sample in site_names:
				cur_col.append(extracted_obj[a_sample][environmental_param_[idx]])



			# Single out these values
			unique_values = numpy.unique(cur_col)

			# Iterate over these single values to create the dictionary
			dict_hash = {}
			reverse_dict = {} # I'm pretty sure there is much better way to do that...
			for a_class, a_value in enumerate(unique_values):
				# Do some stuff
				if a_value is not None:
					dict_hash[a_class] = a_value # Not sure that's the best way to do it
					reverse_dict[a_value] = a_class


			# Convert the actual data in the data matrix to classes
			for a_sample in site_names:
				if extracted_obj[a_sample][environmental_param_[idx]] is not None:
					extracted_obj[a_sample][environmental_param_[idx]] = reverse_dict[extracted_obj[a_sample][environmental_param_[idx]]]
			
			# Add this new map in the environment dictionary
			hashtable_env[environmental_param_[idx]] = dict_hash


	# convert into a numpy array:
	data_matrix = numpy.array(extracted_obj.toLists())

	return data_matrix.transpose(), site_names, environmental_param, is_continuous, hashtable_env


def create_fused_data_matrix(biom_fname, map_fname, foi_fname, is_fasta=False):
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
	if not is_fasta:
		biom_table = parse_biom_table(open(biom_fname, 'U'))
		# Function to convert a biom dictionary to a numpy array
		data_matrix, site_names, variable_names = biom_table_to_array(biom_table)
		env_table, site_names_env, environmental_param, is_continuous, hashtable_env = read_environment_table(open(map_fname,'U'), open(foi_fname, 'U'))
	else: # We read the k-mers from the fasta file instead of the biom file
		# Call the count-k-mer function
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
	is_none = numpy.equal(data_matrix,None) # & numpy.equal(data_matrix,none)
	is_nan = numpy.isnan(data_matrix)
	is_none = numpy.logical_or(is_none, is_nan)
	itemidx = numpy.where(is_none == True)
	return itemidx[0], itemidx[1], is_none


def complete_missing_data(data_matrix, row_location, col_location, method):
	"""
		complete_missing_data(data_matrix, row_location, col_location, method)
		@data_matrix - the matrix organised such that data_matrix[r][c] corresponds to feature r in site c
		@row_location - the list of rows where to find the missing data
		@col_location - the list of col number where the missing data are (i.e. data_matrix[row_location[j][col_locationj]] is missing)
		@method - a string containing the method for completion (should be one of the output of get_learning_method)
		@completed_matrix (return) - a numpy matrix filled witht the estimation based on the data_matrix
		@completed_value (return) - a vector containing the values that have been estimated at the given locations
	"""	
	return []

def get_learning_methods():
	# Return the different methods we want to implement
	return []

def create_random_unknown_mask(data_matrix, from_rows, missing_ratio):
	"""
		create_random_unknown_,ask(data_matrix, from_rows, missing_ratio)
		@data_matrix - the matrix organised such that data_matrix[r][c] corresponds to feature r in site c - It should correspond to the "complete" matrix (i.e. with both -omics information and environment)
		@from_rows  - subset of rows from which the missing data should be created (the idea of that, is that, when dealin with k-mer counts, we can assume these counts to be complete -> No missing data in the first rows)
		@missing_ratio - the ration of number of unknown data by the number of total data 
		@row_none (return) - a list of row numbers where to find a (artificial) missing measure
		@col_none (return) - the same list as above but containing the column numbers
		@unknown_mask (return) - a numpy ndarray containing True if at the location of unknown data
		@new_matrix (return) - a copy of the data matrix with nones at the random locations given by the unkwnon mask
	"""	
	# Allocate memoray for a matrix
	unknown_mask = numpy.empty_like(data_matrix, bool) # Creates a matrix of the size of the data matrix filled with false
	unknown_mask.fill(False)
	r,c = data_matrix.shape
	nb_interesting_rows = len(from_rows)

	# Create random submatrix
	unknown_submask = numpy.random.uniform(0, 1, nb_interesting_rows*c).reshape(nb_interesting_rows,c) <= missing_ratio
	
	# Copy original data_matrix
	new_matrix = data_matrix.copy()

	row_none = []
	col_none = []

	# Go over all the rows of interest and replace the unknwon mask by these values and create the row_non and col_none vectors
	for row_idx, one_row in enumerate(from_rows):
		# First replace the row of the larger mask with the corresponding one in the smaller 		
		unknown_mask[one_row,:] = unknown_submask[row_idx,:]

		# Now create/fill the row_none, col_none vectors
		for one_col in xrange(0,c-1):
			if unknown_mask[one_row,one_col]:
				row_none.append(one_row)
				col_none.append(one_col)
				new_matrix[one_row,one_col] = None

	return row_none, col_none, unknown_mask, new_matrix
