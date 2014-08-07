#!/usr/bin/env python 
import json 
import numpy
import scipy.sparse as sp
from optparse import OptionParser

__author__ = "Gregory Ditzler"
__copyright__ = "Copyright 2013, EESI Laboratory (Drexel University)"
__credits__ = ["Gregory Ditzler"]
__license__ = "GPL"
__version__ = "0.3.0"
__maintainer__ = "Gregory Ditzler"
__email__ = "gregory.ditzler@gmail.com"
__status__ = "development"

def build_parser():
  """
  build a command line option parser geared for running experiments with a 
  biom and map file with the NPFS-\beta algorithm. run the experiment script
  with either -h or --help for the latest options. 
  """
  parser = OptionParser()
  parser.add_option("-a", "--biom",
      dest = "biom",
      help = "Path to metahit biom.")
  parser.add_option("-b", "--map",
      dest = "map",
      help = "Path to metahit map file")
  parser.add_option("-c", "--n_select", 
      dest = "n_select",
      default = 50,
      type = int, 
      help = "Number of base algorithm to select.")
  parser.add_option("-d", "--method", 
      dest = "method",
      default = "MIM",
      help = "Base feature selection algorithm.")
  parser.add_option("-e", "--n_bootstraps",
      dest = "n_bootstraps",
      default = 100, 
      type = int,
      help = "Number of boot straps.")
  parser.add_option("-f", "--processes",
      dest = "processes",
      default = None, 
      type = int, 
      help = "Number of process to launch in the bootstraps.")
  parser.add_option("-g", "--bias",
      dest = "bias",
      default = 0.0, 
      type = float, 
      help = "Bias term for NPFS's null hypothesis.")
  parser.add_option("-i", "--alpha",
      dest = "alpha",
      default = 0.01,
      type = float, 
      help = "Size of the hypothesis test.")
  parser.add_option("-j", "--output",
      dest = "output",
      help = "Path to output.")
  parser.add_option("-k", "--caporaso",
      dest = "caporaso",
      default = "COMMON_SAMPLE_SITE",
      help = "Labels to use with Caporaso's data " \
          + "[BODY_SITE,BODY_HABITAT,SEX,COMMON_SAMPLE_SITE].")
  return parser


def load_biom(fname):
  """
  load a biom file and return a dense matrix 
  :fname - string containing the path to the biom file
  :data - numpy array containing the OTU matrix
  :samples - list containing the sample IDs (important for knowing 
    the labels in the data matrix)
  :features - list containing the feature names
  """
  o = json.loads(open(fname,"U").read())
  if o["matrix_type"] == "sparse":
    data = load_sparse(o)
  else:
    data = load_dense(o)

  samples = []
  for sid in o["columns"]:
    samples.append(sid["id"])
  features = []
  otu_phylo = {}
  for sid in o["rows"]:
    # check to see if the taxonomy is listed, this will generally lead to more 
    # descriptive names for the taxonomies. 
    if sid.has_key("metadata") and sid["metadata"] != None:
      if sid["metadata"].has_key("taxonomy"):
        #features.append(str( \
        #    sid["metadata"]["taxonomy"]).strip( \
        #    "[]").replace(",",";").replace("u'","").replace("'",""))
		features.append(sid["id"])
		otu_phylo[sid["id"]] = json.dumps(sid["metadata"]["taxonomy"])
      else:
        features.append(sid["id"])
        otu_phylo[sid["id"]] = 'Unknown'
    else:
      features.append(sid["id"])
      otu_phylo[sid["id"]] = 'Unknown'
  return data, samples, features, otu_phylo
  
  
  
def load_dense(obj):
  """
  load a biom file in dense format
  :obj - json dictionary from biom file
  :data - dense data matrix
  """
  n_feat,n_sample = obj["shape"]
  data = np.array(obj["data"])
  return data.transpose()

def load_sparse(obj):
  """
  load a biom file in sparse format
  :obj - json dictionary from biom file
  :data - dense data matrix
  """
  n_feat,n_sample = obj["shape"] 
  data = numpy.zeros((n_feat, n_sample))
  for val in obj["data"]:
    data[val[0], val[1]] = val[2]
  data = data.transpose() 
  return data

def load_map(fname):
  """
  load a map file. this function does not have any dependecies on qiime's
  tools. the returned object is a dictionary of dictionaries. the dictionary 
  is indexed by the sample_ID and there is an added field for the the 
  available meta-data. each element in the dictionary is a dictionary with 
  the keys of the meta-data. 
  :fname - string containing the map file path
  :meta_data - dictionary containin the mapping file information  
  """
  f = open(fname, "U")
  mfile = []
  for line in f: 
    mfile.append(line.replace("\n","").replace("#","").split("\t"))
  meta_data_header = mfile.pop(0)

  meta_data = {}
  for sample in mfile:
    sample_id = sample[0]
    meta_data[sample_id] = {}
    for identifier, value in map(None, meta_data_header, sample):
      meta_data[sample_id][identifier] = value 
  return meta_data

def index2feature(indices, features):
  """
  return the feature names corresponding to a set on indices
  :indices - list of integers
  :features - list of feature names
  :out - sub list of features indicated by indices 
  """
  return [features[index] for index in indices]

def write_output(fname, selected_features):
  """
  write a text file of feature names
  :fname - string containing the output file path 
  :selected_features - list of features to write
  """
  f = open(fname, "w")
  for feat in selected_features:
    f.write(feat + "\n")
  f.close()
  return None 

def discretize(labels_in):
  """
  convert a list on labels to a set of integers
  :labels_in - list on non-integer labels 
  :labels_out - list of integer labels
  :keys - dictionary translation from integer to key
  """
  keys = {}
  labels_out = []
  for label_set in enumerate(numpy.unique(labels_in)):
    keys[label_set[1]] = label_set[0]
  for lab in labels_in:
    labels_out.append(keys[lab])
  return numpy.array(labels_out), keys

