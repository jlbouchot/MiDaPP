#!/usr/bin/env python 

import sys
from qiime.parse import parse_mapping_file_to_dict

__author__ = "Gregory Ditzler"
__copyright__ = "Copyright 2013, EESI Lab"
__credits__ = ["Gregory Ditzler"]
__license__ = "GPL"
__version__ = "0.1.0-dev"
__maintainer__ = "Gregory Ditzler"
__email__ = "gregory.ditzler@gmail.com"
__status__ = "Development"


def print_info(opts):
  """
  @opt: options from the command line. see --help
  """
  # import the the data using qiimes utilities. the o variable 
  # contains the dictionary with the map file contents. 
  try:
    o,c = parse_mapping_file_to_dict(open(opts.map))
  except IOError:  
    raise("File may not exist. ")

  # ---- print the metadata
  if opts.print_meta: 
    s = set()
    print "Meta-Data: "
    for k in o.keys():
      s = s.union(set(o[k].keys()))
    for r in s: 
      print "  " + r 

  # ---- print the sample IDs
  if opts.print_ids: 
    print "SampleIDs: "
    for k in o.keys():
      print "  " + k

  # ---- print a column from the map file
  if opts.print_col != None: 
    s = set()
    print "Column: " + opts.print_col + ": "
    for k in o.keys():
      s = s.union(set([o[k][opts.print_col]]))
    for r in s:
      print "  " + r
  return None 


if __name__ == "__main__":
  from optparse import OptionParser
  parser = OptionParser()
  parser.add_option("-a", "--map",
    dest = "map", 
    default = None, 
    help = "Map file."
  )
  parser.add_option("-b", "--print_meta", 
    dest = "print_meta", 
    action = "store_true", 
    help = "Print available metadata. "
  )
  parser.add_option("-c", "--print_ids",
    dest = "print_ids", 
    action = "store_true", 
    help = "Print sample IDs."
  )
  parser.add_option("-d", "--print_col", 
    dest = "print_col", 
    default = None, 
    help = "Print a column from the map file."
  )
  (options, args) = parser.parse_args()
  if options.map == None:
    print "Map file not specified."
    sys.exit(1)
  print_info(options)
  
