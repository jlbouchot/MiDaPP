#!/usr/bin/env python 

import sys
from qiime.parse import parse_mapping_file_to_dict

def print_info(opts):
  try:
    o,c = parse_mapping_file_to_dict(open(opts.map))
  except IOError:  
    raise("File may not exist. ")

  if opts.print_meta: 
    s = set()
    print "Meta-Data: "
    for k in o.keys():
      s = s.union(set(o[k].keys()))
    for r in s: 
      print "  " + r 
  if opts.print_ids: 
    print "SampleIDs: "
    for k in o.keys():
      print "  " + k

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

  (options, args) = parser.parse_args()
  if options.map == None:
    print "Map file not specified."
    sys.exit(1)
  print_info(options)
  
