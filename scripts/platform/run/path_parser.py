#!/usr/bin/env python

#########################################################################
##
## Takes the path to a file and parses out the directory
##
#########################################################################

from __future__ import print_function

from optparse import OptionParser
import os

## Command line options
usage = "usage: %prog [options] <file_path>"
parser = OptionParser(usage=usage)
parser.add_option("-p", dest="getParent", default=False,
                  action="store_true",
                  help="Get the parent, parent directory of the fileName. "
                  "For example if the file is /tmp/foo/test.txt this "
                  "will return /tmp.")

(options, args) = parser.parse_args()
getParent = options.getParent

if len(args) == 0:
    raise LookupError("SCRIPT ERROR: Need to pass filename with complete path to path_parser")

fileName = args[0]


## get canonical path
fileName = os.path.expanduser(fileName)

## get directory
dirName = os.path.dirname(fileName)

## if want parent direcotry
if getParent:
    dirName = os.path.dirname(dirName)
print(dirName)



