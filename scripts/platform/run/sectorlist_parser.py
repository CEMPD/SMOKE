#!/usr/bin/env python
####################################################################
##
##
## Author:
##      Alexis Zubrow,  UNC I.E. ,   January, 2009  initial version
####################################################################

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import csv, string, os


## function for parsing the SECTORLIST file
def parseSectorlist(fileName, sectorName, colName = "sectorcase", testColFlag = False):
    """
    Parses the sector list file and returns a value for a specific sector
    and column.

    Ignores additional white space and differences in capitalization.

    Inputs:
        fileName - sector list file name

        sectorName - the sector to parse out

        colName - the column to select, default is 'sectorcase'

        testColFlag - don't return value, test if column exist or has a real value

    Output:
       colVal - the value for that sector and col.  Typically this would
                be the merge approach to match against the merge
                dates file (ie. to pick appropriate col in merge dates).
                If testColFlag = true, this is True or False if column
                has a real value
    """
    
    ## open sector list file and read into a buffer
    f = open(os.path.expanduser(fileName), "rt")
    buffLst = f.readlines()
    f.close()
    
    ## parse the sector list file based on commas
    ## creates a reader Dictionary that can loop over
    readerDct = csv.DictReader(buffLst, skipinitialspace=True)

    firstTime = True
    
    ## loop over the reader dictionary, find the sector appropriate line
    for rowDct in readerDct:
        if firstTime:
            ## figure out what is the key for sectors
            keys = rowDct.keys()
            sectorKey = None
            colKey = None
            for key in keys:
                if key.upper().strip() == "SECTOR":
                    sectorKey = key
                if key.upper().strip() == colName.upper():
                    colKey = key
            if sectorKey is None:
                raise LookupError("Sector column not in sectorlist: %s" %fileName)
            if colKey is None:
                ## if only trying to test column, return false, not in sectorlist
                if testColFlag:
                    return False
                raise LookupError("Column name (%s) not in sectorlist: %s" \
                      %(colName, fileName))
            firstTime = False

        ## Check sector name, need to ignore case and spaces
        if rowDct[sectorKey].upper().strip() == sectorName.upper().strip():
            tmpVal = rowDct[colKey].strip()
            if testColFlag:
                ## test if the column has some value
                if tmpVal == "": return False
                else: return True
                
            ## return the col value, strip off excess white space
            return tmpVal
        
    ## if get to this point, sector is not in sectorlist:
    raise LookupError("Sector (%s) is not in sector list file: %s" \
          %(sectorName, fileName))


        
##-----------------------------------------------------------------------##
## Actual script


if __name__ == "__main__":
    from optparse import OptionParser

    ## Command line options
    usage = "usage: %prog [options] sectorlist"
    parser = OptionParser(usage=usage)

    parser.add_option("-s", metavar=" sector_name", dest="sectorName", default="",
                      help="Sector name.  Overrides environmental variable SECTOR")
    parser.add_option("-c", metavar=" column", dest="columnName", default="sectorcase",
                      help="Column name to parse out.  Default is 'sectorcase'.")
    parser.add_option("--test_col",  dest="testColFlag", default="False",
                      action="store_true",
                      help="Tests if column exists or has a valid value "
                      "(i.e not '').  Returns 'Y' or 'N'")
    (options, args) = parser.parse_args()

    ## get args
    if len(args) < 1:
        raise LookupError("Must pass SECTORLIST file")
    sectorlistFile = os.path.expanduser(args[0])
    
    ## get options
    sectorName = options.sectorName
    columnName = options.columnName
    testColFlag = options.testColFlag

    ## If any of the options are not set, try to get from the
    ## environmental variables
    if sectorName == "":
        try:
            sectorName = os.environ["SECTOR"]
        except:
            raise LookupError("Sector name must be set.  Not found in environment or passed as an option.")


    if testColFlag == True:
        colExists = parseSectorlist(sectorlistFile, sectorName, columnName, testColFlag)
        if colExists:
            print('Y')
        else:
            print('N')
    else:
        ## get the specific column from the SECTORLIST file
        colVal = parseSectorlist(sectorlistFile, sectorName, columnName)
        print(colVal)
