#!/bin/csh -f
##
## Wrapper script to clean up any unzipped intermediary files.
## Tests if both an unzipped and zipped version exists before
## removing the unzipped version. Sector independent.  
##
## Script created by: A. Zubrow, UNC,  Jan 2009 
##
##*********************************************************************

set exitstat = 0  ## successful run 

echo "Cleaning unzipped intermediary files input file(s) for $ESDATE"

## Depending on source category [A|P|M] check for different files
if ( $MRG_SOURCE == A ) then
    echo "No cleanup of unzipped intermediary files for area" 
else if ( $MRG_SOURCE == M ) then
    echo "No cleanup of unzipped intermediary files for mobile"
else if ( $MRG_SOURCE == P ) then
    if ( -e $PLAY  && -e ${PLAY}.gz ) then
	echo "Removing unzipped PLAY file, $PLAY"
	rm $PLAY
    endif
else
      echo "SCRIPT ERROR: Unknown merge source $MRG_SOURCE"
      exit (1)
endif

## return exit status
exit($exitstat)
