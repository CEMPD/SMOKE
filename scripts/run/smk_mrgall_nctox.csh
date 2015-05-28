#!/bin/csh -f
#
# Version @(#)$Id$
# Path    $Source$
# Date    $Date$
#
# This script sets up needed environment variables for merging
# emissions in SMOKE and calls the scripts that run the SMOKE programs. 
#
# Script created by : M. Houyoux, CEP Environmental Modeling Center 
#
#*********************************************************************

## Set optional customized SMKMERGE output file names
## setenv SMKMERGE_CUSTOM_OUTPUT  N  # Y define your own output file names from SMKMERGE

## Set Assigns file name
setenv ASSIGNS_FILE $SMKROOT/assigns/ASSIGNS.nctox.cmaq.cb05_soa.us12-nc

## Set source category
setenv SMK_SOURCE    E          # source category to process

## Set programs to run...

## For merging from previously generated gridded Smkmerge outputs
setenv RUN_MRGGRID   Y          # run gridded file merge program
#      NOTE: in sample script, Mrggrid used to create merged model-ready CMAQ files

## Program-specific controls...

## For Mrggrid
setenv MRG_DIFF_DAYS        N   # Y allows data from different days to be merged
#      MRGFILES         see "Script settings" below

## Script settings
setenv MRGFILES "AGTS_L NGTS_L BGTS_L PGTS_L RPDGTS_L RPVGTS_L RPPGTS_L RPHGTS_L" # logical file names to merge
setenv MRGGRID_MOLE         Y   # Y outputs mole-based file, musy be consistent with MRGFILES
setenv SRCABBR              abmp # abbreviation for naming log files
setenv PROMPTFLAG           N   # Y prompts for user input
setenv SMK_MAXERROR         100 # maximum number of error messages in log file
setenv SMK_MAXWARNING       100 # maximum number of warning messages in log file
setenv AUTO_DELETE          Y   # Y automatically deletes I/O API NetCDF output files
setenv AUTO_DELETE_LOG      Y   # Y automatically deletes log files
setenv DEBUGMODE            N   # Y runs program in debugger
setenv DEBUG_EXE            pgdbg # debugger to use when DEBUGMODE = Y
##############################################################################

## Loop through days to run Smkmerge and Mrggrid
#
setenv RUN_PART2 Y
setenv RUN_PART4 Y
set cnt = 0
set g_stdate_sav = $G_STDATE
while ( $cnt < $EPI_NDAY )

   @ cnt = $cnt + $NDAYS
   source $ASSIGNS_FILE   # Invoke Assigns file to set new dates

   if ( $MRGGRID_MOLE == Y ) then
      setenv OUTFILE $EGTS_L
   else 
      setenv OUTFILE $EGTS_S
   endif

   source smk_run.csh     # Run programs
   setenv G_STDATE_ADVANCE $cnt

end
setenv RUN_PART2 N
setenv RUN_PART4 N
unsetenv G_STDATE_ADVANCE

#
## Ending of script
#
exit( 0 )
