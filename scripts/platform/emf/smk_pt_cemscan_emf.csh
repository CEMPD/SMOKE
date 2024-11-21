#!/bin/tcsh -fx

# Version @(#)$Id$
# Path    $Source$
# Date    $Date$

# This script sets up needed environment variables for running SMOKE point
# sources using hourly emissions for only the CEMscan program.
#
# This script is intended to be used with the EMF
# source emissions in SMOKE for the EPA 2002 modeling platform, and 
# calls the scripts that runs the SMOKE programs. 
#
# Script created by : M. Houyoux, Environmental Protection Agency
# December, 2007
#
#*********************************************************************

## log w/ EMF server that script is running
$EMF_CLIENT -k $EMF_JOBKEY -s "Running" 

# set source category
setenv SMK_SOURCE P           # source category to process

## time-independent programs
setenv RUN_CEMSCAN   Y        #  run hourly CEM data scanner

setenv PROMPTFLAG         N       # Y (never set to Y for batch processing)
setenv AUTO_DELETE        Y       # Y deletes SMOKE I/O API output files (recommended)
setenv AUTO_DELETE_LOG    Y       # Y automatically deletes logs without asking
setenv DEBUGMODE          N       # Y changes script to use debugger
setenv DEBUG_EXE    totalview     # Sets the debugger to use when DEBUGMODE = Y

##############################################################################

## source the ASSIGN file
source $ASSIGNS_FILE

## List of all the helper scripts that are run in this script
set timetracker  = $SCRIPTS/run/timetracker_v2.csh
set smk_run      = $SCRIPTS/run/smk_run_v9.csh
set scriptname   = $SCRIPTS/emf/smk_pt_cemscan_emf.csh

## If TIMELOG_YN = N, don't do timetracker (default is Y)
if (! $?TIMELOG_YN) setenv TIMELOG_YN Y

# If timelog turned off, unset TIMELOG variable, which will turn timetracker off in the rest of the scripts
if ($TIMELOG_YN != N) setenv TIMELOG $LOGS/timelog_$namelabel.txt

## Set Time Log filename and initialize file
if ( $DEBUGMODE != Y && $DEBUGMODE != y && $TIMELOG_YN != N) then

   setenv TIMELOG $LOGS/timelog_${SECTOR}_${CASE}_cemscan.txt

   $EMF_CLIENT -k $EMF_JOBKEY -m "Running timetracker" -x $timetracker  ## log w/ EMF server
   $timetracker Y $TIMELOG
   if ( $status != 0 ) then
	echo "ERROR: running timetracker"
	$EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: running timetracker" -t "e"
	exit ( 1 )
   endif
endif

## Set up scripting environment variables prior to calling the Assigns file
setenv SUBSECT $SECTOR     # set variable for input/output names
setenv SRCABBR $SUBSECT    # set abbreviation for naming log files
setenv EISECTOR $SECTOR
setenv EMF_PERIOD $YEAR    # Set the EMF_PERIOD to the year


## Determine input CEM directory using EMISHOUR environment variable
set sc = ( `echo $EMISHOUR | grep -o "/" | wc -l` )
set hourpath = ( `echo $EMISHOUR | cut -d"/" -f1-$sc` )

## Set up input and output file names.
setenv FILELIST $hourpath/pthour.lst  # This is a CEMSCAN input file
setenv OUTFILE $CEMSUM
setenv REPFILE $CEMREPORT

##  Populate input list file, using all files in the directory
$EMF_CLIENT -k $EMF_JOBKEY -m "Creating list file" -x $scriptname -p $EMF_PERIOD   ## log w/ EMF server

/bin/ls $hourpath/$HOURLY_PREFIX* > $FILELIST
if ( $status == 1 ) then
    echo "SCRIPT ERROR: Cannot create list file - check directory permissions of"
    echo "              $hourpath"
    echo "              and check that CEM data files are present using prefix $HOURLY_PREFIX"
    $EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: creating list file for Cemscan, see stdout log" -t "e"
    exit ( 1 )
endif

## Run Cemscan
#
setenv RUN_PART0 Y

# Run programs
$EMF_CLIENT -k $EMF_JOBKEY -m "Running Cemscan" -x $smk_run -p $EMF_PERIOD   ## log w/ EMF server
source $smk_run 
if ( $status != 0 ) then
    echo "ERROR: Running Cemscan"
    $EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: Running Cemcsan" -t "e"
    exit ( 1 )
endif

setenv RUN_PART0 N
#
## Ending of script
#
exit( 0 )

 
