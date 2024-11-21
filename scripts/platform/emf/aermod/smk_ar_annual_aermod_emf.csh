#!/bin/tcsh -f

# Version @(#)$Id$
# Path    $Source$
# Date    $Date$

# This script sets up needed environment variables for running SMOKE area
# sources using annual emissions for 1 or more months.
#
# This script is intended to be used with the EMF
# source emissions in SMOKE for the EPA 2002 modeling platform, and 
# calls the scripts that runs the SMOKE programs. 
#
# Script created by : M. Houyoux, Environmental Protection Agency
# July, 2007
# Modified to work w/ EMF: A. Zubrow, UNC - IE,  August, 2007
# Modified to allow county/scc/tprof and cell/county reports: R. Mason,
#     EPA, September 2008
# Modified to support concatenation for GSCNV and GSPRO_COMBO, M. Houyoux Dec, 2008
#
#*********************************************************************

## log w/ EMF server that script is running
$EMF_CLIENT -k $EMF_JOBKEY -s "Running" 

setenv SMK_SOURCE A           # source category to process
setenv MRG_SOURCE $SMK_SOURCE # source category to merge

## month-dependent programs
setenv RUN_SMKINVEN  Y        #  run inventory import program
setenv RUN_SPCMAT    N        #  run speciation matrix program
setenv RUN_GRDMAT    Y        #  run gridding matrix program
setenv RUN_TEMPORAL  Y 

setenv PROMPTFLAG         N       # Y (never set to Y for batch processing)
setenv AUTO_DELETE        Y       # Y deletes SMOKE I/O API output files (recommended)
setenv AUTO_DELETE_LOG    Y       # Y automatically deletes logs without asking
setenv DEBUGMODE          N 
##############################################################################

switch ( $#argv )
   case 0:
   case 1:
   case 2:
      echo "SCRIPT ERROR: Script requires arguments for a grid name"
      echo " "
      echo "  This script expects to be called using one of the following argument lists:"
      echo "     <grid abbrv> <I/O API gridname>"
      echo " "
      echo "  You can either use one approach or the other (differing by the -m or -q options)."
      echo " "
      echo "  In the above list, the arguments are defined as follows:"
      echo "     <grid abbrv>       : Grid abbreviation (e.g., 36US1)"
      echo "     <I/O API gridname> : I/O API gridname that needs to match entry in the"
      echo "                          GRIDDESC input file"
      echo " "
      echo "  Examples:"
      echo "     <script name> 36US1 36US1_148X112"
      echo "              This example runs the script"
      echo "              for the 36US1 grid, with no spinup days and"
      echo "              gives a label to the TIMELOG file of jan-sep." 
      $EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: smoke script did not receive more than 3 arguments" -t "e"
      exit( 1 )
endsw

# Get the first two options for the grid abbreviation and I/O API grid
setenv GRID "$argv[1]"
setenv IOAPI_GRIDNAME_1 "$argv[2]"

## source the ASSIGNS file
source $ASSIGNS_FILE

## List of all the helper scripts that are run in this script
set emf_cleanup  = $SCRIPTS/run/emf_cleanup.csh
set set_months   = $SCRIPTS/run/set_months_v4.csh
set timetracker  = $SCRIPTS/run/timetracker_v2.csh
set combine_data = $SCRIPTS/run/combine_data_v6.csh
set smk_run      = $SCRIPTS/run/smk_run_v9.csh
set qa_run       = $SCRIPTS/run/qa_run_v10.csh
set m3stat       = $SCRIPTS/run/m3stat_chk_v6.csh
set set_days     = $SCRIPTS/run/set_days_v5.csh
set log_analyzer = $SCRIPTS/log_analyzer/log_analyzer.py
set msg_list     = $SCRIPTS/log_analyzer/known_messages.txt
set duplic_chk   = $SCRIPTS/run/duplicate_check.csh
set path_parser  = $SCRIPTS/run/path_parser.py
set aermod_np    = $SCRIPTS/aermod/nonpt.pl
set aermod_og    = $SCRIPTS/aermod/np_oilgas.pl
set aermod_rwc   = $SCRIPTS/aermod/rwc.pl

## If running from EMF, move old EMF-created scripts to "old"
if ( $?EMF_JOBID ) then
   source $emf_cleanup
   if ( $status != 0 ) then
	echo "ERROR: running EMF script/log cleanup script"
	$EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: running EMF script/log cleanup script" -t "e" -x $emf_cleanup
	exit( 1 )
   endif
endif

## Set naming label
set namelabel = ${SECTOR}_${CASE}_${GRID}
if ( $?TLABEL ) then
  set namelabel = ${namelabel}_$TLABEL
endif

## Record the helper scripts being used
set suffix = _$namelabel.txt
echo "# Helper scripts used for $SECTOR" > $LOGS/helper_scripts_list$suffix
echo $emf_cleanup >> $LOGS/helper_scripts_list$suffix
echo $set_months >> $LOGS/helper_scripts_list$suffix
echo $timetracker >> $LOGS/helper_scripts_list$suffix
echo $combine_data >> $LOGS/helper_scripts_list$suffix
echo $smk_run >> $LOGS/helper_scripts_list$suffix
echo $qa_run >> $LOGS/helper_scripts_list$suffix
echo $m3stat >> $LOGS/helper_scripts_list$suffix
echo $set_days >> $LOGS/helper_scripts_list$suffix
echo $log_analyzer >> $LOGS/helper_scripts_list$suffix
echo $msg_list >> $LOGS/helper_scripts_list$suffix
echo $duplic_chk >> $LOGS/helper_scripts_list$suffix
echo $path_parser >> $LOGS/helper_scripts_list$suffix

## Set Time Log filename and initialize file
setenv TIMELOG $LOGS/timelog_$namelabel.txt

# Only initialize TIMELOG if it doesn't already exist, since the timeracker 
#   can now delete/add entries to prevent duplicates
if ( ! -e $TIMELOG ) then 
   $EMF_CLIENT -k $EMF_JOBKEY -m "Initializing Time Log" -x $timetracker  ## log w/ EMF server
   $timetracker Y $TIMELOG
   if ( $status != 0 ) then
	echo "ERROR: running timetracker"
	$EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: running timetracker to initialize time log" -t "e" -x $timetracker
	exit ( 1 )
   endif
endif

## Set up scripting environment variables prior to calling the Assigns file
setenv SUBSECT $SECTOR                   # set variable for input/output names
setenv SRCABBR $SUBSECT                  # set abbreviation for naming log files

## Run Smkinven, Grdmat, and Spcmat

setenv RUN_PART1 Y
source $ASSIGNS_FILE               # Invoke Assigns file

## Construct inventory list (NOTE: For other sectors, the last argument
#      on the combine_data.csh script will need to include the month
#      number) Needs to have a unique env PREFIX that is not used by
#      intermediary or output files.
$combine_data EMISINV $LISTFILE list
if ( $status != 0 ) then
    echo "ERROR: Could not run combine_data.csh to create area file list:"
    echo "       $LISTFILE"
    $EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: Could not run combine_data to create area file list" -t "e" -x $combine_data
    exit ( 1 )
endif
setenv ARINV $LISTFILE  #  Set the area inventory list file to the generic list file

## If DAY_SPECIFIC_YN = Y, then set up list file for daily FF10
if ( $?DAY_SPECIFIC_YN ) then
   if ( $DAY_SPECIFIC_YN == Y ) then
      $combine_data EMISDAY $DAYLISTFILE list
      if ( $status != 0 ) then
          echo "ERROR: Could not run combine_data.csh to create area daily inventory file list:"
          echo "       $DAYLISTFILE"
          $EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: Could not run combine_data to create area file list" -t "e" -x $combine_data
          exit ( 1 )
      endif
      setenv ARDAY $DAYLISTFILE  #  For FF10 daily nonpoint. Not sure if this is correct? (C. Allen - 15 Dec 2011)
   endif
endif

## Check speciation, gridding, and temporal cross-reference files for duplicates

echo "SCRIPT NOTE: Scanning AGREF for duplicate records"
$duplic_chk $AGREF G
if ( $status != 0 ) then
    echo "ERROR: Duplicate records detected in AGREF by duplic_chk.csh, or other script problem"
    $EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: Duplicate records detected in AGREF by duplic_chk.csh, or other script problem" -t "e" -x $duplic_chk
    exit ( 1 )
endif

echo "SCRIPT NOTE: Scanning ATREF for duplicate records"
$duplic_chk $ATREF T
if ( $status != 0 ) then
    echo "ERROR: Duplicate records detected in ATREF by duplic_chk.csh, or other script problem"
    $EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: Duplicate records detected in ATREF by duplic_chk.csh, or other script problem" -t "e" -x $duplic_chk
    exit ( 1 )
endif

## Set the EMF_PERIOD to the year
setenv EMF_PERIOD $YEAR

# Call EMF Client for one-time steps 
$EMF_CLIENT -k $EMF_JOBKEY -m "Running SMOKE steps for one-time steps" -p $EMF_PERIOD   ## log w/ EMF server

# Run programs for "part 1"
source $smk_run 

# Check status of QA run to see if it worked. Give error if failed
if ( $status != 0 ) then
    echo "ERROR: Running smk_run for one-time steps"
    $EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: Running smk_run for one-time steps" -t "e" -x $smk_run -p $EMF_PERIOD
    set exitstat = 1
    goto end_of_script
endif

setenv RUN_PART1 N

## Run programs for "part 2" 

setenv RUN_PART2 Y
setenv MONTH jan           # set variable for month name

## Determine dates to run in this month
setenv MONTH_ARRAY  1     # MONTH_ARRAY can have as many months listed as needed

# Source assigns file BEFORE set_days_v2.csh to set PROCDATES and SMK_RUN_DATES
source $ASSIGNS_FILE                  # Invoke Assigns file to set new dates
setenv SMK_RUN_DATES $SMK_RUN_DATES_1 

## Set the EMF_PERIOD for this month and year
setenv EMF_PERIOD "${MONTH}_${YEAR}"

# Call EMF Client for current period 
$EMF_CLIENT -k $EMF_JOBKEY -m "Running SMOKE steps for month $MONTH" -p $EMF_PERIOD   ## log w/ EMF server

setenv T_TYPE $L_TYPE                 # Set temporal type to type for temporal

# Source assigns file BEFORE set_days_v3.csh to set PROCDATES and SMK_RUN_DATES settings
source $ASSIGNS_FILE                  # Invoke Assigns file to set new dates
setenv SMK_RUN_DATES $SMK_RUN_DATES_1 

setenv SPINUP_ARRAY 0

source $set_days   # Call script to set dates for run
if ( $status != 0 ) then
    echo "ERROR: Running set_days"
    $EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: Running set_days for $MONTH" -t "e" -x $set_days -p $EMF_PERIOD
    exit (1)
endif

setenv RUN_TEMPORAL Y 

# Source assigns file AFTER set_days_v3.csh so that it can use PROCDATES file
source $ASSIGNS_FILE                  # Invoke Assigns file to set new dates

## Run Temporal (using PROCDATES file - run for all days needed in month)
source $smk_run       # Run programs
if ( $status != 0 ) then
    echo "ERROR: Running smk_run for part 2 in $EMF_PERIOD" 
    $EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: Running smk_run for part 2 in $EMF_PERIOD" -t "e" -x $smk_run -p $EMF_PERIOD
    set exitstat = 1
    goto end_of_script
endif

## Generate SMOKE reports for AERMOD
setenv QA_TYPE aermod 
setenv QA_LABEL $SUBSECT           # Used to name the report inputs and outputs
setenv REPLABEL $SUBSECT           # Used internally by Smkreport
setenv RUN_SMKREPORT Y

## Advance date to first representative day so that correct ATSUP is referenced
set g_stdate_sav = $G_STDATE
set line = `head -1 $SMK_RUN_DATES`
@ diff = $line[1] - $g_stdate_sav
setenv G_STDATE_ADVANCE $diff
source $ASSIGNS_FILE   # Invoke Assigns file to set new dates

source $qa_run
# Check status of QA run to see if it worked. Give error if failed
if ( $status != 0 ) then
   echo "ERROR: Running qa_run for $QA_TYPE" 
   $EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: Running qa_run for $QA_TYPE" -t "e" -x $qa_run -p $EMF_PERIOD 
   set exitstat = 1
   goto end_of_script
endif  # If qa script failure or not

# Set path of REPORT for AERMOD helper file scripts
# C. Allen (7 Apr 2017): If Smkreport is turned off via RUNSET (i.e. if the user only wants to run the post-processing step),
#   then REPORT25 is never defined.
if ($?REPORT25) then
   setenv REPORT $REPORT25
else
   setenv REPORT $REPOUT/$QA_TYPE/rep_${QA_LABEL}_${CASE}_${ESDATE}_${QA_TYPE}_${GRID}.txt
endif

unsetenv LOGFILE

setenv RUN_PART2 N

setenv OUTPUT_DIR $PROJECT_ROOT/$CASE/aermod
if ( $?OUTPUT_DIR ) then
    mkdir -p $OUTPUT_DIR
    chmod ug+rwx $OUTPUT_DIR
    foreach subdir (locations emis temporal parameters xwalk)
        mkdir -p $OUTPUT_DIR/$subdir
        chmod ug+rwx $OUTPUT_DIR/$subdir
    end
endif

setenv PERL5LIB $SCRIPTS/aermod
if ($SECTOR == "np_oilgas") then
    echo "SCRIPT NOTE: Generating np_oilgas style AERMOD helper files"
    $EMF_CLIENT -k $EMF_JOBKEY -m "Generating np_oilgas style AERMOD helper files"   ## log w/ EMF server
    $aermod_og
else if ($?ATPRO_DAILY || $SECTOR == "ag") then
    echo "SCRIPT NOTE: Generating RWC style AERMOD helper files"
    $EMF_CLIENT -k $EMF_JOBKEY -m "Generating RWC style AERMOD helper files"   ## log w/ EMF server
    $aermod_rwc
else
## Generate AERMOD helper files
    echo "SCRIPT NOTE: Generating nonpoint style AERMOD helper files"
    $EMF_CLIENT -k $EMF_JOBKEY -m "Generating nonpoint style AERMOD helper files"   ## log w/ EMF server
    $aermod_np
endif

# Label for the end of the script, used during script abort
end_of_script:

## Register time log
echo "SCRIPT NOTE: Registering time log"
$EMF_CLIENT -k $EMF_JOBKEY -F $TIMELOG -T "SMOKE time log (External)" -N "SMOKE timelog $namelabel" -O "Timelog $namelabel (External)"

## Run log file analyzer
$log_analyzer -k $msg_list --list_unknowns -l 3 -f $REPOUT/log_analyzer/rep_logs_${namelabel}_level3.csv -d $LOGS
if ( $status != 0 ) then
   echo "ERROR: running log_analyzer, level 3"
   $EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: running log_analyzer, level3" -t "e" -x $log_analyzer
   set exitstat = 1
endif

## If log analyzer returns exit status 10, then the log analyzer ran to completion, but found some errors or unrecognized
#  warnings. In this case, we want this script to exit with "Failed" status, but we still want to register the level 1 report in the EMF.
#  If log analyzer returns exit status 1, it didn't work at all, so bomb out immediately.

$log_analyzer -k $msg_list --list_unknowns -e 1 -f $REPOUT/log_analyzer/rep_logs_${namelabel}_level1.csv -d $LOGS
if ( $status == 10 ) then
   echo "SCRIPT NOTE: Registering log summary, level1"
   $EMF_CLIENT -k $EMF_JOBKEY -F $REPOUT/log_analyzer/rep_logs_${namelabel}_level1.csv \
        -T "Log summary level 1" -O "Level 1 log summary ${namelabel}"
	
   echo "ERROR: Level 1 errors or Level 0 warnings found. Check log files, fix inputs, and rerun."
   $EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: Level 1 errors or Level 0 warnings found. Check log files, fix inputs, and rerun." -t "e" -x $log_analyzer
   set exitstat = 1

## Register log analyzer output
else if ( $status != 0 ) then
   echo "ERROR: running log_analyzer, level 1"
   $EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: running log_analyzer, level 1" -t "e" -x $log_analyzer
   set exitstat = 1
else
   echo "SCRIPT NOTE: Registering log summary, level1"
   $EMF_CLIENT -k $EMF_JOBKEY -F $REPOUT/log_analyzer/rep_logs_${namelabel}_level1.csv \
        -T "Log summary level 1" -O "Level 1 log summary ${namelabel}"
endif

## Ending of script
#
exit( $exitstat )
