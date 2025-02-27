#!/bin/csh -f

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
# 10/5/2007
# Modified to support concatenation for GSCNV and GSPRO_COMBO, M. Houyoux Dec, 2008
# Modified: Added SPCMAT_PERIOD for GSPRO_COMBO, A Zubrow Mar, 2011
#
#*********************************************************************

## log w/ EMF server that script is running
$EMF_CLIENT -k $EMF_JOBKEY -s "Running" 

# set source category
setenv SMK_SOURCE A           # source category to process
setenv MRG_SOURCE $SMK_SOURCE # source category to merge

setenv PROMPTFLAG         N       # Y (never set to Y for batch processing)
setenv AUTO_DELETE        Y       # Y deletes SMOKE I/O API output files (recommended)
setenv AUTO_DELETE_LOG    Y       # Y automatically deletes logs without asking
if ( ! $?DEBUGMODE ) then
   setenv DEBUGMODE          N       # Y changes script to use debugger
endif
if ( ! $?DEBUG_EXE ) then
   setenv DEBUG_EXE     idb    # Sets the debugger to use when DEBUGMODE = Y
endif

##############################################################################

switch ( $#argv )
   case 0:
   case 1:
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
      echo "  Examples:"
      echo "     <script name> 36US1 36US1_148X112"
      echo "              This example runs the script "
      echo "              for the 36US1 grid, with no spinup days and"
      echo "              gives a label to the TIMELOG file of jan-sep." 
      $EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: smoke script did not receive more than 2 arguments" -t "e"
      exit( 1 )
endsw

# Get the first two options for the grid abbreviation and I/O API grid
setenv GRID "$argv[1]"
setenv IOAPI_GRIDNAME_1 "$argv[2]"

setenv RUN_SMKINVEN Y 
setenv RUN_SPCMAT   N 
setenv RUN_GRDMAT   Y 
setenv RUN_TEMPORAL Y

## source the ASSIGN file
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
set aermod_nr   = $SCRIPTS/aermod/nonroad.pl

## Invoke script to interpret calling arguments
$EMF_CLIENT -k $EMF_JOBKEY -m "Running set_months" -x $set_months  ## log w/ EMF server
set exitstat = 0
source $set_months -m "1 2 3 4 5 6 7 8 9 10 11 12" 0
if ( $exitstat != 0 ) then
    echo "ERROR: setting months"
    $EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: setting months" -t "e" -x $set_months
    exit (1)
endif

setenv SPINUP_DURATION 0 

# Save spinup array from set_months
set spinup_array = ( $SPINUP_LIST )

## Set naming label
set namelabel = ${SECTOR}_${CASE}_${GRID}

## Record the helper scripts being used
set suffix = _$namelabel.txt
echo Helper scripts used for $SECTOR > $LOGS/helper_scripts_list$suffix
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

## If running from EMF, move old EMF-created scripts to "old"
if ( $?EMF_JOBID ) then
   source $emf_cleanup
   if ( $status != 0 ) then
	echo "ERROR: running EMF script/log cleanup script"
	$EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: running EMF script/log cleanup script" -t "e" -x $emf_cleanup
	exit ( 1 )
   endif
endif

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

set monname = ( jan feb mar apr may jun jul aug sep oct nov dec )
set moncode = ( 01 02 03 04 05 06 07 08 09 10 11 12)

# Loop through months as determined from calling arguments.
set mc = 0
set diff = 0
set g_stdate_all = $G_STDATE

foreach m ( $MONTHS_LIST )   # MONTHS_LIST set by set_months.csh

   @ mc = $mc + 1   # month array counter

   setenv SUBSECT ${SECTOR}_${monname[$m]} # set variable for input/output names
   setenv SRCABBR $SUBSECT                  # set abbreviation for naming log files
   setenv MONTH   ${monname[$m]}           # set variable for month name

   ## Run Smkinven, Grdmat, and Spcmat

   setenv RUN_PART1 Y
   source $ASSIGNS_FILE                    # Invoke Assigns file to set new dates

   ## Set the EMF_PERIOD to the year
   setenv EMF_PERIOD "${MONTH}_${YEAR}"

   # Call EMF Client for current period 
   $EMF_CLIENT -k $EMF_JOBKEY -m "Running SMOKE steps for month $MONTH" -p $EMF_PERIOD   ## log w/ EMF server

   ## Construct inventory list (NOTE: For other sectors, the last argument
   #      on the combine_data.csh script will need to include the month
   #      number) Needs to have a unique env PREFIX that is not used by
   #      intermediary or output files.
   $combine_data EMISINV $LISTFILE list $moncode[$m]
   if ( $status != 0 ) then
       echo "ERROR: Could not run combine_data.csh in $EMF_JOBNAME"
       $EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: Could not run combine_data to create area list: $LISTFILE" -t "e" -x $combine_data
       exit ( 1 )
   endif
   setenv ARINV $LISTFILE  #  Set the area inventory list file to the generic list file

   ## Set period for spcmat for month specific profile ratios in GSPRO_COMBO
   setenv SPCMAT_PERIOD $m
   $EMF_CLIENT -k $EMF_JOBKEY -m "Setting Spcmat period to $SPCMAT_PERIOD" -p $EMF_PERIOD   ## log w/ EMF server

   ## Set period for Smkinven so that the correct emissions column of FF10-format inventories is used
   setenv SMKINVEN_MONTH $m
   echo "SMKINVEN_MONTH set to $SMKINVEN_MONTH"
   $EMF_CLIENT -k $EMF_JOBKEY -m "Setting Smkinven month to $SMKINVEN_MONTH" -p $EMF_PERIOD   ## log w/ EMF server

   # Run programs for "part 1"
   source $smk_run 
   if ( $status != 0 ) then
       echo "ERROR: Running smk_run for part 1 in $EMF_PERIOD"
       $EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: Running smk_run for part 1 in $EMF_PERIOD" -t "e" -x $smk_run -p $EMF_PERIOD
       set exitstat = 1
   endif

   # 6/27/12 C. Allen - nonpt Smkinven logs have millions of warnings in them, and this is basically crashing the log analyzer
   #   and causing all kinds of problems. workaround: filter those warnings out ahead of time.
   
   set inlog = $LOGS/smkinven_${SUBSECT}_${CASE}.log
   echo "inlog = $inlog"
   
   if (-e $inlog) then
     echo "Removing the many 'Monthly inventory is missing' warnings from the Smkinven log for mixed annual/monthly sectors."
     echo "This warning appears millions of times and basically crashes the log analyzer."
     grep 'WARNING: Monthly inventory is missing: Annual inventory will be used' $inlog > /dev/null # check to see if warnings exist
     if ($status == 0) then
       echo "Warning messages found, removing"
       grep -v 'WARNING: Monthly inventory is missing: Annual inventory will be used' $inlog > $inlog.temp
       echo "Note: multiple 'WARNING: Monthly inventory is missing: Annual inventory will be used' warnings to avert log analyzer crash" >> $inlog.temp
       mv -f $inlog.temp $inlog
     endif
   endif

   if ($exitstat != 0) goto end_of_script
   
   setenv RUN_PART1 N

   ## Determine dates to run in this month
   setenv RUN_PART2 Y
   setenv MONTH_ARRAY  $m             # MONTH_ARRAY can have as many months listed as needed
   setenv SPINUP_ARRAY $spinup_array[$mc]

   # Source assigns file BEFORE set_days_v3.csh to set PROCDATES and SMK_RUN_DATES settings
   source $ASSIGNS_FILE                  # Invoke Assigns file to set new dates
   setenv SMK_RUN_DATES $SMK_RUN_DATES_1 

   setenv T_TYPE $L_TYPE                # Set temporal type to type for temporal

   source $set_days   # Call script to set dates for run
   if ( $status != 0 ) then
	echo "ERROR: Running set_days"
	$EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: Running set_days for $MONTH" -t "e" -x $set_days -p $EMF_PERIOD
        exit (1)
   endif

   # Settings so that ESDATE will be first date in procdates file, after Assigns file
   #    is sourced
   set line = `head -1 $SMK_RUN_DATES`
   @ diff = $line[1] - $g_stdate_all
   setenv G_STDATE_ADVANCE $diff

   # Source assigns file AFTER set_days_v3.csh so that it can use PROCDATES file
   source $ASSIGNS_FILE                 # Invoke Assigns file to set new dates

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

   source $qa_run
   # Check status of QA run to see if it worked. Give error if failed
   if ( $status != 0 ) then
      echo "ERROR: Running qa_run for $QA_TYPE" 
      $EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: Running qa_run for $QA_TYPE" -t "e" -x $qa_run -p $EMF_PERIOD 
      set exitstat = 1
      goto end_of_script
   endif  # If qa script failure or not

   # Set path of REPORT for AERMOD helper file scripts
   set var_name = REPORT_`echo $MONTH | tr '[a-z]' '[A-Z]'`
   # C. Allen (7 Apr 2017): If Smkreport is turned off via RUNSET (i.e. if the user only wants to run the post-processing step),
   #   then REPORT25 is never defined.
   if ($?REPORT25) then
      setenv $var_name $REPORT25
   else
      setenv $var_name $REPOUT/$QA_TYPE/rep_${QA_LABEL}_${CASE}_${ESDATE}_${QA_TYPE}_${GRID}.txt
   endif

   unsetenv LOGFILE

   setenv RUN_PART2 N

   unsetenv G_STDATE_ADVANCE

end  # End loop over parts

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
echo "SCRIPT NOTE: Generating nonroad style AERMOD helper files"
$EMF_CLIENT -k $EMF_JOBKEY -m "Generating nonroad style AERMOD helper files"   ## log w/ EMF server
$aermod_nr

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
