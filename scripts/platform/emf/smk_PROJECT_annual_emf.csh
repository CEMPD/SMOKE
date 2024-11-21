#!/bin/tcsh -fx

# Version @(#)$Id$
# Path    $Source$
# Date    $Date$

# This script sets up needed environment variables for applying PROJECTION
# and/or CONTROL packets to annual SMOKE (point or area) sources.
# This process requires running SMOKE through the following modules:
# smkinven, spcmat, smkreport (pre-projections QA), cntlmat, and grwinven.
# This script is part of "job 1".  The "job 2" script will run smkinven,
# spcmat, and smkreport on the "projected" and/or "controlled" inventory
# created from this "job 1" script.
#
# This script is intended to be used with the EMF
# source emissions in SMOKE for the EPA 2002 modeling platform, and 
# calls the scripts that runs the SMOKE programs. 
#
# Script created by : M. Houyoux, Environmental Protection Agency
# July, 2007
# Modified to work w/ EMF: A. Zubrow, UNC - IE,  August, 2007
# Modified to apply PROJECTION and/or CONTROL packets by R. Mason, EPA, March, 2008
#
#ISSUES:  do not think we need $USE_CASE_CON with 2 separate jobs now
#*********************************************************************

## log w/ EMF server that script is running
# setenv EMF_CLIENT false (in script that calls this).... to test interactively
# also, build EMF case, which will create a script that calls this (with parameters filled in)
# THEN, setenv EMF_CLIENT false and try executing this!!
$EMF_CLIENT -k $EMF_JOBKEY -s "Running" 

#NOTE: make this script independent of SMK_SOURCE
#setenv SMK_SOURCE A           # source category to process
#setenv MRG_SOURCE $SMK_SOURCE # source category to merge

## month-dependent programs
setenv RUN_SMKINVEN  N        #  run inventory import program
setenv RUN_SPCMAT    N        #  run speciation matrix program
setenv RUN_CNTLMAT   Y        #  run control matrix program
setenv RUN_GRWINVEN  Y        #  run control application program

## quality assurance
set run_smkreport = N      # Y runs reporting for state reports

setenv REGISTER_REPOUT    Y       # Imports Smkreport and Smkmerge reports into EMF
setenv REGISTER_GC_OUT    Y       # Imports Grwinven ORL inventory into EMF
setenv PROMPTFLAG         N       # Y (never set to Y for batch processing)
setenv AUTO_DELETE        Y       # Y deletes SMOKE I/O API output files (recommended)
#setenv AUTO_DELETE        N       # Y deletes SMOKE I/O API output files (recommended)
setenv AUTO_DELETE_LOG    Y       # Y automatically deletes logs without asking
if ( ! $?DEBUGMODE ) then
   setenv DEBUGMODE          N       # Y changes script to use debugger
endif
if ( ! $?DEBUG_EXE ) then
   setenv DEBUG_EXE    totalview     # Sets the debugger to use when DEBUGMODE = Y
endif

##############################################################################

#MODIFY FOR PROJECTIONS-SPECIFIC NUMBER OF ARGUMENTS:
# grid abbrv, I/O API gridname, and spinup are not needed

switch ( $#argv )
   case 0:
   case 1:
      echo "SCRIPT ERROR: Script requires arguments for a SMK_SOURCE"
      echo "              and the -m or -q option with 2 settings."
      echo " "
      echo "  This script expects to be called using one of the following argument lists:"
      echo "     <SMK_SOURCE> -m <monthlist> <label>"
      echo "     <SMK_SOURCE> -q <quarters> <label>"
      echo " "
      echo "  You can either use one approach or the other (differing by the -m or -q options)."
      echo " "
      echo "  In the above list, the arguments are defined as follows:"
      echo "     <SMK_SOURCE>       : SMOKE Source Category: A P or M"
      echo "     <monthlist>        : list of months to run when using the -m option"
      echo "     <quarters>         : list of quarters to run when using the -q option"
      echo "     <label>            : label to put on TIMELOG file and helper-scripts list"
      echo " "
      echo "  Examples:"
      echo "     <script name> A -m '1 2 3' jan-sep"
      echo "              This example runs the script for Jan, Feb, & Mar"
      echo "              for SMOKE area sources, and"
      echo "              gives a label to the TIMELOG file of jan-sep." 
      echo " "
      echo "     <script name> P -q '2 3' apr-sep:"
      echo "               This example runs the script for the 2nd & 3rd quarters,"
      echo "               for SMOKE point sources, and gives"
      echo "               a label to the TIMELOG file of apr-sep."
      $EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: smoke script did not receive more than 1 arguments" -t "e"
      exit( 1 )
endsw

###ram: edit argument calls
## Get the first two options for the grid abbreviation and I/O API grid
setenv SMK_SOURCE "$argv[1]"           # source category to process
setenv MRG_SOURCE $SMK_SOURCE # source category to merge

###############################################################################
# The following environment variables need to be set in the EMF parameters tab:
# these variables should be assigned into the EMF-generated script that calls this script
###############################################################################
## For Cntlmat
#      SMK_AVEDAY_YN (see below)  # Y uses seas emis in assessing cutoff, etc.
#      REPORT_DEFAULTS            # See multi-program controls
#setenv XREF_SICOVERSCC      Y     # not needed for PM NAAQS: most growth factors removed
#setenv COMPARE_REPLACE_CONTROL  Y # Applies CONTROL packet cntls only when Replacement control
#                                   (repflag) > baseline controls
# For Grwinven
#setenv SMK_NUM_CTLMAT       1     # number of control/projection matrices
#setenv SMK_GRWSMKOUT_YN     N     # Y outputs a SMOKE-formatted inventory
#setenv SMK_GRWIDAOUT_YN     N     # Y outputs an IDA-formatted inventory
#setenv SMK_GRWORLOUT_YN     Y     # Y outputs an ORL-formatted inventory
#setenv ORL_NONROAD_OUT      N     # Y outputs NONROAD ORL format
#setenv USE_CASE_CON      N  # set to Y to use $CASE_CON in ASSIGNS file (instead of $CASE)
#setenv SMK_FUTURE_YN Y
#setenv SMK_CONTROL_YN Y
#setenv FYEAR              2020    # future year (in this case ne base)
#setenv CNTLCASE BASE
#setenv FULLSCC_ONLY             Y # assign PROJECTION/CONTROL packet explicitly
#setenv area_sav $AREA
#source $ASSIGNS_FILE   # Invoke Assigns file
#setenv AREA $area_sav
#setenv ACMAT01 $APMAT
#setenv ACMAT02 $ACMAT #ACMAT01= $APMAT when running projections or projections only
###############################################################################

#MAYBE ADD THIS FOR PROJECTION:
setenv RUN_SMKREPORT $run_smkreport

setenv SUBSECT $SECTOR                   # set variable for input/output names
setenv SRCABBR $SUBSECT                  # set abbreviation for naming log files

## source the ASSIGNS file: this needs to be defined for too many steps: simply initialize
if ( ! $?GRID) then
   setenv GRID 36US1        # ASSIGNS.emf demands $GRID be defined
endif
source $ASSIGNS_FILE

#ram: ADD FOR PROJECTIONS, BUT USE ONES YOU NEED + smk_cntl_?.csh
## List of all the helper scripts that are run in this script
set emf_cleanup  = $SCRIPTS/run/emf_cleanup.csh  # For cleanup of EMF-created scripts and stdout logs
set set_months   = $SCRIPTS/run/set_months_v4.csh  # Maybe needed for monthly projections?
set timetracker  = $SCRIPTS/run/timetracker_v2.csh  # Needed - creates TIMELOG file
set combine_data = $SCRIPTS/run/combine_data_v6.csh  # Creates the list files automatically
set smk_run      = $SCRIPTS/run/smk_run_v9.csh
set cntl_run     = $SCRIPTS/run/cntl_run_v6.csh      # Runs CNTLMAT and GRWINVEN

#setenv QA_TYPE  custom                # Used to name the report inputs and outputs
setenv QA_TYPE  proj                  # Used to name the report inputs and outputs
setenv QA_LABEL ${SUBSECT}_preGC           # Used to name the report inputs and outputs
setenv REPLABEL $SUBSECT           # Used internally by Smkreport

mkdir -p $REPOUT/custom/

set qa_run       = $SCRIPTS/run/qa_run_v10.csh

#set set_days     = $SCRIPTS/run/set_days_v3.csh  # Not needed
set log_analyzer = $SCRIPTS/log_analyzer/log_analyzer.py  # checks log files for errors/warnings
set msg_list     = $SCRIPTS/log_analyzer/known_messages.txt  # list of warnings with wildcards and rank

#ADD THIS FOR PROJECTIONS
## If running from EMF, move old EMF-created scripts to "old"
if ( $?EMF_JOBID ) then
   source $emf_cleanup
   if ( $status != 0 ) then
	echo "ERROR: running EMF script/log cleanup script"
	$EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: running EMF script/log cleanup script" -t "e" -x $emf_cleanup
	exit( 1 )
   endif
endif

#IF NEED set_months, use this:

## Invoke script to interpret calling arguments
$EMF_CLIENT -k $EMF_JOBKEY -m "Running set_months" -x $set_months  ## log w/ EMF server
set exitstat = 0
switch ( $#argv )
   case 3:
      source $set_months $argv[2] "$argv[3]"
      if ( $status != 0 ) set exitstat = 1
   breaksw
   case 4:
      source $set_months $argv[2] "$argv[3]"
      if ( $status != 0 ) set exitstat = 1
      setenv TLABEL $argv[4]
   breaksw
endsw
if ( $exitstat != 0 ) then
    echo "ERROR: setting months"
    $EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: setting months" -t "e" -x $set_months  
    exit (1)
endif

#USE THIS FOR PROJECTIONS

## Set naming label: ram, changed to remove GRID
#set namelabel = ${SECTOR}_${CASE}_${GRID}
set namelabel = ${SECTOR}_${CASE}
if ( $?TLABEL ) then
  set namelabel = ${namelabel}_$TLABEL
endif

#ADD SIMILIAR FOR PROJECTIONS:

## Record the helper scripts being used
set suffix = _$namelabel.txt
echo "# Helper scripts used for $SECTOR" > $LOGS/helper_scripts_list$suffix
echo $emf_cleanup >> $LOGS/helper_scripts_list$suffix
echo $set_months >> $LOGS/helper_scripts_list$suffix
echo $timetracker >> $LOGS/helper_scripts_list$suffix
echo $combine_data >> $LOGS/helper_scripts_list$suffix
echo $smk_run >> $LOGS/helper_scripts_list$suffix
echo $qa_run >> $LOGS/helper_scripts_list$suffix
#echo $m3stat >> $LOGS/helper_scripts_list$suffix
#echo $set_days >> $LOGS/helper_scripts_list$suffix
echo $log_analyzer >> $LOGS/helper_scripts_list$suffix
echo $msg_list >> $LOGS/helper_scripts_list$suffix

#ADD THIS FOR PROJECTIONS
## Set Time Log filename and initialize file
setenv TIMELOG $LOGS/timelog_$namelabel.txt

## If TIMELOG_YN = N, don't do timetracker (default is Y)
if (! $?TIMELOG_YN) setenv TIMELOG_YN Y

# Only initialize TIMELOG if it doesn't already exist, since the timeracker 
#   can now delete/add entries to prevent duplicates
if ( ! -e $TIMELOG && $TIMELOG_YN != N ) then 
   $EMF_CLIENT -k $EMF_JOBKEY -m "Initializing Time Log" -x $timetracker  ## log w/ EMF server
   $timetracker Y $TIMELOG
   if ( $status != 0 ) then
	echo "ERROR: running timetracker"
	$EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: running timetracker to initialize time log" -t "e" -x $timetracker
	exit ( 1 )
   endif
endif

# If timelog turned off, unset TIMELOG variable, which will turn timetracker off in the rest of the scripts
if ($TIMELOG_YN == N) unsetenv TIMELOG

#WILL WANT THIS FOR REPORTING. Case WILL NEED TO MATCH APPROACH FOR GSREFTMP_A, GSPROTMP_A, etc.
# Setup GSPRO and GSREF file names
# Set the output file name dates using the most current file matching the prefix
set tmplist = (`env| grep GSPROTMP | cut -d\= -f2` )
set gspro_date = `/bin/ls -1 -t $tmplist | head -1 | sed 's/_/\xa/g' | tail -2 | head -1`
setenv GSPRO $GE_DAT/speciation/$CASE/gspro_${SPC}_${CASE}_${gspro_date}.txt

set tmplist = (`env| grep GSREFTMP | cut -d\= -f2` )
set gsref_date = `/bin/ls -1 -t $tmplist | head -1 | sed 's/_/\xa/g' | tail -2 | head -1`
setenv GSREF $GE_DAT/speciation/$CASE/gsref_${SPC}_${CASE}_${gsref_date}.txt

# Create the GSPRO file from the pieces
$combine_data GSPROTMP $GSPRO concat
if ( $status != 0 ) then
    echo "ERROR: Could not run combine_data.csh to create speciation profile file"
    $EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: Could not run combine_data to create speciation profile file" -t "e" -x $combine_data
    exit ( 1 )
endif

# Create the GSREF file from the pieces
$combine_data GSREFTMP $GSREF concat
if ( $status != 0 ) then
    echo "ERROR: Could not run combine_data.csh to create speciation cross-reference file"
    $EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: Could not run combine_data to create speciation cross-reference file" -t "e" -x $combine_data
    exit ( 1 )
endif

#THIS WOULD BE IN EXSITING SCRIPT IF YOU NEED IT.
## Set up scripting environment variables prior to calling the Assigns file
setenv SUBSECT $SECTOR                   # set variable for input/output names
setenv SRCABBR $SUBSECT                  # set abbreviation for naming log files

## Run Smkinven, Grdmat, and Spcmat

setenv RUN_PART1 Y
source $ASSIGNS_FILE               # Invoke Assigns file

#THIS IS CREATING THE .lst FILE FOR SMKINVEN:

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

echo "script echo: SMK_SOURCE= $SMK_SOURCE INVDIR= $INVDIR SECTOR= $SECTOR SUBSECT= $SUBSECT CASE= $CASE"

if ( $SMK_SOURCE == A ) then
  setenv ARINV $LISTFILE  #  Set the area inventory list file to the generic list file
else
  setenv PTINV $LISTFILE  #  Set the area inventory list file to the generic list file
endif

#YOU'LL WANT TO ADD EMF_PERIOD FOR MESSAGES

## Set the EMF_PERIOD to the year
setenv EMF_PERIOD $YEAR

#ADD MESSAGES LIKE THIS:

# Call EMF Client for one-time steps 
$EMF_CLIENT -k $EMF_JOBKEY -m "Running SMOKE steps for one-time steps" -p $EMF_PERIOD   ## log w/ EMF server

# Run programs for "part 1"
source $smk_run 

echo "OK through this part: immediately after sourcing smk_run"
#ADD STATUS CHECKS FOR ALL HELPER SCRIPTS, and sent ERROR message back.
# Check status of QA run to see if it worked. Give error if failed
if ( $status != 0 ) then
    echo "ERROR: Running smk_run for one-time steps"
    $EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: Running smk_run for one-time steps" -t "e" -x $smk_run -p $EMF_PERIOD
    set exitstat = 1
    goto end_of_script
endif

##############################################################################
##ram: only want to run following reports: 1,2,4,21 ALSO, NOT DEFINING GRID1!!
# question:  what happens if try to run state/MACT report for nonroad or mobile sector?
##############################################################################
## Run QA for part 1, non-gridded reports
## Run only if Smkreport still set to Y (may have been turned off by RUNSET)
## Run only for primary grid
#if ( $RUN_SMKREPORT == Y && $GRID == $GRID_1 ) then
if ( $RUN_SMKREPORT == Y ) then
   setenv QA_TYPE  proj               # Used to name the report inputs and outputs
   setenv QA_LABEL $SUBSECT           # Used to name the report inputs and outputs
   setenv REPLABEL $SUBSECT           # Used internally by Smkreport

   source $qa_run
echo "made it to this part: post-qa_run"

   # Check status of QA run to see if it worked. Give error if failed
   if ( $status != 0 ) then
       echo "ERROR: Running qa_run for $QA_TYPE" 
       $EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: Running qa_run for $QA_TYPE" -t "e" -x $qa_run -p $EMF_PERIOD 
       set exitstat = 1
       goto end_of_script

   ## Otherwise, register internal and external reports with the EMF
   else

      if ( $REGISTER_REPOUT == Y ) then
#         echo "SCRIPT NOTE: Registering internal inventory Smkreport reports"
#         if ( -e $REPORT1  ) $EMF_CLIENT -k $EMF_JOBKEY -F $REPORT1 -T "SMOKE Report" -O "state Smkreport $SECTOR"
#         if ( -e $REPORT2  ) $EMF_CLIENT -k $EMF_JOBKEY -F $REPORT2 -T "SMOKE Report" -O "state/SCC Smkreport $SECTOR"
#         if ( -e $REPORT4  ) $EMF_CLIENT -k $EMF_JOBKEY -F $REPORT4 -T "SMOKE Report" -O "county Smkreport $SECTOR"
#         if ( -e $REPORT9  ) $EMF_CLIENT -k $EMF_JOBKEY -F $REPORT9 -T "SMOKE Report" -O "state/SIC Smkreport $SECTOR"
#         if ( -e $REPORT16 ) $EMF_CLIENT -k $EMF_JOBKEY -F $REPORT16 -T "SMOKE Report" -O "state/SCC/PM25 profile Smkreport $SECTOR"
#         if ( -e $REPORT19 ) $EMF_CLIENT -k $EMF_JOBKEY -F $REPORT19 -T "SMOKE Report" -O "state/NAICS Smkreport $SECTOR"
#         if ( -e $REPORT21 ) $EMF_CLIENT -k $EMF_JOBKEY -F $REPORT21 -T "SMOKE Report" -O "state/MACT Smkreport $SECTOR"

      endif  # if register report outputs
   endif  # If qa script failure or not

endif  # End if primary grid or not

setenv RUN_PART1 N
# PART1B runs cntlmat and grwinven
setenv RUN_PART1B Y

##################################################################################################
# Explicitly define GCNTL file from EMF input SHOULD NOT BE NEEDED
# ram: source cntl_run
# MOVE GRWINVEN OUTPUT TO LOCAL SPACE ON AMBER AFTER MODIFYING HEADER INFORMATION
##################################################################################################
set echo on

if ( $SMK_SOURCE == A ) then
  setenv area_sav $AREA
else
  setenv pnts_sav $PNTS
endif

source $ASSIGNS_FILE   # Invoke Assigns file

##################################################################################################
# the MAT01 variable (from EMF parameters) indicates if first packet is (P)rojection or (C)ontrol
##################################################################################################
if ( $SMK_SOURCE == A ) then
  setenv AREA $area_sav
  if ( $MAT01 == P ) then
    setenv ACMAT01 $APMAT
  else
    setenv ACMAT01 $ACMAT
  endif
  setenv ACMAT02 $ACMAT #ACMAT01= $APMAT when running projections or projections only
else
  setenv PNTS $pnts_sav
  if ( $MAT01 == P ) then
    setenv PCMAT01 $PPMAT
    setenv PCMAT02 $PCMAT
  else
    setenv PCMAT01 $PCMAT
    setenv PCMAT02 $PPMAT
  endif
# MRH insert for testing
#  setenv PCMAT01 $PCMAT
#  setenv PCMAT02 $PPMAT
endif

echo "Beginning cntl_run";
#ram: I think this needs to be defined before grwinven is called
   mkdir -p $SMK_HOME/inputs/$CASE_CON/$SECTOR/work
echo "MAKE PATH WHERE GRWINVEN INVENTORY IS FINALLY PLACED"
source $cntl_run
echo "Finished cntl_run";

set dird = `pwd`
# Move output file to new file name and location
#if ( -e $ARINV_O ) then
if ( $?ARINV_O ) then
  head -5  $ARINV_O > header
  set dlines = (`cat $ARINV_O | wc -l `)
  echo "dlines first= " $dlines
    @ dlines -= 5
  echo "dlines= " $dlines
  tail -$dlines  $ARINV_O > x
  rm  $ARINV_O
endif
if ( $?PTINV_O ) then
  head -5  $PTINV_O > header
  set dlines = (`cat $PTINV_O | wc -l `)
  echo "dlines first= " $dlines
    @ dlines -= 5
  echo "dlines= " $dlines
  tail -$dlines  $PTINV_O > x
  rm  $PTINV_O
endif

echo "#DESC SECTOR is: " $SECTOR > grwinvhead
echo "#DESC CNTLCASE is: " $CASE_CON >> grwinvhead
echo "#DESC $CASE inventory for this $CASE_CON control scenario" >> grwinvhead
echo "#DESC Created $nicedate by $user" >> grwinvhead

if ( $SMK_SOURCE == A ) then
  set tmpfile = $SMK_HOME/inputs/$CASE_CON/$SECTOR/work/arinv_${SECTOR}_${CASE_CON}_${CNTLCASE}_${nicedate}_orl.txt
  cat header grwinvhead x > $tmpfile
else
  set tmpfile = $SMK_HOME/inputs/$CASE_CON/$SECTOR/work/ptinv_${SECTOR}_${CASE_CON}_${CNTLCASE}_${nicedate}_orl.txt
  cat header grwinvhead x > $tmpfile
endif

rm grwinvhead x header

if ( $REGISTER_GC_OUT == Y ) then
  if ( $SMK_SOURCE == A ) then
    if ( $ORL_NONROAD_OUT == Y ) then
      if ( -e $tmpfile) $EMF_CLIENT -k $EMF_JOBKEY -F $tmpfile -N ${SECTOR}_${CASE_CON}_${CNTLCASE} -T "ORL Nonroad Inventory (ARINV)" -O "$SECTOR ${CASE_CON}_${CNTLCASE} Inventory"
    else
      if ( -e $tmpfile) $EMF_CLIENT -k $EMF_JOBKEY -F $tmpfile -N ${SECTOR}_${CASE_CON}_${CNTLCASE} -T "ORL Nonpoint Inventory (ARINV)" -O "$SECTOR ${CASE_CON}_${CNTLCASE} Inventory"
    endif
  else 
    if ( -e $tmpfile) $EMF_CLIENT -k $EMF_JOBKEY -F $tmpfile -N ${SECTOR}_${CASE_CON}_${CNTLCASE} -T "ORL Point Inventory (PTINV)" -O "$SECTOR ${CASE_CON}_${CNTLCASE} Inventory"
  endif
endif


endif
unset echo 
##################################################################################################
#ram : I think it is safe to let the rest stay in this script
##################################################################################################

#THESE ARE EXAMPLES FOR EXTERNAL REPORT REGISTRATION (e.g., large month-specific reports):

# Label for the end of the script, used during script abort
end_of_script:

## Register time log
#echo "SCRIPT NOTE: Registering time log"
#$EMF_CLIENT -k $EMF_JOBKEY -F $TIMELOG -T "SMOKE time log" -N "SMOKE timelog $namelabel" -O "Timelog $namelabel"

## Run log file analyzer
$log_analyzer -k $msg_list --list_unknowns -l 3 -f $REPOUT/log_analyzer/rep_logs_${namelabel}_level3.csv -d $LOGS
if ( $status != 0 ) then
   echo "ERROR: running log_analyzer, level 3"
   $EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: running log_analyzer, level3" -t "e" -x $log_analyzer
   set exitstat = 1

### Register log analyzer output
#else
#   echo "SCRIPT NOTE: Registering log summary, level3"
#   $EMF_CLIENT -k $EMF_JOBKEY -F $REPOUT/log_analyzer/rep_logs_${namelabel}_level3.csv \
#        -T "Log summary level 3" -O "Level 3 log summary ${namelabel}"
endif

$log_analyzer -k $msg_list --list_unknowns -e 1 -f $REPOUT/log_analyzer/rep_logs_${namelabel}_level1.csv -d $LOGS
if ( $status != 0 ) then
   echo "ERROR: The log analyzer has detected warning or error messages in your SMOKE logs which may indicate a problem."
   echo "*****  Please review the priority 0 and 1 messages listed in this report:"
   echo "*****  $REPOUT/log_analyzer/rep_logs_${namelabel}_level1.csv"
   $EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: Level 1 errors found. Check log files, fix inputs, and rerun." -t "e" -x $log_analyzer
   set exitstat = 1

## Register log analyzer output
else
#   echo "SCRIPT NOTE: Registering log summary, level1"
#   $EMF_CLIENT -k $EMF_JOBKEY -F $REPOUT/log_analyzer/rep_logs_${namelabel}_level1.csv \
#        -T "Log summary level 1" -O "Level 1 log summary ${namelabel}"
endif

## Ending of script
#
exit( $exitstat )
