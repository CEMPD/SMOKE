#!/bin/tcsh -f

# Version @(#)$Id$
# Path    $Source$
# Date    $Date$

# This script sets up needed environment variables for running SMOKE point
# sources for the steps necessary to generate AERMOD helper files.
# The external helper scripts use SMOKE reports and inputs to generate the
# final AERMOD helper files.
#
# This script is intended to be used with the EMF
# source emissions in SMOKE for the EPA 2014 modeling platform, and 
# calls the scripts that runs the SMOKE programs. 
#
# Script created by : M. Houyoux, Environmental Protection Agency
# October, 2007
# Modified for support of AERMOD by: J. Beidler, CSRA March 2017
#
#*********************************************************************

## log w/ EMF server that script is running
$EMF_CLIENT -k $EMF_JOBKEY -s "Running" 

# set source category
setenv SMK_SOURCE P           # source category to process
setenv MRG_SOURCE $SMK_SOURCE # source category to merge

## time-independent programs (except Smkinven run monthly for CEM import)
setenv RUN_SMKINVEN  Y        #  run inventory import program
setenv RUN_SPCMAT    N
setenv RUN_GRDMAT    Y        #  run gridding matrix program

if (! $?USE_FF10_DAILY_POINT) then
   setenv USE_FF10_DAILY_POINT N
endif

## If PTHOUR_WHOLE_YEAR=Y, all CEMs are processed at once in a single Smkinven 
#  instance, instead of separate Smkinven instances by month. This will produce 
#  a single non-month-specific PHOUR file.
#  C.Allen 3 Aug 2018: This script is for processing one month of Smkinven at a time.
#    This setting is used by combine_data script.
if (! $?PTHOUR_WHOLE_YEAR) then
   setenv PTHOUR_WHOLE_YEAR N
endif

setenv USE_FF10_HOURLY_POINT N

if (! $?DAY_SPECIFIC_YN) setenv DAY_SPECIFIC_YN N
if (! $?HOUR_SPECIFIC_YN) setenv HOUR_SPECIFIC_YN N
set day_specific_yn_save = $DAY_SPECIFIC_YN
set hour_specific_yn_save = $HOUR_SPECIFIC_YN
setenv DAY_SPECIFIC_YN N
setenv HOUR_SPECIFIC_YN N

setenv PROMPTFLAG         N       # Y (never set to Y for batch processing)
setenv AUTO_DELETE        Y       # Y deletes SMOKE I/O API output files (recommended)
setenv AUTO_DELETE_LOG    Y       # Y automatically deletes logs without asking
if ( ! $?DEBUGMODE ) then
   setenv DEBUGMODE          N       # Y changes script to use debugger
endif
# setenv DEBUG_EXE    totalview     # Sets the debugger to use when DEBUGMODE = Y

##############################################################################

switch ( $#argv )
   case 0:
   case 1:
   case 2:
   case 4:
   case 5:
   case 6:
      echo "SCRIPT ERROR: Script requires arguments for grid name and"
      echo "              a <label> argument for labeling the TIMELOG"
      echo " "
      echo "  This script expects to be called using the following argument list:"
      echo "     <grid abbrv> <I/O API gridname> <label>"
      echo " "
      echo "  In the above list, the arguments are defined as follows:"
      echo "     <grid abbrv>       : Grid abbreviation (e.g., 36km)"
      echo "     <I/O API gridname> : I/O API gridname that needs to match entry in the"
      echo "                          GRIDDESC input file"
      echo "     <label>            : label to put on TIMELOG file and helper-scripts list"
      echo " "
      echo "  Example:"
      echo "     <script name> 36km US36KM_148X112 onetime"
      echo "              This example runs the script for the 36km (CONUS) grid, and"
      echo "              gives a label to the TIMELOG file of onetime"
      $EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: smoke script did not receive 3 arguments" -t "e"
      exit( 1 )
endsw

# Get the first two options for the grid abbreviation and I/O API grid
setenv GRID "$argv[1]"
setenv IOAPI_GRIDNAME_1 "$argv[2]"

## source the ASSIGN file
source $ASSIGNS_FILE

## List of all the helper scripts that are run in this script
echo "Script path $SCRIPTS"
set emf_cleanup  = $SCRIPTS/run/emf_cleanup.csh
set set_months   = $SCRIPTS/run/set_months_v4.csh
set set_days     = $SCRIPTS/run/set_days_v5.csh
set timetracker  = $SCRIPTS/run/timetracker_v2.csh
set combine_data = $SCRIPTS/run/combine_data_v6.csh
set smk_run      = $SCRIPTS/run/smk_run_v9.csh
set qa_run       = $SCRIPTS/run/qa_run_v10.csh
set log_analyzer = $SCRIPTS/log_analyzer/log_analyzer.py
set msg_list     = $SCRIPTS/log_analyzer/known_messages.txt
set duplic_chk   = $SCRIPTS/run/duplicate_check.csh
set aermod_egu   = $SCRIPTS/aermod/ptegu.pl
set phour_convert = $SCRIPTS/aermod/convert_phour
set aermod_pt    = $SCRIPTS/aermod/ptnonipm.pl
set aermod_ap    = $SCRIPTS/aermod/ptairport.pl
set aermod_cmv   = $SCRIPTS/aermod/cmv.pl
set aermod_cmvout   = $SCRIPTS/aermod/cmv_output_hourly

## If running from EMF, move old EMF-created scripts to "old"
if ( $?EMF_JOBID ) then
   source $emf_cleanup
   if ( $status != 0 ) then
	echo "ERROR: running EMF script/log cleanup script"
	$EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: running EMF script/log cleanup script" -t "e" -x $emf_cleanup
	exit( 1 )
   endif
endif

set exitstat = 0
if ( $#argv == 3 ) then
   setenv TLABEL $argv[3]
endif

## Set naming label
set namelabel = ${SECTOR}_${CASE}_${GRID}
if ( $?TLABEL ) set namelabel = ${namelabel}_$TLABEL

## Record the helper scripts being used
set suffix = _$namelabel.txt
echo "# Helper scripts used for $SECTOR" > $LOGS/helper_scripts_list$suffix
echo $emf_cleanup >> $LOGS/helper_scripts_list$suffix
echo $set_months >> $LOGS/helper_scripts_list$suffix
echo $timetracker >> $LOGS/helper_scripts_list$suffix
echo $combine_data >> $LOGS/helper_scripts_list$suffix
echo $smk_run >> $LOGS/helper_scripts_list$suffix
echo $qa_run >> $LOGS/helper_scripts_list$suffix
echo $log_analyzer >> $LOGS/helper_scripts_list$suffix
echo $msg_list >> $LOGS/helper_scripts_list$suffix
echo $duplic_chk >> $LOGS/helper_scripts_list$suffix

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
setenv EISECTOR $SECTOR

## Run Smkinven (annual inventory only), Grdmat, and Spcmat
#
setenv RUN_PART1 Y
source $ASSIGNS_FILE                # Invoke Assigns file to set new dates

## Construct inventory list (NOTE: For other sectors, the last argument
#      on the combine_data.csh script will need to include the month
#      number) Needs to have a unique env PREFIX that is not used by
#      intermediary or output files.
$combine_data EMISINV $LISTFILE list
if ( $status != 0 ) then
    echo "ERROR: Could not run combine_data.csh to create point file list:"
    echo "       $LISTFILE"
    $EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: Could not run combine_data to create point file list" -t "e" -x $combine_data
    exit ( 1 )
endif
setenv PTINV $LISTFILE  #  Set the point inventory list file to the generic list file

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

setenv DAY_SPECIFIC_YN $day_specific_yn_save
setenv HOUR_SPECIFIC_YN $hour_specific_yn_save

if ($HOUR_SPECIFIC_YN == Y) then
   setenv RUN_SPCMAT N
   setenv RUN_GRDMAT N

   ## Don't import annual inventory, that's already done
   setenv IMPORT_AVEINV_YN  N           
 
   source $ASSIGNS_FILE                  # Invoke Assigns file to set for current month

   # Rename original smkinven log so that it doesn't get overwritten
   mv -f $LOGS/smkinven_${SRCABBR}_$CASE.log $LOGS/smkinven_${SRCABBR}_annual_$CASE.log

   setenv PTDAY $DAYLISTFILE  #  Set the point daily-specific inventory list file to the generic list file

   # C.Allen 3 Aug 2018: Run one month at a time via SMKINVEN_MONTH parameter
   # Script is currently hardwired to run all months; if necessary in the future we can add command line month selection capability like other scripts have (e.g. -m 1 2 3)
   set monname = ( jan feb mar apr may jun jul aug sep oct nov dec )
   set moncode = ( 01 02 03 04 05 06 07 08 09 10 11 12)
   foreach mon (01 02 03 04 05 06 07 08 09 10 11 12)
     
      setenv MONTH   ${monname[$mon]}         # set variable for month name
      setenv SUBSECT ${SECTOR}_$MONTH
      setenv SRCABBR $SUBSECT
      setenv EMF_PERIOD "${MONTH}_${YEAR}"
      setenv SMKINVEN_MONTH $mon
     
      source $ASSIGNS_FILE                  # Invoke Assigns file to set for current month
      ## Create the month-specific list file of hourly data for import to Smkinven
      if ( $HOUR_SPECIFIC_YN == Y ) then
         set combine_style = Hourly
         if ( $?NAMEBREAK_HOURLY ) then
	     if ($NAMEBREAK_HOURLY == cmv) then
	       set combine_style = Hourly_CMV
	       setenv MULTIFILE_NAMEBREAK 5 # dummy numeric value
	     else if ($NAMEBREAK_HOURLY == cmv19) then
	       set combine_style = Hourly_CMV19
	       setenv MULTIFILE_NAMEBREAK 5 # dummy numeric value
	     else if ($NAMEBREAK_HOURLY == ptegu) then
	       set combine_style = Hourly_egu
	       setenv MULTIFILE_NAMEBREAK 5 # dummy numeric value
	     else
               setenv MULTIFILE_NAMEBREAK $NAMEBREAK_HOURLY
	     endif
         endif
         $combine_data EMISHOUR $HOURLISTFILE list $moncode[$mon] $combine_style
         if ( $status != 0 ) then
             echo "ERROR: Could not run combine_data.csh for hourly list file in $EMF_JOBNAME"
             $EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: Could not run combine_data for hourly list file in $EMF_JOBNAME for $MONTH" -t "e" -p $EMF_PERIOD
             exit ( 1 )
         endif
         setenv PTHOUR $HOURLISTFILE  #  Set the point hour-specific inventory list file to the generic list file
      endif
      ## Run Smkinven to import day-specific and hour-specific inventories
      source $smk_run    # Run program
      if ( $status != 0 ) then
          echo "ERROR: Running Smkinven for daily data import in $EMF_PERIOD"
          $EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: Running smk_run for part 1 in $EMF_PERIOD" -t "e" -p $EMF_PERIOD
          goto end_of_script
      endif

      ## Filter out millions of "all emissions set to 0" warnings from the Smkinven log
      #  to prevent the log analyzer from using gobs of memory
      #$ Rename smkinven log so that filename makes sense
      #  Also, filter out millions of "all emissions set to 0" warnings from the Smkinven log
      #  to prevent the log analyzer from using gobs of memory
   
      if ( -e $LOGS/smkinven_${SRCABBR}_$CASE.log ) then
        grep -v "WARNING: All emissions set to 0. for ORIS:" $LOGS/smkinven_${SRCABBR}_$CASE.log | grep -v "since hourly heat input," | grep -v "gross load, and steam load are missing" > $LOGS/smkinven_${SRCABBR}_daily_$CASE.log
        rm $LOGS/smkinven_${SRCABBR}_$CASE.log
      endif

      setenv PHOUR${mon} $PHOUR

   end
   
endif 

setenv RUN_PART1 N

if ($?RUNWAYS) then
   set cntl_run     = $SCRIPTS/run/cntl_run_v6.csh      # Runs CNTLMAT and GRWINVEN
   setenv RUN_CNTLMAT   Y        #  run control matrix program
   # PART1B runs cntlmat and grwinven for airports
   setenv RUN_PART1B Y
   source $ASSIGNS_FILE   # Invoke Assigns file
   echo "**TEST: $PCMAT"
   echo "Beginning cntl_run";
   $cntl_run
   if ( $status != 0 ) then
       echo "ERROR: running control run to create Projection matrix"
       $EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: running control run to create Projection matrix" \
           -t "e" -x $cntl_run -p $EMF_PERIOD
       set exitstat = 1
       goto end_of_script
   endif

   $EMF_CLIENT -k $EMF_JOBKEY -m "Finished running control run to create Projection matrix" \
       -x $cntl_run -p $EMF_PERIOD 
   echo "Finished cntl_run, $EMF_PERIOD";

   # PART1B runs cntlmat and grwinven, 
   # set to N so doesn't remove projection matrix next time source ASSIGNS
   setenv RUN_PART1B N
endif

# Loop through months as determined from calling arguments
set mc = 0
set diff = 0
set g_stdate_all = $G_STDATE
set anymerge  = N
set anytmprep = N


setenv MONTH   jan         # set variable for month name
setenv SUBSECT ${SECTOR}_$MONTH
setenv SRCABBR $SUBSECT
setenv EMF_PERIOD "${MONTH}_${YEAR}"

setenv RUN_PART2 Y
setenv MONTH_ARRAY  1              # MONTH_ARRAY can have as many months listed as needed

#setenv SPINUP_ARRAY 0
setenv T_TYPE $L_TYPE                 # Set temporal type to type for temporal
setenv SPINUP_ARRAY 0

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

setenv RUN_TEMPORAL Y 

# Source assigns file AFTER set_days_v3.csh so that it can use PROCDATES file to delete Temporal outputs
source $ASSIGNS_FILE                  # Invoke Assigns file to delete Temporal outputs

## Run Temporal (using PROCDATES file - run for all days in month)
source $smk_run       # Run programs
if ( $status != 0 ) then
    echo "ERROR: Running smk_run for part 2 in $EMF_PERIOD" 
    $EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: Running smk_run for part 2 in $EMF_PERIOD" -t "e" -x $smk_run -p $EMF_PERIOD
    set exitstat = 1
    goto end_of_script
endif

setenv RUN_PART2 N

## Generate SMOKE reports for AERMOD
setenv QA_TYPE aermod 
setenv QA_LABEL $SUBSECT           # Used to name the report inputs and outputs
setenv REPLABEL $SUBSECT           # Used internally by Smkreport
setenv RUN_SMKREPORT Y

## Advance date to first representative day so that correct ATSUP is referenced
set line = `head -1 $SMK_RUN_DATES`
@ diff = $line[1] - $g_stdate_all
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
   setenv REPORT $REPOUT/$QA_TYPE/rep_${QA_LABEL}_${CASE}_${QA_TYPE}_${GRID}.txt
endif

setenv REP_SRC ${REPORT26}
setenv REP_XWALK ${REPORT}_src_crosswalk.txt

unsetenv LOGFILE

setenv OUTPUT_DIR $PROJECT_ROOT/$CASE/aermod
if ( $?OUTPUT_DIR ) then
    mkdir -p $OUTPUT_DIR
    chmod ug+rwx $OUTPUT_DIR
    foreach subdir (locations temporal parameters xwalk emis)
        mkdir -p $OUTPUT_DIR/$subdir
        chmod ug+rwx $OUTPUT_DIR/$subdir
    end
endif

setenv PERL5LIB $SCRIPTS/aermod
## Generate AERMOD helper files
if ($hour_specific_yn_save == Y) then
    if (($combine_style == Hourly_CMV) || ($combine_style == Hourly_CMV19)) then
        # CMV
        echo "SCRIPT NOTE: Generating CMV style AERMOD helper files"
        $EMF_CLIENT -k $EMF_JOBKEY -m "Generating CMV style AERMOD helper files"   ## log w/ EMF server
        setenv REPORT $REPORT25
        $aermod_cmv
        # Simple run group region labeling
        #
        if (`echo $REGION_ABBREV | cut -c-4` == "12US") then
            set gl = "12"
        else
            set gl = `echo $REGION_ABBREV | cut -c-3`
        endif
        setenv CMV_SRC_LIST $OUTPUT_DIR/temporal/CMV_${gl}_source_list.csv
        rm -f $OUTPUT_DIR/temporal/CMVU_${gl}_*_hourly.csv
        rm -f $OUTPUT_DIR/temporal/CMVP_${gl}_*_hourly.csv
        $aermod_cmvout
    else
        # Assume all sectors with daily FF10 inventories are EGU
        echo "SCRIPT NOTE: Generating EGU style AERMOD helper files"
        $EMF_CLIENT -k $EMF_JOBKEY -m "Generating EGU style AERMOD helper files"   ## log w/ EMF server
        setenv PHOUR_OUT $IMD_ROOT/$SECTOR/phour_out_${SECTOR}_${CASE}.txt
        $phour_convert
        $aermod_egu
        $EMF_CLIENT -k $EMF_JOBKEY -F $OUTPUT_DIR/locations/ptegu_location.csv \
             -T "AERMOD Helper Point Location (CSV)" -O "EGU Point locations ${namelabel}"
        $EMF_CLIENT -k $EMF_JOBKEY -F $OUTPUT_DIR/locations/ptegu_point_srcparam.csv \
             -T "AERMOD Helper Point Parameters (CSV)" -O "EGU Point parameters ${namelabel}"
        $EMF_CLIENT -k $EMF_JOBKEY -F $OUTPUT_DIR/locations/ptegu_fug_srcparam.csv \
             -T "AERMOD Helper Fugitive Parameters (CSV)" -O "EGU Fugitive parameters ${namelabel}"
        $EMF_CLIENT -k $EMF_JOBKEY -F $OUTPUT_DIR/locations/ptegu_srcid_emis.csv \
             -T "AERMOD Helper Point Emissions (CSV)" -O "EGU source emissions ${namelabel}"
        $EMF_CLIENT -k $EMF_JOBKEY -F $OUTPUT_DIR/locations/ptegu_srcid_xwalk.csv \
             -T "AERMOD Helper Point Crosswalk (CSV)" -O "EGU source crosswalk ${namelabel}"
    endif
else if ($?RUNWAYS) then
    # Assume all sectors run with $RUNWAYS set are airports
    echo "SCRIPT NOTE: Generating Airport style AERMOD helper files"
    $EMF_CLIENT -k $EMF_JOBKEY -m "Generating airport style AERMOD helper files"   ## log w/ EMF server
    $aermod_ap
else if (`echo $SECTOR|cut -c1-6` == "ptegu_") then
     echo "SCRIPT NOTE: Generating EGU style AERMOD helper files"
    $EMF_CLIENT -k $EMF_JOBKEY -m "Generating EGU style AERMOD helper files"   ## log w/ EMF server
    $aermod_egu
    $EMF_CLIENT -k $EMF_JOBKEY -F $OUTPUT_DIR/locations/ptegu_location.csv \
         -T "AERMOD Helper Point Location (CSV)" -O "EGU Point locations ${namelabel}"
    $EMF_CLIENT -k $EMF_JOBKEY -F $OUTPUT_DIR/locations/ptegu_point_srcparam.csv \
         -T "AERMOD Helper Point Parameters (CSV)" -O "EGU Point parameters ${namelabel}"
    $EMF_CLIENT -k $EMF_JOBKEY -F $OUTPUT_DIR/locations/ptegu_fug_srcparam.csv \
         -T "AERMOD Helper Fugitive Parameters (CSV)" -O "EGU Fugitive parameters ${namelabel}"
    $EMF_CLIENT -k $EMF_JOBKEY -F $OUTPUT_DIR/locations/ptegu_srcid_emis.csv \
         -T "AERMOD Helper Point Emissions (CSV)" -O "EGU source emissions ${namelabel}"
    $EMF_CLIENT -k $EMF_JOBKEY -F $OUTPUT_DIR/locations/ptegu_srcid_xwalk.csv \
         -T "AERMOD Helper Point Crosswalk (CSV)" -O "EGU source crosswalk ${namelabel}"
else
    echo "SCRIPT NOTE: Generating Point style AERMOD helper files"
    $EMF_CLIENT -k $EMF_JOBKEY -m "Generating point style AERMOD helper files"   ## log w/ EMF server
    $aermod_pt
    $EMF_CLIENT -k $EMF_JOBKEY -F $OUTPUT_DIR/locations/point_location.csv \
         -T "AERMOD Helper Point Location (CSV)" -O "Point Point locations ${namelabel}"
    $EMF_CLIENT -k $EMF_JOBKEY -F $OUTPUT_DIR/locations/point_point_srcparam.csv \
         -T "AERMOD Helper Point Parameters (CSV)" -O "Point Point parameters ${namelabel}"
    $EMF_CLIENT -k $EMF_JOBKEY -F $OUTPUT_DIR/locations/point_fug_srcparam.csv \
         -T "AERMOD Helper Fugitive Parameters (CSV)" -O "Point Fugitive parameters ${namelabel}"
    $EMF_CLIENT -k $EMF_JOBKEY -F $OUTPUT_DIR/locations/point_srcid_emis.csv \
         -T "AERMOD Helper Point Emissions (CSV)" -O "Point source emissions ${namelabel}"
    $EMF_CLIENT -k $EMF_JOBKEY -F $OUTPUT_DIR/locations/point_srcid_xwalk.csv \
         -T "AERMOD Helper Point Crosswalk (CSV)" -O "Point source crosswalk ${namelabel}"
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
