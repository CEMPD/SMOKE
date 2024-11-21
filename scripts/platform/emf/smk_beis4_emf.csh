#!/bin/tcsh -f
#
# Version @(#)$Id$
# Path    $Source$
# Date    $Date$
#
# This script sets up needed environment variables for running BEIS3 from SMOKE
# for an annual simulation that supports either the in-line or traditional 
# CMAQ inputs or CAMx inputs.
#
# Script created by : M. Houyoux, EPA, December 2007
#
#*********************************************************************

## log w/ EMF server that script is running
$EMF_CLIENT -k $EMF_JOBKEY -s "Running" 

## Set source category
setenv SMK_SOURCE    B           # source category to process
setenv MRG_SOURCE    $SMK_SOURCE # source category to merge

## Set programs to run...

## Time-independent programs
setenv RUN_NORMBEIS3 N          # run normalized biogenic emissions program
setenv RUN_NORMBEIS4 Y          # run normalized biogenic emissions program

## Time-dependent programs
set run_tmpbeis3 = N          # run temporal adjustments and speciation program
set run_tmpbeis4 = Y          # run temporal adjustments and speciation program
set run_smkmerge = N          # run merge program
set run_m3stat   = Y

## Program-specific controls (note: mostly defined by EMF calling script)

## For Tmpbeis4
setenv BG_CLOUD_TYPE        1   # method used to calculate PAR (no other options yet available)
setenv BIOMET_SAME          Y   # Y indicates temperature and radiation data in same file   
if (! $?BIOSW_YN) setenv BIOSW_YN             Y   # Y uses BIOSEASON file to set winter or summer factors
if (! $?SUMMER_YN) setenv SUMMER_YN            N   # Y assumes summer factors
if (! $?COMPUTE_SEASON_FUNCTION) setenv COMPUTE_SEASON_FUNCTION N # new setting in BEIS4, don't know what it does, but it must be set to N if BIOSW_YN=Y

## Script settings

# allow the user to turn these things off in his/her case
if (! $?REGISTER_REPOUT) then
   setenv REGISTER_REPOUT    Y       # Imports Smkreport and Smkmerge reports into EMF
endif
if (! $?REGISTER_AQMOUT) then
   setenv REGISTER_AQMOUT    Y       # Imports Smkmerge I/O API outputs in EMF
endif
setenv PROMPTFLAG           N   # Y prompts for user input
setenv AUTO_DELETE          Y   # Y automatically deletes I/O API NetCDF output files
setenv AUTO_DELETE_LOG      Y   # Y automatically deletes log files
setenv DEBUGMODE            N   # Y runs program in debugger
setenv DEBUG_EXE     totalview  # debugger to use when DEBUGMODE = Y

##############################################################################

switch ( $#argv )
   case 0:
   case 1:
   case 2:
   case 3:
      echo "SCRIPT ERROR: Script requires arguments for a grid name"
      echo "              and the -m or -q option with 3 settings."
      echo " "
      echo "  This script expects to be called using one of the following argument lists:"
      echo "     <grid abbrv> <I/O API gridname> -m <monthlist> <spinup> <label>"
      echo "     <grid abbrv> <I/O API gridname> -q <quarters> <spinup> <label>"
      echo " "
      echo "  You can either use one approach or the other (differing by the -m or -q options)."
      echo " "
      echo "  In the above list, the arguments are defined as follows:"
      echo "     <grid abbrv>       : Grid abbreviation (e.g., 36US1)"
      echo "     <I/O API gridname> : I/O API gridname that needs to match entry in the"
      echo "                          GRIDDESC input file"
      echo "     <monthlist>        : list of months to run when using the -m option"
      echo "     <quarters>         : list of quarters to run when using the -q option"
      echo "     <spinup>           : set to number of days between 1 and 20 to run a spinup" 
      echo "                          period (value sets number of days), and N otherwise"
      echo "     <label>            : label to put on TIMELOG file and helper-scripts list"
      echo " "
      echo "  Examples:"
      echo "     <script name> 36US1 36US1_148X112 -m '1 2 3' 0 jan-sep"
      echo "              This example runs the script for Jan, Feb, & Mar"
      echo "              for the 36US1 grid, with no spinup days and"
      echo "              gives a label to the TIMELOG file of jan-sep." 
      echo " "
      echo "     <script name> 12EUS1 12EUS1_279X240 -q '2 3' 10 apr-sep:"
      echo "               This example runs the script for the 2nd & 3rd quarters,"
      echo "               for the 12EUS1 grid, including 10 spin-up days, and gives"
      echo "               a label to the TIMELOG file of apr-sep."
      $EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: smoke script did not receive more than 3 arguments" -t "e"
      exit( 1 )
endsw



## Set Smkmerge settings depending on script setting for output format
#     (Note: these must be set before calling Assigns file)
set aqmodel = `echo $OUTPUT_FORMAT | cut -d" " -f1`
switch ($aqmodel)
   case cmaq:
   case CMAQ:
      setenv MRG_TEMPORAL_YN      Y          # Y merges with hourly emissions
      setenv MRG_SPCMAT_YN        Y          # Y merges with speciation matrix
      setenv MRG_GRDOUT_YN        Y          # Y outputs gridded file
      setenv MRG_GRDOUT_UNIT      moles/s    # units for gridded output file
      setenv MRG_TOTOUT_UNIT      moles/day  # units for state and/or county
      setenv OUT_UNITS            2          # molar output units (1 = moles/hr, 2 = moles/s) for Tmpbeis4
      breaksw
   case camx:
   case CAMX:
   case CAMx:
      setenv MRG_TEMPORAL_YN      Y          # Y merges with hourly emissions
      setenv MRG_SPCMAT_YN        Y          # Y merges with speciation matrix
      setenv MRG_GRDOUT_YN        Y          # Y outputs gridded file
      setenv MRG_GRDOUT_UNIT      moles/hr   # units for gridded output file
      setenv MRG_TOTOUT_UNIT      moles/day  # units for state and/or county
      setenv OUT_UNITS            1          # molar output units (1 = moles/hr, 2 = moles/s) for Tmpbeis4
      setenv RUN_SMK2EMIS         Y          #  run conversion of 2-d to UAM binary
      breaksw
   default:
      echo  ERROR: OUTPUT_FORMAT of \"$OUTPUT_FORMAT\" is not known
      echo "      by main run script $EMF_JOBNAME."
      $EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: OUTPUT_FORMAT = $OUTPUT_FORMAT is not recognized" -t "e"
      exit ( 1 )
endsw

# Get the first two options for the grid abbreviation and I/O API grid
setenv GRID "$argv[1]"
setenv IOAPI_GRIDNAME_1 "$argv[2]"

## [**cvy**] Additional ASSIGNS_FILE call added because otherwise,
## $LOGS (used in TIMELOG definition) is undefined
source $ASSIGNS_FILE

## [**cvy**] Moved these definitions earlier in the script, because
## timetracker is referenced a few lines down from here
## List of all the helper scripts that are run in this script
set emf_cleanup  = $SCRIPTS/run/emf_cleanup.csh
set timetracker  = $SCRIPTS/run/timetracker_v2.csh
set smk_run      = $SCRIPTS/run/smk_run_beis4.csh
set m3stat       = $SCRIPTS/run/m3stat_chk_v6.csh
set set_months   = $SCRIPTS/run/set_months_v4.csh
set set_days     = $SCRIPTS/run/set_days_v5.csh
set combine_data = $SCRIPTS/run/combine_data_v6.csh
set log_analyzer = $SCRIPTS/log_analyzer/log_analyzer.py
set msg_list     = $SCRIPTS/log_analyzer/known_messages.txt
set duplic_chk   = $SCRIPTS/run/duplicate_check.csh

## If running from EMF, move old EMF-created scripts to "old"
if ( $?EMF_JOBID ) then
   source $emf_cleanup
   if ( $status != 0 ) then
	echo "ERROR: running EMF script/log cleanup script"
	$EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: running EMF script/log cleanup script" -t "e" -x $emf_cleanup
	exit( 1 )
   endif
endif

## Invoke script to interpret calling arguments
$EMF_CLIENT -k $EMF_JOBKEY -m "Running set_months" -x $set_months  ## log w/ EMF server
set exitstat = 0
switch ( $#argv )
   case 4:
      source $set_months $argv[3] "$argv[4]"
      if ( $status != 0 ) set exitstat = 1
   breaksw
   case 5: 
      source $set_months $argv[3] "$argv[4]" $argv[5]
      if ( $status != 0 ) set exitstat = 1
   breaksw
   case 6:
      source $set_months $argv[3] "$argv[4]" $argv[5]
      if ( $status != 0 ) set exitstat = 1
      setenv TLABEL $argv[6]
   breaksw
endsw
if ( $exitstat != 0 ) then
    echo "ERROR: setting months"
    $EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: setting months" -t "e" -x $set_months  
    exit (1)
endif

# Set spinup duration - the set_months will have QA'f the $argv[5] value
if ( $#argv > 4 ) setenv SPINUP_DURATION $argv[5]

# Save spinup array from set_months
set spinup_array = ( $SPINUP_LIST )

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
echo $m3stat >> $LOGS/helper_scripts_list$suffix
echo $set_days >> $LOGS/helper_scripts_list$suffix
echo $log_analyzer >> $LOGS/helper_scripts_list$suffix
echo $msg_list >> $LOGS/helper_scripts_list$suffix
echo $duplic_chk >> $LOGS/helper_scripts_list$suffix

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

# Setup GSPRO and GSREF file names
   # Set the output file name dates using the most current file matching the prefix
set tmplist = (`env| grep GSPROTMP | cut -d\= -f2` )
set gspro_date = `/bin/ls -1 -t $tmplist | head -1 | sed 's/_/\xa/g' | tail -2 | head -1`
setenv GSPRO $GE_DAT/speciation/$CASE/gspro_${SECTOR}_${SPC}_${CASE}_${gspro_date}.txt

set tmplist = (`env| grep GSREFTMP | cut -d\= -f2` )
set gsref_date = `/bin/ls -1 -t $tmplist | head -1 | sed 's/_/\xa/g' | tail -2 | head -1`
setenv GSREF $GE_DAT/speciation/$CASE/gsref_${SECTOR}_${SPC}_${CASE}_${gsref_date}.txt

if (! -e $GE_DAT/speciation/$CASE) mkdir $GE_DAT/speciation/$CASE

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

## Check speciation, gridding, and temporal cross-reference files for duplicates
echo "SCRIPT NOTE: Scanning GSREF for duplicate records"
$duplic_chk $GSREF S
if ( $status != 0 ) then
    echo "ERROR: Duplicate records detected in GSREF by duplic_chk.csh, or other script problem"
    $EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: Duplicate records detected in GSREF by duplic_chk.csh, or other script problem" -t "e" -x $duplic_chk
    exit ( 1 )
endif

# Check to see if biogenic county sums variable is defined, if not, assume not wanted
if ( ! $?BIO_COUNTY_SUMS ) then
   echo "SCRIPT WARNING: Assuming biogenic county sums not wanted"
   echo "                since BIO_COUNTY_SUMS is not set to Y."
   setenv BIO_COUNTY_SUMS N
endif

# Check to see if biogenic state sums variable is defined, if not, assume not wanted
if ( ! $?BIO_STATE_SUMS ) then
   echo "SCRIPT WARNING: Assuming biogenic state sums not wanted"
   echo "                since BIO_STATE_SUMS is not set to Y."
   setenv BIO_STATE_SUMS N
endif

# source the ASSIGNS file
source $ASSIGNS_FILE

## Set up scripting environment variables prior to calling the Assigns file
setenv SUBSECT $SECTOR                   # set variable for input/output names
setenv SRCABBR $SUBSECT                  # set abbreviation for naming log files

## Set the EMF_PERIOD to the year
setenv EMF_PERIOD $YEAR

## Set WRF_VERSION variable, used to determine which version of tmpbeis4 to use
#  WRF v3 has different SLTYP definitions and does not have the WSAT_PX variable
#  So there is a different version of tmpbeis4 for WRF v3.x applications
#   This alternate executable (tmpbeis4_wrf3) is activated in the smk_run script if WRF_VERSION=3
#  WRF_VERSION can be explicitly set in the EMF case, but if it is undefined, the script attempts to figure it out from the met path
if (! $?WRF_VERSION) then
   echo "WRF_VERSION undefined, attempting to set based on MET_ROOT path"
   if ($MET_ROOT =~ *WRFv4*) then
      setenv WRF_VERSION 4
   else if ($MET_ROOT =~ *WRFv3*) then
      setenv WRF_VERSION 3
   else if ($MET_ROOT =~ *from_Ramboll*) then
      setenv WRF_VERSION 3
   else if ($MET_ROOT =~ *MYR_*) then
      setenv WRF_VERSION 4
   else
      echo "ERROR: WRF_VERSION must be set so we know which version of TMPBEIS to use, the default version of the WRF v3 version"
      echo "If using WRFv4.1 or later, set WRF_VERSION = 4"
      echo "If using an earlier version of WRF, set WRF_VERSION = 3"
      $EMF_CLIENT -k $EMF_JOBKEY -m "WRF_VERSION must be set so we know which version of TMPBEIS to use" -t "e"
      exit (1)
   endif
   echo "WRF_VERSION = $WRF_VERSION"
endif

## Run Normbeis
#
setenv RUN_PART1 Y
source $ASSIGNS_FILE   # Invoke Assigns file

$EMF_CLIENT -k $EMF_JOBKEY -m "Running smk_run -- Normbeis4" -x $smk_run -p $EMF_PERIOD   ## log w/ EMF server
source $smk_run        # Run programs
if ( $status != 0 ) then
    echo "ERROR: Running smk_run for Normbeis4"
        $EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: Running smk_run for Normbeis4" -t "e"
    set exitstat = 1
    goto end_of_script
endif

setenv RUN_PART1 N

## Loop through days to run Tmpbeis4 and Smkmerge
$EMF_CLIENT -k $EMF_JOBKEY -m "Running smk_run -- Tmpbeis4 for all days" -x $smk_run -p $EMF_PERIOD   ## log w/ EMF server

set mrg_repsta_save = $MRG_REPSTA_YN
set mrg_custom = $SMKMERGE_CUSTOM_OUTPUT

set monname = ( jan feb mar apr may jun jul aug sep oct nov dec )

# Loop through months as determined from calling arguments.
set mc = 0
set diff = 0
set g_stdate_all = $G_STDATE

setenv INITIAL_RUN  Y   # Y: running first day of scenario, N for subsequent days

foreach mon ( $MONTHS_LIST ) 

   @ mc = $mc + 1

   setenv MONTH ${monname[$mon]}           # set variable for month name

   ## Determine dates to run in this month
   setenv MONTH_ARRAY  $mon     # MONTH_ARRAY can have as many months listed as needed
   setenv SPINUP_ARRAY $spinup_array[$mc]
   
   setenv T_TYPE $M_TYPE                 # Set temporal type to type for merge
   setenv SMK_RUN_DATES $SMK_RUN_DATES_2 # Set with file name using MTYPE

   ## Set the EMF_PERIOD for this month and year
   setenv EMF_PERIOD "${MONTH}_${YEAR}"

   # Call EMF Client for current period 
   $EMF_CLIENT -k $EMF_JOBKEY -m "Running SMOKE steps for month $MONTH" -p $EMF_PERIOD   ## log w/ EMF server

   source $set_days   # Call script to set dates for run
   # Check status of run to see if it worked. Give error if failed
   if ( $status != 0 ) then
      echo "ERROR: Running set_days in $EMF_PERIOD"
      $EMF_CLIENT - k $EMF_JOBKEY -m "ERROR: Running set_days" -t "e" -x $set_days -p $EMF_PERIOD
      exit (1)
   endif

   ## Determine the number of days to run
   set ndays = `cat $SMK_RUN_DATES | wc -l`
   echo "Number of dates to run for $MONTH" : $ndays

   set n = 0
   set diff = 0
   set g_stdate_sav = $g_stdate_all

   # Loop through days to run during the month.
   while ( $n < $ndays )

      @ n = $n + 1   

      set line = `head -$n $SMK_RUN_DATES | tail -1`
      @ diff = $line[1] - $g_stdate_sav

      setenv G_STDATE_ADVANCE $diff

      setenv RUN_PART2              Y
      setenv RUN_PART4              Y
      setenv MRG_REPSTA_YN          $mrg_repsta_save  # Restore original setting
      setenv SMKMERGE_CUSTOM_OUTPUT $mrg_custom       # Restore original setting

      setenv RUN_TMPBEIS4 $run_tmpbeis4
      setenv RUN_SMKMERGE $run_smkmerge
      setenv RUN_M3STAT   $run_m3stat
      
#      @ cnt = $cnt + $NDAYS
      source $ASSIGNS_FILE   # Invoke Assigns file to set new dates 

      ## Set EMF_PERIOD for this day
      setenv EMF_PERIOD $ESDATE

      # Run Tmpbeis4 and Smkmerge to generate gridded data
      source $smk_run                # Run programs
      if ( $status != 0 ) then
         echo "ERROR: running Tmpbeis4 or Smkmerge for $ESDATE"
         $EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: running Tmpbeis4 or Smkmerge for $ESDATE" -t "e" -p $EMF_PERIOD
         set exitstat = 1
	 goto end_of_script
      endif

      # Run m3stat script on Smkmerge output file
      if ( $RUN_M3STAT == Y ) then
         $m3stat $B3GTS_L
         if ( $status != 0 ) then
            echo "ERROR: running m3stat_chk for $ESDATE"
            $EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: running m3stat_chk for $ESDATE" -t "e"
            set exitstat = 1
	    goto end_of_script
         endif
      endif

      # Run Smkmerge to generate daily totals
      # Separate Smkmerge runs for county and state totals

      ## County reports, if Tmpbeis4 has already run
      if ( $BIO_COUNTY_SUMS == Y && -e $B3GTS_L ) then
         ## RUN_SMKMERGE needs to be Y for county or state reports to run
         set run_smkmerge_save = $RUN_SMKMERGE
      
         setenv RUN_PART2 N 
         setenv RUN_SMKMERGE           Y
         setenv MRG_GRDOUT_YN          N        # Y produces a gridded output file
         setenv MRG_REPCNY_YN          Y        # Y produces a report of emission totals by county
         setenv MRG_REPSTA_YN          N        # Y produces a report of emission totals by state
         setenv MRG_GRDOUT_UNIT        tons/day # units for the gridded output file
         setenv MRG_TOTOUT_UNIT        tons/day # units for the state and county reports
#         setenv SMKMERGE_CUSTOM_OUTPUT N        # Y uses "BOUT" name instead of more complicated one.
      
         ## Sets the filename for REPBGTS_S to "county"
         setenv BEIS_REP_TYPE county
         source $ASSIGNS_FILE

         source $smk_run                # Run programs
         if ( $status != 0 ) then
	    echo "ERROR: running Smkmerge for county sums for $ESDATE"
   	    $EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: running Smkmerge for county sums for $ESDATE" -t "e" -p $EMF_PERIOD
            set exitstat = 1
	    goto end_of_script
         endif
      
         setenv RUN_SMKMERGE $run_smkmerge_save
      
      endif
      
      ## State reports, if Tmpbeis4 has already run
      if ( $BIO_STATE_SUMS == Y && -e $B3GTS_L ) then
         ## RUN_SMKMERGE needs to be Y for county or state reports to run
         set run_smkmerge_save = $RUN_SMKMERGE
      
         setenv RUN_PART2 N 
         setenv RUN_SMKMERGE           Y
         setenv MRG_GRDOUT_YN          N        # Y produces a gridded output file
         setenv MRG_REPCNY_YN          N        # Y produces a report of emission totals by county
         setenv MRG_REPSTA_YN          Y        # Y produces a report of emission totals by state
         setenv MRG_GRDOUT_UNIT        tons/day # units for the gridded output file
         setenv MRG_TOTOUT_UNIT        tons/day # units for the state and county reports
#         setenv SMKMERGE_CUSTOM_OUTPUT N        # Y uses "BOUT" name instead of more complicated one.
      
         ## Sets the filename for REPBGTS_S to "state"
         setenv BEIS_REP_TYPE state
         source $ASSIGNS_FILE

         source $smk_run                # Run programs
         if ( $status != 0 ) then
	    echo "ERROR: running Smkmerge for state sums for $ESDATE"
   	    $EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: running Smkmerge for state sums for $ESDATE" -t "e" -p $EMF_PERIOD
            set exitstat = 1
	    goto end_of_script
	 endif
      
         setenv RUN_SMKMERGE $run_smkmerge_save
      
      endif
       
      if ($RUN_TMPBEIS4 == Y) setenv INITIAL_RUN   N
      
   end # End loop over days
   
   unsetenv G_STDATE_ADVANCE
   
end

setenv RUN_PART2 N 
setenv RUN_PART4 N 

# Register model-ready Tmpbeis4 outputs in EMF, if RUN_TMPBEIS4 wasn't overridden to N by RUNSET
if ( $RUN_TMPBEIS4 == Y && $REGISTER_AQMOUT == Y ) then
#   echo "SCRIPT NOTE: Registering Tmpbeis4 outputs"
#   $EMF_CLIENT -k $EMF_JOBKEY \
#               -D $PREMERGED \
#               -P "emis_mole_${SECTOR}_*_${GRID}_*.ncf" \
#               -T "CMAQ Model Ready Emissions: Sector-specific (External)" \
#               -N "Model-ready CMAQ $SECTOR $CASE" \
#               -O "AQM-ready data $SECTOR"
endif
	
# Register Smkmerge county reports in EMF
if ( $BIO_COUNTY_SUMS == Y && $REGISTER_REPOUT == Y ) then
#   echo "SCRIPT NOTE: Registering Smkmerge county reports"
#   $EMF_CLIENT -k $EMF_JOBKEY \
#               -D $REPOUT/smkmerge/$SECTOR \
#               -P "rep_mass_county_${SECTOR}_*_${GRID}_*.txt" \
#               -T "SMOKE Report (External)" \
#               -N "Smkmerge $SECTOR $CASE county reports" \
#               -O "county Smkmerge $SECTOR (External)"
endif
	
# Register Smkmerge state reports in EMF
if ( $BIO_STATE_SUMS == Y && $REGISTER_REPOUT == Y ) then
#   echo "SCRIPT NOTE: Registering Smkmerge state reports"
#   $EMF_CLIENT -k $EMF_JOBKEY \
#               -D $REPOUT/smkmerge/$SECTOR \
#               -P "rep_mass_state_${SECTOR}_*_${GRID}_*.txt" \
#               -T "Smkmerge Report state (External Multifile)" \
#               -N "Smkmerge $SECTOR $CASE state reports" \
#               -O "state Smkmerge $SECTOR (External)"
endif
      
# Label for the end of the script, used during script abort
end_of_script:

## Register time log
#echo "SCRIPT NOTE: Registering time log"
#$EMF_CLIENT -k $EMF_JOBKEY -F $TIMELOG -T "SMOKE time log (External)" -N "SMOKE timelog $namelabel" -O "Timelog $namelabel (External)"

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

## If log analyzer returns exit status 10, then the log analyzer ran to completion, but found some errors or unrecognized
#  warnings. In this case, we want this script to exit with "Failed" status, but we still want to register the level 1 report in the EMF.
#  If log analyzer returns exit status 1, it didn't work at all, so bomb out immediately.

$log_analyzer -k $msg_list --list_unknowns -e 1 -f $REPOUT/log_analyzer/rep_logs_${namelabel}_level1.csv -d $LOGS
if ( $status == 10 ) then
#   echo "SCRIPT NOTE: Registering log summary, level1"
#   $EMF_CLIENT -k $EMF_JOBKEY -F $REPOUT/log_analyzer/rep_logs_${namelabel}_level1.csv \
#        -T "Log summary level 1" -O "Level 1 log summary ${namelabel}"
	
   echo "ERROR: The log analyzer has detected warning or error messages in your SMOKE logs which may indicate a problem."
   echo "*****  Please review the priority 0 and 1 messages listed in this report:"
   echo "*****  $REPOUT/log_analyzer/rep_logs_${namelabel}_level1.csv"
   $EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: Level 1 errors or Level 0 warnings found. Check log files, fix inputs, and rerun." -t "e" -x $log_analyzer
   set exitstat = 1

## Register log analyzer output
else if ( $status != 0 ) then
   echo "ERROR: running log_analyzer, level 1"
   $EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: running log_analyzer, level 1" -t "e" -x $log_analyzer
   set exitstat = 1
else
#   echo "SCRIPT NOTE: Registering log summary, level1"
#   $EMF_CLIENT -k $EMF_JOBKEY -F $REPOUT/log_analyzer/rep_logs_${namelabel}_level1.csv \
#        -T "Log summary level 1" -O "Level 1 log summary ${namelabel}"
endif

## Ending of script
#
exit( $exitstat )
