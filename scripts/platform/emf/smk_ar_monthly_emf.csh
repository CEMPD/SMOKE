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

## month-dependent programs
set run_smkinven = Y        #  run inventory import program
set run_spcmat   = Y        #  run speciation matrix program
set run_grdmat   = Y        #  run gridding matrix program

## time-dependent programs
set run_temporal = Y       #  run temporal allocation program
set run_smkmerge = Y       #  run merge program
set run_smk2emis = N       #  run conversion of 2-d to UAM binary

## quality assurance
set run_smkreport = Y      # Y runs reporting for state reports
set run_m3stat    = Y      # Y runs the m3stat program on the Smkmerge outputs

# allow the user to turn these things off in his/her case
if (! $?REGISTER_REPOUT) then
   setenv REGISTER_REPOUT    Y       # Imports Smkreport and Smkmerge reports into EMF
endif
if (! $?REGISTER_AQMOUT) then
   setenv REGISTER_AQMOUT    Y       # Imports Smkmerge I/O API outputs in EMF
endif
setenv PROMPTFLAG         N       # Y (never set to Y for batch processing)
setenv AUTO_DELETE        Y       # Y deletes SMOKE I/O API output files (recommended)
setenv AUTO_DELETE_LOG    Y       # Y automatically deletes logs without asking
if ( ! $?DEBUGMODE ) then
   setenv DEBUGMODE          N       # Y changes script to use debugger
endif
if ( ! $?DEBUG_EXE ) then
   setenv DEBUG_EXE     idb    # Sets the debugger to use when DEBUGMODE = Y
endif

# Preserve original setting for DAY_SPECIFIC_YN
set day_specific_sav = N
if ($?DAY_SPECIFIC_YN) then
  if ($DAY_SPECIFIC_YN == Y || $DAY_SPECIFIC_YN == TRUE) then
    set day_specific_sav = Y
  endif
endif
    

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
      breaksw
   case camx:
   case CAMX:
   case CAMx:
      setenv MRG_TEMPORAL_YN      Y          # Y merges with hourly emissions
      setenv MRG_SPCMAT_YN        Y          # Y merges with speciation matrix
      setenv MRG_GRDOUT_YN        Y          # Y outputs gridded file
      setenv MRG_GRDOUT_UNIT      moles/hr   # units for gridded output file
      setenv MRG_TOTOUT_UNIT      moles/day  # units for state and/or county
      set    run_smk2emis =       Y        
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

setenv RUN_SMKINVEN  $run_smkinven
setenv RUN_SPCMAT    $run_spcmat
setenv RUN_GRDMAT    $run_grdmat
setenv RUN_TEMPORAL  $run_temporal
setenv RUN_SMKMERGE  $run_smkmerge
if (! $?RUN_SMKREPORT) then
#  if not set already, set it from the variable.
   setenv RUN_SMKREPORT $run_smkreport
else
#  this allows for someone to override the run report setting in their case
   set run_smkreport = $RUN_SMKREPORT
   if ($RUN_SMKREPORT == 'N') then
       setenv REGISTER_REPOUT    N
#      don't register the reports if you are not creating them
   endif
endif
setenv RUN_SMK2EMIS  $run_smk2emis

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
set python_annl  = $SCRIPTS/annual_report/annual_report_v2.py # can read county or state

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
if ( $?TLABEL ) set namelabel = ${namelabel}_$TLABEL

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
echo $python_annl >> $LOGS/helper_scripts_list$suffix
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

set tmplist = (`env| grep GSCNVTMP | cut -d\= -f2` )
set merge_gscnv = N
if ($#tmplist > 0 ) then
   set gscnv_date = `/bin/ls -1 -t $tmplist | head -1 | sed 's/_/\xa/g' | tail -2 | head -1`
   setenv GSCNV $GE_DAT/speciation/$CASE/gscnv_${SECTOR}_${SPC}_${CASE}_${gscnv_date}.txt
   set merge_gscnv = Y
endif

set tmplist = (`env| grep GSPRO_COMBOTMP | cut -d\= -f2` )
set merge_combo = N
if ($#tmplist > 0 ) then
   set combo_date = `/bin/ls -1 -t $tmplist | head -1 | sed 's/_/\xa/g' | tail -2 | head -1`
   setenv GSPRO_COMBO $GE_DAT/speciation/$CASE/gspro_combo_${SECTOR}_${SPC}_${CASE}_${combo_date}.txt
   set merge_combo = Y
endif

if (! -e $GE_DAT/speciation/$CASE) mkdir $GE_DAT/speciation/$CASE

if (! $?GSPRO) then # Create the GSPRO file from the pieces, if full GSPRO not already defined
   set tmplist = (`env| grep GSPROTMP | cut -d\= -f2` )
   set gspro_date = `/bin/ls -1 -t $tmplist | head -1 | sed 's/_/\xa/g' | tail -2 | head -1`
   setenv GSPRO $GE_DAT/speciation/$CASE/gspro_${SECTOR}_${SPC}_${CASE}_${gspro_date}.txt

   # Create the GSPRO file from the pieces
   $combine_data GSPROTMP $GSPRO concat
   if ( $status != 0 ) then
       echo "ERROR: Could not run combine_data.csh to create speciation profile file"
       $EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: Could not run combine_data to create speciation profile file" -t "e" -x $combine_data
       exit ( 1 )
   endif
endif

# Use GSREF parameter if it exists, otherwise build from GSREFTMP components (updated 26 Jul 2018)
if (! $?GSREF) then
  set tmplist = (`env| grep GSREFTMP | cut -d\= -f2` )
  set gsref_date = `/bin/ls -1 -t $tmplist | head -1 | sed 's/_/\xa/g' | tail -2 | head -1`
  setenv GSREF $GE_DAT/speciation/$CASE/gsref_${SECTOR}_${SPC}_${CASE}_${gsref_date}.txt

  # Create the GSREF file from the pieces
  $combine_data GSREFTMP $GSREF concat
  if ( $status != 0 ) then
      echo "ERROR: Could not run combine_data.csh to create speciation cross-reference file"
      $EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: Could not run combine_data to create speciation cross-reference file" -t "e" -x $combine_data
      exit ( 1 )
  endif
endif

# Create the GSCNV file from the pieces
if ( $merge_gscnv == Y ) then
   $combine_data GSCNVTMP $GSCNV concat
   if ( $status != 0 ) then
       echo "ERROR: Could not run combine_data.csh to create voc-to-tog conversion file (GSCNV)"
       $EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: Could not run combine_data to create voc-to-tog conversion file (GSCNV)" -t "e" -x $combine_data
       exit ( 1 )
   endif
endif

# Create the GSPRO_COMBO file from the pieces
if ( $merge_combo == Y ) then
   $combine_data GSPRO_COMBOTMP $GSPRO_COMBO concat
   if ( $status != 0 ) then
       echo "ERROR: Could not run combine_data.csh to create profiles combinations file (GSPRO_COMBO)"
       $EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: Could not run combine_data to create profiles combinations file (GSPRO_COMBO)" -t "e" -x $combine_data
       exit ( 1 )
   endif
endif

## Check speciation, gridding, and temporal cross-reference files for duplicates
echo "SCRIPT NOTE: Scanning GSREF for duplicate records"
$duplic_chk $GSREF S
if ( $status != 0 ) then
    echo "ERROR: Duplicate records detected in GSREF by duplic_chk.csh, or other script problem"
    $EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: Duplicate records detected in GSREF by duplic_chk.csh, or other script problem" -t "e" -x $duplic_chk
    exit ( 1 )
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

   if ($m == 0) then # For creating annual Smkreports, support added by C.Allen (14 Feb 2019)
      setenv SUBSECT ${SECTOR}_annual # set variable for input/output names
      setenv MONTH   annual           # set variable for month name
      setenv DAY_SPECIFIC_YN N        # Don't do daily processing for the annual reports
   else
      setenv SUBSECT ${SECTOR}_${monname[$m]} # set variable for input/output names
      setenv MONTH   ${monname[$m]}           # set variable for month name
      setenv DAY_SPECIFIC_YN $day_specific_sav # reset to original parameter in case annual reports were already generated in this run
   endif
   setenv SRCABBR $SUBSECT                  # set abbreviation for naming log files

   ## Run Smkinven, Grdmat, and Spcmat

   setenv RUN_SMKINVEN  $run_smkinven
   setenv RUN_SPCMAT    $run_spcmat
   setenv RUN_GRDMAT    $run_grdmat
   setenv RUN_SMKREPORT $run_smkreport

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

   # Christine Allen added daily inventory support 29 Jan 2024. Previously, only annual area script had this
   ## If DAY_SPECIFIC_YN = Y, then set up list file for daily FF10
   if ( $?DAY_SPECIFIC_YN ) then
      if ( $DAY_SPECIFIC_YN == Y ) then
         $combine_data EMISDAY $DAYLISTFILE list $moncode[$m]
         if ( $status != 0 ) then
             echo "ERROR: Could not run combine_data.csh to create area daily inventory file list:"
             echo "       $DAYLISTFILE"
             $EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: Could not run combine_data to create area file list" -t "e" -x $combine_data
             exit ( 1 )
         endif
         setenv ARDAY $DAYLISTFILE  #  For FF10 daily nonpoint
      endif
   endif

   ## If HOUR_SPECIFIC_YN = Y, then set up list file for hourly FF10
   if ( $?HOUR_SPECIFIC_YN ) then
      if ( $HOUR_SPECIFIC_YN == Y ) then
         $combine_data EMISHOUR $HOURLISTFILE list $moncode[$m]
         if ( $status != 0 ) then
             echo "ERROR: Could not run combine_data.csh to create area hourly inventory file list:"
             echo "       $HOURLISTFILE"
             $EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: Could not run combine_data to create area file list" -t "e" -x $combine_data
             exit ( 1 )
         endif
         setenv ARHOUR $HOURLISTFILE  #  For FF10 hourly nonpoint
      endif
   endif

   ## Set period for spcmat for month specific profile ratios in GSPRO_COMBO
   setenv SPCMAT_PERIOD $m
   $EMF_CLIENT -k $EMF_JOBKEY -m "Setting Spcmat period to $SPCMAT_PERIOD" -p $EMF_PERIOD   ## log w/ EMF server

   ## Set period for Smkinven so that the correct emissions column of FF10-format inventories is used
   setenv SMKINVEN_MONTH $m
   echo "SMKINVEN_MONTH set to $SMKINVEN_MONTH"
   $EMF_CLIENT -k $EMF_JOBKEY -m "Setting Smkinven month to $SMKINVEN_MONTH" -p $EMF_PERIOD   ## log w/ EMF server

#   echo "ADAY $ADAY" #test
#   echo "AHOUR $AHOUR" #test

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
   
   ## Run QA for part 1, non-gridded reports
   ## Run only if Smkreport still set to Y (may have been turned off by RUNSET)
    if ( $RUN_SMKREPORT == Y ) then

      setenv QA_TYPE  inv                # Used to name the report inputs and outputs
      setenv QA_LABEL $SUBSECT           # Used to name the report inputs and outputs
      setenv REPLABEL $SUBSECT           # Used internally by Smkreport
      source $qa_run  
   
      # Check status of QA run to see if it worked. Give error if failed
      if ( $status != 0 ) then
          echo "ERROR: Running qa_run, part 1 in $EMF_PERIOD" 
          $EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: Running qa_run, part 1 in $EMF_PERIOD" -t "e" -x $qa_run -p $EMF_PERIOD
          set exitstat = 1
          goto end_of_script

      # Otherwise, register internal and external reports with the EMF
      else

         if ( $REGISTER_REPOUT == Y ) then
#            echo "SCRIPT NOTE: Registering internal inventory Smkreport reports"
#            echo "REGISTER_REPOUT = " $REGISTER_REPOUT
#            if ( -e $REPORT1  ) $EMF_CLIENT -k $EMF_JOBKEY -F $REPORT1 -T "Smkreport state" -O "state Smkreport $SUBSECT"
#            if ( -e $REPORT2  ) $EMF_CLIENT -k $EMF_JOBKEY -F $REPORT2 -T "Smkreport state-SCC" -O "state/SCC Smkreport $SUBSECT"
#            if ( -e $REPORT4  ) $EMF_CLIENT -k $EMF_JOBKEY -F $REPORT4 -T "Smkreport county" -O "county Smkreport $SUBSECT"
#            if ( -e $REPORT16 ) $EMF_CLIENT -k $EMF_JOBKEY -F $REPORT16 -T "Smkreport state-SCC-spec_profile" -O "state/SCC/PM25 profile Smkreport $SUBSECT"
#
#           echo "SCRIPT NOTE: Registering external Smkreport reports"
#            if ( -e $REPORT8 ) $EMF_CLIENT -k $EMF_JOBKEY -F $REPORT8 -T "Smkreport county-SCC (External)" \
#                        -N "rep_${SUBSECT}_inv_county_scc_$CASE (External)" \
#                        -O "county/SCC Smkreport $SUBSECT (External)"
#            if ( -e $REPORT9  ) $EMF_CLIENT -k $EMF_JOBKEY -F $REPORT9 -T "Smkreport state-SIC (External)" \
#                        -N "rep_${SUBSECT}_inv_state_sic_$CASE (External)" \
#	                -O "state/SIC Smkreport $SUBSECT (External)"
#            if ( -e $REPORT19 ) $EMF_CLIENT -k $EMF_JOBKEY -F $REPORT19 -T "Smkreport state-NAICS (External)" \
#                        -N "rep_${SUBSECT}_inv_state_naics_$CASE (External)" \
#	                -O "state/NAICS Smkreport $SUBSECT (External)"
#            if ( -e $REPORT21 ) $EMF_CLIENT -k $EMF_JOBKEY -F $REPORT21 -T "Smkreport state-MACT (External)" \
#                        -N "rep_${SUBSECT}_inv_state_mact_$CASE (External)" \
#                        -O "state/MACT Smkreport $SUBSECT (External)"
			
	 endif # register_repout

      endif # if qa script failure or not

     ## Continue QA, but for state/SCC VOC profiles (PM2.5 profiles were in the first run, and can't 
     #      do both in same Smkreport run)
     #  This section was edited on 7 Feb 2019 in support of emissions modeling workgroup needs
     #  If REPCONFIG_INV2 exists, run Smkreport again
     #  If REPCONFIG_INV3 exists, run Smkreport once more
     #  Same for 2A, 2B, and so on
     if ( $?REPCONFIG_INV2 ) then
       setenv QA_TYPE  inv2               # Used to name the report inputs and outputs
       setenv QA_LABEL $SUBSECT           # Used to name the report inputs and outputs
       setenv REPLABEL $SUBSECT           # Used internally by Smkreport
       source $qa_run  

       # Check status of QA run to see if it worked. Give error if failed
       if ( $status != 0 ) then
           echo "ERROR: Running qa_run for $QA_TYPE" 
           $EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: Running qa_run for $QA_TYPE" -t "e" -x $qa_run -p $EMF_PERIOD 
           set exitstat = 1
           goto end_of_script

       endif  # If qa script failure or not
     endif      # If REPCONFIG_INV2

     if ( $?REPCONFIG_INV3 ) then
       setenv QA_TYPE  inv3               # Used to name the report inputs and outputs
       setenv QA_LABEL $SUBSECT           # Used to name the report inputs and outputs
       setenv REPLABEL $SUBSECT           # Used internally by Smkreport
       source $qa_run  

       # Check status of QA run to see if it worked. Give error if failed
       if ( $status != 0 ) then
           echo "ERROR: Running qa_run for $QA_TYPE" 
           $EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: Running qa_run for $QA_TYPE" -t "e" -x $qa_run -p $EMF_PERIOD 
           set exitstat = 1
           goto end_of_script

       endif  # If qa script failure or not
     endif      # If REPCONFIG_INV3
   
     if ( $?REPCONFIG_INV2_A ) then
       setenv QA_TYPE  inv2a              # Used to name the report inputs and outputs
       setenv QA_LABEL $SUBSECT           # Used to name the report inputs and outputs
       setenv REPLABEL $SUBSECT           # Used internally by Smkreport
       source $qa_run  

       # Check status of QA run to see if it worked. Give error if failed
       if ( $status != 0 ) then
           echo "ERROR: Running qa_run for $QA_TYPE" 
           $EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: Running qa_run for $QA_TYPE" -t "e" -x $qa_run -p $EMF_PERIOD 
           set exitstat = 1
           goto end_of_script

       endif  # If qa script failure or not
     endif      # If REPCONFIG_INV2

     if ( $?REPCONFIG_INV3_A ) then
       setenv QA_TYPE  inv3a               # Used to name the report inputs and outputs
       setenv QA_LABEL $SUBSECT           # Used to name the report inputs and outputs
       setenv REPLABEL $SUBSECT           # Used internally by Smkreport
       source $qa_run  

       # Check status of QA run to see if it worked. Give error if failed
       if ( $status != 0 ) then
           echo "ERROR: Running qa_run for $QA_TYPE" 
           $EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: Running qa_run for $QA_TYPE" -t "e" -x $qa_run -p $EMF_PERIOD 
           set exitstat = 1
           goto end_of_script

       endif  # If qa script failure or not
     endif      # If REPCONFIG_INV3

     if ( $?REPCONFIG_INV2_B ) then
       setenv QA_TYPE  inv2b              # Used to name the report inputs and outputs
       setenv QA_LABEL $SUBSECT           # Used to name the report inputs and outputs
       setenv REPLABEL $SUBSECT           # Used internally by Smkreport
       source $qa_run  

       # Check status of QA run to see if it worked. Give error if failed
       if ( $status != 0 ) then
           echo "ERROR: Running qa_run for $QA_TYPE" 
           $EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: Running qa_run for $QA_TYPE" -t "e" -x $qa_run -p $EMF_PERIOD 
           set exitstat = 1
           goto end_of_script

       endif  # If qa script failure or not
     endif      # If REPCONFIG_INV2

     if ( $?REPCONFIG_INV3_B ) then
       setenv QA_TYPE  inv3b               # Used to name the report inputs and outputs
       setenv QA_LABEL $SUBSECT           # Used to name the report inputs and outputs
       setenv REPLABEL $SUBSECT           # Used internally by Smkreport
       source $qa_run  

       # Check status of QA run to see if it worked. Give error if failed
       if ( $status != 0 ) then
           echo "ERROR: Running qa_run for $QA_TYPE" 
           $EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: Running qa_run for $QA_TYPE" -t "e" -x $qa_run -p $EMF_PERIOD 
           set exitstat = 1
           goto end_of_script

       endif  # If qa script failure or not
     endif      # If REPCONFIG_INV3
     if ( $?REPCONFIG_INV2_C ) then
       setenv QA_TYPE  inv2c              # Used to name the report inputs and outputs
       setenv QA_LABEL $SUBSECT           # Used to name the report inputs and outputs
       setenv REPLABEL $SUBSECT           # Used internally by Smkreport
       source $qa_run  

       # Check status of QA run to see if it worked. Give error if failed
       if ( $status != 0 ) then
           echo "ERROR: Running qa_run for $QA_TYPE" 
           $EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: Running qa_run for $QA_TYPE" -t "e" -x $qa_run -p $EMF_PERIOD 
           set exitstat = 1
           goto end_of_script

       endif  # If qa script failure or not
     endif      # If REPCONFIG_INV2

     if ( $?REPCONFIG_INV3_C ) then
       setenv QA_TYPE  inv3c               # Used to name the report inputs and outputs
       setenv QA_LABEL $SUBSECT           # Used to name the report inputs and outputs
       setenv REPLABEL $SUBSECT           # Used internally by Smkreport
       source $qa_run  

       # Check status of QA run to see if it worked. Give error if failed
       if ( $status != 0 ) then
           echo "ERROR: Running qa_run for $QA_TYPE" 
           $EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: Running qa_run for $QA_TYPE" -t "e" -x $qa_run -p $EMF_PERIOD 
           set exitstat = 1
           goto end_of_script

       endif  # If qa script failure or not
     endif      # If REPCONFIG_INV3
 
   endif # End if primary grid or not

   ## Run QA for part 1, gridded reports
   if ( $RUN_SMKREPORT == Y ) then
      setenv QA_TYPE  invgrid            # Used to name the report inputs and outputs
      setenv QA_LABEL $SUBSECT           # Used to name the report inputs and outputs
      setenv REPLABEL $SUBSECT           # Used internally by Smkreport
      source $qa_run

      # Check status of QA run to see if it worked. Give error if failed
      if ( $status != 0 ) then
          echo "ERROR: Running qa_run for $QA_TYPE" 
          $EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: Running qa_run for $QA_TYPE" -t "e" -x $qa_run -p $EMF_PERIOD 
          set exitstat = 1
          goto end_of_script

      ## Otherwise, register internal reports with the EMF
      else
         if ( $REGISTER_REPOUT == Y ) then
#            echo "SCRIPT NOTE: Registering internal grid-specific Smkreport reports"
#            if ( -e $REPORT14B) $EMF_CLIENT -k $EMF_JOBKEY -F $REPORT14B -T "Smkreport state-SCC-srgcode" -O "state/SCC/Srgcode Smkreport $SUBSECT"
#	 
#            echo "SCRIPT NOTE: Registering external grid-specific Smkreport reports"
#            if ( -e $REPORT7  ) $EMF_CLIENT -k $EMF_JOBKEY -F $REPORT7 -T "Smkreport state-cell (External)" \
#                                -N "rep_${SUBSECT}_invgrid_state_grid_${GRID}_$CASE (External)" \
#      	                        -O "state w/gridding Smkreport $SUBSECT (External)"
#            if ( -e $REPORT17 ) $EMF_CLIENT -k $EMF_JOBKEY -F $REPORT17 -T "Smkreport grid cell (External)" \
#                                -N "rep_${SUBSECT}_invgrid_cell_${GRID}_$CASE (External)" \
# 	                        -O "grid cell Smkreport $SUBSECT (External)"
         endif # register_repout

      endif  # If qa script failure or not

   endif

   setenv RUN_PART1 N
   
   if ($m == 0) continue # don't run Temporal or anything afterwards for the annual run

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

   setenv RUN_TEMPORAL  $run_temporal
   
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

   setenv RUN_PART2 N

   ## Run Smkmerge and Smkreport, and optionally, Smk2emis (for CAMx)

   setenv T_TYPE $M_TYPE                 # Set temporal type to type for merge
   setenv SMK_RUN_DATES $SMK_RUN_DATES_2 # Set with file name using MTYPE

   source $set_days   # Call script to set dates for run
   if ( $status != 0 ) then
	echo "ERROR: Running set_days"
        $EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: Running set_days for $MONTH" -t "e" -x $set_days -p $EMF_PERIOD
        exit (1)
   endif

   ## Determine the number of days to run
   set ndays = `cat $SMK_RUN_DATES | wc -l`
   echo "Number of dates to run for month ($m)" : $ndays

   set n = 0
   set diff = 0
   set g_stdate_sav = $g_stdate_all

   while ( $n < $ndays )

      @ n = $n + 1   

      setenv RUN_SMKREPORT $run_smkreport
      
      ## This is necessary so that the assigns file doesn't delete the Temporal files
      setenv RUN_TEMPORAL N
      
      set line = `head -$n $SMK_RUN_DATES | tail -1`
      @ diff = $line[1] - $g_stdate_sav

      setenv G_STDATE_ADVANCE $diff
#      setenv RUN_PART2 Y   # Commenting this out fixes the "Temporal files are deleted after they are created if RUNSET = Y" bug
      source $ASSIGNS_FILE               # Invoke Assigns file to set new dates

      ## Set EMF_PERIOD for this day
      setenv EMF_PERIOD $ESDATE

      # Run QA for part 2
      # C. Allen (26 May 2017): We almost always turn these reports off anyway, and this occasionally
      #   causes other users problems, so we are now turning the temporal Smkreports off for good.
#      setenv QA_TYPE  temporal           # Used to name the report inputs and outputs
#      setenv QA_LABEL $SUBSECT           # Used to name the report inputs and outputs
#      setenv REPLABEL $SUBSECT           # Used internally by Smkreport
#
#     source $qa_run  
#      if ( $status != 0 ) then
#         echo "ERROR: running QA & Smkreport for $ESDATE"
#         $EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: running QA & Smkreport for $ESDATE" -t "e" -x $qa_run -p $EMF_PERIOD
#         set exitstat = 1
#         goto end_of_script
#
#      endif

      # Run QA for part 2: one day per month to get the temporal profile assignments, if REPCONFIG_TEMP is set
      if ($?REPCONFIG_TEMP && $n == 1) then
         setenv QA_TYPE  temporal           # Used to name the report inputs and outputs
         setenv QA_LABEL $SUBSECT           # Used to name the report inputs and outputs
         setenv REPLABEL $SUBSECT           # Used internally by Smkreport
         source $qa_run  
         if ( $status != 0 ) then
            echo "ERROR: running QA & Smkreport for Temporal"
            $EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: running QA & Smkreport for Temporal" -t "e" -x $qa_run -p $EMF_PERIOD
            set exitstat = 1
            goto end_of_script
         endif
      endif
      
      setenv RUN_PART2 N

      setenv RUN_SMKMERGE  $run_smkmerge
      setenv RUN_SMK2EMIS  $run_smk2emis
      setenv RUN_M3STAT    $run_m3stat  
      
      setenv RUN_PART4 Y
      source $ASSIGNS_FILE

      # Run programs part 4
      source $smk_run        
      if ( $status != 0 ) then	    
	 echo "ERROR: Running smk_run for part 4 for $ESDATE"
         $EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: Running smk_run for part 4" -t "e" -x $smk_run -p $EMF_PERIOD
         set exitstat = 1
	 goto end_of_script

      endif

      ## C. Allen configured for multiple files per day (All of my comments contain [c])
      if ( $?AOUT && $RUN_M3STAT == Y ) then
      
         # [c] If $AOUT exists, this implies there is only one Smkmerge file per day. In this case, run m3stat once and get out
	 if ( -e $AOUT ) then
            $m3stat $AOUT
            if ( $status != 0 ) then
   	       echo "ERROR: running m3stat_chk for $ESDATE"
	       $EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: running m3stat_chk for $ESDATE" -t "e" -x $m3stat -p $EMF_PERIOD
               set exitstat = 1
	       goto end_of_script
            endif
	 # [c] Otherwise, check to see how many files there are, and run m3stat on each one
	 else
	    # [c] "aout_prefix" is prefix of $AOUT, everything before the file number
	    set aout_prefix = $PREMERGED/emis_${MM}_${SECTOR}_${ESDATE}_${GRID}_${SPC}_${CASE}
	    set num_smkmerge_files = `/bin/ls -1 $aout_prefix.*.ncf | wc -l`
	    echo "SCRIPT NOTE: Running m3stat on $num_smkmerge_files model-ready files"
	    set fc = 0
	    # [c] This allows for case where number of files is 0; in that case the while loop is skipped altogether
	    # [c] The run_m3stat script has an optional 3rd command-line parameter that I originally added to handle inline approach;
	    # [c]   this is a "file stamp" appended to the file names so that the files don't get overwritten each time
	    while ( $fc < $num_smkmerge_files )
 	       @ fc = $fc + 1
	       $m3stat $aout_prefix.$fc.ncf $fc
               if ( $status != 0 ) then
      	          echo "ERROR: running m3stat_chk for $ESDATE"
	          $EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: running m3stat_chk for $ESDATE" -t "e" -x $m3stat -p $EMF_PERIOD
                  set exitstat = 1
                  goto end_of_script
               endif
	    end # while
	 endif # -e $AOUT
      endif # $?AOUT

      setenv RUN_PART4 N

   end  # End loop over days

   unsetenv G_STDATE_ADVANCE

end  # End loop over parts

## Set day count depending on leap year and spinup period
if ( $BASE_YEAR == 2000 || $BASE_YEAR == 2004 ||  $BASE_YEAR == 2008 || $BASE_YEAR == 2012 || $BASE_YEAR == 2016  || $BASE_YEAR == 2020) then
   set alldays = 366
else
   set alldays = 365
endif

@ alldays = $alldays + $SPINUP_DURATION

## Set the expected file counts needed for Smkreport output registration
switch ( $L_TYPE ) 
  case all:
     set tmprep_expected = $alldays
  breaksw

  case week:
     set tmprep_expected = 84
     if ( $RUN_HOLIDAYS == Y ) set tmprep_expected = 100
  breaksw

  case mwdss:
     set tmprep_expected = 48
     if ( $RUN_HOLIDAYS == Y ) set tmprep_expected = 64
  breaksw

  case aveday:
     set tmprep_expected = 12
     if ( $RUN_HOLIDAYS == Y ) set tmprep_expected = 28
  breaksw
endsw

## Set the expected file counts needed for Smkmerge output registration
switch ( $M_TYPE ) 
  case all:
     set mrg_expected = $alldays
  breaksw

  case week:
     set mrg_expected = 84
     if ( $RUN_HOLIDAYS == Y ) set mrg_expected = 100
     if ( $RUN_HOLIDAYS == Y && $BASE_YEAR == 2021 ) set mrg_expected = 101
     if ( $RUN_HOLIDAYS == Y && $BASE_YEAR == 2020 ) set mrg_expected = 101
     if ( $RUN_HOLIDAYS == Y && $BASE_YEAR == 2019 ) set mrg_expected = 101
     if ( $RUN_HOLIDAYS == Y && $BASE_YEAR == 2016 ) set mrg_expected = 101
     if ( $RUN_HOLIDAYS == Y && $BASE_YEAR == 2015 ) set mrg_expected = 101
  breaksw

  case mwdss:
     set mrg_expected = 48
     if ( $RUN_HOLIDAYS == Y ) set mrg_expected = 64
     if ( $RUN_HOLIDAYS == Y && $BASE_YEAR == 2021 ) set mrg_expected = 65
     if ( $RUN_HOLIDAYS == Y && $BASE_YEAR == 2020 ) set mrg_expected = 65
     if ( $RUN_HOLIDAYS == Y && $BASE_YEAR == 2019 ) set mrg_expected = 65
     if ( $RUN_HOLIDAYS == Y && $BASE_YEAR == 2016 ) set mrg_expected = 65
     if ( $RUN_HOLIDAYS == Y && $BASE_YEAR == 2015 ) set mrg_expected = 65
  breaksw

  case aveday:
     set mrg_expected = 12
     if ( $RUN_HOLIDAYS == Y ) set mrg_expected = 28
  breaksw
endsw

## NOTE: 2004, 2009, 2010 have an extra holiday for July 4th
if ( $RUN_HOLIDAYS == Y && ($BASE_YEAR == 2004 || $BASE_YEAR == 2009 || $BASE_YEAR == 2010) ) then
  if ( $L_TYPE != all ) @ tmprep_expected = $tmprep_expected + 1
  if ( $M_TYPE != all ) @ mrg_expected = $mrg_expected + 1
endif

## Determine if all temporal report files are available to import
##   Use state temporal report as sample.
set tmprepcnt = `/bin/ls -1 $REPOUT/temporal/rep_${SECTOR}_*_${CASE}_*_temporal_state.txt | wc -l`

## Determine if all smkmerge report files are available to import
set mrgrepcnt = `/bin/ls -1 $REPOUT/smkmerge/$SECTOR/rep_mole_${SECTOR}_*_${GRID}_*.txt | wc -l`

## Determine if all smkmerge model-ready files are available to import
set mrgcnt = `/bin/ls -1 $PREMERGED/emis_mole_${SECTOR}_*_${GRID}_*.ncf* | wc -l`

# Register Smkreport temporal outputs with EMF for current month
# Make sure Smkreport was run for any days thus far (don't want to re-register data if it wasn't run)
if ( $REGISTER_REPOUT == Y ) then
   if ( 0 ) then # if ( $tmprepcnt == $tmprep_expected ) then
      echo "SCRIPT NOTE: Registering external Smkreport temporal reports"
      $EMF_CLIENT -k $EMF_JOBKEY \
           -D $REPOUT/temporal \
           -P "rep_${SECTOR}_*_${CASE}_*_temporal_county.txt" \
           -T "Smkreport county daily (External Multifile)" \
           -N "rep_${SECTOR}_temporal_county_$CASE" \
           -O "county daily Smkreport $SECTOR (External)"

      $EMF_CLIENT -k $EMF_JOBKEY \
           -D $REPOUT/temporal \
           -P "rep_${SECTOR}_*_${CASE}_*_temporal_county_moncode.txt" \
           -T "Smkreport county-moncode daily (External Multifile)" \
           -N "rep_${SECTOR}_temporal_county_moncode_$CASE" \
           -O "county/moncode daily Smkreport $SECTOR (External)"

      $EMF_CLIENT -k $EMF_JOBKEY \
           -D $REPOUT/temporal \
           -P "rep_${SECTOR}_*_${CASE}_*_temporal_state.txt" \
           -T "Smkreport state daily (External Multifile)" \
           -N "rep_${SECTOR}_temporal_state_$CASE" \
           -O "state daily Smkreport $SECTOR (External)"

      $EMF_CLIENT -k $EMF_JOBKEY \
           -D $REPOUT/temporal \
           -P "rep_${SECTOR}_*_${CASE}_*_temporal_state_scc.txt" \
           -T "Smkreport state-SCC daily (External Multifile)" \
           -N "rep_${SECTOR}_temporal_state_scc_$CASE" \
           -O "state/SCC daily Smkreport $SECTOR (External)"
   else
      echo "SCRIPT NOTE: Registration of Smkreport temporal reports is currently automatically turned off"
      echo "             because these reports are no longer being generated."
#      echo "SCRIPT NOTE: Not registering Smkreport temporal reports because only found"
#      echo "             $tmprepcnt state temporal reports out of expected $tmprep_expected" 

   endif

   ## Register Smkmerge reports
   if ( $mrgcnt == $mrg_expected ) then
#      echo "SCRIPT NOTE: Registering Smkmerge reports"
#      $EMF_CLIENT -k $EMF_JOBKEY \
#           -D $REPOUT/smkmerge/$SECTOR \
#           -P "rep_mole_${SECTOR}_*_${GRID}_*.txt" \
#           -T "Smkmerge Report state (External Multifile)" \
#           -N "Smkmerge $SECTOR $CASE state reports" \
#           -O "state Smkmerge $SECTOR (External)"

      ## Run python annual report script, if all Smkmerge reports are there,
      #  and if RUN_PYTHON_ANNUAL parameter is Y (default is N)
      if ( ! $?RUN_PYTHON_ANNUAL ) then
         setenv RUN_PYTHON_ANNUAL N
      endif	 
      if ( $RUN_PYTHON_ANNUAL != N  && $?SECTORLIST ) then

         echo "SCRIPT NOTE: Running python annual report script"

         ## Get root of the merge dates directory from smk merge date file
         if ( $?MRGDATE_FILES) then
	    ## get path from first merge date file
	    echo "Getting merge dates path from MRGDATE_FILES"
	    set mrgdate_root = `$path_parser $MRGDATE_FILES`
         else 
	    ## no merge dates file, so get from parent case's inputs/mrggrid
	    echo "MRGDATE_FILES not defined, so merge dates path defaults to $INVDIR/mrggrid"
     	    set mrgdate_root = $INVDIR/mrggrid
         endif

         if (! -e $REPOUT/annual_report) mkdir -p $REPOUT/annual_report

         # ANNUAL_REPORT_OUTPUT determines which reports annual_report.py
         #   creates, either "state", "county", or "both". Default is "both".
         if (! $?ANNUAL_REPORT_OUTPUT) setenv ANNUAL_REPORT_OUTPUT "both"
	 
         if ($MRG_REPCNY_YN == "Y") then # county Smkmerge reports


           $python_annl -R $REPOUT/smkmerge -D $mrgdate_root \
                        -O $REPOUT/annual_report \
	    	        -o annual_${CASE}_${SECTOR}_${GRID}_${SPC} \
			-r county -a $ANNUAL_REPORT_OUTPUT \
	                $SCRIPTS/annual_report/parameter_file_$SPC.txt
			
	 else # state Smkmerge reports
	  
	   # if Smkmerge creates state reports, ANNUAL_REPORT_OUTPUT must = "state"
	   setenv ANNUAL_REPORT_OUTPUT "state"

           $python_annl -R $REPOUT/smkmerge -D $mrgdate_root \
                        -O $REPOUT/annual_report \
	  	        -o annual_${CASE}_${SECTOR}_${GRID}_${SPC} \
			-r state -a $ANNUAL_REPORT_OUTPUT \
	                $SCRIPTS/annual_report/parameter_file_$SPC.txt

         endif

         if ($status == 0) then
            ## Register the output w/ the EMF
	    if ($ANNUAL_REPORT_OUTPUT == "both" || $ANNUAL_REPORT_OUTPUT == "state") then
#              echo "SCRIPT NOTE: Registering smkmerge state totals"
#              $EMF_CLIENT -k $EMF_JOBKEY \
#   	        -F  $REPOUT/annual_report/annual_${CASE}_${SECTOR}_${GRID}_${SPC}_state_emf.csv \
#                -T "Smkmerge report state annual summary (CSV)" -O "state totals ${namelabel}"    
	    endif
	    if ($ANNUAL_REPORT_OUTPUT == "both" || $ANNUAL_REPORT_OUTPUT == "county") then
#              echo "SCRIPT NOTE: Registering smkmerge county totals"
#              $EMF_CLIENT -k $EMF_JOBKEY \
#   	        -F  $REPOUT/annual_report/annual_${CASE}_${SECTOR}_${GRID}_${SPC}_county_emf.csv \
#                -T "Smkmerge report county annual summary (CSV)" -O "county totals ${namelabel}"    
	    endif
         endif

      endif

   else
      echo "SCRIPT NOTE: Not registering Smkmerge reports because only found"
      echo "             $mrgcnt Smkmerge reports out of expected $mrg_expected" 

   endif

endif

# Register Smkmerge AQM-ready outputs with EMF for current month
#if ( $REGISTER_AQMOUT == Y ) then
#   if ( $mrgcnt == $mrg_expected ) then
#      echo "SCRIPT NOTE: Registering Smkmerge outputs"
#      $EMF_CLIENT -k $EMF_JOBKEY \
#           -D $PREMERGED \
#           -P "emis_mole_${SECTOR}_*_${GRID}_*.ncf" \
#           -T "CMAQ Model Ready Emissions: Sector-specific (External)" \
#           -N "Model-ready CMAQ $SECTOR $CASE" \
#           -O "AQM-ready data $SECTOR"
#   else
#      echo "SCRIPT NOTE: Not registering Smkmerge model-ready files because only found"
#      echo "             $mrgcnt Smkmerge outputs out of expected $mrg_expected" 
#
#   endif
#endif

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

## Register log analyzer output
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