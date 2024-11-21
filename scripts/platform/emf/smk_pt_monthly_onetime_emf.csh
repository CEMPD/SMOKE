#!/bin/tcsh -f

# Version @(#)$Id$
# Path    $Source$
# Date    $Date$

# This script sets up needed environment variables for running SMOKE point
# sources using annual emissions for only the "onetime" steps (Smkinven,
# Spcmat, Grdmat, Elevpoint and inventory QA).
#
# This script is intended to be used with the EMF
# source emissions in SMOKE for the EPA 2002 modeling platform, and 
# calls the scripts that runs the SMOKE programs. 
#
# Script created by : M. Houyoux, Environmental Protection Agency
# October, 2007
# Modified to support concatenation for GSCNV and GSPRO_COMBO, M. Houyoux Dec, 2008
#
#*********************************************************************

## log w/ EMF server that script is running
$EMF_CLIENT -k $EMF_JOBKEY -s "Running" 

# set source category
setenv SMK_SOURCE P           # source category to process
setenv MRG_SOURCE $SMK_SOURCE # source category to merge

## time-independent programs (except Smkinven run monthly for CEM import)
setenv RUN_SMKINVEN  Y        #  run inventory import program
setenv RUN_SPCMAT    Y        #  run speciation matrix program

## run GRDMAT only if not PART1ONLY
if ( $?SMK_PART1ONLY ) then
   if ( $SMK_PART1ONLY == Y ) then
      setenv RUN_GRDMAT    N        #  run gridding matrix program
   else
      setenv RUN_GRDMAT    Y        #  run gridding matrix program
   endif
else
   setenv RUN_GRDMAT    Y        #  run gridding matrix program
endif

## Set whether or not to run Elevpoint based on other settings
if ( $SMK_SPECELEV_YN == Y ) then
   setenv RUN_ELEVPOINT Y        #  run elevated source selection program (once only since PELVCONFIG not based on emissions)
endif

## (C. Allen) ELEVPOINT_DAILY means Elevpoint is run once for every day. This is used for ptfire.
#  If ELEVPOINT_DAILY = Y, this means we do NOT want to run Elevpoint "annually", so set RUN_ELEVPOINT to N.
#  If ELEVPOINT_DAILY = N or is undefined, leave RUN_ELEVPOINT as is (set above based on SMK_SPECELEV_YN).
if ( $?ELEVPOINT_DAILY ) then
   if ( $ELEVPOINT_DAILY == Y ) then
      setenv RUN_ELEVPOINT N
   endif
endif

## (C. Allen, 21 Jun 2013) New parameter USE_FF10_DAILY_POINT to be used with FF10 daily point inventories.
#  If Y, then after Smkreport, run Smkinven once for the daily/hourly processing, since we only need to run that
#  once for the year now rather than separately for every month. Defaults to N for backwards compatibility.
if (! $?USE_FF10_DAILY_POINT) then
   setenv USE_FF10_DAILY_POINT N
endif

## (J. Beidler 1 Nov 2016) New parameter USE_FF10_HOURLY_POINT to be used with FF10 hourly point inventories.
#  If Y, then after Smkreport, run Smkinven once for the daily/hourly processing, since we only need to run that
#  once for the year now rather than separately for every month. Defaults to N for backwards compatibility.
if (! $?USE_FF10_HOURLY_POINT) then
   setenv USE_FF10_HOURLY_POINT N
endif

## (C. Allen, 21 Jun 2013) In concert with USE_FF10_DAILY_POINT, preserve DAY_SPECIFIC_YN and HOUR_SPECIFIC_YN
#  settings for later. Set to N for now for the first Smkinven run.
#  I haven't yet tested whether SMOKE can handle FF10 daily and FF10 annual in a single Smkinven run.
if (! $?DAY_SPECIFIC_YN) setenv DAY_SPECIFIC_YN N
if (! $?HOUR_SPECIFIC_YN) setenv HOUR_SPECIFIC_YN N
set day_specific_yn_save = $DAY_SPECIFIC_YN
set hour_specific_yn_save = $HOUR_SPECIFIC_YN
setenv DAY_SPECIFIC_YN N
setenv HOUR_SPECIFIC_YN N

## quality assurance

# allow the user to turn these things off in his/her case
#if (! $?REGISTER_REPOUT) then
   setenv REGISTER_REPOUT    N       # Imports Smkreport and Smkmerge reports into EMF
#endif
#if (! $?REGISTER_AQMOUT) then
   setenv REGISTER_AQMOUT    N       # Imports Smkmerge I/O API outputs in EMF
#endif
if (! $?RUN_SMKREPORT) then
#  if not set already, set it from the variable.
   setenv RUN_SMKREPORT Y # Y runs reporting for state reports
else
#  this allows for someone to override the run report setting in their case
   set run_smkreport = $RUN_SMKREPORT
   if ($RUN_SMKREPORT == 'N') then
       setenv REGISTER_REPOUT    N
#      don't register the reports if you are not creating them
   endif
endif
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
   case 3:
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

## Set EVs accordingly depending on whether we want inline approach or not
## If INLINE_MODE isn't set, default to "off"
if (! $?INLINE_MODE ) then
  setenv INLINE_MODE      off
endif

if ( $INLINE_MODE == only || $INLINE_MODE == both ) then
   setenv SMK_ELEV_METHOD 2
else
   # Don't overwrite SMK_ELEV_METHOD if it's already defined, for backwards compatibility
   if (! $?SMK_ELEV_METHOD ) setenv SMK_ELEV_METHOD 1
endif

# Get the first two options for the grid abbreviation and I/O API grid
setenv GRID "$argv[1]"
setenv IOAPI_GRIDNAME_1 "$argv[2]"

## source the ASSIGN file
source $ASSIGNS_FILE

## List of all the helper scripts that are run in this script
set emf_cleanup  = $SCRIPTS/run/emf_cleanup.csh
set set_months   = $SCRIPTS/run/set_months_v4.csh
set timetracker  = $SCRIPTS/run/timetracker_v2.csh
set combine_data = $SCRIPTS/run/combine_data_v6.csh
set smk_run      = $SCRIPTS/run/smk_run_v9.csh
set qa_run       = $SCRIPTS/run/qa_run_v10.csh
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

# Setup GSPRO, GSREF, GSCNV, and GSPRO_COMBO file names
# Set the output file name dates using the most current file matching the prefix
set tmplist = (`env| grep GSPROTMP | cut -d\= -f2` )
if ($#tmplist > 0 ) then
   set gspro_date = `/bin/ls -1 -t $tmplist | head -1 | sed 's/_/\xa/g' | tail -2 | head -1`
   setenv GSPRO $GE_DAT/speciation/$CASE/gspro_${SECTOR}_${SPC}_${CASE}_${gspro_date}.txt
endif

set tmplist = (`env| grep GSREFTMP | cut -d\= -f2` )
if ($#tmplist > 0 ) then
   set gsref_date = `/bin/ls -1 -t $tmplist | head -1 | sed 's/_/\xa/g' | tail -2 | head -1`
   setenv GSREF $GE_DAT/speciation/$CASE/gsref_${SECTOR}_${SPC}_${CASE}_${gsref_date}.txt
endif

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
   else
      setenv SUBSECT ${SECTOR}_${monname[$m]} # set variable for input/output names
      setenv MONTH   ${monname[$m]}           # set variable for month name
   endif
   setenv SRCABBR $SUBSECT                  # set abbreviation for naming log files
   setenv EISECTOR $SUBSECT

   ## Run Smkinven (annual inventory only), Grdmat, and Spcmat
   #
   setenv RUN_PART1 Y
   source $ASSIGNS_FILE                # Invoke Assigns file to set new dates

   ## Set the EMF_PERIOD to the year
   setenv EMF_PERIOD "${MONTH}_${YEAR}"

   # Call EMF Client for current period 
   $EMF_CLIENT -k $EMF_JOBKEY -m "Running SMOKE one-time steps for month $MONTH" -p $EMF_PERIOD   ## log w/ EMF server

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
   setenv PTINV $LISTFILE  #  Set the area inventory list file to the generic list file

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

   ## Run QA for part 1, non-gridded reports
   ## Run only if Smkreport still set to Y (may have been turned off by RUNSET)
   ## Run only for primary grid
   #if ( $RUN_SMKREPORT == Y && $GRID == $GRID_1 ) then
   # Primary grid no longer set with Regions 4/14/11 J. Beidler
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
            echo "SCRIPT NOTE: Registering internal inventory Smkreport reports"
            echo "REGISTER_REPOUT = " $REGISTER_REPOUT
            if ( -e $REPORT1  ) $EMF_CLIENT -k $EMF_JOBKEY -F $REPORT1 -T "Smkreport state" -O "state Smkreport $SUBSECT"
            if ( -e $REPORT2  ) $EMF_CLIENT -k $EMF_JOBKEY -F $REPORT2 -T "Smkreport state-SCC" -O "state/SCC Smkreport $SUBSECT"
            if ( -e $REPORT4  ) $EMF_CLIENT -k $EMF_JOBKEY -F $REPORT4 -T "Smkreport county" -O "county Smkreport $SUBSECT"
            if ( -e $REPORT16 ) $EMF_CLIENT -k $EMF_JOBKEY -F $REPORT16 -T "Smkreport state-SCC-spec_profile" -O "state/SCC/PM25 profile Smkreport $SUBSECT"

            echo "SCRIPT NOTE: Registering external Smkreport reports"
            if ( -e $REPORT8 ) $EMF_CLIENT -k $EMF_JOBKEY -F $REPORT8 -T "Smkreport county-SCC (External)" \
                        -N "rep_${SUBSECT}_inv_county_scc_$CASE (External)" \
                        -O "county/SCC Smkreport $SUBSECT (External)"
            if ( -e $REPORT9  ) $EMF_CLIENT -k $EMF_JOBKEY -F $REPORT9 -T "Smkreport state-SIC (External)" \
                        -N "rep_${SUBSECT}_inv_state_sic_$CASE (External)" \
	                -O "state/SIC Smkreport $SUBSECT (External)"
            if ( -e $REPORT19 ) $EMF_CLIENT -k $EMF_JOBKEY -F $REPORT19 -T "Smkreport state-NAICS (External)" \
                        -N "rep_${SUBSECT}_inv_state_naics_$CASE (External)" \
	                -O "state/NAICS Smkreport $SUBSECT (External)"
            if ( -e $REPORT21 ) $EMF_CLIENT -k $EMF_JOBKEY -F $REPORT21 -T "Smkreport state-MACT (External)" \
                        -N "rep_${SUBSECT}_inv_state_mact_$CASE (External)" \
                        -O "state/MACT Smkreport $SUBSECT (External)"
			
	 endif # register_repout

      endif # if qa script failure or not

      # Continue QA, but for state/SCC VOC profiles (PM2.5 profiles were in the first run, and can't do both in same Smkreport run)
      #  This section was edited on 7 Feb 2019 in support of emissions modeling workgroup needs
      #  If REPCONFIG_INV2 exists, run Smkreport again
      #  If REPCONFIG_INV3 exists, run Smkreport once more
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
   
   endif # run smkreport

   ## Run QA for part 1, gridded reports (run for all grids - i.e., do not check that $GRID = $GRID1)
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
            echo "SCRIPT NOTE: Registering internal grid-specific Smkreport reports"
            if ( -e $REPORT14B) $EMF_CLIENT -k $EMF_JOBKEY -F $REPORT14B -T "Smkreport state-SCC-srgcode" -O "state/SCC/Srgcode Smkreport $SUBSECT"
	 
            echo "SCRIPT NOTE: Registering external grid-specific Smkreport reports"
            if ( -e $REPORT7  ) $EMF_CLIENT -k $EMF_JOBKEY -F $REPORT7 -T "Smkreport state-cell (External)" \
                                -N "rep_${SUBSECT}_invgrid_state_grid_${GRID}_$CASE (External)" \
      	                        -O "state w/gridding Smkreport $SUBSECT (External)"
            if ( -e $REPORT17 ) $EMF_CLIENT -k $EMF_JOBKEY -F $REPORT17 -T "Smkreport grid cell (External)" \
                                -N "rep_${SUBSECT}_invgrid_cell_${GRID}_$CASE (External)" \
 	                        -O "grid cell Smkreport $SUBSECT (External)"
         endif # register_repout

      endif  # If qa script failure or not

   endif

   setenv RUN_PART1 N

## New parameters USE_FF10_DAILY_POINT and USE_FF10_HOURLY_POINT are not relevant for monthly point processing. 
#  Keep Smkinven daily/hourly processing in the daily job.

   ## Run Elevpoint (only once, since elevated source criteria is not emissions-based)
   #
   setenv RUN_PART3 Y
   source $ASSIGNS_FILE                 # Invoke Assigns file to set new dates
   
   # Call EMF Client for elevpoint
   $EMF_CLIENT -k $EMF_JOBKEY -m "Running Elevpoint" -p $EMF_PERIOD   ## log w/ EMF server

   source $smk_run 
   if ( $status != 0 ) then
       echo "ERROR: Running smk_run for Elevpoint"
       $EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: Running smk_run for Elevpoint" -t "e" -x $smk_run -p $EMF_PERIOD 
       set exitstat = 1
       goto end_of_script
   endif

   setenv RUN_PART3 N

end  # End loop over months

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
