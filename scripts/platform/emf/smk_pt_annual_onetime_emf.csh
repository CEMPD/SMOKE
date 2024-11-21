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

if ($?GCNTL) then  # if GCNTL exists, run Cntlmat; otherwise don't run Cntlmat
   setenv RUN_CNTLMAT Y        #  run control matrix program
   setenv MRG_CTLMAT_MULT $MRG_SOURCE  # 'A' merges with area projection matrix produced from CNTLMAT
   setenv MRG_REPCTL_YN   Y            # Y Report control totals 
   setenv CUSTOM_GCNTL    Y
else
   setenv RUN_CNTLMAT N
endif

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
if (! $?REGISTER_REPOUT) then
   setenv REGISTER_REPOUT    Y       # Imports Smkreport and Smkmerge reports into EMF
endif
if (! $?REGISTER_AQMOUT) then
   setenv REGISTER_AQMOUT    Y       # Imports Smkmerge I/O API outputs in EMF
endif
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
set cntl_run     = $SCRIPTS/run/cntl_run_v6.csh      # Runs CNTLMAT and GRWINVEN
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
echo $cntl_run >> $LOGS/helper_scripts_list$suffix
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

## Run QA for part 1, non-gridded reports
## Run only if Smkreport still set to Y (may have been turned off by RUNSET)
## Run only for primary grid
#if ( $RUN_SMKREPORT == Y && $GRID == $GRID_1 ) then
# Primary grid no longer set with Regions 4/14/11 J. Beidler
if ( $RUN_SMKREPORT == Y ) then

   if ( $?CASE_CON ) then
      setenv QA_TYPE proj
   else
      setenv QA_TYPE  inv                # Used to name the report inputs and outputs
   endif

   setenv QA_LABEL $SUBSECT           # Used to name the report inputs and outputs
   setenv REPLABEL $SUBSECT           # Used internally by Smkreport
   source $qa_run

   # Check status of QA run to see if it worked. Give error if failed
   if ( $status != 0 ) then
       echo "ERROR: Running qa_run for $QA_TYPE" 
       $EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: Running qa_run for $QA_TYPE" -t "e" -x $qa_run -p $EMF_PERIOD 
       set exitstat = 1
       goto end_of_script

   ## Otherwise, register internal and external reports with the EMF
   else if ( $REGISTER_REPOUT == Y ) then

#      echo "SCRIPT NOTE: Registering internal inventory Smkreport reports"
#      if ( -e $REPORT1  ) $EMF_CLIENT -k $EMF_JOBKEY -F $REPORT1 -T "Smkreport state" -O "state Smkreport $SECTOR"
#      if ( -e $REPORT2  ) $EMF_CLIENT -k $EMF_JOBKEY -F $REPORT2 -T "Smkreport state-SCC" -O "state/SCC Smkreport $SECTOR"
#      if ( -e $REPORT4  ) $EMF_CLIENT -k $EMF_JOBKEY -F $REPORT4 -T "Smkreport county" -O "county Smkreport $SECTOR"
#      if ( -e $REPORT16 ) $EMF_CLIENT -k $EMF_JOBKEY -F $REPORT16 -T "Smkreport state-SCC-spec_profile" -O "state/SCC/PM25 profile Smkreport $SECTOR"
#
#      echo "SCRIPT NOTE: Registering external Smkreport reports"
#      if ( -e $REPORT8 ) $EMF_CLIENT -k $EMF_JOBKEY -F $REPORT8 -T "Smkreport county-SCC (External)" \
#                  -N "rep_${SECTOR}_inv_county_scc_$CASE (External)" \
#                  -O "county/SCC Smkreport $SECTOR (External)"
#      if ( -e $REPORT9  ) $EMF_CLIENT -k $EMF_JOBKEY -F $REPORT9 -T "Smkreport state-SIC (External)" \
#                  -N "rep_${SECTOR}_inv_state_sic_$CASE (External)" \
#                  -O "state/SIC Smkreport $SECTOR (External)"
#      if ( -e $REPORT19 ) $EMF_CLIENT -k $EMF_JOBKEY -F $REPORT19 -T "Smkreport state-NAICS (External)" \
#                  -N "rep_${SECTOR}_inv_state_naics_$CASE (External)" \
#	          -O "state/NAICS Smkreport $SECTOR (External)"
#      if ( -e $REPORT20 ) $EMF_CLIENT -k $EMF_JOBKEY -F $REPORT20 -T "Smkreport plant-county-ORIS (External)" \
#                  -N "rep_${SECTOR}_inv_plant_oris_$CASE (External)" \
#                  -O "county/plant/ORIS Smkreport $SECTOR (External)"
#      if ( -e $REPORT21 ) $EMF_CLIENT -k $EMF_JOBKEY -F $REPORT21 -T "Smkreport state-MACT (External)" \
#                  -N "rep_${SECTOR}_inv_state_mact_$CASE (External)" \
#	          -O "state/MACT Smkreport $SECTOR (External)"

   endif  # If qa script failure or not

   if ( $?SMK_PART1ONLY ) then
      if ( $SMK_PART1ONLY == Y ) then
         goto end_of_script
      endif
   endif

   ## Continue QA, but for state/SCC VOC profiles (PM2.5 profiles were in the first run, and can't 
   #      do both in same Smkreport run)
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

endif  # End if primary grid or not

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
   else if ( $REGISTER_REPOUT == Y ) then
#      echo "SCRIPT NOTE: Registering external grid-specific Smkreport reports"
#      if ( -e $REPORT7  ) $EMF_CLIENT -k $EMF_JOBKEY -F $REPORT7 -T "Smkreport state-cell (External)" \
#                          -N "rep_${SECTOR}_invgrid_state_grid_${GRID}_$CASE (External)" \
#                          -O "state w/gridding Smkreport $SECTOR (External)"
#      if ( -e $REPORT17 ) $EMF_CLIENT -k $EMF_JOBKEY -F $REPORT17 -T "Smkreport grid cell (External)" \
#                          -N "rep_${SECTOR}_invgrid_cell_${GRID}_$CASE (External)" \
#	                  -O "grid cell Smkreport $SECTOR (External)"
#      if ( -e $REPORT11 ) $EMF_CLIENT -k $EMF_JOBKEY -F $REPORT11 -T "Smkreport plant-cell (External)" \
#                          -N "rep_${SECTOR}_invgrid_plant_cell_${GRID}_$CASE (External)" \
#                          -O "plant/cell Smkreport $SECTOR (External)"
#      if ( -e $REPORT23 ) $EMF_CLIENT -k $EMF_JOBKEY -F $REPORT23 -T "Smkreport grid cell-county (External)" \
#                          -N "rep_${SECTOR}_invgrid_cell_county_${GRID}_$CASE (External)" \
#                          -O "grid cell w/county Smkreport $SECTOR (External)"
   endif  # If qa script failure or not

endif

## (C. Allen, 21 Jun 2013) New parameter USE_FF10_DAILY_POINT to be used with FF10 daily point inventories.
#  If Y, then after Smkreport, run Smkinven once for the daily/hourly processing, since we only need to run that
#  once for the year now rather than separately for every month. Defaults to N for backwards compatibility.

if ($USE_FF10_DAILY_POINT == Y) then

   ## This script is linear with no loops, so we can turn Spcmat and Grdmat off at this point without worrying about turning them
   #  back on later.
   setenv RUN_SPCMAT N
   setenv RUN_GRDMAT N

   ## Set DAY_SPECIFIC_YN and HOUR_SPECIFIC_YN back to their original settings.
   setenv DAY_SPECIFIC_YN $day_specific_yn_save
   setenv HOUR_SPECIFIC_YN $hour_specific_yn_save

   ## Don't import annual inventory, that's already done
   setenv IMPORT_AVEINV_YN  N           
 
   source $ASSIGNS_FILE                  # Invoke Assigns file to set for current month

   # Rename original smkinven log so that it doesn't get overwritten
   mv -f $LOGS/smkinven_${SRCABBR}_$CASE.log $LOGS/smkinven_${SRCABBR}_annual_$CASE.log

   ## Setup environment variables for using in the PTDAY list files
   #  FF10 daily point should always be internal, one for the full year, so no need
   #  for the old NAMEBREAK_DAILY nonsense. Instead, anything input whose e.v. starts
   #  with EMISDAY, such as EMISDAY_A, goes in the ptday .lst file.

   if ( $?MULTIFILE_NAMEBREAK ) unsetenv MULTIFILE_NAMEBREAK
   $combine_data EMISDAY $DAYLISTFILE list
   if ( $status != 0 ) then
       echo "ERROR: Could not run combine_data.csh for daily list file in $EMF_JOBNAME"
       $EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: Could not run combine_data for daily list file in $EMF_JOBNAME" -t "e" -p $EMF_PERIOD
       exit ( 1 )
   endif

   setenv PTDAY $DAYLISTFILE  #  Set the point daily-specific inventory list file to the generic list file

   ## Now, this is where it gets tricky. Used to be, for the CEM data, we'd have separate pthour.lst files
   #  for each month. Now we will need to put all PTHOUR files for the full year into a single pthour.lst.
   #  We still need NAMEBREAK_HOURLY for this.

   ## Create the month-specific list file of hourly data for import to Smkinven
   if ( $HOUR_SPECIFIC_YN == Y ) then
      if ( $?NAMEBREAK_HOURLY ) then
          setenv MULTIFILE_NAMEBREAK $NAMEBREAK_HOURLY
      endif
      $combine_data EMISHOUR $HOURLISTFILE list 01 Hourly  # the 01 for the month is sort of a dummy value
      if ( $status != 0 ) then
          echo "ERROR: Could not run combine_data.csh for hourly list file in $EMF_JOBNAME"
          $EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: Could not run combine_data for hourly list file in $EMF_JOBNAME" -t "e" -p $EMF_PERIOD
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

endif # USE_FF10_DAILY_POINT = Y

## (J. Beidler 1 Nov 2016) New parameter USE_FF10_HOURLY_POINT to be used with FF10 hourly point inventories.
#  If Y, then after Smkreport, run Smkinven once for the hourly processing, since we only need to run that
#  once for the year now rather than separately for every month. Defaults to N for backwards compatibility.

if ($USE_FF10_HOURLY_POINT == Y) then

   ## This script is linear with no loops, so we can turn Spcmat and Grdmat off at this point without worrying about turning them
   #  back on later.
   setenv RUN_SPCMAT N
   setenv RUN_GRDMAT N

   ## Set HOUR_SPECIFIC_YN back to their original settings.
   setenv HOUR_SPECIFIC_YN $hour_specific_yn_save

   ## Don't import annual inventory, that's already done
   setenv IMPORT_AVEINV_YN  N           
 
   source $ASSIGNS_FILE                  # Invoke Assigns file to set for current month

   # Rename original smkinven log so that it doesn't get overwritten
   mv -f $LOGS/smkinven_${SRCABBR}_$CASE.log $LOGS/smkinven_${SRCABBR}_annual_$CASE.log

   setenv PTHOUR $HOURLISTFILE

   ## Now, this is where it gets tricky. Used to be, for the CEM data, we'd have separate pthour.lst files
   #  for each month. Now we will need to put all PTHOUR files for the full year into a single pthour.lst.
   #  We still need NAMEBREAK_HOURLY for this.

   ## Create the month-specific list file of hourly data for import to Smkinven
   if ( $HOUR_SPECIFIC_YN == Y ) then
      if ( $?NAMEBREAK_HOURLY ) then
          setenv MULTIFILE_NAMEBREAK $NAMEBREAK_HOURLY
      endif
      $combine_data EMISHOUR $HOURLISTFILE list 01 Hourly  # the 01 for the month is sort of a dummy value
      if ( $status != 0 ) then
          echo "ERROR: Could not run combine_data.csh for hourly list file in $EMF_JOBNAME"
          $EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: Could not run combine_data for hourly list file in $EMF_JOBNAME" -t "e" -p $EMF_PERIOD
          exit ( 1 )
      endif
      setenv PTHOUR $HOURLISTFILE  #  Set the point hour-specific inventory list file to the generic list file
   endif

   ## Run Smkinven to import day-specific and hour-specific inventories
   source $smk_run    # Run program
   if ( $status != 0 ) then
       echo "ERROR: Running Smkinven for hourly data import in $EMF_PERIOD"
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

endif # USE_FF10_HOURLY_POINT = Y

setenv RUN_PART1 N

if ($RUN_CNTLMAT == Y) then
   # PART1B runs cntlmat and grwinven
   setenv RUN_PART1B Y

   source $ASSIGNS_FILE   # Invoke Assigns file

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
