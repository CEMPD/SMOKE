#!/bin/tcsh -f

# Version @(#)$Id$
# Path    $Source$
# Date    $Date$

# This script merges the three onroad components together, RPP/RPD/RPV/RPH, for
# SMOKE-MOVES 2014. Based loosely on smk_sectormerge_emf.csh.
#
# This script is intended to be used with the EMF
# source emissions in SMOKE for the EPA 2011 v6 modeling platform, and 
# calls the scripts that runs the SMOKE programs. 
#
# The SECTORLIST does NOT need to be defined for this script to work.
# The modes RPP/RPD/RPV/RPH are hard-wired into the script. The sector name
# "onroad" is not. Daily emissions are assumed. This kind of hardwiring
# can be changed if necessary, with the help of a SECTORLIST-like text file.
#
# Original script created by : M. Houyoux, Environmental Protection Agency
# February, 2008
#
# Modified by: C. Allen, CSC (May 2011, Jul 2012)
# Modified by: J. Beidler, CSC (Sep 2014)
# C.Allen added support for RPHO on 13 Dec 2019

## log w/ EMF server that script is running
$EMF_CLIENT -k $EMF_JOBKEY -s "Running" 

setenv SRCABBR $SECTOR

set run_m3stat    = Y      # Y runs the m3stat program on the Smkmerge outputs

setenv PROMPTFLAG         N       # Y (never set to Y for batch processing)
setenv AUTO_DELETE        Y       # Y deletes SMOKE I/O API output files (recommended)
setenv AUTO_DELETE_LOG    Y       # Y automatically deletes logs without asking

if ( ! $?DEBUGMODE ) then
   setenv DEBUGMODE          N       # Y changes script to use debugger
endif
setenv DEBUG_EXE    totalview     # Sets the debugger to use when DEBUGMODE = Y

set debug_script = N
##############################################################################

switch ( $#argv )
   case 0:
   case 1:
   case 2:
      echo "SCRIPT ERROR: Script requires 1 argument for a grid name"
      echo "              and the -m or -q option with 3 settings."
      echo " "
      echo "  This script expects to be called using one of the following argument lists:"
      echo "     <grid abbrv> -m <monthlist> <spinup> <label>"
      echo "     <grid abbrv> -q <quarters> <spinup> <label>"
      echo " "
      echo "  You can either use one approach or the other (differing by the -m or -q options)."
      echo " "
      echo "  In the above list, the arguments are defined as follows:"
      echo "     <grid abbrv>       : Grid abbreviation (e.g., 36US1)"
      echo "     <monthlist>        : list of months to run when using the -m option"
      echo "     <quarters>         : list of quarters to run when using the -q option"
      echo "     <spinup>           : set to number of days between 1 and 20 to run a spinup" 
      echo "                          period (value sets number of days), and N otherwise"
      echo "     <label>            : label to put on TIMELOG file and helper-scripts list"
      echo " "
      echo "  Examples:"
      echo "     <script name> 36US1 -m '1 2 3' 0 jan-sep"
      echo "              This example runs the script for Jan, Feb, & Mar"
      echo "              for the 36US1 grid, with no spinup days and"
      echo "              gives a label to the TIMELOG file of jan-sep." 
      echo " "
      echo "     <script name> 12EUS1 -q '2 3' 10 apr-sep:"
      echo "               This example runs the script for the 2nd & 3rd quarters,"
      echo "               for the 12EUS1 grid, including 10 spin-up days, and gives"
      echo "               a label to the TIMELOG file of apr-sep."
      $EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: smoke script did not receive more than 2 arguments" -t "e"
      exit( 1 )
endsw

# Get the first two options for the grid abbreviation and I/O API grid
setenv GRID "$argv[1]"

set exitstat = 0

## Initialize zip option to N, if not set
if ( ! $?GZIP_OUTPUTS ) then
   setenv GZIP_OUTPUTS N
endif

## source the ASSIGNS file to get EVs: MM, OUTPUT, SCRIPTS, SMK_BIN
##   ASSIGNS, and IOAPIDIR, LOGS
source $ASSIGNS_FILE

## List of all the helper scripts that are run in this script
set emf_cleanup   = $SCRIPTS/run/emf_cleanup.csh
set set_months    = $SCRIPTS/run/set_months_v4.csh
set timetracker   = $SCRIPTS/run/timetracker_v2.csh
set smk_run       = $SCRIPTS/run/smk_run_v9.csh
set m3stat        = $SCRIPTS/run/m3stat_chk_v6.csh
set set_days      = $SCRIPTS/run/set_days_v5.csh
set log_analyzer  = $SCRIPTS/log_analyzer/log_analyzer.py
set msg_list      = $SCRIPTS/log_analyzer/known_messages.txt
set domain_totals = $SCRIPTS/run/domaintotal_new.x
set path_parser   = $SCRIPTS/run/path_parser.py
set aggregator    = /garnet/work/bte/WO109.1/pyagg/agg.py
set windower      = /garnet/work/bte/WO109.1/pyagg/subset.py
#set sectorlist_parser  = $SCRIPTS/run/sectorlist_parser.py

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
## Define MONTHS_LIST and SPINUP_LIST
$EMF_CLIENT -k $EMF_JOBKEY -m "Running set_months" -x $set_months  ## log w/ EMF server
switch ( $#argv )
   case 3:
      source $set_months $argv[2] "$argv[3]"
      if ( $status != 0 ) set exitstat = 1
   breaksw
   case 4: 
      source $set_months $argv[2] "$argv[3]" $argv[4]
      if ( $status != 0 ) set exitstat = 1
   breaksw
   case 5:
      source $set_months $argv[2] "$argv[3]" $argv[4]
      if ( $status != 0 ) set exitstat = 1
      setenv TLABEL $argv[5]
   breaksw
endsw
if ( $exitstat != 0 ) then
    echo "ERROR: setting months"
    $EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: setting months" -t "e" -x $set_months  
    exit (1)
endif

# Set spinup duration - the set_months will have QA'f the $argv[4] value
if ( $#argv > 3 ) setenv SPINUP_DURATION $argv[4]

## Set naming label
set namelabel = ${SECTOR}_${CASE}_${GRID}
if ( $?TLABEL != 0 ) then
    set namelabel = ${namelabel}_$TLABEL
endif
if ( $?JOB_GROUP != 0 ) then
    set namelabel = ${namelabel}_$JOB_GROUP
endif

## Set spinup list to local variable, so can use as array
set spinup_list = ( $SPINUP_LIST )

## Record the helper scripts being used
set suffix = _$namelabel.txt
echo "# Helper scripts used for mrggrid" > $LOGS/helper_scripts_list$suffix
echo $emf_cleanup >> $LOGS/helper_scripts_list$suffix
echo $set_months >> $LOGS/helper_scripts_list$suffix
echo $timetracker >> $LOGS/helper_scripts_list$suffix
echo $smk_run >> $LOGS/helper_scripts_list$suffix
echo $m3stat >> $LOGS/helper_scripts_list$suffix
echo $set_days >> $LOGS/helper_scripts_list$suffix
echo $log_analyzer >> $LOGS/helper_scripts_list$suffix
echo $msg_list >> $LOGS/helper_scripts_list$suffix
echo $path_parser >> $LOGS/helper_scripts_list$suffix
echo $aggregator >> $LOGS/helper_scripts_list$suffix
echo $windower >> $LOGS/helper_scripts_list$suffix
#echo $sectorlist_parser >> $LOGS/helper_scripts_list$suffix

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

# Initialize sector arrays, assuming maximum number of sectors of 60 - upgraded from 40 by C. Allen (CSC) on 17feb2011.  Added project_root column for new platform. by J. Beidler(CSC) 19Apr2011
set sector_name        = (0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0)
set sector_case        = (0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0)
set sector_mtype       = (0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0)
set sector_spin        = (0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0)
set sector_col         = (0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0)
set sector_year        = (0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0)
set sector_datefile    = (0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0)
set sector_file        = (0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0)
set sector_ev          = (N N N N N N N N N N N N N N N N N N N N N N N N N N N N N N N N N N N N N N N N N N N N N N N N N N N N N N N N N N N N)
set sector_fileset_cnt = (0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0)
set sector_orig        = (0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0)
set sector_spc         = (0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0)
set sector_mrgg        = (0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0)
set sector_project     = (0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0)

# Initialize leap years that we'll handle
set leap_years = (1992 1996 2000 2004 2008 2012 2016 2020 2024 2028 2032 2036)
set daysinmonth = ( 31 28 31 30 31 30 31 31 30 31 30 31)

## Adjust daysinmonth for leap years
echo $leap_years | /bin/grep -q $BASE_YEAR 
if ( $status == 0 ) then
   set daysinmonth[2] = 29
endif

# Other initializations
setenv RUN_MRGGRID Y    # Run the mrggrid program
set exitstat = 0

## Define the mrggrid reports if using adjust factor or tagging
setenv REPMERGE_ADJ $REPOUT/mrggrid/$GRID/$SPC/rep_adj_${SECTOR}_{$namelabel}.csv
setenv REPMERGE_SUM $REPOUT/mrggrid/$GRID/$SPC/rep_adjsummary_${SECTOR}_{$namelabel}.csv
setenv REPMERGE_TAG $REPOUT/mrggrid/$GRID/$SPC/rep_tag_${SECTOR}_{$namelabel}.csv

## remove mrggrid reports if they exist so doesn't keep appending
if ( -e $REPMERGE_ADJ ) then
    rm $REPMERGE_ADJ
    rm $REPMERGE_SUM
endif
if ( -e $REPMERGE_TAG ) then
    rm $REPMERGE_TAG
endif

## Check if top of output directory exists, if not, raise an error
if ( ! -e $OUT_ROOT ) then
    echo "ERROR: Output root doesn't exist, $OUT_ROOT"
    $EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: Output root doesn't exist, $OUT_ROOT" -t "e"
    set exitstat = 1
    goto end_of_script
endif

if ( ! $?RUN_SOURCESENS ) then
    setenv RUN_SOURCESENS N
endif

### If source sensitivity, need source sector override file
#if ( $RUN_SOURCESENS == Y ) then
#    echo "Running merge for source sensitivity"
#    set sector_override = $IMD_ROOT/mrggrid/source_sector_override_${CASE}_${GRID}.txt
#    if ( ! -e $sector_override ) then
#	echo "SCRIPT ERROR: source sector override file doesn't exist "
#	$EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: " -t "e"
#	set exitstat = 1
#	goto end_of_script
#    endif
#
#    ## Unique list of sectors to override (ie. set to source sens case over parent case)
#    set source_sense_sectors = `sort -u $sector_override`
#    echo "SCRIPT NOTE: Sectors to merge from source sensitivity case: $source_sense_sectors"
#    $EMF_CLIENT -k $EMF_JOBKEY -m "Sectors to merge from source sensitivity case: $source_sense_sectors"
#endif

## Get root of the merge dates directory from smk merge date files, if exists
## smk merge date files should be in $mrgdate_root (including december of previous year)
if ( $?MRGDATE_FILES) then
    echo "SCRIPT NOTE: Getting representative merge dates path from MRGDATE_FILES"
    set mrgdate_root = `$path_parser $MRGDATE_FILES`
endif    

    
# Loop through months specified in arguments, incl spinup period
   # Use same approach as is used in main run scripts
set mc = 0
foreach mon ( $MONTHS_LIST )    # Note that some months can appear more than once with different spin-up group

   @ mc = $mc + 1 

   # Initialise the days in the month of interest
   #    and set the year based on the spinup group for current month
   set daysmonth = $daysinmonth[$mon]
   set date_counter = 0
   set year = $BASE_YEAR

   # Adjust the start day and the days in month based on the "spinup" status of current month
   switch ( $spinup_list[$mc] )

      case a:                         # All days is same as initialized values
         echo "SCRIPT NOTE: Using all dates for month $mon during merge"
         set year = $BASE_YEAR
      breaksw
 
      case p:                         # All days but spinup days at end of month
         echo "SCRIPT NOTE: Using non-spinup dates for month $mon during merge"
         @ daysmonth = $daysmonth - $SPINUP_DURATION
         set year = $BASE_YEAR
      breaksw

      case s:                         # Spinup days at end of month
         echo "SCRIPT NOTE: Using spinup dates for month $mon during merge"
         @ date_counter = $daysmonth - $SPINUP_DURATION
         set year = $BASE_YEAR
      breaksw

      case sp:                        # Spinup days at end of month, previous year
         echo "SCRIPT NOTE: Using spinup dates for month $mon of previous year during merge"
         @ date_counter = $daysmonth - $SPINUP_DURATION
         @ year = $year - 1
      breaksw

   endsw

   # Set 2-digit month
   set mon2 = $mon
   if ( $mon2 < 10 ) set mon2 = 0$mon2

   # Ensure director for FILELIST files is created and has proper permissions 
   if ( ! -e $IMD_ROOT/merge/$SECTOR ) then
      mkdir -p $IMD_ROOT/merge/$SECTOR
      chmod ug+rwx $IMD_ROOT/merge/$SECTOR
   endif

   if ( ! -e $PREMERGED/../mrggrid ) then
      mkdir -p $PREMERGED/../mrggrid
      chmod ug+rwx $PREMERGED/../mrggrid
   endif

   # Define FILELIST name for current month
   setenv FILELIST $IMD_ROOT/merge/$SECTOR/filelist_${CASE}_${GRID}_$mon.lst

   # Initialize MRG_DIFF_DAYS to N as default (is reset based on merge type from sectors file)
   setenv MRG_DIFF_DAYS N

   # (J. Beidler, 3 Sep 2014) Remove all refueling separation. No longer a separate sector.
   # Add RPH as an onroad sector for MOVES 2014.
   set nsector = 4  # RPD/RPH/RPV/RPP   

   if (! $?NO_RPH) set NO_RPH = ""
   if (! $?NO_RPP) set NO_RPP = ""
   if (! $?NO_RPV) set NO_RPV = ""
   if (! $?INCLUDE_RPHO) set INCLUDE_RPHO = ""
   if (! $?INCLUDE_RPS) set INCLUDE_RPS = ""

   # (C. Allen, 12 Feb 2016) Add support for merges without either RPH or RPP.
   # (13 Dec 2019) Add support for RPHO
   # (21 Oct 2020) Add support for RPS, and also add a NO_RPV flag
   if ($NO_RPH == Y) then
      @ nsector = $nsector - 1
   endif
   if ($NO_RPP == Y) then
      @ nsector = $nsector - 1
   endif
   if ($NO_RPV == Y) then
      @ nsector = $nsector - 1
   endif
   if ($INCLUDE_RPHO == Y) then
      @ nsector = $nsector + 1
   endif
   if ($INCLUDE_RPS == Y) then
      @ nsector = $nsector + 1
   endif

   echo "nsector = $nsector"

   # Check for maximum of 60 sectors and abort if count is higher, upgraded from 40 by C. Allen (CSC) on 17feb2011
   if ( $nsector > 60 ) then
      echo "SCRIPT ERROR: smk_sectormerge_emf.csh is limited to a maximum"
      echo "              of sixty sectors, but $nsector sectors attempted"
      $EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: limit of 60 sectors, but $nsector sectors requested" -t "e"
      exit( 1 )
   endif

   # Preprocess and QA sector list and settings
   # Loop through sectors
   set nl = 0
   set seccnt = 0
   set outcnt = 0
   
   ## test the number of environment variables
   if ( $debug_script == Y ) then
      ## test the environment
      env
      set nenv_vars = `env |wc -l`
      echo "SCRIPT NOTE: Number of environment variables in env, $nenv_vars"

   endif
   
   set rpxs = ("RPD")
   if ($NO_RPH != Y) set rpxs = ($rpxs " RPH")
   if ($NO_RPP != Y) set rpxs = ($rpxs " RPP")
   if ($NO_RPV != Y) set rpxs = ($rpxs " RPV")
   if ($INCLUDE_RPHO == Y) set rpxs = ($rpxs " RPHO")
   if ($INCLUDE_RPS == Y) set rpxs = ($rpxs " RPS")
   
   foreach moves_type ( $rpxs )

         @ nl = $nl + 1
         @ seccnt = $seccnt + 1
	 
         set tmpsector              = $moves_type
	 set sector_name[$seccnt]   = $tmpsector
         set tmpcase                = $CASE
         set sector_year[$seccnt]   = $BASE_YEAR
         set sector_mtype[$seccnt]  = "all"
         set sector_spin[$seccnt]   = "actualMet"
	 
	 # Check for optional parameters, 

#	 ## speciation; if it's not there, then set to default $SPC
	 set sector_spc[$seccnt]    = $SPC
	 
#	 ## merge sector, if not there set to default Y
         set sector_mrgg[$seccnt]   = Y	 

#	 ## project root; if it's not there, then set to default $PROJECT_ROOT
	 set sector_project[$seccnt]    = $PROJECT_ROOT

         set sector_case[$seccnt] = $tmpcase

         # If previous-year spinup period for output, adjust sector's year depending 
         #      on sector's approach for the previous-year spinup
         set tmpyear = $sector_year[$seccnt]
         set skipsector = 0
         if ( $spinup_list[$mc] == sp ) then
            if ( $sector_spin[$seccnt] == actualMet || ( $sector_spin[$seccnt] == SectBaseYr && $sector_mtype[$seccnt] != all ) ) then
                @ tmpyear = $tmpyear - 1
            else
               if ( $sector_spin[$seccnt] == none ) then
                  set skipsector = 1
	       endif  
            endif
         endif

         if ( $skipsector == 1 ) goto end_of_preproc

	 # Define environment variable for sector
         # Get length of sector name.
         set charcnt = `echo $sector_name[$seccnt] | wc -c` 

#        # set environment variable as uppercase sector name, appended with counter if necessary to avoid duplicates
         set sector_ev[$seccnt] = `echo $sector_name[$seccnt] | tr '[:lower:]' '[:upper:]'`

#         # Set MRG_DIFF_DAYS for current day if any sectors not using "all"
         if ( $sector_mtype[$seccnt] != all ) setenv MRG_DIFF_DAYS Y

         ## Set datesfile name 
	 if ( $?MRGDATE_FILES) then
            ## have representative dates file, therefore generate from it's path
            set sector_datefile[$seccnt] = ${mrgdate_root}/smk_merge_dates_$tmpyear$mon2.txt
	 else 
	    ## generate from  current sector/case/year (note SECTOR should = mrggrid)
            set sector_datefile[$seccnt] = ${DAT_ROOT}/inputs/$sector_case[$seccnt]/$SECTOR/smk_merge_dates_$tmpyear$mon2.txt
	 endif

	 ## ensure it exists
         if ( ! -e $sector_datefile[$seccnt] ) then
            echo "SCRIPT ERROR: Date file name based on SectorCase and SectBaseYr is not found:"
            echo "              "$sector_datefile[$seccnt]
            $EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: Date file $sector_datefile[$seccnt] not found" -t "e"
            set exitstat = 1
	    goto end_of_script
         endif

         # Determine acceptable merge types from the month-specific dates file
         set mrgtypes = ( `head -1 $sector_datefile[$seccnt] | cut -d, -f2- | sed 's/,//g'` )

         # Check that merge type is one of the types available in the month-specific date file
         # Set column number from $sector_datefile[$seccnt] (col initialized to 1 on purpose to skip merge date column)
         set found = N    # initialize
         set col = 1      # initialize
         foreach t ( $mrgtypes )
            @ col = $col + 1
            if ( $t == $sector_mtype[$seccnt] ) then
               set found = Y
               set sector_col[$seccnt] = $col
            endif
         end

         if ( $found == N ) then
            echo SCRIPT ERROR: Merge type \"$sector_mtype[$seccnt]\" is not found in header of the merge dates file:
            echo "              "$sector_datefile[$seccnt]
            $EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: merge type $sector_mtype[$seccnt] not found in date file" -t "e"
            exit ( 1 )
         endif

         # If got this far (i.e., not a spinup period with a sector using SpinupApproach = none)
         #     then iterate output counter
         @ outcnt = $outcnt + 1

         # Write this sector's variable to FILELIST for input to Mrggrid
         if ( $outcnt == 1 ) then
            echo $sector_ev[$seccnt] > $FILELIST
         else
            echo $sector_ev[$seccnt] >> $FILELIST
         endif

         # To skip parts of loop when in previous-year spinup and PrevYrSpinup=none
         end_of_preproc:

   end      # End loop through sectors

   if ( $debug_script == Y ) then
   ## debug
       echo "sectors: $sector_name"
       echo "sector cases: $sector_case"
       echo "sector mtypes: $sector_mtype"
       echo "sectors speciation: $sector_spc"
       echo "sectors merge: $sector_mrgg"
   endif

   # Loop through days in month
   while ( $date_counter < $daysmonth )

      @ date_counter = $date_counter + 1     

      # Add leading zeros to day and month 
      set day   = $date_counter
      if ( $day < 10 ) set day = 0$day

      # Set current output date (YYYYMMDD)
      set date = $year$mon2$day

      ## Set EMF_PERIOD for this day
      setenv EMF_PERIOD $date

      setenv RUN_PART4 Y

      ## Set dates now, so that we can reference RUNSET and skip this date if it's not supposed to be run
      ## Since not using Assigns file to handle day changes, must set dates
      #   manually for log file naming and controlling the run
      setenv ESDATE   $date

      set juldate = `$IOAPIDIR/juldate $mon $day $year | /bin/grep $year | cut -d, -f2`
      setenv G_STDATE $juldate

      ## Reset RUN_MRGGRID and RUN_M3STAT before calling runsettings.csh
      setenv RUN_MRGGRID Y
      setenv RUN_M3STAT $run_m3stat
#      setenv RUN_DOMAIN $run_domain

      if ( $?RUNSET ) then
         if ( -e $RUNSET ) then
            source $ASSIGNS/runsettings.csh $SECTOR $GRID PART4 $G_STDATE $RUNSET
	 endif
      endif
      
      ## If runsettings.csh set RUN_MRGGRID to N, then skip to the next date in the loop immediately,
      ## so that we don't do any unnecessary file unzips
      if ( $RUN_MRGGRID == N ) then
         echo "## Skipping merge for $mon2/$day/$year, based on RUNSET"
         continue
      endif
   
      ## Print message to stdout log
      echo "## Running merge script for $mon2/$day/$year, using data files:"

      ## Loop through sectors
      set seccnt = 0
      while ( $seccnt < $nsector )

         @ seccnt = $seccnt + 1

	 ## check if merging in this sector, if not skip to next sector
	 if ( $sector_mrgg[$seccnt] == N ) goto end_of_sectorloop
    
         if ( $sector_case[$seccnt] == none ) then

            ## For SectorCase = "none", interpret environment variables in string provided
            set namevars_tmp = ( `echo $sector_mtype[$seccnt] | sed 's/\$/ /g'` )

            ## Loop through EVs provided and set the file name (previously they have been checked that they are set)
            set sector_file[$seccnt] = ""
            foreach tmpvar ( $namevars_tmp )
               set nextpart = `env | /bin/grep ^"$tmpvar=" | cut -d\= -f2`
               set sector_file[$seccnt] = $sector_file[$seccnt]$nextpart
            end

            ## Define sector's environment variable (in filelist) to the actual file name
            setenv $sector_ev[$seccnt] $sector_file[$seccnt]
    
         ## Otherwise, for most SectorCases (not equal to "none")
         else
            # If previous-year spinup period for output, adjust sector's year depending 
            #      on sector's approach for the previous-year spinup
            set tmpyear = $sector_year[$seccnt]
            set skipsector = 0
            if ( $spinup_list[$mc] == sp ) then
               if ( $sector_spin[$seccnt] == actualMet || ( $sector_spin[$seccnt] == SectBaseYr && $sector_mtype[$seccnt] != all ) ) then
                   @ tmpyear = $tmpyear - 1
               else
                  if ( $sector_spin[$seccnt] == none ) then
                     set skipsector = 1
                  endif
               endif
            endif

            if ( $skipsector == 1 ) goto end_of_sectorloop

            # Set date to look for in the datefile, by appending sector year with 
            #   month and day.
            set searchdate = $tmpyear$mon2$day

            # Handle case when current year is a leap year, but referenced year is
            #   not, by assigning different date that matches correct day of week.
            # Assume that subtracting 7 days will not land on a holiday in Feb.
            if ( $mon2 == '02' && $day == 29 ) then
               # Search for sector year in list of leap years
               echo $leap_years | /bin/grep -q $sector_year[$seccnt]

               # If it's not a leap year, then change the searchdate
               if ( $status == 1 ) then
                  set searchdate = $tmpyear${mon2}'22'
                  echo "SCRIPT WARNING: Using 02/22/$tmpyear to instead of 02/29/$tmpyear"
                  echo "                since $sector_name[$seccnt] sector year $tmpyear is not a,"
                  echo "                leap year, but run year $year is a leap year."
               endif
            endif

            #  Retrieve date to use for current sector/year/temporal-type (set by sector_col)
            set sector = $sector_name[$seccnt]                                     # Set tmp sector name
            set case   = $sector_case[$seccnt]                                     # Set tmp case for this sector
	    set speciation = $sector_spc[$seccnt]                                  # Set tmp speciation for this sector
            set case_imd = $sector_project[$seccnt]/$case/premerged              # Set tmp project/case intermediate path for this sector
            set esdate = `/bin/grep ^$searchdate $sector_datefile[$seccnt] | cut -d\, -f$sector_col[$seccnt]` # extract column from sector_datefile for current sector

            set filebase = $case_imd/$SECTOR/$sector/emis_${MM}_${sector}_${SECTOR}_${esdate}_${GRID}_${speciation}_${case}
	    set filebasetmp = $case_imd/mrggrid/emis_${MM}_${sector}_${SECTOR}_${esdate}_${GRID}_${speciation}_${case}

            # Determine if current day's file is a fileset or not based on the name
            if ( -e $filebase.1.ncf ||  \
                 -e $filebase.1.ncf.gz ) then
               set sector_fileset_cnt[$seccnt] = `ls -1 $filebase.*.ncf* | wc -l`
            endif

            # Define environment variable for current sector/day to actual file name
            # Note: use same name even for filesets, since FilesetAPI deals with that
            ###################### Set input file name #########################
            set    sector_orig[$seccnt] = $filebase.ncf
	    set    sector_file[$seccnt] = $filebasetmp.${EMF_JOBID}.tmp.ncf
            setenv $sector_ev[$seccnt] $sector_file[$seccnt]

            # This next code segment is inside the "sectorcase != none" wrapper, so that it is not
	    # performed when sectorcase == none.

            # Unzipping and file checks for non-fileset files
            if ( $sector_fileset_cnt[$seccnt] == 0 ) then

               # If non-fileset file is zipped, unzip and set flag for later zipping for this sector
               if ( -e $sector_orig[$seccnt].gz ) then
                  echo "      " $sector_name[$seccnt] unzipping
                  /bin/zcat -f $sector_orig[$seccnt].gz > $sector_file[$seccnt]

               # If it is unzipped, then set a link to it in the temporary directory
   	       else  
   	          ln -fs $sector_orig[$seccnt] $sector_file[$seccnt]   
               endif
         
               # If file does not exist, report error and continue
               if ( ! -e $sector_file[$seccnt] ) then
                  echo "SCRIPT ERROR: data file expected but not found:"
                  echo "              "$sector_file[$seccnt]
                  $EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: Data file expected but not found. See run log." -t "e"
                  set exitstat = 1

               ## Otherwise, print message to stdout log
               else
                  echo "    ##" $sector_name[$seccnt]  $sector_file[$seccnt]
               endif

            # Unzipping and file checks for fileset files
            else
               # Loop through fileset
               set tmpcnt = 0
               while ( $tmpcnt < $sector_fileset_cnt[$seccnt] )
                  @ tmpcnt = $tmpcnt + 1
                  set tmpfile = $filebase.$tmpcnt.ncf

                  # If file in fileset is zipped, unzip it.
                  if ( -e $tmpfile.gz ) then
                     echo "      " $sector_name[$seccnt] unzipping fileset file number $tmpcnt
	             /bin/zcat -f $tmpfile.gz > $filebasetmp.${EMF_JOBID}.tmp.$tmpcnt.ncf

                  # If it is unzipped, then set a link to it in the temporary directory
                   else  
                     ln -fs $tmpfile $filebasetmp.${EMF_JOBID}.tmp.$tmpcnt.ncf
                  endif

                  # Double check to make sure zip worked and file is present (and not deleted, moved, or migrated)
                  if ( ! -e $filebasetmp.${EMF_JOBID}.tmp.$tmpcnt.ncf ) then
                     echo "SCRIPT ERROR: data file expected but not found:"
                     echo "              "$filebasetmp.${EMF_JOBID}.tmp.$tmpcnt.ncf
                     $EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: Data file expected but not found. See run log." -t "e"
                     set exitstat = 1
                  else
                     echo "    ##" $sector_name[$seccnt] $filebasetmp.${EMF_JOBID}.tmp.$tmpcnt.ncf 
                  endif

               end  # End loop over files in fileset

            endif  # if ( $sector_fileset_cnt[$seccnt] == 0 )

         endif   # if SectorCase is "none" or not

         #  Skip to here when previous-year spinup and Sector has approach "none"
         end_of_sectorloop:

      end      # End loop on sectors

      ## Abort now if exitstatus is non-zero (because one or more files don't exist)
      if ( $exitstat > 0 ) then
         echo "SCRIPT ERROR: One or more input files expected by script"
         echo "              does not exist."   
         exit ( $exitstat ) 

      endif

      ## Set output file names
      setenv OUTFILE    $case_imd/$SECTOR/emis_${MM}_${SECTOR}_${date}_${GRID}_${SPC}_${CASE}.ncf
      setenv REPMRGGRID $REPOUT/mrggrid/$GRID/$SPC/mrggrid_dates_${date}_${SECTOR}.rpt
      if ( ! -e $REPOUT/mrggrid/$GRID/$SPC ) then
         mkdir -p $REPOUT/mrggrid/$GRID/$SPC
         chmod ug+rwx $REPOUT/mrggrid/$GRID/$SPC
      endif
      
      ###################### Run Mrggrid program #########################
      # Run mrggrid program
      $smk_run
      if ( $status != 0 ) then
          echo "ERROR: Running smk_run for one-time steps"
          $EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: Running smk_run for one-time steps" -t "e" -x $smk_run -p $EMF_PERIOD
          set exitstat = 1
          goto end_of_script
      endif

      ## Loop through sectors for zipping and removing temporary files
      set seccnt = 0
      while ( $seccnt < $nsector )

         @ seccnt = $seccnt + 1

         # Remove the temporary files and links from the temporary directory

         # For non-fileset files
	 # Don't delete file if sectorcase == none, because then sector_file is the original file, not a softlink
         if ( $sector_fileset_cnt[$seccnt] == 0 && $sector_case[$seccnt] != none ) then
            /bin/rm -f $sector_file[$seccnt]

         # For fileset files
         else
            set sector = $sector_name[$seccnt]
            set case   = $sector_case[$seccnt]
	    set speciation = $sector_spc[$seccnt]
            set case_imd = $sector_project[$seccnt]/$case/premerged              # Set project/case intermediate path for this sector
            set esdate =  `/bin/grep ^$searchdate $sector_datefile[$seccnt] | cut -d\, -f$sector_col[$seccnt]`
            set tmpcnt = 0
            while ( $tmpcnt < $sector_fileset_cnt[$seccnt] )
               @ tmpcnt = $tmpcnt + 1
               set tmpfile = $filebasetmp.${EMF_JOBID}.tmp.$tmpcnt.ncf
               /bin/rm -f $tmpfile

            end  # End loop over fileset

         endif

      end        # End loop on sectors for zipping  

      ## Run M3STAT on output file from Mrggrid, depending on runset settings
      if ( -e $OUTFILE && $RUN_M3STAT == Y ) then
         $m3stat $OUTFILE
         if ( $status != 0 ) then
	    echo "ERROR: running m3stat_chk for $date"
	    $EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: running m3stat_chk for $date" -t "e" -x $m3stat -p $EMF_PERIOD
            set exitstat = 1
	    goto end_of_script
         endif
      endif

      ## If override exists, set RUN_DOMAIN to override value.
      ## An alternative method of turning on or off RUN_DOMAIN
      if ( $?OVERRIDE_RUN_DOMAIN != 0 ) then
	    echo "Overriding RUN_DOMAIN, setting to $OVERRIDE_RUN_DOMAIN"
	    setenv RUN_DOMAIN $OVERRIDE_RUN_DOMAIN
      endif

      # Run domain totals, but not for spinup
      # Commented out for SMOKE-MOVES
#      if ( -e $OUTFILE && $RUN_DOMAIN == Y ) then

#         ## Define "permanent" setting for OUTPUT, since that EV is used by the domain totals program
#         set output_sav = $OUTPUT

#         setenv FILENAME $OUTFILE
#         setenv OUTPUT $REPOUT/$SECTOR/$GRID/$SPC/domain_totals_$date.rpt
#        unsetenv LOGFILE
#         unsetenv LOG_FILE
#	 echo "Running domaintotals for $date"
#         $domain_totals
#         if ( $status != 0 ) then
#            echo "ERROR: running domain totals for $date"
#            $EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: running domain totals for $date" -t "e" -x $domain_totals -p $EMF_PERIOD
#            set exitstat = 1
#	    goto end_of_script
#         endif

#	 setenv OUTPUT $output_sav

 #     endif

      # For SMOKE-MOVES, run the windower to create 36km emissions, and the aggregator to create E12/W12 emissions
      # Currently, these pythons scripts - $aggregator and $windower as defined above - are hard-wired to ALWAYS
      # create emissions for all three grids based on the 12MERGEUS1 emissions. This can be changed
      # down the road if necessary. At some point I may add flags which can be set from the EMF that turn off
      # either the aggregator or the windower, or both.
      
      # If the grid isn't 12MERGEUS1, skip this step
      if ($GRID == 12MERGEUS1) then

         echo "Running aggregator to create 36US1 emissions from 12MERGEUS1 emissions..."      
         /usr/local/bin/python2.6 $aggregator
         if ( $status != 0 ) then
	    echo "ERROR: running 36US1 aggregator for $date"
	    $EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: running 36US1 aggregator for $date" -t "e" -x $aggregator -p $EMF_PERIOD
            set exitstat = 1
	    goto end_of_script
         endif

         echo "Running windower to create 12EUS1 and 12WUS1 emissions from 12MERGEUS1 emissions..."      
         /usr/local/bin/python2.6 $windower
         if ( $status != 0 ) then
	    echo "ERROR: running 12km windower for $date"
	    $EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: running 12km windower for $date" -t "e" -x $windower -p $EMF_PERIOD
            set exitstat = 1
	    goto end_of_script
         endif
      
      endif

      ## Zip Mrggrid outputs if setting requests this
      if ( $GZIP_OUTPUTS == Y ) then
         /bin/gzip -f $OUTFILE
      endif
 
   end     # End loop through days in month

   ## send EMF message back that completed this month
   echo "Completed sector merge for $mon2/$year"
   $EMF_CLIENT -k $EMF_JOBKEY -m  "Completed sector merge for $mon2/$year" 

end        # End loop through months

## Set file extension for counting and registering Mrggrid outputs,
#    which depends on zipping option set before calling script.
set fileext = ncf
if ( $GZIP_OUTPUTS == Y ) set fileext = ncf.gz

## Determine if all smkmerge model-ready files are available to import
set mrgcnt = `/bin/ls -1 $case_imd/$SECTOR/emis_${MM}_$SECTOR_*_${GRID}*.$fileext | wc -l`

## Set day count depending on leap year and spinup period
set alldays = 365
echo $leap_years | /bin/grep -q $BASE_YEAR
if ( $status == 0 ) set alldays = 366

@ alldays = $alldays + $SPINUP_DURATION

## Register Smkmerge AQM-ready outputs with EMF for entire year
## or if set to always register outputs
#if ( $REGISTER_AQMOUT == Y ) then
#   if ( $mrgcnt == $alldays || $REGISTER_AQMOUT_ALWAYS == Y) then
#      echo "SCRIPT NOTE: Registering Mrggrid outputs, $mrgcnt model-ready onroad files"
#      $EMF_CLIENT -k $EMF_JOBKEY \
#           -D $OUTPUT \
#           -P "emis_${MM}_all_*_${GRID}_*$fileext" \
#           -T "CMAQ Model Ready Emissions: Merged" \
#           -N "Model-ready CMAQ $CASE $GRID" \
#           -O "AQM-ready data $GRID"
#   else
#      echo "SCRIPT NOTE: Not registering Mrggrid model-ready files because only found"
#      echo "             $mrgcnt Mrggrid outputs out of expected $alldays" 
#
#   endif
#endif

# Label for the end of the script, used during script abort
end_of_script:

## Register time log
#if ( -e $TIMELOG ) then
#   echo "SCRIPT NOTE: Registering time log"
#   $EMF_CLIENT -k $EMF_JOBKEY -F $TIMELOG -T "SMOKE time log (External)" -N "SMOKE timelog $namelabel" -O "Timelog $namelabel (External)"
#endif

### register mrggrid adjusted report
#if ( -e $REPMERGE_ADJ ) then
#    echo "SCRIPT NOTE: Registering MRGGRID adjustment reports"
#    $EMF_CLIENT -k $EMF_JOBKEY -F $REPMERGE_ADJ -T "Smkreport mrggrid adjustment" -O "mrggrid Smkreport adjustment sectors $namelabel" 
#    $EMF_CLIENT -k $EMF_JOBKEY -F $REPMERGE_SUM -T "Smkreport mrggrid adjustment" -O "mrggrid Smkreport adjustment summary $namelabel" 
#endif

### register mrggrid tagging report
#if ( -e $REPMERGE_TAG ) then
#    echo "SCRIPT NOTE: Registering MRGGRID tagging reports"
#    $EMF_CLIENT -k $EMF_JOBKEY -F $REPMERGE_TAG -T "Smkreport mrggrid adjustment" -O "mrggrid Smkreport tagging sectors $namelabel" 
#endif

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

## If source sensitivity, register sector override file
if ( $RUN_SOURCESENS == Y ) then
#   $EMF_CLIENT -k $EMF_JOBKEY \
#	-F  $sector_override -T "Sector List Sensitivity Override (External)" \
#	-O "sector override ${namelabel}"    
endif

## Ending of script
#
exit( $exitstat )

