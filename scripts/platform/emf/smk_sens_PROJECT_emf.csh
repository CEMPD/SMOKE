#!/bin/tcsh -f

# Version @(#)$Id$
# Path    $Source$
# Date    $Date$

# This script sets up needed environment variables for applying
# PROJECTION packets to annual SMOKE (point or area) sources.  This
# process requires running SMOKE through the following modules:
# cntlmat, grwinven (optional), and smkmerge.  

# This script is intended to be used with the EMF source emissions in
# SMOKE for the EPA 2002 modeling platform, and calls the scripts that
# runs the SMOKE programs.

#
# Script created by : M. Houyoux, Environmental Protection Agency July, 2007 
# Modified to work w/ EMF: A. Zubrow, UNC - IE, August, 2007 
# Modified to apply PROJECTION and/or CONTROL packets by R. Mason, EPA, March, 2008 
# Modified for PROJECTION packets using PARENT CASE: A. Zubrow, UNC Dec, 2008
# Modified to fix variable name for reruns: J. Beidler, CSC, Jul, 2009 
##*********************************************************************

## log w/ EMF server that script is running
# setenv EMF_CLIENT false (in script that calls this).... to test interactively
# also, build EMF case, which will create a script that calls this (with parameters filled in)
# THEN, setenv EMF_CLIENT false and try executing this!!
$EMF_CLIENT -k $EMF_JOBKEY -s "Running" 

## month-dependent programs
setenv RUN_SMKINVEN  N        #  run inventory import program
setenv RUN_SPCMAT    N        #  run speciation matrix program
setenv RUN_CNTLMAT   Y        #  run control matrix program
setenv RUN_GRWINVEN  N        #  run control application program

## time-dependent programs
set run_smkmerge = Y       #  run merge program
set run_smk2emis = N       #  run conversion of 2-d to UAM binary
setenv RUN_TEMPORAL N      #  run temoral allocation program

## quality assurance
set run_smkreport = Y      # Y runs reporting for state reports
set run_m3stat    = Y      # Y runs the m3stat program on the Smkmerge outputs

if (! $?REGISTER_REPOUT) then
   setenv REGISTER_REPOUT    Y       # Imports Smkreport and Smkmerge reports into EMF
endif
setenv REGISTER_GC_OUT    Y       # Imports Grwinven ORL inventory into EMF
setenv PROMPTFLAG         N       # Y (never set to Y for batch processing)
setenv AUTO_DELETE        Y       # Y deletes SMOKE I/O API output files (recommended)
setenv AUTO_DELETE_LOG    Y       # Y automatically deletes logs without asking
if ( ! $?DEBUGMODE ) then
   setenv DEBUGMODE          N       # Y changes script to use debugger
endif
setenv DEBUG_EXE    totalview     # Sets the debugger to use when DEBUGMODE = Y

##############################################################################

## Get arguments
switch ( $#argv )
   case 0:
   case 1:
   case 2:
   case 3:
   case 4:
   case 5:
      echo "SCRIPT ERROR: Script requires arguments for a emission period, "
      echo "              source category, grid name, I/O API grid name and the "
      echo "              -m or -q option with 2-3 settings."
      echo " "
      echo "  This script expects to be called using one of the following argument lists:"
      echo "     <emiss_period><source categ> <grid abbrv> <I/O API gridname> -m <monthlist> <spinup> <label>"
      echo "     <emiss_period><source categ> <grid abbrv> <I/O API gridname> -q <quarters> <spinup> <label>"
      echo " "
      echo "  You can either use one approach or the other (differing by the -m or -q options)."
      echo " "
      echo "  In the above list, the arguments are defined as follows:"
      echo "     <emiss_period>     : Emission period, either monthly or annual"
      echo "     <source categ>     : Source category, e.g. A, M, or P"
      echo "     <grid abbrv>       : Grid abbreviation (e.g., 36US1)"
      echo "     <I/O API gridname> : I/O API gridname that needs to match entry in the"
      echo "                          GRIDDESC input file"
      echo "     <monthlist>        : list of months to run when using the -m option"
      echo "     <quarters>         : list of quarters to run when using the -q option"
      echo "     <spinup>           : set to number of days between 1 and 20 to run a spinup" 
      echo "                          period (value sets number of days), and N otherwise"
      echo "     <label>            : label to put on TIMELOG file and helper-scripts list (optional)" 
      echo " "
      echo "  Examples:"
      echo "     <script name> annual A 36US1 36US1_148X112 -m '1 2 3' 0 jan-sep"
      echo "              This example runs the annual script for area source category, "
      echo "              for Jan, Feb, & Mar on the 36US1 grid, with no spinup days and"
      echo "              gives a label to the TIMELOG file of jan-sep." 
      echo " "
      $EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: smoke script did not receive more than 4 arguments" -t "e"
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
      exit( 1 )
endsw

## Determine if running based on annual or monthly emissions
if ( "$argv[1]" == monthly ) then
    set run_monthly = Y
else if ( "$argv[1]" == annual ) then
    set run_monthly = N
else
    echo "SCRIPT ERROR: First argument must equal 'annual' or 'monthly'"
      $EMF_CLIENT -k $EMF_JOBKEY -m "SCRIPT ERROR: First argument must equal 'annual' or 'monthly'" -t "e"
    exit( 1 )

endif 

# Get the next three options for the source category, grid abbreviation and I/O API grid
setenv SMK_SOURCE "$argv[2]"
setenv GRID "$argv[3]"
setenv IOAPI_GRIDNAME_1 "$argv[4]"

# category to process during merge
setenv MRG_SOURCE $SMK_SOURCE # source category to merge


## SMKMERGE use control matrix
setenv MRG_CTLMAT_MULT $MRG_SOURCE  # 'A' merges with area projection matrix produced from CNTLMAT
setenv MRG_REPCTL_YN   Y            # Y Report control totals 

## If point category, setup some additional parameters
if ( $SMK_SOURCE == P ) then

    setenv EISECTOR $SECTOR

    ## default INLINE_MODE to off
    if ( $INLINE_MODE == "" ) then
	setenv INLINE_MODE "off"
    endif
    if ( $INLINE_MODE == "off" ) then
	setenv MRG_LAYERS_YN Y
	setenv SMK_ELEV_METHOD 1
    else if ( $INLINE_MODE == "both" ) then
	setenv MRG_LAYERS_YN N
	setenv SMK_ELEV_METHOD 2
    else if ( $INLINE_MODE == "only" ) then
	setenv MRG_LAYERS_YN N
	setenv SMK_ELEV_METHOD 2
	setenv MRG_GRIDOUT_YN N
    else
	echo "SCRIPT ERROR: Unkown INLINE_MODE, $INLINE_MODE"
	$EMF_CLIENT -k $EMF_JOBKEY -m "SCRIPT ERROR: Unkown INLINE_MODE, $INLINE_MODE" -t "e" 
	exit( 1 )
    endif
endif

## test if PARENT_CASE exists, if not fail
if ( ! $?PARENT_CASE) then
   echo "ERROR: must set PARENT_CASE to run this type of sensitivity script"
   exit(1)
else
   echo "Source level sensitivity PARENT CASE = $PARENT_CASE"
endif

## store original CASE and PARENT_CASE and set CASE_CON to original CASE
set case_orig=$CASE
setenv CASE_CON $case_orig
set parent_case_orig=$PARENT_CASE

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

setenv SUBSECT $SECTOR                   # set variable for input/output names
setenv SRCABBR $SUBSECT                  # set abbreviation for naming log files

source $ASSIGNS_FILE

# List of all the helper scripts that are run in this script
set emf_cleanup  = $SCRIPTS/run/emf_cleanup.csh  # For cleanup of EMF-created scripts and stdout logs
set set_months   = $SCRIPTS/run/set_months_v4.csh  # Maybe needed for monthly projections?
set timetracker  = $SCRIPTS/run/timetracker_v2.csh  # Needed - creates TIMELOG file
set combine_data = $SCRIPTS/run/combine_data_v6.csh  # Creates the list files automatically
set smk_run      = $SCRIPTS/run/smk_run_v9.csh
set cntl_run     = $SCRIPTS/run/cntl_run_v6.csh      # Runs CNTLMAT and GRWINVEN
set qa_run       = $SCRIPTS/run/qa_run_v10.csh
set m3stat       = $SCRIPTS/run/m3stat_chk_v6.csh
set m3stat_merge = $SCRIPTS/run/m3stat_smkmerge.csh
set set_days     = $SCRIPTS/run/set_days_v5.csh  # Not needed
set log_analyzer = $SCRIPTS/log_analyzer/log_analyzer.py  # checks log files for errors/warnings
set msg_list     = $SCRIPTS/log_analyzer/known_messages.txt  # list of warnings with wildcards and rank
set path_parser  = $SCRIPTS/run/path_parser.py
set chk_smkmerge = $SCRIPTS/run/chk_smkmerge.csh
set clean_smkmerge = $SCRIPTS/run/clean_smkmerge.csh
set python_annl  = $SCRIPTS/annual_report/annual_report_v2.py
set sectorlist_parser  = $SCRIPTS/run/sectorlist_parser.py

## Using CUSTOM GCNTL (i.e. as an explicit input
setenv CUSTOM_GCNTL Y

mkdir -p $REPOUT/custom/


## If running from EMF, move old EMF-created scripts to "old"
if ( $?EMF_JOBID ) then
   source $emf_cleanup
   if ( $status != 0 ) then
	echo "ERROR: running EMF script/log cleanup script"
	$EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: running EMF script/log cleanup script" -t "e" -x $emf_cleanup
	exit( 1 )
   endif
endif


## Invoke script to interpret rest of calling arguments
$EMF_CLIENT -k $EMF_JOBKEY -m "Running set_months" -x $set_months  ## log w/ EMF server
set exitstat = 0
switch ( $#argv )
   case 6:
      source $set_months $argv[5] "$argv[6]"
      if ( $status != 0 ) set exitstat = 1
   breaksw
   case 7: 
      source $set_months $argv[5] "$argv[6]" $argv[7]
      if ( $status != 0 ) set exitstat = 1
   breaksw
   case 8:
      source $set_months $argv[5] "$argv[6]" $argv[7]
      if ( $status != 0 ) set exitstat = 1
      setenv TLABEL $argv[8]
   breaksw
endsw
if ( $exitstat != 0 ) then
    echo "ERROR: setting months"
    $EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: setting months" -t "e" -x $set_months  
    exit (1)
endif

## store set months arguments
set setmonths_flag = $argv[5]  # -m or -q
set setmonths_monlist = "$argv[6]"  # month list or quarter list

# Set spinup duration - the set_months will have QA'f the $argv[7] value
if ( $#argv > 6 ) setenv SPINUP_DURATION $argv[7]

# Save spinup array from set_months
set spinup_array = ( $SPINUP_LIST )


## Set naming label
set namelabel = ${SECTOR}_${CASE}_${GRID}
if ( $?TLABEL ) then
  set namelabel = ${namelabel}_$TLABEL
endif
if ( $?JOB_GROUP != 0 ) then
    set namelabel = ${namelabel}_$JOB_GROUP
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
echo $m3stat_merge >> $LOGS/helper_scripts_list$suffix
echo $set_days >> $LOGS/helper_scripts_list$suffix
echo $log_analyzer >> $LOGS/helper_scripts_list$suffix
echo $msg_list >> $LOGS/helper_scripts_list$suffix
echo $path_parser >> $LOGS/helper_scripts_list$suffix
echo $chk_smkmerge >> $LOGS/helper_scripts_list$suffix
echo $clean_smkmerge >> $LOGS/helper_scripts_list$suffix
echo $python_annl >> $LOGS/helper_scripts_list$suffix
echo $sectorlist_parser >> $LOGS/helper_scripts_list$suffix


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

## Check if using MRG_BYDAY flag, typically for annual pt sectors
## set a flag and save the original setting
if (! $?MRG_BYDAY ) then
    set mrgbyday_flag = N
else
    set mrgbyday_flag = Y
    set mrg_byday_sav = $MRG_BYDAY
endif

####################################################################
## Run smkmerge with intermediary from PARENT_CASE and
## include new projection matrix
####################################################################


## Check that environment variables for required inputs are set and exist
if ( ! $?SECTORLIST ) then
   echo "SCRIPT ERROR: Environment variable SECTORLIST is not defined"
   echo "              but is required by smk_sectormerge_emf.csh"
   set exitstat = 1
 
   if ( ! -e $SECTORLIST ) then
      echo "SCRIPT ERROR: SECTORLIST file is not found:"
      echo "              "$SECTORLIST
      $EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: SECTORLIST file is not found" -t "e"
      exit( 1 )
   endif
endif

## Parse SECTORLIST, determine the case for this particular sector
set sector_case = `$sectorlist_parser $SECTORLIST`
if ( $status != 0 ) then
    echo "ERROR: running sectorlist_parser"
    $EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: running sectorlist_parser on $SECTORLIST" -t "e" -x $sectorlist_parser
    exit ( 1 )
endif

## if case from SECTORLIST is same as PARENT_CASE or is 
## set to environment variable $CASE, retain PARENT_CASE
## otherwise, replace PARENT_CASE with case from SECTORLIST
if ( $sector_case == $PARENT_CASE || $sector_case == '$CASE' ) then
    echo "Retaining PARENT_CASE as $PARENT_CASE"
    $EMF_CLIENT -k $EMF_JOBKEY -m "Retaining PARENT_CASE as $PARENT_CASE"
else
    setenv PARENT_CASE $sector_case
    echo "Resetting PARENT_CASE to $PARENT_CASE"
    $EMF_CLIENT -k $EMF_JOBKEY -m "Resetting PARENT_CASE to $PARENT_CASE"

    ## reset the merge type and holiday
    set sect_mrgapproach =  `$sectorlist_parser -c mrgapproach $SECTORLIST`
    setenv M_TYPE `echo $sect_mrgapproach |cut -d "_" -f 1`
    setenv RUN_HOLIDAYS `echo $sect_mrgapproach |cut -d "_" -f 2`
    echo "Resetting M_TYPE to $M_TYPE"
    echo "Resetting Run Holidays to $RUN_HOLIDAYS"


endif

set monname = ( jan feb mar apr may jun jul aug sep oct nov dec )

# Loop through months as determined from calling arguments.
set mc = 0
set diff = 0
set g_stdate_all = $G_STDATE

## set what merge programs and related programs to run
setenv RUN_SMKMERGE  $run_smkmerge
setenv RUN_SMK2EMIS  $run_smk2emis
setenv RUN_M3STAT    $run_m3stat  
setenv RUN_PART4 Y
setenv RUN_PART2 N
#setenv RUN_SMKREPORT $run_smkreport

## This is necessary so that the assigns file doesn't delete the Temporal files
setenv RUN_TEMPORAL N

#set run_monthly = N
set first_time = Y
foreach mon ( $MONTHS_LIST ) 

  echo "test 407"

   @ mc = $mc + 1   # month array counter

   setenv MONTH ${monname[$mon]}           # set variable for month name

   ## Determine dates to run in this month
   setenv MONTH_ARRAY  $mon     # MONTH_ARRAY can have as many months listed as needed
   setenv SPINUP_ARRAY $spinup_array[$mc]

   echo "test 417"

   ## First time through and for monthly runs
   if ( $first_time == Y || $run_monthly == Y ) then
	if ( $run_monthly == Y ) then
	    ## Append month to SUBSECT, used in all intermediate file names
	    setenv SUBSECT ${SECTOR}_${MONTH}         # set variable for input/output names
	    setenv EMF_PERIOD ${MONTH}_$YEAR
	    
	    ## add month to namelabel
	    set namelabel2 = ${MONTH}_$namelabel
	else
	## Set up scripting environment variables prior to calling the Assigns file
	    setenv SUBSECT $SECTOR                   # set variable for input/output names
	    setenv EMF_PERIOD $YEAR

	    ## set namelabel2 to same as namelabel
	    set namelabel2 = $namelabel
	endif
	
	setenv SRCABBR $SUBSECT                  # set abbreviation for naming log files

	#setenv QA_TYPE  custom                # Used to name the report inputs and outputs
	setenv QA_TYPE  proj                  # Used to name the report inputs and outputs
	setenv QA_LABEL ${SUBSECT}_preGC           # Used to name the report inputs and outputs
	setenv REPLABEL $SUBSECT           # Used internally by Smkreport


	####################################################################
	## Run cntlmat to generate Projection matrix. Optionaly run grwinven
	####################################################################
	setenv RUN_PART1 N

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

        echo "test 475"

	## Register the output w/ the EMF
	echo "SCRIPT NOTE: Registering control projection report"
	## determine environmental variable based on source
	if ( $SMK_SOURCE == A ) then
	    set projrep = $APROJREP
	else if ( $SMK_SOURCE == P ) then
	    set projrep = $PPROJREP
	else if ( $SMK_SOURCE == M ) then
	    set projrep = $MPROJREP
	else
	    echo "ERROR: Unknown smoke source $SMK_SOURCE"
	    set exitstat = 1
	    goto end_of_script
	endif

#	$EMF_CLIENT -k $EMF_JOBKEY \
#	    -F  $projrep -T "Cntlmat Projection Report (External)" \
#	    -O "control projection report ${namelabel2}"    

	## turn off first time flag
	set first_time = N
   endif    ## end of cntlmat

   echo "test 500"

   ## Now setup for smkmerge from PARENT_CASE
   ## set what merge programs and related programs to run
   setenv RUN_SMKMERGE  $run_smkmerge
   setenv RUN_SMK2EMIS  $run_smk2emis
   setenv RUN_M3STAT    $run_m3stat  
   setenv RUN_PART4 Y
   setenv RUN_PART2 N

   ## Set the EMF_PERIOD for this month and year
   setenv EMF_PERIOD "${MONTH}_${YEAR}"

   # Call EMF Client for current period 
   $EMF_CLIENT -k $EMF_JOBKEY -m "Running SMKMERGE for month $MONTH" -p $EMF_PERIOD   ## log w/ EMF server

   ## Run Smkmerge and Smkreport, and optionally, Smk2emis (for CAMx)

   setenv T_TYPE $M_TYPE                 # Set temporal type to type for merge
   setenv SMK_RUN_DATES $SMK_RUN_DATES_2 # Set with file name using MTYPE

   source $set_days   # Call script to set dates for run

   # Check status of run to see if it worked. Give error if failed
   if ( $status != 0 ) then
      echo "ERROR: Running set_days in $EMF_PERIOD"
      $EMF_CLIENT - k $EMF_JOBKEY -m "ERROR: Running set_days" -t "e" -x $set_days -p $EMF_PERIOD
      exit (1)
   endif

   echo "test 530"

   ## Determine the number of days to run
   set ndays = `cat $SMK_RUN_DATES | wc -l`
   echo "Number of dates to run for $MONTH" : $ndays

   set n = 0
   set diff = 0
   set g_stdate_sav = $g_stdate_all

   # Loop through days to run during the month.
   while ( $n < $ndays )

     echo "test 543"

      @ n = $n + 1   
      
      set line = `head -$n $SMK_RUN_DATES | tail -1`
      @ diff = $line[1] - $g_stdate_sav

      setenv G_STDATE_ADVANCE $diff
      source $ASSIGNS_FILE   # Invoke Assigns file to set new dates

      ## Set EMF_PERIOD for this day
      setenv EMF_PERIOD $ESDATE
      
      echo "test 556"

      echo "Running smkmerge for $ESDATE"
      
      ## if MRG_BYDAY, have to setup representative dates
      if ( $mrgbyday_flag == Y ) then
	echo "Running Skmerge w/ merge-by-dat flag"
	## Determine day to use for merging
	set day_line = `grep ^$G_STDATE $SMK_RUN_DATES`
	set holiday = `echo $day_line | cut -c18`  

	## Setup for iterating the date in the Assigns file.
	set mon2 = `echo $day_line[2] | cut -c 5-6`
      
	# For holidays, run that holiday specifically with the prepared PTMP file
	if ( $holiday == H && $RUN_HOLIDAYS == Y ) then
	    # For December 2001 holidays, use 2002 temporals. Smkmerge ordinarily crashes because
	    # the PLAY and PTMP times are different. However, if we fake "MRG_BYDAY", it ignores that issue
	    if ( $G_STDATE < ${YEAR}000 ) then # Use 2002 temporals for 2001 Christmas holidays
		set g_stdate_temp = $G_STDATE
		@ g_stdate_temp = $g_stdate_temp + 1000 # Change 2001 to 2002
		## Account for leap years in calculaton of julian date
		if ( ${YEAR} % 4 == 1 ) then  # spinup year is leap year, but modeling year is not
		    @ g_stdate_temp = $g_stdate_temp - 1
		endif
		if ( ${YEAR} % 4 == 0 ) then  # modeling year is leap year, but spinup year is not
		    @ g_stdate_temp = $g_stdate_temp + 1
		endif
		setenv PTMP_MON ${PTMPNAME}${g_stdate_temp}.ncf
		setenv PTMP_TUE  $PTMP_MON
		setenv PTMP_WED  $PTMP_MON
		setenv PTMP_THU  $PTMP_MON
		setenv PTMP_FRI  $PTMP_MON
		setenv PTMP_SAT  $PTMP_MON
		setenv PTMP_SUN  $PTMP_MON
	    else # All 2002 holidays
		set mrg_byday_sav = $MRG_BYDAY
		unsetenv MRG_BYDAY 
	    endif

	# For non-holidays, run using the merge-by-day approach
	else
	    # Determine which days to use for each day of the week, depending on L_TYPE setting
	    if ( $L_TYPE == aveday ) then # For aveday, every day is Tuesday
		set mon_date = `grep -v " H " $SMK_RUN_DATES_1 | grep Tuesday  | grep " $YEAR$mon2" | cut -c1-7`
		set tue_date = `grep -v " H " $SMK_RUN_DATES_1 | grep Tuesday  | grep " $YEAR$mon2" | cut -c1-7`
		set sat_date = `grep -v " H " $SMK_RUN_DATES_1 | grep Tuesday  | grep " $YEAR$mon2" | cut -c1-7`
		set sun_date = `grep -v " H " $SMK_RUN_DATES_1 | grep Tuesday  | grep " $YEAR$mon2" | cut -c1-7`
	    else # default is mwdss
		set mon_date = `grep -v " H " $SMK_RUN_DATES_1 | grep Monday   | grep " $YEAR$mon2" | cut -c1-7`
		set tue_date = `grep -v " H " $SMK_RUN_DATES_1 | grep Tuesday  | grep " $YEAR$mon2" | cut -c1-7`
		set sat_date = `grep -v " H " $SMK_RUN_DATES_1 | grep Saturday | grep " $YEAR$mon2" | cut -c1-7`
		set sun_date = `grep -v " H " $SMK_RUN_DATES_1 | grep Sunday   | grep " $YEAR$mon2" | cut -c1-7`
	    endif
	 
	    echo "MRG_BYDAY mondate is $mon_date"
	    echo "MRG_BYDAY tuedate is $tue_date"
	    echo "MRG_BYDAY satdate is $sat_date"
	    echo "MRG_BYDAY sundate is $sun_date"

	    # Set PTMP files for BYDAY approach
	    setenv PTMP_MON  ${PTMPNAME}$mon_date.ncf
	    setenv PTMP_TUE  ${PTMPNAME}$tue_date.ncf
	    setenv PTMP_WED  $PTMP_TUE
	    setenv PTMP_THU  $PTMP_TUE
	    setenv PTMP_FRI  $PTMP_TUE
	    setenv PTMP_SAT  ${PTMPNAME}$sat_date.ncf
	    setenv PTMP_SUN  ${PTMPNAME}$sun_date.ncf

	endif


      endif  ## end MRG_BYDAY flag


      ## set appropriate control matrix to the projection matrix
      ## depends on the merge source
      if ( $MRG_SOURCE == A ) then
	    setenv ACMAT $APMAT
      else if ( $MRG_SOURCE == P ) then
	    setenv PCMAT $PPMAT
      else if ( $MRG_SOURCE == M ) then
	    setenv MCMAT $MPMAT
      else
	    echo "ERROR: Unknown merge source $MRG_SOURCE"
	    set exitstat = 1
	    goto end_of_script
      endif

      ## Test for necessary intermediary files
      #  C.Allen: Skip this since it doesn't work when intermediate files are split into multiple components
#      $chk_smkmerge        
#      set unzip_status = $status
#      ## status = 99 means unzipped some intermediary files (still successful return)
#      if ( $unzip_status != 0 && $unzip_status != 99 ) then	    
#	    echo "ERROR: Missing intermediary files for smkmerge for $ESDATE"
#	    $EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: Missing intermediary files for smkmerge for $ESDATE" -t "e" -x $chk_smkmerge -p $EMF_PERIOD 
#            set exitstat = 1
#	    goto end_of_script
#      endif

      # Run programs part 4
      source $smk_run        
      if ( $status != 0 ) then	    
	    echo "ERROR: Running smk_run for part 4 for $ESDATE"
	    $EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: Running smk_run for part 4 for $ESDATE" -t "e" -x $smk_run -p $EMF_PERIOD 
            set exitstat = 1
	    goto end_of_script
      endif

      ## Run M3STAT on smkmerge file(s)
      $m3stat_merge
      if ( $status != 0 ) then	    
	    echo "ERROR: Running m3stat_smkmerge for for $ESDATE"
	    $EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: Running m3stat_smkmerge for $ESDATE" -t "e" -x $m3stat_merge -p $EMF_PERIOD 
            set exitstat = 1
	    goto end_of_script
      endif

      ## Check if need to remove unzipped intermediary files
#      if ( $unzip_status == 99 ) then
#	    $clean_smkmerge
#	    if ( $status != 0 ) then
#		echo "ERROR: Running clean_smkmerge for for $ESDATE"
#		$EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: Running clean_smkmerge for $ESDATE" -t "e" -x $clean_smkmerge -p $EMF_PERIOD 
#		set exitstat = 1
#		goto end_of_script
#	    endif
#      endif
      
      # Zip output files if ZIPOUT=Y
      # Added by C.Allen, 13 Sep 2018
      # This will only work if Smkmerge is run, otherwise the EOUT and INLN parameters aren't set
      if ( $?ZIPOUT && $?EOUT ) then
         if ( $ZIPOUT == Y && -e $EOUT ) then
            gzip -vf $EOUT
         endif
      endif
      if ( $?ZIPOUT && $?INLN ) then
         if ( $ZIPOUT == Y && -e $INLN ) then
            gzip -vf $INLN
         endif
      endif
      
      ## if merge-by-day flag, reset MRG_BYDAY
      if ( $mrgbyday_flag == Y ) then
	    setenv MRG_BYDAY $mrg_byday_sav
      endif
      
   end # End loop over days

   unsetenv G_STDATE_ADVANCE

   # 5/13/09: C. Allen moved this inside the months script; otherwise the January 1st model-ready file is deleted after it is created
   setenv RUN_PART4 N

   echo "test 683"

end  # End loop over months

### Register Smkmerge AQM-ready outputs with EMF 
#echo "SCRIPT NOTE: Registering Smkmerge outputs"
#$EMF_CLIENT -k $EMF_JOBKEY \
#    -D $INTERMED \
#    -P "emis_mole_${SECTOR}_*_${GRID}_*.ncf" \
#    -T "CMAQ Model Ready Emissions: Sector-specific (External)" \
#    -N "Model-ready CMAQ $SECTOR $CASE" \
#    -O "AQM-ready data $SECTOR"
#
### Register smkmerge state reports with EMF
#echo "SCRIPT NOTE: Registering Smkmerge state reports"
#$EMF_CLIENT -k $EMF_JOBKEY \
#    -D $REPOUT/smkmerge/$SECTOR \
#    -P "rep_mole_${SECTOR}_*_${GRID}_*.txt" \
#    -T "Smkmerge Report state (External Multifile)" \
#    -N "Smkmerge $SECTOR $CASE state reports" \
#    -O "state Smkmerge $SECTOR (External)"

## Run python annual summary of smkreports, even if less than a year
## register output w/ EMF
if ( ! $?RUN_PYTHON_ANNUAL ) then
    setenv RUN_PYTHON_ANNUAL N
endif	 


if ( $RUN_PYTHON_ANNUAL != N  && $?SECTORLIST ) then
    ## resource set months w/ arguments, but turn off any spinup
    source $set_months $setmonths_flag "$setmonths_monlist"
    if ( $status != 0 ) then
	    echo "ERROR: Running set_months for annual_report.py "
            set exitstat = 1
	    goto end_of_script
    endif

    ## Get first and last month to run
    set months_list = ($MONTHS_LIST)
    set monS = $months_list[1]
    set monE = $months_list[$#months_list]

    echo "SCRIPT NOTE: Running python annual report script for months $monS - $monE"

    ## Get root of the merge dates directory from smk merge date file
    if ( $?MRGDATE_FILES) then
	## get path from first merge date file
	echo "Getting merge dates path from MRGDATE_FILES"
	set mrgdate_root = `$path_parser $MRGDATE_FILES`
    else 
	## no merge dates file, so get from parent case's inputs/mrggrid
	echo "Getting merge dates path from PARENT_CASE"
	set mrgdate_root = $PARENT_INVDIR/mrggrid
    endif

    if (! -e $REPOUT/annual_report) mkdir -p $REPOUT/annual_report

    set repregion = "county"
    if ($MRG_REPSTA_YN == Y && $MRG_REPCNY_YN == N) set repregion = "state"

    ## run python annual report
    $python_annl -R $REPOUT/smkmerge -D $mrgdate_root \
	    --start_month=$monS --end_month=$monE \
	    -O $REPOUT/annual_report -r $repregion \
	    -o annual_${CASE}_${SECTOR}_${GRID}_${SPC} \
	    $SCRIPTS/annual_report/parameter_file_$SPC.txt
    if ( $status != 0 ) then
	    echo "ERROR: Running annual_report.py "
	    $EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: Running running annual_report" \
		-t e -t "e" -x $SCRIPTS/annual_report/annual_report.py
            set exitstat = 1
	    goto end_of_script
    endif

    ## Register the output w/ the EMF
#    echo "SCRIPT NOTE: Registering smkmerge state totals"
#    $EMF_CLIENT -k $EMF_JOBKEY \
#	-F  $REPOUT/annual_report/annual_${CASE}_${SECTOR}_${GRID}_${SPC}_emf.csv \
#        -T "Smkmerge report state annual summary (CSV)" -O "state totals ${namelabel}"    
endif

###########################################################
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
if ( $status == 10 ) then
#   echo "SCRIPT NOTE: Registering log summary, level1"
#   $EMF_CLIENT -k $EMF_JOBKEY -F $REPOUT/log_analyzer/rep_logs_${namelabel}_level1.csv \
#        -T "Log summary level 1" -O "Level 1 log summary ${namelabel}"
	
   echo "ERROR: The log analyzer has detected warning or error messages in your SMOKE logs which may indicate a problem."
   echo "*****  Please review the priority 0 and 1 messages listed in this report:"
   echo "*****  $REPOUT/log_analyzer/rep_logs_${namelabel}_level1.csv"
   $EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: Level 1 errors or Level 0 warnings found. Check log files, fix inputs, and rerun." -t "e" -x $log_analyzer
   set exitstat = 1

else if ( $status != 0 ) then
   echo "ERROR: running log_analyzer, level 1."
   $EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: running log_analyzer, level 1" -t "e" -x $log_analyzer
   set exitstat = 1

## Register log analyzer output
else
#   echo "SCRIPT NOTE: Registering log summary, level1"
#   $EMF_CLIENT -k $EMF_JOBKEY -F $REPOUT/log_analyzer/rep_logs_${namelabel}_level1.csv \
#        -T "Log summary level 1" -O "Level 1 log summary ${namelabel}"
endif

## If ran successfully, add this sector to source sector override file
## used by sectormerge
#  C.Allen 12 Mar 2018: Removed the "ran successfully" criteria because almost every run returns status=1 from log analyzer
#if ( $exitstat == 0 ) then
   echo "Adding $SECTOR to source sector override file"
   # Ensure mrrgirddirector for source sector override file is created
   if ( ! -e $IMD_ROOT/mrggrid ) then
      mkdir -p $IMD_ROOT/mrggrid
      chmod ug+rwx $IMD_ROOT/mrggrid
   endif
   
   ## append this sector to the source sector override file
   echo $SECTOR >> $IMD_ROOT/mrggrid/source_sector_override_${CASE}_${GRID}.txt
#endif
## Ending of script
#
exit( $exitstat )
