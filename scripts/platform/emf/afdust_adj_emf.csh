#!/bin/tcsh -f
# This script will do the complete afdust adjustment
# Updated 28/2/2012 J. Beidler
# Updated 26 Mar 2013 J. Beidler - moved to new version of QA script

## log w/ EMF server that script is running
$EMF_CLIENT -k $EMF_JOBKEY -s "Running"

set exitstat = 0

switch ( $#argv )
   case 0:
   case 1:
   case 2:
      echo "SCRIPT ERROR: Script requires arguments for a grid name"
      echo "              and the -m or -q option with 3 settings."
      echo " "
      echo "  This script expects to be called using one of the following argument lists:"
      echo "     <grid abbrv> -m <monthlist> <spinup>"
      echo "     <grid abbrv> -q <quarters> <spinup>"
      echo " "
      echo "  You can either use one approach or the other (differing by the -m or -q options)."
      echo " "
      echo "  In the above list, the arguments are defined as follows:"
      echo "     <grid abbrv>       : Grid abbreviation (e.g., 36US1)"
      echo "     <monthlist>        : list of months to run when using the -m option"
      echo "     <quarters>         : list of quarters to run when using the -q option"
      echo "     <spinup>           : set to number of days between 1 and 20 to run a spinup"
      echo "                          period (value sets number of days), and N otherwise"
      echo " "
      echo "  Examples:"
      echo "     <script name> 36US1 -m '1 2 3' 0"
      echo "              This example runs the script for Jan, Feb, & Mar"
      echo "              for the 36US1 grid, with no spinup days and"
      echo " "
      echo "     <script name> 12EUS1 -q '2 3' 10"
      echo "               This example runs the script for the 2nd & 3rd quarters,"
      echo "               for the 12EUS1 grid, including 10 spin-up days, and gives"
      $EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: smoke script did not receive more than 2 arguments" -t "e"
      exit( 1 )
endsw

setenv GRID "$argv[1]"

########################
## Standard EMF script calls
## Modified from M. Houyoux's EMF scripts
########################

## source the ASSIGN file
source $ASSIGNS_FILE

## List of all the helper scripts that are run in this script
set emf_cleanup  = $SCRIPTS/run/emf_cleanup.csh
set set_months   = $SCRIPTS/run/set_months_v4.csh
set timetracker  = $SCRIPTS/run/timetracker_v2.csh
set set_days     = $SCRIPTS/run/set_days_v5.csh
set log_analyzer = $SCRIPTS/log_analyzer/log_analyzer.py
set msg_list     = $SCRIPTS/log_analyzer/known_messages.txt
set duplic_chk   = $SCRIPTS/run/duplicate_check.csh
set path_parser  = $SCRIPTS/run/path_parser.py
set mult         = $SCRIPTS/afdust_adj/mult.x
set precip_adj   = $SCRIPTS/afdust_adj/apply_precip_adj_wrf.x
set ann_sum      = $SCRIPTS/afdust_adj/afdust_ann_report.py

## If running from EMF, move old EMF-created scripts to "old"
if ( $?EMF_JOBID ) then
   source $emf_cleanup 
   if ( $status != 0 ) then
        echo "ERROR: running EMF script/log cleanup script"
        $EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: running EMF script/log cleanup script" -t "e" -x $emf_cleanup
        exit ( 1 )
   endif
endif

## Invoke script to interpret calling arguments
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
endsw
if ( $exitstat != 0 ) then
    echo "ERROR: setting months"
    $EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: setting months" -t "e" -x $set_months
    exit ( 1 )
endif

# Set spinup duration - the set_months will have QA'f the $argv[4] value
if ( $#argv > 3 ) setenv SPINUP_DURATION $argv[4]

# Save spinup array from set_months
set spinup_array = ( $SPINUP_LIST )

## Set naming label
set namelabel = ${SECTOR}_${CASE}_${GRID}

## Set the EMF_PERIOD to the year
setenv EMF_PERIOD $YEAR

## Set WRF_VERSION variable, used by the precip_adj program to determine
#  the appropriate soil types (new soil type added in WRF v4.1).
#  If WRF_VERSION = 4, use the newer soil types.
#  If WRF_VERSION = anything else, use the older soil types.
#  If WRF_VERSION undefined, attempt to determine the WRF version from the path
#    of the MET_ROOT. 
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
      echo "ERROR: WRF_VERSION must be set so the precip adjustment script knows which set of SLTYP values to use"
      echo "If using WRFv4.1 or later, set WRF_VERSION = 4"
      echo "If using an earlier version of WRF, set WRF_VERSION = 3"
      $EMF_CLIENT -k $EMF_JOBKEY -m "WRF_VERSION must be set so the precip adjustment script knows which set of SLTYP values to use" -t "e"
      exit (1)
   endif
   echo "WRF_VERSION = $WRF_VERSION"
endif

## Record the helper scripts being used
set suffix = _$namelabel.txt
echo "# Helper scripts used for $SECTOR" > $LOGS/helper_scripts_list$suffix
echo $emf_cleanup >> $LOGS/helper_scripts_list$suffix
echo $set_months >> $LOGS/helper_scripts_list$suffix
echo $timetracker >> $LOGS/helper_scripts_list$suffix
echo $set_days >> $LOGS/helper_scripts_list$suffix
echo $log_analyzer >> $LOGS/helper_scripts_list$suffix
echo $msg_list >> $LOGS/helper_scripts_list$suffix
echo $duplic_chk >> $LOGS/helper_scripts_list$suffix
echo $path_parser >> $LOGS/helper_scripts_list$suffix
echo $mult        >> $LOGS/helper_scripts_list$suffix
echo $precip_adj  >> $LOGS/helper_scripts_list$suffix
echo $ann_sum     >> $LOGS/helper_scripts_list$suffix
#echo $adj_rep     >> $LOGS/helper_scripts_list$suffix

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

###################################################

# Set the path to the smoke merge dates directory
set mrgdate_root = `$path_parser $MRGDATE_FILES` 
setenv MRGDATE_ROOT $mrgdate_root

set orig_path = $PREMERGED
set adj_path = ${PREMERGED}_adj

# Create adjustment path if it does not exist
if ( ! -d ${adj_path} ) then
	mkdir ${adj_path}
endif

set months = ($MONTHS_LIST)

echo "months = $MONTHS_LIST"

## Set up the spinups for this script
# C. Allen (30 Aug 2012): Rewrote this section. If SPINUP_DURATION > 0,
#   then the spinup month is already first in MONTHS_LIST. Basically,
#   we know to use a spinup year if the first month = 12 and SPINUP_DURATION > 0.
set s_year = $BASE_YEAR
set spin = $SPINUP_DURATION
if ($spin > 0 && $months[1] == 12) then
        @ s_year = $BASE_YEAR - 1
endif

set firstiteration = 1

# Loop over the months
foreach mon ( $months )



	# For the spinup month set the days from that month and only with the number of spin days 
	# C. Allen (30 Aug 2012): Rewrote this section. If this is the first month ($firstiteration = 1)
	#   and December spinup was detected above ($s_year != $BASE_YEAR) then load the spinup dates file
	#   this time. Then set firstiteration = 0 so that we don't use the spinup year for the rest of the run.
	#   Can't use "month = months[1]" as a criterion instead of $firstiteration because December could 
	#   come up twice, once at the start of the run for spinup, and then again at the end of the base year.
	if ($s_year != $BASE_YEAR && $firstiteration) then
		if ($mon < 10) set mon = 0${mon}
		set smk_dat = ${mrgdate_root}/smk_merge_dates_${s_year}${mon}.txt
		set days = $spin
		set firstiteration = 0

	# Set smk_merge_dates file and days in the month being run
	else
		if ($mon < 10) set mon = 0${mon}
		set smk_dat = ${mrgdate_root}/smk_merge_dates_${BASE_YEAR}${mon}.txt
		set days = `wc -l ${smk_dat} | cut -d' ' -f1`
		@ days = $days - 1
	endif

        echo "smk_dat = $smk_dat"
	echo "days = $days"

	#### this script will apply the xportfraction to a sector-processed set of emissions
	#   files using IOAPI programs
	#
	#   Script written by G. Pouliot:  /home/pou/afdust/new_meth/apply_xportfrac.csh
	#   Script edited by R. Cleary 05/20/2010
	#
	#####################################################################################

	setenv INFILE1 $XPORTFRAC

	# C. Allen (25 Nov 2014): Script was hardwired to process holidays, which we don't want to do for othafdust.
	set datecolumn = 7
	if ($RUN_HOLIDAYS == N) set datecolumn = 6

	# Get the date information for the last x days in the file.  Typically this is the number of days in the month.
	foreach date (`tail -${days} ${smk_dat} | cut -d',' -f${datecolumn}`)
	
	        echo "date = $date"

		setenv INFILE2 ${orig_path}/emis_mole_${SECTOR}_${date}_${GRID}_${EMF_SPC}_${CASE}.ncf
		setenv OUTFILE ${orig_path}/emis_mole_${SECTOR}_${date}_${GRID}_${EMF_SPC}_${CASE}_xportfrac.ncf
		setenv LOGFILE ${INTERMED}/logs/xportfrac_${date}_${GRID}_${EMF_SPC}_${CASE}.log

		# 7 Jan 2022: Christine Allen edited so that OUTFILE and LOGFILE are deleted first if they already exist
		if ( -e $LOGFILE ) rm -f $LOGFILE
		if ( -e $OUTFILE ) rm -f $OUTFILE

		# Call EMF Client for current period
		$EMF_CLIENT -k $EMF_JOBKEY -m "Running xportfrac adjustments"
		if (-e $LOGFILE) rm -f $LOGFILE

		$mult

		if ( $status != 0 ) then
		        echo "ERROR: Running xportfrac adjustments "
		        $EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: Running xportfrac adjustments" -t "e"
		        exit ( 1 )
		endif

	end

	#### 
	#   this script will apply a precipitation adjustment to the afdust emissions
	#
	#   Script written by G. Pouliot:  /home/pou/codes/afdust/src/apply_adjustment_jun_2002.csh
	#   Script edited by R. Cleary 12/17/2010
	#   Script modified for terra 03/07/2011
	#   Script modified for D31 project
	#
	#####################################################################################

	foreach date (`tail -${days} ${smk_dat} | cut -d',' -f1`)

		set outdate = ${date}
		set indate = `grep ^${date} ${smk_dat}| cut -d',' -f${datecolumn}`

		set metdate = `echo ${date} | cut -c3-8`
		setenv METCRO2D ${MET_ROOT}/METCRO2D_${metdate}

#		if ( ! -d ${adj_path}/logs ) mkdir ${adj_path}/logs
		setenv LOGFILE ${INTERMED}/logs/apply.precip.adj.${EMF_SPC}.${GRID}.${CASE}.${outdate}.log

		setenv INFILE ${orig_path}/emis_mole_${SECTOR}_${indate}_${GRID}_${EMF_SPC}_${CASE}_xportfrac.ncf
		setenv OUTFILE ${adj_path}/emis_mole_${SECTOR}_adj_${outdate}_${GRID}_${EMF_SPC}_${CASE}.ncf

		if ( -e $LOGFILE ) rm -f $LOGFILE
		if ( -e $OUTFILE ) rm -f $OUTFILE

		# Call EMF Client for current period
		$EMF_CLIENT -k $EMF_JOBKEY -m "Running precip adjustments"

		$precip_adj

		if ( $status != 0 ) then
		        echo "ERROR: Running precip adjustments "
		        $EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: Running precip adjustments" -t "e"
		        exit ( 1 )
		endif
	end

end # months loop 

#########################
## QA
#########################
### Run the onroad annual summary for temperature unadjusted onroad emissions.
#setenv RUN_SECT "UNADJ"

# Call EMF Client for current period
#$EMF_CLIENT -k $EMF_JOBKEY -m "Running annual summary for unadjusted emissions"    ## log w/ EMF server

#/usr/bin/python /garnet/oaqps/smoke/test/smoke3.0/scripts/afdust_adj/afdust_ann_sum.py
#if ( $status != 0 ) then
#	echo "ERROR: Running annual summary for unadjusted emissions"
#	$EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: Running annual summary for unadjusted emissions" -t "e" -p $EMF_PERIOD
#	exit ( 1 )
#endif

### Run the onroad annual summary for temperature adjusted onroad emissions.
setenv RUN_SECT "ADJ"

# Call EMF Client for current period
$EMF_CLIENT -k $EMF_JOBKEY -m "Running annual summary for adjusted emissions"    ## log w/ EMF server

$ann_sum

if ( $status != 0 ) then
	echo "ERROR: Running annual summary for adjusted emissions"
	$EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: Running annual summary for adjusted emissions" -t "e" -p $EMF_PERIOD
	exit ( 1 )
endif

#########################
## EMF script clean-up
#########################

## Check if RUN_SOURCESENS exists, if not set to N
if ( $?RUN_SOURCESENS == 0 ) then
    setenv RUN_SOURCESENS N
endif

## If running source sensitivity
## If ran successfully, add this sector to source sector override file
## used by sectormerge
if ( $exitstat == 0  && $RUN_SOURCESENS == Y ) then
   echo "Adding $SECTOR to source sector override file"
   # Ensure mrrgirddirector for source sector override file is created
   if ( ! -e $IMD_ROOT/mrggrid ) then
      mkdir -p $IMD_ROOT/mrggrid
      chmod ug+rwx $IMD_ROOT/mrggrid
   endif

   ## append this sector /w "adj" appended to the source sector override file
   echo ${SECTOR}_adj >> $IMD_ROOT/mrggrid/source_sector_override_${CASE}_${GRID}.txt
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
