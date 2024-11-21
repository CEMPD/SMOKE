#!/bin/tcsh -f
# This script allows one to use the combine utility, which is most commonly used for CB6-CAMx to CB6-CMAQ conversions, via an EMF job.

## log w/ EMF server that script is running
$EMF_CLIENT -k $EMF_JOBKEY -s "Running"

set exitstat = 0

switch ( $#argv )
   case 0:
   case 1:
   case 2:
      echo "SCRIPT ERROR: Script requires arguments for the grid name"
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

if (! $?SPECIES_DEF) then
   echo "ERROR: SPECIES_DEF input is required to run combine utility."
   $EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: SPECIES_DEF input not defined" -t "e"
   exit ( 1 )
endif

# This controls the length of the m3xtract extractions.
# Default is 250000, which works on RHEL7.
# For RHEL8, this should be set to blank in the EMF case.
if (! $?M3XTRACT_LENGTH) then
   setenv M3XTRACT_LENGTH 250000
endif

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
set combine_exec = $SCRIPTS/combine/combine-terra.exe

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
echo $combine_exec >> $LOGS/helper_scripts_list$suffix

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

set outfile_created = 0

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

	if (! -e $PREMERGED/pre_combine) mkdir $PREMERGED/pre_combine

	foreach date (`tail -${days} ${smk_dat} | cut -d',' -f1`)

		setenv OUTFILE $PREMERGED/emis_mole_${SECTOR}_${date}_${GRID}_${EMF_SPC}_${CASE}.ncf
		if (-e $PREMERGED/pre_combine/emis_mole_${SECTOR}_${date}_${GRID}_${EMF_SPC}_${CASE}.ncf || -e $PREMERGED/pre_combine/emis_mole_${SECTOR}_${date}_${GRID}_${EMF_SPC}_${CASE}.ncf.gz ) then
	    	   echo "ERROR: Emissions file emis_mole_${SECTOR}_${date}_${GRID}_${EMF_SPC}_${CASE}.ncf already exists in pre_combine directory, suggesting combine utility has already been run."
		   echo "       If you are trying to rerun the utility, move the original files from the pre_combine directory to the original premerged directory and resubmit."
	    	   $EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: Running combine utility" -t "e"
		   exit ( 1 )
		endif

		if (! -e $OUTFILE && -e $OUTFILE.gz) gunzip -fv $OUTFILE
		if (-e $OUTFILE) then
		   mv $OUTFILE $PREMERGED/pre_combine/
		   setenv INFILE1 $PREMERGED/pre_combine/emis_mole_${SECTOR}_${date}_${GRID}_${EMF_SPC}_${CASE}.ncf
		else
		   echo "$OUTFILE not found for this date; moving on to next date"
		   continue
		endif

		# Call EMF Client for current period
		$EMF_CLIENT -k $EMF_JOBKEY -m "Running combine utility for $date"

		$combine_exec
		if ( $status != 0 ) then
		        echo "ERROR: Running combine utility"
		        $EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: Running combine utility for $date" -t "e"
		        exit ( 1 )
		endif
		
		# Some tools, such as PAVE, won't work on raw outputs created by combine, unless we pass them through m3xtract.
		
	        # temporary file with m3xtract instructions
	        set m3xtract_in = $INTERMED/.m3xtract.${MONTH}_$$.in
 	        if (-e $m3xtract_in) rm -f $m3xtract_in
	        setenv LOGFILE $INTERMED/.m3xtract.${MONTH}_$$.log
	        if (-e $LOGFILE) rm -f $LOGFILE
	 	
		setenv INFILE $OUTFILE
		setenv OUTFILE $INFILE.xt
		if (-e $OUTFILE) rm -f $OUTFILE
		
	        echo INFILE >> $m3xtract_in
	        echo 0 >> $m3xtract_in
	        echo -1 >> $m3xtract_in
	        echo "" >> $m3xtract_in
	        echo 0 >> $m3xtract_in
	        echo $M3XTRACT_LENGTH >> $m3xtract_in
	        echo OUTFILE >> $m3xtract_in
		
	        setenv PROMPTFLAG Y # Y needed for m3xtract
	       
	        $IOAPIDIR/m3xtract < $m3xtract_in
                if ( $status != 0 ) then
      	           echo "ERROR: running m3xtract for $date"
	           $EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: running m3xtract for $date" -t "e" -x m3xtract -p $EMF_PERIOD
                   set exitstat = 1
                   goto end_of_script
                endif
		
		if (-e $OUTFILE) set outfile_created = 1

		rm -f $INFILE $m3xtract_in $LOGFILE
		mv -fv $OUTFILE $INFILE

	end

end # months loop 

# Since script skips over days with no input file, we don't get an error if it just skips over every day.
# Check to see if ANY output files were created; if not, that's obviously bad, so trigger an error.
if ( $outfile_created == 0 ) then
        echo "ERROR: Running combine utility - no output files were created. Check standard output for error messages."
        $EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: Running combine utility" -t "e"
        exit ( 1 )
endif

#########################
## EMF script clean-up
#########################

# I don't think we want to add to the source_sector_override here, but I'm keeping the code - just commented out - in case we do need it

### Check if RUN_SOURCESENS exists, if not set to N
#if ( $?RUN_SOURCESENS == 0 ) then
#    setenv RUN_SOURCESENS N
#endif
#
### If running source sensitivity
### If ran successfully, add this sector to source sector override file
### used by sectormerge
#if ( $exitstat == 0  && $RUN_SOURCESENS == Y ) then
#   echo "Adding $SECTOR to source sector override file"
#   # Ensure mrrgirddirector for source sector override file is created
#   if ( ! -e $IMD_ROOT/mrggrid ) then
#      mkdir -p $IMD_ROOT/mrggrid
#      chmod ug+rwx $IMD_ROOT/mrggrid
#   endif
#
#   ## append this sector to the source sector override file
#   echo ${SECTOR} >> $IMD_ROOT/mrggrid/source_sector_override_${CASE}_${GRID}.txt
#endif

# Label for the end of the script, used during script abort
# No log analyzer here
end_of_script:

## Ending of script
#
exit( $exitstat )
