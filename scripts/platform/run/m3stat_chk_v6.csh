#!/bin/tcsh -f
#
##### run m3stat on SMKMERGE output to check for NaN and negative values
# Jonathan Petters 5/11/04
# Updated by M. Houyoux to use a script argument on 11/22/2006
#
# Modified for integrating w/ EMF:  A. Zubrow - IE UNC  August, 2007

echo 'Now running M3STAT'

switch ( $# )
   case 0:
      echo "SCRIPT ERROR: m3stat_chk_v2.csh script requires an argument that"
      echo "    gives the file name on which run the script."
      exit( 1 )
endsw

## Make sure that the EMF time period is set
if ( $?EMF_PERIOD == 0 ) then
    set emf_period = ""
else
    set emf_period = $EMF_PERIOD
endif

## If EMF_CLIENT is not defined, assume that the user has not setup
#    the assigns file to use the new EMF-ready helper scripts
if ( ! $?EMF_CLIENT ) then
   setenv EMF_CLIENT false
   echo 'NOTE: EMF_CLIENT setting assumed to be "false" in' !$
endif

## helper scripts used by this program
set timetracker = $SCRIPTS/run/timetracker_v2.csh
set checklog = $SCRIPTS/run/checklogfile.csh

setenv TOOLS $IOAPIDIR

setenv INFILE $argv[1]

## For inline approach, m3stat gets run twice (once for POUT and once for INLN).
#  For these cases, a second argument is sent that allows this script to produce
#  separate log and report files for each.
if ( $#argv >= 2 ) then
  setenv LOG_LABEL $argv[2]
endif

setenv LOGDIR  $M3STAT_LOGS
if ( ! -e $LOGDIR ) then
   mkdir -p $LOGDIR
endif

if ( $?LOG_LABEL ) then
  setenv LOGFILE $LOGDIR/m3stat_${SUBSECT}_${ESDATE}_${GRID}_${SPC}_${CASE}_${LOG_LABEL}.log
else
  setenv LOGFILE $LOGDIR/m3stat_${SUBSECT}_${ESDATE}_${GRID}_${SPC}_$CASE.log
endif
if ( -e $LOGFILE ) then
  /bin/rm $LOGFILE
endif

if ( $?LOG_LABEL ) then
   setenv REPORT $LOGDIR/m3stat_${SUBSECT}_${ESDATE}_${GRID}_${SPC}_${CASE}_${LOG_LABEL}.rpt
else
   setenv REPORT $LOGDIR/m3stat_${SUBSECT}_${ESDATE}_${GRID}_${SPC}_$CASE.rpt
endif
if ( -e $REPORT ) then
   /bin/rm -f $REPORT
endif 

set m3stat_in = $INTERMED/.m3stat_${MONTH}_$$.in 
if ( -e $m3stat_in ) then
   /bin/rm -f $m3stat_in
endif

# Make sure that the timelog is set.  Whether the output file is
#    set or not determines whether the timetracking file is created.
# Also, time tracking is not written if using DEBUGMODE = Y
set timelog_yn = N
if ( $?TIMELOG && $DEBUGMODE != Y ) then
   set timelog_yn = Y
endif

echo INFILE >> $m3stat_in
echo REPORT >> $m3stat_in
echo Y >> $m3stat_in
echo $G_STDATE >> $m3stat_in
echo 0 >> $m3stat_in
echo 250000 >> $m3stat_in

set startdt = `date +%m/%d/%Y,%T`

$TOOLS/m3stat < $m3stat_in
if ( $timelog_yn == Y ) then
   $timetracker N $TIMELOG $startdt m3stat $ESDATE
   if ( $status != 0 ) then
       echo "ERROR: Problem calling timetracker from m3stat script"
       $EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: Problem calling timetracker from m3stat script" -x $timetracker  -t "e" -p $emf_period ## log w/ EMF server
       exit ( 1 )
   endif  
endif 

setenv LOG_FILE $LOGFILE

# C. Allen update 22 Oct 2013: If a bad number is found, we want the script to error out immediately.
# So, I added "exit (1)" to each of the if/then blocks below.

# Made case sensitive so that T_NONANAL species doesn't trigger an error
grep -E "NaN|Infinity" $REPORT 
if ( $status == 0 ) then
   echo '*******************************' >> $LOGFILE
   echo '* ERROR detected in M3STAT report:' >> $LOGFILE
   echo "* NaN or Infinity found in "$REPORT >> $LOGFILE
   echo '*******************************' >> $LOGFILE
   $EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: NaN or Infinity found in $REPORT" -p $emf_period  -t "e" ## Log to EMF Server
   exit ( 1 )
endif

grep ' -' $REPORT 
if ( $status == 0 ) then
   echo '*******************************' >> $LOGFILE
   echo '* ERROR detected in M3STAT report:' >> $LOGFILE
   echo "* Negative # found in "$REPORT >> $LOGFILE
   echo '*******************************' >> $LOGFILE
   $EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: Negative # found in $REPORT" -p $emf_period -t "e" ## Log to EMF Server
   exit ( 1 )
endif

# Check log file
$checklog
if ( $status != 0 ) then
  echo SCRIPT ERROR: ERROR detected in M3STAT
  $EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: detected in M3STAT" -p $emf_period -t "e"  ## Log to EMF Server
  exit ( 1 )
endif

/bin/rm -rf $m3stat_in
   
  

    
