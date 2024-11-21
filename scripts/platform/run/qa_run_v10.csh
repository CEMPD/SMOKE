#!/bin/csh -f
#
# $Version$
# $Path$
# $Date: 2003/03/14 00:30:52 $
#
# This script runs the SMOKE QA processors.  
#
# Script created by : M. Houyoux
# Last edited : January, 2008
# Modified March 2017 for AERMOD support.
#*********************************************************************

# If smkreport is not to be run, abort
if ( $RUN_SMKREPORT != Y ) then
   exit( 0 )
endif

# Create directory for output program logs

setenv OUTLOG $LOGS
if ( ! -e $OUTLOG ) then
   mkdir -p $OUTLOG
   chmod ug+w $OUTLOG
endif

# Initialize exit status
set exitstat = 0

# Make sure that debug mode and debug executable are set
if ( $?DEBUGMODE ) then
   set debugmode = $DEBUGMODE
else
   set debugmode = N
endif

if ( $?DEBUG_LOCATION ) then
   set debug_exe = $DEBUG_LOCATION
else
   if ( $debugmode == Y ) then
      echo "*** ERROR: DEBUG_LOCATION must be defined and point to debug executable location if DEBUGMODE=Y ***"
      $EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: DEBUG_LOCATION must be defined and point to debug executable location if DEBUGMODE=Y" -t "e" -p $emf_period ## log w/ EMF server
      exit ( 1 )
   endif
endif

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
set movelog = $SCRIPTS/run/movelog.csh

# Make sure that the timelog is set.  Whether the output file is
#    set or not determines whether the timetracking file is created.
# Also, time tracking is not written if using DEBUGMODE = Y
set timelog_yn = N
if ( $?TIMELOG && $DEBUGMODE != Y ) then
   set timelog_yn = Y
endif

# Ensure that SMK_SOURCE is defined
if ( $?SMK_SOURCE ) then

# SMK_SOURCE is not defined   
else
   echo 'SCRIPT ERROR: Environment variable SMK_SOURCE is not set,'
   echo '              but it is needed to use the qa_run.csh script!'
   set exitstat = 1

endif

# Check if QA_TYPE variable is set
if ( $?QA_TYPE ) then
else

   echo 'SCRIPT ERROR: Environment variable QA_TYPE is not set,'
   echo '              but it is needed to use the qa_run.csh script!'
   echo '              Valid values are "inventory" or "monthly".'
   set exitstat = 1
   
endif

# Abort if already had an error
if ( $exitstat == 1 ) then
   $EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: detected in qa_run" -p $emf_period -t "e" ## log w/ EMF server
   exit( $exitstat )
endif

# Check if QA label is set, and initialize file
# If QA_LABEL is set, use it for naming reports using "label" variable
if ( $?QA_LABEL ) then

   # For any part of processing that is day-specific, add date in label
   # 7 Feb 2019: Disabled adding the date to the label for PART2 to suppoer emissions modeling workgroup report by temporal profile
   if ( $RUN_PART3 == Y || $RUN_PART4 == Y ) then
      set label = ${QA_LABEL}_${CASE}_${ESDATE}
      set runperiod = $ESDATE
   else
      set label = ${QA_LABEL}_${CASE}
      set runperiod = inv
   endif

# If QA_LABEL is not set, set report "label" variable to defaults
else

   # For any part of processing that is day-specific, add date in label
   # 7 Feb 2019: Disabled adding the date to the label for PART2 to suppoer emissions modeling workgroup report by temporal profile
   if ( $RUN_PART3 == Y || $RUN_PART4 == Y ) then
      set label = ${SECTOR}_${CASE}_${ESDATE}
      set runperiod = $ESDATE
   else
      set label = ${SECTOR}_${CASE}
      set runperiod = inv
   endif

endif

if ( $?QA_TYPE ) then
   set qa_type = $QA_TYPE

   if ( $qa_type == none ) then 
      echo 'SCRIPT NOTE: "QA_TYPE" is set to "none", so no QA will be run.'
      $EMF_CLIENT -k $EMF_JOBKEY -m "QA_TYPE is set to none, so no QA will be run" -p $emf_period  ## Log to EMF Server
      exit( $exitstat )
   endif

else
   set qa_type = none
   echo 'SCRIPT NOTE: "QA_TYPE" is not set, so no QA will be run.'
   $EMF_CLIENT -k $EMF_JOBKEY -m "QA_TYPE is not set, so no QA will be run" -p $emf_period   ## Log to EMF Server
    
endif

# Abort if already had an error
if ( $exitstat == 1 ) then
   $EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: detected in qa_run" -p $emf_period -t "e" ## Log to EMF Server
   exit( $exitstat )

endif

# Set output directory and make sure that it's there & set permissions
if ( $qa_type != custom ) then
   set outdir = $REPOUT/$qa_type
endif

## based on qa_type, set up the repconfig file for reports
if ( $qa_type == inv ) then
    if ( $?REPCONFIG_INV ) then
	setenv REPCONFIG $REPCONFIG_INV
    endif
endif
if ( $qa_type == proj ) then
    set outdir = $REPOUT/proj
    if ( $?REPCONFIG_INV ) then
        setenv REPCONFIG $REPCONFIG_INV
    endif
endif
if ( $qa_type == inv2 ) then
    set outdir = $REPOUT/inv
    set qa_type = inv
    if ( $?REPCONFIG_INV2 ) then
	setenv REPCONFIG $REPCONFIG_INV2
    endif
endif
if ( $qa_type == inv3 ) then
    set outdir = $REPOUT/inv
    set qa_type = inv
    if ( $?REPCONFIG_INV3 ) then
	setenv REPCONFIG $REPCONFIG_INV3
    endif
endif
if ( $qa_type == inv2a ) then
    set outdir = $REPOUT/inv
    set qa_type = inv
    if ( $?REPCONFIG_INV2_A ) then
	setenv REPCONFIG $REPCONFIG_INV2_A
    endif
endif
if ( $qa_type == inv2b ) then
    set outdir = $REPOUT/inv
    set qa_type = inv
    if ( $?REPCONFIG_INV2_B ) then
	setenv REPCONFIG $REPCONFIG_INV2_B
    endif
endif
if ( $qa_type == inv2c ) then
    set outdir = $REPOUT/inv
    set qa_type = inv
    if ( $?REPCONFIG_INV2_C ) then
	setenv REPCONFIG $REPCONFIG_INV2_C
    endif
endif
if ( $qa_type == inv3a ) then
    set outdir = $REPOUT/inv
    set qa_type = inv
    if ( $?REPCONFIG_INV3_A ) then
	setenv REPCONFIG $REPCONFIG_INV3_A
    endif
endif
if ( $qa_type == inv3b ) then
    set outdir = $REPOUT/inv
    set qa_type = inv
    if ( $?REPCONFIG_INV3_B ) then
	setenv REPCONFIG $REPCONFIG_INV3_B
    endif
endif
if ( $qa_type == inv3c ) then
    set outdir = $REPOUT/inv
    set qa_type = inv
    if ( $?REPCONFIG_INV3_C ) then
	setenv REPCONFIG $REPCONFIG_INV3_C
    endif
endif
if ( $qa_type == invgrid ) then
    set outdir = $REPOUT/inv
    if ( $?REPCONFIG_GRID ) then
	setenv REPCONFIG $REPCONFIG_GRID
    endif
endif
if ( $qa_type == temporal ) then
    if ( $?REPCONFIG_TEMP ) then
	setenv REPCONFIG $REPCONFIG_TEMP
    endif
endif
if ( $qa_type == aermod ) then
    set outdir = $REPOUT/aermod
    if ( $?REPCONFIG_AERMOD ) then
        setenv REPCONFIG $REPCONFIG_AERMOD
    endif
endif
if ( $?REPCONFIG == 0 ) then
    echo "ERROR: REPCONFIG is not set"
    $EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: REPCONFIG is not set" -p $emf_period -t "e" ##Log w/ EMF server
    exit ( 1 )
endif

if ( $?outdir ) then
   mkdir -p $outdir
   chmod ug+rwx $outdir
endif

# Set up input file for Smkreport depending on settings

# Reports used for all sectors and source categories

if ( $qa_type != custom ) then
   setenv REPORT1   $outdir/rep_${label}_${qa_type}_state.txt
   setenv REPORT2   $outdir/rep_${label}_${qa_type}_state_scc.txt
   setenv REPORT3   $outdir/rep_${label}_${qa_type}_state_scctier2.txt
   setenv REPORT4   $outdir/rep_${label}_${qa_type}_county.txt
   setenv REPORT5   $outdir/rep_${label}_${qa_type}_county_scc_NAareas.txt
   setenv REPORT6   $outdir/rep_${label}_${qa_type}_scc_pm25prof.txt
   setenv REPORT7   $outdir/rep_${label}_${qa_type}_state_grid_${GRID}.txt
   setenv REPORT8   $outdir/rep_${label}_${qa_type}_county_scc.txt
   setenv REPORT8A  $outdir/rep_${label}_${qa_type}_county_scc7.txt # added by C. Allen (CSC), 07feb2012, for 2007 platform
   setenv REPORT9   $outdir/rep_${label}_${qa_type}_state_sic.txt
   setenv REPORT10  $outdir/rep_${label}_${qa_type}_facility_scc.txt
   setenv REPORT14  $outdir/rep_${label}_${qa_type}_scc_srgid_${GRID}.txt
   setenv REPORT14B $outdir/rep_${label}_${qa_type}_state_scc_srgid_${GRID}.txt
   setenv REPORT14N $outdir/rep_${label}_${qa_type}_national_srgid_${GRID}.txt
   setenv REPORT16  $outdir/rep_${label}_${qa_type}_state_scc_pm25prof.txt
   setenv REPORT17  $outdir/rep_${label}_${qa_type}_cell_${GRID}.txt
   setenv REPORT18  $outdir/rep_${label}_${qa_type}_county_moncode.txt
   setenv REPORT19  $outdir/rep_${label}_${qa_type}_state_naics.txt
   setenv REPORT20  $outdir/rep_${label}_${qa_type}_facility_oris.txt
   setenv REPORT21  $outdir/rep_${label}_${qa_type}_state_mact.txt
   setenv REPORT22  $outdir/rep_${label}_${qa_type}_state_scc_vocprof.txt
   setenv REPORT22N $outdir/rep_${label}_${qa_type}_state_scc_nonhapvocprof.txt
   setenv REPORT22A $outdir/rep_${label}_${qa_type}_state_scc_exhvocprof.txt
   setenv REPORT22B $outdir/rep_${label}_${qa_type}_state_scc_evpvocprof.txt
   setenv REPORT22C $outdir/rep_${label}_${qa_type}_state_scc_rflvocprof.txt
   setenv REPORT23  $outdir/rep_${label}_${qa_type}_cell_county_${GRID}.txt
   # REPORT24 added by C. Allen (CSC), 18may2011, for SMOKE-MOVES. County/SCC7 grid cell reports 
   # are needed in order to generate 12EUS1 and 12WUS1 reports based on 12MERGEUS1 emissions.
   # SCC7 used instead of SCC10, and without descrptions, to reduce size of reports.
   setenv REPORT24  $outdir/rep_${label}_${qa_type}_cell_county_scc7_${GRID}.txt 
   setenv REPORT25  $outdir/rep_${label}_${qa_type}_${GRID}.txt
   setenv REPORT26  $outdir/rep_src_${label}_${qa_type}_${GRID}.txt
   # These were added for the needs of the emissions modeling workgroup (7 Feb 2019)
   setenv REPORT31S  $outdir/rep_${label}_${qa_type}_state_srgid.txt
   setenv REPORT31C  $outdir/rep_${label}_${qa_type}_county_srgid.txt
   setenv REPORT32S  $outdir/rep_${label}_${qa_type}_state_pm25prof.txt
   setenv REPORT32C  $outdir/rep_${label}_${qa_type}_county_pm25prof.txt
   setenv REPORT32XS  $outdir/rep_${label}_${qa_type}_state_exhpm25prof.txt
   setenv REPORT32XC  $outdir/rep_${label}_${qa_type}_county_exhpm25prof.txt
   setenv REPORT32BS  $outdir/rep_${label}_${qa_type}_state_brkpm25prof.txt
   setenv REPORT32BC  $outdir/rep_${label}_${qa_type}_county_brkpm25prof.txt
   setenv REPORT32TS  $outdir/rep_${label}_${qa_type}_state_tirpm25prof.txt
   setenv REPORT32TC  $outdir/rep_${label}_${qa_type}_county_tirpm25prof.txt
   setenv REPORT33S  $outdir/rep_${label}_${qa_type}_state_vocprof.txt
   setenv REPORT33C  $outdir/rep_${label}_${qa_type}_county_vocprof.txt
   setenv REPORT33XS  $outdir/rep_${label}_${qa_type}_state_exhvocprof.txt
   setenv REPORT33XC  $outdir/rep_${label}_${qa_type}_county_exhvocprof.txt
   setenv REPORT33ES  $outdir/rep_${label}_${qa_type}_state_evpvocprof.txt
   setenv REPORT33EC  $outdir/rep_${label}_${qa_type}_county_evpvocprof.txt
   setenv REPORT33RS  $outdir/rep_${label}_${qa_type}_state_rflvocprof.txt
   setenv REPORT33RC  $outdir/rep_${label}_${qa_type}_county_rflvocprof.txt
   setenv REPORT34S  $outdir/rep_${label}_${qa_type}_state_nonhapvocprof.txt
   setenv REPORT34C  $outdir/rep_${label}_${qa_type}_county_nonhapvocprof.txt
   setenv REPORT34XS  $outdir/rep_${label}_${qa_type}_state_exhnonhapvocprof.txt
   setenv REPORT34XC  $outdir/rep_${label}_${qa_type}_county_exhnonhapvocprof.txt
   setenv REPORT34ES  $outdir/rep_${label}_${qa_type}_state_evpnonhapvocprof.txt
   setenv REPORT34EC  $outdir/rep_${label}_${qa_type}_county_evpnonhapvocprof.txt
   setenv REPORT34RS  $outdir/rep_${label}_${qa_type}_state_rflnonhapvocprof.txt
   setenv REPORT34RC  $outdir/rep_${label}_${qa_type}_county_rflnonhapvocprof.txt
   setenv REPORT35S  $outdir/rep_${label}_${qa_type}_state_temporalprof.txt
   setenv REPORT35C  $outdir/rep_${label}_${qa_type}_county_temporalprof.txt
   setenv REPORT35SC  $outdir/rep_${label}_${qa_type}_state_scc_temporalprof.txt
   setenv REPORT35N  $outdir/rep_${label}_${qa_type}_national_temporalprof.txt
   setenv REPORT40   $outdir/rep_${label}_${qa_type}_national_scc.txt
   setenv REPORT41  $outdir/rep_${label}_${qa_type}_national_vocprof.txt
   setenv REPORT41N $outdir/rep_${label}_${qa_type}_national_nonhapvocprof.txt
   setenv REPORT41A $outdir/rep_${label}_${qa_type}_national_exhvocprof.txt
   setenv REPORT41B $outdir/rep_${label}_${qa_type}_national_evpvocprof.txt
   setenv REPORT41C $outdir/rep_${label}_${qa_type}_national_rflvocprof.txt
   setenv REPORT34XN  $outdir/rep_${label}_${qa_type}_national_exhnonhapvocprof.txt
   setenv REPORT34EN  $outdir/rep_${label}_${qa_type}_national_evpnonhapvocprof.txt
   setenv REPORT34RN  $outdir/rep_${label}_${qa_type}_national_rflnonhapvocprof.txt
   setenv REPORT32N  $outdir/rep_${label}_${qa_type}_national_pm25prof.txt
   setenv REPORT32XN  $outdir/rep_${label}_${qa_type}_national_exhpm25prof.txt
   setenv REPORT32BN  $outdir/rep_${label}_${qa_type}_national_brkpm25prof.txt
   setenv REPORT32TN  $outdir/rep_${label}_${qa_type}_national_tirpm25prof.txt

   if ( $SMK_SOURCE == P ) then

      setenv REPORT10 $outdir/rep_${label}_${qa_type}_facility_scc.txt
      setenv REPORT11 $outdir/rep_${label}_${qa_type}_facility_cell_${GRID}.txt
      setenv REPORT12 $outdir/rep_${label}_${qa_type}_stackparm.txt
      setenv REPORT13 $outdir/rep_${label}_${qa_type}_facility.txt

   endif

endif

# Abort if already had an error
if ( $exitstat == 1 ) then
   $EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: detected in qa_run" -p $emf_period -t "e"  ## Log to EMF Server
   exit( $exitstat )
endif
  
#
### Smkreport processing for area, mobile, or point sources
#
set debugexestat = 0
set exestat = 0
setenv TMPLOG   $OUTLOG/smkreport_${label}_${QA_TYPE}.log   # Use original here for inv2
if ( $?RUN_SMKREPORT ) then
   if ( $RUN_SMKREPORT == Y ) then

      if ( -e $TMPLOG ) then
#	 $EMF_CLIENT -k $EMF_JOBKEY -x $movelog -p $emf_period ## log w/ EMF server
	 source $movelog
      endif

      if ( $exitstat == 0 ) then         # Run program
         setenv LOGFILE $TMPLOG
         if ( $debugmode == Y ) then
            if (-e $debug_exe/smkreport) then
  	       $debug_exe/smkreport
               if ( $status != 0 ) then
                  echo ERROR detected in smkreport
                  $EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: running smkreport debug version" -p $emf_period -t "e" ## Log to EMF Server
                  exit( 1 )
               endif
	    else
	       echo 'SCRIPT ERROR: smkreport program does not exist in:'
  	       echo '              '$debug_exe
               set exitstat = 1
	    endif
         else
            if ( -e $SMK_BIN/smkreport ) then

               set startdt = `date +%m/%d/%Y,%T`
#	       $EMF_CLIENT -k $EMF_JOBKEY -x $SMK_BIN/smkreport -p $emf_period ## log w/ EMF server
               time $SMK_BIN/smkreport
               if ( $timelog_yn == Y ) then
#		  $EMF_CLIENT -k $EMF_JOBKEY -x $timetracker -p $emf_period ## log w/ EMF server
                  $timetracker N $TIMELOG $startdt smkreport $runperiod
                  if ( $status != 0 ) then
                      echo "ERROR: Problem calling timetracker from qa_run script"
                      $EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: Problem calling timetracker from qa_run script" -x $timetracker  -t "e" -p $emf_period ## log w/ EMF server
                      exit ( 1 )
                  endif  
               endif 
	       
           endif
	
#	       $EMF_CLIENT -k $EMF_JOBKEY -x $checklog -m "Checking log $LOGFILE" -p $emf_period ## log w/ EMF server
               $checklog
               if ( $status != 0 ) then
                 echo ERROR detected in Smkreport
		 $EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: detected in Smkreport"  -p $emf_period -t "e" ## Log to EMF Server
                 exit( 1 )
               endif
            else
               set exestat = 1 
            endif
         endif
      endif

      if ( -e $SCRIPTS/fort.99 ) then
         mv $LOGFILE $LOGFILE.tmp
         cat $LOGFILE.tmp $SCRIPTS/fort.99 > $LOGFILE
         /bin/rm -rf $LOGFILE.tmp
         /bin/rm -rf $SCRIPTS/fort.99
      endif

      if ( $exestat == 1 ) then
	 echo 'SCRIPT ERROR: smkreport program does not exist in:'
	 echo '              '$SMK_BIN
         set exitstat = 1
      endif

      if ( $debugexestat == 1 ) then
	 echo 'SCRIPT ERROR: smkreport.debug program does not exist in:'
	 echo '              '$QA_SRC
         set exitstat = 1
      endif
      if ( $exitstat == 1 ) then
	 $EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: detected in Smkreport" -p $emf_period  -t "e" ## Log to EMF Server
      endif

   endif
endif

#
## Ending of script with exit status
#
exit( $exitstat )

