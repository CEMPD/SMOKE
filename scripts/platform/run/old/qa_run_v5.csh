#!/bin/csh -f
#
# $Version$
# $Path$
# $Date: 2003/03/14 00:30:52 $
#
# This script runs the SMOKE QA processors.  
#
# Script created by : M. Houyoux, North Carolina 
#                     Supercomputing Center
# Last edited : September, 2000
#
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

if ( $?DEBUG_EXE ) then
   set debug_exe = $DEBUG_EXE
else
   set debug_exe = dbx
   if ( $debugmode == Y ) then
      echo 'NOTE: DEBUG_EXE setting assumed to be "dbx" in' !$
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
set timetracker = $SCRIPTS/run/timetracker.csh
set checklog = $SCRIPTS/run/checklogfile.csh
set movelog = $SCRIPTS/run/movelog.csh

# Make sure that the timelog is set.  Whether the output file is
#    set or not determines whether the timetracking file is created.
# Also, time tracking is not written if using DEBUGMODE = Y
set timelog_yn = N
if ( $?TIMELOG && $DEBUGMODE != Y ) then
   set timelog_yn = Y
endif

# Make sure that debug mode and debug executable are set
if ( $?DEBUGMODE ) then
   set debugmode = $DEBUGMODE
else
   set debugmode = N
endif

if ( $?DEBUG_EXE ) then
   set debug_exe = $DEBUG_EXE
else
   set debug_exe = dbx
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
   if ( $RUN_PART2 == Y || $RUN_PART3 == Y || $RUN_PART4 == Y ) then
      set label = ${CASE}_${ESDATE}_${QA_LABEL}
   else
      set label = ${CASE}_$QA_LABEL
   endif

# If QA_LABEL is not set, set report "label" variable to defaults
else

   # For any part of processing that is day-specific, add date in label
   if ( $RUN_PART2 == Y || $RUN_PART3 == Y || $RUN_PART4 == Y ) then
      set label = ${CASE}_${ESDATE}
   else
      set label = ${CASE}
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

# Adjust RUN_ID to consider that it might have things appended to it
# Check if RUN_ID has an underscore included
#set leftrun_id = ( `echo $SECTOR |  cut -d"_" -f1` )
#if ( $SECTOR == $leftrun_id ) then
   set run_id = $SECTOR 
#else
#   set run_id = $leftrun_id
#endif

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
if ( $qa_type == inv2 ) then
    set outdir = $REPOUT/inv
    set qa_type = inv
    if ( $?REPCONFIG_INV2 ) then
	setenv REPCONFIG $REPCONFIG_INV2
    endif
endif
if ( $qa_type == temporal ) then
    if ( $?REPCONFIG_TEMP ) then
	setenv REPCONFIG $REPCONFIG_TEMP
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
   setenv REPORT1   $outdir/rep_${run_id}_${qa_type}_state_$label.txt
   setenv REPORT2   $outdir/rep_${run_id}_${qa_type}_state_scc_$label.txt
   setenv REPORT3   $outdir/rep_${run_id}_${qa_type}_state_scctier2_$label.txt
   setenv REPORT4   $outdir/rep_${run_id}_${qa_type}_county_$label.txt
   setenv REPORT5   $outdir/rep_${run_id}_${qa_type}_county_scc_NAareas_$label.txt
   setenv REPORT6   $outdir/rep_${run_id}_${qa_type}_scc_pm25prof_$label.txt
   setenv REPORT7   $outdir/rep_${run_id}_${qa_type}_state_grid_${GRID}_$label.txt
   setenv REPORT8   $outdir/rep_${run_id}_${qa_type}_county_scc_$label.txt
   setenv REPORT9   $outdir/rep_${run_id}_${qa_type}_state_sic_$label.txt
   setenv REPORT10  $outdir/rep_${run_id}_${qa_type}_plant_scc_$label.txt
   setenv REPORT14  $outdir/rep_${run_id}_${qa_type}_scc_srgid_$label.txt
   setenv REPORT14B $outdir/rep_${run_id}_${qa_type}_state_scc_srgid_$label.txt
   setenv REPORT15  $outdir/rep_${run_id}_${qa_type}_state_pm25prof_$label.txt
   setenv REPORT16  $outdir/rep_${run_id}_${qa_type}_state_scc_pm25prof_$label.txt
   setenv REPORT17  $outdir/rep_${run_id}_${qa_type}_cell_$label.txt
   setenv REPORT18  $outdir/rep_${run_id}_${qa_type}_county_moncode_$label.txt
   setenv REPORT19  $outdir/rep_${run_id}_${qa_type}_plant_oris_$label.txt
   setenv REPORT20  $outdir/rep_${run_id}_${qa_type}_state_naics_$label.txt
   setenv REPORT21  $outdir/rep_${run_id}_${qa_type}_state_mact_$label.txt
   setenv REPORT22  $outdir/rep_${run_id}_${qa_type}_state_scc_vocprof_$label.txt

   if ( $SMK_SOURCE == P ) then

      setenv REPORT10 $outdir/rep_${run_id}_${qa_type}_plant_scc_$label.txt
      setenv REPORT11 $outdir/rep_${run_id}_${qa_type}_plant_cell_${GRID}_$label.txt
      setenv REPORT12 $outdir/rep_${run_id}_${qa_type}_stackparm_$label.txt
      setenv REPORT13 $outdir/rep_${run_id}_${qa_type}_plant_$label.txt

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
setenv TMPLOG   $OUTLOG/smkreport_${run_id}_${qa_type}_$label.log
if ( $?RUN_SMKREPORT ) then
   if ( $RUN_SMKREPORT == Y ) then

      if ( -e $TMPLOG ) then
#	 $EMF_CLIENT -k $EMF_JOBKEY -x $movelog -p $emf_period ## log w/ EMF server
	 source $movelog
      endif

      if ( $exitstat == 0 ) then         # Run program
         setenv LOGFILE $TMPLOG
         if ( $debugmode == Y ) then
            if ( -e $QA_SRC/smkreport.debug ) then
               $debug_exe $QA_SRC/smkreport.debug
            else
                set debugexestat = 1
            endif
         else
            if ( -e $SMK_BIN/smkreport ) then

               set startdt = `date +%m/%d/%Y,%T`
#	       $EMF_CLIENT -k $EMF_JOBKEY -x $SMK_BIN/smkreport -p $emf_period ## log w/ EMF server
               time $SMK_BIN/smkreport
               if ( $timelog_yn == Y ) then
#		  $EMF_CLIENT -k $EMF_JOBKEY -x $timetracker -p $emf_period ## log w/ EMF server
                  $timetracker N $TIMELOG $startdt smkreport $ESDATE
                  if ( $status != 0 ) then
                      echo "ERROR: Problem calling timetracker from qa_run script"
                      $EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: Problem calling timetracker from qa_run script" -x $timetracker  -t "e" -p $emf_period ## log w/ EMF server
                      exit ( 1 )
                  endif  
               endif 
	       
               if ( $?BUD_REPORTS ) then 
                  if ( $BUD_REPORTS == Y ) then
	             echo "Copying Reports to bud"
		     setenv BUD_REPORT  $BUD_INSTALL/2002/smoke/reports/$platform/$CASE/$QA_TYPE

		     if ( $qa_type != custom ) then
	              scp $outdir/rep_${run_id}_${qa_type}_state_$label.txt bud.rtp.epa.gov:$BUD_REPORT >&! /dev/null
                      scp $outdir/rep_${run_id}_${qa_type}_state_scc_$label.txt bud.rtp.epa.gov:$BUD_REPORT >&! /dev/null
                      scp $outdir/rep_${run_id}_${qa_type}_state_scctier2_$label.txt bud.rtp.epa.gov:$BUD_REPORT >&! /dev/null
                      scp $outdir/rep_${run_id}_${qa_type}_county_$label.txt bud.rtp.epa.gov:$BUD_REPORT >&! /dev/null
                      scp $outdir/rep_${run_id}_${qa_type}_county_scc_NAareas_$label.txt bud.rtp.epa.gov:$BUD_REPORT >&! /dev/null
                      scp $outdir/rep_${run_id}_${qa_type}_scc_pm25prof_$label.txt bud.rtp.epa.gov:$BUD_REPORT >&! /dev/null
                      scp $outdir/rep_${run_id}_${qa_type}_state_grid_${GRID}_$label.txt bud.rtp.epa.gov:$BUD_REPORT >&! /dev/null
                      scp $outdir/rep_${run_id}_${qa_type}_county_scc_all_areas_$label.txt bud.rtp.epa.gov:$BUD_REPORT >&! /dev/null
                      scp $outdir/rep_${run_id}_${qa_type}_state_pm25prof_$label.txt bud.rtp.epa.gov:$BUD_REPORT >&! /dev/null
                      scp $outdir/rep_${run_id}_${qa_type}_state_scc_pm25prof_$label.txt bud.rtp.epa.gov:$BUD_REPORT >&! /dev/null
		            if ( $SMK_SOURCE == P ) then
                             scp $outdir/rep_${run_id}_${qa_type}_state_sic_$label.txt bud.rtp.epa.gov:$BUD_REPORT >&! /dev/null
                             scp $outdir/rep_${run_id}_${qa_type}_plant_scc_$label.txt bud.rtp.epa.gov:$BUD_REPORT >&! /dev/null
                             scp $outdir/rep_${run_id}_${qa_type}_plant_cell_${GRID}_$label.txt bud.rtp.epa.gov:$BUD_REPORT >&! /dev/null
                             scp $outdir/rep_${run_id}_${qa_type}_stackparm_$label.txt bud.rtp.epa.gov:$BUD_REPORT >&! /dev/null
                             scp $outdir/rep_${run_id}_${qa_type}_plant_$label.txt bud.rtp.epa.gov:$BUD_REPORT >&! /dev/null
			    else
                             scp $outdir/rep_${run_id}_${qa_type}_scc_srgid_$label.txt bud.rtp.epa.gov:$BUD_REPORT >&! /dev/null
                            endif
		    endif
		     if ( $status != 0 ) then
                        echo ERROR detected in copying to bud - run scp -v to debug
                     endif
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

