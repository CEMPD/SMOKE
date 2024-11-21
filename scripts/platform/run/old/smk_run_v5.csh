#!/bin/csh -f

# Version @(#)$Id: smk_run.csh,v 1.1.2.15 2003/09/22 18:59:21 cseppan Exp $
# Path    $Source: /afs/isis/depts/cep/emc/apps/archive/smoke/smoke/scripts/run/Attic/smk_run.csh,v $
# Date    $Date: 2003/09/22 18:59:21 $

# This script runs the SMOKE processors.  
#
# Time independent programs only need to be run once.
# Time dependent programs need to be processed once
# for each day needed for the air quality simulation.  
#
# Script created by : M. Houyoux and J. Vukovich, MCNC
#                     Environmental Modeling Center 
# Last edited : November 2006
#
# Modified to work w/ EMF: A Zubrow, UNC - IE August, 2007
#
#*********************************************************************

# Create directory for output program logs

setenv OUTLOG $LOGS
if ( ! -e $OUTLOG ) then
   mkdir -p $OUTLOG
   chmod ug+w $OUTLOG
endif

# Initialize exit status
set exitstat = 0

# Make sure that the timelog is set.  Whether the output file is
#    set or not determines whether the timetracking file is created.
# Also, time tracking is not written if using DEBUGMODE = Y
set timelog_yn = N
if ( $?TIMELOG && $DEBUGMODE != Y ) then
   set timelog_yn = Y
else
   echo "SCRIPT WARNING: Variable TIMELOG is not set, so time tracking file"
   echo "                will not be created."
endif

## Make sure that the EMF time period is set
if ( $?EMF_PERIOD == 0 ) then
    set emf_period = ""
else
    set emf_period = $EMF_PERIOD
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
   if ( $debugmode == Y ) then
      echo 'NOTE: DEBUG_EXE setting assumed to be "dbx" in smk_run.csh'
   endif
endif

## If EMF_CLIENT is not defined, assume that the user has not setup
#    the assigns file to use the new EMF-ready helper scripts
if ( ! $?EMF_CLIENT ) then
   setenv EMF_CLIENT false
   echo 'NOTE: EMF_CLIENT setting assumed to be "false" in smk_run.csh'
endif

## helper scripts
set movelog = $SCRIPTS/run/movelog.csh
set make_invdir = $SCRIPTS/run/make_invdir.csh
set timetracker = $SCRIPTS/run/timetracker_v2.csh 
set checklog = $SCRIPTS/run/checklogfile.csh

### Set version of BEIS3 executables to use.
set nb = normbeis3
set tb = tmpbeis3

### Ensure new controller variables are set
if ( $?RUN_PART1 ) then
   if ( $RUN_PART1 == Y || $RUN_PART1 == y ) then
      setenv RUN_PART1 Y
      echo 'Running part 1...'
      set runperiod = inv
   endif
else
   setenv RUN_PART1 N
endif
if ( $?RUN_PART2 ) then
   if ( $RUN_PART2 == Y || $RUN_PART2 == y ) then
      setenv RUN_PART2 Y
      echo "Running part 2, for $ESDATE ..."
      set runperiod = $ESDATE
   endif
else
   setenv RUN_PART2 N 
endif
if ( $?RUN_PART3 ) then
   if ( $RUN_PART3 == Y || $RUN_PART3 == y ) then
      setenv RUN_PART3 Y
      echo 'Running part 3 ...'
      set runperiod = $ESDATE
   endif
else
   setenv RUN_PART3 N 
endif
if ( $?RUN_PART4 ) then
   if ( $RUN_PART4 == Y || $RUN_PART4 == y ) then
      setenv RUN_PART4 Y
      echo "Running part 4, for $ESDATE..."
      set runperiod = $ESDATE
   endif
else
   setenv RUN_PART4 N 
endif

#
### Scan CEM Data
#
set debugexestat = 0
set exestat = 0
setenv TMPLOG   $OUTLOG/cemscan_${SRCABBR}_$CASE.log
if ( $?RUN_CEMSCAN ) then
   if ( $RUN_CEMSCAN == 'Y' && $RUN_PART0 == Y ) then

      if ( -e $TMPLOG ) then
	 source $movelog
      endif

      ##  Create output directories, if needed
      source $make_invdir
      set exitstat = $status

      if ( $exitstat == 0 ) then         # Run program
         setenv LOGFILE $TMPLOG

         if ( $debugmode == Y ) then
            if ( -e $IV_SRC/cemscan.debug ) then
               $debug_exe $IV_SRC/cemscan.debug
            else
                set debugexestat = 1
            endif
         else
            if ( -e $SMK_BIN/cemscan ) then

               set startdt = `date +%m/%d/%Y,%T`
#	       $EMF_CLIENT -k $EMF_JOBKEY -x $SMK_BIN/cemscan -m "Running Cemscan -- $emf_period" -p $emf_period ## Log to EMF server
               time $SMK_BIN/cemscan
               if ( $timelog_yn == Y ) then
                  $timetracker N $TIMELOG $startdt cemscan $runperiod
                  if ( $status != 0 ) then
                      echo "ERROR: Problem calling timetracker from smk_run script"
                      $EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: Problem calling timetracker from smk_run script" -x $timetracker  -t "e" -p $emf_period ## log w/ EMF server
                      exit ( 1 )
                  endif  
               endif 

#	       $EMF_CLIENT -k $EMF_JOBKEY -x $checklog -m "Checking log $LOGFILE" -p $emf_period ## log w/ EMF server
               $checklog
               if ( $status != 0 ) then
                 echo ERROR detected in Cemscan
		 $EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: detected in Cemscan" -p $emf_period -t "e" ## Log to EMF Server
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
	 echo 'SCRIPT ERROR: cemscan program does not exist in:'
	 echo '              '$SMK_BIN
         set exitstat = 1
      endif

      if ( $debugexestat == 1 ) then
	 echo 'SCRIPT ERROR: cemscan.debug program does not exist in:'
	 echo '              '$UT_SRC
         set exitstat = 1
      endif
      if ( $exitstat == 1 ) then
	 $EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: running Cemscan" -p $emf_period -t "e"  ## Log to EMF Server
	 exit ( 1 )
      endif

   endif
endif

#
### Raw Inventory processing
#
set debugexestat = 0
set exestat = 0
setenv TMPLOG   $OUTLOG/smkinven_${SRCABBR}_$CASE.log
if ( $?RUN_SMKINVEN ) then
   if ( $RUN_SMKINVEN == 'Y' && $RUN_PART1 == Y ) then

      if ( -e $TMPLOG ) then
	 source $movelog
      endif

      ##  Create output directories, if needed
      source $make_invdir
      set exitstat = $status

      if ( $exitstat == 0 ) then         # Run program
         setenv LOGFILE $TMPLOG

         if ( $debugmode == Y ) then
            if ( -e $IV_SRC/smkinven.debug ) then
               $debug_exe $IV_SRC/smkinven.debug
            else
                set debugexestat = 1
            endif
         else
            if ( -e $SMK_BIN/smkinven ) then

               set startdt = `date +%m/%d/%Y,%T`
#	       $EMF_CLIENT -k $EMF_JOBKEY -x $SMK_BIN/smkinven -m "Running Smkinven -- $emf_period" -p $emf_period ## Log to EMF server
               time $SMK_BIN/smkinven
               if ( $timelog_yn == Y ) then
                  $timetracker N $TIMELOG $startdt smkinven $runperiod
                  if ( $status != 0 ) then
                      echo "ERROR: Problem calling timetracker from smk_run script"
                      $EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: Problem calling timetracker from smk_run script" -x $timetracker  -t "e" -p $emf_period ## log w/ EMF server
                      exit ( 1 )
                  endif  
               endif 

#	       $EMF_CLIENT -k $EMF_JOBKEY -x $checklog -m "Checking log $LOGFILE" -p $emf_period ## log w/ EMF server
               $checklog
               if ( $status != 0 ) then
                 echo "ERROR: detected in Smkinven"
		 $EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: detected in Smkinven" -p $emf_period -t "e" ## Log to EMF Server
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

      if ( -e $SMK_TMPDIR/import_tmp ) then
         /bin/rm -rf $SMK_TMPDIR/import_tmp
      endif

      if ( $exestat == 1 ) then
	 echo 'SCRIPT ERROR: smkinven program does not exist in:'
	 echo '              '$SMK_BIN
         set exitstat = 1
      endif

      if ( $debugexestat == 1 ) then
	 echo 'SCRIPT ERROR: smkinven.debug program does not exist in:'
	 echo '              '$IV_SRC
         set exitstat = 1
      endif
      if ( $exitstat == 1 ) then
	 $EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: running Smkinven" -p $emf_period -t "e" ## Log to EMF Server
	 exit ( 1 )
      endif

   endif
endif

#
### County or gridded landuse import
#
set debugexestat = 0
set exestat = 0
setenv TMPLOG   $OUTLOG/rawbio_${SRCABBR}_${CASE}_$GRID.log
if ( $?RUN_RAWBIO ) then

   if ( $RUN_RAWBIO == 'Y' && $RUN_PART1 == Y ) then

      # Use summer emission factors, if they are set for first rawbio run
      if ( $?S_BFAC ) then
          setenv BFAC $S_BFAC

      else

         if ( $?BFAC ) then
            echo 'SCRIPT WARNING: default biogenic emission factors, BFAC file:'
            echo '               '$BFAC
            echo '                are used directly based on user setting.'

         else
            echo 'SCRIPT ERROR: neither BFAC nor S_BFAC environment variables are set'
            set exitstat = 1

         endif

      endif

      set season = n
      # Use summer-specific processing if season-switch option in use
      if ( $?BIOSW_YN ) then
         if ( $BIOSW_YN == Y ) then
            set season = y
            setenv TMPLOG $OUTLOG/rawbio_${SRCABBR}_summer_${CASE}_$GRID.log
         endif  
      endif  

      if ( -e $TMPLOG ) then
	 source $movelog
      endif

      if ( $exitstat == 0 ) then         # Run program for standard or summer
         setenv LOGFILE $TMPLOG

         if ( $debugmode == Y ) then
            if ( -e $BG_SRC/rawbio.debug ) then
               $debug_exe $BG_SRC/rawbio.debug
            else
                set debugexestat = 1
            endif
         else
            if ( -e $SMK_BIN/rawbio ) then

               set startdt = `date +%m/%d/%Y,%T`
#	       $EMF_CLIENT -k $EMF_JOBKEY -x $SMK_BIN/rawbio -m "Running Rawbio -- $emf_period" -p $emf_period ## Log to EMF server
               time $SMK_BIN/rawbio
               if ( $timelog_yn == Y ) then
                  $timetracker N $TIMELOG $startdt rawbio $runperiod
                  if ( $status != 0 ) then
                      echo "ERROR: Problem calling timetracker from smk_run script"
                      $EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: Problem calling timetracker from smk_run script" -x $timetracker  -t "e" -p $emf_period ## log w/ EMF server
                      exit ( 1 )
                  endif  
               endif

#	       $EMF_CLIENT -k $EMF_JOBKEY -x $checklog -m "Checking log $LOGFILE" -p $emf_period ## log w/ EMF server
               $checklog
               if ( $status != 0 ) then
                 echo ERROR detected in Rawbio
		 $EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: detected in Rawbio" -p $emf_period -t "e" ## Log to EMF Server
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

      # Now do winter-specific processing, if season-switch option in use
      if ( $season == y ) then

         if ( $?W_BFAC && $?BGRDW ) then
            # Set input file 
            setenv BFAC $W_BFAC
            setenv BGRD $BGRDW
            setenv TMPLOG $OUTLOG/rawbio${SRCABBR}_winter_${CASE}_$GRID.log

            if ( -e $TMPLOG ) then
	       source $movelog
            endif

            if ( $exitstat == 0 ) then         # Run program for winter
               setenv LOGFILE $TMPLOG

               if ( $debugmode == Y ) then
                  if ( -e $BG_SRC/rawbio.debug ) then
                     $debug_exe $BG_SRC/rawbio.debug
                  else
                      set debugexestat = 1
                  endif
               else
                  if ( -e $SMK_BIN/rawbio ) then

                     set startdt = `date +%m/%d/%Y,%T`
#		     $EMF_CLIENT -k $EMF_JOBKEY -x $SMK_BIN/rawbio -m "Running Rawbio -- $emf_period" -p $emf_period ## Log to EMF server
                     time $SMK_BIN/rawbio
                     if ( $timelog_yn == Y ) then
                        $timetracker N $TIMELOG $startdt rawbio $runperiod
                        if ( $status != 0 ) then
                            echo "ERROR: Problem calling timetracker from smk_run script"
                            $EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: Problem calling timetracker from smk_run script" -x $timetracker  -t "e" -p $emf_period ## log w/ EMF server
                            exit ( 1 )
                        endif  
                     endif

#		     $EMF_CLIENT -k $EMF_JOBKEY -x $checklog -m "Checking log $LOGFILE" -p $emf_period ## log w/ EMF server
        	     $checklog
        	     if ( $status != 0 ) then
                       echo ERROR detected in Rawbio
		       $EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: detected in Rawbio" -p $emf_period -t "e" ## Log to EMF Server
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

         else
      	    echo 'SCRIPT ERROR: BIOSW_YN (biogenics seasonality switch) set to'
            echo '              Y, but W_BFAC (winter biogenic factors input) and/or'
            echo '              BGRDW (winter normalized emissions) environment'
            echo '              variables are not defined.'
	    set exitstat = 1
         endif
    
      endif
   endif

   if ( $exestat == 1 ) then
      echo 'SCRIPT ERROR: rawbio program does not exist in:'
      echo '              '$SMK_BIN
      set exitstat = 1
   endif

   if ( $debugexestat == 1 ) then
      echo 'SCRIPT ERROR: rawbio.debug program does not exist in:'
      echo '              '$BG_SRC
      set exitstat = 1
   endif
   if ( $exitstat == 1 ) then
      $EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: running Rawbio" -p $emf_period -t "e" ## Log to EMF Server
      exit ( 1 )
   endif

endif


#
# BELD3 gridded landuse import for use in BEIS3
#

set debugexestat = 0
set exestat = 0
setenv TMPLOG   $OUTLOG/${nb}_${SRCABBR}_${CASE}_$GRID.log
if ( $?RUN_NORMBEIS3 ) then

   if ( $RUN_NORMBEIS3 == 'Y' && $RUN_PART1 == Y ) then

      set season = n
      # Use summer-specific processing if season-switch option in use
      if ( $?BIOSW_YN ) then
         if ( $BIOSW_YN == Y ) then
            set season = y
            setenv TMPLOG $OUTLOG/${nb}_${SRCABBR}_summer_${CASE}_$GRID.log
         endif  
      endif  

      if ( -e $TMPLOG ) then
	 source $movelog
      endif

      if ( $exitstat == 0 ) then         # Run program for standard or summer
         setenv LOGFILE $TMPLOG

         if ( $debugmode == Y ) then
            if ( -e $BG_SRC/$nb.debug ) then
               $debug_exe $BG_SRC/$nb.debug
            else
                set debugexestat = 1
            endif
         else
            if ( -e $SMK_BIN/$nb ) then

               set startdt = `date +%m/%d/%Y,%T`
#	       $EMF_CLIENT -k $EMF_JOBKEY -x $SMK_BIN/$nb -m "Running Normbeis -- $emf_period" -p $emf_period ## Log to EMF server
               time $SMK_BIN/$nb
               if ( $timelog_yn == Y ) then
                  $timetracker N $TIMELOG $startdt $nb $runperiod
                  if ( $status != 0 ) then
                      echo "ERROR: Problem calling timetracker from smk_run script"
                      $EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: Problem calling timetracker from smk_run script" -x $timetracker  -t "e" -p $emf_period ## log w/ EMF server
                      exit ( 1 )
                  endif  
               endif

#	       $EMF_CLIENT -k $EMF_JOBKEY -x $checklog -m "Checking log $LOGFILE" -p $emf_period ## log w/ EMF server
               $checklog
               if ( $status != 0 ) then
                 echo ERROR detected in $nb
		 $EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: detected in $nb" -p $emf_period -t "e" ## Log to EMF Server
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

   endif

   if ( $exestat == 1 ) then
      echo 'SCRIPT ERROR: '$nb' program does not exist in:'
      echo '              '$SMK_BIN
      set exitstat = 1
   endif

   if ( $debugexestat == 1 ) then
      echo 'SCRIPT ERROR: '$nb'.debug program does not exist in:'
      echo '              '$BG_SRC
      set exitstat = 1
   endif
   if ( $exitstat == 1 ) then
      $EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: running Normbeis" -p $emf_period -t "e" ## Log to EMF Server
      exit ( 1 )
   endif

endif

#
### Speciation Matrix generation
#
set debugexestat = 0
set exestat = 0
setenv TMPLOG   $OUTLOG/spcmat_${SRCABBR}_${CASE}_$SPC.log
if ( $?RUN_SPCMAT ) then
   if ( $RUN_SPCMAT == 'Y' && $RUN_PART1 == Y ) then

      if ( -e $TMPLOG ) then
	 source $movelog
      endif

      if ( $exitstat == 0 ) then         # Run program
         setenv LOGFILE $TMPLOG

         if ( $debugmode == Y ) then
            if ( -e $SP_SRC/spcmat.debug ) then
               $debug_exe $SP_SRC/spcmat.debug
            else
                set debugexestat = 1
            endif
         else
            if ( -e $SMK_BIN/spcmat ) then

               set startdt = `date +%m/%d/%Y,%T`
#	       $EMF_CLIENT -k $EMF_JOBKEY -x $SMK_BIN/spcmat -m "Running Spcmat -- $emf_period" -p $emf_period ## Log to EMF server
               time $SMK_BIN/spcmat
               if ( $timelog_yn == Y ) then
                  $timetracker N $TIMELOG $startdt spcmat $runperiod
                  if ( $status != 0 ) then
                      echo "ERROR: Problem calling timetracker from smk_run script"
                      $EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: Problem calling timetracker from smk_run script" -x $timetracker  -t "e" -p $emf_period ## log w/ EMF server
                      exit ( 1 )
                  endif  
               endif

#	       $EMF_CLIENT -k $EMF_JOBKEY -x $checklog -m "Checking log $LOGFILE" -p $emf_period ## log w/ EMF server
               $checklog
               if ( $status != 0 ) then
                 echo "ERROR: detected in Spcmat"
		 $EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: detected in Spcmat" -p $emf_period -t "e" ## Log to EMF Server
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
	 echo 'SCRIPT ERROR: spcmat program does not exist in:'
	 echo '              '$SMK_BIN
         set exitstat = 1
      endif

      if ( $debugexestat == 1 ) then
	 echo 'SCRIPT ERROR: spcmat.debug program does not exist in:'
	 echo '              '$SP_SRC
         set exitstat = 1
      endif
      if ( $exitstat == 1 ) then
	 $EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: running Spcmat" -p $emf_period -t "e" ## Log to EMF Server
	 exit ( 1 )
      endif

   endif
endif

#
### Gridding Matrix generation
#
set debugexestat = 0
set exestat = 0
setenv TMPLOG   $OUTLOG/grdmat_${SRCABBR}_${CASE}_$GRID.log
if ( $?RUN_GRDMAT ) then
   if ( $RUN_GRDMAT == 'Y' && $RUN_PART1 == Y ) then

      if ( -e $TMPLOG ) then
	 source $movelog
      endif

      if ( $exitstat == 0 ) then         # Run program
         setenv LOGFILE $TMPLOG

         if ( $debugmode == Y ) then
            if ( -e $GD_SRC/grdmat.debug ) then
               $debug_exe $GD_SRC/grdmat.debug
            else
                set debugexestat = 1
            endif
         else
            if ( -e $SMK_BIN/grdmat ) then

               set startdt = `date +%m/%d/%Y,%T`
#	       $EMF_CLIENT -k $EMF_JOBKEY -x $SMK_BIN/grdmat -m "Running Grdmat -- $emf_period" -p $emf_period ## Log to EMF server
               time $SMK_BIN/grdmat
               if ( $timelog_yn == Y ) then
                  $timetracker N $TIMELOG $startdt grdmat $runperiod
                  if ( $status != 0 ) then
                      echo "ERROR: Problem calling timetracker from smk_run script"
                      $EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: Problem calling timetracker from smk_run script" -x $timetracker  -t "e" -p $emf_period ## log w/ EMF server
                      exit ( 1 )
                  endif  
               endif

#	       $EMF_CLIENT -k $EMF_JOBKEY -x $checklog -m "Checking log $LOGFILE" -p $emf_period ## log w/ EMF server
               $checklog
               if ( $status != 0 ) then
                 echo "ERROR: detected in Grdmat"
		 $EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: detected in Grdmat" -p $emf_period -t "e" ## Log to EMF Server
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
	 echo 'SCRIPT ERROR: grdmat program does not exist in:'
	 echo '              '$SMK_BIN
         set exitstat = 1
      endif

      if ( $debugexestat == 1 ) then
	 echo 'SCRIPT ERROR: grdmat.debug program does not exist in:'
	 echo '              '$GD_SRC
         set exitstat = 1
      endif
      if ( $exitstat == 1 ) then
	 $EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: running Grdmat" -p $emf_period -t "e" ## Log to EMF Server
	 exit ( 1 )
      endif
   endif
endif

#
### Mobile setup for MOBILE6 runs
#
set debugexestat = 0
set exestat = 0
setenv TMPLOG   $OUTLOG/mbsetup_${SRCABBR}_${CASE}.log
if ( $?RUN_MBSETUP ) then
   if ( $RUN_MBSETUP == 'Y' && $RUN_PART1 == Y ) then

      if ( -e $TMPLOG ) then
         source $movelog
      endif

      if ( $exitstat == 0 ) then         # Run program
         setenv LOGFILE $TMPLOG
         if ( $debugmode == Y ) then
            if ( -e $MB_SRC/mbsetup.debug ) then
               $debug_exe $MB_SRC/mbsetup.debug
            else
                set debugexestat = 1
            endif
         else
            if ( -e $SMK_BIN/mbsetup ) then

               set startdt = `date +%m/%d/%Y,%T`
#	       $EMF_CLIENT -k $EMF_JOBKEY -x $SMK_BIN/mbsetup -m "Running Mbsetup -- $emf_period" -p $emf_period ## Log to EMF server
               time $SMK_BIN/mbsetup
               if ( $timelog_yn == Y ) then
                  $timetracker N $TIMELOG $startdt mbsetup $runperiod
                  if ( $status != 0 ) then
                      echo "ERROR: Problem calling timetracker from smk_run script"
                      $EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: Problem calling timetracker from smk_run script" -x $timetracker  -t "e" -p $emf_period ## log w/ EMF server
                      exit ( 1 )
                  endif  
               endif
 
#	       $EMF_CLIENT -k $EMF_JOBKEY -x $checklog -m "Checking log $LOGFILE" -p $emf_period ## log w/ EMF server
               $checklog
               if ( $status != 0 ) then
                 echo ERROR detected in Mbsetup
		 $EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: detected in Mbsetup" -p $emf_period -t "e" ## Log to EMF Server
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
	 echo 'SCRIPT ERROR: mbsetup program does not exist in:'
	 echo '              '$SMK_BIN
         set exitstat = 1
      endif

      if ( $debugexestat == 1 ) then
	 echo 'SCRIPT ERROR: mbsetup.debug program does not exist in:'
	 echo '              '$MB_SRC
         set exitstat = 1
      endif
      if ( $exitstat == 1 ) then
	 $EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: running Mbsetup" -p $emf_period -t "e" ## Log to EMF Server
	 exit ( 1 )
      endif

   endif
endif

#
### Pre-mobile-source processing
#
set debugexestat = 0
set exestat = 0
setenv TMPLOG   $OUTLOG/premobl_${SRCABBR}_${CASE}_$GRID.log
if ( $?RUN_PREMOBL ) then
   if ( $RUN_PREMOBL == 'Y' && $RUN_PART1 == Y ) then

      if ( -e $TMPLOG ) then
	 source $movelog
      endif

      if ( -e $METLIST ) /bin/rm -rf $METLIST

      ls $METDAT/METCRO2D* > $METLIST

      if ( $exitstat == 0 ) then         # Run program
         setenv LOGFILE $TMPLOG

         if ( $debugmode == Y ) then
            if ( -e $MB_SRC/premobl.debug ) then
               $debug_exe $MB_SRC/premobl.debug
            else
                set debugexestat = 1
            endif
         else
            if ( -e $SMK_BIN/premobl ) then

               set startdt = `date +%m/%d/%Y,%T`
#	       $EMF_CLIENT -k $EMF_JOBKEY -x $SMK_BIN/premobl -m "Running Premobl -- $emf_period" -p $emf_period ## Log to EMF server
               time $SMK_BIN/premobl
               if ( $timelog_yn == Y ) then
                  $timetracker N $TIMELOG $startdt premobl $runperiod
                  if ( $status != 0 ) then
                      echo "ERROR: Problem calling timetracker from smk_run script"
                      $EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: Problem calling timetracker from smk_run script" -x $timetracker  -t "e" -p $emf_period ## log w/ EMF server
                      exit ( 1 )
                  endif  
               endif

#	       $EMF_CLIENT -k $EMF_JOBKEY -x $checklog -m "Checking log $LOGFILE" -p $emf_period ## log w/ EMF server
               $checklog
               if ( $status != 0 ) then
                 echo ERROR detected in Premobl
		 $EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: detected in Premobl" -p $emf_period -t "e" ## Log to EMF Server
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
	 echo 'SCRIPT ERROR: premobl program does not exist in:'
	 echo '              '$SMK_BIN
         set exitstat = 1
      endif

      if ( $debugexestat == 1 ) then
	 echo 'SCRIPT ERROR: premobl.debug program does not exist in:'
	 echo '              '$MB_SRC
         set exitstat = 1
      endif
      if ( $exitstat == 1 ) then
         $EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: running Premobl" -p $emf_period -t "e" ## Log to EMF Server
	 exit ( 1 )
      endif

   endif
endif

#
### Meteorology scan for bioseason file
# 
set debugexestat = 0
set exestat = 0
setenv TMPLOG   $OUTLOG/metscan_${SRCABBR}_${CASE}_${ESDATE}_$GRID.log
if ( $?RUN_METSCAN ) then
   if ( $RUN_METSCAN == 'Y' && $RUN_PART1 == Y ) then

      if ( -e $TMPLOG ) then
	 source $movelog
      endif

      if ( $exitstat == 0 ) then         # Run program
        setenv LOGFILE $TMPLOG
        if ( $debugmode == Y ) then
            if ( -e $UT_SRC/metscan.debug ) then
               $debug_exe $UT_SRC/metscan.debug
            else
                set debugexestat = 1
            endif
         else
            if ( -e $SMK_BIN/metscan ) then

               set startdt = `date +%m/%d/%Y,%T`
#	       $EMF_CLIENT -k $EMF_JOBKEY -x $SMK_BIN/metscan -m "Running Metscan -- $emf_period" -p $emf_period ## Log to EMF server
               time $SMK_BIN/metscan
               if ( $timelog_yn == Y ) then
                  $timetracker N $TIMELOG $startdt metscan $runperiod
                  if ( $status != 0 ) then
                      echo "ERROR: Problem calling timetracker from smk_run script"
                      $EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: Problem calling timetracker from smk_run script" -x $timetracker  -t "e" -p $emf_period ## log w/ EMF server
                      exit ( 1 )
                  endif  
               endif

#	       $EMF_CLIENT -k $EMF_JOBKEY -x $checklog -m "Checking log $LOGFILE" -p $emf_period ## log w/ EMF server
               $checklog
               if ( $status != 0 ) then
                 echo ERROR detected in Metscan
		 $EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: detected in Metscan" -p $emf_period -t "e" ## Log to EMF Server
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
	 echo 'SCRIPT ERROR: metscan program does not exist in:'
	 echo '              '$SMK_BIN
         set exitstat = 1
      endif

      if ( $debugexestat == 1 ) then
	 echo 'SCRIPT ERROR: metscan.debug program does not exist in:'
	 echo '              '$UT_SRC
         set exitstat = 1
      endif
      if ( $exitstat == 1 ) then
         $EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: running Metscan" -p $emf_period -t "e" ## Log to EMF Server
	 exit ( 1 )
      endif

   endif
endif

#
### Temporal Allocation - anthropogenic sources
# 
set debugexestat = 0
set exestat = 0
if ( $SMK_SOURCE == 'M' ) then
   setenv TMPLOG   $OUTLOG/temporal_${SRCABBR}_${CASE}_${ESDATE}_$GRID.log
else
   setenv TMPLOG   $OUTLOG/temporal_${SRCABBR}_${CASE}_${ESDATE}.log
endif
if ( $?RUN_TEMPORAL ) then
   if ( $RUN_TEMPORAL == 'Y' && $RUN_PART2 == Y ) then

      if ( -e $TMPLOG ) then
	 source $movelog
      endif

      ## For mobile sources, create MEFLIST file if there are any files
      #    to put into the list
      if ( $SMK_SOURCE == M ) then
         setenv MEFLIST $SMK_EMISPATH/meflist.txt
         set ef_cnt = `ls $SMK_EMISPATH/emisfacs*ncf | wc -l`
         if ( $ef_cnt > 0 ) then

            if ( -e $SMK_EMISPATH/meflist.txt ) /bin/rm -rf $SMK_EMISPATH/meflist.txt

            set extension = `ls $SMK_EMISPATH/emisfacs*ncf | cut -d. -f4`
            if ( $extension[1] == 'ncf' ) then
               ls $SMK_EMISPATH/emisfacs*ncf > $SMK_EMISPATH/meflist.txt
            else 
               ls $SMK_EMISPATH/emisfacs*1.ncf | cut -d. -f1,2,3,5 > $SMK_EMISPATH/meflist.txt
            endif
         endif
      endif

      if ( $exitstat == 0 ) then         # Run program
         setenv LOGFILE $TMPLOG
        if ( $debugmode == Y ) then
            if ( -e $TM_SRC/temporal.debug ) then
               $debug_exe $TM_SRC/temporal.debug
            else
                set debugexestat = 1
            endif
         else
            if ( -e $SMK_BIN/temporal ) then

               set startdt = `date +%m/%d/%Y,%T`
#	       $EMF_CLIENT -k $EMF_JOBKEY -x $SMK_BIN/temporal -m "Running Temporal -- $emf_period" -p $emf_period ## Log to EMF server
               time $SMK_BIN/temporal
               if ( $timelog_yn == Y ) then
                  $timetracker N $TIMELOG $startdt temporal $runperiod
                  if ( $status != 0 ) then
                      echo "ERROR: Problem calling timetracker from smk_run script"
                      $EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: Problem calling timetracker from smk_run script" -x $timetracker  -t "e" -p $emf_period ## log w/ EMF server
                      exit ( 1 )
                  endif  
               endif

#	       $EMF_CLIENT -k $EMF_JOBKEY -x $checklog -m "Checking log $LOGFILE" -p $emf_period ## log w/ EMF server
               $checklog
               if ( $status != 0 ) then
                 echo "ERROR: detected in Temporal"
		 $EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: detected in Temporal" -p $emf_period -t "e" ## Log to EMF Server
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
	 echo 'SCRIPT ERROR: temporal program does not exist in:'
	 echo '              '$SMK_BIN
         set exitstat = 1
      endif

      if ( $debugexestat == 1 ) then
	 echo 'SCRIPT ERROR: temporal.debug program does not exist in:'
	 echo '              '$TM_SRC
         set exitstat = 1
      endif
      if ( $exitstat == 1 ) then
	 $EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: running Temporal" -p $emf_period -t "e" ## Log to EMF Server
	 exit ( 1 )
      endif
   endif
endif

#
### Elevated source selection
#
set debugexestat = 0
set exestat = 0
setenv TMPLOG   $OUTLOG/elevpoint_${SRCABBR}_${CASE}.log
if ( $?RUN_ELEVPOINT ) then
   if ( $RUN_ELEVPOINT == Y && $RUN_PART3 == Y ) then      

      if ( -e $TMPLOG ) then
	 source $movelog
      endif

      # Create PTMPLIST file, in case it is needed.
      if ( -e $PTMP ) then
         setenv PTMPLIST $INVDIR/other/ptmplist.txt
         if ( -e $PTMPLIST ) /bin/rm -rf $PTMPLIST

         ls $SMKDAT/run_$PSCEN/*/ptmp*$PSCEN*ncf > $PTMPLIST
      endif

      if ( $exitstat == 0 ) then         # Run program
         setenv LOGFILE $TMPLOG
        if ( $debugmode == Y ) then
            if ( -e $PT_SRC/elevpoint.debug ) then
               $debug_exe $PT_SRC/elevpoint.debug
            else
                set debugexestat = 1
            endif
         else
            if ( -e $SMK_BIN/elevpoint ) then

               set startdt = `date +%m/%d/%Y,%T`
#	       $EMF_CLIENT -k $EMF_JOBKEY -x $SMK_BIN/elevpoint -m "Running Elevpoint -- $emf_period" -p $emf_period ## Log to EMF server
               time $SMK_BIN/elevpoint
               if ( $timelog_yn == Y ) then
                  $timetracker N $TIMELOG $startdt elevpoint $runperiod
                  if ( $status != 0 ) then
                      echo "ERROR: Problem calling timetracker from smk_run script"
                      $EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: Problem calling timetracker from smk_run script" -x $timetracker  -t "e" -p $emf_period ## log w/ EMF server
                      exit ( 1 )
                  endif  
               endif

#	       $EMF_CLIENT -k $EMF_JOBKEY -x $checklog -m "Checking log $LOGFILE" -p $emf_period ## log w/ EMF server
               $checklog
               if ( $status != 0 ) then
                 echo ERROR detected in Elevpoint
		 $EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: detected in Elevpoint" -p $emf_period -t "e" ## Log to EMF Server
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
	 echo 'SCRIPT ERROR: elevpoint program does not exist in:'
	 echo '              '$SMK_BIN
         set exitstat = 1
      endif

      if ( $debugexestat == 1 ) then
	 echo 'SCRIPT ERROR: elevpoint.debug program does not exist in:'
	 echo '              '$PT_SRC
         set exitstat = 1
      endif
      if ( $exitstat == 1 ) then
         $EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: running Elevpoint" -p $emf_period -t "e" ## Log to EMF Server
	 exit ( 1 )
      endif

   endif
endif

#
### Layer fractions creation
#
set debugexestat = 0
set exestat = 0
setenv TMPLOG   $OUTLOG/laypoint_${SRCABBR}_${CASE}_${ESDATE}_$GRID.log
if ( $?RUN_LAYPOINT ) then
   if ( $RUN_LAYPOINT == 'Y' && $RUN_PART4 == Y ) then

      if ( $?SMK_PING_METHOD ) then
         if ( $SMK_PING_METHOD == 2 ) then
            setenv SMK_PING_YN Y
         endif
      endif

      if ( -e $TMPLOG ) then
	 source $movelog
      endif

      if ( $exitstat == 0 ) then         # Run program
         setenv LOGFILE $TMPLOG
        if ( $debugmode == Y ) then
            if ( -e $PT_SRC/laypoint.debug ) then
               $debug_exe $PT_SRC/laypoint.debug
            else
                set debugexestat = 1
            endif
         else
            if ( -e $SMK_BIN/laypoint ) then

               set startdt = `date +%m/%d/%Y,%T`
#	       $EMF_CLIENT -k $EMF_JOBKEY -x $SMK_BIN/laypoint -m "Running Laypoint -- $emf_period" -p $emf_period ## Log to EMF server
               time $SMK_BIN/laypoint
               if ( $timelog_yn == Y ) then
                  $timetracker N $TIMELOG $startdt laypoint $runperiod
                  if ( $status != 0 ) then
                      echo "ERROR: Problem calling timetracker from smk_run script"
                      $EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: Problem calling timetracker from smk_run script" -x $timetracker  -t "e" -p $emf_period ## log w/ EMF server
                      exit ( 1 )
                  endif  
               endif

#	       $EMF_CLIENT -k $EMF_JOBKEY -x $checklog -m "Checking log $LOGFILE" -p $emf_period ## log w/ EMF server
               $checklog
               if ( $status != 0 ) then
                 echo ERROR detected in Laypoint
		 $EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: detected in Laypoint" -p $emf_period -t "e" ## Log to EMF Server
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
	 echo 'SCRIPT ERROR: laypoint program does not exist in:'
	 echo '              '$SMK_BIN
         set exitstat = 1
      endif

      if ( $debugexestat == 1 ) then
	 echo 'SCRIPT ERROR: laypoint.debug program does not exist in:'
	 echo '              '$PT_SRC
         set exitstat = 1
      endif
      if ( $exitstat == 1 ) then
         $EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: running Laypoint" -p $emf_period -t "e" ## Log to EMF Server
	 exit ( 1 )
      endif

   endif
endif

#
### Temporal Allocation - biogenic sources
# 
set debugexestat = 0
set exestat = 0
setenv TMPLOG   $OUTLOG/tmpbio_${SRCABBR}_${CASE}_${ESDATE}_$GRID.log
if ( $?RUN_TMPBIO ) then
   if ( $RUN_TMPBIO == 'Y' && $RUN_PART2 == Y ) then

      if ( -e $TMPLOG ) then
	 source $movelog
      endif

      if ( -e $METLIST ) /bin/rm -rf $METLIST
      if ( -e $RADLIST ) /bin/rm -rf $RADLIST

      ls $MET_FILE1 > $METLIST
      ls $MET_FILE2 > $RADLIST

      if ( $exitstat == 0 ) then         # Run program
        setenv LOGFILE $TMPLOG
        if ( $debugmode == Y ) then
            if ( -e $BG_SRC/tmpbio.debug ) then
               $debug_exe $BG_SRC/tmpbio.debug
            else
                set debugexestat = 1
            endif
         else
            if ( -e $SMK_BIN/tmpbio ) then

               set startdt = `date +%m/%d/%Y,%T`
#	       $EMF_CLIENT -k $EMF_JOBKEY -x $SMK_BIN/tmpbio -m "Running Tmpbio -- $emf_period" -p $emf_period ## Log to EMF server
               time $SMK_BIN/tmpbio
               if ( $timelog_yn == Y ) then
                  $timetracker N $TIMELOG $startdt tmpbio $runperiod
                  if ( $status != 0 ) then
                      echo "ERROR: Problem calling timetracker from smk_run script"
                      $EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: Problem calling timetracker from smk_run script" -x $timetracker  -t "e" -p $emf_period ## log w/ EMF server
                      exit ( 1 )
                  endif  
               endif

#	       $EMF_CLIENT -k $EMF_JOBKEY -x $checklog -m "Checking log $LOGFILE" -p $emf_period ## log w/ EMF server
               $checklog
               if ( $status != 0 ) then
                 echo ERROR detected in Tmpbio
		 $EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: detected in Tmpbio" -p $emf_period -t "e" ## Log to EMF Server 
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
	 echo 'SCRIPT ERROR: tmpbio program does not exist in:'
	 echo '              '$SMK_BIN
         set exitstat = 1
      endif

      if ( $debugexestat == 1 ) then
	 echo 'SCRIPT ERROR: tmpbio.debug program does not exist in:'
	 echo '              '$BG_SRC
         set exitstat = 1
      endif
      if ( $exitstat == 1 ) then
         $EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: running Tmpbio" -p $emf_period -t "e" ## Log to EMF Server
	 exit ( 1 )
      endif

   endif
endif

#
# Temporal allocation for biogenics for BEIS3
#

set debugexestat = 0
set exestat = 0
setenv TMPLOG   $OUTLOG/${tb}_${SRCABBR}_${CASE}_${ESDATE}_$GRID.log
if ( $?RUN_TMPBEIS3 ) then
   if ( $RUN_TMPBEIS3 == 'Y' && $RUN_PART2 == Y ) then

      if ( -e $TMPLOG ) then
	 source $movelog
      endif

      if ( $exitstat == 0 ) then         # Run program
        setenv LOGFILE $TMPLOG
        if ( $debugmode == Y ) then
            if ( -e $BG_SRC/$tb.debug ) then
               $debug_exe $BG_SRC/$tb.debug
            else
                set debugexestat = 1
            endif
         else
            if ( -e $SMK_BIN/$tb ) then

               set startdt = `date +%m/%d/%Y,%T`
#	       $EMF_CLIENT -k $EMF_JOBKEY -x $SMK_BIN/$tb -m "Running Beis -- $emf_period" -p $emf_period ## Log to EMF server
               time $SMK_BIN/$tb
               if ( $timelog_yn == Y ) then
                  $timetracker N $TIMELOG $startdt $tb $runperiod
                  if ( $status != 0 ) then
                      echo "ERROR: Problem calling timetracker from smk_run script"
                      $EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: Problem calling timetracker from smk_run script" -x $timetracker  -t "e" -p $emf_period ## log w/ EMF server
                      exit ( 1 )
                  endif  
               endif

#	       $EMF_CLIENT -k $EMF_JOBKEY -x $checklog -m "Checking log $LOGFILE" -p $emf_period ## log w/ EMF server
               $checklog
               if ( $status != 0 ) then
                 echo ERROR detected in $tb
		 $EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: detected in $tb" -p $emf_period -t "e" ## Log to EMF Server
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
	 echo "SCRIPT ERROR: "$tb" program does not exist in:"
	 echo "              "$SMK_BIN
         set exitstat = 1
      endif

      if ( $debugexestat == 1 ) then
	 echo "SCRIPT ERROR: "$tb".debug program does not exist in:"
	 echo "              "$BG_SRC
         set exitstat = 1
      endif
      if ( $exitstat == 1 ) then
         $EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: running $tb" -p $emf_period -t "e" ## Log to EMF Server
	 exit ( 1 )
      endif

   endif
endif

#

#
### Merging
#
set debugexestat = 0
set exestat = 0
setenv TMPLOG   $OUTLOG/smkmerge_${SRCABBR}_${CASE}_${ESDATE}_${GRID}_$SPC.log
if ( $?RUN_SMKMERGE ) then
   if ( $RUN_SMKMERGE == 'Y' && $RUN_PART4 == Y ) then

      # Set mole/mass-based speciation matrices.
      set unit = tons
      if( $?MRG_GRDOUT_UNIT ) then 
         set unit = `echo $MRG_GRDOUT_UNIT | cut -c1-4`
      endif
      if ( $unit == mole ) then   # mole
         if ( $?ASMAT_L ) then
            setenv ASMAT $ASMAT_L
         endif
         if ( $?BGTS_L ) then
            setenv BGTS $BGTS_L
         endif
         if ( $?MSMAT_L ) then
            setenv MSMAT $MSMAT_L
         endif
         if ( $?PSMAT_L ) then
            setenv PSMAT $PSMAT_L
         endif

      else                       # mass
         if ( $?ASMAT_S ) then
            setenv ASMAT $ASMAT_S
         endif
         if ( $?BGTS_S ) then
            setenv BGTS $BGTS_S
         endif
         if ( $?MSMAT_S ) then
            setenv MSMAT $MSMAT_S
         endif
         if ( $?PSMAT_S ) then
            setenv PSMAT $PSMAT_S
         endif
      endif

      if ( -e $TMPLOG ) then
         source $movelog
      endif

      if ( $exitstat == 0 ) then         # Run program
         setenv LOGFILE $TMPLOG
         if ( $debugmode == Y ) then
            if ( -e $MG_SRC/smkmerge.debug ) then
               $debug_exe $MG_SRC/smkmerge.debug
            else
                set debugexestat = 1
            endif
         else
            if ( -e $SMK_BIN/smkmerge ) then

               set startdt = `date +%m/%d/%Y,%T`
#	       $EMF_CLIENT -k $EMF_JOBKEY -x $SMK_BIN/smkmerge -m "Running Smkmerge -- $emf_period" -p $emf_period ## Log to EMF server
               time $SMK_BIN/smkmerge
               if ( $timelog_yn == Y ) then
                  $timetracker N $TIMELOG $startdt smkmerge $runperiod
                  if ( $status != 0 ) then
                      echo "ERROR: Problem calling timetracker from smk_run script"
                      $EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: Problem calling timetracker from smk_run script" -x $timetracker  -t "e" -p $emf_period ## log w/ EMF server
                      exit ( 1 )
                  endif  
               endif

#	       $EMF_CLIENT -k $EMF_JOBKEY -x $checklog -m "Checking log $LOGFILE" -p $emf_period ## log w/ EMF server
               $checklog
               if ( $status != 0 ) then
                 echo "ERROR: detected in Smkmerge"
		 $EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: detected in Smkmerge" -p $emf_period -t "e" ## Log to EMF Server
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
	 echo 'SCRIPT ERROR: smkmerge program does not exist in:'
	 echo '              '$SMK_BIN
         set exitstat = 1
      endif

      if ( $debugexestat == 1 ) then
	 echo 'SCRIPT ERROR: smkmerge.debug program does not exist in:'
	 echo '              '$MG_SRC
         set exitstat = 1
      endif
      if ( $exitstat == 1 ) then
	 $EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: running Smkmerge" -p $emf_period -t "e" ## Log to EMF Server
	 exit ( 1 )
      endif
   endif
endif

set debugexestat = 0
set exestat = 0
setenv TMPLOG   $OUTLOG/mrggrid_${SRCABBR}_${CASE}_${ESDATE}_${GRID}_$SPC.log
if ( $?RUN_MRGGRID ) then
   if ( $RUN_MRGGRID == 'Y' && $RUN_PART4 == Y ) then

      if ( -e $TMPLOG ) then
	 source $movelog
      endif

      ## If Mrggrid file list defined, create list file
      if ( $?MRGFILES ) then

         ## Remove existing file list, if it is there
         if ( -e $FILELIST ) then
            /bin/rm -rf $FILELIST
         endif
         set mrg_cnt = 0
         foreach f ( $MRGFILES )
            @ mrg_cnt = $mrg_cnt + 1
            if ( $mrg_cnt == 1 ) then
               echo $f > $FILELIST
            else
               echo $f >> $FILELIST
            endif
         end

         if ( $mrg_cnt == 0 ) then
             echo "SCRIPT ERROR: MRGFILES defined, but no logical file names included."
             echo "              MRGFILES script variable = "${MRGFILES}.
             echo "              Please reset and rerun script"
             set exitstat = 1
         echo
             echo "SCRIPT NOTE: File FILELIST created with logical files:"
             echo "             "$MRGFILES
         endif
      endif

      if ( $exitstat == 0 ) then         # Run program
         setenv LOGFILE $TMPLOG
         if ( $debugmode == Y ) then
            if ( -e $MG_SRC/mrggrid.debug ) then
               $debug_exe $MG_SRC/mrggrid.debug
            else
                set debugexestat = 1
            endif
         else
            if ( -e $SMK_BIN/mrggrid ) then

               set startdt = `date +%m/%d/%Y,%T`
#	       $EMF_CLIENT -k $EMF_JOBKEY -x $SMK_BIN/mrggrid -m "Running Mrggrid -- $emf_period" -p $emf_period ## Log to EMF server
               time $SMK_BIN/mrggrid
               if ( $timelog_yn == Y ) then
                  $timetracker N $TIMELOG $startdt mrggrid $runperiod
                  if ( $status != 0 ) then
                      echo "ERROR: Problem calling timetracker from smk_run script"
                      $EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: Problem calling timetracker from smk_run script" -x $timetracker  -t "e" -p $emf_period ## log w/ EMF server
                      exit ( 1 )
                  endif  
               endif

#	       $EMF_CLIENT -k $EMF_JOBKEY -x $checklog -m "Checking log $LOGFILE" -p $emf_period ## log w/ EMF server
               $checklog
               if ( $status != 0 ) then
                 echo "ERROR: detected in Mrggrid"
		 $EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: detected in Mrggrid" -p $emf_period -t "e" ## Log to EMF Server
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
	 echo 'SCRIPT ERROR: mrggrid program does not exist in:'
	 echo '              '$SMK_BIN
         set exitstat = 1
      endif

      if ( $debugexestat == 1 ) then
	 echo 'SCRIPT ERROR: mrggrid.debug program does not exist in:'
	 echo '              '$MG_SRC
         set exitstat = 1
      endif
      if ( $exitstat == 1 ) then
	 $EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: running Mrggrid" -p $emf_period -t "e" ## Log to EMF Server
	 exit ( 1 )
      endif
   endif
endif

#
### Formats conversion
#
set debugexestat = 0
set exestat = 0
setenv TMPLOG   $OUTLOG/smk2emis_${SRCABBR}_${CASE}_${ESDATE}_$GRID.log
if ( $?RUN_SMK2EMIS ) then
   if ( $RUN_SMK2EMIS == 'Y' && $RUN_PART4 == Y ) then

      if ( -e $TMPLOG ) then
	 source $movelog
      endif

      if ( $exitstat == 0 ) then         # Run program
         setenv LOGFILE $TMPLOG
         if ( $debugmode == Y ) then
            if ( -e $UT_SRC/smk2emis.debug ) then
               $debug_exe $UT_SRC/smk2emis.debug
            else
                set debugexestat = 1
            endif
         else
            if ( -e $SMK_BIN/smk2emis ) then

               set startdt = `date +%m/%d/%Y,%T`
#	       $EMF_CLIENT -k $EMF_JOBKEY -x $SMK_BIN/smk2emis -m "Running Smk2emis -- $emf_period" -p $emf_period ## Log to EMF server
               time $SMK_BIN/smk2emis
               if ( $timelog_yn == Y ) then
                  $timetracker N $TIMELOG $startdt smk2emis $runperiod
                  if ( $status != 0 ) then
                      echo "ERROR: Problem calling timetracker from smk_run script"
                      $EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: Problem calling timetracker from smk_run script" -x $timetracker  -t "e" -p $emf_period ## log w/ EMF server
                      exit ( 1 )
                  endif  
               endif

#	       $EMF_CLIENT -k $EMF_JOBKEY -x $checklog -m "Checking log $LOGFILE" -p $emf_period ## log w/ EMF server
               $checklog
               if ( $status != 0 ) then
                 echo ERROR detected in Smk2emis
		 $EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: detected in Smk2emis" -p $emf_period -t "e" ## Log to EMF Server
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
	 echo 'SCRIPT ERROR: smk2emis program does not exist in:'
	 echo '              '$SMK_BIN
         set exitstat = 1
      endif

      if ( $debugexestat == 1 ) then
	 echo 'SCRIPT ERROR: smk2emis.debug program does not exist in:'
	 echo '              '$UT_SRC
         set exitstat = 1
      endif
      if ( $exitstat == 1 ) then
         $EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: running Smk2emis" -p $emf_period -t "e" ## Log to EMF Server
	 exit ( 1 )
      endif

   endif
endif

#
## Ending of script with exit status
#
exit( $exitstat )

