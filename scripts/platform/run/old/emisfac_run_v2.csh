#!/bin/csh -f

# Version @(#)$Id: emisfac_run.csh,v 1.1.2.7 2003/08/14 13:57:32 cseppan Exp $
# Path    $Source: /afs/isis/depts/cep/emc/apps/archive/smoke/smoke/scripts/run/Attic/emisfac_run.csh,v $
# Date    $Date: 2003/08/14 13:57:32 $

# This script runs the SMOKE Emisfac processor for MOBILE6 emission factors.  
#
# Script created by : M. Houyoux 
#                     Carolina Environmental Program
# Last edited : February, 2003
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
set movelog = $SCRIPTS/run/movelog.csh

#
### Emission factor creation
#
set debugexestat = 0
set exestat = 0
setenv TMPLOG   $OUTLOG/emisfac.$SRCABBR.$GROUP_TYPE.$GRID.log
if ( $?RUN_EMISFAC ) then
   if ( $RUN_EMISFAC == 'Y' ) then

      if ( -e $TMPLOG ) then
	 source $movelog
      endif

      # Create M6LIST file...
      # First figure out if there are files in scenario directory
      setenv M6LIST $MBDAT/m6list.$MSCEN.$EF_YEAR.txt
      set m6_input_dir = $MBDAT/m6_$EF_YEAR
      if ( $?CNTLCASE ) then
         setenv M6LIST $MBDAT/m6list.$MSCEN.${EF_YEAR}_$CNTLCASE.txt
         set m6_input_dir = $MBDAT/mb_${EF_YEAR}_$CNTLCASE
      endif
      set ef_cnt = `ls -1 $m6_input_dir/*.in | wc -l`
      if ( $ef_cnt > 0 ) then
         if ( -e $M6LIST ) /bin/rm -rf $M6LIST
         ls $m6_input_dir/*.in > $M6LIST
      else
         echo "SCRIPT ERROR: No MOBILE6 *.in scenario files found in the directory:"
         echo "              $m6_input_dir"
         echo "              Please put these files in the correct place and"
         echo "              try again."
         exit ( 1 )
      endif

      ## Set HOURLYT input file
      setenv HOURLYT `ls -1 $SMK_METPATH/*${GROUP_TYPE}* | head -1`
      if ( $HOURLYT == ' ' ) then
         echo 'SCRIPT ERROR: Could not find Premobl output files'
         echo "              for GROUP_TYPE $GROUP_TYPE in SMK_METPATH"
         echo "              directory:"
         echo "              $METPATH"
         set exitstat == 1
      endif

      if ( $exitstat == 0 ) then         # Run program
         setenv LOGFILE $TMPLOG
         if ( $debugmode == Y ) then
            if ( -e $MB_SRC/emisfac.debug ) then
               $debug_exe $MB_SRC/emisfac.debug
            else
                set debugexestat = 1
            endif
         else
            if ( -e $SMK_BIN/emisfac ) then
	       $EMF_CLIENT -k $EMF_JOBKEY -x $SMK_BIN/emisfac -m "Running Emisfac -- $emf_period" -p $emf_period ## Log to EMF server
               time $SMK_BIN/emisfac
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
	 echo 'SCRIPT ERROR: emisfac program does not exist in:'
	 echo '              '$SMK_BIN
         set exitstat = 1
      endif

      if ( $debugexestat == 1 ) then
	 echo 'SCRIPT ERROR: emisfac.debug program does not exist in:'
	 echo '              '$MB_SRC
         set exitstat = 1
      endif
      if ( $exitstat == 1 ) then
         $EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: script error in emisfac_run" -p $emf_period -t "e" ## Log to EMF Server
      endif

   endif
endif

exit( $exitstat )
