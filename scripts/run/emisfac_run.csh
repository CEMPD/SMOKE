#!/bin/csh -f

# Version @(#)$Id$
# Path    $Source$
# Date    $Date$

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
endif

#
### Emission factor creation
#
set debugexestat = 0
set exestat = 0
setenv TMPLOG   $OUTLOG/emisfac.$SRCABBR.$GROUP_TYPE.$GRID.log
if ( $?RUN_EMISFAC ) then
   if ( $RUN_EMISFAC == 'Y' ) then

      if ( -e $TMPLOG ) then
	 source $SCRIPTS/run/movelog.csh
      endif

      # Create M6LIST file...
      # First figure out if there are files in scenario directory
      setenv M6LIST $MBDAT/m6list.$MSCEN.$EF_YEAR.txt
#      set m6_input_dir = $MBDAT/m6_$EF_YEAR
#      if ( $?CNTLCASE ) then
#         setenv M6LIST $MBDAT/m6list.$MSCEN.${EF_YEAR}_$CNTLCASE.txt
#         set m6_input_dir = $MBDAT/mb_${EF_YEAR}_$CNTLCASE
#      endif
#      set ef_cnt = `ls -1 $m6_input_dir/*.in | wc -l`
#      if ( $ef_cnt > 0 ) then
#         if ( -e $M6LIST ) /bin/rm -rf $M6LIST
#         ls $m6_input_dir/*.in > $M6LIST
#      else
#         echo "SCRIPT ERROR: No MOBILE6 *.in scenario files found in the directory:"
#         echo "              $m6_input_dir"
#         echo "              Please put these files in the correct place and"
#         echo "              try again."
#         exit ( 1 )
#      endif

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

   endif
endif

exit( $exitstat )
