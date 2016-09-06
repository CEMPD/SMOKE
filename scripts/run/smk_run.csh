#!/bin/csh -f

# Version @(#)$Id$
# Path    $Source$
# Date    $Date$

# This script runs the SMOKE processors.  
#
# Time independent programs only need to be run once.
# Time dependent programs need to be processed once
# for each day needed for the air quality simulation.  
#
# Script created by : M. Houyoux and J. Vukovich, MCNC
#                     Environmental Modeling Center 
# Last edited : February 2003
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
   set debug_exe = idb
endif

### Ensure new controller variables are set
if ( $?RUN_PART1 ) then
   if ( $RUN_PART1 == Y || $RUN_PART1 == y ) then
      setenv RUN_PART1 Y
      echo 'Running part 1...'
   endif
else
   setenv RUN_PART1 N
endif
if ( $?RUN_PART2 ) then
   if ( $RUN_PART2 == Y || $RUN_PART2 == y ) then
      setenv RUN_PART2 Y
      echo "Running part 2, for $ESDATE ..."
   endif
else
   setenv RUN_PART2 N 
endif
if ( $?RUN_PART3 ) then
   if ( $RUN_PART3 == Y || $RUN_PART3 == y ) then
      setenv RUN_PART3 Y
      echo 'Running part 3 ...'
   endif
else
   setenv RUN_PART3 N 
endif
if ( $?RUN_PART4 ) then
   if ( $RUN_PART4 == Y || $RUN_PART4 == y ) then
      setenv RUN_PART4 Y
      echo "Running part 4, for $ESDATE..."
   endif
else
   setenv RUN_PART4 N 
endif

#
### Raw Inventory processing
#
set debugexestat = 0
set exestat = 0
setenv TMPLOG   $OUTLOG/smkinven.$SRCABBR.$INVEN.log
if ( $?RUN_SMKINVEN ) then
   if ( $RUN_SMKINVEN == 'Y' && $RUN_PART1 == Y ) then

      if ( -e $TMPLOG ) then
	 source $SCRIPTS/run/movelog.csh
      endif

      ##  Create output directories, if needed
      source $SCRIPTS/run/make_invdir.csh
      set exitstat = $status

      if ( $exitstat == 0 ) then         # Run program
         setenv LOGFILE $TMPLOG

         if ( $debugmode == Y ) then
            if ( -e $SMK_BIN/smkinven ) then
               $debug_exe $SMK_BIN/smkinven
            else
                set debugexestat = 1
            endif
         else
            if ( -e $SMK_BIN/smkinven ) then
               time $SMK_BIN/smkinven
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
	 echo 'SCRIPT ERROR: smkinven program does not exist in:'
	 echo '              '$SMK_BIN
         set exitstat = 1
      endif

   endif
endif

#
# BELD3 gridded landuse import for use in BEIS3
#

set debugexestat = 0
set exestat = 0
setenv TMPLOG   $OUTLOG/normbeis3.$SRCABBR.$INVEN.$GRID.log
if ( $?RUN_NORMBEIS3 ) then

   if ( $RUN_NORMBEIS3 == 'Y' && $RUN_PART1 == Y ) then

      set season = n
      # Use summer-specific processing if season-switch option in use
      if ( $?BIOSW_YN ) then
         if ( $BIOSW_YN == Y ) then
            set season = y
            setenv TMPLOG $OUTLOG/normbeis3.$SRCABBR.sumr.$INVEN.$GRID.log
         endif  
      endif  

      if ( -e $TMPLOG ) then
	 source $SCRIPTS/run/movelog.csh
      endif

      if ( $exitstat == 0 ) then         # Run program for standard or summer
         setenv LOGFILE $TMPLOG

         if ( $debugmode == Y ) then
            if ( -e $SMK_BIN/normbeis3 ) then
               $debug_exe $SMK_BIN/normbeis3
            else
                set debugexestat = 1
            endif
         else
            if ( -e $SMK_BIN/normbeis3 ) then
               time $SMK_BIN/normbeis3
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
      echo 'SCRIPT ERROR: normbeis3 program does not exist in:'
      echo '              '$SMK_BIN
      set exitstat = 1
   endif

endif

#
### Speciation Matrix generation
#
set debugexestat = 0
set exestat = 0
setenv TMPLOG   $OUTLOG/spcmat.$SRCABBR.$INVEN.$SPC.log
if ( $?RUN_SPCMAT ) then
   if ( $RUN_SPCMAT == 'Y' && $RUN_PART1 == Y ) then

      if ( -e $TMPLOG ) then
	 source $SCRIPTS/run/movelog.csh
      endif

      if ( $exitstat == 0 ) then         # Run program
         setenv LOGFILE $TMPLOG

         if ( $debugmode == Y ) then
            if ( -e $SMK_BIN/spcmat ) then
               $debug_exe $SMK_BIN/spcmat
            else
                set debugexestat = 1
            endif
         else
            if ( -e $SMK_BIN/spcmat ) then
               time $SMK_BIN/spcmat
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

   endif
endif

#
### Gridding Matrix generation
#
set debugexestat = 0
set exestat = 0
setenv TMPLOG   $OUTLOG/grdmat.$SRCABBR.$INVEN.$GRID.log
if ( $?RUN_GRDMAT ) then
   if ( $RUN_GRDMAT == 'Y' && $RUN_PART1 == Y ) then

      if ( -e $TMPLOG ) then
	 source $SCRIPTS/run/movelog.csh
      endif

      if ( $exitstat == 0 ) then         # Run program
         setenv LOGFILE $TMPLOG

         if ( $debugmode == Y ) then
            if ( -e $SMK_BIN/grdmat ) then
               $debug_exe $SMK_BIN/grdmat
            else
                set debugexestat = 1
            endif
         else
            if ( -e $SMK_BIN/grdmat ) then
               time $SMK_BIN/grdmat
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

   endif
endif

#
### Meteorology scan for bioseason file
# 
set debugexestat = 0
set exestat = 0
setenv TMPLOG   $OUTLOG/metscan.$SRCABBR.$INVEN.$ESDATE.$GRID.log
if ( $?RUN_METSCAN ) then
   if ( $RUN_METSCAN == 'Y' && $RUN_PART1 == Y ) then

      if ( -e $TMPLOG ) then
	 source $SCRIPTS/run/movelog.csh
      endif

      if ( $exitstat == 0 ) then         # Run program
        setenv LOGFILE $TMPLOG
        if ( $debugmode == Y ) then
            if ( -e $SMK_BIN/metscan ) then
               $debug_exe $SMK_BIN/metscan
            else
                set debugexestat = 1
            endif
         else
            if ( -e $SMK_BIN/metscan ) then
               time $SMK_BIN/metscan
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

   endif
endif

#
### Temporal Allocation - anthropogenic sources
# 
set debugexestat = 0
set exestat = 0
if ( $SMK_SOURCE == 'M' ) then
   setenv TMPLOG   $OUTLOG/temporal.$SRCABBR.$INVEN.$ESDATE.$GRID.log
else
   setenv TMPLOG   $OUTLOG/temporal.$SRCABBR.$INVEN.$ESDATE.log
endif
if ( $?RUN_TEMPORAL ) then
   if ( $RUN_TEMPORAL == 'Y' && $RUN_PART2 == Y ) then

      if ( -e $TMPLOG ) then
	 source $SCRIPTS/run/movelog.csh
      endif

      if ( $exitstat == 0 ) then         # Run program
         setenv LOGFILE $TMPLOG
        if ( $debugmode == Y ) then
            if ( -e $SMK_BIN/temporal ) then
               $debug_exe $SMK_BIN/temporal
            else
                set debugexestat = 1
            endif
         else
            if ( -e $SMK_BIN/temporal ) then
               time $SMK_BIN/temporal
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

   endif
endif

#
### Elevated source selection
#
set debugexestat = 0
set exestat = 0
setenv TMPLOG   $OUTLOG/elevpoint.$SRCABBR.$INVEN.log
if ( $?RUN_ELEVPOINT ) then
   if ( $RUN_ELEVPOINT == Y && $RUN_PART3 == Y ) then      

      if ( -e $TMPLOG ) then
	 source $SCRIPTS/run/movelog.csh
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
            if ( -e $SMK_BIN/elevpoint ) then
               $debug_exe $SMK_BIN/elevpoint
            else
                set debugexestat = 1
            endif
         else
            if ( -e $SMK_BIN/elevpoint ) then
               time $SMK_BIN/elevpoint
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

   endif
endif

#
### Layer fractions creation
#
set debugexestat = 0
set exestat = 0
setenv TMPLOG   $OUTLOG/laypoint.$SRCABBR.$INVEN.$ESDATE.$GRID.log
if ( $?RUN_LAYPOINT ) then
   if ( $RUN_LAYPOINT == 'Y' && $RUN_PART4 == Y ) then

      if ( $?SMK_PING_METHOD ) then
         if ( $SMK_PING_METHOD == 2 ) then
            setenv SMK_PING_YN Y
         endif
      endif

      if ( -e $TMPLOG ) then
	 source $SCRIPTS/run/movelog.csh
      endif

      if ( $exitstat == 0 ) then         # Run program
         setenv LOGFILE $TMPLOG
        if ( $debugmode == Y ) then
            if ( -e $SMK_BIN/laypoint ) then
               $debug_exe $SMK_BIN/laypoint
            else
                set debugexestat = 1
            endif
         else
            if ( -e $SMK_BIN/laypoint ) then
               time $SMK_BIN/laypoint
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

   endif
endif

#
# Temporal allocation for biogenics for BEIS3
#

set debugexestat = 0
set exestat = 0
setenv TMPLOG   $OUTLOG/tmpbeis3.$SRCABBR.$INVEN.$ESDATE.$GRID.log
if ( $?RUN_TMPBEIS3 ) then
   if ( $RUN_TMPBEIS3 == 'Y' && $RUN_PART2 == Y ) then

      if ( -e $TMPLOG ) then
	 source $SCRIPTS/run/movelog.csh
      endif

      if ( $exitstat == 0 ) then         # Run program
        setenv LOGFILE $TMPLOG
        if ( $debugmode == Y ) then
            if ( -e $SMK_BIN/tmpbeis3 ) then
               $debug_exe $SMK_BIN/tmpbeis3
            else
                set debugexestat = 1
            endif
         else
            if ( -e $SMK_BIN/tmpbeis3 ) then
               time $SMK_BIN/tmpbeis3
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
	 echo 'SCRIPT ERROR: tmpbeis3 program does not exist in:'
	 echo '              '$SMK_BIN
         set exitstat = 1
      endif

   endif
endif

#

#
### Merging
#
set debugexestat = 0
set exestat = 0
setenv TMPLOG   $OUTLOG/smkmerge.$SRCABBR.$INVEN.$ESDATE.$GRID.log
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
         source $SCRIPTS/run/movelog.csh
      endif

      if ( $exitstat == 0 ) then         # Run program
         setenv LOGFILE $TMPLOG
         if ( $debugmode == Y ) then
            if ( -e $SMK_BIN/smkmerge ) then
               $debug_exe $SMK_BIN/smkmerge
            else
                set debugexestat = 1
            endif
         else
            if ( -e $SMK_BIN/smkmerge ) then
               time $SMK_BIN/smkmerge
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

   endif
endif

#
### Merging for MOVES output
#
set debugexestat = 0
set exestat = 0
setenv TMPLOG   $OUTLOG/movesmrg.$SRCABBR.$INVEN.$ESDATE.$GRID.log
if ( $?RUN_MOVESMRG ) then
   if ( $RUN_MOVESMRG == 'Y' && $RUN_PART4 == Y ) then

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
         source $SCRIPTS/run/movelog.csh
      endif

      if ( $exitstat == 0 ) then         # Run program
         setenv LOGFILE $TMPLOG
         if ( $debugmode == Y ) then
            if ( -e $SMK_BIN/movesmrg ) then
               $debug_exe $SMK_BIN/movesmrg
            else
                set debugexestat = 1
            endif
         else
            if ( -e $SMK_BIN/movesmrg ) then
               time $SMK_BIN/movesmrg
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
	 echo 'SCRIPT ERROR: movesmrg program does not exist in:'
	 echo '              '$SMK_BIN
         set exitstat = 1
      endif

   endif
endif

# Merging for multiple sector output files

set debugexestat = 0
set exestat = 0
setenv TMPLOG   $OUTLOG/mrggrid.$SRCABBR.$INVEN.$ESDATE.$GRID.log
if ( $?RUN_MRGGRID ) then
   if ( $RUN_MRGGRID == 'Y' && $RUN_PART4 == Y ) then

      if ( -e $TMPLOG ) then
	 source $SCRIPTS/run/movelog.csh
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
            if ( -e $SMK_BIN/mrggrid ) then
               $debug_exe $SMK_BIN/mrggrid
            else
                set debugexestat = 1
            endif
         else
            if ( -e $SMK_BIN/mrggrid ) then
               time $SMK_BIN/mrggrid
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

   endif
endif

#
### Formats conversion
#
set debugexestat = 0
set exestat = 0
setenv TMPLOG   $OUTLOG/smk2emis.$SRCABBR.$INVEN.$ESDATE.$GRID.log
if ( $?RUN_SMK2EMIS ) then
   if ( $RUN_SMK2EMIS == 'Y' && $RUN_PART4 == Y ) then

      if ( -e $TMPLOG ) then
	 source $SCRIPTS/run/movelog.csh
      endif

      if ( $exitstat == 0 ) then         # Run program
         setenv LOGFILE $TMPLOG
         if ( $debugmode == Y ) then
            if ( -e $SMK_BIN/smk2emis ) then
               $debug_exe $SMK_BIN/smk2emis
            else
                set debugexestat = 1
            endif
         else
            if ( -e $SMK_BIN/smk2emis ) then
               time $SMK_BIN/smk2emis
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

   endif
endif

#
## Ending of script with exit status
#
exit( $exitstat )

