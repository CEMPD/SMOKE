#!/bin/csh -fx
#
# $Version$
# $Path$
# $Date$
#
# This script runs the SMOKE QA processors.  
#
# Script created by : M. Houyoux, North Carolina 
#                     Supercomputing Center
# Last edited : September, 2000
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
   exit( $exitstat )
endif

# Check if QA label is set, and initialize file
if ( $?QA_LABEL ) then

   set ilabl = $FYIOP.$QA_LABEL
   set slabl = $ESCEN.$QA_LABEL

else

   set ilabl = $INVOP
   set slabl = $ESCEN

endif

switch ( $QA_TYPE ) 
   case custom:
   case none:
      set qa_type = $QA_TYPE
      breaksw

   case part1:
   case PART1:
      if ( $?RUN_PART1 ) then
         if ( $RUN_PART1 == Y ) then
            set qa_type = inventory
         else
            set qa_type = none
         endif
      else
         set qa_type = none
      endif
      breaksw

   case part2:
   case PART2:
      if ( $?RUN_PART2 ) then
         if ( $RUN_PART2 == Y ) then
            set qa_type = temporal
         else
            set qa_type = none
         endif
      else
         set qa_type = none
      endif
      breaksw

   case part3:
   case PART3:
      if ( $?RUN_PART3 ) then
         if ( $RUN_PART3 == Y ) then
            set qa_type = elevpoint
         else
            set qa_type = none
         endif
      else
         set qa_type = none
      endif
      breaksw

   case part4:
   case PART4:
      if ( $?RUN_PART4 ) then
         if ( $RUN_PART4 == Y ) then
            set qa_type = laypoint
         else
            set qa_type = none
         endif
      else
         set qa_type = none
      endif

   case all:
      set qa_type = inventory
      if ( $?RUN_PART1 ) then
         if ( $RUN_PART1 == Y ) then
            set qa_type = inventory
            breaksw
         else
            set qa_type = none
         endif
      else
         set qa_type = none
         breaksw
      endif
      if ( $?RUN_PART2 ) then
         if ( $RUN_PART2 == Y ) then
            set qa_type = temporal
            breaksw
         else
            set qa_type = none
         endif
      endif
      if ( $?RUN_PART3 ) then
         if ( $RUN_PART3 == Y && $SMK_SOURCE == P ) then
            set qa_type = elevpoint
            breaksw
         else
            set qa_type = none
         endif
      endif
      if ( $?RUN_PART4 ) then
         if ( $RUN_PART4 == Y && $SMK_SOURCE == P ) then
            set qa_type = laypoint
            breaksw
         else
            set qa_type = none
         endif
      endif
      breaksw

endsw

# Set up input file for Smkreport depending on settings

switch ( $qa_type ) 

case inventory:

   switch ( $SMK_SOURCE )
   case A:
      setenv REPCONFIG $INVDIR/other/repconfig.ar.inv.txt
      setenv REPORT1 $REPSTAT/a.state.$ilabl.rpt
      setenv REPORT2 $REPSTAT/a.county.$ilabl.rpt
      setenv REPORT3 $REPSTAT/a.scc.$ilabl.rpt
      setenv REPORT4 $REPSTAT/a.state_scc.$ilabl.rpt
      setenv REPORT5 $REPSTAT/ag.state.$GRID.$ilabl.rpt
      setenv REPORT6 $REPSTAT/ag.scc.$GRID.$ilabl.rpt
      breaksw
      
   case M:
      setenv REPCONFIG $INVDIR/other/repconfig.mb.inv.txt
      setenv REPORT1 $REPSTAT/m.state.$ilabl.rpt
      setenv REPORT2 $REPSTAT/m.county.$ilabl.rpt
      setenv REPORT3 $REPSTAT/m.scc.$ilabl.rpt
      setenv REPORT4 $REPSTAT/m.state_scc.$ilabl.rpt
      setenv REPORT5 $REPSTAT/m.state_rclas.$ilabl.rpt
      setenv REPORT6 $REPSTAT/mg.state.$GRID.$ilabl.rpt
      setenv REPORT7 $REPSTAT/mg.scc.$GRID.$ilabl.rpt
      breaksw
      
   case P:
      setenv REPCONFIG $INVDIR/other/repconfig.pt.inv.txt
      setenv REPORT1 $REPSTAT/p.state.$ilabl.rpt
      setenv REPORT2 $REPSTAT/p.county.$ilabl.rpt
      setenv REPORT3 $REPSTAT/p.scc.$ilabl.rpt
      setenv REPORT4 $REPSTAT/p.stackparm.$ilabl.rpt
      setenv REPORT5 $REPSTAT/p.state_scc.$ilabl.rpt
      setenv REPORT6 $REPSTAT/pg.state.$GRID.$ilabl.rpt
      setenv REPORT7 $REPSTAT/pg.scc.$GRID.$ilabl.rpt
      breaksw
      
   endsw
   
   set logabbr = inv.$ilabl
   
   breaksw
   
case temporal:

   switch ( $SMK_SOURCE )
   case A:
      setenv REPCONFIG $INVDIR/other/repconfig.ar.temporal.txt
      setenv REPORT1  $REPSCEN/at.state.$slabl.rpt
      setenv REPORT2  $REPSCEN/at.county.$slabl.rpt
      setenv REPORT3  $REPSCEN/at.scc.$slabl.rpt
      setenv REPORT4  $REPSCEN/at.state_scc.$slabl.rpt
      setenv REPORT5  $REPSCEN/at.hour_scc.$slabl.rpt
      setenv REPORT6  $REPSCEN/agt.state.$GRID.$slabl.rpt
      setenv REPORT7  $REPSCEN/ats.scc.$slabl.rpt
      setenv REPORT8  $REPSCEN/ats.state_scc.$slabl.rpt
      breaksw
      
   case M:
      setenv REPCONFIG $INVDIR/other/repconfig.mb.temporal.txt
      setenv REPORT1  $REPSCEN/mt.state.$slabl.rpt
      setenv REPORT2  $REPSCEN/mt.county.$slabl.rpt
      setenv REPORT3  $REPSCEN/mt.scc.$slabl.rpt
      setenv REPORT4  $REPSCEN/mt.state_scc.$slabl.rpt
      setenv REPORT5  $REPSCEN/mt.hour_scc.$slabl.rpt
      setenv REPORT6  $REPSCEN/mt.state_rclas.$slabl.rpt
      setenv REPORT7  $REPSCEN/mt.county_rclas.$slabl.rpt
      setenv REPORT8  $REPSCEN/mgt.state.$GRID.$slabl.rpt
      setenv REPORT9  $REPSCEN/mts.scc.$slabl.rpt
      setenv REPORT10 $REPSCEN/mts.state_scc.$slabl.rpt
      breaksw
      
   case P:
      setenv REPCONFIG $INVDIR/other/repconfig.pt.temporal.txt
      setenv REPORT1  $REPSCEN/pt.state.$slabl.rpt
      setenv REPORT2  $REPSCEN/pt.county.$slabl.rpt
      setenv REPORT3  $REPSCEN/pt.scc.$slabl.rpt
      setenv REPORT4  $REPSCEN/pt.state_scc.$slabl.rpt
      setenv REPORT5  $REPSCEN/pt.hour_scc.$slabl.rpt
      setenv REPORT6  $REPSCEN/pgt.state.$GRID.$slabl.rpt
      setenv REPORT7  $REPSCEN/pgt.source.$GRID.$slabl.rpt
      setenv REPORT8  $REPSCEN/pts.scc.$slabl.rpt
      setenv REPORT9  $REPSCEN/pts.state_scc.$slabl.rpt
      breaksw
      
   endsw
   
   set logabbr = temporal.$ESDATE.$slabl
   
   breaksw

case elevpoint:
   setenv REPCONFIG $INVDIR/other/repconfig.pt.elev.txt
   setenv REPORT1  $REPSTAT/p.src_elev.$ilabl.rpt

   set logabbr = elevpoint.$ESDATE.$slabl
   
   breaksw

case laypoint:
   setenv REPCONFIG $INVDIR/other/repconfig.pt.lfrac.txt
   setenv REPORT1   $REPSCEN/pt.state_layers.$ilabl.rpt

   set logabbr = laypoint.$ESDATE.$slabl
   
   breaksw

case custom:
   set logabbr = custom.$ilabl
   breaksw

case none:
   echo 'SCRIPT NOTE: exiting qa_run.csh without running QA.'
   echo "             QA_TYPE  = "$QA_TYPE
   echo "             RUN_PART1= "$RUN_PART1
   echo "             RUN_PART2= "$RUN_PART2
   echo "             RUN_PART3= "$RUN_PART3
   echo "             RUN_PART4= "$RUN_PART4
   exit( 0 )

default
   echo 'SCRIPT ERROR: Environment variable QA_TYPE is set to'
   echo '              an unknown setting.  Valid values are'
   echo '              "inventory", "temporal", "elevpoint", or "laypoint".'
   set exitstat = 1
   breaksw
      
endsw

# Abort if already had an error
if ( $exitstat == 1 ) then
   exit( $exitstat )
endif
  
#
### Smkreport processing for area, mobile, or point sources
#
set debugexestat = 0
set exestat = 0
setenv TMPLOG   $OUTLOG/smkreport.$SRCABBR.$logabbr.log
if ( $?RUN_SMKREPORT ) then
   if ( $RUN_SMKREPORT == 'Y' && -e $SMK_BIN/smkreport ) then

      if ( -e $TMPLOG ) then
	 source $SCRIPTS/run/movelog.csh
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
               time $SMK_BIN/smkreport
            else
               set exestat = 1 
            endif
         endif
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

   endif
endif

#
## Ending of script with exit status
#
exit( $exitstat )

