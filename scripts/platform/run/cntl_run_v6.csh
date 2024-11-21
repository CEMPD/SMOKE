#!/bin/csh -f

# Version @(#)$Id: cntl_run.csh,v 1.1.2.6 2003/03/14 14:12:43 mhouyoux Exp $
# Path    $Source: /afs/isis/depts/cep/emc/apps/archive/smoke/smoke/scripts/run/Attic/cntl_run.csh,v $
# Date    $Date: 2003/03/14 14:12:43 $

# This script runs the SMOKE control/projection and inventory growth program
#
# Script created by : M. Houyoux , CEP Environmental Modeling Center 
# Modified: A. Zubrow, UNC, for custom pre-defined GCNTL  Dec, 2008
# 
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
      echo 'NOTE: DEBUG_EXE setting assumed to be "dbx" in cntl_run.csh'
   endif
endif

## Make sure that the EMF time period is set
if ( ! $?EMF_PERIOD ) then
    set emf_period = ""
else
    set emf_period = $EMF_PERIOD
endif

## If EMF_CLIENT is not defined, assume that the user has not setup
#    the assigns file to use the new EMF-ready helper scripts
if ( ! $?EMF_CLIENT ) then
   setenv EMF_CLIENT false
   echo 'NOTE: EMF_CLIENT setting assumed to be "false" in cntl_run.csh'
endif

## Check the CUSTOM GCNTL flag is set
if ( ! $?CUSTOM_GCNTL ) then
    setenv CUSTOM_GCNTL N
endif

## helper scripts
set movelog = $SCRIPTS/run/movelog.csh
set make_invdir = $SCRIPTS/run/make_invdir.csh
set timetracker = $SCRIPTS/run/timetracker_v2.csh 
set checklog = $SCRIPTS/run/checklogfile.csh


### Ensure new controller variables are set
if ( $?RUN_PART1B ) then
   if ( $RUN_PART1B == Y || $RUN_PART1B == y ) then
      setenv RUN_PART1B Y
      echo 'Running part 1b...'
   endif
else
   setenv RUN_PART1B N
endif

#
### Control matrix generation
#
set debugexestat = 0
set exestat = 0
setenv TMPLOG   $OUTLOG/cntlmat.$SRCABBR.$CASE.log
if ( $?RUN_CNTLMAT ) then
   if ( $RUN_CNTLMAT == Y && $RUN_PART1B == Y ) then

      ## Move log
      if ( -e $TMPLOG ) then
         source $movelog
      endif

      # Make temporary control file directory and make sure that it is
      #    writable by the user
      if ( ! -e $SMK_TMPDIR ) then
          mkdir -p $SMK_TMPDIR
          chmod ug+w $SMK_TMPDIR

      else
	 set line   = ( `/bin/ls -ld $SMK_TMPDIR` )
	 set owner  = $line[3]
	 set check  = ( `echo $line | grep $user` )
	 if ( $status == 0 ) then

             set permis = ( `echo $line | grep drwxrw` )
             if( $status != 0 ) then
        	 chmod ug+w $SMK_TMPDIR
             endif

	 else

             set permis = ( `echo $line | grep drwxrw` )
             if ( $status != 0 ) then
        	 echo "NOTE: Do not have write permission for temporary control"
                 echo "      file directory:"
        	 echo "      $SMK_TMPDIR"
        	 echo "      Check with user $owner for write permissions."
                 set exitstat = 1
             endif

	 endif
      endif

      ## if custom GCNTL, use predefined input file
      if ($CUSTOM_GCNTL == Y) then
	 ## test if GCNTL is defined
	 if (! $?GCNTL) then
	    echo "SCRIPT ERROR: GCNTL must be pre-defined if using CUSTOM_GCNTL"
	    exit (1)
	 endif
	 echo "Using CUSTOM GCNTL = $GCNTL"
         setenv PCMAT    $INTERMED/pcmat_${SECTOR}_${CASE}.ncf
         setenv PRMAT_L  $INTERMED/prmat_l_${SECTOR}_${SPC}_$CASE.ncf  # Reactivity matrix, mole based
         setenv PRMAT_S  $INTERMED/prmat_s_${SECTOR}_${SPC}_$CASE.ncf  # Reactivity matrix, mass based
         setenv PRSUP    $INTERMED/prsup_${SECTOR}_${CASE}.txt         # Reactivity ASCII part
         setenv PPMAT    $INTERMED/ppmat_${SECTOR}_${CASE}.ncf         # Growth matrix
      else
      ## Define name for input file
	 set cntl_indir = $INVDIR/$SECTOR

	 if ( $?FYEAR ) then
	    if ( $?CNTLCASE ) then
		if (! $?GCNTL) then
		    setenv GCNTL $cntl_indir/gcntl_${YEAR}_${FYEAR}_$CNTLCASE.txt
		endif
	    else
		setenv GCNTL $cntl_indir/gcntl_${YEAR}_$FYEAR.txt
	    endif
	 else
	    if ( $?CNTLCASE ) then
		setenv GCNTL $cntl_indir/gcntl_$CNTLCASE.txt
	    else
		echo "SCRIPT ERROR: At least FYEAR or CNTLCASE variables must be"
		echo "              set to set GCNTL file name for Cntlmat program."
		exit( 1 )
	    endif
	 endif
      endif
      
      ## Remove any existing tmp files
      if ( -e $SMK_TMPDIR/filelist.txt ) /bin/rm -rf $SMK_TMPDIR/filelist.txt
      ls -1 $SMK_TMPDIR > $SMK_TMPDIR/filelist.txt
      set list = ( `cat $SMK_TMPDIR/filelist.txt | grep cntlmat` )
      if ( $status == 0 ) then
         /bin/rm -rf $SMK_TMPDIR/cntlmat*
      endif

      ## Run program
      if ( $exitstat == 0 ) then
         setenv LOGFILE $TMPLOG

         if ( $debugmode == Y ) then
            if ( -e $CL_SRC/cntlmat.debug ) then
               $debug_exe $CL_SRC/cntlmat.debug
            else
                set debugexestat = 1
            endif
         else
            if ( -e $SMK_BIN/cntlmat ) then

               set startdt = `date +%m/%d/%Y,%T`
#	       $EMF_CLIENT -k $EMF_JOBKEY -x $SMK_BIN/cntlmat -m "Running Cntlmat -- $emf_period" -p $emf_period ## Log to EMF server
               time $SMK_BIN/cntlmat
               if ( $timelog_yn == Y ) then
                  $timetracker N $TIMELOG $startdt cntlmat $ESDATE
                  if ( $status != 0 ) then
                      echo "ERROR: Problem calling timetracker from cntl_run script"
                      $EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: Problem calling timetracker from cntl_run script" -x $timetracker  -t "e" -p $emf_period ## log w/ EMF server
                      exit ( 1 )
                  endif  
               endif 

               $checklog
               if ( $status != 0 ) then
                 echo ERROR detected in Cntlmat
		 $EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: detected in Cntlmat" -p $emf_period -t "e" ## Log to EMF Server
                 exit( 1 )
               endif
            else
               set exestat = 1
            endif
         endif
      endif

      ## Move any info from fort.99 that may have been written by I/O API
      if ( -e $SCRIPTS/fort.99 ) then
         mv $LOGFILE $LOGFILE.tmp
         cat $LOGFILE.tmp $SCRIPTS/fort.99 > $LOGFILE
         /bin/rm -rf $LOGFILE.tmp
         /bin/rm -rf $SCRIPTS/fort.99
      endif

      if ( $exestat == 1 ) then
         echo 'SCRIPT ERROR: cntlmat program does not exist in:'
         echo '              '$SMK_BIN
         set exitstat = 1
      endif

      if ( $debugexestat == 1 ) then
         echo 'SCRIPT ERROR: cntlmat.debug program does not exist in:'
         echo '              '$CL_SRC
         set exitstat = 1
      endif

      ## Remove any existing tmp files
      if ( -e $SMK_TMPDIR/filelist.txt )  /bin/rm -rf $SMK_TMPDIR/filelist.txt
      ls -1 $SMK_TMPDIR > $SMK_TMPDIR/filelist.txt
      set list = ( `cat $SMK_TMPDIR/filelist.txt | grep cntlmat` )
      if ( $status == 0 ) then
#        /bin/rm -rf $SMK_TMPDIR/cntlmat*
      endif
      if ( $exitstat == 1 ) then
         $EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: script error running Cntlmat" -p $emf_period -t "e" ## Log to EMF Server
      endif

   endif
endif

#
### Inventory growth and control
#
setenv TMPLOG   $OUTLOG/grwinven.$SRCABBR.$CASE.log
if ( $?RUN_GRWINVEN ) then
   if ( $RUN_GRWINVEN == Y && $RUN_PART1B == Y ) then

      ## Move log
      if ( -e $TMPLOG ) then
         source $movelog
      endif

      # Make temporary control file directory and make sure that it is
      #    writable by the user
      if ( ! -e $SMK_TMPDIR ) then
          mkdir -p $SMK_TMPDIR
          chmod ug+w $SMK_TMPDIR

      else
	 set line   = ( `/bin/ls -ld $SMK_TMPDIR` )
	 set owner  = $line[3]
	 set check  = ( `echo $line | grep $user` )
	 if ( $status == 0 ) then

             set permis = ( `echo $line | grep drwxrw` )
             if( $status != 0 ) then
        	 chmod ug+w $SMK_TMPDIR
             endif

	 else

             set permis = ( `echo $line | grep drwxrw` )
             if ( $status != 0 ) then
        	 echo "NOTE: Do not have write permission for temporary control"
                 echo "      file directory:"
        	 echo "      $SMK_TMPDIR"
        	 echo "      Check with user $owner for write permissions."
                 set exitstat = 1
             endif

	 endif
      endif

      ##  Create output directories, if needed
      source $make_invdir
      set exitstat = $status

      ## Remove any existing tmp files
      if ( -e $SMK_TMPDIR/filelist.txt )  /bin/rm -rf $SMK_TMPDIR/filelist.txt
      ls -1 $SMK_TMPDIR > $SMK_TMPDIR/filelist.txt
      set list = ( `cat $SMK_TMPDIR/filelist.txt | grep grwinven` )
      if ( $status == 0 ) then
         /bin/rm -rf $SMK_TMPDIR/grwinven*
      endif

      ## Run program
      if ( $exitstat == 0 ) then
         setenv LOGFILE $TMPLOG

         if ( $debugmode == Y ) then
            if ( -e $IV_SRC/grwinven.debug ) then
               $debug_exe $IV_SRC/grwinven.debug
            else
                set debugexestat = 1
            endif
         else
            if ( -e $SMK_BIN/grwinven ) then

               set startdt = `date +%m/%d/%Y,%T`
#	       $EMF_CLIENT -k $EMF_JOBKEY -x $SMK_BIN/grwinven -m "Running Grwinven -- $emf_period" -p $emf_period ## Log to EMF server
               time $SMK_BIN/grwinven
               if ( $timelog_yn == Y ) then
                  $timetracker N $TIMELOG $startdt grwinven $ESDATE
                  if ( $status != 0 ) then
                      echo "ERROR: Problem calling timetracker from cntl_run script"
                      $EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: Problem calling timetracker from cntl_run script" -x $timetracker  -t "e" -p $emf_period ## log w/ EMF server
                      exit ( 1 )
                  endif  
               endif 

               $checklog
               if ( $status != 0 ) then
                 echo ERROR detected in Grwinven
		 $EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: detected in Grwinven" -p $emf_period -t "e" ## Log to EMF Server
                 exit( 1 )
               endif
            else
               set exestat = 1
            endif
         endif
      endif

      ## Move any info from fort.99 that may have been written by I/O API
      if ( -e $SCRIPTS/fort.99 ) then
         mv $LOGFILE $LOGFILE.tmp
         cat $LOGFILE.tmp $SCRIPTS/fort.99 > $LOGFILE
         /bin/rm -rf $LOGFILE.tmp
         /bin/rm -rf $SCRIPTS/fort.99
      endif

      if ( $exestat == 1 ) then
         echo 'SCRIPT ERROR: grwinven program does not exist in:'
         echo '              '$SMK_BIN
         set exitstat = 1
      endif

      if ( $debugexestat == 1 ) then
         echo 'SCRIPT ERROR: grwinven.debug program does not exist in:'
         echo '              '$IV_SRC
         set exitstat = 1
      endif

      ## Remove any existing tmp files
      if ( -e $SMK_TMPDIR/filelist.txt )  /bin/rm -rf $SMK_TMPDIR/filelist.txt
      ls -1 $SMK_TMPDIR > $SMK_TMPDIR/filelist.txt
      set list = ( `cat $SMK_TMPDIR/filelist.txt | grep grwinven` )
      if ( $status == 0 ) then
         /bin/rm -rf $SMK_TMPDIR/grwinven*
      endif
      if ( $exitstat == 1 ) then
         $EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: running Grwinven" -p $emf_period -t "e" ## Log to EMF Server
      endif

   endif
endif

#
## Ending of script with exit status
#
exit( $exitstat )

