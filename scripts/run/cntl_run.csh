#!/bin/csh -fx

# Version @(#)$Id$
# Path    $Source$
# Date    $Date$

# This script runs the SMOKE control/projection and inventory growth program
#
# Script created by : M. Houyoux , CEP Environmental Modeling Center 
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
   set debug_exe = dbx
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

#
### Control matrix generation
#
set debugexestat = 0
set exestat = 0
setenv TMPLOG   $OUTLOG/cntlmat.$SRCABBR.$INVEN.log
if ( $?RUN_CNTLMAT ) then
   if ( $RUN_CNTLMAT == Y && $RUN_PART1 == Y ) then

      ## Move log
      if ( -e $TMPLOG ) then
         source $SCRIPTS/run/movelog.csh
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

      ## Define name for input file
      switch ( $SMK_SOURCE )
         case A:
            set cntl_indir = $ARDAT
            if ( $NONROAD == Y ) set cntl_indir = $NRDAT
            breaksw
         case M:
            set cntl_indir = $MBDAT
            breaksw
         case P:
            set cntl_indir = $PTDAT
            breaksw
      endsw

      if ( $?FYEAR ) then
         if ( $?CNTLCASE ) then
            setenv GCNTL $cntl_indir/gcntl.${YEAR}_${FYEAR}_$CNTLCASE.txt
         else
            setenv GCNTL $cntl_indir/gcntl.${YEAR}_$FYEAR.txt
         endif
      else
         if ( $?CNTLCASE ) then
            setenv GCNTL $cntl_indir/gcntl.$CNTLCASE.txt
         else
            echo "SCRIPT ERROR: At least FYEAR or CNTLCASE variables must be"
            echo "              set to set GCNTL file name for Cntlmat program."
            exit( 1 )
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
               time $SMK_BIN/cntlmat
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

   endif
endif

#
### Inventory growth and control
#
setenv TMPLOG   $OUTLOG/grwinven.$SRCABBR.$INVEN.log
if ( $?RUN_GRWINVEN ) then
   if ( $RUN_GRWINVEN == Y && $RUN_PART1 == Y ) then

      ## Move log
      if ( -e $TMPLOG ) then
         source $SCRIPTS/run/movelog.csh
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
      source $SCRIPTS/run/make_invdir.csh
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
               time $SMK_BIN/grwinven
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

   endif
endif

#
## Ending of script with exit status
#
exit( $exitstat )

