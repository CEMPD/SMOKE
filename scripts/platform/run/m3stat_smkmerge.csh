#!/bin/tcsh -f
##
## Wrapper script to run m3stat on smkmerge output files.
## Assumes that is being called from a script that has just run
## smkmerge. Sector independent.  
##
## These code blocks were copied from the relevant sector scripts.
##
## Script created by: A. Zubrow, UNC,  Jan 2009 
##
##*********************************************************************

set exitstat = 0  ## successful run 

## actual script to call M3STAT function
set m3stat       = $SCRIPTS/run/m3stat_chk_v6.csh

## Make sure that the EMF time period is set
if ( $?EMF_PERIOD == 0 ) then
    setenv EMF_PERIOD "" 
endif

## If EMF_CLIENT is not defined, assume that the user has not setup
#    the assigns file to use the new EMF-ready helper scripts
if ( ! $?EMF_CLIENT ) then
   setenv EMF_CLIENT false
   setenv EMF_JOBKEY  ""
   echo 'NOTE: EMF_CLIENT setting assumed to be "false" in' !$
endif


echo "Running m3stat on smkmerge file(s) for $ESDATE"

## Depending on source category [A|P|M] run m3stat on different files
if ( $MRG_SOURCE == A ) then
      # Run m3stat script on Smkmerge output file
      ## C. Allen configured for multiple files per day (All of my comments contain [c])
      if ( $?AOUT && $RUN_M3STAT == Y ) then
      
         # [c] If $AOUT exists, this implies there is only one Smkmerge file per day. In this case, run m3stat once and get out
	 if ( -e $AOUT ) then
            $m3stat $AOUT
            if ( $status != 0 ) then
   	       echo "ERROR: running m3stat_chk for $ESDATE"
	       $EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: running m3stat_chk for $ESDATE" -t "e" -x $m3stat -p $EMF_PERIOD
               set exitstat = 1
	       goto end_of_script
            endif
	 # [c] Otherwise, check to see how many files there are, and run m3stat on each one
	 else
	    # [c] "aout_prefix" is prefix of $AOUT, everything before the file number
	    set aout_prefix = $INTERMED/emis_${MM}_${SECTOR}_${ESDATE}_${GRID}_${SPC}_${CASE}
	    set num_smkmerge_files = `/bin/ls -1 $aout_prefix.*.ncf | wc -l`
	    echo "SCRIPT NOTE: Running m3stat on $num_smkmerge_files model-ready files"
	    set fc = 0
	    # [c] This allows for case where number of files is 0; in that case the while loop is skipped altogether
	    # [c] The run_m3stat script has an optional 3rd command-line parameter that I originally added to handle inline approach;
	    # [c]   this is a "file stamp" appended to the file names so that the files don't get overwritten each time
	    while ( $fc < $num_smkmerge_files )
 	       @ fc = $fc + 1
	       $m3stat $aout_prefix.$fc.ncf $fc
               if ( $status != 0 ) then
      	          echo "ERROR: running m3stat_chk for $ESDATE"
	          $EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: running m3stat_chk for $ESDATE" -t "e" -x $m3stat -p $EMF_PERIOD
		  set exitstat = 1
		  goto end_of_script
               endif
	    end # while
	 endif # -e $AOUT
      endif # $?AOUT

else if ( $MRG_SOURCE == M ) then
      # Run m3stat script on Smkmerge output file
      ## C. Allen configured for multiple files per day (All of my comments contain [c])
      if ( $?MOUT && $RUN_M3STAT == Y ) then
      
         # [c] If $MOUT exists, this implies there is only one Smkmerge file per day. In this case, run m3stat once and get out
	 if ( -e $MOUT ) then
            $m3stat $MOUT
            if ( $status != 0 ) then
   	       echo "ERROR: running m3stat_chk for $ESDATE"
	       $EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: running m3stat_chk for $ESDATE" -t "e" -x $m3stat -p $EMF_PERIOD
               set exitstat = 1
	       goto end_of_script
            endif
	 # [c] Otherwise, check to see how many files there are, and run m3stat on each one
	 else
	    # [c] "mout_prefix" is prefix of $MOUT, everything before the file number
	    set mout_prefix = $INTERMED/emis_${MM}_${SECTOR}_${ESDATE}_${GRID}_${SPC}_${CASE}
	    set num_smkmerge_files = `/bin/ls -1 $mout_prefix.*.ncf | wc -l`
	    echo "SCRIPT NOTE: Running m3stat on $num_smkmerge_files model-ready files"
	    set fc = 0
	    # [c] This allows for case where number of files is 0; in that case the while loop is skipped altogether
	    # [c] The run_m3stat script has an optional 3rd command-line parameter that I originally added to handle inline approach;
	    # [c]   this is a "file stamp" appended to the file names so that the files don't get overwritten each time
	    while ( $fc < $num_smkmerge_files )
 	       @ fc = $fc + 1
	       $m3stat $mout_prefix.$fc.ncf $fc
               if ( $status != 0 ) then
      	          echo "ERROR: running m3stat_chk for $ESDATE"
	          $EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: running m3stat_chk for $ESDATE" -t "e" -x $m3stat -p $EMF_PERIOD
                  set exitstat = 1
                  goto end_of_script
               endif
	    end # while
	 endif # -e $MOUT
      endif # $?MOUT

else if ( $MRG_SOURCE == P ) then
      # Run m3stat script on POUT, unless we're only creating INLN file
      if ( $INLINE_MODE != only ) then
         # Run m3stat script on Smkmerge output file
         ## C. Allen configured for multiple files per day (All of my comments contain [c])
         if ( $?POUT && $RUN_M3STAT == Y ) then
      
            # [c] If $POUT exists, this implies there is only one Smkmerge file per day. In this case, run m3stat once and get out
	    if ( -e $POUT ) then
               $m3stat $POUT
               if ( $status != 0 ) then
   	          echo "ERROR: running m3stat_chk for $ESDATE"
	          $EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: running m3stat_chk for $ESDATE" -t "e" -x $m3stat -p $EMF_PERIOD
                  set exitstat = 1
	          goto end_of_script
               endif
	    # [c] Otherwise, check to see how many files there are, and run m3stat on each one
	    else
	       # [c] "pout_prefix" is prefix of $POUT, everything before the file number
	       set pout_prefix = $INTERMED/emis_${MM}_${SECTOR}_${ESDATE}_${GRID}_${SPC}_${CASE}
	       set num_smkmerge_files = `/bin/ls -1 $pout_prefix.*.ncf | wc -l`
	       echo "SCRIPT NOTE: Running m3stat on $num_smkmerge_files model-ready files"
	       set fc = 0
	       # [c] This allows for case where number of files is 0; in that case the while loop is skipped altogether
	       # [c] The run_m3stat script has an optional 3rd command-line parameter that I originally added to handle inline approach;
   	       # [c]   this is a "file stamp" appended to the file names so that the files don't get overwritten each time
               while ( $fc < $num_smkmerge_files )
 	          @ fc = $fc + 1
   	          $m3stat $pout_prefix.$fc.ncf $fc
                  if ( $status != 0 ) then
                     echo "ERROR: running m3stat_chk for $ESDATE"
	             $EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: running m3stat_chk for $ESDATE" -t "e" -x $m3stat -p $EMF_PERIOD
                     set exitstat = 1
                     goto end_of_script
                  endif
   	       end # while
	    endif # -e $POUT
         endif # $?POUT
      endif

      # Run m3stat script on INLN, unless we're not running inline approach
      if ( $INLINE_MODE != off ) then
         if ( $?INLN && $RUN_M3STAT == Y ) then
            $m3stat $INLN inln
            if ( $status != 0 ) then
	       echo "ERROR: running $m3stat for $ESDATE, INLN"
	       $EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: running m3stat_chk for $ESDATE, INLN" -t "e" -x $m3stat -p $EMF_PERIOD
               set exitstat = 1
	       goto end_of_script
            endif
         endif
      endif


else
      echo "SCRIPT ERROR: Unknown merge source $MRG_SOURCE"
      exit (1)
endif

###########################################################
# Label for the end of the script, used during script abort
end_of_script:

## successfully completed
exit($exitstat)
