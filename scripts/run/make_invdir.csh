#!/bin/csh -f

# Script for automatically setting E.V.s for Smkinven to build output
#   file names and directories.

      switch( $SMK_SOURCE )
      case A:
         setenv INVNAME1 area
         setenv INVNAME2 asrc

         if ( $?NONROAD ) then
            if ( $NONROAD == Y ) then
               setenv INVNAME1 nroad
               setenv INVNAME2 nrsrc
            endif
         endif
      breaksw

      case M:
         setenv INVNAME1 mobl
         setenv INVNAME2 msrc
      breaksw

      case P:
         setenv INVNAME1 pnts
         setenv INVNAME2 psrc
      breaksw

      default:
         echo "SCRIPT ERROR: SMK_SOURCE value of $SMK_SOURCE not"
         echo "              recognized in smk_run.csh"
         set exitstat = 1
      endsw

      if ( $?FYEAR ) then
         set fyr2 = `echo $FYEAR | cut -c3-4`
         setenv INVNAME1 ${INVNAME1}_$fyr2
      endif

      if ( $?CNTLCASE ) then
         setenv INVNAME1 ${INVNAME1}_$CNTLCASE
      endif

      if ( -e $INVOPD ) then
         mkdir -p $INVOPD/${INVNAME1}_dat
      endif

exit( 0 )
