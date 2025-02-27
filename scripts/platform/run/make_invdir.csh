#!/bin/tcsh -f

# Script for automatically setting E.V.s for Smkinven to build output
#   file names and directories.

      switch( $SMK_SOURCE )
      case A:
         setenv INVNAME1 area_${SUBSECT}
         setenv INVNAME2 asrc_${SUBSECT}

      breaksw

      case M:
         setenv INVNAME1 mobl_${SUBSECT}
         setenv INVNAME2 msrc_${SUBSECT}
      breaksw

      case P:
         setenv INVNAME1 pnts_${SUBSECT}
         setenv INVNAME2 psrc_${SUBSECT}
      breaksw

      default:
         echo "SCRIPT ERROR: SMK_SOURCE value of $SMK_SOURCE not"
         echo "              recognized in smk_run.csh"
         set exitstat = 1
      endsw

      if ( -e $INTERMED ) then
         mkdir -p $INTERMED/${INVNAME1}_dat
      endif

exit( 0 )
