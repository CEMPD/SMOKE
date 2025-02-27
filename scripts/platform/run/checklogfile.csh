#!/bin/tcsh
#
##### check for Normal completion in a log file
#
echo now checking log file $LOGFILE
grep 'Normal Completion of program' $LOGFILE > /dev/null
if ( $status != 0 ) then
   echo '*******************************'
   echo '* ERROR detected in logfile:' 
   echo "* "$LOGFILE
   echo '*******************************'
   exit( 1 )
endif
