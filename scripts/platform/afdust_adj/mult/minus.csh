#!/bin/csh

setenv DATADIR    /ptmp/pou/



setenv INFILE $DATADIR/aqm_pollutant_20070526_agts.CONUS.nofires.ncf
setenv LOGFILE     minus.aqm_pollutant_20070526_agts.CONUS.nofires.ncf.log
setenv OUTFILE   $DATADIR/minus.aqm_pollutant_20070526_agts.CONUS.nofires.ncf




if ( -e $LOGFILE) then
   rm -f $LOGFILE
endif
if ( -e $OUTFILE) then
   rm $OUTFILE
endif

minus.x
