#!/bin/tcsh -f

# Version @(#)$Id$
# Path    $Source$
# Date    $Date$

# This script sets up needed environment variables for running SMOKE point
# sources for the steps necessary to generate AERMOD helper files.
# The external helper scripts use SMOKE reports and inputs to generate the
# final AERMOD helper files.
#
# This script is intended to be used with the EMF
# source emissions in SMOKE for the EPA 2014 modeling platform, and 
# calls the scripts that runs the SMOKE programs. 
#
# Script created by : M. Houyoux, Environmental Protection Agency
# October, 2007
# Modified for support of AERMOD by: J. Beidler, CSRA March 2017
#
#*********************************************************************

## log w/ EMF server that script is running
$EMF_CLIENT -k $EMF_JOBKEY -s "Running" 
setenv GRID "12US2"
## source the ASSIGN file
source $ASSIGNS_FILE

set emf_cleanup  = $SCRIPTS/run/emf_cleanup.csh
## If running from EMF, move old EMF-created scripts to "old"
if ( $?EMF_JOBID ) then
   source $emf_cleanup
   if ( $status != 0 ) then
	echo "ERROR: running EMF script/log cleanup script"
	$EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: running EMF script/log cleanup script" -t "e" -x $emf_cleanup
	exit( 1 )
   endif
endif

set exitstat = 0

setenv OUTPUT_DIR $PROJECT_ROOT/$CASE/aermod
if ( $?OUTPUT_DIR ) then
    mkdir -p $OUTPUT_DIR
    chmod ug+rwx $OUTPUT_DIR
    foreach subdir (locations temporal parameters xwalk emis)
        mkdir -p $OUTPUT_DIR/$subdir
        chmod ug+rwx $OUTPUT_DIR/$subdir
    end
endif

## Set naming label
set namelabel = ${SECTOR}_${CASE}_${GRID}

# Combine point files
set locout = $OUTPUT_DIR/locations/point_combined_location.csv
rm -f $locout
head -1 $OUTPUT_DIR/locations/point_combined_12US2_location.csv > $locout
grep -hv ^state $OUTPUT_DIR/locations/point_combined_*_location.csv >> $locout
set ppout = $OUTPUT_DIR/parameters/point_combined_point_srcparam.csv
rm -f $ppout
head -1 $OUTPUT_DIR/parameters/point_combined_12US2_point_srcparam.csv > $ppout
grep -hv ^facility_id $OUTPUT_DIR/parameters/point_combined_*_point_srcparam.csv >> $ppout
set fpout = $OUTPUT_DIR/parameters/point_combined_fug_srcparam.csv
rm -f $fpout
head -1 $OUTPUT_DIR/parameters/point_combined_12US2_fug_srcparam.csv > $fpout
grep -hv ^facility_id $OUTPUT_DIR/parameters/point_combined_*_fug_srcparam.csv >> $fpout
set tmpout = $OUTPUT_DIR/temporal/point_combined_temporal.csv
rm -f $tmpout
head -1 $OUTPUT_DIR/temporal/point_combined_12US2_temporal.csv > $tmpout
grep -hv ^facility_id $OUTPUT_DIR/temporal/point_combined_*_temporal.csv >> $tmpout
set xout = $OUTPUT_DIR/xwalk/point_combined_process_releasept_emis.csv
rm -f $xout
head -1 $OUTPUT_DIR/xwalk/point_combined_12US2_process_releasept_emis.csv > $xout
grep -hv ^state $OUTPUT_DIR/xwalk/point_combined_*_process_releasept_emis.csv >> $xout
set srcout = $OUTPUT_DIR/xwalk/point_combined_srcid_xwalk.csv
rm -f $srcout
head -1 $OUTPUT_DIR/xwalk/point_combined_12US2_srcid_xwalk.csv > $srcout
grep -hv ^state $OUTPUT_DIR/xwalk/point_combined_*_srcid_xwalk.csv >> $srcout
echo "SCRIPT NOTE: Generating Point style AERMOD helper files"
$EMF_CLIENT -k $EMF_JOBKEY -m "Generating point style AERMOD helper files"   ## log w/ EMF server
$EMF_CLIENT -k $EMF_JOBKEY -F $OUTPUT_DIR/locations/point_combined_location.csv \
     -T "AERMOD Helper Point Location (CSV)" -O "Point Point locations ${namelabel}"
$EMF_CLIENT -k $EMF_JOBKEY -F $OUTPUT_DIR/parameters/point_combined_point_srcparam.csv \
     -T "AERMOD Helper Point Parameters (CSV)" -O "Point Point parameters ${namelabel}"
$EMF_CLIENT -k $EMF_JOBKEY -F $OUTPUT_DIR/parameters/point_combined_fug_srcparam.csv \
     -T "AERMOD Helper Fugitive Parameters (CSV)" -O "Point Fugitive parameters ${namelabel}"
$EMF_CLIENT -k $EMF_JOBKEY -F $OUTPUT_DIR/temporal/point_combined_temporal.csv \
     -T "AERMOD Helper Point Temporal (CSV)" -O "Point Temporal ${namelabel}"
$EMF_CLIENT -k $EMF_JOBKEY -F $OUTPUT_DIR/xwalk/point_combined_process_releasept_emis.csv \
     -T "AERMOD Helper Point Emissions (CSV)" -O "Point source emissions ${namelabel}"
$EMF_CLIENT -k $EMF_JOBKEY -F $OUTPUT_DIR/xwalk/point_combined_srcid_xwalk.csv \
     -T "AERMOD Helper Point Crosswalk (CSV)" -O "Point source crosswalk ${namelabel}"
   
# Label for the end of the script, used during script abort
end_of_script:
## Ending of script
#
exit( $exitstat )
