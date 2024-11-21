#!/bin/csh -f
## log w/ EMF server that script is running
$EMF_CLIENT -k $EMF_JOBKEY -s "Running" 

setenv GRID "12US2"
## source the ASSIGN file
source $ASSIGNS_FILE

setenv OUTPUT_DIR "$PROJECT_ROOT/$CASE/aermod"
setenv WORK_PATH $OUTPUT_DIR

if (! -d $OUTPUT_DIR/qa/spreadsheets) then
    mkdir -p $OUTPUT_DIR/qa/spreadsheets
endif

## List of all the helper scripts that are run in this script
set emf_cleanup  = $SCRIPTS/run/emf_cleanup.csh
set point_qa = $SCRIPTS/emf/aermod/qa/point_qa.py
set point_loc = $SCRIPTS/emf/aermod/qa/point_loc.py
set xlsx = $SCRIPTS/emf/aermod/qa/gen_xlsx.py

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

echo "SCRIPT NOTE: Generating Point style AERMOD helper QA"
## Set naming label
set namelabel = ${SECTOR}_${CASE}_${GRID}

$point_qa
$point_loc
$xlsx $SECTOR $CASE $OUTPUT_DIR 

## Ending of script
#
exit( $exitstat )
