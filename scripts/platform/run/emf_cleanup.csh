#!/bin/tcsh -f
# Script for automatically moving old EMF scripts for the same job as the
# one being run when this helper script is called

## Check for all needed environment variables. If any aren't available,
## print a warning message and abort

set exitstat = 0
if ( ! $?EMF_JOBNAME ) then
   echo "SCRIPT WARNING: Could not move old scripts because EMF_JOBNAME was not set"
   set exitstat = 1
endif

if ( ! $?EMF_SCRIPTDIR ) then
   echo "SCRIPT WARNING: Could not move old scripts because EMF_SCRIPTDIR was not set"
   set exitstat = 1
endif

if ( ! $?EMF_SCRIPTNAME ) then
   echo "SCRIPT WARNING: Could not move old scripts because EMF_SCRIPTNAME was not set"
   set exitstat = 1
endif

if ( ! $?EMF_LOGNAME ) then
   echo "SCRIPT WARNING: Could not move old scripts because EMF_LOGNAME was not set"
   set exitstat = 1
endif

if ( ! $?SECTOR ) then
   echo "SCRIPT WARNING: Could not move old scripts because SECTOR was not set"
   set exitstat = 1
endif

if ( ! $?CASE ) then
   echo "SCRIPT WARNING: Could not move old scripts because CASE was not set"
   set exitstat = 1
endif

if ( ! $?GRID ) then
   echo "SCRIPT WARNING: Could not move old scripts because GRID was not set"
   set exitstat = 1
endif

if ( $exitstat > 0 ) then
   exit ( $exitstat )
endif

## Create the directory to move the old scripts to and change the permissions
## to group-writable
mkdir -p $EMF_SCRIPTDIR/old/$SECTOR
chmod g+rwx $EMF_SCRIPTDIR/old 
chmod g+rwx $EMF_SCRIPTDIR/old/$SECTOR

## Create the directory to move the old logs to and change the permissions 
## to group writable
mkdir -p $EMF_SCRIPTDIR/old_logs/$SECTOR
chmod g+rwx $EMF_SCRIPTDIR/old_logs 
chmod g+rwx $EMF_SCRIPTDIR/old_logs/$SECTOR

## Move the old scripts
## Old job naming (unique on name)
set num = `/bin/ls -1 "$EMF_SCRIPTDIR/${EMF_JOBNAME}_${CASE}"* | wc -l`
if ( $num > 0 ) then
   foreach f ( "$EMF_SCRIPTDIR/${EMF_JOBNAME}_${CASE}"* )
      if ( "$f" != $EMF_SCRIPTNAME ) then
         /bin/mv $f $EMF_SCRIPTDIR/old/$SECTOR
         if ( $status != 0 ) set exitstat = 1
      endif
   end

else
    ## New job uniqueness (unique on name, sector, grid), the script names are EMF_JOBNAME_SECTOR_GRID_CASE
    set num = `/bin/ls -1 "$EMF_SCRIPTDIR/${EMF_JOBNAME}_${SECTOR}_${GRID}_${CASE}"* | wc -l`
    if ( $num > 0 ) then
	foreach f ( "$EMF_SCRIPTDIR/${EMF_JOBNAME}_${SECTOR}_${GRID}_${CASE}"* )
	    if ( "$f" != $EMF_SCRIPTNAME ) then
		/bin/mv $f $EMF_SCRIPTDIR/old/$SECTOR
		if ( $status != 0 ) set exitstat = 1
	    endif
	end
    endif

endif

## Move the old logs
## Old job log naming (unique on name)
set num = `/bin/ls -1 "$EMF_SCRIPTDIR/logs/${EMF_JOBNAME}_${CASE}"* | wc -l`
if ( $num > 0 ) then
   foreach f ( "$EMF_SCRIPTDIR/logs/${EMF_JOBNAME}_${CASE}"* )
      if ( "$f" != $EMF_LOGNAME ) then
         /bin/mv $f $EMF_SCRIPTDIR/old_logs/$SECTOR
         if ( $status != 0 ) set exitstat = 1
      endif
   end
else
    ## New job log uniqueness (unique on name, sector, grid), the script names are EMF_JOBNAME_SECTOR_GRID_CASE
    set num = `/bin/ls -1 "$EMF_SCRIPTDIR/logs/${EMF_JOBNAME}_${SECTOR}_${GRID}_${CASE}"* | wc -l`
    if ( $num > 0 ) then
	foreach f ( "$EMF_SCRIPTDIR/logs/${EMF_JOBNAME}_${SECTOR}_${GRID}_${CASE}"* )
	    if ( "$f" != $EMF_LOGNAME ) then
		/bin/mv $f $EMF_SCRIPTDIR/old_logs/$SECTOR
		if ( $status != 0 ) set exitstat = 1
	    endif
	end
    endif

endif
  
exit ( $exitstat )
