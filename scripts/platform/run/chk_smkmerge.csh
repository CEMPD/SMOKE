#!/bin/tcsh -f
##
## Wrapper script to test for existence of smkmerge input files.
## Assumes that is being called from a script that is about to run
## smkmerge (i.e. the dates have been correctly set). Sector independent.  
##
## Script created by: A. Zubrow, UNC,  Jan 2009 
##
##*********************************************************************

set exitstat = 0  ## successful run 

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


echo "Checking for smkmerge input file(s) for $ESDATE"

## Depending on source category [A|P|M] check for different files
if ( $MRG_SOURCE == A ) then
    if (! -e $AREA ) then
	echo "ERROR: Necessary intermediary file does not exist: $AREA"
	set exitstat = 1
	$EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: Necessary intermediary file does not exist: $AREA" -t "e"
    endif
    if (! -e $ATMP ) then
	echo "ERROR: Necessary intermediary file does not exist: $ATMP"
	set exitstat = 1
	$EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: Necessary intermediary file does not exist: $ATMP" -t "e"
    endif
    if (! -e $ASMAT_L ) then
	echo "ERROR: Necessary intermediary file does not exist: $ASMAT_L"
	set exitstat = 1
	$EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: Necessary intermediary file does not exist: $ASMAT_L" -t "e"
    endif
    if (! -e $ASMAT_S ) then
	echo "ERROR: Necessary intermediary file does not exist: $ASMAT_S"
	set exitstat = 1
	$EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: Necessary intermediary file does not exist: $ASMAT_S" -t "e"
    endif
    if (! -e $AGMAT ) then
	echo "ERROR: Necessary intermediary file does not exist: $AGMAT"
	set exitstat = 1
	$EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: Necessary intermediary file does not exist: $AGMAT" -t "e"
    endif

else if ( $MRG_SOURCE == M ) then
    if (! -e $MOBL ) then
	echo "ERROR: Necessary intermediary file does not exist: $MOBL"
	set exitstat = 1
	$EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: Necessary intermediary file does not exist: $MOBL" -t "e"
    endif
    if (! -e $MTMP ) then
	echo "ERROR: Necessary intermediary file does not exist: $MTMP"
	set exitstat = 1
	$EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: Necessary intermediary file does not exist: $MTMP" -t "e"
    endif
    if (! -e $MSMAT_L ) then
	echo "ERROR: Necessary intermediary file does not exist: $MSMAT_L"
	set exitstat = 1
	$EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: Necessary intermediary file does not exist: $MSMAT_L" -t "e"
    endif
    if (! -e $MSMAT_S ) then
	echo "ERROR: Necessary intermediary file does not exist: $MSMAT_S"
	set exitstat = 1
	$EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: Necessary intermediary file does not exist: $MSMAT_S" -t "e"
    endif
    if (! -e $MGMAT ) then
	echo "ERROR: Necessary intermediary file does not exist: $MGMAT"
	set exitstat = 1
	$EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: Necessary intermediary file does not exist: $MGMAT" -t "e"
    endif

else if ( $MRG_SOURCE == P ) then
    if (! -e $PNTS ) then
	echo "ERROR: Necessary intermediary file does not exist: $PNTS"
	set exitstat = 1
	$EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: Necessary intermediary file does not exist: $PNTS" -t "e"
    endif
    if (! -e $PTMP ) then
	echo "ERROR: Necessary intermediary file does not exist: $PTMP"
	set exitstat = 1
	$EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: Necessary intermediary file does not exist: $PTMP" -t "e"
    endif
    if (! -e $PSMAT_L ) then
	echo "ERROR: Necessary intermediary file does not exist: $PSMAT_L"
	set exitstat = 1
	$EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: Necessary intermediary file does not exist: $PSMAT_L" -t "e"
    endif
    if (! -e $PSMAT_S ) then
	echo "ERROR: Necessary intermediary file does not exist: $PSMAT_S"
	set exitstat = 1
	$EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: Necessary intermediary file does not exist: $PSMAT_S" -t "e"
    endif
    if (! -e $PGMAT ) then
	echo "ERROR: Necessary intermediary file does not exist: $PGMAT"
	set exitstat = 1
	$EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: Necessary intermediary file does not exist: $PGMAT" -t "e"
    endif

    ## depending 2-D or 3-D emission test for PLAY or PELV
    ## default INLINE_MODE to off
    if (! $?INLINE_MODE ) then
	setenv INLINE_MODE "off"
    endif
    if ( $INLINE_MODE == "off" ) then
	## need PLAY
	if (! -e $PLAY ) then
	    ## check if zip version of file exists
	    if (! -e ${PLAY}.gz ) then
		echo "ERROR: Necessary intermediary file does not exist: $PLAY"
		set exitstat = 1
		$EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: Necessary intermediary file does not exist: $PLAY" -t "e"
	    else
		## gzip version of file exists, unzip to new file
		## retain zip file
		echo "Unzipping PLAY file"
		gunzip -c ${PLAY}.gz > $PLAY

		## set special exit status indicating unzipped
		set exitstat = 99
	    endif

	endif
    else
	## need PELV
	if (! -e $PELV ) then
	    echo "ERROR: Necessary intermediary file does not exist: $PELV"
	    set exitstat = 1
	    $EMF_CLIENT -k $EMF_JOBKEY -m "ERROR: Necessary intermediary file does not exist: $PELV" -t "e"
	endif
    endif

else
      echo "SCRIPT ERROR: Unknown merge source $MRG_SOURCE"
      exit (1)
endif

## return exit status
exit($exitstat)
