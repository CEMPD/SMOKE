#!/bin/csh -f

# Version @(#)$Id$
# Path    $Source$
# Date    $Path$

# This script sets up needed environment variables for creating future-year
# point-source emissions and calls the script that runs the SMOKE programs. 
#
# Script created by : M. Houyoux, MCNC Environmental Modeling Center 
#                     April, 2002
#*********************************************************************

echo " "
echo "NOTE: THIS SCRIPT HAS NOT BEEN INTEGRATED FOR THE SMOKE v1.5 RELEASE"
echo " "

exit( 1 )


# set Assigns file names
setenv ASSIGNS_FILE $SMKROOT/assigns/ASSIGNS.m3demo.cb-iv.ed36

# set future year
setenv FYEAR 2007              # year of future case

# set source category
setenv SMK_SOURCE  P           # source category to process
setenv MRG_SOURCE  P           # source category to merge
setenv SRCABBR     pt.$FYEAR   # abbreviation for naming log files

# time independent programs
setenv RUN_CNTLMAT   Y            # Y runs control matrix program
setenv RUN_GRWINVEN  Y            # Y runs control matrix program

# time-dependent programs
setenv RUN_TEMPORAL  Y            # Y runs temporal allocation program
setenv RUN_SMKMERGE  Y            # Y runs merge program

# quality assurance
setenv RUN_SMKREPORT Y         # Y runs reporting for state reports

# Program-specific controls...

# For Cntlmat
setenv TMP_CTL_PATH         $SMKDAT/tmp # directory for temporary program files
setenv PROJECTION_YR_SPEC   Y           # projection packet entries in year-specific format 
setenv CONTROL_REPORT       Y     # Y outputs a report of factors by source

# For Grwinven
setenv SMK_NUM_CTLMAT       1     # number of control/projection matrices
setenv SMK_GRWSMKOUT_YN     Y     # Y outputs a SMOKE-formatted inventory
setenv SMK_GRWIDAOUT_YN     Y     # Y outputs an IDA-formatted inventory

# For Temporal
setenv RENORM_TPROF         Y     # Y renormalizes temporal profiles
setenv UNIFORM_TPROF_YN     N     # Y makes all temporal profiles uniform
setenv ZONE4WM              Y     # Y uses time zones for start of day & month

# For Smkmerge
setenv MRG_TEMPORAL_YN      Y          # Y merges with hourly emissions
setenv MRG_SPCMAT_YN        Y          # Y merges with speciation matrix
setenv MRG_LAYERS_YN        Y          # Y merges with layer fractions
setenv MRG_GRDOUT_YN        Y          # Y outputs gridded file
setenv MRG_REPSTA_YN        N          # Y outputs state totals
setenv MRG_REPCNY_YN        N          # Y outputs county totals
setenv SMK_ASCIIELEV_YN     N          # Y outputs ASCII elevated file
setenv MRG_GRDOUT_UNIT      moles/s    # units for gridded output file
setenv MRG_TOTOUT_UNIT      moles/day  # units for state and/or county totals
setenv MRG_REPORT_TIME      230000     # hour in OUTZONE for reporting emissions
setenv MRG_MARKETPEN_YN     N          # apply reac. controls market penetration
setenv SPC_INPUT            mole       # [mass|mole] controls speciation merge
setenv SMK_ASCIIELEV_FMT    UAM        # [UAM] format of ASCII elevated file

# For Smkreport
setenv QA_TYPE              state      # run state reports at each proc stage

# Multiple-program controls
setenv PROMPTFLAG           N     # Y (never set to Y for batch processing)
setenv REPORT_DEFAULTS      N     # Y reports default profile application
setenv SMK_O3SEASON_YN      N     # Y uses O3-season emissions instead of annual
setenv SMK_PING_METHOD      0     # 1 processes and outputs for PinG, 0 no PING used
setenv WEST_HSPHERE         Y     # Y converts stack coords to wst hemisphere
setenv SMK_USE_ACTVNAMS     ' '   # Y force use of activities list
setenv SMK_USE_SIPOLS       ' '   # Y force use if pollutants list
setenv SMK_EMLAYS           12    # number of emissions layers
setenv SMK_MAXWARNING       100   # maximum number of warnings in log file
setenv SMK_MAXERROR         100   # maximum number of errors in log file
setenv OUTZONE              0     # output time zone of emissions, GMT = 0
setenv SMK_DEFAULT_TZONE    5     # time zone to fix in missing COSTCY file
setenv SMK_TMPPATH          $SMKDAT/tmp # directory for temporary files
setenv AUTO_DELETE_LOG      N     # Y automatically deletes logs without asking
setenv IOAPI_GRIDNAME_1     T2_36_C  # Grid name defined in the GRIDDESC file

# Invoke Assigns file
source $ASSIGNS_FILE

##############################################################################

## File name override
#
setenv GCNTL   $PPROJ
if ( $?PPMAT ) then
   setenv PCMAT01 $PPMAT
endif

## Run control and growth
#
source cntl_run.csh

## Setup for future-year processing
#
setenv SMK_FUTURE_YN Y
source $ASSIGNS_FILE

## Run programs
#
source smk_run.csh

## Run QA
#
source qa_run.csh

#
## Ending of script
#
exit( $status )

