#!/bin/csh -f

# Version @(#)$Id$
# Path    $Source$
# Date    $Date$

# This script sets up needed environment variables for running point source
# emissions in SMOKE, and calls the script that runs the SMOKE programs. 
#
# Script created by : M. Houyoux, MCNC Environmental Modeling Center 
#
#*********************************************************************

# set Assigns file name
setenv ASSIGNS_FILE $SMKROOT/assigns/ASSIGNS.net96.us36

# set source category
setenv SMK_SOURCE P           # source category to process
setenv MRG_SOURCE P           # source category to merge
setenv SRCABBR    pt          # abbreviation for naming log files
setenv MRG_CTLMAT_MULT ' '    # [A|P|AP] for merging with multiplier controls
setenv MRG_CTLMAT_ADD  ' '    # [A|P|AP] for merging with additive controls
setenv MRG_CTLMAT_REAC ' '    # [A|M|P|AMP] for merging with reactivity controls

# time independent programs
setenv RUN_SMKINVEN  Y        #  run inventory import program
setenv RUN_SPCMAT    Y        #  run speciation matrix program
setenv RUN_GRDMAT    Y        #  run gridding matrix program
setenv RUN_CNTLMAT   N        #  run control matrix program
setenv RUN_ELEVPOINT N        #  run elevated/PinG sources selection program

# time-dependent programs

setenv RUN_LAYPOINT  Y        #  run layer fractions program
setenv RUN_TEMPORAL  Y        #  run temporal allocation program
setenv RUN_SMKMERGE  Y        #  run merge program

# quality assurance
setenv RUN_SMKREPORT Y         # Y runs reporting for state reports

# Program-specific controls...

# For Smkinven
setenv HOURLY_TO_DAILY      N # Y reads daily total only from hourly file
setenv HOURLY_TO_PROFILE    N # Y converts hourly data to source-specific profs
setenv IMPORT_AVEINV_YN     Y # Y then import annual/average inventory
setenv RAW_DUP_CHECK        N # Y errors on duplicate records
setenv RAW_SRC_CHECK        N # Y errors if missing pollutants or activities
setenv VELOC_RECALC         N # Y recalculates velocity from diam and flow
setenv WKDAY_NORMALIZE      Y # Y normalizes weekly profiles by weekdays
setenv SMKINVEN_FORMULA     "PMC=PM10-PM2_5" # Internal PMC calculation

# For Spcmat
setenv POLLUTANT_CONVERSION Y     # Y uses ROG to TOG file, for example
setenv SMK_GSREF_FIXED      N     # Y uses fixed-format speciation x-ref
setenv SPEC_OUTPUT          ALL   # [ALL|MASS|MOLE] controls program outputs

# For Cntlmat
setenv TMP_CTL_PATH $EDSS_ROOT/data/tmp # directory for temporary program files
setenv CONTROL_REPORT       Y     # Y outputs a report of factors by source

# For Elevpoint
setenv SMK_CUTOFF_HT        0.    # Elev source cutoff ht [m] or 0. if not used
setenv SMK_ENG2METRIC_YN    N     # Y converts stack parameters to metric units

# For Temporal
setenv RENORM_TPROF         Y     # Y renormalizes temporal profiles
setenv UNIFORM_TPROF_YN     N     # Y makes all temporal profiles uniform
setenv ZONE4WM              Y     # Y uses time zones for start of day & month

# For Laypoint
setenv REP_LAYER_MAX        ' '   # Layer no. for reporting high plume rise

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
setenv HOUR_SPECIFIC_YN     N     # Y imports and uses hour-specific inventory
setenv DAY_SPECIFIC_YN      N     # Y imports and uses day-specific inventory
setenv SMK_O3SEASON_YN      N     # Y uses O3-season emissions instead of annual
setenv SMK_PING_YN          N     # Y processes and outputs for PinG 
setenv SMK_SPECELEV_YN      N     # Y uses the indicator for major/minor sources
setenv WEST_HSPHERE         Y     # Y converts stack coords to wst hemisphere
setenv SMK_USE_ACTVNAMS     ' '   # Y force use of activities list
setenv SMK_USE_SIPOLS       ' '   # Y force use if pollutants list
setenv SMK_EMLAYS           12    # number of emissions layers
setenv SMK_MAXWARNING       100   # maximum number of warnings in log file
setenv SMK_MAXERROR         100   # maximum number of errors in log file
setenv OUTZONE              0     # output time zone of emissions
setenv SMK_DEFAULT_TZONE    5     # time zone to fix in missing COSTCY file
setenv AUTO_DELETE_LOG      Y     # Y automatically deletes logs without asking

# Invoke Assigns file
source $ASSIGNS_FILE

##############################################################################

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

