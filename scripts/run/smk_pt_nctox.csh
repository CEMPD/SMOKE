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
setenv ASSIGNS_FILE $SMKROOT/assigns/ASSIGNS.nctox.cmaq.cb4p25_wtox.us36-nc

# set source category
setenv SMK_SOURCE P           # source category to process
setenv MRG_SOURCE P           # source category to merge
setenv MRG_CTLMAT_MULT ' '    # [A|P|AP] for merging with multiplier controls
setenv MRG_CTLMAT_ADD  ' '    # [A|P|AP] for merging with additive controls
setenv MRG_CTLMAT_REAC ' '    # [A|M|P|AMP] for merging with reactivity controls

## time independent programs
setenv RUN_SMKINVEN  Y        #  run inventory import program
setenv RUN_SPCMAT    Y        #  run speciation matrix program
setenv RUN_GRDMAT    Y        #  run gridding matrix program
setenv RUN_CNTLMAT   N        #  run control matrix program

## time-dependent programs
setenv RUN_TEMPORAL  Y        #  run temporal allocation program
setenv RUN_ELEVPOINT Y        #  run elevated/PinG sources selection program
setenv RUN_LAYPOINT  Y        #  run layer fractions program
setenv RUN_SMKMERGE  Y        #  run merge program
setenv RUN_SMK2EMIS  N        #  run conversion of 2-d to UAM binary

## quality assurance
setenv RUN_SMKREPORT N        # Y runs reporting for state reports

## Program-specific controls...

## For Smkinven
setenv FILL_ANN_WSEAS       N # Y fills annual value when only seasonal is provided
setenv HOURLY_TO_DAILY      N # Y reads daily total only from hourly file
setenv HOURLY_TO_PROFILE    N # Y converts hourly data to source-specific profs
setenv IMPORT_AVEINV_YN     Y # Y then import annual/average inventory
setenv RAW_DUP_CHECK        N # Y errors on duplicate records
setenv SMK_BASEYR_OVERRIDE  0 # Enter year of the base year when future-year inven provided
setenv SMKINVEN_FORMULA     "PMC=PM10-PM2_5" # Internal PMC calculation
setenv SMK_SWITCH_EPSXY     N # Y corrects x-y location in file for for EPS format
setenv WEST_HSPHERE         Y # Y converts ALL stack coords to western hemisphere
setenv WKDAY_NORMALIZE      Y # Y normalizes weekly profiles by weekdays
#     DAY_SPECIFIC_YN     # see multiple-program controls, below
#     HOUR_SPECIFIC_YN    # see multiple-program controls, below
#     OUTZONE             # see multiple-program controls, below
#     REPORT_DEFAULTS     # see multiple-program controls, below
#     VELOC_RECALC        # See multiple-program controls, below

## For Spcmat
setenv POLLUTANT_CONVERSION Y     # Y uses ROG to TOG file, for example
#     REPORT_DEFAULTS         # see multiple-program controls, below

## For Cntlmat
setenv REACTIVITY_POL       ' '   # Set to VOC or ROG (only for reactivity controls) 
#     REPORT_DEFAULTS         # see multiple-program controls, below
#     SMK_O3SEASON_YN         # see multiple-program controls

# For Elevpoint
setenv SMK_ELEV_METHOD      1     # 0=Laypoint sets elev srcs; 1=use PELVCONFIG
#     SMK_PING_METHOD         # see multiple-program controls, below

# For Temporal
setenv RENORM_TPROF         Y     # Y renormalizes temporal profiles
setenv UNIFORM_TPROF_YN     N     # Y makes all temporal profiles uniform
setenv ZONE4WM              Y     # Y uses time zones for start of day & month
#     DAY_SPECIFIC_YN         # see multiple-program controls, below
#     HOUR_SPECIFIC_YN        # see multiple-program controls, below
#     OUTZONE                 # see multiple-program controls, below
#     REPORT_DEFAULTS         # see multiple-program controls, below
#     SMK_O3SEASON_YN         # see multiple-program controls, below

# For Laypoint
setenv REP_LAYER_MAX        ' '   # Layer no. for reporting high plume rise
setenv SMK_SPECELEV_YN      N     # Y: Laypoint uses Elevpoint outputs to pick elevated
#     SMK_EMLAYS              # see multiple-program controls, below
#     EXPLICIT_PLUME_YN       # see multiple-program controls, below
#     VELOC_RECALC            # eee multiple-program controls, below

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
#     EXPLICIT_PLUME_YN           # see multiple-program controls
#     SMK_EMLAYS                  # see multiple-program controls
#     SMK_O3SEASON_YN             # see multiple-program controls
#     SMK_PING_METHOD             # see multiple-program controls, below

# Multiple-program controls
setenv DAY_SPECIFIC_YN      N     # Y imports and uses day-specific inventory
setenv EXPLICIT_PLUME_YN    N     # Y for special wildfire processing for UAM/REMSAD/CAMx
setenv HOUR_SPECIFIC_YN     N     # Y imports and uses hour-specific inventory
setenv OUTZONE              0     # output time zone of emissions
setenv REPORT_DEFAULTS      N     # Y reports default profile application
setenv SMK_EMLAYS           12    # number of emissions layers
setenv SMK_DEFAULT_TZONE    5     # time zone to fix in missing COSTCY file
setenv SMK_O3SEASON_YN      N     # Y uses O3-season emissions instead of annual
setenv SMK_MAXWARNING       100   # maximum number of warnings in log file
setenv SMK_MAXERROR         100   # maximum number of errors in log file
setenv SMK_PING_METHOD      1     # 1 outputs for PinG (using Elevpoint outputs), 0 no PING
setenv SMK_SPECELEV_YN      N     # Y uses the indicator for major/minor sources
setenv VELOC_RECALC         N     # Y recalculates velocity from diam and flow

# Script settings
setenv SRCABBR            pt      # abbreviation for naming log files
setenv QA_TYPE            state   # run state reports at each proc stage
setenv PROMPTFLAG         N       # Y (never set to Y for batch processing)
setenv AUTO_DELETE        Y       # Y deletes SMOKE I/O API output files (recommended)
setenv AUTO_DELETE_LOG    Y       # Y automatically deletes logs without asking
setenv DEBUGMODE          N       # Y changes script to use debugger
setenv DEBUG_EXE          dbx     # Sets the debugger to use when DEBUGMODE = Y

##############################################################################

## Set controls to run processing in parts
setenv RUN_PART1 N    # Smkinven, Spcmat, Grdmat, Cntlmat
setenv RUN_PART2 N    # Temporal
setenv RUN_PART3 N    # Elevpoint
setenv RUN_PART4 N    # Laypoint, Smkmerge, Smk2emis

## Run Smkinven, Spcmat, Grdmat, Cntlmat, if needed
#
setenv RUN_PART1 Y
source $ASSIGNS_FILE   # Invoke Assigns file
source smk_run.csh     # Run programs
source qa_run.csh      # Run QA for part 1
setenv RUN_PART1 N

## Loop through days to run Temporal
#
setenv RUN_PART2 Y
set cnt = 0
set g_stdate_sav = $G_STDATE
while ( $cnt < $EPI_NDAY )

   @ cnt = $cnt + $NDAYS
   source $ASSIGNS_FILE   # Invoke Assigns file to set new dates
   source smk_run.csh     # Run programs
   source qa_run.csh      # Run QA for part 2

   setenv G_STDATE_ADVANCE $NDAYS

end
setenv RUN_PART2 N
unsetenv G_STDATE_ADVANCE

## Run Elevpoint
#
setenv RUN_PART3 Y
setenv G_STDATE  $g_stdate_sav
setenv ESDATE `$IOAPIDIR/datshift $G_STDATE 0`
source $ASSIGNS_FILE   # Invoke Assigns file to set new dates
source smk_run.csh     # Run programs
source qa_run.csh      # Run QA for part 3
setenv RUN_PART3 N

## Loop through dats to run Laypoint and Smkmerge
setenv RUN_PART4 Y
set cnt = 0
while ( $cnt < $EPI_NDAY )

   @ cnt = $cnt + $NDAYS
   source $ASSIGNS_FILE   # Invoke Assigns file to set new dates
   source smk_run.csh     # Run programs

   setenv G_STDATE_ADVANCE $NDAYS

end
setenv RUN_PART4 N

#
## Ending of script
#
exit( $outstat )

