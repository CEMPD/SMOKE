#!/bin/csh -f
#BSUB -c 0:20

# Version @(#)$Id$
# Path    $Source$
# Date    $Date$

# This script sets up needed environment variables for merging
# emissions in SMOKE, and calls the script that runs the SMOKE programs. 
#
# Script created by : M. Houyoux, CEP Environmental Modeling Center 
#
#*********************************************************************

# set Assigns file name
setenv ASSIGNS_FILE $SMKROOT/assigns/ASSIGNS.nctox.cmaq.cb4p25_wtox.us36-nc

# set source category
setenv MRG_SOURCE ABMP        # source category to merge
setenv MRG_CTLMAT_MULT ' '    # [A|P|AP] for merging with multiplier controls
setenv MRG_CTLMAT_ADD  ' '    # [A|P|AP] for merging with additive controls
setenv MRG_CTLMAT_REAC ' '    # [A|M|P|AMP] for merging with reactivity controls

## For merging from matrices
setenv RUN_SMKMERGE  N        #  run merge program
#      NOTE: in sample script, not run b/c nonroad is separate category

## For merging from previously generated gridded Smkmerge outputs
setenv RUN_MRGGRID   Y        #  run pre-gridded merge program
#      MRGFILES           #  see script settings, below
#      NOTE: in nctox script, mrggrid used to create merged model-ready CMAQ files

## For converting to UAM binary format from either
setenv RUN_SMK2EMIS  N        #  run conversion of 2-d to UAM binary

## Program-specific controls...

# For Smkmerge
setenv MRG_TEMPORAL_YN      Y          # Y merges with hourly emissions
setenv MRG_SPCMAT_YN        N          # Y merges with speciation matrix
setenv MRG_LAYERS_YN        N          # Y merges with layer fractions
setenv MRG_GRDOUT_YN        Y          # Y outputs gridded file
setenv MRG_REPSTA_YN        N          # Y outputs state totals
setenv MRG_REPCNY_YN        N          # Y outputs county totals
setenv SMK_ASCIIELEV_YN     N          # Y outputs ASCII elevated file
setenv MRG_GRDOUT_UNIT      tons/hr    # units for gridded output file
setenv MRG_TOTOUT_UNIT      tons/day   # units for state and/or county totals
setenv MRG_REPORT_TIME      230000     # hour in OUTZONE for reporting emissions
setenv MRG_MARKETPEN_YN     N          # apply reac. controls market penetration
#     EXPLICIT_PLUME_YN           # see multiple-program controls
#     SMK_EMLAYS                  # see multiple-program controls
#     SMK_O3SEASON_YN             # see multiple-program controls
#     SMK_PING_METHOD             # see multiple-program controls, below

# Multiple-program controls
setenv DAY_SPECIFIC_YN      N     # Y imports and uses day-specific inventory
setenv EXPLICIT_PLUME_YN    N     # Y for special wildfire processing for UAM/REMSAD/CAMx
setenv HOUR_SPECIFIC_YN     Y     # Y imports and uses hour-specific inventory
setenv REPORT_DEFAULTS      N     # Y reports default profile application
setenv SMK_EMLAYS           12    # number of emissions layers
setenv SMK_DEFAULT_TZONE    5     # time zone to fix in missing COSTCY file
setenv SMK_O3SEASON_YN      N     # Y uses O3-season emissions instead of annual
setenv SMK_MAXWARNING       100   # maximum number of warnings in log file
setenv SMK_MAXERROR         100   # maximum number of errors in log file
setenv SMK_PING_METHOD      1     # 1 outputs for PinG (using Elevpoint outputs), 0 no PING

# Script settings
setenv MRGFILES   "AGTS_L NGTS_L BGTS_L MGTS_L PGTS3D_L"  # Logical files for Mrggrid
setenv SRCABBR            abmp    # abbreviation for naming log files
setenv PROMPTFLAG         N       # Y (never set to Y for batch processing)
setenv AUTO_DELETE        Y       # Y deletes SMOKE I/O API output files (recommended)
setenv AUTO_DELETE_LOG    Y       # Y automatically deletes logs without asking
setenv DEBUGMODE          N       # Y changes script to use debugger
setenv DEBUG_EXE          dbx     # Sets the debugger to use when DEBUGMODE = Y

# Override settings
# setenv SPC_OVERRIDE  cmaq.cb4p25  # Chemical mechanism override
# setenv INVTABLE_OVERRIDE          # Inventory table override
setenv A_SPC_OVERRIDE cmaq.cb4p25
setenv N_SPC_OVERRIDE cmaq.cb4p25_wtox
setenv B_SPC_OVERRIDE cmaq.cb4p25
setenv M_SPC_OVERRIDE cmaq.cb4p25_wtox
setenv P_SPC_OVERRIDE cmaq.cb4p25

##############################################################################

## Loop through days to run Smkmerge, Mrggrid and Smk2emis
#
setenv RUN_PART2 Y
setenv RUN_PART4 Y
set cnt = 0
set g_stdate_sav = $G_STDATE
while ( $cnt < $EPI_NDAY )

   @ cnt = $cnt + $NDAYS
   source $ASSIGNS_FILE   # Invoke Assigns file to set new dates
   source smk_run.csh     # Run programs
   setenv G_STDATE_ADVANCE $NDAYS

end
setenv RUN_PART2 N
setenv RUN_PART4 N
unsetenv G_STDATE_ADVANCE

#
## Ending of script
#
exit( 0 )

