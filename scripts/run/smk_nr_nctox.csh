#!/bin/csh -f
#BSUB 0:30

# Version @(#)$Id$
# Path    $Source$
# Date    $Date$

# This script sets up needed environment variables for running area source
# emissions in SMOKE, and calls the script that runs the SMOKE programs. 
#
# Script created by : M. Houyoux, CEP Environmental Modeling Center 
#
#*********************************************************************

# set Assigns file name
setenv ASSIGNS_FILE $SMKROOT/assigns/ASSIGNS.nctox.cmaq.cb4p25_wtox.us36-nc

# set source category
setenv SMK_SOURCE A           # source category to process
setenv MRG_SOURCE A           # source category to merge
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
setenv RUN_SMKMERGE  Y        #  run merge program
setenv RUN_SMK2EMIS  N        #  run conversion of 2-d to UAM binary

## quality assurance
setenv RUN_SMKREPORT Y        # Y runs reporting for state reports

## Program-specific controls...

## For Smkinven
setenv FILL_ANNUAL          N  # Y fills annual value when only average day is provided
setenv IMPORT_GRDIOAPI_YN   N  # Y imported gridded I/O API inventory
setenv RAW_DUP_CHECK        N  # Y errors on duplicate records
setenv SMK_ARTOPNT_YN       Y  # Y uses area-to-point conversions
setenv SMK_BASEYR_OVERRIDE  0  # Enter year of the base year when future-year inven provided
setenv SMK_NHAPEXCLUDE_YN   Y  # Y uses NonHAP exclusions file
setenv SMKINVEN_FORMULA     "PMC=PM10-PM2_5" # Internal PMC calculation
setenv WEST_HSPHERE         Y  # Y converts ALL stack coords to western hemisphere
setenv WKDAY_NORMALIZE      N  # Y normalizes weekly profiles by weekdays
#     REPORT_DEFAULTS      # see multiple-program controls, below

## For Grdmat
setenv SMK_DEFAULT_SRGID    8     # default surrogate code (8=popl'n)

## For Spcmat
setenv POLLUTANT_CONVERSION Y     # Y uses ROG to TOG file, for example
#     REPORT_DEFAULTS         # see multiple-program controls, below

## For Cntlmat
setenv REACTIVITY_POL       ' '   # Set to VOC or ROG (only for reactivity controls) 
#     REPORT_DEFAULTS         # see multiple-program controls, below
#     SMK_AVEDAY_YN           # see multiple-program controls

# For Temporal
setenv RENORM_TPROF         Y     # Y renormalizes temporal profiles
setenv UNIFORM_TPROF_YN     N     # Y makes all temporal profiles uniform
setenv ZONE4WM              Y     # Y uses time zones for start of day & month
#     OUTZONE                 # see multiple-program controls, below
#     REPORT_DEFAULTS         # see multiple-program controls, below
#     SMK_AVEDAY_YN           # see multiple-program controls, below
#     Date/time settings      # in Assigns file

# For Smkmerge
setenv MRG_TEMPORAL_YN      Y          # Y merges with hourly emissions
setenv MRG_SPCMAT_YN        Y          # Y merges with speciation matrix
setenv MRG_GRDOUT_YN        Y          # Y outputs gridded file
setenv MRG_REPSTA_YN        N          # Y outputs state totals
setenv MRG_REPCNY_YN        N          # Y outputs county totals
setenv MRG_GRDOUT_UNIT      moles/s    # units for gridded output file
setenv MRG_TOTOUT_UNIT      moles/day  # units for state and/or county totals
setenv MRG_REPORT_TIME      230000     # hour in OUTZONE for reporting emissions
setenv MRG_MARKETPEN_YN     N          # apply reac. controls market penetration
#     SMK_AVEDAY_YN               # see multiple-program controls

# For Smk2emis
setenv SMK2EMIS_VMAP_YN     N     # Y uses name remapping file

# Multiple-program controls
setenv OUTZONE              0     # output time zone of emissions
setenv REPORT_DEFAULTS      N     # Y reports default profile application
setenv SMK_DEFAULT_TZONE    5     # time zone to fix in missing COSTCY file
setenv SMK_AVEDAY_YN        N     # Y uses average day emissions instead of annual
setenv SMK_MAXWARNING       100   # maximum number of warnings in log file
setenv SMK_MAXERROR         100   # maximum number of errors in log file

# Script settings
setenv SRCABBR            nr      # abbreviation for naming log files
setenv QA_TYPE            all     # [none, all, part1, part2, or custom]
setenv NONROAD            Y       # Y resets area files to nonroad files
setenv PROMPTFLAG         N       # Y (never set to Y for batch processing)
setenv AUTO_DELETE        Y       # Y deletes SMOKE I/O API output files (recommended)
setenv AUTO_DELETE_LOG    Y       # Y automatically deletes logs without asking
setenv DEBUGMODE          N       # Y changes script to use debugger
setenv DEBUG_EXE          dbx     # Sets the debugger to use when DEBUGMODE = Y

# Override settings
# setenv SPC_OVERRIDE  cmaq.cb4p25  # Chemical mechanism override
# setenv YEAR_OVERRIDE              # Overrides YEAR (base) in Assigns file
# setenv INVTABLE_OVERRIDE          # Inventory table override

##############################################################################

## Run Smkinven, Spcmat, Grdmat, Cntlmat, if needed
#
setenv RUN_PART1 Y
source $ASSIGNS_FILE   # Invoke Assigns file
source smk_run.csh     # Run programs
source qa_run.csh      # Run QA for part 1
setenv RUN_PART1 N

## Loop through days to run Temporal, Smkmerge, and Smk2emis
#
setenv RUN_PART2 Y
setenv RUN_PART4 Y
set cnt = 0
set g_stdate_sav = $G_STDATE
while ( $cnt < $EPI_NDAY )

   @ cnt = $cnt + $NDAYS
   source $ASSIGNS_FILE   # Invoke Assigns file to set new dates
   source smk_run.csh     # Run programs
   source qa_run.csh      # Run QA for part 2
   setenv G_STDATE_ADVANCE $cnt

end
setenv RUN_PART2 N
setenv RUN_PART4 N
unsetenv G_STDATE_ADVANCE

#
## Ending of script
#
exit( 0 )

