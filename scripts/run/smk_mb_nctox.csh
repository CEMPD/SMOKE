#!/bin/csh -f

# Version @(#)$Id$
# Path    $Source$
# Date    $Date$

# This script sets up needed environment variables for running on-road mobile
# source emissions in SMOKE, and calls the script that runs the SMOKE
# programs. 
#
# Script created by : M. Houyoux, CEP Environmental Modeling Center 
#
#*********************************************************************

# set Assigns file name
setenv ASSIGNS_FILE $SMKROOT/assigns/ASSIGNS.nctox.cmaq.cb4p25_wtox.mb.us36-nc

# set source category
setenv SMK_SOURCE M           # source category to process
setenv MRG_SOURCE M           # source category to merge
setenv MRG_CTLMAT_MULT ' '    # [A|P|AP] for merging with multiplier controls
setenv MRG_CTLMAT_ADD  ' '    # [A|P|AP] for merging with additive controls
setenv MRG_CTLMAT_REAC ' '    # [A|M|P|AMP] for merging with reactivity controls

## time independent programs
setenv RUN_SMKINVEN  Y        #  run inventory import program
setenv RUN_SPCMAT    N        #  run speciation matrix program
setenv RUN_GRDMAT    N        #  run gridding matrix program
setenv RUN_CNTLMAT   N        #  run control matrix program
setenv RUN_MBSETUP   N        #  run speed/temperature setup program        

## episode dependent programs
setenv RUN_PREMOBL   N        #  Y runs temperature preprocessing program
setenv RUN_EMISFAC   N        #  Y runs emission factors program

## time-dependent programs
setenv RUN_TEMPORAL  N        #  run temporal allocation program
setenv RUN_SMKMERGE  N        #  run merge program
setenv RUN_SMK2EMIS  N        #  run conversion of 2-d to UAM binary

## quality assurance
setenv RUN_SMKREPORT N        # Y runs reporting for state reports

## Program-specific controls...

## For Smkinven
setenv FILL_ANN_WSEAS       N  # Y fills annual value when only seasonal is provided
setenv IMPORT_VMTMIX_YN     N  # Y imports VMT mix - needed for EMS-95 inputs
setenv RAW_DUP_CHECK        N  # Y errors on duplicate records
setenv SMK_BASEYR_OVERRIDE  0  # Enter year of the base year when future-year inven provided
setenv SMK_NHAPEXLCUDE_YN   Y  # Y uses NonHAP exclusions file
setenv SMKINVEN_FORMULA    " " # Internal PMC calculation
setenv SMK_VMTMIX_FIXFMT    N  # Y uses fixed-format VMT Mix file (EMS-95 input format only)
setenv WEST_HSPHERE         Y  # Y converts ALL stack coords to western hemisphere
setenv WKDAY_NORMALIZE      Y  # Y normalizes weekly profiles by weekdays
#     OUTZONE              # see multiple-program controls, below
#     REPORT_DEFAULTS      # see multiple-program controls, below

## For Grdmat
setenv SMK_DEFAULT_SRGID    8     # default surrogate code (8=popl'n)

## For Spcmat
setenv POLLUTANT_CONVERSION Y     # Y uses ROG to TOG file, for example
#     REPORT_DEFAULTS         # see multiple-program controls, below

## For Cntlmat
setenv REACTIVITY_POL       ' '   # Set to VOC or ROG (only for reactivity controls) 
#     REPORT_DEFAULTS         # see multiple-program controls, below
#     SMK_O3SEASON_YN         # see multiple-program controls

# For Mbsetup (none)
#     Temporary file paths    # set automatically

# For Prediur
setenv SMK_MINTEMP           0.   # MOBILE6 minimum allowed daily temperature [deg F]
setenv SMK_MAXTEMP         120.   # MOBILE6 Maximum allowed daily temperature [deg F]
setenv UNIFORM_STIME        -1    # -1 or HHMMSS for uniform start hour for "days"
#     OUTZONE                 # see multiple-program controls, below
#     Episode settings        # in Assigns file
#     Temporary file paths    # set automatically

# For Emisfac
setenv MB_HC_TYPE           VOC   # hydrocarbon type from MOBILE6
setenv REPLACE_TEMPERATURES   Y   # Y: use gridded hourly tmprs; N: use tmprs in M6 files
#     EF_YEAR                 # SMK_FUTURE_YN=N: set to $YEAR from Assigns
#                             # SMK_FUTURE_YN=Y: set to $FYEAR from scripts
#     GROUP_TYPE              # set to [daily|weekly|monthly|episode], see below
#     Temporary file paths    # set automatically

# For Temporal
setenv RENORM_TPROF         Y     # Y renormalizes temporal profiles
setenv UNIFORM_TPROF_YN     N     # Y makes all temporal profiles uniform
setenv ZONE4WM              Y     # Y uses time zones for start of day & month
#     OUTZONE                 # see multiple-program controls, below
#     REPORT_DEFAULTS         # see multiple-program controls, below
#     SMK_O3SEASON_YN         # see multiple-program controls, below
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
#     SMK_O3SEASON_YN             # see multiple-program controls

# Multiple-program controls
setenv OUTZONE              0     # output time zone of emissions
setenv REPORT_DEFAULTS      N     # Y reports default profile application
setenv SMK_DEFAULT_TZONE    5     # time zone to fix in missing COSTCY file
setenv SMK_O3SEASON_YN      N     # Y uses O3-season emissions instead of annual
setenv SMK_MAXWARNING       100   # maximum number of warnings in log file
setenv SMK_MAXERROR         100   # maximum number of errors in log file

# Script settings
setenv SRCABBR            mb      # abbreviation for naming log files
setenv QA_TYPE            all     # [none, all, part1-part4, or custom]
setenv PROMPTFLAG         N       # Y (never set to Y for batch processing)
setenv AUTO_DELETE        Y       # Y deletes SMOKE I/O API output files (recommended)
setenv AUTO_DELETE_LOG    Y       # Y automatically deletes logs without asking
setenv DEBUGMODE          Y       # Y changes script to use debugger
setenv DEBUG_EXE          dbx     # Sets the debugger to use when DEBUGMODE = Y

##############################################################################

## Run Smkinven, Spcmat, Grdmat, Cntlmat, Mbsetup, and Premobl if needed
#
setenv RUN_PART1 Y
source $ASSIGNS_FILE   # Invoke Assigns file
source smk_run.csh     # Run programs
source qa_run.csh      # Run QA for part 1
setenv RUN_PART1 N

## Run emisfac for cases needed by temperature choices in MVREF file
#
foreach GROUP_TYPE  ( episode )  # can include ( daily weekly monthly episode )
  source emisfac_run.csh     # Run programs
endif

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
   setenv G_STDATE_ADVANCE $NDAYS

end
setenv RUN_PART2 N
setenv RUN_PART4 N
unsetenv G_STDATE_ADVANCE

#
## Ending of script
#
exit( 0 )

