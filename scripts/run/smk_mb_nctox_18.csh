#!/bin/csh -f

# Version @(#)$Id$
# Path    $Source$
# Date    $Date$

# This script sets up needed environment variables for creating future-year
# mobile-source emissions and calls the script that runs the SMOKE programs. 
#
# Script created by : M. Houyoux, CEP Environmental Modeling Center 
#                     February 2003
#*********************************************************************

# set Assigns file names
setenv ASSIGNS_FILE $SMKROOT/assigns/ASSIGNS.nctox.cmaq.cb4p25_wtox.us36-nc

# set future year
setenv FYEAR 2018              # year of future case

# set source category
setenv SMK_SOURCE  M           # source category to process
setenv MRG_SOURCE  M           # source category to merge

## time independent programs
setenv RUN_CNTLMAT   Y         # Y runs control matrix program
setenv RUN_GRWINVEN  Y         # Y runs control matrix program
setenv RUN_MBSETUP   Y        #  run speed/temperature setup program        

## episode dependent programs
setenv RUN_PREMOBL   Y        #  Y runs temperature preprocessing program
setenv RUN_EMISFAC   Y        #  Y runs emission factors program

## time-dependent programs
setenv RUN_TEMPORAL  Y         # Y runs temporal allocation program
setenv RUN_SMKMERGE  Y         # Y runs merge program
setenv RUN_SMK2EMIS  N         # run conversion of 2-d to UAM binary

## quality assurance
setenv RUN_SMKREPORT Y         # Y runs reporting for state reports

# Program-specific controls...

# For Cntlmat
#      SMK_O3SEASON_YN (see below)# Y uses seas emis in assessing cutoff, etc.
#      REPORT_DEFAULTS            # See multi-program controls

# For Grwinven
setenv SMK_NUM_CTLMAT       1     # number of control/projection matrices
setenv SMK_GRWSMKOUT_YN     Y     # Y outputs a SMOKE-formatted inventory
setenv SMK_GRWIDAOUT_YN     Y     # Y outputs an IDA-formatted inventory

# For Mbsetup 
#     USE_SPEED_PROFILES      # see multiple-program controls, below
#     Temporary file paths    # set automatically

# For Prediur
setenv SMK_MINTEMP           0.   # MOBILE6 minimum allowed daily temperature [deg F]
setenv SMK_MAXTEMP         120.   # MOBILE6 Maximum allowed daily temperature [deg F]
setenv UNIFORM_STIME        -1    # -1 or HHMMSS for uniform start hour for "days"
#     OUTZONE                 # see multiple-program controls, below
#     Episode settings        # in Assigns file
#     Temporary file paths    # set automatically

# For Emisfac
setenv MB_HC_TYPE           TOG   # hydrocarbon type from MOBILE6
setenv REPLACE_TEMPERATURES   Y   # Y: use gridded hourly tmprs; N: use tmprs in M6 files
#     EF_YEAR                 # SMK_FUTURE_YN=N: set to $YEAR from Assigns
#                             # SMK_FUTURE_YN=Y: set to $FYEAR from scripts
#     GROUP_TYPE              # set to [daily|weekly|monthly|episode], see below
#     USE_SPEED_PROFILES      # see multiple-program controls, below
#     Temporary file paths    # set automatically

# For Temporal
setenv RENORM_TPROF         Y     # Y renormalizes temporal profiles
setenv UNIFORM_TPROF_YN     N     # Y makes all temporal profiles uniform
setenv ZONE4WM              Y     # Y uses time zones for start of day & month
#     OUTZONE                 # see multiple-program controls, below
#     REPORT_DEFAULTS         # see multiple-program controls, below
#     SMK_O3SEASON_YN         # see multiple-program controls, below

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
setenv USE_SPEED_PROFILES   N     # Y uses speed profiles instead of inventory speeds

# Script settings
setenv SRCABBR            mb.$FYEAR # abbreviation for naming log files
setenv QA_TYPE            all     # [none, all, part1, part2, or custom]
setenv PROMPTFLAG         N       # Y (never set to Y for batch processing)
setenv AUTO_DELETE        Y       # Y deletes SMOKE I/O API output files (recommended)
setenv AUTO_DELETE_LOG    Y       # Y automatically deletes logs without asking
setenv DEBUGMODE          N       # Y changes script to use debugger
setenv DEBUG_EXE          dbx     # Sets the debugger to use when DEBUGMODE = Y

# Override settings (comment out if not used)
setenv SPC_OVERRIDE       cmaq.cb4p25_wtox.m  # Chemical mechanism override
setenv YEAR_OVERRIDE      1999                # Overrides YEAR (base) in Assigns file
setenv INVTABLE_OVERRIDE  invtable_onroad.cb4.120202.txt   # Inventory table override
# setenv CNTLCASE               # Control case

##############################################################################

## NOTE: GCNTL file for mobile sources must be in the $MBDAT directory for the
#        base case simulation and should have the names, as follows:
#          Projection only     :  gcntl.$YEAR_$FYEAR.txt,
#          Control only        :  gcntl.$CNTLCASE.txt,
#          Projection & control:  gcntl.$YEAR_$FYEAR_$CNTLCASE.txt,
#        where YEAR is set in the Assigns file and FYEAR and CNTLCASE
#        are set in this script.
#
#        In addition, the MOBILE6 inputs for the case you want to run
#        should all have file names with a .in extension and should
#        be placed in a *directory* named as follows:
#          No control name    : $INVDIR/mobile/m6_$EF_YEAR/
#          With control name  : $INVDIR/mobile/m6_$EFYEAR_$CNTLCASE/
#        where EF_YEAR is either YEAR if SMK_FUTURE_YN = N and is
#        FYEAR if SMK_FUTURE_YN = Y.  The CNTLCASE variable will only
#        be used in the directory name if it is defined.

setenv RUN_PART1 Y
source $ASSIGNS_FILE   # Invoke Assigns file

## Set projection and control matrices to use in Grwinven
#      NOTE: Grwinven setting SMK_NUM_CTLMAT > 1 to use MCMAT02-04
setenv MCMAT01 $MPMAT
# setenv MCMAT02 
# setenv MCMAT03 
# setenv MCMAT04 

## Run Cntlmat and Grwinven
#
if ( $?CNTLCASE ) then
    setenv SRCABBR $SRCABBR.$CNTLCASE 
endif
source cntl_run.csh    # Run programs

## Set up for future-year and/or control processing
#
setenv SMK_FUTURE_YN Y
setenv SMK_CONTROL_YN N
setenv RUN_PART1 N     # Temporarily end part 1 to prevent file deletion
source $ASSIGNS_FILE

## Run QA for new grown inventory
#
setenv RUN_PART1 Y     # Restart part 1
source qa_run.csh      # Run QA for part 1

## Run Mbsetup and Premobl
#
source smk_run.csh     # Run remaining programs for part 1
setenv RUN_PART1 N     # End part 1

## Run emisfac for cases needed by temperature choices in MVREF file
#
foreach group ( episode )  # can include ( daily weekly monthly episode )
  setenv GROUP_TYPE $group
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
   setenv G_STDATE_ADVANCE $cnt

end
setenv RUN_PART2 N
setenv RUN_PART4 N
unsetenv G_STDATE_ADVANCE

## Ending of script
#
exit( $status )

