#!/bin/csh -f

# Version @(#)$Id$
# Path    $Source$
# Date    $Path$

# This script sets up needed environment variables for creating future-year
# area-source emissions and calls the script that runs the SMOKE programs. 
#
# Script created by : M. Houyoux, CEP Environmental Modeling Center 
#                     February 2003
#*********************************************************************

# set Assigns file names
setenv ASSIGNS_FILE $SMKROOT/assigns/ASSIGNS.nctox.cmaq.cb4p25_wtox.us36-nc

# set future year
setenv FYEAR 2018              # year of future case

# set source category
setenv SMK_SOURCE  A           # source category to process
setenv MRG_SOURCE  A           # source category to merge
setenv SRCABBR     ar.$FYEAR   # abbreviation for naming log files

# time independent programs
setenv RUN_CNTLMAT   Y         # Y runs control matrix program
setenv RUN_GRWINVEN  Y         # Y runs control matrix program

# time-dependent programs
setenv RUN_TEMPORAL  Y         # Y runs temporal allocation program
setenv RUN_SMKMERGE  Y         # Y runs merge program
setenv RUN_SMK2EMIS  N         # run conversion of 2-d to UAM binary

# quality assurance
setenv RUN_SMKREPORT Y         # Y runs reporting for state reports

# Program-specific controls...

# For Cntlmat
#      SMK_O3SEASON_YN (see below)# Y uses seas emis in assessing cutoff, etc.
#      REPORT_DEFAULTS            # See multi-program controls

# For Grwinven
setenv SMK_NUM_CTLMAT       1     # number of control/projection matrices
setenv SMK_GRWSMKOUT_YN     Y     # Y outputs a SMOKE-formatted inventory
setenv SMK_GRWIDAOUT_YN     Y     # Y outputs an IDA-formatted inventory

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

# Script settings
setenv SRCABBR            ar      # abbreviation for naming log files
setenv QA_TYPE            all     # [none, all, part1, part2, or custom]
setenv PROMPTFLAG         N       # Y (never set to Y for batch processing)
setenv AUTO_DELETE        Y       # Y deletes SMOKE I/O API output files (recommended)
setenv AUTO_DELETE_LOG    Y       # Y automatically deletes logs without asking
setenv DEBUGMODE          N       # Y changes script to use debugger
setenv DEBUG_EXE          ldb     # Sets the debugger to use when DEBUGMODE = Y

# Override settings (comment out if not used)
setenv SPC_OVERRIDE  cmaq.cb4p25  # Chemical mechanism override
# setenv CNTLCASE                   # Control case
# setenv INVTABLE_OVERRIDE          # Inventory table override

##############################################################################

## NOTE: GCNTL file for area sources must be in the $ARDAT directory for the
#        base case simulation and should have the names, as follows:
#            Projection only       :  gcntl.$YEAR_$FYEAR.txt,
#            Control only          :  gcntl.$CNTLCASE.txt,
#            Projection and control:  gcntl.$YEAR_$FYEAR_$CNTLCASE.txt,
#        where YEAR is set in the Assigns file and FYEAR and CNTLCASE
#        are set in this script.

setenv RUN_PART1 Y
source $ASSIGNS_FILE   # Invoke Assigns file

## Set projection and control matrices to use in Grwinven
setenv ACMAT01 $APMAT
# setenv ACMAT02 
# setenv ACMAT03 
# setenv ACMAT04 

## Run Cntlmat and Grwinven
#
source cntl_run.csh    # Run programs

## Set up for future-year and/or control processing
#
setenv SMK_FUTURE_YN Y
setenv RUN_PART1 N     # Temporarily end part 1
source $ASSIGNS_FILE

setenv RUN_PART1 Y     # Restart part 1
source qa_run.csh      # Run QA for part 1
setenv RUN_PART1 N     # End part 1

## Run Temporal, Smkmerge, and Smk2emis programs
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

## Ending of script
#
exit( $status )

