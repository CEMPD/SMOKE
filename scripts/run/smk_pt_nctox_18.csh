#!/bin/csh -f

# Version @(#)$Id$
# Path    $Source$
# Date    $Path$

# This script sets up needed environment variables for creating future-year
# point-source emissions and calls the script that runs the SMOKE programs. 
#
# Script created by : M. Houyoux, CEP Environmental Modeling Center 
#                     February 2003
#*********************************************************************

# set Assigns file names
setenv ASSIGNS_FILE $SMKROOT/assigns/ASSIGNS.nctox.cmaq.cb4p25_wtox.us36-nc

# set future year
setenv FYEAR 2018              # year of future case

# set source category
setenv SMK_SOURCE  P           # source category to process
setenv MRG_SOURCE  P           # source category to merge

# time independent programs
setenv RUN_CNTLMAT   Y         # Y runs control matrix program
setenv RUN_GRWINVEN  Y         # Y runs control application program

# time-dependent programs
setenv RUN_TEMPORAL  Y         # Y runs temporal allocation program
setenv RUN_ELEVPOINT Y        #  run elevated/PinG sources selection program
setenv RUN_SMKMERGE  Y         # Y runs merge program
setenv RUN_SMK2EMIS  N         # run conversion of 2-d to UAM binary

# quality assurance
setenv RUN_SMKREPORT Y         # Y runs reporting for state reports

# Program-specific controls...

# For Cntlmat
#      SMK_AVEDAY_YN (see below)  # Y uses seas emis in assessing cutoff, etc.
#      REPORT_DEFAULTS            # See multi-program controls

# For Grwinven
setenv SMK_NUM_CTLMAT       1     # number of control/projection matrices
setenv SMK_GRWSMKOUT_YN     Y     # Y outputs a SMOKE-formatted inventory
setenv SMK_GRWIDAOUT_YN     Y     # Y outputs an IDA-formatted inventory

# For Elevpoint
setenv SMK_ELEV_METHOD      1     # 0=Laypoint sets elev srcs; 1=use PELVCONFIG
setenv UNIFORM_STIME        -1    # -1 or HHMMSS for uniform start hour for daily emissions days
#     SMK_PING_METHOD         # see multiple-program controls, below

# For Temporal
setenv RENORM_TPROF         Y     # Y renormalizes temporal profiles
setenv UNIFORM_TPROF_YN     N     # Y makes all temporal profiles uniform
setenv ZONE4WM              Y     # Y uses time zones for start of day & month
#     DAY_SPECIFIC_YN         # see multiple-program controls, below
#     HOUR_SPECIFIC_YN        # see multiple-program controls, below
#     OUTZONE                 # see multiple-program controls, below
#     REPORT_DEFAULTS         # see multiple-program controls, below
#     SMK_AVEDAY_YN           # see multiple-program controls, below
#     Date/time settings      # in Assigns file

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
#     SMK_AVEDAY_YN               # see multiple-program controls
#     SMK_PING_METHOD             # see multiple-program controls, below

# For Smk2emis
setenv SMK2EMIS_VMAP_YN     N     # Y uses name remapping file

# Multiple-program controls
setenv DAY_SPECIFIC_YN      N     # Y imports and uses day-specific inventory
setenv EXPLICIT_PLUME_YN    N     # Y for special wildfire processing for UAM/REMSAD/CAMx
setenv HOUR_SPECIFIC_YN     N     # Y imports and uses hour-specific inventory
setenv OUTZONE              0     # output time zone of emissions
setenv REPORT_DEFAULTS      N     # Y reports default profile application
setenv SMK_EMLAYS           12    # number of emissions layers
setenv SMK_DEFAULT_TZONE    5     # time zone to fix in missing COSTCY file
setenv SMK_AVEDAY_YN        N     # Y uses average day emissions instead of annual
setenv SMK_MAXWARNING       100   # maximum number of warnings in log file
setenv SMK_MAXERROR         100   # maximum number of errors in log file
setenv SMK_PING_METHOD      1     # 1 outputs for PinG (using Elevpoint outputs), 0 no PING
setenv SMK_SPECELEV_YN      N     # Y uses the indicator for major/minor sources
setenv VELOC_RECALC         N     # Y recalculates velocity from diam and flow

# Script settings
setenv SRCABBR        pt.$FYEAR   # abbreviation for naming log files
setenv QA_TYPE            all     # [none, all, part1-part4, or custom]
setenv PROMPTFLAG         N       # Y (never set to Y for batch processing)
setenv AUTO_DELETE        Y       # Y deletes SMOKE I/O API output files (recommended)
setenv AUTO_DELETE_LOG    Y       # Y automatically deletes logs without asking
setenv DEBUGMODE          N       # Y changes script to use debugger
setenv DEBUG_EXE          idb     # Sets the debugger to use when DEBUGMODE = Y

# Override settings (comment out if not used)
# setenv SPC_OVERRIDE  cmaq.cb4p25  # Chemical mechanism override
# setenv YEAR_OVERRIDE              # Overrides YEAR (base) in Assigns file
# setenv INVTABLE_OVERRIDE          # Inventory table override
# setenv CNTLCASE                   # Control case

##############################################################################

## NOTE: GCNTL file for point sources must be in the $PRDAT directory for the
#        base case simulation and should have the names, as follows:
#            Projection only       :  gcntl.$YEAR_$FYEAR.txt,
#            Control only          :  gcntl.$CNTLCASE.txt,
#            Projection and control:  gcntl.$YEAR_$FYEAR_$CNTLCASE.txt,
#        where YEAR is set in the Assigns file and FYEAR and CNTLCASE
#        are set in this script.

setenv RUN_PART1 Y
source $ASSIGNS_FILE   # Invoke Assigns file

## Set projection and control matrices to use in Grwinven
#      NOTE: Grwinven setting SMK_NUM_CTLMAT > 1 to use PCMAT02-04
setenv PCMAT01 $PPMAT
# setenv PCMAT02 
# setenv PCMAT03 
# setenv PCMAT04 

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
setenv RUN_PART1 N     # Temporarily end part 1
source $ASSIGNS_FILE

setenv RUN_PART1 Y     # Restart part 1
source qa_run.csh      # Run QA for part 1
setenv RUN_PART1 N     # End part 1

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
   setenv G_STDATE_ADVANCE $cnt

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

## Loop through days to run Smkmerge and Smk2emis 
#
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
setenv RUN_PART4 N
unsetenv G_STDATE_ADVANCE

## Ending of script
#
exit( $status )

