#!/bin/csh -f
#
# Version @(#)$Id$
# Path    $Source$
# Date    $Path$
#
# This script sets up needed environment variables for creating future-year
# point-source emissions and calls the scripts that run the SMOKE programs. 
#
# Script created by : M. Houyoux, CEP Environmental Modeling Center 
#                     February 2003
#*********************************************************************

## Set Assigns file name
setenv ASSIGNS_FILE $SMKROOT/assigns/ASSIGNS.nctox.cmaq.cb4p25_wtox.us12-nc

## Set future year
setenv FYEAR 2018               # year of future case

## Set source category
setenv SMK_SOURCE    P          # source category to process
setenv MRG_SOURCE    P          # source category to merge
setenv MRG_CTLMAT_MULT ' '      # source category to merge multiplicative controls
setenv MRG_CTLMAT_REAC ' '      # source category to merge reactivity controls

## Set programs to run...

## Time-independent programs
setenv RUN_CNTLMAT   Y          # run control and growth matrix program
setenv RUN_GRWINVEN  Y          # run control and growth application program

## Time-dependent programs
setenv RUN_TEMPORAL  Y          # run temporal allocation program
setenv RUN_ELEVPOINT Y          # run elevated/PinG sources selection program
setenv RUN_SMKMERGE  Y          # run merge program

## Quality assurance
setenv RUN_SMKREPORT Y          # run emissions reporting program

# Program-specific controls...

## Program-specific controls...

# For Cntlmat
setenv PROJECTION_YN_SPEC   Y   # Y uses year-specific /PROJECTION/ packet data
setenv REACTIVITY_POL       ' ' # pollutant name to use in computing reactivity matrix
setenv XREF_SICOVERSCC      Y   # Y matches sources by SIC before SCC
#      REPORT_DEFAULTS  see "Multiple-program controls" below
#      SMK_AVEDAY_YN    see "Multiple-program controls" below
#      SMK_TMPDIR       set by assigns/set_dirs.scr script

## For Grwinven
setenv SMK_GRWIDAOUT_YN     N   # Y outputs an IDA-formatted inventory
setenv SMK_GRWORLOUT_YN     Y   # Y outputs an ORL-formatted inventory
setenv SMK_GRWSMKOUT_YN     Y   # Y outputs a SMOKE-formatted inventory
setenv SMK_NUM_CTLMAT       1   # number of control/growth matrices

## For Elevpoint
setenv SMK_ELEV_METHOD      1   # 1 uses PELVCONFIG file to determine elevated sources
setenv UNIFORM_STIME        -1  # indicates day start time; -1 uses time zones
#      IOAPI_ISPH       set by Assigns file
#      SMK_PING_METHOD  see "Multiple-program controls" below

## For Temporal
setenv RENORM_TPROF         Y   # Y normalizes the temporal profiles
setenv UNIFORM_TPROF_YN     N   # Y uses uniform temporal profiles for all sources
setenv ZONE4WM              Y   # Y applies temporal profiles using time zones
#      DAY_SPECIFIC_YN  see "Multiple-program controls" below
#      HOUR_SPECIFIC_YN see "Multiple-program controls" below
#      OUTZONE          see "Multiple-program controls" below
#      REPORT_DEFAULTS  see "Multiple-program controls" below
#      SMK_AVEDAY_YN    see "Multiple-program controls" below
#      SMK_MAXERROR     see "Multiple-program controls" below
#      SMK_MAXWARNING   see "Multiple-program controls" below

## For Smkmerge
setenv MRG_LAYERS_YN        Y   # Y produces layered output
setenv MRG_SPCMAT_YN        Y   # Y produces speciated output 
setenv MRG_TEMPORAL_YN      Y   # Y produces temporally allocated output
setenv MRG_GRDOUT_YN        Y   # Y produces a gridded output file
setenv MRG_REPCNY_YN        Y   # Y produces a report of emission totals by county
setenv MRG_REPSTA_YN        Y   # Y produces a report of emission totals by state
setenv MRG_REPCTL_YN        N   # Y separately reports controlled emissions
setenv MRG_GRDOUT_UNIT      moles/s   # units for the gridded output file
setenv MRG_TOTOUT_UNIT      moles/day # units for the state and county reports
setenv MRG_MARKETPEN_YN     N   # Y uses market penetration from reactivity matrices
setenv SMK_ASCIIELEV_YN     N   # Y creates an ASCII elevated point sources file
setenv SMK_REPORT_TIME      230000    # hour for reporting daily emissions
#      EXPLICIT_PLUMES_YN see "Multiple-program controls" below
#      SMK_AVEDAY_YN    see "Multiple-program controls" below
#      SMK_EMLAYS       see "Multiple-program controls" below
#      SMK_PING_METHOD  see "Multiple-program controls" below

## For Smkreport
setenv REPORT_ZERO_VALUES   N   # Y outputs entries with all zero values

# Multiple-program controls
setenv DAY_SPECIFIC_YN      N   # Y imports and uses day-specific inventory data
setenv EXPLICIT_PLUME_YN    N   # Y processes only sources using explicit plume rise
setenv HOUR_SPECIFIC_YN     N   # Y imports and uses hour-specific inventory data
setenv OUTZONE              0   # time zone of output emissions
setenv REPORT_DEFAULTS      N   # Y reports sources that use default cross-reference
setenv SMK_EMLAYS           12  # number of emissions layers
setenv SMK_AVEDAY_YN        N   # Y uses average-day emissions instead of annual emissions
setenv SMK_MAXERROR         100 # maximum number of error messages in log file
setenv SMK_MAXWARNING       100 # maximum number of warning messages in log file
setenv SMK_PING_METHOD      1   # 1 processes and outputs PinG sources
setenv VELOC_RECALC         N   # Y recalculates stack velocity using flow and diameter

## Script settings
setenv SRCABBR              pt.$FYEAR # abbreviation for naming log files
setenv QA_TYPE              all # type of QA to perform [none, all, part1, part2, or custom]
setenv PROMPTFLAG           N   # Y prompts for user input
setenv AUTO_DELETE          Y   # Y automatically deletes I/O API NetCDF output files
setenv AUTO_DELETE_LOG      Y   # Y automatically deletes log files
setenv DEBUGMODE            N   # Y runs program in debugger
setenv DEBUG_EXE            pgdbg # debugger to use when DEBUGMODE = Y

## Assigns file override settings
# setenv SPC_OVERRIDE  cmaq.cb4p25  # chemical mechanism override
# setenv YEAR_OVERRIDE              # base year override
# setenv INVTABLE_OVERRIDE          # inventory table override
# setenv CNTLCASE                   # control case name

##############################################################################

## NOTE: GCNTL file for point sources must be in the $PTDAT directory for the
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

## Loop through days to run Smkmerge 
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
