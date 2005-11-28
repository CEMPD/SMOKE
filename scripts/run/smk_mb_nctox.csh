#!/bin/csh -f
#
# Version @(#)$Id$
# Path    $Source$
# Date    $Date$
#
# This script sets up needed environment variables for processing on-road mobile
# source emissions in SMOKE and calls the scripts that run the SMOKE
# programs. 
#
# Script created by : M. Houyoux, CEP Environmental Modeling Center 
#
#*********************************************************************

## Set Assigns file name
setenv ASSIGNS_FILE $SMKROOT/assigns/ASSIGNS.nctox.cmaq.cb4p25_wtox.us12-nc

## Set source category
setenv SMK_SOURCE    M          # source category to process
setenv MRG_SOURCE    M          # source category to merge

## Set programs to run...

## Time-independent programs
setenv RUN_SMKINVEN  Y          # run inventory import program
setenv RUN_SPCMAT    Y          # run speciation matrix program
setenv RUN_GRDMAT    Y          # run gridding matrix program
setenv RUN_MBSETUP   Y          # run speed/temperature setup program

## Episode-dependent programs
setenv RUN_PREMOBL   Y          # run temperature processing program
setenv RUN_EMISFAC   Y          # run emission factors program

## Time-dependent programs
setenv RUN_TEMPORAL  Y          # run temporal allocation program
setenv RUN_SMKMERGE  Y          # run merge program

## Quality assurance
setenv RUN_SMKREPORT Y          # run emissions reporting program

## Program-specific controls...

## For Smkinven
setenv FILL_ANNUAL          N   # Y fills annual data with average-day data
setenv IMPORT_VMTMIX_YN     N   # Y imports VMT mix data for use with EMS-95 inventory
setenv RAW_DUP_CHECK        N   # Y checks for duplicate records
setenv SMK_BASEYR_OVERRIDE  0   # year to override the base year of the inventory
setenv SMK_DEFAULT_TZONE    5   # default time zone for sources not in the COSTCY file
setenv SMK_EMS95_FIXFMT     N   # Y indicates that EMS-95 inventory is fixed format
setenv SMK_NHAPEXCLUDE_YN   Y   # Y uses NHAPEXCLUDE file when integrating toxic sources
setenv SMKINVEN_FORMULA     " " # formula for computing emissions value
setenv WEST_HSPHERE         Y   # Y converts longitudes to negative values
setenv WKDAY_NORMALIZE      N   # Y treats average-day emissions as weekday only
setenv WRITE_ANN_ZERO       N   # Y writes zero emission values to intermediate inventory
#      INVNAME1         set by make_invdir.csh script
#      INVNAME2         set by make_invdir.csh scripts
#      OUTZONE          see "Multiple-program controls" below
#      REPORT_DEFAULTS  see "Multiple-program controls" below
#      SMK_MAXERROR     see "Multiple-program controls" below
#      SMK_MAXWARNING   see "Multiple-program controls" below
#      SMK_TMPDIR       set by assigns/set_dirs.scr script

## For Grdmat
setenv SMK_DEFAULT_SRGID    100 # surrogate code number to use as fallback
#      IOAPI_ISPH       set by Assigns file
#      REPORT_DEFAULTS  see "Multiple-program controls" below

## For Spcmat
setenv POLLUTANT_CONVERSION Y   # Y uses the GSCNV pollutant conversion file
#      REPORT_DEFAULTS  see "Multiple-program controls" below

## For Mbsetup
setenv USE_INVSPD_DEFAULT   Y   # Y uses inventory speeds when missing hourly speeds
#      USE_SPEED_PROFILES see "Multiple-program controls" below

## For Premobl
setenv SMK_MAXTEMP          120. # maximum allowable hourly temperature in deg F
setenv SMK_MINTEMP          0.   # minimum allowable hourly temperature in deg F
setenv TVARNAME             TEMPG # temperature variable name in meteorology files
setenv UNIFORM_STIME        -1   # indicates day start time; -1 uses time zones
#      OUTZONE          see "Multiple-program controls" below
#      SMK_METPATH      set by assigns/set_dirs.scr script

## For Emisfac
setenv ADJUST_HR_SPEED      Y   # Y adjusts hourly speeds for freeway ramps
setenv ADJUST_INV_SPEED     Y   # Y adjusts inventory speeds for freeway ramps
setenv M6_FLAT_VMT          Y   # Y uses flat hourly VMT profile when running MOBILE6
setenv MB_HC_TYPE           TOG # name of hydrocarbon pollutant to generate factors for
setenv REPLACE_HUMIDITY     Y   # Y replaces humidity data in MOBILE6 inputs
setenv REPLACE_TEMPERATURES Y   # Y replaces temperature data in MOBILE6 inputs
#      EF_YEAR          set by assigns/set_case.scr script
#      GROUP_TYPE       set to [daily|weekly|monthly|episode] when running Emisfac below
#      SMK_EMISPATH     set by assigns/set_dirs.scr script
#      SMK_M6PATH       set by assigns/set_dirs.scr script
#      SMK_METPATH      set by assigns/set_dirs.scr script
#      SMK_SPDPATH      set by assigns/set_dirs.scr script
#      USE_SPEED_PROFILES see "Multiple-program controls" below

## For Temporal
setenv RENORM_TPROF         Y   # Y normalizes the temporal profiles
setenv UNIFORM_TPROF_YN     N   # Y uses uniform temporal profiles for all sources
setenv ZONE4WM              Y   # Y applies temporal profiles using time zones
#      OUTZONE          see "Multiple-program controls" below
#      REPORT_DEFAULTS  see "Multiple-program controls" below
#      SMK_AVEDAY_YN    see "Multiple-program controls" below
#      SMK_MAXERROR     see "Multiple-program controls" below
#      SMK_MAXWARNING   see "Multiple-program controls" below

## For Smkmerge
setenv MRG_SPCMAT_YN        Y   # Y produces speciated output 
setenv MRG_TEMPORAL_YN      Y   # Y produces temporally allocated output
setenv MRG_GRDOUT_YN        Y   # Y produces a gridded output file
setenv MRG_REPCNY_YN        Y   # Y produces a report of emission totals by county
setenv MRG_REPSTA_YN        Y   # Y produces a report of emission totals by state
setenv MRG_GRDOUT_UNIT      moles/s   # units for the gridded output file
setenv MRG_TOTOUT_UNIT      moles/day # units for the state and county reports
setenv SMK_REPORT_TIME      230000    # hour for reporting daily emissions
#      SMK_AVEDAY_YN    see "Multiple-program controls" below

## For Smkreport
setenv REPORT_ZERO_VALUES   N   # Y outputs entries with all zero values

## Multiple-program controls
setenv OUTZONE              0   # time zone of output emissions
setenv REPORT_DEFAULTS      N   # Y reports sources that use default cross-reference
setenv SMK_AVEDAY_YN        N   # Y uses average-day emissions instead of annual emissions
setenv SMK_MAXERROR         100 # maximum number of error messages in log file
setenv SMK_MAXWARNING       100 # maximum number of warning messages in log file
setenv USE_SPEED_PROFILES   N   # Y uses hourly speed profiles instead of inventory speeds

## Script settings
setenv SRCABBR              mb  # abbreviation for naming log files
setenv QA_TYPE              all # type of QA to perform [none, all, part1, part2, or custom]
setenv PROMPTFLAG           N   # Y prompts for user input
setenv AUTO_DELETE          Y   # Y automatically deletes I/O API NetCDF output files
setenv AUTO_DELETE_LOG      Y   # Y automatically deletes log files
setenv DEBUGMODE            N   # Y runs program in debugger
setenv DEBUG_EXE            pgdbg # debugger to use when DEBUGMODE = Y

## Assigns file override settings
# setenv SPC_OVERRIDE  cmaq.cb4p25  # chemical mechanism override
setenv SPC_SRC_SPECIFIC     Y       # Y uses source specific speciation files
setenv YEAR_OVERRIDE        1999    # base year override
setenv INVTABLE_OVERRIDE    invtable_onroad.cb4.010804.txt # inventory table override

##############################################################################

# NOTE:  The MOBILE6 inputs for the case you want to run
#        should all have file names with a .in extension and should
#        be placed in a *directory* named as follows:
#          No control name    : $INVDIR/mobile/m6_$EF_YEAR/
#          With control name  : $INVDIR/mobile/m6_$EFYEAR_$CNTLCASE/
#        where EF_YEAR is either YEAR if SMK_FUTURE_YN = N and is
#        FYEAR if SMK_FUTURE_YN = Y.  The CNTLCASE variable will only
#        be used in the directory name if it is defined.

## Run Smkinven, Spcmat, Grdmat, Mbsetup, and Premobl
#
setenv RUN_PART1 Y
source $ASSIGNS_FILE   # Invoke Assigns file
source smk_run.csh     # Run programs
source qa_run.csh      # Run QA for part 1
setenv RUN_PART1 N

## Run Emisfac for cases needed by temperature choices in MVREF file
#
foreach group ( episode )  # can include ( daily weekly monthly episode )
  setenv GROUP_TYPE $group
  source emisfac_run.csh     # Run programs
endif

## Loop through days to run Temporal and Smkmerge
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
