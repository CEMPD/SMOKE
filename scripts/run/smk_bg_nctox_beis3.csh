#!/bin/csh -f

# Version @(#)$Id$
# Path    $Source$
# Date    $Date$

# This script sets up needed environment variables for running biogenic source
# emissions in SMOKE using BEIS3, and calls the script that runs the SMOKE programs. 
#
# Script created by : M. Houyoux, CEP Environmental Modeling Center 
#
#
#*********************************************************************

# set Assigns file name
setenv ASSIGNS_FILE $SMKROOT/assigns/ASSIGNS.nctox.cmaq.cb4p25_wtox.us36-nc

# set source category
setenv SMK_SOURCE B            # source category to process
setenv MRG_SOURCE B            # source category to merge

# time independent programs
setenv RUN_NORMBEIS3    Y        # Y runs normalized emissions program

# time-dependent programs
setenv RUN_TMPBEIS3     Y        # Y runs temporal adjustments and speciation program
setenv RUN_SMKMERGE     Y        # Y runs merge program

# For Tmpbeis3
setenv BG_CLOUD_TYPE    1        # method used to calculate PAR
setenv BIOG_SPRO        BV309    # speciation profile code to use for biogenics
setenv BIOMET_SAME      N        # Y: temperature and rad/cld data in same file   
setenv BIOSW_YN         Y        # Y for using seasons file in Tmpbeis3 
setenv OUTZONE          0        # Output time zone
setenv RAD_VAR          RGRND    # name of radiation/cloud variable
setenv TMPR_VAR         TA       # name of temperature variable
setenv PRES_VAR         PRES     # name of surface pressure variable

# For Smkmerge - getting county reports
setenv MRG_TEMPORAL_YN    Y         # Y merges with hourly emissions
setenv MRG_SPCMAT_YN      Y         # Y merges with speciation matrix
setenv MRG_REPSTA_YN      Y         # Y outputs state totals
setenv MRG_REPCNY_YN      N         # Y outputs county totals
setenv MRG_GRDOUT_UNIT    tons/day  # units for gridded output file
setenv MRG_TOTOUT_UNIT    tons/day  # units for state and/or county totals
setenv AREA_SURROGATE_NUM 3         # number for land-area surrogate

# Script settings
setenv SRCABBR          bg      # abbreviation for naming log files
setenv PROMPTFLAG       N       # Y (never set to Y for batch processing)
setenv AUTO_DELETE      Y       # Y deletes SMOKE I/O API output files (recommended)
setenv AUTO_DELETE_LOG  Y       # Y automatically deletes logs without asking
setenv DEBUGMODE        N       # Y changes script to use debugger
setenv DEBUG_EXE        dbx     # Sets the debugger to use when DEBUGMODE = Y

# Override settings
# setenv SPC_OVERRIDE  cmaq.cb4p25  # Chemical mechanism override 
# setenv INVTABLE_OVERRIDE          # Inventory table override

##############################################################################

## Run Normbeis3
#
setenv RUN_PART1 Y
source $ASSIGNS_FILE   # Invoke Assigns file
source smk_run.csh     # Run programs
setenv RUN_PART1 N

## Loop through days to run Tmpbeis3 and Smkmerge mass reports
setenv RUN_PART2 Y
setenv RUN_PART4 Y
set cnt = 0
set g_stdate_sav = $G_STDATE
while ( $cnt < $EPI_NDAY )

   @ cnt = $cnt + $NDAYS
   source $ASSIGNS_FILE   # Invoke Assigns file to set new dates

   # Override file names to use Tmpbeis3 outputs
   setenv BGTS_L     $B3GTS_L
   setenv BGTS_S     $B3GTS_S
   setenv BGTS_L_O   $B3GTS_L_O
   setenv BGTS_S_O   $B3GTS_S_O
   setenv REPBGTS_L  $REPB3GTS_L
   setenv REPBGTS_S  $REPB3GTS_S

   source smk_run.csh     # Run programs
   setenv G_STDATE_ADVANCE $cnt

end
setenv RUN_PART2 N 
setenv RUN_PART4 N 
unsetenv G_STDATE_ADVANCE

## Additional Smkmerge run for CMAQ
## Loop through days to run Smkmerge units coversion
setenv RUN_PART4 Y 

setenv MRG_GRDOUT_YN      Y         # Y outputs gridded data
setenv MRG_REPSTA_YN      N         # Y outputs state totals
setenv MRG_REPCNY_YN      N         # Y outputs county totals
setenv MRG_GRDOUT_UNIT    moles/s   # units for gridded output file
setenv MRG_TOTOUT_UNIT    moles/hr  # units for state and/or county totals

set cnt = 0
set g_stdate_sav = $G_STDATE
while ( $cnt < $EPI_NDAY )

   @ cnt = $cnt + $NDAYS
   source $ASSIGNS_FILE   # Invoke Assigns file to set new dates

   # Override file names to use Tmpbeis3 outputs
   setenv BGTS_L     $B3GTS_L
   setenv BGTS_S     $B3GTS_S
   setenv BGTS_L_O   $B3GTS_L_O
   setenv BGTS_S_O   $B3GTS_S_O

   source smk_run.csh     # Run programs
   setenv G_STDATE_ADVANCE $cnt

end
setenv RUN_PART2 N
setenv RUN_PART4 N
unsetenv G_STDATE_ADVANCE

#
#
## Ending of script
#
exit( 0 )

