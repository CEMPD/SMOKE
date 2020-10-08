#!/bin/csh -f
#
# Version @(#)$Id$
# Path    $Source$
# Date    $Date$
#
# This script sets up needed environment variables for processing biogenic source
# emissions in SMOKE using BEIS3 and calls the scripts that run the SMOKE programs. 
#
# Script created by : M. Houyoux, CEP Environmental Modeling Center 
#
#*********************************************************************

## Set optional customized SMKMERGE output file names
## setenv SMKMERGE_CUSTOM_OUTPUT  N  # Y define your own output file names from SMKMERGE

## Set Assigns file name
setenv ASSIGNS_FILE $SMKROOT/assigns/ASSIGNS.nctox.cmaq.cb05_soa.us12-nc

## Set source category
setenv SMK_SOURCE    B          # source category to process
setenv MRG_SOURCE    B          # source category to merge

setenv BEIS_VERSION  3.7        # version of BEIS3 to use (currently 3.14, 3.61 or 3.7)

## Set programs to run...

## Time-independent programs
setenv RUN_NORMBEIS3 Y          # run normalized biogenic emissions program

## Time-dependent programs
setenv RUN_TMPBEIS3  Y          # run temporal adjustments and speciation program
setenv RUN_SMKMERGE  Y          # run merge program

## Program-specific controls...

## For Normbeis3
#      BEIS_VERSION     already set above

## For Tmpbeis3.60 and 3.14
setenv BG_CLOUD_TYPE        1     # method used to calculate PAR
setenv BIOG_SPRO            B10C5 # speciation profile code to use for biogenics
setenv BIOMET_SAME          Y     # Y indicates temperature and radiation data in same file
setenv BIOSW_YN             N     # Y uses BIOSEASON file to set winter or summer factors (for annual simulations)
setenv SUMMER_YN            Y     # Y assumes summer factors (for short episodic simulations)
setenv OUTZONE              0     # time zone of output emissions
setenv TMPR_VAR             TEMP2 # name of temperature variable
setenv PRES_VAR             PRSFC # name of surface pressure variable
setenv OUT_UNITS            2     # molar output units (1 = moles/hr, 2 = moles/s)
setenv PX_VERSION           Y     # Y indicates that met data is from PX version of MM5 or WRF
setenv SOILT_VAR            SOIT1 # name of soil temperature variable if using PX version
setenv ISLTYP_VAR           SLTYP # name of soil type variable if using PX version
setenv SOILM_VAR            SOIM1 # name of soil moisture variable if using PX version
setenv RN_VAR               RN    # name of non-convective rainfall variable
setenv RC_VAR               RC    # name of convective rainfall variable
setenv INITIAL_RUN          Y     # Y: running first day of scenario, N for subsequent days

## For BEIS3.60 only
setenv RGRND_VAR           RGRND  # solar radiation reaching the surface variable (WATTS/M**2)
setenv RADYNI_VAR          RADYNI # aerodynamic resistance variable     (s/m)
setenv LAI_VAR             LAI    # leaf area index variable from met model (non dimensional)
setenv Q2_VAR              Q2     # 2 meter water vapor mixing ratio variable  (kg/kg)
setenv RSTOMI_VAR          RSTOMI # inverse of the stomatal resistance variable (M/S)
setenv USTAR_VAR           USTAR  # surface friction velocity variable  (M/S)
setenv TEMPG_VAR           TEMPG  # skin temperature at ground variable (K)

## For BEIS3.14 only
setenv RAD_VAR              RGRND # name of radiation/cloud variable for BEI3.14 (see RGRND_VAR for BEIS 3.60)

## For Smkmerge
#      NOTE: Smkmerge run to create state and county emission total reports
setenv AREA_SURROGATE_NUM   340 # surrogate code number for land-area surrogate
setenv MRG_SPCMAT_YN        Y   # Y produces speciated output 
setenv MRG_TEMPORAL_YN      Y   # Y produces temporally allocated output
setenv MRG_GRDOUT_YN        Y   # Y produces a gridded output file
setenv MRG_REPCNY_YN        Y   # Y produces a report of emission totals by county
setenv MRG_REPSTA_YN        Y   # Y produces a report of emission totals by state
setenv MRG_GRDOUT_UNIT      tons/day # units for the gridded output file
setenv MRG_TOTOUT_UNIT      tons/day # units for the state and county reports
setenv SMK_REPORT_TIME      230000   # hour for reporting daily emissions

## Script settings
setenv SRCABBR              bg  # abbreviation for naming log files
setenv PROMPTFLAG           N   # Y prompts for user input
setenv AUTO_DELETE          Y   # Y automatically deletes I/O API NetCDF output files
setenv AUTO_DELETE_LOG      Y   # Y automatically deletes log files
setenv DEBUGMODE            N   # Y runs program in debugger
setenv DEBUG_EXE            pgdbg # debugger to use when DEBUGMODE = Y

## Assigns file override settings
# setenv SPC_OVERRIDE  cmaq.cb4p25  # chemical mechanism override
# setenv INVTABLE_OVERRIDE          # inventory table override

##############################################################################

## Run Normbeis3
#
setenv RUN_PART1 Y
source $ASSIGNS_FILE   # Invoke Assigns file
source smk_run.csh     # Run programs
setenv RUN_PART1 N

## Loop through days to run Tmpbeis3 and Smkmerge
setenv RUN_PART2 Y
setenv RUN_PART4 Y
set cnt = 0
set g_stdate_sav = $G_STDATE
while ( $cnt < $EPI_NDAY )

   if ( $cnt > 0 ) then
      setenv INITIAL_RUN N
   endif

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

#
#
## Ending of script
#
exit( 0 )
