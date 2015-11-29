#! /bin/csh -f

# This script sets up environment variables and runs the Metmoves utility 
# Metmoves processes meteorology files preprocessed by Metcombine utility
# to create custom met files needed by MOVES Driver Scripts and SMOKE-Movesmrg. 
#
# This script assumes that the meteorology files are preprocessed and named with
# the following convention:
#
#     METCOMBO_<YYYYDDD>
#     METCOMBO_<YYYYDDD>
#
# Script created by B.H. Baek, CEMPD (9/20/2011)
#
#*****************************************************************************

# set Assigns file name
setenv ASSIGNS_FILE $SMKROOT/assigns/ASSIGNS.nctox.cmaq.cb05_soa.us12-nc

# Metcombine settings
setenv SRG_LIST  "100, 240"          # list of surrogates to retrieve a list of counties

# Script settings
setenv PROMPTFLAG N                  # Y prompts for input

# Source category
setenv SMK_SOURCE M                  # M : mobile sources
setenv SRCABBR   rateperprofile 

# Define averaging method for RH and Min/Max temperatures
setenv TVARNAME          TEMP2       # chosen temperature variable name

# Define the modeling period
setenv STDATE 2005182                # Starting date in Julian
setenv ENDATE 2005213                # Ending date in Julian 

##############################################################################
# Source your Assign file

source $ASSIGNS_FILE

# Define output file names for MOVES and SMOKE models
setenv SMOKE_OUTFILE $METDAT/SMOKE_${GRID}_${STDATE}-${ENDATE}.txt
setenv MOVES_OUTFILE $METDAT/MOVES_${GRID}_${STDATE}-${ENDATE}.txt

$SMK_BIN/met4moves

# remove temporary output file
rm -f TMP_COMBINED_SRG.txt
