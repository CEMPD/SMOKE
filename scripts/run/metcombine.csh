#!/bin/csh -f

# This script sets up environment variables and runs the Metcombine utility.
# Metcombine extracts data from 2D and the 1st layer of 3D meteorology files
# to create custom met files needed by Premobl.
#
# This script assumes that the meteorology files are named with the following
# convention:
#     METCRO2D_<YYYYDDD>
#     METCRO3D_<YYYYDDD>
#
# Script created by C. Seppanen, CEP 
#
#*****************************************************************************

# set Assigns file name
setenv ASSIGNS_FILE $SMKROOT/assigns/ASSIGNS.nctox.cmaq.cb4p25_wtox.us36-nc

# Metcombine settings
setenv VARLIST "TEMPG, PRES, QV"   # list of variables to extract

# Script settings
setenv PROMPTFLAG N                # Y prompts for input

##############################################################################

source $ASSIGNS_FILE

# loop through METCRO2D files in METDAT directory
foreach met2d ( $METDAT/METCRO2D* )

   set date_end = `echo $met2d | wc -c`
   set date_start = $date_end
   @ date_start = $date_start - 7

   set date = `echo $met2d | cut -c${date_start}-${date_end}`

   set met3d = $METDAT/METCRO3D_${date}

   if ( -e $met3d ) then
       setenv LOGFILE  $LOGS/metcombine.$date.log
       setenv METFILE1 $met2d
       setenv METFILE2 $met3d
       setenv OUTFILE  $METDAT/METCOMBO_${date} 

       $SMK_BIN/metcombine
   endif

end

