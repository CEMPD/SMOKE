#!/bin/tcsh -fx
#  Script takes GRIDCRO landuse file for any day and generate gridded transport
#  fractions.   The FORTRAN code is located at:
#  /home/jvukovic/apps/smoke_v45/subsys/smoke/src/jmv_tools/gen_afdust_tfrac.f
#
#  Input files are the following:
#    GRIDCRO2D file (env var is BELD4 cause we initially that is what we tried to use)
#    The gridcro2d_lu_types.12US2.csv file maps GRIDCRO2D landuse types to capture 
#      classes 
#    The captureclass_fractions.txt file assigns transport fractions to each capture 
#      class 
#
#    LU_PCT_YN is input environment variable set to N since GRIDCRO2D landuse is in 
#      fractions instead of percent (BELD4 is in percent)
#
#    NOTE: for MCIP 4.5 and earlier: GRIDCRO2D contains LUFRAC variables, therefore only GRIDCRO2D is required;
#          MCIP v5.0 and later: LUFRAC variables are not in GRIDCRO2D but written out to LUFRAC_CRO file, 
#          thus both LUFRACRO and GRIDCRO2D are required
#
#  Output file is:
#    Gridded transport fractions on same grid as GRICRO2D file
#
#  Created by J. Vukovich  Sept 2018
#  Updated by H. Tran (UNC-IE) Nov 2024: Improvised for working with MCIP v5.0 and later output files
#
########################################################################

set wkdir = /proj/ie/proj/SMOKE/htran/smoke_training_Oct2024/2018gg_18j/scripts/utilities/afdust_xportfrac

setenv PROMPTFLAG N

setenv YYMMDD 210121  # Specific any day that has valid GRIDCRO2D (and LUFRAC_CRO, if applicable) 

setenv INPDIR /proj/ie/proj/EDF-OG/SMOKE/platform_2021hb/met
setenv OUTDIR /proj/ie/proj/EDF-OG/SMOKE/platform_2021hb/met/xportfrac

setenv BDIR            $wkdir

setenv LUFRACRO        $INPDIR/LUFRAC_CRO.12US1.35L.$YYMMDD # Required for MCIP 5.0 and later; If missing or undefined; expects LUFRAC variables to be in GRIDCRO2D
setenv GRIDCRO2D       $INPDIR/GRIDCRO2D.12US1.35L.$YYMMDD  # Required

setenv BELD4TOCAPTURE  $wkdir/gridcro2d_lu_types.12US2.csv  # can be used for any grid
setenv CAPFRACS        $wkdir/captureclass_fractions.txt
setenv LU_PCT_YN       N   # Y = use LUFRAC in percentage (%); N = use LUFRAC as fractions

setenv OUTFILE         $OUTDIR/xportfrac.12US1.35L.$YYMMDD

$wkdir/gen_afdust_tfrac_v2

exit
