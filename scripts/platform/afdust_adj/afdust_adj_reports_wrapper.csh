#!/bin/csh -f
#SBATCH --export=NONE

setenv SECTOR canada_afdust
setenv BASE_YEAR 2021
setenv GRID 12US1
setenv EMF_SPC cmaq_cb6ae7
setenv CASE 2021hb_cb6_21k
setenv PROJECT_ROOT /work/EMIS/em_v9/hapcap
setenv IMD_ROOT ${PROJECT_ROOT}/$CASE/premerged
setenv PREMERGED $IMD_ROOT/$SECTOR
setenv RUN_SECT "ADJ"
setenv MRGDATE_ROOT /work/EMIS/em_v6.3/ge_dat/smk_dates/$BASE_YEAR
setenv MONTHS_LIST "1 2 3 4 5 6 7 8 9 10 11 12"
setenv SCRIPTS /work/EMIS/smoke/smoke4.6/scripts
#setenv MONTHS_LIST "2" #5 6 7 8 9"

#/work/EMIS/smoke/smoke4.5/scripts/afdust_adj/afdust_ann_report_test.py
/work/EMIS/smoke/smoke4.6/scripts/afdust_adj/afdust_ann_report.py
