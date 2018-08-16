#!/bin/bash

# locations of AERMOD-source reports created by Smkreport
export REPORT_JAN=$REPSCEN/rep_nonroad_jan_2011ek_cb6v2_v6_11g_aermod.txt
export REPORT_FEB=$REPSCEN/rep_nonroad_feb_2011ek_cb6v2_v6_11g_aermod.txt
export REPORT_MAR=$REPSCEN/rep_nonroad_mar_2011ek_cb6v2_v6_11g_aermod.txt
export REPORT_APR=$REPSCEN/rep_nonroad_apr_2011ek_cb6v2_v6_11g_aermod.txt
export REPORT_MAY=$REPSCEN/rep_nonroad_may_2011ek_cb6v2_v6_11g_aermod.txt
export REPORT_JUN=$REPSCEN/rep_nonroad_jun_2011ek_cb6v2_v6_11g_aermod.txt
export REPORT_JUL=$REPSCEN/rep_nonroad_jul_2011ek_cb6v2_v6_11g_aermod.txt
export REPORT_AUG=$REPSCEN/rep_nonroad_aug_2011ek_cb6v2_v6_11g_aermod.txt
export REPORT_SEP=$REPSCEN/rep_nonroad_sep_2011ek_cb6v2_v6_11g_aermod.txt
export REPORT_OCT=$REPSCEN/rep_nonroad_oct_2011ek_cb6v2_v6_11g_aermod.txt
export REPORT_NOV=$REPSCEN/rep_nonroad_nov_2011ek_cb6v2_v6_11g_aermod.txt
export REPORT_DEC=$REPSCEN/rep_nonroad_dec_2011ek_cb6v2_v6_11g_aermod.txt

# locations of source groups data files
export GROUP_PARAMS=$SCRIPTS/aermod/nonpoint_rungroup_parameters.csv
export SOURCE_GROUPS=$SCRIPTS/aermod/aermod_groups_2014_25jan2018.csv

# locations of temporal profile files
export ATPRO_MONTHLY=$GE_DAT/amptpro_general_2011platform_tpro_monthly_20nov2015_v2
export ATPRO_WEEKLY=$GE_DAT/amptpro_general_2011platform_tpro_weekly_13nov2014_v1
export ATPRO_HOURLY=$GE_DAT/amptpro_general_2011platform_tpro_hourly_13nov2014_v1

# directory where output files should be created
export OUTPUT_DIR=$OUTPUT/aermod

$SCRIPTS/aermod/nonroad.pl
