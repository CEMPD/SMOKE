#!/bin/bash

# location of AERMOD-source report created by Smkreport
export REPORT=$REPSCEN/rep_np_oilgas_2011ek_cb6v2_v6_11g_aermod.txt

# locations of source groups data files
export GROUP_PARAMS=$SCRIPTS/aermod/nonpoint_rungroup_parameters.csv
export SOURCE_GROUPS=$SCRIPTS/aermod/aermod_groups_2014_25jan2018.csv

# locations of temporal profile files
export ATPRO_MONTHLY=$GE_DAT/amptpro_general_2011platform_tpro_monthly_20nov2015_v2
export ATPRO_WEEKLY=$GE_DAT/amptpro_general_2011platform_tpro_weekly_13nov2014_v1
export ATPRO_HOURLY=$GE_DAT/amptpro_general_2011platform_tpro_hourly_13nov2014_v1

# directory where output files should be created
export OUTPUT_DIR=$OUTPUT/aermod

$SCRIPTS/aermod/np_oilgas.pl
