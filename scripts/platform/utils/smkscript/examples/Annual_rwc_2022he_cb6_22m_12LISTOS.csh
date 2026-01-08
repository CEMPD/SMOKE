#!/bin/csh -f
#SBATCH --export=NONE

limit stacksize unlimited

setenv SECTOR "rwc"

if ($?SLURM_SUBMIT_DIR) then
  cd $SLURM_SUBMIT_DIR
endif

## Definitions for case name, directory structures, etc, that are used
#  by every sector in the case
#  Anything defined in directory_definitions can be overridden here
#  if desired
source /proj/ie/proj/SMOKE/htran/12LISTOS/2022he_cb6_22m/scripts/directory_definitions.csh

## Months for emissions processing, and spinup duration
#  In the EPA emissions modeling platforms, the only sectors that use
#    SPINUP_DURATION are biogenics and the final sector merge (Mrggrid).
#  Elsewhere, SPINUP_DURATION = 0, and when Mrggrid runs for spinup days,
#    base year emissions are used for the spinup year for all sectors except
#    biogenics.
#  Effective Jan 2019, SPINUP_DURATION now should work for all months.
#  SPINUP_MONTH_END (new for Jan 2019) specifies whether the last $SPINUP_DURATION
#    days of quarter 2/3/4 should be run at the end of a quarter (Y), or at the start
#    of the next quarter (N). For example, if runningwith SPINUP_DURATION = 10:
#    When N (old behavior), Q1 will include 10 day spinup and end on 3/21; Q2 will
#    cover 3/22 through 6/20. When Y, Q1 will include 10 day spinup and end on 3/31
#    (including all of March), remaining quarters will function as if spinup = 0.
setenv RUN_MONTHS "1 2 3 4 5 6 7 8 9 10 11 12"
setenv SPINUP_DURATION "0"
setenv SPINUP_MONTH_END "Y"

## Emissions modeling year
#  (i.e. meteorological year, not necessarily the inventory year"
setenv BASE_YEAR "2022"
setenv EPI_STDATE_TIME "${BASE_YEAR}-01-01 00:00:00.0"
setenv EPI_ENDATE_TIME "${BASE_YEAR}-12-31 23:59:00.0"

## Inventory case name, if inventories are coming from a different case (they usually aren't)
#  CASEINPUTS is defined in directory_definitions and optionally overridden here
#setenv INVENTORY_CASE "2011ek_cb6v2_v6_11g"
#setenv CASEINPUTS "$INSTALL_DIR/$INVENTORY_CASE/inputs"


## Inputs for all sectors
# Spatial Gridding Surrogates
setenv SRGDESC "${GE_DAT}/gridding/srgdesc_CONUS12_2022v2_04apr2025_12sep2025_v2.txt" # VALID PATH
setenv AGREF "${GE_DAT}/gridding/agref_us_2020platform_13aug2025_nf_v8.txt" # VALID PATH
setenv SRGPRO "${GE_DAT}/gridding/surrogates/CONUS12_2022v2_04apr2025/USA_100_NOFILL.txt" # VALID PATH
# Temporal Profile Cross-reference
setenv ATREF "${GE_DAT}/temporal/tref_general_2022platform_09oct2024_v1" # VALID PATH
setenv PTREF "${GE_DAT}/temporal/tref_general_2022platform_31may2025_nf_v4" # VALID PATH
setenv MTREF "${GE_DAT}/temporal/mtref_onroad_MOVES_2017NEI_03sep2020_v1" # VALID PATH
setenv MTREF_1 "${GE_DAT}/temporal/mtref_2022v2_onroad_MOVES_month01_04apr2025_v0" # VALID PATH
setenv MTREF_12 "${GE_DAT}/temporal/mtref_2022v2_onroad_MOVES_month12_04apr2025_v0" # VALID PATH
setenv MTREF_5 "${GE_DAT}/temporal/mtref_2022v2_onroad_MOVES_month05_04apr2025_v0" # VALID PATH
setenv MTREF_9 "${GE_DAT}/temporal/mtref_2022v2_onroad_MOVES_month09_04apr2025_v0" # VALID PATH
setenv MTREF_2 "${GE_DAT}/temporal/mtref_2022v2_onroad_MOVES_month02_04apr2025_v0" # VALID PATH
setenv MTREF_3 "${GE_DAT}/temporal/mtref_2022v2_onroad_MOVES_month03_04apr2025_v0" # VALID PATH
setenv MTREF_4 "${GE_DAT}/temporal/mtref_2022v2_onroad_MOVES_month04_04apr2025_v0" # VALID PATH
setenv MTREF_6 "${GE_DAT}/temporal/mtref_2022v2_onroad_MOVES_month06_04apr2025_v0" # VALID PATH
setenv MTREF_10 "${GE_DAT}/temporal/mtref_2022v2_onroad_MOVES_month10_04apr2025_v0" # VALID PATH
setenv MTREF_11 "${GE_DAT}/temporal/mtref_2022v2_onroad_MOVES_month11_04apr2025_v0" # VALID PATH
setenv MTREF_8 "${GE_DAT}/temporal/mtref_2022v2_onroad_MOVES_month08_04apr2025_v0" # VALID PATH
setenv MTREF_7 "${GE_DAT}/temporal/mtref_2022v2_onroad_MOVES_month07_04apr2025_v0" # VALID PATH
# Temporal Profile Definitions and Split Factors
setenv MTPRO_MONTHLY "${GE_DAT}/temporal/tpro_monthly_general_2022platform_13jun2024_nf_v1" # VALID PATH
setenv ATPRO_MONTHLY "${GE_DAT}/temporal/tpro_monthly_general_2022platform_13jun2024_nf_v1" # VALID PATH
setenv ATPRO_WEEKLY "${GE_DAT}/temporal/tpro_weekly_general_2022platform_26nov2024_v1" # VALID PATH
setenv ATPRO_HOURLY "${GE_DAT}/temporal/tpro_hourly_general_2022platform_26nov2024_v1" # VALID PATH
setenv PTPRO_WEEKLY "${GE_DAT}/temporal/tpro_weekly_general_2022platform_26nov2024_v1" # VALID PATH
setenv PTPRO_MONTHLY "${GE_DAT}/temporal/tpro_monthly_general_2022platform_13jun2024_nf_v1" # VALID PATH
setenv PTPRO_HOURLY "${GE_DAT}/temporal/tpro_hourly_general_2022platform_26nov2024_v1" # VALID PATH
setenv MTPRO_WEEKLY "${GE_DAT}/temporal/mtpro_weekly_onroad_MOVES_2022v2_04apr2025_v0" # VALID PATH
setenv MTPRO_HOURLY "${GE_DAT}/temporal/mtpro_hourly_onroad_MOVES_2022v2_01apr2025_v0" # VALID PATH
setenv ATPRO_HOURLY_NCF "${GE_DAT}/temporal/Gentpro_TPRO_HOUR_HOURLY_BASH_NH3.agNH3_bash_2022m_12US1.ncf" # VALID PATH
setenv ATPRO_DAILY "${GE_DAT}/temporal/atpro_daily_Gentpro_RWC_2022hc_19jun2024_nf_v1" # VALID PATH
# Point Sources Parameters
setenv ARTOPNT "${GE_DAT}/artopnt_2002detroit_20aug2019_v2.txt" # VALID PATH
setenv PELVCONFIG "${GE_DAT}/point/pelvconfig_elevate_everything_17apr2020_v0.txt" # VALID PATH
setenv ORISDESC "${GE_DAT}/smkreport/orisdesc_04dec2006_v0.txt" # VALID PATH
setenv PSTK "${GE_DAT}/point/pstk_13nov2018_v1.txt" # VALID PATH
# Grid Descriptions
# (https://www.cmascenter.org/ioapi/documentation/all_versions/html/GRIDDESC.html)
setenv GRIDDESC "${GE_DAT}/gridding/griddesc_lambertonly_18jan2019_v7.txt" # VALID PATH
# Merge Date Settings
# (https://github.com/CEMPD/SMOKE/wiki/A.-Overall-Instructions-on-Running-SMOKE-using-EPA's-Emissions-Modeling-Platforms#mrgdate_files)
setenv MRGDATE_FILES "${INSTALL_DIR}/smoke5.2/scripts/smk_dates/2022/smk_merge_dates_202112.txt" # VALID PATH
# Other settings
setenv COSTCY "${GE_DAT}/costcy_for_2017platform_20dec2023_v8.txt" # VALID PATH
setenv SCCDESC "${GE_DAT}/smkreport/sccdesc_2021platform_29aug2025_nf_v1.txt" # VALID PATH
setenv HOLIDAYS "${GE_DAT}/temporal/holidays_19dec2023_v5.txt" # VALID PATH
setenv GSCNV "${GE_DAT}/speciation/gscnv_SPECIATE_5_3_fromS2S_20230911_16jun2025_v1.txt" # VALID PATH
setenv NAICSDESC "${GE_DAT}/smkreport/naicsdesc_20220222_29feb2024_nf_v1.txt" # VALID PATH
setenv EFTABLES "${CASEINPUTS}/onroad/eftables/rateperdistance_smoke_aq_cb6_saprc_20250402_rates2022v2-2022-20250306_10003_1.csv" # MISSING PATH
setenv INVTABLE "${GE_DAT}/invtable_2022platform_integrate_15sep2025_nf_v3.txt" # VALID PATH

## Inputs specific to this sector and/or job
# Speciation Profile Cross-reference
setenv GSREFTMP_X1A "${GE_DAT}/speciation/gsref_minnesota_metals_05dec2017_v0.txt" # VALID PATH
setenv GSREFTMP_X5B "${GE_DAT}/speciation/gsref_hg_2017platform_geothermal_06apr2020_v0.txt" # VALID PATH
setenv GSREF_SUMMER "${GE_DAT}/speciation/gsref_MOVES5_custom_speciation_2022hd_all_nostarts_summer_15apr2025_v0.txt" # VALID PATH
setenv GSREF_WINTER "${GE_DAT}/speciation/gsref_MOVES5_custom_speciation_2022hd_all_nostarts_winter_15apr2025_v0.txt" # VALID PATH
setenv GSREFTMP_A "${GE_DAT}/speciation/gsref_static_cap_pf4_2014v1_platform_23jan2017_v0.txt" # VALID PATH
setenv GSREFTMP_C "${GE_DAT}/speciation/gsref_static_nox_hono_pf4_2014v1_platform_04jan2024_nf_v4.txt" # VALID PATH
setenv GSREFTMP_B "${GE_DAT}/speciation/gsref_sulf__2014v1_platform_23jan2017_v0.txt" # VALID PATH
setenv GSREFTMP_E "${GE_DAT}/speciation/gsref_pm25_2022platform_04jun2025_v4.txt" # VALID PATH
setenv GSREFTMP_F "${GE_DAT}/speciation/gsref_voc_2022platform_16jun2025_v3.txt" # VALID PATH
setenv GSREFTMP_H "${GE_DAT}/speciation/gsref_nonhapvoc_2022platform_10apr2025_nf_v2.txt" # VALID PATH
setenv GSREFTMP_X1 "${GE_DAT}/speciation/gsref_metals_chromium_stationary_2014v1_platform_17jan2017_v0.txt" # VALID PATH
setenv GSREFTMP_X4 "${GE_DAT}/speciation/gsref_dieselpm_2014v1_platform_15apr2020_nf_v2.txt" # VALID PATH
setenv GSREFTMP_X5 "${GE_DAT}/speciation/gsref_hg_2017_2021_platforms_02jun2025_nf_v6.txt" # VALID PATH
setenv GSREFTMP_H1 "${GE_DAT}/speciation/gsref_nonhapvoc_2022hd_oilgas_basin_specific_ramboll_16apr2025_v0.txt" # VALID PATH
setenv GSREFTMP_F1 "${GE_DAT}/speciation/gsref_voc_2022hd_oilgas_basin_specific_ramboll_16apr2025_v0.txt" # VALID PATH
setenv GSREFTMP_H2 "${GE_DAT}/speciation/gsref_nonhapvoc_2022hd_oilgas_combo_ERG_16apr2025_v0.txt" # VALID PATH
setenv GSREFTMP_F2 "${GE_DAT}/speciation/gsref_voc_2022hd_oilgas_combo_ERG_16apr2025_v0.txt" # VALID PATH
setenv GSREF_STARTS_1 "${GE_DAT}/speciation/gsref_MOVES5_custom_speciation_2022hd_all_starts_m1_21apr2025_v0.txt" # VALID PATH
setenv GSREF_STARTS_2 "${GE_DAT}/speciation/gsref_MOVES5_custom_speciation_2022hd_all_starts_m2_21apr2025_v0.txt" # VALID PATH
setenv GSREF_STARTS_4 "${GE_DAT}/speciation/gsref_MOVES5_custom_speciation_2022hd_all_starts_m4_21apr2025_v0.txt" # VALID PATH
setenv GSREF_STARTS_3 "${GE_DAT}/speciation/gsref_MOVES5_custom_speciation_2022hd_all_starts_m3_21apr2025_v0.txt" # VALID PATH
setenv GSREF_STARTS_5 "${GE_DAT}/speciation/gsref_MOVES5_custom_speciation_2022hd_all_starts_m5_21apr2025_v0.txt" # VALID PATH
setenv GSREF_STARTS_6 "${GE_DAT}/speciation/gsref_MOVES5_custom_speciation_2022hd_all_starts_m6_21apr2025_v0.txt" # VALID PATH
setenv GSREF_STARTS_7 "${GE_DAT}/speciation/gsref_MOVES5_custom_speciation_2022hd_all_starts_m7_21apr2025_v0.txt" # VALID PATH
setenv GSREF_STARTS_8 "${GE_DAT}/speciation/gsref_MOVES5_custom_speciation_2022hd_all_starts_m8_21apr2025_v0.txt" # VALID PATH
setenv GSREF_STARTS_9 "${GE_DAT}/speciation/gsref_MOVES5_custom_speciation_2022hd_all_starts_m9_21apr2025_v0.txt" # VALID PATH
setenv GSREF_STARTS_10 "${GE_DAT}/speciation/gsref_MOVES5_custom_speciation_2022hd_all_starts_m10_21apr2025_v0.txt" # VALID PATH
setenv GSREF_STARTS_11 "${GE_DAT}/speciation/gsref_MOVES5_custom_speciation_2022hd_all_starts_m11_21apr2025_v0.txt" # VALID PATH
setenv GSREF_STARTS_12 "${GE_DAT}/speciation/gsref_MOVES5_custom_speciation_2022hd_all_starts_m12_21apr2025_v0.txt" # VALID PATH
# Speciation Profiles Definitons with Split Factors
setenv GSPROTMP_A "${GE_DAT}/speciation/gspro_static_cmaq_21feb2012_v13.txt" # VALID PATH
setenv GSPROTMP_C "${GE_DAT}/speciation/gspro_PM2_5_AE6_Spec_5_3_S2S_20230911_01nov2023_v0.txt" # VALID PATH
setenv GSPROTMP_B "${GE_DAT}/speciation/gspro_sulf_29jun2007_v1.txt" # VALID PATH
setenv GSPROTMP_D "${GE_DAT}/speciation/gspro_CB6R3_AE7_nointegrate_Spec_5_3_S2S_20230911_01nov2023_nf_v1.txt" # VALID PATH
setenv GSPROTMP_E "${GE_DAT}/speciation/gspro_CB6R3_AE7_integrate_Spec_5_3_S2S_20230911_01nov2023_v1.txt" # VALID PATH
setenv GSPROTMP_I "${GE_DAT}/speciation/gspro_speciated_voc_cb6cmaq_21jul2017_v1.txt" # VALID PATH
setenv GSPROTMP_G "${GE_DAT}/speciation/gspro_nox_hono_pf4_06aug2008_v0.txt" # VALID PATH
setenv GSPROTMP_F "${GE_DAT}/speciation/gspro_haps_2022platform_03jun2024_v0.txt" # VALID PATH
setenv GSPROTMP_M "${GE_DAT}/speciation/gspro_NMOG_HAPs_only_22dec2022_v1.txt" # VALID PATH
setenv GSPROTMP_E1 "${GE_DAT}/speciation/gspro_CB6R3_AE7_integrate_tracer_Spec_5_3_S2S_20230911_01nov2023_v0.txt" # VALID PATH
setenv GSPROTMP_X1 "${GE_DAT}/speciation/gspro_hapmetals_05dec2017_nf_v2.txt" # VALID PATH
setenv GSPROTMP_X2 "${GE_DAT}/speciation/gspro_chromium_05dec2017_nf_v2.txt" # VALID PATH
setenv GSPROTMP_X3 "${GE_DAT}/speciation/gspro_hg_2017platform_16apr2020_nf_v1.txt" # VALID PATH
setenv GSPROTMP_X4 "${GE_DAT}/speciation/gspro_pahs_2014v1_nata_16jul2024_v1.txt" # VALID PATH
setenv GSPROTMP_X5 "${GE_DAT}/speciation/gspro_diesel_pm_2014v1_nata_21oct2016_nf_v1.txt" # VALID PATH
setenv GSPROTMP_D1 "${GE_DAT}/speciation/gspro_CB6R3_AE7_criteria_tracer_Spec_5_3_S2S_20230911_01nov2023_v1.txt" # VALID PATH
# Temporal Profile Cross-reference
setenv ATREF "${GE_DAT}/temporal/atref_2020platform_rwc_10apr2025_nf_v4" # VALID PATH ; DUPLICATED SETTING
# Temporal Profile Definitions and Split Factors
setenv ATPRO_MONTHLY "${GE_DAT}/temporal/atpro_monthly_Gentpro_RWC_2022hc_19jun2024_nf_v1" # VALID PATH ; DUPLICATED SETTING
# SMOKE Report Configurations
setenv REPCONFIG_INV "${GE_DAT}/smkreport/repconfig/repconfig_area_inv_2016beta_15nov2021_nf_v3.txt" # VALID PATH
setenv REPCONFIG_GRID "${GE_DAT}/smkreport/repconfig/repconfig_area_inv_grid_2016beta_15nov2021_nf_v2.txt" # VALID PATH
setenv REPCONFIG_TEMP "${GE_DAT}/smkreport/repconfig/repconfig_area_temporal_2016beta_15nov2021_nf_v1.txt" # VALID PATH
setenv REPCONFIG_INV2 "${GE_DAT}/smkreport/repconfig/repconfig_area_inv_nonhapvocprof_2016beta_15nov2021_nf_v1.txt" # VALID PATH
# Emission Inputs
setenv EMISINV_A "${CASEINPUTS}/rwc/rwc_2022v2_platform_EPA_09apr2025_nf_v1_12LISTOS.csv" # VALID PATH
setenv EMISINV_B "${CASEINPUTS}/rwc/2022hd_from_rwc_2020NEI_NONPOINT_20230222_SLT_05aug2025_nf_v1_12LISTOS.csv" # VALID PATH


## Parameters for all sectors
# L_TYPE and M_TYPE
# (https://github.com/CEMPD/SMOKE/wiki/A.-Overall-Instructions-on-Running-SMOKE-using-EPA's-Emissions-Modeling-Platforms#l_type-and-m_type)
setenv L_TYPE "mwdss" # # Temporal type
setenv M_TYPE "mwdss" # # Merge type
# Other settings
setenv M3XTRACT_LENGTH "240000" # # Length of m3xtract extractions
setenv PLATFORM "v9.2" # # Platform name
setenv SPC "$EMF_SPC" # # Speciation type name
setenv FILL_ANNUAL "N" # # Fill annual values
setenv SMKINVEN_FORMULA "PMC=PM10-PM2_5" # # Formula for Smkinven
setenv SMK_DEFAULT_SRGID "100" # # Default surrogate code
setenv REPORT_DEFAULTS "Y" # # Report default profiles used
setenv SMK_AVEDAY_YN "N" # # Use average day emissions
setenv SMK_MAXWARNING "10" # # Maximum warnings printed
setenv SMK_MAXERROR "10000" # # Maximum errors printed
setenv FULLSCC_ONLY "Y" # # Match full SCCs
setenv RUN_HOLIDAYS "Y" # # Run holidays
setenv RAW_DUP_CHECK "N" # # Check for duplicate sources
setenv SMK_SPECELEV_YN "Y" # # Laypoint uses Elevpoint to set sources for plume rise calc
setenv SMK_PING_METHOD "0" # # Plume-in-grid method
setenv INLINE_MODE "only" # # Run in "inline" mode
setenv RUN_PYTHON_ANNUAL "Y" # # Run script for Smkmerge annual totals
setenv NO_SPC_ZERO_EMIS "Y" # # Don't speciate zero emission SCCs
setenv WRF_VERSION "4" # # WRF version (affects afdust adjustments)
setenv DEFAULT_CONV_FAC_ERROR "Y" # # Error if profile missing from GSCNV
setenv OUTPUT_FORMAT "$EMF_AQM" # # Model output format
setenv SMKMERGE_CUSTOM_OUTPUT "Y" # # Custom merge output
setenv MRG_REPSTA_YN "N" # # Output state totals
setenv IOAPI_ISPH "20" # # I/O API Sphere type
setenv POLLUTANT_CONVERSION "Y" # # Use pollutant conversion
setenv RENORM_TPROF "Y" # # Renormalize temporal profiles
setenv OUTZONE "0" # # Output time zone
setenv MRG_REPCNY_YN "Y" # # Output county totals
setenv WEST_HSPHERE "Y" # # Western hemisphere?
setenv MRG_MARKETPEN_YN "N" # # Include market penetration


## Parameters specific to this sector and/or job
# L_TYPE and M_TYPE
# (https://github.com/CEMPD/SMOKE/wiki/A.-Overall-Instructions-on-Running-SMOKE-using-EPA's-Emissions-Modeling-Platforms#l_type-and-m_type)
setenv L_TYPE "all" # # Temporal type ; DUPLICATED SETTING
setenv M_TYPE "all" # # Merge type ; DUPLICATED SETTING
# Other settings
setenv NONHAP_TYPE "VOC" # # Nonhap Type
setenv SMK_PROCESS_HAPS "ALL" # # Use NHAPEXCLUDE file


$RUNSCRIPTS/emf/smk_ar_annual_emf.csh $REGION_ABBREV $REGION_IOAPI_GRIDNAME -m "$RUN_MONTHS" $SPINUP_DURATION all
