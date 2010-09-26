#!/usr/bin/perl
#
# Filename   : runspec_generator_v0.1.pl
# Author     : Michele Jimenez, ENVIRON International Corp.
# Version    : 0.1
# Description: Generate RunSpec control files input to MOVES2010
#
#  In order to optimize computing times in MOVES2010 generate the
#  fewest number of runs that will produce all the necessary
#  emission factors.
#
# 	RunSpec files generated for:
#	a) on-road operating; rate per distance table
#       b) off-network processes; rate per vehicle table
#	c) vapor venting off-network; rate per profile table
#
#======================================================================
#= Runspec Generator - a MOVES preprocessor utility
#=
#= Copyright (C) 2010 ENVIRON International Corporation
#=
#= The Runspec Generator is free software; you can redistribute it 
#= and/or modify it under the terms of the GNU General Public License
#= as published by the Free Software Foundation; either version 3
#= of the License, or (at your option) any later version.
#=
#= This utility is distributed in the hope that it will be useful,
#= but WITHOUT ANY WARRANTY; without even the implied warranty of
#= MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#= GNU General Public License for more details.
#=
#= You should have received a copy of the GNU General Public License
#= along with this program.  If not, see <http://www.gnu.org/licenses/>.
#=======================================================================

use strict;
use FileHandle;

($#ARGV >= 1) or die "Usage: runspec_generator.pl RunControlfile RepCounty [-csh]\n";
# input file variables
my ($RunControlFile, $RepCntyFile, $GenCSH, $batchFile, $importBatchFile, $dbListFile, $fileout);

$RunControlFile = $ARGV[0];
$RepCntyFile = $ARGV[1];
$GenCSH = ($#ARGV > 1 && $ARGV[2] eq "-csh") ? 1 : 0;

# run control variables
my ($line, $i, $j, $ip, $ic, $jj, @line);
my ($dbhost, $batchrun, $outdir, $moveshome, $modelyear, $User_polls, $dayofweek, $MetFile);
my (@User_polls, @dayofweek);
my ($pollsFlg, $WeekDayFlag, $WeekEndFlag);

# repcounty variables
my ($cntRepCnty, @repFips, @repAge, @repIM, @repFuelSup, @repFuelForm, @repPop, @repVMT);

my ($default_dummy);
$default_dummy = "filename_dummy_holder.csv";

#=========================================================================================================
# Set Parameters
#=========================================================================================================
#
# MOVES2010 speed bins --------------------------------------------------------------

# MOVES2010 pollutant list ----------------------------------------------------------
#  poll_ID   poll_name Dependence  Pollutant
#     1      THC                   Total Gaseous Hydrocarbons
#     79     NMHC      THC, CH4    Non-Methane Hydrocarbons
#     80     NMOG      NMHC        Non-Methane Organic Gases
#     86     TOG       NMOG, CH4   Total Organic Gases
#     87     VOC       NMHC        Volatile Organic Compounds
#     2      CO                    Carbon Monoxide (CO)
#     3      NOX                   Oxides of Nitrogen
#     30     NH3                   Ammonia (NH3)
#     32     NO        NOX         Nitrogen Oxide
#     33     NO2       NOX         Nitrogen Dioxide
#     31     SO2       Tot Energy  Sulfur Dioxide (SO2)
#     100    TOTPM10               Primary Exhaust PM10  - Total
#     101    OCARB10   OCARB2_5    Primary PM10 - Organic Carbon
#     102    ECARB10   ECARB2_5    Primary PM10 - Elemental Carbon
#     105    SO4_10    Tot Energy  Primary PM10 - Sulfate Particulate
#     106    BRAKE10               Primary PM10 - Brakewear Particulate
#     107    TIRE10                Primary PM10 - Tirewear Particulate
#     110    TOTPM2_5              Primary Exhaust PM2.5 - Total
#     111    OCARB2_5              Primary PM2.5 - Organic Carbon
#     112    ECARB2_5              Primary PM2.5 - Elemental Carbon
#     115    SO4_2_5   Tot Energy  Primary PM2.5 - Sulfate Particulate
#     116    BRAKE2_5              Primary PM2.5 - Brakewear Particulate
#     117    TIRE2_5               Primary PM2.5 - Tirewear Particulate
#     91                           Total Energy Consumption
#     92                           Petroleum Energy Consumption
#     93                           Fossil Fuel Energy Consumption
#     5      CH4                   Methane (CH4)
#     6      N20                   Nitrous Oxide (N2O)
#     90     CO2       Tot Energy  Atmospheric CO2
#     98     CO2EQ                 CO2 Equivalent
#     20     BENZ      VOC         Benzene
#     21     ETHA      VOC         Ethanol
#     22     MTBE      VOC         MTBE
#     23     NAPH      THC,TOTPM10 Naphthalene
#     24     BUTA      VOC         1,3-Butadiene
#     25     FORM      VOC         Formaldehyde
#     26     ACET      VOC         Acetaldehyde
#     27     ACRO      VOC         Acrolein

my ( @pollOptions, @pollsList, @pollsListID, @pollsListName);
my ( @pollsByOptionOZONE, @pollsByOptionTOXICS, @pollsByOptionPM, @pollsByOptionGHG );
my ( @pollsOutList );

@pollsList = ( "THC", "NMHC", "NMOG", "TOG", "VOC", "CO", "NOX", "NH3", "NO", "NO2",
               "SO2", "TOTPM10", "OCARB10", "ECARB10", "SO4_10", "BRAKE10", "TIRE10", "TOTPM2_5", "OCARB2_5", "ECARB2_5",
               "SO4_2_5", "BRAKE2_5", "TIRE2_5", "CH4", "N20", "CO2", "CO2EQ", "BENZ", "ETHA", "MTBE",
               "NAPH", "BUTA", "FORM", "ACET", "ACRO", "TENERGY");

@pollsListID = ( 1, 79, 80, 86, 87, 2, 3, 30, 32, 33, 31, 100, 101, 102, 105, 106, 107, 110, 111,
               112, 115, 116, 117, 91, 92, 93, 5, 6, 90, 98, 20, 21, 22, 23, 24, 25, 26, 27);

@pollsListName = ("Total Gaseous Hydrocarbons",
                 "Non-Methane Hydrocarbons",
                 "Non-Methane Organic Gases",
                 "Total Organic Gases",
                 "Volatile Organic Compounds",
                 "Carbon Monoxide (CO)",
                 "Oxides of Nitrogen",
                 "Ammonia (NH3)",
                 "Nitrogen Oxide",
                 "Nitrogen Dioxide",
                 "Sulfur Dioxide (SO2)",
                 "Primary Exhaust PM10  - Total",
                 "Primary PM10 - Organic Carbon",
                 "Primary PM10 - Elemental Carbon",
                 "Primary PM10 - Sulfate Particulate",
                 "Primary PM10 - Brakewear Particulate",
                 "Primary PM10 - Tirewear Particulate",
                 "Primary Exhaust PM2.5 - Total",
                 "Primary PM2.5 - Organic Carbon",
                 "Primary PM2.5 - Elemental Carbon",
                 "Primary PM2.5 - Sulfate Particulate",
                 "Primary PM2.5 - Brakewear Particulate",
                 "Primary PM2.5 - Tirewear Particulate",
                 "Total Energy Consumption",
                 "Petroleum Energy Consumption",
                 "Fossil Fuel Energy Consumption",
                 "Methane (CH4)",
                 "Nitrous Oxide (N2O)",
                 "Atmospheric CO2",
                 "CO2 Equivalent",
                 "Benzene",
                 "Ethanol",
                 "MTBE",
                 "Naphthalene",
                 "1,3-Butadiene",
                 "Formaldehyde",
                 "Acetaldehyde",
                 "Acrolein");

@pollOptions = ("OZONE", "PM", "TOXICS", "GHG");

#  A subset of MOVES2010 pollutants are generated for each user option specified.
#  Taken from the design document of Task4, Table 4.
@pollsByOptionOZONE = (1,5,79,80,86,87,2,3,32,33);
@pollsByOptionTOXICS = (1,5,79,80,86,87,20,22,23,24,25,26,27,100,101,102,105,91);
@pollsByOptionPM = (1,5,79,80,86,87,3,30,32,33,31,100,101,102,105,106,107,110,111,112,115,116,117,91,20);
@pollsByOptionGHG = (90,91,92,93,5,6,98);

# Process Types ------------------------------------------------------------------
my (@processID, %processName, %PollProc_tablemap);

@processID = (1,2,9,10,11,12,13,15,16,17,90);    # 18,19 refueling not included
$processName{"01"} = "Running Exhaust";
$processName{"02"} = "Start Exhaust";
$processName{"09"} = "Brakewear";
$processName{"10"} = "Tirewear"; 
$processName{"11"} = "Evap Permeation";
$processName{"12"} = "Evap Fuel Vapor Venting"; 
$processName{"13"} = "Evap Fuel Leaks"; 
$processName{"15"} = "Crankcase Running Exhaust";
$processName{"16"} = "Crankcase Start Exhaust"; 
$processName{"17"} = "Crankcase Extended Idle Exhaust"; 
#$processName{"18"} = "Refueling Displacement Vapor Loss";   #exclude refueling
#$processName{"19"} = "Refueling Spillage Loss";   #exclude refueling
$processName{"90"} = "Extended Idle Exhaust";

# Process Types ------------------------------------------------------------------
# Table to map pollutant-process IDs to the output runspec files by EF table type
# Pollutant.Process reference - process is the last two characters of the code
# Binary: 100 = rate/distance   10 = rate/vehicle        1 = rate/profile
#

$PollProc_tablemap{"101"} = "100";
$PollProc_tablemap{"102"} = "010";
$PollProc_tablemap{"111"} = "110";
$PollProc_tablemap{"112"} = "101";
$PollProc_tablemap{"113"} = "110";
$PollProc_tablemap{"115"} = "100";
$PollProc_tablemap{"116"} = "010";
$PollProc_tablemap{"117"} = "010";
#$PollProc_tablemap{"118"} = "100";  #exclude refueling
#$PollProc_tablemap{"119"} = "100";  #exclude refueling
$PollProc_tablemap{"190"} = "010";
$PollProc_tablemap{"201"} = "100";
$PollProc_tablemap{"202"} = "010";
$PollProc_tablemap{"215"} = "100";
$PollProc_tablemap{"216"} = "010";
$PollProc_tablemap{"217"} = "010";
$PollProc_tablemap{"290"} = "010";
$PollProc_tablemap{"301"} = "100";
$PollProc_tablemap{"302"} = "010";
$PollProc_tablemap{"315"} = "100";
$PollProc_tablemap{"316"} = "010";
$PollProc_tablemap{"317"} = "010";
$PollProc_tablemap{"390"} = "010";
$PollProc_tablemap{"501"} = "100";
$PollProc_tablemap{"502"} = "010";
$PollProc_tablemap{"511"} = "110";
$PollProc_tablemap{"512"} = "101";
$PollProc_tablemap{"513"} = "110";
$PollProc_tablemap{"515"} = "100";
$PollProc_tablemap{"516"} = "010";
$PollProc_tablemap{"517"} = "010";
#$PollProc_tablemap{"518"} = "100";  #exclude refueling
#$PollProc_tablemap{"519"} = "100";  #exclude refueling
$PollProc_tablemap{"590"} = "010";
$PollProc_tablemap{"601"} = "100";
$PollProc_tablemap{"602"} = "010";
$PollProc_tablemap{"615"} = "100";
$PollProc_tablemap{"616"} = "010";
$PollProc_tablemap{"2001"} = "100";
$PollProc_tablemap{"2002"} = "010";
$PollProc_tablemap{"2011"} = "110";
$PollProc_tablemap{"2012"} = "101";
$PollProc_tablemap{"2013"} = "110";
$PollProc_tablemap{"2015"} = "100";
$PollProc_tablemap{"2016"} = "010";
$PollProc_tablemap{"2017"} = "010";
#$PollProc_tablemap{"2018"} = "100";  #exclude refueling
#$PollProc_tablemap{"2019"} = "100";  #exclude refueling
$PollProc_tablemap{"2090"} = "010";
$PollProc_tablemap{"2101"} = "100";
$PollProc_tablemap{"2102"} = "010";
$PollProc_tablemap{"2111"} = "110";
$PollProc_tablemap{"2112"} = "101";
$PollProc_tablemap{"2113"} = "110";
$PollProc_tablemap{"2115"} = "100";
$PollProc_tablemap{"2116"} = "010";
$PollProc_tablemap{"2117"} = "010";
#$PollProc_tablemap{"2118"} = "100";  #exclude refueling
#$PollProc_tablemap{"2119"} = "100";  #exclude refueling
$PollProc_tablemap{"2190"} = "010";
$PollProc_tablemap{"2201"} = "100";
$PollProc_tablemap{"2202"} = "010";
$PollProc_tablemap{"2211"} = "110";
$PollProc_tablemap{"2212"} = "101";
$PollProc_tablemap{"2213"} = "110";
$PollProc_tablemap{"2215"} = "100";
$PollProc_tablemap{"2216"} = "010";
$PollProc_tablemap{"2217"} = "010";
#$PollProc_tablemap{"2218"} = "100";  #exclude refueling
#$PollProc_tablemap{"2219"} = "100";  #exclude refueling
$PollProc_tablemap{"2290"} = "010";
$PollProc_tablemap{"2301"} = "100";
$PollProc_tablemap{"2302"} = "010";
$PollProc_tablemap{"2311"} = "110";
$PollProc_tablemap{"2312"} = "101";
$PollProc_tablemap{"2313"} = "110";
$PollProc_tablemap{"2315"} = "100";
$PollProc_tablemap{"2316"} = "010";
$PollProc_tablemap{"2317"} = "010";
#$PollProc_tablemap{"2318"} = "100";  #exclude refueling
#$PollProc_tablemap{"2319"} = "100";  #exclude refueling
$PollProc_tablemap{"2390"} = "010";
$PollProc_tablemap{"2401"} = "100";
$PollProc_tablemap{"2402"} = "010";
$PollProc_tablemap{"2415"} = "100";
$PollProc_tablemap{"2416"} = "010";
$PollProc_tablemap{"2417"} = "010";
$PollProc_tablemap{"2490"} = "010";
$PollProc_tablemap{"2501"} = "100";
$PollProc_tablemap{"2502"} = "010";
$PollProc_tablemap{"2515"} = "100";
$PollProc_tablemap{"2516"} = "010";
$PollProc_tablemap{"2517"} = "010";
$PollProc_tablemap{"2590"} = "010";
$PollProc_tablemap{"2601"} = "100";
$PollProc_tablemap{"2602"} = "010";
$PollProc_tablemap{"2615"} = "100";
$PollProc_tablemap{"2616"} = "010";
$PollProc_tablemap{"2617"} = "010";
$PollProc_tablemap{"2690"} = "010";
$PollProc_tablemap{"2701"} = "100";
$PollProc_tablemap{"2702"} = "010";
$PollProc_tablemap{"2715"} = "100";
$PollProc_tablemap{"2716"} = "010";
$PollProc_tablemap{"2717"} = "010";
$PollProc_tablemap{"2790"} = "010";
$PollProc_tablemap{"3001"} = "100";
$PollProc_tablemap{"3002"} = "010";
$PollProc_tablemap{"3015"} = "100";
$PollProc_tablemap{"3016"} = "010";
$PollProc_tablemap{"3017"} = "010";
$PollProc_tablemap{"3090"} = "010";
$PollProc_tablemap{"3101"} = "100";
$PollProc_tablemap{"3102"} = "010";
$PollProc_tablemap{"3115"} = "100";
$PollProc_tablemap{"3116"} = "010";
$PollProc_tablemap{"3117"} = "010";
$PollProc_tablemap{"3190"} = "010";
$PollProc_tablemap{"3201"} = "100";
$PollProc_tablemap{"3202"} = "010";
$PollProc_tablemap{"3215"} = "100";
$PollProc_tablemap{"3216"} = "010";
$PollProc_tablemap{"3217"} = "010";
$PollProc_tablemap{"3290"} = "010";
$PollProc_tablemap{"3301"} = "100";
$PollProc_tablemap{"3302"} = "010";
$PollProc_tablemap{"3315"} = "100";
$PollProc_tablemap{"3316"} = "010";
$PollProc_tablemap{"3317"} = "010";
$PollProc_tablemap{"3390"} = "010";
$PollProc_tablemap{"7901"} = "100";
$PollProc_tablemap{"7902"} = "010";
$PollProc_tablemap{"7911"} = "110";
$PollProc_tablemap{"7912"} = "101";
$PollProc_tablemap{"7913"} = "110";
$PollProc_tablemap{"7915"} = "100";
$PollProc_tablemap{"7916"} = "010";
$PollProc_tablemap{"7917"} = "010";
#$PollProc_tablemap{"7918"} = "100";  #exclude refueling
#$PollProc_tablemap{"7919"} = "100";  #exclude refueling
$PollProc_tablemap{"7990"} = "010";
$PollProc_tablemap{"8001"} = "100";
$PollProc_tablemap{"8002"} = "010";
$PollProc_tablemap{"8011"} = "110";
$PollProc_tablemap{"8012"} = "101";
$PollProc_tablemap{"8013"} = "110";
$PollProc_tablemap{"8015"} = "100";
$PollProc_tablemap{"8016"} = "010";
$PollProc_tablemap{"8017"} = "010";
#$PollProc_tablemap{"8018"} = "100";  #exclude refueling
#$PollProc_tablemap{"8019"} = "100";  #exclude refueling
$PollProc_tablemap{"8090"} = "010";
$PollProc_tablemap{"8601"} = "100";
$PollProc_tablemap{"8602"} = "010";
$PollProc_tablemap{"8611"} = "110";
$PollProc_tablemap{"8612"} = "101";
$PollProc_tablemap{"8613"} = "110";
$PollProc_tablemap{"8615"} = "100";
$PollProc_tablemap{"8616"} = "010";
$PollProc_tablemap{"8617"} = "010";
#$PollProc_tablemap{"8618"} = "100";  #exclude refueling
#$PollProc_tablemap{"8619"} = "100";  #exclude refueling
$PollProc_tablemap{"8690"} = "010";
$PollProc_tablemap{"8701"} = "100";
$PollProc_tablemap{"8702"} = "010";
$PollProc_tablemap{"8711"} = "110";
$PollProc_tablemap{"8712"} = "101";
$PollProc_tablemap{"8713"} = "110";
$PollProc_tablemap{"8715"} = "100";
$PollProc_tablemap{"8716"} = "010";
$PollProc_tablemap{"8717"} = "010";
#$PollProc_tablemap{"8718"} = "100";  #exclude refueling
#$PollProc_tablemap{"8719"} = "100";  #exclude refueling
$PollProc_tablemap{"8790"} = "010";
$PollProc_tablemap{"9001"} = "100";
$PollProc_tablemap{"9002"} = "010";
$PollProc_tablemap{"9090"} = "010";
$PollProc_tablemap{"9101"} = "100";
$PollProc_tablemap{"9102"} = "010";
$PollProc_tablemap{"9190"} = "010";
$PollProc_tablemap{"9201"} = "100";
$PollProc_tablemap{"9202"} = "010";
$PollProc_tablemap{"9290"} = "010";
$PollProc_tablemap{"9301"} = "100";
$PollProc_tablemap{"9302"} = "010";
$PollProc_tablemap{"9390"} = "010";
$PollProc_tablemap{"9801"} = "100";
$PollProc_tablemap{"9802"} = "010";
$PollProc_tablemap{"9890"} = "010";
$PollProc_tablemap{"10001"} = "100";
$PollProc_tablemap{"10002"} = "010";
$PollProc_tablemap{"10015"} = "100";
$PollProc_tablemap{"10016"} = "010";
$PollProc_tablemap{"10017"} = "010";
$PollProc_tablemap{"10090"} = "010";
$PollProc_tablemap{"10101"} = "100";
$PollProc_tablemap{"10102"} = "010";
$PollProc_tablemap{"10115"} = "100";
$PollProc_tablemap{"10116"} = "010";
$PollProc_tablemap{"10117"} = "010";
$PollProc_tablemap{"10190"} = "010";
$PollProc_tablemap{"10201"} = "100";
$PollProc_tablemap{"10202"} = "010";
$PollProc_tablemap{"10215"} = "100";
$PollProc_tablemap{"10216"} = "010";
$PollProc_tablemap{"10217"} = "010";
$PollProc_tablemap{"10290"} = "010";
$PollProc_tablemap{"10501"} = "100";
$PollProc_tablemap{"10502"} = "010";
$PollProc_tablemap{"10515"} = "100";
$PollProc_tablemap{"10516"} = "010";
$PollProc_tablemap{"10517"} = "010";
$PollProc_tablemap{"10590"} = "010";
$PollProc_tablemap{"10609"} = "100";
$PollProc_tablemap{"10710"} = "100";
$PollProc_tablemap{"11001"} = "100";
$PollProc_tablemap{"11002"} = "010";
$PollProc_tablemap{"11015"} = "100";
$PollProc_tablemap{"11016"} = "010";
$PollProc_tablemap{"11017"} = "010";
$PollProc_tablemap{"11090"} = "010";
$PollProc_tablemap{"11101"} = "100";
$PollProc_tablemap{"11102"} = "010";
$PollProc_tablemap{"11115"} = "100";
$PollProc_tablemap{"11116"} = "010";
$PollProc_tablemap{"11117"} = "010";
$PollProc_tablemap{"11190"} = "010";
$PollProc_tablemap{"11201"} = "100";
$PollProc_tablemap{"11202"} = "010";
$PollProc_tablemap{"11215"} = "100";
$PollProc_tablemap{"11216"} = "010";
$PollProc_tablemap{"11217"} = "010";
$PollProc_tablemap{"11290"} = "010";
$PollProc_tablemap{"11501"} = "100";
$PollProc_tablemap{"11502"} = "010";
$PollProc_tablemap{"11515"} = "100";
$PollProc_tablemap{"11516"} = "010";
$PollProc_tablemap{"11517"} = "010";
$PollProc_tablemap{"11590"} = "010";
$PollProc_tablemap{"11609"} = "100";
$PollProc_tablemap{"11710"} = "100";

#=========================================================================================================
# Read  RUN CONTROL file
#=========================================================================================================
#
# Example of run control file.  Format: KEYWORD = input parameters
#
#          DBHOST	= hostname   
#          BATCHRUN	= CENRAP
#          OUTDIR       = C:\Program Files\MOVES20091214\runspec_files\tests\
#          MODELYEAR	= 2005
#          POLLUTANTS 	= OZONE, TOXICS, PM, GHG
#          DAYOFWEEK	= WEEKDAY, WEEKEND
#          METFILE	= c:\movesdata\cenrap\2005\2005_repcounty_met.in

# open run control file
open(CONTROLFILE, "$RunControlFile") or die "Unable to open run Control file: $RunControlFile\n";


while (<CONTROLFILE>)
{
    chomp;
    $line = trim($_);

    if (($line eq "") || (substr($line, 0, 1) eq "#"))
    {
        next;
    }

    @line = split(/=/, $_);
    NXTLINE:
    {
        if (uc trim($line[0]) eq "DBHOST")           { $dbhost           = trim($line[1]); last NXTLINE; }
        if (uc trim($line[0]) eq "BATCHRUN")         { $batchrun         = trim($line[1]); last NXTLINE; }
        if (uc trim($line[0]) eq "OUTDIR")           { $outdir           = trim($line[1]); last NXTLINE; }
        if (uc trim($line[0]) eq "MOVESHOME")        { $moveshome        = trim($line[1]); last NXTLINE; }
        if (uc trim($line[0]) eq "MODELYEAR")        { $modelyear        = trim($line[1]); last NXTLINE; }
        if (uc trim($line[0]) eq "POLLUTANTS")       { $User_polls       = uc trim($line[1]); last NXTLINE; }
        if (uc trim($line[0]) eq "DAYOFWEEK")        { $dayofweek        = uc trim($line[1]); last NXTLINE; }
        if (uc trim($line[0]) eq "METFILE")          { $MetFile          = trim($line[1]); last NXTLINE; }
    }
}

# close the file
close CONTROLFILE;

# split up the multiple-value inputs -------------
@User_polls = split(/,/, $User_polls);
@dayofweek  = split(/,/, $dayofweek);

# Verify user run control file contains valid parameters ----------------

$moveshome = "C:\\Program Files\\MOVES20091214" if ($moveshome eq '');
printf "MOVES home directory is - %s\n",$moveshome;

die "ERROR - invalid value for MODELYEAR ('$modelyear'). Valid values are 1990 and 1999 - 2050 inclusive."
if !(($modelyear == 1990) || ($modelyear >= 1999 && $modelyear <= 2050));

for($i=0;$i<=$#User_polls;++$i) {
	for($j=0;$j<=$#pollOptions;++$j) {
		goto NXTPOLL if ( uc trim($User_polls[$i]) eq $pollOptions[$j] );
	}
	die "ERROR - invalid value for POLLUTANTS ('$User_polls[$i]'). Valid values are 'OZONE','TOXICS','HC','PM','GHG'.";
NXTPOLL:
}
$pollsFlg = "_";
for($i=0;$i<=$#User_polls;++$i) {
	$pollsFlg = $pollsFlg . (uc trim ($User_polls[$i]));
}
for($i=0;$i<=$#dayofweek;++$i) {
	die "ERROR - invalid value for DAYOFWEEK ('$dayofweek[$i]'). Valid values are 'WEEKDAY' and 'WEEKEND'."
	if (trim($dayofweek[$i]) ne 'WEEKDAY' && trim($dayofweek[$i]) ne 'WEEKEND');
	$WeekDayFlag = 1 if (trim($dayofweek[$i]) eq 'WEEKDAY');
	$WeekEndFlag = 1 if (trim($dayofweek[$i]) eq 'WEEKEND');
}
#  if day of week not specified, assume both weekday and weekend day run
if ( $WeekDayFlag == 0 && $WeekEndFlag == 0 ) {
	$WeekDayFlag = 1;
	$WeekEndFlag = 1;
}

# check last character of output directory name - it must include the slash
die "ERROR - invalid pathname OUTDIR ('$outdir'). Directory name must end in a slash."
if !( (substr($outdir,-1,1) eq "/") || (substr($outdir,-1,1) eq "\\") );

# end of run control file reading and parsing
 
#=========================================================================================================
#   open the batch files for list of runspec and data importer filenames
#=========================================================================================================
my $olen = length($outdir); 
my $slash = substr($outdir,-1,1);

$batchFile = $outdir . $batchrun . "_" . $modelyear . "runspec.";
$batchFile .= $GenCSH ? "csh" : "bat";
open(BATFILE, ">$batchFile") or die "Unable to open batch file for list of run spec filenames: $batchFile\n";
printf BATFILE "#!/bin/csh -xf\n" if $GenCSH;
printf BATFILE "cd \"%s\"\n",$moveshome;
printf BATFILE "source setenv.csh\n" if $GenCSH;
printf BATFILE "call setenv.bat\n" unless $GenCSH;
printf BATFILE "\@echo on\n" unless $GenCSH;
printf BATFILE "echo \"\"" if $GenCSH;
printf BATFILE "type null" unless $GenCSH;
printf BATFILE " > \"%s\"%s%s\n", substr($outdir,0,$olen-1),$slash,"runlog_".$batchrun."_".$modelyear.".txt";

$importBatchFile = $outdir . "/" . $batchrun . "_" . $modelyear . "importer.";
$importBatchFile .= $GenCSH ? "csh" : "bat";
open(IMPFILE, ">$importBatchFile") or die "Unable to open batch file for list of data importer filenames: $importBatchFile\n";
printf IMPFILE "#!/bin/csh -xf\n" if $GenCSH;
printf IMPFILE "cd \"%s\"\n",$moveshome;
printf IMPFILE "source setenv.csh\n" if $GenCSH;
printf IMPFILE "call setenv.bat\n" unless $GenCSH;
printf IMPFILE "\@echo on\n" unless $GenCSH;
printf IMPFILE "echo \"\"" if $GenCSH;
printf IMPFILE "type null" unless $GenCSH;
printf IMPFILE " > \"%s\"%s%s\n", substr($outdir,0,$olen-1),$slash,"importlog_".$batchrun."_".$modelyear.".txt";
 
#=========================================================================================================
# Generate complete list of pollutants and processes 

for($i=0;$i<=$#User_polls;++$i) {
	for($j=0;$j<=$#pollOptions;++$j) {
		&setPollsList(@pollsByOptionOZONE)  if ( uc trim($User_polls[$i]) eq "OZONE" );
		&setPollsList(@pollsByOptionTOXICS) if ( uc trim($User_polls[$i]) eq "TOXICS" );
		&setPollsList(@pollsByOptionPM)     if ( uc trim($User_polls[$i]) eq "PM" );
		&setPollsList(@pollsByOptionGHG)    if ( uc trim($User_polls[$i]) eq "GHG" );
	}
}

#qa
#for($i=0;$i<=$#pollsListID;++$i) {
#printf "%2.2d %d %d\n",$i,$pollsListID[$i], $pollsOutList[$i];
#}

#=========================================================================================================
# Read the RepCounty file
#=========================================================================================================

# open rep county input file

$cntRepCnty = 0;

open(REPFILE, "$RepCntyFile") or die "Unable to open representative county file: $RepCntyFile\n";

while (<REPFILE>)
{
    chomp;
    $line = trim($_);

    @line = split(/=/, $_);
    NXTREP:
    {
    if (($line eq "") || (substr($line, 0, 1) eq "#"))    { last NXTREP; }
    if (uc substr($line,0,10) eq "<REPCOUNTY")            { ++$cntRepCnty; last NXTREP; }
    if (uc trim($line[0]) eq "FIPS")                      { $repFips[$cntRepCnty] = trim($line[1]); last NXTREP; }
    if (uc trim($line[0]) eq "AGE")                       { $repAge[$cntRepCnty] = trim($line[1]); last NXTREP; }
    if (uc trim($line[0]) eq "IM")                        { $repIM[$cntRepCnty] = trim($line[1]); last NXTREP; }
    if (uc trim($line[0]) eq "FUELSUPPLY")                { $repFuelSup[$cntRepCnty] = trim($line[1]); last NXTREP; }
    if (uc trim($line[0]) eq "FUELFORM")                  { $repFuelForm[$cntRepCnty] = trim($line[1]); last NXTREP; }
    if (uc trim($line[0]) eq "POP")                       { $repPop[$cntRepCnty] = trim($line[1]); last NXTREP; }
    if (uc trim($line[0]) eq "HPMSVMT")                   { $repVMT[$cntRepCnty] = trim($line[1]); last NXTREP; }
    }
}
#qa printf "check repcounty in repc %s\n",$repFips[$cntRepCnty];

#=========================================================================================================
# Read the met4moves input met file and call routines to generate RunSpec files and DataImport files
#=========================================================================================================

# open the met data
open(METFILE, "$MetFile") or die "Unable to open met file: $MetFile\n";
my ($VV_TempInc, $RD_TempInc, $RV_TempInc, $temperature);
my ($RDminT, $RDmaxT, $RVminT, $RVmaxT);
my ($t, $numRDruns, $cntRDtemps, $cntRVtemps, @RDtemps, @RVtemps);
my ($cntMet, $MetRep, $MetMonth, $MetProfId, $MetRH,%fipsList);
my ($MetMin, $MetMax, @MetTempProf);
my ($outputDB, $scenarioID, @scenarioIDList);
my ($cntyidx,$flg,$ref,$pp,$process,$code,$tID,$codeOut,$pollRef,$cnty);

# --- read met data records ==============================
#
$cntMet = 0;
while (<METFILE>)
{
    chomp;
    $line = trim($_);
    @line = split(/\s+/, trim($_));

# read header record for temperature bin increments
    NXTMET:
    {
    if (($line eq "") || (substr($line, 0, 1) eq "#"))    { last NXTMET; }
    if (uc trim($line[0]) eq "PP_TEMP_INCREMENT")         { $VV_TempInc = trim($line[1]); last NXTMET;  }
    if (uc trim($line[0]) eq "PD_TEMP_INCREMENT")         { $RD_TempInc = trim($line[1]); last NXTMET;  }
    if (uc trim($line[0]) eq "PV_TEMP_INCREMENT")         { $RV_TempInc = trim($line[1]); last NXTMET;  }
    ++$cntMet;

    # split data fields, determine record types, and parse each appropriately
    # min/max record:  fips, fuelMonth, key, rh, tmin, tmax
    # diurnal record:  fips, fuelMonth, key, rh, temp1, temp2, ...temp24

    @line = split(/\s+/, $_);
    $MetRep = substr(trim($line[0]),1,5);
    $cntyidx = &getRepCnty();
	if ($cntyidx <= 0) {die "ERROR : Repcounty input file incomplete.  Missing county:  $MetRep\nREPCOUNTY packets must exist for all FIPS in input met file.\n";}
    $MetMonth = trim($line[1]);
    $MetProfId = trim($line[2]);
    $MetRH = $line[3];
    $outputDB = $MetRep . "_" . $modelyear . $pollsFlg;
    $fipsList{$MetRep} = 1;

    if ($MetProfId eq "min_max" )
    {
	$MetMin = $line[4];
	$MetMax = $line[5];

#       --- determine the temperature bins between min/max temperatures
        &setMinMax();
#
#       --- rate per distance runs --- ======================================================================
        for($t=1;$t<=$numRDruns;++$t)
	{


	   $scenarioID = "RD_".$MetRep . "_" . $modelyear."_".$MetMonth ."_T". int($RDtemps[$t][1])."_". int($RDtemps[$t][24]);
#	   --- write the meteorology MOVES input csv formatted file
	   $fileout = $outdir .  $scenarioID . "_zmh.csv";
           open (OUTFL,">$fileout") || die "Cannot open file: $fileout\n";
	   printf OUTFL "monthID,zoneID,hourID,temperature,relHumidity\n";

	   for($i=0;$i<=23;++$i) 
	   { 
              printf OUTFL "%d,%s0,%d,%5.1f,%5.1f\n",$MetMonth,$MetRep,$i+1, $RDtemps[$t][$i+1], $MetRH;
	   }
           close(OUTFL);

#          --- write the data importer for this runspec 
	   $fileout = $outdir . $scenarioID . "_imp.xml";
           open (OUTFL,">$fileout") || die "Cannot open file: $fileout\n";
	   printf IMPFILE "java gov.epa.otaq.moves.master.commandline.MOVESCommandLine -i \"%s\"%s%s", 
                                                                    substr($outdir,0,$olen-1),$slash,$scenarioID."_imp.xml";
	   printf IMPFILE " >> \"%s\"%s%s\n", substr($outdir,0,$olen-1),$slash,"importlog_".$batchrun."_".$modelyear.".txt";
	   RD_writeDataImporter();
	   close(OUTFL);

#          --- write the runspec file
	   $fileout = $outdir . "/".  $scenarioID . "_mrs.xml";
           open (OUTFL,">$fileout") || die "Cannot open file: $fileout\n";
	   printf BATFILE "java gov.epa.otaq.moves.master.commandline.MOVESCommandLine -r \"%s\"%s%s", 
                                                                    substr($outdir,0,$olen-1),$slash,$scenarioID."_mrs.xml";
	   printf BATFILE " >> \"%s\"%s%s\n", substr($outdir,0,$olen-1),$slash,"runlog_".$batchrun."_".$modelyear.".txt";
	   RD_writeRunSpec();
	   close(OUTFL);
	}

#        --- rate per vehicle runs --- ========================================================================
        for($t=1;$t<=$cntRVtemps;++$t)
	{
	   $scenarioID = "RV_" . $MetRep . "_" . $modelyear . "_" . $MetMonth . "_T" . int($RVtemps[$t]);

#	   --- write the meteorology MOVES input csv formatted file
	   $fileout = $outdir . "/".  $scenarioID . "_zmh.csv";
           open (OUTFL,">$fileout") || die "Cannot open file: $fileout\n";
	   printf OUTFL "monthID,zoneID,hourID,temperature,relHumidity\n";
	   for($i=0;$i<=23;++$i) 
	   { 
              printf OUTFL "%d,%s0,%d,%5.1f,%5.1f\n",$MetMonth,$MetRep,$i+1, $RVtemps[$t], $MetRH;
	   }
           close(OUTFL);

#          --- write the data importer for this runspec 
	   $fileout = $outdir . "/".  $scenarioID . "_imp.xml";
           open (OUTFL,">$fileout") || die "Cannot open file: $fileout\n";
	   printf IMPFILE "java gov.epa.otaq.moves.master.commandline.MOVESCommandLine -i \"%s\"%s%s", 
                                                                    substr($outdir,0,$olen-1),$slash,$scenarioID."_imp.xml";
	   printf IMPFILE " >> \"%s\"%s%s\n", substr($outdir,0,$olen-1),$slash,"importlog_".$batchrun."_".$modelyear.".txt";
	   RV_writeDataImporter();
	   close(OUTFL);

#          --- write the runspec file
	   $fileout = $outdir . "/".  $scenarioID . "_mrs.xml";
           open (OUTFL,">$fileout") || die "Cannot open file: $fileout\n";
	   printf BATFILE "java gov.epa.otaq.moves.master.commandline.MOVESCommandLine -r \"%s\"%s%s", 
                                                                    substr($outdir,0,$olen-1),$slash,$scenarioID."_mrs.xml";
	   printf BATFILE " >> \"%s\"%s%s\n", substr($outdir,0,$olen-1),$slash,"runlog_".$batchrun."_".$modelyear.".txt";
	   RV_writeRunSpec();
	   close(OUTFL);
	}
    }
    else
#   24-hour temperature profiles required for vapor venting emissions mode
    {
	$scenarioID = "RP_" . $MetRep . "_" . $modelyear . "_" . $MetMonth . "_prof" . $MetProfId;

#	--- write the meteorology MOVES input csv formatted file
	$fileout = $outdir . "/".  $scenarioID . "_zmh.csv";
        open (OUTFL,">$fileout") || die "Cannot open file: $fileout\n";
	   printf OUTFL "monthID,zoneID,hourID,temperature,relHumidity\n";
	for($i=0;$i<=23;++$i) 
	{ 
           printf OUTFL "%d,%s0,%d,%5.1f,%5.1f\n",$MetMonth,$MetRep,$i+1, $line[4+$i], $MetRH;
	}
        close(OUTFL);

#       --- write the data importer for this runspec 
	$fileout = $outdir . "/".  $scenarioID . "_imp.xml";
        open (OUTFL,">$fileout") || die "Cannot open file: $fileout\n";
	printf IMPFILE "java gov.epa.otaq.moves.master.commandline.MOVESCommandLine -i \"%s\"%s%s", 
                                                                 substr($outdir,0,$olen-1),$slash,$scenarioID."_imp.xml";
	printf IMPFILE " >> \"%s\"%s%s\n", substr($outdir,0,$olen-1),$slash,"importlog_".$batchrun."_".$modelyear.".txt";
	VV_writeDataImporter(); 
	close(OUTFL);

#       --- write the runspec file
	$fileout = $outdir . "/".  $scenarioID."_mrs.xml";
        open (OUTFL,">$fileout") || die "Cannot open file: $fileout\n";
	printf BATFILE "java gov.epa.otaq.moves.master.commandline.MOVESCommandLine -r \"%s\"%s%s", 
                                                                 substr($outdir,0,$olen-1),$slash,$scenarioID."_mrs.xml";
	printf BATFILE " >> \"%s\"%s%s\n", substr($outdir,0,$olen-1),$slash,"runlog_".$batchrun."_".$modelyear.".txt";
	VV_writeRunSpec();
	close(OUTFL);

    }  # end if for profileid types
    }  # end of NXTMET

}  # end read met file

$dbListFile = $outdir . "/" . $batchrun . "_" . $modelyear . "outputDBs.txt";
open(DBFILE, ">$dbListFile") or die "Unable to open file for list of output DB names: $dbListFile\n";

printf DBFILE "%s\n",$dbhost;
printf DBFILE "%s\n",$outdir;
foreach $cnty (sort(keys(%fipsList)))
{
	printf DBFILE "%5.5d_%4.4d%s\n",$cnty, $modelyear, $pollsFlg;
}
#=========================================================================================================
#               SUBROUTINES
#=========================================================================================================

#=========================================================================================================
# Set the output pollutant list flag for each User pollutant option
#=========================================================================================================
sub setPollsList()
{
	foreach $_ (@_)
		{
			for($ip=0;$ip<=$#pollsListID;++$ip) {
				$pollsOutList[$ip] = 1 if ($pollsListID[$ip] == $_);
			}
		}
}  #end subroutine setPollsList

#=========================================================================================================
#  Determine the repcounty file county index
#=========================================================================================================
sub getRepCnty
{
   for($ic=1;$ic<=$cntRepCnty;++$ic) 
   {
	if ($repFips[$ic] eq $MetRep)
	{
	   return $ic;
	}
   }

}  # end subroutine getRepCnty

#=========================================================================================================
#  Set the min/max temperature bins for each runspec file
#=========================================================================================================
sub setMinMax
{
   $numRDruns = 1;
   $cntRDtemps = 0;
   $RDminT = (int($MetMin/$RD_TempInc) ) * $RD_TempInc;
   $RDminT = $RDminT - $RD_TempInc if ($MetMin <= 0.0);
   $RDmaxT = (int($MetMax/$RD_TempInc) ) * $RD_TempInc + $RD_TempInc;
   $RDmaxT = $RDmaxT - $RD_TempInc if ($MetMax <= 0.0);
   $temperature = $RDminT;
   while ($temperature <= $RDmaxT)
	{
	   ++$cntRDtemps;
	   if ( $cntRDtemps > 24 ) {$cntRDtemps = 1; ++$numRDruns;}
	   $RDtemps[$numRDruns][$cntRDtemps] = $temperature;
	   $temperature = $temperature + $RD_TempInc;
	}
#  if the number of temps < 24 then pad with last temp through hour 24
   if ($cntRDtemps < 24)
   {
	   for($i=$cntRDtemps+1;$i<=24;++$i){
		$RDtemps[$numRDruns][$i] = $temperature - $RD_TempInc; }
   }

   $cntRVtemps = 0;
   $RVminT = (int($MetMin/$RV_TempInc) ) * $RV_TempInc;
   $RVminT = $RVminT - $RV_TempInc if ($MetMin <= 0.0);
   $RVmaxT = (int($MetMax/$RV_TempInc) ) * $RV_TempInc + $RV_TempInc;
   $RVmaxT = $RVmaxT - $RV_TempInc if ($MetMax <= 0.0);
   $temperature = $RVminT;
   while ($temperature <= $RVmaxT)
	{
	   ++$cntRVtemps;
	   $RVtemps[$cntRVtemps] = $temperature;
	   $temperature = $temperature + $RV_TempInc;
	}

} # end subroutine setMinMax

#=========================================================================================================
#  Generate the RunSpec files for categoryA; on-network operating mode
#=========================================================================================================
sub RD_writeRunSpec
{
   printf OUTFL "\t<runspec>\n";
   printf OUTFL "\t<description><![CDATA[RunSpec Generator for MOVES2010 - %s]]></description>\n",$scenarioID;
   printf OUTFL "\t<modelscale value=\"Rates\"/>\n";
   printf OUTFL "\t<modeldomain value=\"SINGLE\"/>\n";

   &geoselect();
   &timespan(1);
   &vehsel();

   printf OUTFL "\t<roadtypes>\n";
   printf OUTFL "\t\t<roadtype roadtypeid=\"1\" roadtypename=\"Off-Network\"/>\n";
   printf OUTFL "\t\t<roadtype roadtypeid=\"2\" roadtypename=\"Rural Restricted Access\"/>\n";
   printf OUTFL "\t\t<roadtype roadtypeid=\"3\" roadtypename=\"Rural Unrestricted Access\"/>\n";
   printf OUTFL "\t\t<roadtype roadtypeid=\"4\" roadtypename=\"Urban Restricted Access\"/>\n";
   printf OUTFL "\t\t<roadtype roadtypeid=\"5\" roadtypename=\"Urban Unrestricted Access\"/>\n";
   printf OUTFL "\t</roadtypes>\n";

   &pollProc (1);   # pass the column of interest for this runspec type

   &rspend();
}  # end subroutine RD_writeRunSpec

#=========================================================================================================
#  Generate the Data Importer files for categoryA; on-network operating mode
#=========================================================================================================
sub RD_writeDataImporter
{
   printf OUTFL "\t<moves>\n";
   printf OUTFL "\t\t<importer mode=\"county\">\n";
   printf OUTFL "\t\t<filters>\n";

   &geoselect();
   &timespan(1);
   &vehsel();

   printf OUTFL "\t<roadtypes>\n";
   printf OUTFL "\t\t<roadtype roadtypeid=\"1\" roadtypename=\"Off-Network\"/>\n";
   printf OUTFL "\t\t<roadtype roadtypeid=\"2\" roadtypename=\"Rural Restricted Access\"/>\n";
   printf OUTFL "\t\t<roadtype roadtypeid=\"3\" roadtypename=\"Rural Unrestricted Access\"/>\n";
   printf OUTFL "\t\t<roadtype roadtypeid=\"4\" roadtypename=\"Urban Restricted Access\"/>\n";
   printf OUTFL "\t\t<roadtype roadtypeid=\"5\" roadtypename=\"Urban Unrestricted Access\"/>\n";
   printf OUTFL "\t</roadtypes>\n";

   &pollProc (1);   # pass the column of interest for this runspec type
   printf OUTFL "\t\t</filters>\n";

   &dataImporter();
}  # end subroutine RD_writeDataImporter

#=========================================================================================================
#  Generate the RunSpec files for categoryB; offnetwork processes
#=========================================================================================================
sub RV_writeRunSpec
{
   printf OUTFL "\t<runspec>\n";
   printf OUTFL "\t<description><![CDATA[RunSpec Generator for MOVES2010 - %s]]></description>\n",$scenarioID;
   printf OUTFL "\t<modelscale value=\"Rates\"/>\n";
   printf OUTFL "\t<modeldomain value=\"SINGLE\"/>\n";

   &geoselect();
   &timespan(0);
   &vehsel();

   printf OUTFL "\t<roadtypes>\n";
   printf OUTFL "\t\t<roadtype roadtypeid=\"1\" roadtypename=\"Off-Network\"/>\n";
   printf OUTFL "\t\t<roadtype roadtypeid=\"2\" roadtypename=\"Rural Restricted Access\"/>\n";
   printf OUTFL "\t\t<roadtype roadtypeid=\"3\" roadtypename=\"Rural Unrestricted Access\"/>\n";
   printf OUTFL "\t\t<roadtype roadtypeid=\"4\" roadtypename=\"Urban Restricted Access\"/>\n";
   printf OUTFL "\t\t<roadtype roadtypeid=\"5\" roadtypename=\"Urban Unrestricted Access\"/>\n";
   printf OUTFL "\t</roadtypes>\n";

   &pollProc (2);   # pass the column of interest for this runspec type

   &rspend();
}  # end subroutine RV_writeRunSpec

#=========================================================================================================
#  Generate the Data Importer files for categoryB; offnetwork processes
#=========================================================================================================
sub RV_writeDataImporter
{
   printf OUTFL "\t<moves>\n";
   printf OUTFL "\t\t<importer mode=\"county\">\n";
   printf OUTFL "\t\t<filters>\n";

   &geoselect();
   &timespan(0);
   &vehsel();

   printf OUTFL "\t<roadtypes>\n";
   printf OUTFL "\t\t<roadtype roadtypeid=\"1\" roadtypename=\"Off-Network\"/>\n";
   printf OUTFL "\t\t<roadtype roadtypeid=\"2\" roadtypename=\"Rural Restricted Access\"/>\n";
   printf OUTFL "\t\t<roadtype roadtypeid=\"3\" roadtypename=\"Rural Unrestricted Access\"/>\n";
   printf OUTFL "\t\t<roadtype roadtypeid=\"4\" roadtypename=\"Urban Restricted Access\"/>\n";
   printf OUTFL "\t\t<roadtype roadtypeid=\"5\" roadtypename=\"Urban Unrestricted Access\"/>\n";
   printf OUTFL "\t</roadtypes>\n";

   &pollProc (2);   # pass the column of interest for this runspec type
   printf OUTFL "\t\t</filters>\n";

   &dataImporter();
}  # end subroutine RV_writeDataImporter

#=========================================================================================================
#  Generate the RunSpec files for categoryC; vapor venting
#=========================================================================================================
sub VV_writeRunSpec
{
   printf OUTFL "\t<runspec>\n";
   printf OUTFL "\t<description><![CDATA[RunSpec Generator for MOVES2010 - %s]]></description>\n",$scenarioID;
   printf OUTFL "\t<modelscale value=\"Rates\"/>\n";
   printf OUTFL "\t<modeldomain value=\"SINGLE\"/>\n";

   &geoselect();
   &timespan(0);
   &vehsel();

   printf OUTFL "\t<roadtypes>\n";
   printf OUTFL "\t\t<roadtype roadtypeid=\"1\" roadtypename=\"Off-Network\"/>\n";
   printf OUTFL "\t\t<roadtype roadtypeid=\"2\" roadtypename=\"Rural Restricted Access\"/>\n";
   printf OUTFL "\t\t<roadtype roadtypeid=\"3\" roadtypename=\"Rural Unrestricted Access\"/>\n";
   printf OUTFL "\t\t<roadtype roadtypeid=\"4\" roadtypename=\"Urban Restricted Access\"/>\n";
   printf OUTFL "\t\t<roadtype roadtypeid=\"5\" roadtypename=\"Urban Unrestricted Access\"/>\n";
   printf OUTFL "\t</roadtypes>\n";

   &pollProc (3);   # pass the column of interest for this runspec type

   &rspend();
}  # end subroutine VV_writeRunSpec

#=========================================================================================================
#  Generate the Data Importer files for categoryC; vapor venting
#=========================================================================================================
sub VV_writeDataImporter
{
   printf OUTFL "\t<moves>\n";
   printf OUTFL "\t\t<importer mode=\"county\">\n";
   printf OUTFL "\t\t<filters>\n";

   &geoselect();
   &timespan(0);
   &vehsel();

   printf OUTFL "\t<roadtypes>\n";
   printf OUTFL "\t\t<roadtype roadtypeid=\"1\" roadtypename=\"Off-Network\"/>\n";
   printf OUTFL "\t\t<roadtype roadtypeid=\"2\" roadtypename=\"Rural Restricted Access\"/>\n";
   printf OUTFL "\t\t<roadtype roadtypeid=\"3\" roadtypename=\"Rural Unrestricted Access\"/>\n";
   printf OUTFL "\t\t<roadtype roadtypeid=\"4\" roadtypename=\"Urban Restricted Access\"/>\n";
   printf OUTFL "\t\t<roadtype roadtypeid=\"5\" roadtypename=\"Urban Unrestricted Access\"/>\n";
   printf OUTFL "\t</roadtypes>\n";

   &pollProc (3);   # pass the column of interest for this runspec type
   printf OUTFL "\t\t</filters>\n";

   &dataImporter();
}  # end subroutine VV_writeDataImporter

# --- geographic selections - ========================================================================
sub geoselect
{
   printf OUTFL "\t<geographicselections>\n";
   printf OUTFL "\t\t<geographicselection type=\"COUNTY\" key=\"%d\" description=\"\"/>", $MetRep;
   printf OUTFL "\t</geographicselections>\n";
        
}  # end geoselect subroutine

# --- timespan - ====================================================================================
sub timespan
{
   $flg = $_[0];

   printf OUTFL "\t<timespan>\n";
   printf OUTFL "\t\t<year key=\"%4d\"/>\n", $modelyear;
   printf OUTFL "\t\t<month id=\"%d\"/>\n", $MetMonth;
   if ( $flg == 1 ) {
       printf OUTFL "\t\t<day id=\"5\"/>\n"; }
   else {
       printf OUTFL "\t\t<day id=\"5\"/>\n" if ($WeekDayFlag);
       printf OUTFL "\t\t<day id=\"2\"/>\n" if ($WeekEndFlag);
   }
   printf OUTFL "\t\t<beginhour id=\"1\"/>\n";
   printf OUTFL "\t\t<endhour id=\"24\"/>\n";
   printf OUTFL "\t\t<aggregateBy key=\"Hour\"/>\n";
   printf OUTFL "\t</timespan>\n";
} # end timespan subroutine

# --- onroad vehicle selections - ====================================================================
sub vehsel
{
   #--Note - gas intercity bus, gas combination long-haul truck, diesel motorcycles not supported in MOVES2010 DB
   printf OUTFL "\t<onroadvehicleselections>\n";
   printf OUTFL "\t\t<onroadvehicleselection fueltypeid=\"1\" fueltypedesc=\"Gasoline\" sourcetypeid=\"11\" sourcetypename=\"Motorcycle\"/>\n";
   printf OUTFL "\t\t<onroadvehicleselection fueltypeid=\"1\" fueltypedesc=\"Gasoline\" sourcetypeid=\"21\" sourcetypename=\"Passenger Car\"/>\n";
   printf OUTFL "\t\t<onroadvehicleselection fueltypeid=\"1\" fueltypedesc=\"Gasoline\" sourcetypeid=\"31\" sourcetypename=\"Passenger Truck\"/>\n";
   printf OUTFL "\t\t<onroadvehicleselection fueltypeid=\"1\" fueltypedesc=\"Gasoline\" sourcetypeid=\"32\" sourcetypename=\"Light Commercial Truck\"/>\n";
   printf OUTFL "\t\t<onroadvehicleselection fueltypeid=\"1\" fueltypedesc=\"Gasoline\" sourcetypeid=\"42\" sourcetypename=\"Transit Bus\"/>\n";
   printf OUTFL "\t\t<onroadvehicleselection fueltypeid=\"1\" fueltypedesc=\"Gasoline\" sourcetypeid=\"43\" sourcetypename=\"School Bus\"/>\n";
   printf OUTFL "\t\t<onroadvehicleselection fueltypeid=\"1\" fueltypedesc=\"Gasoline\" sourcetypeid=\"51\" sourcetypename=\"Refuse Truck\"/>\n";
   printf OUTFL "\t\t<onroadvehicleselection fueltypeid=\"1\" fueltypedesc=\"Gasoline\" sourcetypeid=\"52\" sourcetypename=\"Single Unit Short-haul Truck\"/>\n";
   printf OUTFL "\t\t<onroadvehicleselection fueltypeid=\"1\" fueltypedesc=\"Gasoline\" sourcetypeid=\"53\" sourcetypename=\"Single Unit Long-haul Truck\"/>\n";
   printf OUTFL "\t\t<onroadvehicleselection fueltypeid=\"1\" fueltypedesc=\"Gasoline\" sourcetypeid=\"54\" sourcetypename=\"Motor Home\"/>\n";
   printf OUTFL "\t\t<onroadvehicleselection fueltypeid=\"1\" fueltypedesc=\"Gasoline\" sourcetypeid=\"61\" sourcetypename=\"Combination Short-haul Truck\"/>\n";
   printf OUTFL "\t\t<onroadvehicleselection fueltypeid=\"2\" fueltypedesc=\"Diesel Fuel\" sourcetypeid=\"21\" sourcetypename=\"Passenger Car\"/> \n";
   printf OUTFL "\t\t<onroadvehicleselection fueltypeid=\"2\" fueltypedesc=\"Diesel Fuel\" sourcetypeid=\"31\" sourcetypename=\"Passenger Truck\" />\n";
   printf OUTFL "\t\t<onroadvehicleselection fueltypeid=\"2\" fueltypedesc=\"Diesel Fuel\" sourcetypeid=\"32\" sourcetypename=\"Light Commercial Truck\"/>\n";
   printf OUTFL "\t\t<onroadvehicleselection fueltypeid=\"2\" fueltypedesc=\"Diesel Fuel\" sourcetypeid=\"41\" sourcetypename=\"Intercity Bus\"/> \n";
   printf OUTFL "\t\t<onroadvehicleselection fueltypeid=\"2\" fueltypedesc=\"Diesel Fuel\" sourcetypeid=\"42\" sourcetypename=\"Transit Bus\"/>\n";
   printf OUTFL "\t\t<onroadvehicleselection fueltypeid=\"2\" fueltypedesc=\"Diesel Fuel\" sourcetypeid=\"43\" sourcetypename=\"School Bus\"/>\n";
   printf OUTFL "\t\t<onroadvehicleselection fueltypeid=\"2\" fueltypedesc=\"Diesel Fuel\" sourcetypeid=\"51\" sourcetypename=\"Refuse Truck\"/>\n";
   printf OUTFL "\t\t<onroadvehicleselection fueltypeid=\"2\" fueltypedesc=\"Diesel Fuel\" sourcetypeid=\"52\" sourcetypename=\"Single Unit Short-haul Truck\"/>\n";
   printf OUTFL "\t\t<onroadvehicleselection fueltypeid=\"2\" fueltypedesc=\"Diesel Fuel\" sourcetypeid=\"53\" sourcetypename=\"Single Unit Long-haul Truck\"/>\n";
   printf OUTFL "\t\t<onroadvehicleselection fueltypeid=\"2\" fueltypedesc=\"Diesel Fuel\" sourcetypeid=\"54\" sourcetypename=\"Motor Home\"/>\n";
   printf OUTFL "\t\t<onroadvehicleselection fueltypeid=\"2\" fueltypedesc=\"Diesel Fuel\" sourcetypeid=\"61\" sourcetypename=\"Combination Short-haul Truck\"/>\n";
   printf OUTFL "\t\t<onroadvehicleselection fueltypeid=\"2\" fueltypedesc=\"Diesel Fuel\" sourcetypeid=\"62\" sourcetypename=\"Combination Long-haul Truck\"/>\n";
   printf OUTFL "\t</onroadvehicleselections>\n";
}  # end vehsel subroutine

# --- pollutant process associations - ====================================================================
sub pollProc
{
   $ref = $_[0];
   $ref = $ref - 1;

#qa printf "reference is %d\n", $ref;
my ($qac);
$qac = 0;
   printf OUTFL "\t<pollutantprocessassociations>\n";
   foreach $pp (keys %PollProc_tablemap)
   {
	++$qac;
	$process = substr($pp,-2);
	$code = $PollProc_tablemap{$pp};
	$codeOut = substr($code,$ref,1);
#qa printf "pp %s process %s code %s codeout %s\n", $pp, $process,$code,$codeOut if ($qac < 10);
	if ($codeOut eq "1") 
	{
		$pollRef = substr($pp,0,length($pp)-2);

		for($ip=0;$ip<=$#pollsListID;++$ip) {
			if ( $pollRef eq $pollsListID[$ip] && $pollsOutList[$ip] == 1)
			{
			#  print to output file
				printf OUTFL "\t\t<pollutantprocessassociation pollutantkey=\"%d\" ", $pollsListID[$ip];
				printf OUTFL "pollutantname=\"%s\" ", $pollsListName[$ip];
				printf OUTFL "processkey=\"%d\" ", $process;
				printf OUTFL "processname=\"%s\"/>\n", $processName{$process};
			}
		}

	} # end for this runspec table
   } # end each of pollutant - process mappings
   printf OUTFL "\t</pollutantprocessassociations>\n";
} # end pollProc subroutine


# --- finish up the data importer file - ====================================================================
sub dataImporter
{
   printf OUTFL "\t<databaseselection servername=\"localhost\" databasename=\"%s\"/>\n",$scenarioID."_in";

   printf OUTFL "\t<agedistribution>\n";
   printf OUTFL "\t\t<description><![CDATA[]]></description>\n";
   printf OUTFL "\t\t<parts>\n";
   printf OUTFL "\t\t\t<sourceTypeAgeDistribution>\n";
   printf OUTFL "\t\t\t<filename>%s</filename>\n", $repAge[$cntyidx];
   printf OUTFL "\t\t\t</sourceTypeAgeDistribution>\n";
   printf OUTFL "\t\t</parts>\n";
   printf OUTFL "\t</agedistribution>\n";

   printf OUTFL "\t<avgspeeddistribution>\n";
   printf OUTFL "\t\t<description><![CDATA[]]></description>\n";
   printf OUTFL "\t\t<parts>\n";
   printf OUTFL "\t\t\t<avgSpeedDistribution>\n";
   printf OUTFL "\t\t\t<filename>%sdummy_avgspeeddistribution.csv</filename>\n",$outdir;
   printf OUTFL "\t\t\t</avgSpeedDistribution>\n";
   printf OUTFL "\t\t</parts>\n";
   printf OUTFL "\t</avgspeeddistribution>\n";

   printf OUTFL "\t<fuelsupply>\n";
   printf OUTFL "\t\t<description><![CDATA[]]></description>\n";
   printf OUTFL "\t\t<parts>\n";
   printf OUTFL "\t\t\t<fuelSupply>\n";
   printf OUTFL "\t\t\t<filename>%s</filename>\n", $repFuelSup[$cntyidx];
   printf OUTFL "\t\t\t</fuelSupply>\n";
   printf OUTFL "\t\t</parts>\n";
   printf OUTFL "\t</fuelsupply>\n";

   printf OUTFL "\t<fuelformulation>\n";
   printf OUTFL "\t\t<description><![CDATA[]]></description>\n";
   printf OUTFL "\t\t<parts>\n";
   printf OUTFL "\t\t\t<fuelFormulation>\n";
   printf OUTFL "\t\t\t<filename>%s</filename>\n", $repFuelForm[$cntyidx];
   printf OUTFL "\t\t\t</fuelFormulation>\n";
   printf OUTFL "\t\t</parts>\n";
   printf OUTFL "\t</fuelformulation>\n";

   printf OUTFL "\t<zonemonthhour>\n";
   printf OUTFL "\t\t<description><![CDATA[]]></description>\n";
   printf OUTFL "\t\t<parts>\n";
   printf OUTFL "\t\t\t<zoneMonthHour>\n";
   printf OUTFL "\t\t\t<filename>%s</filename>\n", $outdir.$scenarioID."_zmh.csv";
   printf OUTFL "\t\t\t</zoneMonthHour>\n";
   printf OUTFL "\t\t</parts>\n";
   printf OUTFL "\t</zonemonthhour>\n";

   printf OUTFL "\t<rampfraction>\n";
   printf OUTFL "\t\t<description><![CDATA[]]></description>\n";
   printf OUTFL "\t\t<parts>\n";
   printf OUTFL "\t\t\t<roadType>\n";
   printf OUTFL "\t\t\t<filename></filename>\n";
   printf OUTFL "\t\t\t</roadType>\n";
   printf OUTFL "\t\t</parts>\n";
   printf OUTFL "\t</rampfraction>\n";

   printf OUTFL "\t<roadtypedistribution>\n";
   printf OUTFL "\t\t<description><![CDATA[]]></description>\n";
   printf OUTFL "\t\t<parts>\n";
   printf OUTFL "\t\t\t<roadTypeDistribution>\n";
   printf OUTFL "\t\t\t<filename>%sdummy_roadtypedistribution.csv</filename>\n",$outdir;
   printf OUTFL "\t\t\t</roadTypeDistribution>\n";
   printf OUTFL "\t\t</parts>\n";
   printf OUTFL "\t</roadtypedistribution>\n";

   printf OUTFL "\t<sourcetypepopulation>\n";
   printf OUTFL "\t\t<description><![CDATA[]]></description>\n";
   printf OUTFL "\t\t<parts>\n";
   printf OUTFL "\t\t\t<sourceTypeYear>\n";
   printf OUTFL "\t\t\t<filename>%s</filename>\n",$repPop[$cntyidx];
   printf OUTFL "\t\t\t</sourceTypeYear>\n";
   printf OUTFL "\t\t</parts>\n";
   printf OUTFL "\t</sourcetypepopulation>\n";

   printf OUTFL "\t<vehicletypevmt>\n";
   printf OUTFL "\t\t<description><![CDATA[]]></description>\n";
   printf OUTFL "\t\t<parts>\n";
   printf OUTFL "\t\t\t<HPMSVTypeYear>\n";
   printf OUTFL "\t\t\t<filename>%s</filename>\n",$repVMT[$cntyidx];
   printf OUTFL "\t\t\t</HPMSVTypeYear>\n";
   printf OUTFL "\t\t\t<monthVMTFraction>\n";
   printf OUTFL "\t\t\t<filename>%sdummy_monthvmtfraction.csv</filename>\n",$outdir;
   printf OUTFL "\t\t\t</monthVMTFraction>\n";
   printf OUTFL "\t\t\t<dayVMTFraction>\n";
   printf OUTFL "\t\t\t<filename>%sdummy_dayvmtfraction.csv</filename>\n",$outdir;
   printf OUTFL "\t\t\t</dayVMTFraction>\n";
   printf OUTFL "\t\t\t<hourVMTFraction>\n";
   printf OUTFL "\t\t\t<filename>%sdummy_hourvmtfraction.csv</filename>\n",$outdir;
   printf OUTFL "\t\t\t</hourVMTFraction>\n";
   printf OUTFL "\t\t</parts>\n";
   printf OUTFL "\t</vehicletypevmt>\n";

   printf OUTFL "\t<imcoverage>\n";
   printf OUTFL "\t\t<description><![CDATA[]]></description>\n";
   printf OUTFL "\t\t<parts>\n";
   printf OUTFL "\t\t\t<IMCoverage>\n";
   printf OUTFL "\t\t\t<filename>%s</filename>\n", $repIM[$cntyidx];
   printf OUTFL "\t\t\t</IMCoverage>\n";
   printf OUTFL "\t\t</parts>\n";
   printf OUTFL "\t</imcoverage>\n";

   printf OUTFL "\t</importer>\n";
   printf OUTFL "</moves>\n";

}  # end of dataImporter subroutine


# --- finish up the runspec  file - ====================================================================
sub rspend
{
   printf OUTFL "\t<databaseselections>\n";
   printf OUTFL "\t</databaseselections>\n";
   printf OUTFL "\t<internalcontrolstrategies>\n";

   printf OUTFL "\t<internalcontrolstrategy classname=\"gov.epa.otaq.moves.master.implementation.ghg.internalcontrolstrategies.rateofprogress.RateOfProgressStrategy\"><![CDATA[ useParameters	No ]]></internalcontrolstrategy>\n";
   printf OUTFL "\t</internalcontrolstrategies>\n";
   printf OUTFL "\t<inputdatabase servername=\"\" databasename=\"\" description=\"\"/>\n";
   printf OUTFL "\t<uncertaintyparameters uncertaintymodeenabled=\"false\" numberofrunspersimulation=\"0\" numberofsimulations=\"0\"/>\n";
   printf OUTFL "\t<geographicoutputdetail description=\"LINK\"/>\n";

   printf OUTFL "\t<outputemissionsbreakdownselection>\n";
   printf OUTFL "\t\t<modelyear selected=\"true\"/>\n";
   printf OUTFL "\t\t<fueltype selected=\"true\"/>\n";
   printf OUTFL "\t\t<emissionprocess selected=\"true\"/>\n";
   printf OUTFL "\t\t<onroadoffroad selected=\"true\"/>\n";
   printf OUTFL "\t\t<roadtype selected=\"true\"/>\n";
   printf OUTFL "\t\t<sourceusetype selected=\"true\"/>\n";
   printf OUTFL "\t\t<movesvehicletype selected=\"false\"/>\n";
   printf OUTFL "\t\t<onroadscc selected=\"false\"/>\n";
   printf OUTFL "\t\t<offroadscc selected=\"false\"/>\n";
   printf OUTFL "\t\t<estimateuncertainty selected=\"false\" numberOfIterations=\"2\" keepSampledData=\"false\" keepIterations=\"false\"/>\n";
   printf OUTFL "\t\t<segment selected=\"false\"/>\n";
   printf OUTFL "\t\t<hpclass selected=\"false\"/>\n";
   printf OUTFL "\t</outputemissionsbreakdownselection>\n";

   printf OUTFL "\t<outputdatabase servername=\"%s\" databasename=\"%s\" description=\"\"/>\n",$dbhost,$outputDB;

   printf OUTFL "\t<outputtimestep value=\"Hour\"/>\n";
   printf OUTFL "\t<outputvmtdata value=\"true\"/>\n";
   printf OUTFL "\t<outputsho value=\"true\"/>\n";
   printf OUTFL "\t<outputsh value=\"true\"/>\n";
   printf OUTFL "\t<outputshp value=\"true\"/>\n";
   printf OUTFL "\t<outputshidling value=\"true\"/>\n";
   printf OUTFL "\t<outputstarts value=\"true\"/>\n";
   printf OUTFL "\t<outputpopulation value=\"true\"/>\n";
   printf OUTFL "\t<scaleinputdatabase servername=\"%s\" databasename=\"%s\" description=\"\"/>\n",$dbhost,$scenarioID."_in";
   printf OUTFL "\t<pmsize value=\"0\"/>\n";
   printf OUTFL "\t<outputfactors>\n";
   printf OUTFL "\t\t<timefactors selected=\"true\" units=\"Hours\"/>\n";
   printf OUTFL "\t\t<distancefactors selected=\"true\" units=\"Miles\"/>\n";
   printf OUTFL "\t\t<massfactors selected=\"true\" units=\"Grams\" energyunits=\"Joules\"/>\n";
   printf OUTFL "\t</outputfactors>\n";
   printf OUTFL"\t<savedata>\n";
   printf OUTFL"\t</savedata>\n";
   printf OUTFL "\t<donotexecute>\n";
   printf OUTFL "\t</donotexecute>\n";

   printf OUTFL "\t<generatordatabase shouldsave=\"false\" servername=\"\" databasename=\"\" description=\"\"/>\n";
   printf OUTFL	"\t\t<donotperformfinalaggregation selected=\"false\"/>\n";
   printf OUTFL "\t<lookuptableflags scenarioid=\"%s\" truncateoutput=\"true\" truncateactivity=\"true\"/>\n",$scenarioID;

   printf OUTFL "</runspec>\n";

}  # end of rspend subroutine

sub trim
{
    my $s = shift;
    # remove leading spaces
    $s =~ s/^\s+//;
    # remove trailing spaces
    $s =~ s/\s+$//;
    # remove leading tabs
    $s =~ s/^\t+//;
    # remove trailing tabs
    $s =~ s/\t+$//;
    return $s; 
}

