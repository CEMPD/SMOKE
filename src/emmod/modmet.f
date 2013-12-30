        MODULE MODMET

!***********************************************************************
!  Module body starts at line
!
!  DESCRIPTION:
!     This module contains the derived meteorology data for applying emission
!     factors to activity data.
!
!  PRECONDITIONS REQUIRED:
!
!  SUBROUTINES AND FUNCTIONS CALLED:
!
!  REVISION HISTORY:
!     Created 6/99 by M. Houyoux
!
!***************************************************************************
!
! Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
!                System
! File: @(#)$Id$
!
! COPYRIGHT (C) 2004, Environmental Modeling for Policy Development
! All Rights Reserved
! 
! Carolina Environmental Program
! University of North Carolina at Chapel Hill
! 137 E. Franklin St., CB# 6116
! Chapel Hill, NC 27599-6116
! 
! smoke@unc.edu
!
! Pathname: $Source$
! Last updated: $Date$ 
!
!****************************************************************************

        IMPLICIT NONE
        
        INCLUDE 'EMPRVT3.EXT'

!...........   Setting for range of valid min/max temperatures
        INTEGER, PUBLIC :: NFUEL = 0      ! no of fuelmonth
        REAL,    PUBLIC :: MINTEMP = 0.   ! minimum temperature
        REAL,    PUBLIC :: MAXTEMP = 0.   ! maximum temperature

!...........   Source-based meteorology data (dim: NSRC)
        REAL,    ALLOCATABLE, PUBLIC :: TASRC   ( : )   ! temperature in Kelvin
        REAL,    ALLOCATABLE, PUBLIC :: QVSRC   ( : )   ! water vapor mixing ratio
        REAL,    ALLOCATABLE, PUBLIC :: PRESSRC ( : )   ! pressure in pascals
        REAL,    ALLOCATABLE, PUBLIC :: MAXTSRC ( : )   ! max temp per county
        REAL,    ALLOCATABLE, PUBLIC :: MINTSRC ( : )   ! min temp per county

!...........   Hourly meteorology data
        REAL,    ALLOCATABLE, PUBLIC :: TKHOUR  ( :,: ) ! temps by source per hour
        INTEGER, ALLOCATABLE, PUBLIC :: NTKHOUR ( :,: ) ! no of hours for monthly average for each source 
        REAL,    ALLOCATABLE, PUBLIC :: MAXTDAY ( : )   ! max temp per county per day
        REAL,    ALLOCATABLE, PUBLIC :: MINTDAY ( : )   ! min temp per county per day

        REAL,    ALLOCATABLE, PUBLIC :: RHTBIN  ( :,:,: ) ! relative humidity by temp bins 
        INTEGER, ALLOCATABLE, PUBLIC :: NRHTBIN ( :,:,: ) ! no of hours by tempe bins 

        INTEGER, ALLOCATABLE, PUBLIC :: FUELIDX ( :,: ) ! monthID for fuelmonth ref county
        INTEGER, ALLOCATABLE, PUBLIC :: FUELCNTY( :,: ) ! no of county per ref county per fulemonth
        REAL,    ALLOCATABLE, PUBLIC :: RHFUEL  ( :,: ) ! relative humidity by county per fuelmonth
        REAL,    ALLOCATABLE, PUBLIC :: MAXTFUEL( :,: ) ! max temperature by county per fuelmonth
        REAL,    ALLOCATABLE, PUBLIC :: MINTFUEL( :,: ) ! min temperature by county per fuelmonth
        REAL,    ALLOCATABLE, PUBLIC :: TKFUEL  ( :,:,: ) ! temp profile by county per fuelmonth

        REAL,    ALLOCATABLE, PUBLIC :: TPCNTY( : )     ! daily temps by county
        REAL,    ALLOCATABLE, PUBLIC :: RHCNTY( : )     ! daily RH by county

        REAL,    ALLOCATABLE, PUBLIC :: TDYCNTY ( : )   ! daily temps by county
        REAL,    ALLOCATABLE, PUBLIC :: QVDYCNTY( : )   ! daily mixing ratios by county
        REAL,    ALLOCATABLE, PUBLIC :: BPDYCNTY( : )   ! daily barometric pressure by county
        CHARACTER(FIPLEN3), ALLOCATABLE, PUBLIC :: DYCODES ( : ) ! FIPS codes for daily counties

        REAL,    ALLOCATABLE, PUBLIC :: TWKCNTY ( : )   ! weekly temps by county
        REAL,    ALLOCATABLE, PUBLIC :: QVWKCNTY( : )   ! weekly mixing ratios by county
        REAL,    ALLOCATABLE, PUBLIC :: BPWKCNTY( : )   ! weekly barometric pressure by county
        CHARACTER(FIPLEN3), ALLOCATABLE, PUBLIC :: WKCODES ( : ) ! FIPS codes for weekly counties

        REAL,    ALLOCATABLE, PUBLIC :: TMNCNTY ( : )   ! monthly temps by county
        REAL,    ALLOCATABLE, PUBLIC :: QVMNCNTY( : )   ! monthly mixing ratios by county
        REAL,    ALLOCATABLE, PUBLIC :: BPMNCNTY( : )   ! monthly barometric pressure by county
        CHARACTER(FIPLEN3), ALLOCATABLE, PUBLIC :: MNCODES ( : ) ! FIPS codes for monthly counties

        REAL,    ALLOCATABLE, PUBLIC :: TEPCNTY ( : )   ! episode temps by county
        REAL,    ALLOCATABLE, PUBLIC :: QVEPCNTY( : )   ! episode mixing ratios by county
        REAL,    ALLOCATABLE, PUBLIC :: BPEPCNTY( : )   ! episode barometric pressure by county
        CHARACTER(FIPLEN3), ALLOCATABLE, PUBLIC :: EPCODES ( : ) ! FIPS codes for episode counties

!...........   Daily meteorology data
        REAL,    ALLOCATABLE, PUBLIC :: BPDAY( : )      ! average daily barometric pressure by county

        END MODULE MODMET
