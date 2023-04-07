
        MODULE MODMVSMRG

!***********************************************************************
!  Module body starts at line
!
!  DESCRIPTION:
!     This module contains the public variables and allocatable arrays 
!     used by Movesmrg.
!
!  PRECONDITIONS REQUIRED:
!
!  SUBROUTINES AND FUNCTIONS CALLED:
!
!  REVISION HISTORY:
!     Created 4/10 by C. Seppanen
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

C.........  Program settings
        LOGICAL, PUBLIC :: RPDFLAG  ! mode is rate-per-distance
        LOGICAL, PUBLIC :: RPHFLAG  ! mode is rate-per-hour for extended idle
        LOGICAL, PUBLIC :: RPVFLAG  ! mode is rate-per-vehicle
        LOGICAL, PUBLIC :: RPPFLAG  ! mode is rate-per-profile
        LOGICAL, PUBLIC :: RPSFLAG  ! mode is rate-per-start
        LOGICAL, PUBLIC :: ONIFLAG  ! mode is off-network idling
        LOGICAL, PUBLIC :: MTMP_OUT ! output temporal intermediate MTMP files
        LOGICAL, PUBLIC :: MOPTIMIZE ! memory optimize option
        
        CHARACTER(300), PUBLIC :: MVFILDIR  ! directory for MOVES output files

C.........  Meteorology information
        CHARACTER(IOVLEN3), PUBLIC :: TVARNAME  ! name of temperature variable to read
        CHARACTER(IOVLEN3), PUBLIC :: DATANAME  ! name of activity data name from Temporal
        CHARACTER(16), PUBLIC :: METNAME        ! logical name for meteorology file

C.........  Average min and max temperatures
        REAL, ALLOCATABLE, PUBLIC :: AVGMIN( :,:,: )  ! minimum monthly temperature for each county
        REAL, ALLOCATABLE, PUBLIC :: AVGMAX( :,:,: )  ! maximum monthly temperature for each county

C.........  Text file unit numbers        
        INTEGER, PUBLIC :: XDEV         ! unit no. for county cross-reference file
        INTEGER, PUBLIC :: MDEV         ! unit no. for county fuel month file
        INTEGER, PUBLIC :: FDEV         ! unit no. for reference county file listing
        INTEGER, PUBLIC :: MMDEV        ! unit no. for Met4moves output file

C.........  Source to grid cell lookups
        INTEGER, ALLOCATABLE, PUBLIC :: NSRCCELLS( : )      ! no. cells for each source
        INTEGER, ALLOCATABLE, PUBLIC :: SRCCELLS( :,: )     ! list of cells for each source
        REAL,    ALLOCATABLE, PUBLIC :: SRCCELLFRACS( :,: ) ! grid cell fractions for each source

C.........  Reference county information
        INTEGER, ALLOCATABLE, PUBLIC :: NREFSRCS( : )       ! no. sources for each ref county
        INTEGER, ALLOCATABLE, PUBLIC :: REFSRCS( :,: )      ! list of srcs for each ref county
        CHARACTER(100), ALLOCATABLE, PUBLIC :: MRCLIST( :,: ) ! emfac file for each ref county and month

C.........  NONHAPTOG calculation information
        INTEGER, PUBLIC :: NHAP                                ! number of HAPs to subtract
        CHARACTER(IOVLEN3), ALLOCATABLE, PUBLIC :: HAPNAM( : ) ! list of HAPs

C.........  Emission factors data
        REAL,                 PUBLIC :: TEMPBIN           ! temperature buffer for max/min profiles
        REAL,                 PUBLIC :: MNTEMP            ! temperature for lowest bin for EMIS_TABLE output file
        REAL,                 PUBLIC :: MXTEMP            ! temperature for highest bin for EMIS_TABLE output file
        REAL,                 PUBLIC :: TMPINC            ! temperature increment for EMIS_TABLE output file
        INTEGER,              PUBLIC :: NTBINS = 1        ! no of temperature bins for EMIS TABLE output file (defualt=1)
        LOGICAL, PUBLIC              :: ETABLEFLAG        ! ture: output precomputed gridded hourly emissions by temp bin to EMIS_TABLE
        LOGICAL, PUBLIC              :: NOXADJFLAG        ! true: apply humidity adjustment for NOx emissions
        LOGICAL, PUBLIC              :: NOXADJEQS         ! true: apply older humidity adjustment for NOx emissions
        CHARACTER(FLTLEN3), PUBLIC   :: GASFLTYP          ! gasoline fuel type code in NEI MOVES SCC '01'
        CHARACTER(FLTLEN3), PUBLIC   :: DISFLTYP          ! diesel fuel type code in NEI MOVES SCC '02'
        CHARACTER(FLTLEN3), PUBLIC   :: CNGFLTYP          ! CNG fuel type code in NEI MOVES SCC '03'
        CHARACTER(FLTLEN3), PUBLIC   :: LPGFLTYP          ! LPG fuel type code in NEI MOVES SCC '04'
        CHARACTER(FLTLEN3), PUBLIC   :: ETHFLTYP          ! ethanol fuel type code in NEI MOVES SCC '05'
        LOGICAL, ALLOCATABLE, PUBLIC :: GASFL( : )        ! index for gasoline fuel type SCC
        LOGICAL, ALLOCATABLE, PUBLIC :: DISFL( : )        ! index for diesel fuel type SCC
        LOGICAL, ALLOCATABLE, PUBLIC :: CNGFL( : )        ! index for CNG fuel type SCC
        LOGICAL, ALLOCATABLE, PUBLIC :: LPGFL( : )        ! index for LPG fuel type SCC
        LOGICAL, ALLOCATABLE, PUBLIC :: ETHFL( : )        ! index for ethanol fuel type SCC

        INTEGER, ALLOCATABLE, PUBLIC :: EMPOLIDX( : )     ! index of emission pollutant name
        INTEGER,                        PUBLIC :: NMVSPOLS         ! number of MOVES pollutants/species
        CHARACTER(IOVLEN3),ALLOCATABLE, PUBLIC :: MVSPOLNAMS( : )  ! arry for poll/spc names


        INTEGER, PUBLIC :: NEMTEMPS                       ! no. temperatures for current emision factors
        REAL, ALLOCATABLE, PUBLIC :: EMTEMPS( : )         ! list of temps for emission factors
        REAL, ALLOCATABLE, PUBLIC :: EMXTEMPS( : )        ! list of max. temps in profiles
        INTEGER, ALLOCATABLE, PUBLIC :: EMTEMPIDX( : )    ! index to sorted temperature profiles

        REAL, ALLOCATABLE, PUBLIC :: RPDEMFACS( :,:,:,: )  ! rate-per-distance emission factors
                                                           ! SCC, speed bin, temp, process, pollutant

        REAL, ALLOCATABLE, PUBLIC :: RPHEMFACS( :,:,: )    ! rate-per-hour emission factors
                                                           ! SCC, temp, process, pollutant

        REAL, ALLOCATABLE, PUBLIC :: RPVEMFACS( :,:,:,:,: )  ! rate-per-vehicle emission factors
                                                             ! day, SCC, hour, temp, process, pollutant

        REAL, ALLOCATABLE, PUBLIC :: RPPEMFACS( :,:,:,:,: )  ! rate-per-profile emission factors
                                                             ! day, SCC, hour, temp profile, process, pollutant

C.........  Hourly speed data and control factor data
        LOGICAL, PUBLIC :: SPDPROFLAG                     ! use hourly speed data
        LOGICAL, PUBLIC :: SPDISTFLAG                     ! use hourly speed data
        REAL, ALLOCATABLE, PUBLIC :: SPDPRO( :,:,:,: )    ! indexes: FIP, SCC, weekend/weekday, local hour
        REAL, ALLOCATABLE, PUBLIC :: SPDIST( :,:,:,:,: )  ! indexes: FIP, SCC, weekend/weekday, local hour, spd bins
        LOGICAL, PUBLIC :: CFFLAG                      ! use control factor data
        LOGICAL, PUBLIC :: EXPCFFLAG                   ! use explicit poll/species specific control factor data
        LOGICAL, PUBLIC :: REFCFFLAG                   ! use reference county-specific control factor data
        REAL, ALLOCATABLE, PUBLIC :: CFPRO( :,:,:,: )  ! factor indexes: FIP, SCC, pollutant,month
        REAL, ALLOCATABLE, PUBLIC :: CFITC( :,:,:,: )  ! intercept indexes: FIP, SCC, pollutant,month

C.........  Index from per-source inventory array to INVSCC array (based on MICNY in MODSTCY)
        INTEGER, ALLOCATABLE, PUBLIC :: MISCC( : )     ! dim NMSRC

C.........  Speciation matrixes (mole- and mass-based)
        CHARACTER(16), PUBLIC :: MSNAME_L
        CHARACTER(16), PUBLIC :: PSNAME_L
        CHARACTER(16), PUBLIC :: MSNAME_S
        CHARACTER(16), PUBLIC :: PSNAME_S
        CHARACTER(16), PUBLIC :: GRDENV
        CHARACTER(16), PUBLIC :: TOTENV
        INTEGER, PUBLIC :: MNSMATV_L = 0  ! number of pol-to-species combinations
        INTEGER, PUBLIC :: PNSMATV_L = 0  ! number of pol-to-species combinations
        INTEGER, PUBLIC :: MNSMATV_S = 0
        INTEGER, PUBLIC :: PNSMATV_S = 0
        CHARACTER(PLSLEN3), ALLOCATABLE, PUBLIC :: MSVDESC_L( : )  ! pollutant-to-species names
        CHARACTER(PLSLEN3), ALLOCATABLE, PUBLIC :: PSVDESC_L( : )  ! pollutant-to-species names
        CHARACTER(PLSLEN3), ALLOCATABLE, PUBLIC :: MSVDESC_S( : )
        CHARACTER(PLSLEN3), ALLOCATABLE, PUBLIC :: PSVDESC_S( : )
        CHARACTER(PLSLEN3), ALLOCATABLE, PUBLIC :: MSVUNIT_L( : )  ! pollutant-to-species units
        CHARACTER(PLSLEN3), ALLOCATABLE, PUBLIC :: PSVUNIT_L( : )  ! pollutant-to-species units
        CHARACTER(PLSLEN3), ALLOCATABLE, PUBLIC :: MSVUNIT_S( : )
        CHARACTER(PLSLEN3), ALLOCATABLE, PUBLIC :: PSVUNIT_S( : )
        REAL, ALLOCATABLE, PUBLIC :: MSMATX_L( :,: )  ! speciation matrix, dim nmsrc
        REAL, ALLOCATABLE, PUBLIC :: PSMATX_L( :,: )  ! speciation matrix, dim nmsrc
        REAL, ALLOCATABLE, PUBLIC :: MSMATX_S( :,: )
        REAL, ALLOCATABLE, PUBLIC :: PSMATX_S( :,: )
        CHARACTER(IOULEN3), ALLOCATABLE, PUBLIC :: SPCUNIT_L( : ) ! speciation units
        CHARACTER(IOULEN3), ALLOCATABLE, PUBLIC :: SPCUNIT_S( : )

C.........  Non-speciated emissions reporting
        LOGICAL, ALLOCATABLE, PUBLIC :: EANAMREP( : )  ! indicates when non-speciation emissions should be saved

        END MODULE MODMVSMRG
