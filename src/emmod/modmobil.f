
        MODULE MODMOBIL

!***********************************************************************
!  Module body starts at line
!
!  DESCRIPTION:
!     This module contains the public data that is used for mobile sources
!     only.
!
!  PRECONDITIONS REQUIRED:
!
!  SUBROUTINES AND FUNCTIONS CALLED:
!
!  REVISION HISTORY:
!     Created 2/2000 by M. Houyoux
!
!***************************************************************************
!
! Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
!                System
! File: %W%
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
! Pathname: %P%
! Last updated: %G% %U%
!
!****************************************************************************

        IMPLICIT NONE

        INCLUDE 'EMPRVT3.EXT'

!.........  Mobile source characteristics tables
        INTEGER                                 :: NVTYPE ! no. veh types
        INTEGER,            ALLOCATABLE, PUBLIC :: IVTIDLST( : ) ! IDs
        CHARACTER(VTPLEN3), ALLOCATABLE, PUBLIC :: CVTYPLST( : ) ! names

        INTEGER                                 :: NRCLAS ! no. road classes
        INTEGER,            ALLOCATABLE, PUBLIC :: AMSRDCLS( : ) ! AIRS rclas
        INTEGER,            ALLOCATABLE, PUBLIC :: RDWAYTYP( : ) ! roadway type

!.........  Unsorted VMT Mix table
        REAL,               ALLOCATABLE, PUBLIC :: VMTMIXA ( :,: )

!.........  Unsorted speeds table
        REAL,               ALLOCATABLE, PUBLIC :: SPDTBLA ( : )

!.........  Sorted SCC table for mobile sources
        INTEGER, PUBLIC :: NSCCTBL
        CHARACTER(VTPLEN3+RWTLEN3),ALLOCATABLE, PUBLIC:: SCCRVC( : )
        CHARACTER(SCCLEN3), ALLOCATABLE, PUBLIC :: SCCTBL ( : )

!.........  Sorted speed profiles
        INTEGER                                 :: NSPDPROF        ! no. speed profiles
        INTEGER,            ALLOCATABLE, PUBLIC :: SPDNUMS ( : )   ! profile numbers
        REAL,               ALLOCATABLE, PUBLIC :: SPDPROFS( :,: ) ! 24 hour speed profiles

        END MODULE MODMOBIL
