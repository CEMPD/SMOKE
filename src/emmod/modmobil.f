
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
! COPYRIGHT (C) 2000, MCNC--North Carolina Supercomputing Center
! All Rights Reserved
!
! See file COPYRIGHT for conditions of use.
!
! Environmental Programs Group
! MCNC--North Carolina Supercomputing Center
! P.O. Box 12889
! Research Triangle Park, NC  27709-2889
!
! env_progs@mcnc.org
!
! Pathname: %P%
! Last updated: %G% %U%
!
!****************************************************************************

        INCLUDE 'EMPRVT3.EXT'

!.........  Mobile source characteristics tables
        INTEGER                                     :: NVTYPE ! no. veh types
        INTEGER               , ALLOCATABLE, PUBLIC :: IVTIDLST( : ) ! IDs
        CHARACTER(LEN=VTPLEN3), ALLOCATABLE, PUBLIC :: CVTYPLST( : ) ! names

        INTEGER                                     :: NRCLAS ! no. road classes
        INTEGER               , ALLOCATABLE, PUBLIC :: AMSRDCLS( : ) ! AIRS rclas
        INTEGER               , ALLOCATABLE, PUBLIC :: RDWAYTYP( : ) ! roadway type

!.........  Unsorted VMT Mix table
        REAL                  , ALLOCATABLE, PUBLIC :: VMTMIXA ( :,: )

!.........  Unsorted speeds table
        REAL                  , ALLOCATABLE, PUBLIC :: SPDTBLA ( : )

!.........  Sorted SCC table for mobile sources
        INTEGER, PUBLIC :: NSCCTBL
        CHARACTER(LEN=VTPLEN3+RWTLEN3),ALLOCATABLE, PUBLIC:: SCCRVC( : )
        CHARACTER(LEN=SCCLEN3), ALLOCATABLE, PUBLIC :: SCCTBL ( : )

        END MODULE MODMOBIL
