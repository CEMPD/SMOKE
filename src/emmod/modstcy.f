
        MODULE MODSTCY

!***********************************************************************
!
!  DESCRIPTION:
!     This module contains the public allocatable arrays for doing country, 
!     state and county summaries
!
!  PRECONDITIONS REQUIRED:
!
!  SUBROUTINES AND FUNCTIONS CALLED:
!
!  REVISION HISTORY:
!     
!
!***********************************************************************
!
! Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
!                System
! File: @(#)$Id$
!
! COPYRIGHT (C) 1999, MCNC--North Carolina Supercomputing Center
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
! Pathname: $Source$
! Last updated: $Date$ 
!
!****************************************************************************

        INCLUDE 'EMPRVT3.EXT'   !  emissions private parameters

!.........  Indices from per-source inventory arrays to county array
        INTEGER, ALLOCATABLE, PUBLIC :: AICNY( : )  ! dim NASRC
        INTEGER, ALLOCATABLE, PUBLIC :: MICNY( : )  ! dim NMSRC
        INTEGER, ALLOCATABLE, PUBLIC :: PICNY( : )  ! dim NPSRC

!.........  Codes, names, and their dimensions
        INTEGER, PUBLIC :: NCOUNTRY    ! Number of counties
        INTEGER, PUBLIC :: NSTATE      ! Number of countries/states
        INTEGER, PUBLIC :: NCOUNTY     ! Number of countries/states/counties

        INTEGER, ALLOCATABLE, PUBLIC :: CTRYCOD( : ) ! country codes
        INTEGER, ALLOCATABLE, PUBLIC :: STATCOD( : ) ! country/state codes
        INTEGER, ALLOCATABLE, PUBLIC :: CNTYCOD( : ) ! country/state/county code

        CHARACTER*20, ALLOCATABLE, PUBLIC :: CTRYNAM( : ) ! Country names
        CHARACTER*20, ALLOCATABLE, PUBLIC :: STATNAM( : ) ! State names
        CHARACTER*20, ALLOCATABLE, PUBLIC :: CNTYNAM( : ) ! County names

!.........  Time zone by county array
        INTEGER, ALLOCATABLE, PUBLIC :: CNTYTZON( : )  ! time zone by county

!.........  Daylight time exemptions (false is exempt)
        LOGICAL, ALLOCATABLE, PUBLIC :: USEDAYLT( : )  ! true: use daylight time

!.........  ORIS list with state and county codes
        INTEGER, PUBLIC :: NORIS     ! Number of all ORIS codes

        INTEGER, ALLOCATABLE, PUBLIC :: ORISFIP( : ) ! country/state/county code
        CHARACTER(LEN=ORSLEN3), ALLOCATABLE, PUBLIC :: ORISLST( : ) ! ORIS codes
        CHARACTER(LEN=DSCLEN3), ALLOCATABLE, PUBLIC :: ORISDSC( : ) ! ORIS desc
        

        END MODULE MODSTCY
