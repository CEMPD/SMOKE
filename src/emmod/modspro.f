        MODULE MODSPRO

!***********************************************************************
!  Module body starts at line
!
!  DESCRIPTION:
!     This module contains the public data for the speciation profiles
!
!  PRECONDITIONS REQUIRED:
!
!  SUBROUTINES AND FUNCTIONS CALLED:
!
!  REVISION HISTORY:
!     Created 3/99 by M. Houyoux
!
!***************************************************************************
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

        INCLUDE 'EMPRVT3.EXT'

!.........  NOTE: The unique speciation profiles and speciation tables are
!.........        used multiple times: once for each inventory pollutant that 
!.........        needs speciation

!.........  Unique speciation profiles list and index to sorted tables
        INTEGER,              PUBLIC :: NSPROF         ! Number in unique list

        ! position in INPRF of start of each profile
        INTEGER, ALLOCATABLE, PUBLIC :: IDXSPRO ( : )

        ! number of species in this profile
        INTEGER, ALLOCATABLE, PUBLIC :: NSPECIES( : )
  
        ! index to all-species list for this profile
        INTEGER, ALLOCATABLE, PUBLIC :: IDXSSPEC( :,: )
  
        ! unique list of each profile for searching
        CHARACTER(LEN=SPNLEN3), ALLOCATABLE, PUBLIC :: SPROFN( : )

!.........  Sorted speciation tables
        INTEGER,              PUBLIC :: MXSPFUL   ! Max no. in unprocessed table
        INTEGER,              PUBLIC :: MXSPEC    ! max no. of species per pol
        INTEGER,              PUBLIC :: NSPFUL    ! Number in unprocessed table

        REAL   , ALLOCATABLE, PUBLIC :: MOLEFACT( : ) ! mole-based spec factors
        REAL   , ALLOCATABLE, PUBLIC :: MASSFACT( : ) ! mass-based spec factors

        ! speciation profile codes
	CHARACTER(LEN=SPNLEN3), ALLOCATABLE, PUBLIC :: INPRF( : )
  
        ! names of species
	CHARACTER(LEN=IOVLEN3), ALLOCATABLE, PUBLIC :: SPECID( : )

!.........  Table of species names per inventory pollutant

	CHARACTER(LEN=IOVLEN3), ALLOCATABLE, PUBLIC :: SPCNAMES( :,: )

!.........  Sorted groups of pollutant to pollutant conversion factors

!.........  Default FIPS code=0, SCC=0 (for all pollutants)
        REAL                  , ALLOCATABLE, PUBLIC :: CNVFC00( : )

!.........  FIPS code=0, SCC=all (for all pollutants)
        INTEGER                            , PUBLIC :: NCNV1
        REAL                  , ALLOCATABLE, PUBLIC :: CNVFC01( :,: ) 
        CHARACTER(LEN=SCCLEN3), ALLOCATABLE, PUBLIC :: CNVRT01( : )

!.........  FIPS code=country/state default, SCC=all (for all pollutants)
        INTEGER                            , PUBLIC :: NCNV2
        REAL                  , ALLOCATABLE, PUBLIC :: CNVFC02( :,: ) 
        CHARACTER(LEN=STSLEN3), ALLOCATABLE, PUBLIC :: CNVRT02( : )

!.........  FIPS code=all, SCC=all (for all pollutants)
        INTEGER                            , PUBLIC :: NCNV3
        REAL                  , ALLOCATABLE, PUBLIC :: CNVFC03( :,: ) 
        CHARACTER(LEN=FPSLEN3), ALLOCATABLE, PUBLIC :: CNVRT03( : )

        END MODULE MODSPRO
