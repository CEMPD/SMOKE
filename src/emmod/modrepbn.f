
        MODULE MODREPBN

!***********************************************************************
!  Module body starts at line
!
!  DESCRIPTION:
!     This module contains the public data that are used for reporting. The
!     main purpose of this module is to contain the data structures needed
!     for grouping records into bins for QA and reporting.
!
!  PRECONDITIONS REQUIRED:
!
!  SUBROUTINES AND FUNCTIONS CALLED:
!
!  REVISION HISTORY:
!     Created 7/2000 by M. Houyoux
!
!***************************************************************************
!
! Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
!                System
! File: @(#)$Id$
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
! Pathname: $Source$
! Last updated: $Date$ 
!
!****************************************************************************

        INCLUDE 'EMPRVT3.EXT'

!.........  Module-specific parameters
        INTEGER, PARAMETER, PUBLIC :: LV1 = IOVLEN3 + 2     ! "S:" len = 2
        INTEGER, PARAMETER, PUBLIC :: LV2 = IOVLEN3 * 2 + 4 ! ETJOIN len = 2
        INTEGER, PARAMETER, PUBLIC :: LV3 = IOVLEN3 * 3 + 6 

!.........  Define types needed for module
        TYPE :: BYETYPE

            SEQUENCE

            INTEGER            :: ETP        ! output index for emission types
            INTEGER            :: DAT        ! output index for pol/act
            INTEGER            :: AGG        ! output for either

            LOGICAL            :: SPCYN      ! true: speciation applies

        END TYPE

        TYPE :: BYSPVAR

            SEQUENCE

            INTEGER            :: SPC        ! output index for species
            INTEGER            :: ETPSPC     ! output index for emis type || spc
            INTEGER            :: PRCSPC     ! output index for process || spc
            INTEGER            :: SUMETP     ! output index for summed to etype
            INTEGER            :: SUMPOL     ! output index for summed to pol
            INTEGER            :: AGG        ! output for any

        END TYPE

!.........  Output records to be summed into bins
        INTEGER, PUBLIC :: NOUTREC = 0   ! no. of output records to be summed

        INTEGER, ALLOCATABLE, PUBLIC :: OUTSRC ( : )  ! smoke src ID
        INTEGER, ALLOCATABLE, PUBLIC :: OUTBIN ( : )  ! bin number
        INTEGER, ALLOCATABLE, PUBLIC :: OUTCELL( : )  ! cell number or zero

        REAL   , ALLOCATABLE, PUBLIC :: OUTGFAC( : )  ! gridding factor or 1.

!.........  Scalar values for data bins
        INTEGER, PUBLIC :: NOUTBINS = 0  ! no. of output bins

!.........  Grouped output information and data values
        INTEGER, ALLOCATABLE, PUBLIC :: BINCOIDX ( : )   ! index to country name
        INTEGER, ALLOCATABLE, PUBLIC :: BINCYIDX ( : )   ! index to county name
        INTEGER, ALLOCATABLE, PUBLIC :: BINREGN  ( : )   ! region code
        INTEGER, ALLOCATABLE, PUBLIC :: BINRCL   ( : )   ! roadclass code
        INTEGER, ALLOCATABLE, PUBLIC :: BINSMKID ( : )   ! SMOKE source ID
        INTEGER, ALLOCATABLE, PUBLIC :: BINSNMIDX( : )   ! SCC name index
        INTEGER, ALLOCATABLE, PUBLIC :: BINSTIDX ( : )   ! index to state name
        INTEGER, ALLOCATABLE, PUBLIC :: BINX     ( : )   ! x cell
        INTEGER, ALLOCATABLE, PUBLIC :: BINY     ( : )   ! y cell

        REAL   , ALLOCATABLE, PUBLIC :: BINDATA ( :,: ) ! output data values

        CHARACTER*1, ALLOCATABLE, PUBLIC :: BINELEV( : )! elevated flag
        CHARACTER(LEN=SCCLEN3), ALLOCATABLE, PUBLIC :: BINSCC( : )   ! SCC

!.........  Arrays for determining output from emission types to report columns
!.........  Dimensioned ( NIPPA, NREPORT )
        TYPE( BYETYPE ), ALLOCATABLE, PUBLIC :: TODOUT( :,: )

        CHARACTER(LEN=LV2), ALLOCATABLE, PUBLIC :: ETPNAM( : ) ! emis type names
        CHARACTER(LEN=LV1), ALLOCATABLE, PUBLIC :: DATNAM( : ) ! pol/act names

!.........  Arrays for determining output from species variables to report cols
!.........  Dimensioned ( NSVARS, NREPORT )
        INTEGER, PUBLIC :: NSVARS = 0       ! no. speciation variables

        TYPE( BYSPVAR ), ALLOCATABLE, PUBLIC :: TOSOUT( :,: )          

        CHARACTER(LEN=LV1), ALLOCATABLE, PUBLIC :: SPCNAM   ( : )  ! valid species names
        CHARACTER(LEN=LV3), ALLOCATABLE, PUBLIC :: ETPSPCNAM( : )  ! valid emis type || species
        CHARACTER(LEN=LV2), ALLOCATABLE, PUBLIC :: PRCSPCNAM( : )  ! valid process || species
        CHARACTER(LEN=LV2), ALLOCATABLE, PUBLIC :: SUMETPNAM( : )  ! S: || emis type names
        CHARACTER(LEN=LV1), ALLOCATABLE, PUBLIC :: SUMPOLNAM( : )  ! S: || pollutant names

        CHARACTER(LEN=IOULEN3), ALLOCATABLE, PUBLIC :: SLUNIT( : ) ! spc var units
        CHARACTER(LEN=IOULEN3), ALLOCATABLE, PUBLIC :: SSUNIT( : ) ! spc var units

!.........  Arrays for referencing input data needed across whole program run
        INTEGER, PUBLIC :: NDATIN = 0 ! Actual number of data vars input
        INTEGER, PUBLIC :: NSPCIN = 0 ! Actual number of speciation vars input

        INTEGER, ALLOCATABLE, PUBLIC :: DATIDX( : )  ! index global-to-input
        INTEGER, ALLOCATABLE, PUBLIC :: SPCIDX( : )  ! index global-to-input
        INTEGER, ALLOCATABLE, PUBLIC :: SPCTODAT( : )! index spc-to-data

        LOGICAL, ALLOCATABLE, PUBLIC :: DATOUT( : )  ! true: inven data needed
        LOGICAL, ALLOCATABLE, PUBLIC :: SPCOUT( : )  ! true: spc factors needed

        CHARACTER(LEN=IOVLEN3), ALLOCATABLE, PUBLIC :: INNAMES( : ) ! emis names

!.........  Unique species list and count
        INTEGER, PUBLIC :: NMSPC = 0  ! no. species in whole run
        CHARACTER(LEN=LV1), ALLOCATABLE, PUBLIC :: EMNAM( : ) ! species names

        END MODULE MODREPBN
