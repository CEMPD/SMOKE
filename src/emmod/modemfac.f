
        MODULE MODEMFAC

!***********************************************************************
!  Module body starts at line
!
!  DESCRIPTION:
!     This module contains the public allocatable arrays for using emission
!     factors with activity data
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

!.........  Sizes of the arrays...
        INTEGER, PUBLIC :: NEFS      ! no. EFs

!.........  Unique lists of source characteristics and associated arrays...
        CHARACTER(LEN=IODLEN3), ALLOCATABLE, PUBLIC :: EFSDSC( : )
        CHARACTER(LEN=IOVLEN3), ALLOCATABLE, PUBLIC :: EFSNAM( : )
        CHARACTER(LEN=IOULEN3), ALLOCATABLE, PUBLIC :: EFSUNT( : )
        
!.........  Variables for EF subtraction
        CHARACTER(LEN=IOULEN3), ALLOCATABLE, PUBLIC :: SUBPOLS( : )   ! final subtraction list
        CHARACTER(LEN=IOULEN3), PUBLIC :: INPUTHC   ! input hydrocarbon name
        CHARACTER(LEN=IOULEN3), PUBLIC :: OUTPUTHC  ! output hydrocarbon name
        INTEGER               , PUBLIC :: NSUBPOL   ! no. pols to subtract

!.........  Emission processes and emission types (proccess//pol) and related
!.........  Second dimension is number of activities in the inventory
        INTEGER               , PUBLIC:: MXETYPE     ! Max NETYPE
        INTEGER               , PUBLIC:: NEPOL       ! no. pols from emis types

        INTEGER, ALLOCATABLE, PUBLIC:: NETYPE( : )   ! no. e-types per actvt
        INTEGER, ALLOCATABLE, PUBLIC:: EMTIDX( :,: ) ! index to EMTPOL

        CHARACTER(LEN=IODLEN3), ALLOCATABLE, PUBLIC:: EMTPOL( : )!pols of EMTNAM
        CHARACTER(LEN=IOVLEN3), ALLOCATABLE, PUBLIC:: EMTNAM( :,: )! em type
        CHARACTER(LEN=IOULEN3), ALLOCATABLE, PUBLIC:: EMTUNT( :,: )! et units
        CHARACTER(LEN=IODLEN3), ALLOCATABLE, PUBLIC:: EMTDSC( :,: )! et descript

!.........  Variables for running MOBILE6
        INTEGER, ALLOCATABLE, PUBLIC :: SCENLIST ( :,: )  ! scenario number and local-as-arterial
                                                          ! flag for each source
                                                          
        CHARACTER*300, ALLOCATABLE :: M6LIST( : )  ! contents of M6LIST file
        
        INTEGER, PUBLIC :: NUMSCEN                  ! total no. of M6 scenarios

!.........  Derived data type for emissions pointer array
        TYPE :: EF_PTR
            REAL, DIMENSION( :,:,:,:,: ), POINTER :: PTR  ! dim: scen, poll, veh, road, hour
        END TYPE
        
        TYPE( EF_PTR ), ALLOCATABLE, PUBLIC :: EMISSIONS( : ) ! master emission factors array
                                                              ! size: no. emission processes

!.........  Variables for dealing with user-defined toxic pollutants
        INTEGER, PUBLIC :: NTOTHAPS                             ! total number of HAPS
        INTEGER,           ALLOCATABLE, PUBLIC :: PMHAPS( : )   ! indicates if HAP is PM based (size: NMAP)
        CHARACTER(LEN=16), ALLOCATABLE, PUBLIC :: HAPNAMES( : ) ! pollutant names (size: NTOTHAPS)
        REAL,              ALLOCATABLE, PUBLIC :: HAPEFS( :,:,:,:,:,: ) ! HAP emission factors
                                                                      ! size (SCEN,IV,IP,EF,IH,IFAC)     

        END MODULE MODEMFAC
