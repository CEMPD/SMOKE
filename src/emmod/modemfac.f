
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

!.........  Sizes of the arrays...
        INTEGER, PUBLIC :: NEFS      ! no. EFs

!.........  Unique lists of source characteristics and associated arrays...
        CHARACTER(IODLEN3), ALLOCATABLE, PUBLIC :: EFSDSC( : )
        CHARACTER(IOVLEN3), ALLOCATABLE, PUBLIC :: EFSNAM( : )
        CHARACTER(IOULEN3), ALLOCATABLE, PUBLIC :: EFSUNT( : )
        
!.........  Variables for EF subtraction
        CHARACTER(IOULEN3), ALLOCATABLE, PUBLIC :: SUBPOLS( : )   ! final subtraction list
        CHARACTER(IOULEN3), PUBLIC :: INPUTHC   ! input hydrocarbon name
        CHARACTER(IOULEN3), PUBLIC :: OUTPUTHC  ! output hydrocarbon name
        INTEGER,            PUBLIC :: NSUBPOL   ! no. pols to subtract

!.........  Emission processes and emission types (proccess//pol) and related
!.........  Second dimension is number of activities in the inventory
        INTEGER,            PUBLIC:: MXETYPE     ! Max NETYPE
        INTEGER,            PUBLIC:: NEPOL       ! no. pols from emis types

        INTEGER, ALLOCATABLE, PUBLIC:: NETYPE( : )   ! no. e-types per actvt
        INTEGER, ALLOCATABLE, PUBLIC:: EMTIDX( :,: ) ! index to EMTPOL

        CHARACTER(IODLEN3), ALLOCATABLE, PUBLIC:: EMTPOL( : )!pols of EMTNAM
        CHARACTER(IOVLEN3), ALLOCATABLE, PUBLIC:: EMTNAM( :,: )! em type
        CHARACTER(IOULEN3), ALLOCATABLE, PUBLIC:: EMTUNT( :,: )! et units
        CHARACTER(IODLEN3), ALLOCATABLE, PUBLIC:: EMTDSC( :,: )! et descript

!.........  Variables for running MOBILE6
        INTEGER, ALLOCATABLE, PUBLIC :: SCENLIST ( :,: )  ! scenario number and local-as-arterial
                                                          ! flag for each source
                                                          
        CHARACTER(300), ALLOCATABLE :: M6LIST( : )  ! contents of M6LIST file
        
        INTEGER, PUBLIC :: NUMSCEN                  ! total no. of M6 scenarios

!.........  Derived data type for emissions pointer array
        TYPE :: EF_PTR
            REAL, DIMENSION( :,:,:,:,: ), POINTER :: PTR  ! dim: scen, poll, veh, road, hour
        END TYPE
        
        TYPE( EF_PTR ), ALLOCATABLE, PUBLIC :: EMISSIONS( : ) ! master emission factors array
                                                              ! size: no. emission processes

!.........  Variables for dealing with user-defined toxic pollutants
        INTEGER, PUBLIC :: NTOTHAPS                             ! total number of HAPS
        INTEGER,       ALLOCATABLE, PUBLIC :: PMHAPS( : )   ! indicates if HAP is PM based (size: NMAP)
        CHARACTER(16), ALLOCATABLE, PUBLIC :: HAPNAMES( : ) ! pollutant names (size: NTOTHAPS)
        REAL,          ALLOCATABLE, PUBLIC :: HAPEFS( :,:,:,:,:,: ) ! HAP emission factors
                                                                    ! size (SCEN,IVEH,IP,EF,IH,IFAC)     
!.........  Emission factor arrays        
        CHARACTER(256), ALLOCATABLE, PUBLIC :: EFLIST( : )    ! listing of emission factor file names
        CHARACTER(16) , ALLOCATABLE, PUBLIC :: EFLOGS( : )    ! listing of ef logical file names
        INTEGER       , ALLOCATABLE, PUBLIC :: EFDAYS( :,:,: )! ef file by day for each time period
        REAL          , ALLOCATABLE, PUBLIC :: TEMPEF( : )    ! temporary holding array for efs
        CHARACTER     , ALLOCATABLE, PUBLIC :: EFTYPE( :,: )  ! ef file type (day, week, etc.) for each src
        INTEGER       , ALLOCATABLE, PUBLIC :: EFIDX ( :,: )  ! location of ef in file for each source
        LOGICAL       , ALLOCATABLE, PUBLIC :: USETIME( :,: ) ! true: time period is used

!.........  Variables for mapping MOBILE6-to-SMOKE vehicle types
        INTEGER, PUBLIC :: NVTYPE                        ! number of SMOKE vehicle types
        INTEGER, ALLOCATABLE, PUBLIC :: M6VEHMAP( : )    ! SMOKE vehicle type for each MOBILE6 type
        INTEGER, ALLOCATABLE, PUBLIC :: SMKVEH2EF( :,: ) ! SMOKE vehicle / emission process mapping
        END MODULE MODEMFAC
