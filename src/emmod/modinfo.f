
        MODULE MODINFO

!***********************************************************************
!  Module body starts at line
!
!  DESCRIPTION:
!     This module contains the public data that is source-category-specific
!     such as the CATEGORY variable, the left and right sizes of the SCC
!     the number of sources.  This might not be useful for programs that 
!     combine multiple source categories.
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

!.........  Source-category-specific variables

        INTEGER    , PUBLIC :: BYEAR    ! base inventory year
        INTEGER    , PUBLIC :: INVPIDX=0! annual/average day idx
        INTEGER    , PUBLIC :: INV_MON  ! base inventory month 
        INTEGER    , PUBLIC :: CATLEN   ! length of CATEGORY string
        INTEGER    , PUBLIC :: JSCC  =0 ! position in source chars of SCC (or 0)
        INTEGER    , PUBLIC :: JSTACK=0 ! position in source chars of stack
        INTEGER    , PUBLIC :: LSCCEND  ! end of left-SCC
        INTEGER    , PUBLIC :: MXCHRS   ! max no. of source characteristics
        INTEGER    , PUBLIC :: NCHARS   ! actual no. of source characteristics
        INTEGER    , PUBLIC :: NIACT =0 ! no. unique activities in inventory
        INTEGER    , PUBLIC :: NIPOL =0 ! no. unique pollutants in inventory
        INTEGER    , PUBLIC :: NIPPA =0 ! NIACT + NIPOL
        INTEGER    , PUBLIC :: NIPAS =0 ! NIACT + NIPOL + NSPDAT
        INTEGER    , PUBLIC :: NPACT =0 ! number of variables per activity
        INTEGER    , PUBLIC :: NPPOL =0 ! number of variables per pollutant
        INTEGER    , PUBLIC :: NSPDAT=0 ! number of special data variables
        INTEGER    , PUBLIC :: NSRC  =0 ! number of SMOKE sources
        INTEGER    , PUBLIC :: PLTIDX   ! index of plant (if any) in SC_BEGP
        INTEGER    , PUBLIC :: RSCCBEG  ! beginning of right-SCC
        INTEGER    , PUBLIC :: SCCLEV1  ! right-most position of level-1 of SCC
        INTEGER    , PUBLIC :: SCCLEV2  ! right-most position of level-2 of SCC
        INTEGER    , PUBLIC :: SCCLEV3  ! right-most position of level-3 of SCC
        INTEGER    , PUBLIC :: SCCLEV4  ! right-most position of level-4 of SCC

!.........  Positions of pollutant-specific inventory data in storage array
        INTEGER :: NC1 = 0 !  pos in 2nd dim of POLVLA of primary control code
        INTEGER :: NC2 = 0 !  pos in 2nd dim of POLVLA of secondary cntrl code
        INTEGER :: NCE = 0 !  pos in 2nd dim of POLVLA of control efficiency
        INTEGER :: NEF = 0 !  pos in 2nd dim of POLVLA of emission factors
        INTEGER :: NEM = 0 !  pos in 2nd dim of POLVLA of annual emissions
        INTEGER :: NDY = 0 !  pos in 2nd dim of POLVLA of average day emis
        INTEGER :: NRE = 0 !  pos in 2nd dim of POLVLA of rule effectivenss
        INTEGER :: NRP = 0 !  pos in 2nd dim of POLVLA of rule penetration
        PUBLIC NC1, NC2, NCE, NEF, NEM, NDY, NRE, NRP

        CHARACTER,    PUBLIC :: CRL      = ' ' ! 'A', 'M', or 'P'
        CHARACTER(6), PUBLIC :: CATEGORY = ' ' ! 'AREA', 'MOBILE', or 'POINT'
        CHARACTER(6), PUBLIC :: CATDESC  = ' ' ! 'Area', 'Mobile', or 'Point'

!.........  Allocatable/local variables for formulas
        CHARACTER(512),     PUBLIC      :: VAR_FORMULA = ' ' ! SMKINVEN_FORMULA
        CHARACTER(IOULEN3), PUBLIC      :: NETCDFUNIT        ! NetCDF poll unit
        INTEGER,            PUBLIC      :: NCOMP = 0 ! no of Smkinven formuals  
        LOGICAL,            ALLOCATABLE :: CHKPLUS ( : ) ! true: formula uses a + sign
        LOGICAL,            ALLOCATABLE :: CHKMINUS( : ) ! true: formula uses a - sign
        CHARACTER(128),     ALLOCATABLE :: FORMULAS( : ) ! formulas for emissions
        CHARACTER(IOVLEN3), ALLOCATABLE :: VIN_A( : )    ! first variable in equation
        CHARACTER(IOVLEN3), ALLOCATABLE :: VIN_B( : )    ! second variable in equation
        CHARACTER(IOVLEN3), ALLOCATABLE :: VNAME( : )    ! computed variable name
        PUBLIC CHKPLUS, CHKMINUS, FORMULAS, VIN_A, VIN_B, VNAME

!.........  Arrays of positions for source characteristics in full string of
!           source characteristics (dimension NCHARS)
        INTEGER    , ALLOCATABLE, PUBLIC :: SC_BEGP( : )
        INTEGER    , ALLOCATABLE, PUBLIC :: SC_ENDP( : )

!.........  Inventory pollutants and index to master list dimensioned by NIPOL
        INTEGER,            ALLOCATABLE, PUBLIC :: EIIDX( : )
        CHARACTER(IOVLEN3), ALLOCATABLE, PUBLIC :: EINAM( : )

!.........  Inventory activities and index to master list, dimensioned by NIACT
        INTEGER,            ALLOCATABLE, PUBLIC :: AVIDX ( : )
        CHARACTER(IOVLEN3), ALLOCATABLE, PUBLIC :: ACTVTY( : )

!.........  Inventory pollutants and inventory activies, dimensioned by NIPPA
        CHARACTER(IOVLEN3), ALLOCATABLE, PUBLIC :: EANAM ( : )
        CHARACTER(IOVLEN3), ALLOCATABLE, PUBLIC :: EAREAD( : ) ! for read3

!.........  Inventory pollutants/activies units and descriptions
        CHARACTER(IOULEN3), ALLOCATABLE, PUBLIC :: EAUNIT( : )
        CHARACTER(IODLEN3), ALLOCATABLE, PUBLIC :: EADESC( : )

!.........  Map file pollutants list and physical file names
        INTEGER,                         PUBLIC :: NMAP = 0
        CHARACTER(IOVLEN3), ALLOCATABLE, PUBLIC :: MAPNAM( : )
        CHARACTER(PHYLEN3), ALLOCATABLE, PUBLIC :: MAPFIL( : )

!.........  Units conversions information
        REAL, ALLOCATABLE, PUBLIC :: EACNV( : )     ! units conv factors

!.........  Arrays for reading inventory file headers
        INTEGER,                         PUBLIC :: NPOLID     ! POLID original no of pollutants
        CHARACTER(IOVLEN3), ALLOCATABLE, PUBLIC :: TMPNAM( : )! var names
        INTEGER,            ALLOCATABLE, PUBLIC :: DATPOS( : )! var positn

!.........  Units for source attributes (such as stack parms)
        CHARACTER(IOULEN3), ALLOCATABLE, PUBLIC :: ATTRUNIT( : )

        END MODULE MODINFO
