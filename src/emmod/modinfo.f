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

!.........  Source-category-specific variables

        INTEGER    , PUBLIC :: BYEAR    ! base inventory year
        INTEGER    , PUBLIC :: CATLEN   ! length of CATEGORY string
        INTEGER    , PUBLIC :: JSCC     ! position in source chars of SCC (or 0)
        INTEGER    , PUBLIC :: LSCCEND  ! end of left-SCC
        INTEGER    , PUBLIC :: MXCHRS   ! max no. of source characteristics
        INTEGER    , PUBLIC :: NCHARS   ! actual no. of source characteristics
        INTEGER    , PUBLIC :: NEMSFILE ! no. EMS-95 files
        INTEGER    , PUBLIC :: NIACT    ! no. unique activities in inventory
        INTEGER    , PUBLIC :: NIPOL    ! no. unique pollutants in inventory
        INTEGER    , PUBLIC :: NIPPA    ! NIACT + NIPPA
        INTEGER    , PUBLIC :: NPACT    ! number of variables per activity
        INTEGER    , PUBLIC :: NPPOL    ! number of variables per pollutant
        INTEGER    , PUBLIC :: NSRC     ! number of SMOKE sources
        INTEGER    , PUBLIC :: PLTIDX   ! index of plant (if any) in SC_BEGP
        INTEGER    , PUBLIC :: RSCCBEG  ! beginning of right-SCC

!.........  Positions of pollutant-specific inventory data in storage array
        INTEGER :: NC1 = 0 !  pos in 2nd dim of POLVLA of primary control code
        INTEGER :: NC2 = 0 !  pos in 2nd dim of POLVLA of secondary cntrl code
        INTEGER :: NCE = 0 !  pos in 2nd dim of POLVLA of control efficiency
        INTEGER :: NEF = 0 !  pos in 2nd dim of POLVLA of emission factors
        INTEGER :: NEM = 0 !  pos in 2nd dim of POLVLA of annual emissions
        INTEGER :: NOZ = 0 !  pos in 2nd dim of POLVLA of ozone season emis
        INTEGER :: NRE = 0 !  pos in 2nd dim of POLVLA of rule effectivenss
        INTEGER :: NRP = 0 !  pos in 2nd dim of POLVLA of rule penetration
        PUBLIC NC1, NC2, NCE, NEF, NEM, NOZ, NRE, NRP

        CHARACTER*1, PUBLIC :: CRL      ! 'A', 'M', or 'P'
        CHARACTER*6, PUBLIC :: CATEGORY ! 'AREA', 'MOBILE', or 'POINT'
        CHARACTER*6, PUBLIC :: CATDESC  ! 'Area', 'Mobile', or 'Point'

!.........  Arrays of positions for source characteristics in full string of
!           source characteristics (dimension NCHARS)
        INTEGER    , ALLOCATABLE, PUBLIC :: SC_BEGP( : )
        INTEGER    , ALLOCATABLE, PUBLIC :: SC_ENDP( : )

!.........  Inventory pollutants dimensioned by NIPOL
        CHARACTER(LEN=IOVLEN3), ALLOCATABLE, PUBLIC :: EINAM( : )

!.........  Inventory activities, dimensioned by NIACT
        CHARACTER(LEN=IOVLEN3), ALLOCATABLE, PUBLIC :: ACTVTY( : )

!.........  Inventory pollutants and inventory activies, dimensioned by NIPPA
        CHARACTER(LEN=IOVLEN3), ALLOCATABLE, PUBLIC :: EANAM( : )

!.........  Mobile source characteristics tables
        INTEGER                                     :: NVTYPE ! no. veh types
        INTEGER               , ALLOCATABLE, PUBLIC :: IVTIDLST( : ) ! IDs
        CHARACTER(LEN=VTPLEN3), ALLOCATABLE, PUBLIC :: CVTYPLST( : ) ! names

        INTEGER                                     :: NRCLAS ! no. road classes
        INTEGER               , ALLOCATABLE, PUBLIC :: AMSRDCLS( : ) ! AIRS rclas
        INTEGER               , ALLOCATABLE, PUBLIC :: RDWAYTYP( : ) ! roadway type

        END MODULE MODINFO
