
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

!.........  Source-category-specific variables

        INTEGER    , PUBLIC :: BYEAR    ! base inventory year
        INTEGER    , PUBLIC :: CATLEN   ! length of CATEGORY string
        INTEGER    , PUBLIC :: JSCC  =0 ! position in source chars of SCC (or 0)
        INTEGER    , PUBLIC :: JSTACK=0 ! position in source chars of stack
        INTEGER    , PUBLIC :: LSCCEND  ! end of left-SCC
        INTEGER    , PUBLIC :: MXCHRS   ! max no. of source characteristics
        INTEGER    , PUBLIC :: NCHARS   ! actual no. of source characteristics
        INTEGER    , PUBLIC :: NEMSFILE ! no. EMS-95 files
        INTEGER    , PUBLIC :: NIACT =0 ! no. unique activities in inventory
        INTEGER    , PUBLIC :: NIPOL =0 ! no. unique pollutants in inventory
        INTEGER    , PUBLIC :: NIPPA =0 ! NIACT + NIPPA
        INTEGER    , PUBLIC :: NIPAS =0 ! NIACT + NIPPA + NSPDAT
        INTEGER    , PUBLIC :: NPACT =0 ! number of variables per activity
        INTEGER    , PUBLIC :: NPPOL =0 ! number of variables per pollutant
        INTEGER    , PUBLIC :: NSPDAT=0 ! number of special data variables
        INTEGER    , PUBLIC :: NSRC  =0 ! number of SMOKE sources
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

        CHARACTER*1, PUBLIC :: CRL      = ' ' ! 'A', 'M', or 'P'
        CHARACTER*6, PUBLIC :: CATEGORY = ' ' ! 'AREA', 'MOBILE', or 'POINT'
        CHARACTER*6, PUBLIC :: CATDESC  = ' ' ! 'Area', 'Mobile', or 'Point'

!.........  Arrays of positions for source characteristics in full string of
!           source characteristics (dimension NCHARS)
        INTEGER    , ALLOCATABLE, PUBLIC :: SC_BEGP( : )
        INTEGER    , ALLOCATABLE, PUBLIC :: SC_ENDP( : )

!.........  Inventory pollutants and index to master list dimensioned by NIPOL
        INTEGER               , ALLOCATABLE, PUBLIC :: EIIDX( : )
        CHARACTER(LEN=IOVLEN3), ALLOCATABLE, PUBLIC :: EINAM( : )

!.........  Inventory activities and index to master list, dimensioned by NIACT
        INTEGER               , ALLOCATABLE, PUBLIC :: AVIDX ( : )
        CHARACTER(LEN=IOVLEN3), ALLOCATABLE, PUBLIC :: ACTVTY( : )

!.........  Inventory pollutants and inventory activies, dimensioned by NIPPA
        CHARACTER(LEN=IOVLEN3), ALLOCATABLE, PUBLIC :: EANAM ( : )
        CHARACTER(LEN=IOVLEN3), ALLOCATABLE, PUBLIC :: EAREAD( : ) ! for read3

!.........  Inventory pollutants/activies units and descriptions
        CHARACTER(LEN=IOULEN3), ALLOCATABLE, PUBLIC :: EAUNIT( : )
        CHARACTER(LEN=IODLEN3), ALLOCATABLE, PUBLIC :: EADESC( : )

!.........  Units conversions information
        INTEGER               , PUBLIC :: INVPIDX = 0    ! annual/O3 season idx
        REAL, ALLOCATABLE, PUBLIC :: EACNV( : )          ! units conv factors

!.........  Arrays for reading inventory file headers
        CHARACTER(LEN=IOVLEN3), ALLOCATABLE, PUBLIC :: TMPNAM( : )! var names
        INTEGER               , ALLOCATABLE, PUBLIC :: DATPOS( : )! var positn

!.........  Units for source attributes (such as stack parms)
        CHARACTER(LEN=IOULEN3), ALLOCATABLE, PUBLIC :: ATTRUNIT( : )

        END MODULE MODINFO
