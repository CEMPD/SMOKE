
        MODULE MODLISTS

!***********************************************************************
!  Module body starts at line
!
!  DESCRIPTION:
!     This module contains the public allocatable arrays for various lists
!     from the inventory, such as for SCCs, SICs, inventory pollutants, 
!     and other criteria for selecting the parts of cross-reference files that
!     are relevant.
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
! COPYRIGHT (C) 2002, MCNC Environmental Modeling Center
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

        INTEGER, PUBLIC :: NINVSIC  ! no. unique SICs in inventory
        INTEGER, PUBLIC :: NINVSCC  ! no. unique SCCs in inventory
        INTEGER, PUBLIC :: NINVSCL  ! no. unique left SCCs in inventory
        INTEGER, PUBLIC :: NSCCPSIC ! no. all SCCs for all SICs
        INTEGER, PUBLIC :: NINVIFIP ! no. unique country/state/county codes
        INTEGER, PUBLIC :: NORISBLR ! no. unique ORIS // boilers
        INTEGER, PUBLIC :: NORISPNT ! no. unique ORIS // points
        INTEGER, PUBLIC :: NINVORIS ! no. unique ORIS 

!.........  Controllers
        LOGICAL, PUBLIC :: ORISFLAG  ! true: create ORIS-based arrays

!.........  Unique lists of source characteristics and associated arrays...

!.........  SIC arrays dimensioned by NINVSIC
        INTEGER, ALLOCATABLE, PUBLIC :: INVSIC ( : ) ! SICs
        INTEGER, ALLOCATABLE, PUBLIC :: IBEGSIC( : ) ! start position in SICSCC
        INTEGER, ALLOCATABLE, PUBLIC :: IENDSIC( : ) ! end   position in SICSCC

!.........  SCC arrays dimensioned by NINVSCC
        CHARACTER(LEN=SCCLEN3), ALLOCATABLE, PUBLIC :: INVSCC( : ) ! SCCs
        CHARACTER(LEN=SCCLEN3), ALLOCATABLE, PUBLIC :: INVSCL( : ) ! left SCCs
        CHARACTER(LEN=SDSLEN3), ALLOCATABLE, PUBLIC :: SCCDESC( : ) ! descrptn

!.........  SCCs for all SICs
        CHARACTER(LEN=SCCLEN3), ALLOCATABLE, PUBLIC :: SCCPSIC( : )

!.........  Country/state/county codes dimensioned by NINVFIP
        INTEGER, ALLOCATABLE, PUBLIC :: INVIFIP( : )

!.........  ORIS arrays
        INTEGER               , ALLOCATABLE, PUBLIC :: INVORFP( : ) ! FIPS for ORIS in inventory
        INTEGER               , ALLOCATABLE, PUBLIC :: OBSRCBG( : ) ! 1st source per ORIS/boiler
        INTEGER               , ALLOCATABLE, PUBLIC :: OBSRCNT( : ) ! source count per ORIS/boiler
        INTEGER               , ALLOCATABLE, PUBLIC :: OPSRCBG( : ) ! 1st source per ORIS/point
        INTEGER               , ALLOCATABLE, PUBLIC :: OPSRCNT( : ) ! source count per ORIS/point
        LOGICAL               , ALLOCATABLE, PUBLIC :: IORSMTCH( : ) ! true: inventory ORIS matched to CEM
        CHARACTER(LEN=ORSLEN3), ALLOCATABLE, PUBLIC :: INVORIS( : ) ! unique ORIS
        CHARACTER(LEN=DSCLEN3), ALLOCATABLE, PUBLIC :: INVODSC( : ) ! plant description from inventory
        CHARACTER(LEN=OBRLEN3), ALLOCATABLE, PUBLIC :: ORISBLR( : ) ! ORIS // boiler
        CHARACTER(LEN=OPTLEN3), ALLOCATABLE, PUBLIC :: ORISPNT( : ) ! ORIS // point

!.........  For valid pollutants and activities...

C.........  Full list of inventory pollutants/activities (in output order)
        INTEGER, PUBLIC :: MXIDAT = 0   ! Max no of inv pols & acvtys

        INTEGER, ALLOCATABLE, PUBLIC :: INVDCOD( : ) ! 5-digit pollutant/actvty code
        INTEGER, ALLOCATABLE, PUBLIC :: INVSTAT( : ) ! Status (<0 activity; >0 pol)

        REAL   , ALLOCATABLE, PUBLIC :: INVDCNV( : ) ! local conversion factor

        CHARACTER(LEN=IOVLEN3), ALLOCATABLE, PUBLIC :: INVDNAM( : ) ! name 
        CHARACTER(LEN=IOULEN3), ALLOCATABLE, PUBLIC :: INVDUNT( : ) ! units

        END MODULE MODLISTS
