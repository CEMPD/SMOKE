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
! COPYRIGHT (C) 1998, MCNC--North Carolina Supercomputing Center
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
        INTEGER, PUBLIC :: NINVSCL  ! no. unique left-SCCs in inventory

!.........  Unique lists of source characteristics and associated arrays...

!.........  SIC arrays dimensioned by NINVSIC
        INTEGER, ALLOCATABLE, PUBLIC :: INVSIC ( : ) ! SICs
        INTEGER, ALLOCATABLE, PUBLIC :: IBEGSIC( : ) ! start position in INVSCC
        INTEGER, ALLOCATABLE, PUBLIC :: IENDSIC( : ) ! end   position in INVSCC

!.........  SCC arrays dimensioned by NINVSCC
        CHARACTER(LEN=SCCLEN3), ALLOCATABLE, PUBLIC :: INVSCC( : ) ! SCCs
        CHARACTER(LEN=SCCLEN3), ALLOCATABLE, PUBLIC :: INVSCL( : ) ! left SCCs

        END MODULE MODLISTS
