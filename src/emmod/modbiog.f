
        MODULE MODBIOG

!***********************************************************************
!  Module body starts at line 41
!
!  DESCRIPTION:
!     This module contains the public variables and allocatable arrays 
!     used only in the biogenic emissions module.
!
!  PRECONDITIONS REQUIRED:
!
!  SUBROUTINES AND FUNCTIONS CALLED:
!
!  REVISION HISTORY:
!     Created 3/00 by J. Vukovich
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
!*************************************************************************

        INCLUDE 'EMPRVT3.EXT'

C...........   Emission factor, vegetation types tables:

        INTEGER, PUBLIC ::    NVEG                     !  Number of veg types
        INTEGER, ALLOCATABLE, PUBLIC :: LAI  ( : )     !  Leaf area index
        REAL,    ALLOCATABLE, PUBLIC :: EMFAC( :,: )   !  Emission factors

        CHARACTER(LEN=BVGLEN3), ALLOCATABLE, PUBLIC :: VEGID( : )     !  Veg types

C.......   Gridded normalized emissions arrays

        REAL, ALLOCATABLE, PUBLIC ::  PINE( :, : )  !  pine forest
        REAL, ALLOCATABLE, PUBLIC ::  DECD( :, : )  !  deciduous forest
        REAL, ALLOCATABLE, PUBLIC ::  CONF( :, : )  !  other coniferous forest
        REAL, ALLOCATABLE, PUBLIC ::  AGRC( :, : )  !  grasslands
        REAL, ALLOCATABLE, PUBLIC ::  LEAF( :, : )  !  leaf area
        REAL, ALLOCATABLE, PUBLIC ::  OTHR( :, : )  !  other biogenic area

        REAL, ALLOCATABLE, PUBLIC ::  AVLAI(  : )   ! avg leaf area index

        REAL, ALLOCATABLE, PUBLIC ::  GRASNO( : )   ! Nitric oxide from grass
        REAL, ALLOCATABLE, PUBLIC ::  FORENO( : )   ! Nitric oxide from forest
        REAL, ALLOCATABLE, PUBLIC ::  WETLNO( : )   ! Nitric oxide from wetlands
        REAL, ALLOCATABLE, PUBLIC ::  AGRINO( : )   ! Nitric oxide from agricl


        END MODULE MODBIOG
