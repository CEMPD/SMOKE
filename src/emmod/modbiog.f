
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
!*************************************************************************

        IMPLICIT NONE

        INCLUDE 'EMPRVT3.EXT'

C...........   Emission factor, vegetation types tables:

        INTEGER, PUBLIC ::    NVEG                     !  Number of veg types
        INTEGER, ALLOCATABLE, PUBLIC :: LAI  ( : )     !  Leaf area index
        REAL,    ALLOCATABLE, PUBLIC :: EMFAC( :,: )   !  Emission factors

        CHARACTER(BVGLEN3), ALLOCATABLE, PUBLIC :: VEGID( : )     !  Veg types

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
