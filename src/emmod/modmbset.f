
        MODULE MODMBSET

!***********************************************************************
!  Module body starts at line 42
!
!  DESCRIPTION:
!     This module contains the public allocatable arrays used when
!     setting up MOBILE6 information.
!
!  PRECONDITIONS REQUIRED:
!
!  SUBROUTINES AND FUNCTIONS CALLED:
!
!  REVISION HISTORY:
!     Created 10/01 by C. Seppanen
!
!***************************************************************************
!
! Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
!                System
! File: @(#)$Id$
!
! COPYRIGHT (C) 2002, MCNC Environmental Modeling Center
C All Rights Reserved
C
C See file COPYRIGHT for conditions of use.
C
C Environmental Modeling Center
C MCNC
C P.O. Box 12889
C Research Triangle Park, NC  27709-2889
C
C smoke@emc.mcnc.org
!
! Pathname: $Source$
! Last updated: $Date$ 
!
!****************************************************************************

        IMPLICIT NONE

!.........  Arrays for storing file information
        INTEGER, ALLOCATABLE, PUBLIC :: MCREFSORT ( :,: ) ! sorted MCREF data  
        INTEGER, ALLOCATABLE, PUBLIC :: MCREFIDX ( :,: )  ! index into MCREF by ref. county
        
        INTEGER, ALLOCATABLE, PUBLIC :: MVREFSORT ( :,: ) ! sorted MVREF data
        
        INTEGER, PUBLIC :: NREFC           ! no. of reference counties
        INTEGER, PUBLIC :: NINVC           ! no. of unique counties in inventory
        
        INTEGER, PUBLIC :: NREFFLAGS = 3   ! no. of settings flags in MVREF file
        
!.........  Time period constants
        INTEGER, PARAMETER :: DAILY    = 1    ! daily temperature profiles
        INTEGER, PARAMETER :: WEEKLY   = 2    ! weekly temperature averaging
        INTEGER, PARAMETER :: MONTHLY  = 3    ! monthly temperature averaging
        INTEGER, PARAMETER :: EPISLEN  = 4    ! episode length temp. averaging
        
        END MODULE MODMBSET
        
