
        MODULE MODTMPRL

!***********************************************************************
!  Module body starts at line 41
!
!  DESCRIPTION:
!     This module contains temporal allocation information. 
!
!  PRECONDITIONS REQUIRED:
!
!  SUBROUTINES AND FUNCTIONS CALLED:
!
!  REVISION HISTORY:
!     Created 1/99 by M. Houyoux for temporal profiles only
!     Changed 10/2000 by MRH to include additional temporal info.
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

!.........  Holiday dates arrays
        INTEGER, PUBLIC :: NHOLIDAY  ! number of holidays

        INTEGER, ALLOCATABLE, PUBLIC :: HOLREGN ( : ) ! Region code of holiday
        INTEGER, ALLOCATABLE, PUBLIC :: HOLJDATE( : ) ! Julian date of holidays
	INTEGER, ALLOCATABLE, PUBLIC :: HOLALTDY( : ) ! alternative day of week

!.........  Sorted temporal profiles
        INTEGER, PUBLIC :: NMON   ! number of monthly profiles
        INTEGER, PUBLIC :: NWEK   ! number of weekly profiles
        INTEGER, PUBLIC :: NHRL   ! number of diurnal profiles 

        INTEGER, ALLOCATABLE, PUBLIC :: MONREF( : )   ! Monthly codes
        INTEGER, ALLOCATABLE, PUBLIC :: WEKREF( : )   ! Weekly codes 
        INTEGER, ALLOCATABLE, PUBLIC :: HRLREF( : )   ! Diurnal codes

        REAL   , ALLOCATABLE, PUBLIC :: MONFAC( :,: ) ! Monthly factors
        REAL   , ALLOCATABLE, PUBLIC :: WEKFAC( :,: ) ! Weekly facs (week-norm)
        REAL   , ALLOCATABLE, PUBLIC :: XWKFAC( :,: ) ! Weekly facs (wkday-norm)
        REAL   , ALLOCATABLE, PUBLIC :: HRLFAC( :,:,: ) ! Hourly factors

        END MODULE MODTMPRL
