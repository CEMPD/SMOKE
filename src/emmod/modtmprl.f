        MODULE MODTPRO

!***********************************************************************
!  Module body starts at line 41
!
!  DESCRIPTION:
!     This module contains the public allocatable arrays for temporal profile
!     tables.
!
!  PRECONDITIONS REQUIRED:
!
!  SUBROUTINES AND FUNCTIONS CALLED:
!
!  REVISION HISTORY:
!     Created 1/99 by M. Houyoux
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

!.........  Sorted temporal profiles
        INTEGER, PUBLIC :: NMON   ! number of monthly profiles
        INTEGER, PUBLIC :: NWEK   ! number of weekly profiles
        INTEGER, PUBLIC :: NWKD   ! number of weekday diurnal profiles 
        INTEGER, PUBLIC :: NEND   ! number of weekend diurnal profiles 

        INTEGER, ALLOCATABLE, PUBLIC :: MONREF( : )   ! Monthly codes
        INTEGER, ALLOCATABLE, PUBLIC :: WEKREF( : )   ! Weekly codes 
        INTEGER, ALLOCATABLE, PUBLIC :: WKDREF( : )   ! Weekday-diurnal codes
        INTEGER, ALLOCATABLE, PUBLIC :: ENDREF( : )   ! Weekend-diurnal codes

        REAL   , ALLOCATABLE, PUBLIC :: MONFAC( :,: ) ! Monthly factors
        REAL   , ALLOCATABLE, PUBLIC :: WEKFAC( :,: ) ! Weekly facs (week-norm)
        REAL   , ALLOCATABLE, PUBLIC :: XWKFAC( :,: ) ! Weekly facs (wkday-norm)
        REAL   , ALLOCATABLE, PUBLIC :: WKDFAC( :,: ) ! Weekday-diurnal factors
        REAL   , ALLOCATABLE, PUBLIC :: ENDFAC( :,: ) ! Weekend-diurnal factors        

        END MODULE MODTPRO
