        MODULE MODMET

!***********************************************************************
!  Module body starts at line
!
!  DESCRIPTION:
!     This module contains the derived meteorology data for applying emission
!     factors to activity data.
!
!  PRECONDITIONS REQUIRED:
!
!  SUBROUTINES AND FUNCTIONS CALLED:
!
!  REVISION HISTORY:
!     Created 6/99 by M. Houyoux
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

!...........   Setting for range of valid min/max temperatures
!        REAL   , PUBLIC :: MINT_MIN = 0.  ! min of range of minima
!        REAL   , PUBLIC :: MINT_MAX = 0.  ! max of range of minima
!        REAL   , PUBLIC :: MAXT_MIN = 0.  ! min of range of maxima
!        REAL   , PUBLIC :: MAXT_MAX = 0.  ! max of range of maxima
!        REAL   , PUBLIC :: TMMINVL  = 0.  ! temperature interval
!        REAL   , PUBLIC :: TMXINVL  = 0.  ! maximum temperature interval
!        INTEGER, PUBLIC :: NTMPR    = 0   ! number of valid single temperatures
!        INTEGER, PUBLIC :: NVLDTMM  = 0   ! number of valid combos
        REAL, PUBLIC :: MINTEMP = 0.   ! minimum temperature
        REAL, PUBLIC :: MAXTEMP = 0.   ! maximum temperature

!...........   Valid list of temperatures (dim: NTMPR)
!        REAL   , ALLOCATABLE, PUBLIC :: VLDTMPR( : )   ! valid temperatures

!...........   Valid list of min/max temperature combinations (dim: NVLDTMM)
!        REAL   , ALLOCATABLE, PUBLIC :: VLDTMIN( : )   ! valid mins
!        REAL   , ALLOCATABLE, PUBLIC :: VLDTMAX( : )   ! valid maxs

!...........   Daily min/max temperatures [K] (dim: NSRC)
        REAL   , ALLOCATABLE, PUBLIC :: TASRC   ( : )   ! per-source tmprs
!        REAL   , ALLOCATABLE, PUBLIC :: TKMAX   ( : )   ! working max
!        REAL   , ALLOCATABLE, PUBLIC :: TKMIN   ( : )   ! working min
!        REAL   , ALLOCATABLE, PUBLIC :: TKMAXOUT( :,: ) ! output  max
!        REAL   , ALLOCATABLE, PUBLIC :: TKMINOUT( :,: ) ! output  min
!        INTEGER, ALLOCATABLE, PUBLIC :: METIDX  ( :,: ) ! idx to master min/max

!...........   Hourly temperatures [K] (dim: NSRC)
!...              for Mobile5 processing, index 0 = 12 AM local time
!...              for Mobile6 processing, index 0 = 6 AM local time
        REAL   , ALLOCATABLE, PUBLIC :: TKHOUR  ( :,: ) ! temps by source per hour 
        REAL   , ALLOCATABLE, PUBLIC :: TKCOUNTY( :,: ) ! temps by county per hour

        END MODULE MODMET
