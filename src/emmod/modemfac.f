        MODULE MODEMFAC

!***********************************************************************
!  Module body starts at line
!
!  DESCRIPTION:
!     This module contains the public allocatable arrays for using emission
!     factors with activity data
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

!.........  Sizes of the arrays...

        INTEGER, PUBLIC :: MXXACTV   ! max. no. activities in EFs xref
        INTEGER, PUBLIC :: MXXNPSI   ! max. no. PSIs in EFs xref
        INTEGER, PUBLIC :: MXPPGRP   ! max no. PSIs in a group of scenarios
        INTEGER, PUBLIC :: NDIU      ! no. diurnal EFs
        INTEGER, PUBLIC :: NNDI      ! no. non-diurnal EFs
        INTEGER, PUBLIC :: NMMTEF    ! no. ef-ref/temperature combinations

!.........  Unique lists of source characteristics and associated arrays...

!.........  List of unique PSIs for each activity of interest
        INTEGER, ALLOCATABLE, PUBLIC :: NPSI   ( : )   ! no. PSIs per activity
        INTEGER, ALLOCATABLE, PUBLIC :: PSILIST( :,: ) ! param scheme indices

!.........  Emission factor names - general arrays
        CHARACTER(LEN=IOVLEN3), ALLOCATABLE, PUBLIC :: DIUNAM( : )
        CHARACTER(LEN=IOVLEN3), ALLOCATABLE, PUBLIC :: NDINAM( : )

!.........  Emission factor index & temperature index combination arrays
        INTEGER, ALLOCATABLE, PUBLIC :: MMTEFPSI( : )  ! PSIs of combinations
        INTEGER, ALLOCATABLE, PUBLIC :: MMTEFIDX( : )  ! index to min/max tmprs
        INTEGER, ALLOCATABLE, PUBLIC :: MMTEFPTR( : )  ! pointer from PSIALL to
                                                       !  current indx MMTEFIDX
        INTEGER, ALLOCATABLE, PUBLIC :: MMTEFEND( : )  ! pntr to final tmpr indx
                                                       !  for PSI in MMTEFIDX

!.........  Emission factors. dim: temperatures, vehicle types, emis processes
        REAL, ALLOCATABLE, PUBLIC :: EFACDIU( :,:,: ) ! diurnal emission factors
        REAL, ALLOCATABLE, PUBLIC :: EFACNDI( :,:,: ) ! non-diurnal emis factors

!.........  Emission factors inputs data table
        INTEGER, PUBLIC :: NPSIDAT   ! no. PSIs in EF data file (MOBILE inputs)

        INTEGER, ALLOCATABLE, PUBLIC:: PSIDAT  ( : ) ! sorted PSIs
        INTEGER, ALLOCATABLE, PUBLIC:: PSIDATA ( : ) ! unsorted PSIs
        INTEGER, ALLOCATABLE, PUBLIC:: PDATINDX( : ) ! sorting index
        INTEGER, ALLOCATABLE, PUBLIC:: PDATLINE( : ) ! unsrt no. lines in file
        INTEGER, ALLOCATABLE, PUBLIC:: PDATMCNT( : ) ! unsrt no. scns in group
        INTEGER, ALLOCATABLE, PUBLIC:: PDATPNTR( : ) ! unsrt pointer
        INTEGER, ALLOCATABLE, PUBLIC:: PDATROOT( : ) ! unsrt root PSI of scenrio
        INTEGER, ALLOCATABLE, PUBLIC:: PDATTYPE( : ) ! unsrt 0=pure, >0=no. in combo

!.........  Combination emission factors table
        INTEGER, PUBLIC :: NCMBOPSI  ! no. combination PSIs
        INTEGER, PUBLIC :: MXNCOMBO  ! max no. contributing PSI per combo PSI

        INTEGER, ALLOCATABLE, PUBLIC :: CMBOPSI( :,: ) ! contributing PSIs
        INTEGER, ALLOCATABLE, PUBLIC :: CMBOFAC( :,: ) ! factors

        END MODULE MODEMFAC
