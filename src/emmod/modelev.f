
        MODULE MODELEV

!***********************************************************************
!  Module body starts at line 42
!
!  DESCRIPTION:
!     This module contains the public arrays for processing major point
!     sources and plume-in-grid point sources
!
!  PRECONDITIONS REQUIRED:
!
!  SUBROUTINES AND FUNCTIONS CALLED:
!
!  REVISION HISTORY:
!     Created 8/99 by M. Houyoux
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

!...........   Number of stack groups
        INTEGER, PUBLIC :: NGROUP

!...........   Allocatable arrays for specifying source type
        INTEGER, ALLOCATABLE, PUBLIC:: GROUPID ( : ) ! stack group ID
        REAL   , ALLOCATABLE, PUBLIC:: ELEVFLTR( : ) ! =0. for lower, =1. elev
        LOGICAL, ALLOCATABLE, PUBLIC:: LMAJOR  ( : ) ! true: src is a major src
        LOGICAL, ALLOCATABLE, PUBLIC:: LPING   ( : ) ! true: src is a PinG src

!...........   Allocatable arrays for stack groups
        INTEGER, ALLOCATABLE, PUBLIC:: GRPGID( : ) ! sorted stack group ID
        INTEGER, ALLOCATABLE, PUBLIC:: GRPGIDA( : )! unsorted stack group ID
        INTEGER, ALLOCATABLE, PUBLIC:: GRPCNT( : ) ! no. stacks per group
        INTEGER, ALLOCATABLE, PUBLIC:: GRPCOL( : ) ! col number
        INTEGER, ALLOCATABLE, PUBLIC:: GRPROW( : ) ! row number
        INTEGER, ALLOCATABLE, PUBLIC:: GRPIDX( : ) ! sorting index for GRPGIDA

        REAL   , ALLOCATABLE, PUBLIC:: GRPDM ( : ) ! group intl stack diam [m]
        REAL   , ALLOCATABLE, PUBLIC:: GRPFL ( : ) ! group exit flw rate [m^3/s]
        REAL   , ALLOCATABLE, PUBLIC:: GRPHT ( : ) ! group stack height [m]
        REAL   , ALLOCATABLE, PUBLIC:: GRPLAT( : ) ! group latitude [degrees]
        REAL   , ALLOCATABLE, PUBLIC:: GRPLON( : ) ! group longitude [degrees]
        REAL   , ALLOCATABLE, PUBLIC:: GRPTK ( : ) ! group exit temperature [K]
        REAL   , ALLOCATABLE, PUBLIC:: GRPVE ( : ) ! group exit velocity [m/s]
        REAL   , ALLOCATABLE, PUBLIC:: GRPXL ( : ) ! x-location given projection
        REAL   , ALLOCATABLE, PUBLIC:: GRPYL ( : ) ! y-location given projection

!...........   Allocatable arrays for major and PinG sources indexing
        INTEGER, ALLOCATABLE, PUBLIC:: ELEVSIDX( : ) ! Elev source -> all srcs
        INTEGER, ALLOCATABLE, PUBLIC:: PINGGIDX( : ) ! PinG source -> PinG group

!.........  Elevated source and plume-in-grid emissions arrays
        REAL   , ALLOCATABLE, PUBLIC :: PGRPEMIS( : ) ! PinG group emissions
        REAL   , ALLOCATABLE, PUBLIC :: ELEVEMIS( : ) ! Major elev srcs emis

        END MODULE MODELEV
