        MODULE MODSURG

!***********************************************************************
!
!  DESCRIPTION:
!     This module contains the public allocatable arrays for spatial
!     surrogates files.
!
!  PRECONDITIONS REQUIRED:
!
!  SUBROUTINES AND FUNCTIONS CALLED:
!
!  REVISION HISTORY:
!     
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

        INTEGER, PUBLIC :: NSRGREC      ! Number of surrogate file entries

        INTEGER, ALLOCATABLE, PUBLIC :: IDXSRGA( : )   ! Sorting index
        INTEGER, ALLOCATABLE, PUBLIC :: IDXSRGB( : )   ! Sorting index 
        INTEGER, ALLOCATABLE, PUBLIC :: SCELLA ( : )   ! Surrogate cell numbers
        INTEGER, ALLOCATABLE, PUBLIC :: SFIPSA ( : )   ! Surrogate Country/
                                                       ! State/County codes
        INTEGER, ALLOCATABLE, PUBLIC :: SSRGIDA( : )   ! Surrogate ID codes
        REAL   , ALLOCATABLE, PUBLIC :: SFRACA ( : )   ! Surrogate fractions

C...........   Surrogate-cell::FIPS table
        
        INTEGER,              PUBLIC :: NSRGFIPS     ! no. srgs county codes
        INTEGER,              PUBLIC :: NSRGS        ! no. surrogate ID codes

        INTEGER, ALLOCATABLE, PUBLIC :: NCELLS ( : ) ! no. cells per county code
        INTEGER, ALLOCATABLE, PUBLIC :: SRGFIPS( : ) ! list of cnty codes
        INTEGER, ALLOCATABLE, PUBLIC :: SRGLIST( : ) ! list of surrgate codes
        INTEGER, ALLOCATABLE, PUBLIC :: FIPCELL( :, : )    ! cell numbers
        REAL   , ALLOCATABLE, PUBLIC :: SRGFRAC( :, :, : ) ! surrogate fractions
        REAL   , ALLOCATABLE, PUBLIC :: SRGCSUM( :, : )    ! surg sum by county

        END MODULE MODSURG
