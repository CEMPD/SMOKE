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
!****************************************************************************

        IMPLICIT NONE
        
        INCLUDE 'EMPRVT3.EXT'

        INTEGER, PUBLIC :: NSRGREC      ! Number of surrogate file entries

        INTEGER, ALLOCATABLE, PUBLIC :: IDXSRGA( : )   ! Sorting index
        INTEGER, ALLOCATABLE, PUBLIC :: IDXSRGB( : )   ! Sorting index 
        INTEGER, ALLOCATABLE, PUBLIC :: SCELLA ( : )   ! Surrogate cell numbers
        CHARACTER(FIPLEN3), ALLOCATABLE, PUBLIC :: SFIPSA ( : )   ! Surrogate Country/
                                                       ! State/County codes
        INTEGER, ALLOCATABLE, PUBLIC :: SSRGIDA( : )   ! Surrogate ID codes
        REAL   , ALLOCATABLE, PUBLIC :: SFRACA ( : )   ! Surrogate fractions

C...........   Surrogate-cell::FIPS table
        
        INTEGER,              PUBLIC :: NSRGFIPS     ! no. srgs county codes
        INTEGER,              PUBLIC :: NSRGS        ! no. surrogate ID codes
        INTEGER,              PUBLIC :: MXCFIP       ! Max cells per county code
        INTEGER,              PUBLIC :: SRGNROWS     ! no. rows in surrogates file
        INTEGER,              PUBLIC :: SRGNCOLS     ! no. cols in surrogates file
        CHARACTER(16),        PUBLIC :: SRGFMT       !  surrogates format

        INTEGER, ALLOCATABLE, PUBLIC :: NTLINES( : )       ! no. line buffer for each surrogate
        INTEGER, ALLOCATABLE, PUBLIC :: NCELLS ( : )       ! no. cells per county code
        CHARACTER(FIPLEN3), ALLOCATABLE, PUBLIC :: SRGFIPS( : ) ! list of cnty codes
        INTEGER, ALLOCATABLE, PUBLIC :: SRGLIST( : )       ! list of surrgate codes
        INTEGER, ALLOCATABLE, PUBLIC :: FIPCELL( :, : )    ! cell numbers
        REAL   , ALLOCATABLE, PUBLIC :: SRGFRAC( :, :, : ) ! surrogate fractions
        REAL   , ALLOCATABLE, PUBLIC :: SRGCSUM( :, : )    ! surg sum by county

C...........   Surrogate description files : SRGDESC

        INTEGER,                      PUBLIC :: NTSRGDSC     ! total no of surrogate files in SRGDESC
        INTEGER        , ALLOCATABLE, PUBLIC :: SRGFCOD( : ) ! surrgate code
        CHARACTER( 60 ), ALLOCATABLE, PUBLIC :: SRGFREG( : ) ! surrgate region
        CHARACTER( 60 ), ALLOCATABLE, PUBLIC :: SRGFDES( : ) ! surrgate description
        CHARACTER( 60 ), ALLOCATABLE, PUBLIC :: SRGFNAM( : ) ! surrgate file name

        END MODULE MODSURG
