        MODULE MODAR2PT

!***********************************************************************
!  Module body starts at line 40
!
!  DESCRIPTION:
!     This module contains the public allocatable arrays for the area-to-point
!     adjustment factors
!
!  PRECONDITIONS REQUIRED:
!
!  SUBROUTINES AND FUNCTIONS CALLED:
!
!  REVISION HISTORY:
!     Created 11/02 by M. Houyoux
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

        INCLUDE 'EMPRVT3.EXT'   !  private emissions string widths parameters

!.........  Define types needed for module
        TYPE :: AR2PT
            SEQUENCE
            CHARACTER(FIPLEN3) FIP   ! country, state, and county code
            REAL          LAT   ! latitude
            REAL          LON   ! longitude
            REAL          ALLOC ! allocation factor
            CHARACTER(25) NAME  ! airport name
        END TYPE

        TYPE :: A2PREPTYPE
            INTEGER            STATE     ! state code
            CHARACTER(SCCLEN3) SCC       ! SCC code
            INTEGER            POLL      ! pollutant code
            INTEGER            NFIPS     ! number of FIPS codes
            REAL               ORIGEMIS  ! emissions before processing
            REAL               SUMEMIS   ! summed emissions after
        END TYPE

!.........  Area-to-point table. Second dimension is the number of tables (and is
!           used to dimention the NAR2PT array). First dimension is the
!           maximum number of rows in any table.
        INTEGER                   , PUBLIC :: NTABLA2P        ! number of tables
        INTEGER                   , PUBLIC :: MXROWA2P        ! max. number of rows
        INTEGER      , ALLOCATABLE, PUBLIC :: NAR2PT  ( : )   ! no. entries in each table
        TYPE( AR2PT ), ALLOCATABLE, PUBLIC :: AR2PTABL( :,: ) ! area-to-point table
        
!.........  Area-to-point table SCCs.
        INTEGER                        , PUBLIC :: NA2PSCC         ! no. of area-to-point SCCs
        CHARACTER(SCCLEN3), ALLOCATABLE, PUBLIC :: A2PSCC( : )  ! area-to-point SCCs

!.........  Reporting array (dimension is NCONDSRC)
        INTEGER                        , PUBLIC :: NCONDSRC        ! no. of condensed srcs
        TYPE( A2PREPTYPE ), ALLOCATABLE, PUBLIC :: REPAR2PT( : )   ! area-to-point report

        END MODULE MODAR2PT
