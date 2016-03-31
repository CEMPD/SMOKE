
        MODULE MODGRID

!***********************************************************************
!  Module body starts at line
!
!  DESCRIPTION:
!     This module contains the public variables and allocatable arrays 
!     used for gridding.
!
!  PRECONDITIONS REQUIRED:
!
!  SUBROUTINES AND FUNCTIONS CALLED:
!
!  REVISION HISTORY:
!     Created 10/2000 by M. Houyoux
!
!***************************************************************************
!
! Project Title: EDSS Tools Library
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

C...........   INCLUDES:
        INCLUDE 'IOPRVT3.EXT'   !  I/O API tools private string lengths

!.........  Horizontal grid information
        CHARACTER(IOVLEN3), PUBLIC :: GRDNM = ' '  ! grid name
        CHARACTER(IOVLEN3), PUBLIC :: COORD = ' '  ! coord system name
        INTEGER, PUBLIC :: GDTYP = -1     ! i/o api grid type code
        REAL(8), PUBLIC :: P_ALP = 0.D0   ! projection alpha
        REAL(8), PUBLIC :: P_BET = 0.D0   ! projection beta
        REAL(8), PUBLIC :: P_GAM = 0.D0   ! projection gamma
        REAL(8), PUBLIC :: XCENT = 0.D0   ! x-center of projection
        REAL(8), PUBLIC :: YCENT = 0.D0   ! y-center of projection
        REAL(8), PUBLIC :: XORIG = 0.D0   ! x-origin of grid
        REAL(8), PUBLIC :: YORIG = 0.D0   ! y-origin of grid
        REAL(8), PUBLIC :: XCELL = 0.D0   ! x-dim of cells
        REAL(8), PUBLIC :: YCELL = 0.D0   ! y-dim of cells
        INTEGER, PUBLIC :: NCOLS = 0      ! number of columns in grid
        INTEGER, PUBLIC :: NROWS = 0      ! number of rows in grid
        INTEGER, PUBLIC :: NGRID = 0      ! number of cells in grid
        INTEGER, PUBLIC :: XDIFF = 0      ! subgrid fewer cols (nxsub = nx - xdiff)
        INTEGER, PUBLIC :: YDIFF = 0      ! subgrid fewer cols (nysub = ny - ydiff)
        INTEGER, PUBLIC :: XOFF  = 0      ! subgrid offset (x-sub = x - xoff)
        INTEGER, PUBLIC :: YOFF  = 0      ! subgrid offset
        INTEGER, PUBLIC :: XOFF_A= 0      ! tmp subgrid offset (x-sub = x - xoff)
        INTEGER, PUBLIC :: YOFF_A= 0      ! tmp subgrid offset
        LOGICAL, PUBLIC :: OFFLAG = .FALSE. ! true: subgrid offset has been set

!.........  Vertical structure information
        INTEGER, PUBLIC :: NLAYS = 1       ! number of layers
        INTEGER, PUBLIC :: VGTYP  = -1     ! type of vertical coordinates
        REAL   , PUBLIC :: VGTOP  = 0.0    ! model-top, for sigma coord types
        REAL   , ALLOCATABLE, PUBLIC :: VGLVS( : ) ! vertical coordinate values

        END MODULE MODGRID
