
        MODULE MODBEIS3

C***********************************************************************
C  Module body starts at line 42 
C
C  DESCRIPTION:
C     This module contains the public variables and allocatable arrays 
C     used in calculating biogenic emissions using BEIS v3.09 or v3.12.
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION HISTORY:
C     03/01: prototype by Jeff Vukovich
C
C***************************************************************************
C
C Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
C                System
C File: @(#)$Id$
C
C COPYRIGHT (C) 2004, Environmental Modeling for Policy Development
C All Rights Reserved
C 
C Carolina Environmental Program
C University of North Carolina at Chapel Hill
C 137 E. Franklin St., CB# 6116
C Chapel Hill, NC 27599-6116
C 
C smoke@unc.edu
C
C Pathname: $Source$
C Last updated: $Date$ 
C
C***********************************************************************

        IMPLICIT NONE

C...........   Emission factor, vegetation types tables:

        INTEGER, PUBLIC ::    NVEG                     !  Number of veg types
        INTEGER, ALLOCATABLE, PUBLIC :: LAI( : )       !  Leaf area index
        REAL,    ALLOCATABLE, PUBLIC :: EMFAC( :, : )  !  Emission factors
        REAL,    ALLOCATABLE, PUBLIC :: LFBIO( : )     !  Dry leaf biomass
        REAL,    ALLOCATABLE, PUBLIC :: WFAC( : )      !  Winter biomass factor
        REAL,    ALLOCATABLE, PUBLIC :: SLW( : )       !  Specific leaf weight

        REAL,    ALLOCATABLE, PUBLIC :: AVGISOP( :, :, : )   ! avg isoprene
        REAL,    ALLOCATABLE, PUBLIC :: AVGOVOC( :, :, : )   ! avg other VOCs
        REAL,    ALLOCATABLE, PUBLIC :: AVGMONO( :, :, : )   ! avg monoterpenes
        REAL,    ALLOCATABLE, PUBLIC :: AVGNO( :, :, : )     ! avg nitric oxide
        REAL,    ALLOCATABLE, PUBLIC :: AVGLAI( :, :, :, : ) ! avg leaf index

        REAL,    ALLOCATABLE, PUBLIC :: AVGEMIS( :, :, :, : )  ! avg emissions (3.12)
        REAL,    ALLOCATABLE, PUBLIC ::  NOEMIS( :, :, : )     ! NO emissions (3.12)

        CHARACTER(16), ALLOCATABLE, PUBLIC :: VEGID( : )     !  Veg types

        END MODULE MODBEIS3
