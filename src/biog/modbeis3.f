
        MODULE MODBEIS3

C***********************************************************************
C  Module body starts at line 42 
C
C  DESCRIPTION:
C     This module contains the public variables and allocatable arrays 
C     used only in the biogenic emissions BEIS3v0.9 module.
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
C COPYRIGHT (C) 2001, MCNC--North Carolina Supercomputing Center
C All Rights Reserved
C
C See file COPYRIGHT for conditions of use.
C
C MCNC-Environmental Programs Group
C P.O. Box 12889
C Research Triangle Park, NC  27709-2889
C
C env_progs@mcnc.org
C
C Pathname: $Source$
C Last updated: $Date$ 
C
C***********************************************************************
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
        REAL,    ALLOCATABLE, PUBLIC :: AVGLAI( :, :, : )    ! avg leaf index

        CHARACTER*16, ALLOCATABLE, PUBLIC :: VEGID( : )     !  Veg types


        END MODULE MODBEIS3
