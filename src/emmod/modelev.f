
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

!...........   Number of stack groups
        INTEGER, PUBLIC :: NGROUP

!...........   Allocatable arrays for specifying source type
        INTEGER, ALLOCATABLE, PUBLIC:: GROUPID ( : ) ! stack group ID
        INTEGER, ALLOCATABLE, PUBLIC:: GINDEX  ( : ) ! stack group sorting index
        REAL   , ALLOCATABLE, PUBLIC:: ELEVFLTR( : ) ! =0. for lower, =1. elev
        LOGICAL, ALLOCATABLE, PUBLIC:: LMAJOR  ( : ) ! true: src is a major src
        LOGICAL, ALLOCATABLE, PUBLIC:: LPING   ( : ) ! true: src is a PinG src
        REAL*8 , ALLOCATABLE, PUBLIC:: SRCXL   ( : ) ! x-location given projection
        REAL*8 , ALLOCATABLE, PUBLIC:: SRCYL   ( : ) ! y-location given projection

!...........   Allocatable arrays for stack groups
        INTEGER, ALLOCATABLE, PUBLIC:: GRPGID( : ) ! sorted stack group ID
        INTEGER, ALLOCATABLE, PUBLIC:: GRPGIDA( : )! unsorted stack group ID
        INTEGER, ALLOCATABLE, PUBLIC:: GRPCNT( : ) ! no. stacks per group
        INTEGER, ALLOCATABLE, PUBLIC:: GRPCOL( : ) ! col number
        INTEGER, ALLOCATABLE, PUBLIC:: GRPROW( : ) ! row number
        INTEGER, ALLOCATABLE, PUBLIC:: GRPIDX( : ) ! sorting index for GRPGIDA
        CHARACTER(FIPLEN3), ALLOCATABLE, PUBLIC:: GRPFIP( : ) ! fips code of stack        

        REAL   , ALLOCATABLE, PUBLIC:: GRPDM ( : ) ! group intl stack diam [m]
        REAL   , ALLOCATABLE, PUBLIC:: GRPFL ( : ) ! group exit flw rate [m^3/s]
        REAL   , ALLOCATABLE, PUBLIC:: GRPHT ( : ) ! group stack height [m]
        REAL*8 , ALLOCATABLE, PUBLIC:: GRPLAT( : ) ! group latitude [degrees]
        REAL*8 , ALLOCATABLE, PUBLIC:: GRPLON( : ) ! group longitude [degrees]
        REAL   , ALLOCATABLE, PUBLIC:: GRPTK ( : ) ! group exit temperature [K]
        REAL   , ALLOCATABLE, PUBLIC:: GRPVE ( : ) ! group exit velocity [m/s]
        REAL*8 , ALLOCATABLE, PUBLIC:: GRPXL ( : ) ! x-location given projection
        REAL*8 , ALLOCATABLE, PUBLIC:: GRPYL ( : ) ! y-location given projection
        REAL   , ALLOCATABLE, PUBLIC:: GRPXX ( : ) ! x-location given projection
        REAL   , ALLOCATABLE, PUBLIC:: GRPYY ( : ) ! y-location given projection
        REAL   , ALLOCATABLE, PUBLIC:: GRPACRES( :) ! acres/day for a fire
        INTEGER, ALLOCATABLE, PUBLIC:: GRPLPING( : ) ! flag indicating a ping source group
        INTEGER, ALLOCATABLE, PUBLIC:: GRPLMAJOR( : ) ! flag indicating a major source group
         
!...........   Variables for grouping criteria, major source criteria,
!              and PinG source criteria
        INTEGER, PUBLIC :: NGRPVAR  = 0 ! No. of variables used to set groups
        INTEGER, PUBLIC :: NGRPCRIT = 0 ! No. of OR criteria to set groups
        INTEGER, PUBLIC :: MXGRPCHK = 0 ! Max no. of AND criteria per OR for grps

        INTEGER, PUBLIC :: NELVCRIT = 0 ! No. of OR criteria to set major srcs
        INTEGER, PUBLIC :: MXELVCHK = 0 ! Max no. of AND criteria per OR for grps
        INTEGER, PUBLIC :: NELVCSRC = 0 ! No. specific sources in input file

        INTEGER, PUBLIC :: NPNGCRIT = 0 ! No. of OR criteria to set PinG sources
        INTEGER, PUBLIC :: MXPNGCHK = 0 ! Max no. of AND criteria per OR for PinG
        INTEGER, PUBLIC :: NPNGCSRC = 0 ! No. specific sources in input file

        INTEGER, PUBLIC :: NEVPEMV  = 0 ! No. of emissions vars used for both
        INTEGER, PUBLIC :: NEVPVAR  = 0 ! No. of all vars used for both

        LOGICAL, PUBLIC :: LELVRNK  = .FALSE. ! true: use ranking for elevated criteria
        LOGICAL, PUBLIC :: LPNGRNK  = .FALSE. ! true: use ranking for PinG criteria
        LOGICAL, PUBLIC :: LCUTOFF  = .FALSE. ! true: plume rise cutoff is used

!...........   Allocatable arrays for grouping criteria, major source criteria,
!              and PinG source criteria
        REAL,         ALLOCATABLE, PUBLIC :: GRPVALS ( :,:,: ) ! comprisn value
        CHARACTER(6), ALLOCATABLE, PUBLIC :: GRPTYPES( :,:,: ) ! comprisn type

        REAL,         ALLOCATABLE, PUBLIC :: ELVVALS ( :,:,: ) ! comprisn value
        CHARACTER(6), ALLOCATABLE, PUBLIC :: ELVTYPES( :,:,: ) ! comprisn type
        CHARACTER(PLTLEN3), ALLOCATABLE, PUBLIC :: ELVCHRS( :,:,: ) ! cmpr string
        CHARACTER(ALLLEN3), ALLOCATABLE, PUBLIC :: ELVCSRC( : )

        REAL,         ALLOCATABLE, PUBLIC :: PNGVALS ( :,:,: ) ! comprisn value
        CHARACTER(6), ALLOCATABLE, PUBLIC :: PNGTYPES( :,:,: ) ! comprisn type
        CHARACTER(PLTLEN3), ALLOCATABLE, PUBLIC :: PNGCHRS( :,:,: ) ! cmpr string
        CHARACTER(ALLLEN3), ALLOCATABLE, PUBLIC :: PNGCSRC( : )

        INTEGER, ALLOCATABLE, PUBLIC :: EVPEMIDX( : )  ! indx to EINAM
        LOGICAL, ALLOCATABLE, PUBLIC :: EVPESTAT( : )  ! true: used for elv
        LOGICAL, ALLOCATABLE, PUBLIC :: EVPPSTAT( : )  ! true: used for PinG

!...........   Allocatable arrays for computing data by source
        INTEGER, ALLOCATABLE, PUBLIC :: MXEIDX( :,: )  ! sorting index
        INTEGER, ALLOCATABLE, PUBLIC :: MXRANK( :,: )  ! sorted rank
        REAL,    ALLOCATABLE, PUBLIC :: MXEMIS( :,: )  ! maximum daily emis
        REAL,    ALLOCATABLE, PUBLIC :: RISE  ( : )    ! analytical plume rise

!...........   Allocatable arrays for major and PinG sources indexing
        INTEGER, ALLOCATABLE, PUBLIC:: ELEVSIDX( : ) ! Elev source -> all srcs
        INTEGER, ALLOCATABLE, PUBLIC:: PINGGIDX( : ) ! PinG source -> PinG group

!.........  Elevated source and plume-in-grid emissions arrays
        REAL,    ALLOCATABLE, PUBLIC :: PGRPEMIS( : ) ! PinG group emissions
        REAL,    ALLOCATABLE, PUBLIC :: ELEVEMIS( : ) ! Major elev srcs emis

!.........  Hourly plume rise file values and arrays
        INTEGER,              PUBLIC :: NHRSRC       ! number of explicit srcs
        INTEGER, ALLOCATABLE, PUBLIC :: INDXH  ( : ) ! source list (by hour)
        INTEGER, ALLOCATABLE, PUBLIC :: ELEVSRC( : ) ! static source list
        REAL   , ALLOCATABLE, PUBLIC :: LAY1F  ( : ) ! layer-1 fraction
        REAL   , ALLOCATABLE, PUBLIC :: PLMBOT ( : ) ! plume bottom [m]
        REAL   , ALLOCATABLE, PUBLIC :: PLMTOP ( : ) ! plume top [m]
        REAL   , ALLOCATABLE, PUBLIC :: HRSTKTK( : ) ! stack temperature [K]
        REAL   , ALLOCATABLE, PUBLIC :: HRSTKVE( : ) ! stack exit velocity [m/s]
        REAL   , ALLOCATABLE, PUBLIC :: HRSTKFL( : ) ! stack exit flow rate [m/s]

!.........  Processing wildfire sources
        REAL   , ALLOCATABLE, PUBLIC :: DAY_ACRES(:) ! number of acres per day
        INTEGER, ALLOCATABLE, PUBLIC :: DAY_INDEX(:)
        REAL   , ALLOCATABLE, PUBLIC :: ACRES(:)
        LOGICAL, PUBLIC :: FFLAG    = .TRUE.  ! true if source sector is a fire source

!.........  Source apportionment revised groups
        INTEGER,              PUBLIC :: NELEVGRPS       ! number of revised stack groups
        INTEGER, ALLOCATABLE, PUBLIC :: ELEVGRPID( : )  ! revised stack group ID per src
        INTEGER, ALLOCATABLE, PUBLIC :: ELEVSTKGRP( : ) ! mapping to orig. stack group
        INTEGER, ALLOCATABLE, PUBLIC :: ELEVSRCGRP( : ) ! mapping to orig. src group
        INTEGER, ALLOCATABLE, PUBLIC :: ELEVSTKCNT( : ) ! stack count per group
        REAL,    ALLOCATABLE, PUBLIC :: EMELEVGRP( : )  ! summed group emissions
        LOGICAL, PUBLIC :: SGFIREFLAG = .FALSE. ! true if stack groups has fire data

        END MODULE MODELEV
