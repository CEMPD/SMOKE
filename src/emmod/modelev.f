
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
        REAL       , ALLOCATABLE, PUBLIC :: GRPVALS ( :,:,: ) ! comprisn value
        CHARACTER*6, ALLOCATABLE, PUBLIC :: GRPTYPES( :,:,: ) ! comprisn type

        REAL       , ALLOCATABLE, PUBLIC :: ELVVALS ( :,:,: ) ! comprisn value
        CHARACTER*6, ALLOCATABLE, PUBLIC :: ELVTYPES( :,:,: ) ! comprisn type
        CHARACTER(LEN=PLTLEN3), ALLOCATABLE, PUBLIC :: ELVCHRS( :,:,: ) ! cmpr string
        CHARACTER(LEN=ALLLEN3), ALLOCATABLE, PUBLIC :: ELVCSRC( : )

        REAL       , ALLOCATABLE, PUBLIC :: PNGVALS ( :,:,: ) ! comprisn value
        CHARACTER*6, ALLOCATABLE, PUBLIC :: PNGTYPES( :,:,: ) ! comprisn type
        CHARACTER(LEN=PLTLEN3), ALLOCATABLE, PUBLIC :: PNGCHRS( :,:,: ) ! cmpr string
        CHARACTER(LEN=ALLLEN3), ALLOCATABLE, PUBLIC :: PNGCSRC( : )

        INTEGER    , ALLOCATABLE, PUBLIC :: EVPEMIDX( : )  ! indx to EINAM
        LOGICAL    , ALLOCATABLE, PUBLIC :: EVPESTAT( : )  ! true: used for elv
        LOGICAL    , ALLOCATABLE, PUBLIC :: EVPPSTAT( : )  ! true: used for PinG

!...........   Allocatable arrays for computing data by source
        INTEGER    , ALLOCATABLE, PUBLIC :: MXEIDX( :,: )  ! sorting index
        INTEGER    , ALLOCATABLE, PUBLIC :: MXRANK( :,: )  ! sorted rank
        REAL       , ALLOCATABLE, PUBLIC :: MXEMIS( :,: )  ! maximum daily emis
        REAL       , ALLOCATABLE, PUBLIC :: RISE  ( : )    ! analytical plume rise

!...........   Allocatable arrays for major and PinG sources indexing
        INTEGER, ALLOCATABLE, PUBLIC:: ELEVSIDX( : ) ! Elev source -> all srcs
        INTEGER, ALLOCATABLE, PUBLIC:: PINGGIDX( : ) ! PinG source -> PinG group

!.........  Elevated source and plume-in-grid emissions arrays
        REAL   , ALLOCATABLE, PUBLIC :: PGRPEMIS( : ) ! PinG group emissions
        REAL   , ALLOCATABLE, PUBLIC :: ELEVEMIS( : ) ! Major elev srcs emis

!.........  Hourly plume rise file values and arrays
        INTEGER             , PUBLIC :: NHRSRC       ! number of explicit srcs
        INTEGER, ALLOCATABLE, PUBLIC :: INDXH  ( : ) ! source list (by hour)
        INTEGER, ALLOCATABLE, PUBLIC :: ELEVSRC( : ) ! static source list
        REAL   , ALLOCATABLE, PUBLIC :: LAY1F  ( : ) ! layer-1 fraction
        REAL   , ALLOCATABLE, PUBLIC :: PLMBOT ( : ) ! plume bottom [m]
        REAL   , ALLOCATABLE, PUBLIC :: PLMTOP ( : ) ! plume top [m]
        REAL   , ALLOCATABLE, PUBLIC :: HRSTKTK( : ) ! stack temperature [K]
        REAL   , ALLOCATABLE, PUBLIC :: HRSTKVE( : ) ! stack exit velocity [m/s]
        REAL   , ALLOCATABLE, PUBLIC :: HRSTKFL( : ) ! stack exit flow rate [m/s]

        END MODULE MODELEV
