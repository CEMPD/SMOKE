        MODULE MODXREF

!***********************************************************************
!  Module body starts at line 40
!
!  DESCRIPTION:
!     This module contains the public allocatable arrays for cross-reference
!     tables
!
!  PRECONDITIONS REQUIRED:
!
!  SUBROUTINES AND FUNCTIONS CALLED:
!
!  REVISION HISTORY:
!     Created 1/99 by M. Houyoux
!
!***************************************************************************
!
! Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
!                System
! File: @(#)$Id$
!
! COPYRIGHT (C) 1998, MCNC--North Carolina Supercomputing Center
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

        INCLUDE 'EMPRVT3.EXT'   !  private emissions string widths parameters

!.........  Common per-source arrays with profile numbers
        INTEGER, ALLOCATABLE :: MDEX ( :,: ) ! monthly       profile subscript
        INTEGER, ALLOCATABLE :: WDEX ( :,: ) ! weeky         profile subscript
        INTEGER, ALLOCATABLE :: DDEX ( :,: ) ! wkday diurnal profile subscript
        INTEGER, ALLOCATABLE :: EDEX ( :,: ) ! wkend-diurnal profile subscript
 
!.........  Number to add to monthly profile to indicate pollutant-specific

        INTEGER, PARAMETER :: ADDPS = 90000

!.........  Count of number of entries in each (non-default) table

        INTEGER, PARAMETER, PUBLIC :: NXTYPES = 16
        INTEGER, PUBLIC            :: TXCNT( NXTYPES )

!.........  Sorted groups of point source temporal cross-references 
!.........  Default FIPS code=0, SCC=0

        INTEGER         MPRT01  !  month-of-year
        INTEGER         WPRT01  !  day-of-week
        INTEGER         DPRT01  !  diurnal

!.........  FIPS code=0, SCC=left  ( TXCNT(2) )

        INTEGER               , ALLOCATABLE, PUBLIC :: MPRT02( :,: )
        INTEGER               , ALLOCATABLE, PUBLIC :: WPRT02( :,: )
        INTEGER               , ALLOCATABLE, PUBLIC :: DPRT02( :,: )
        CHARACTER(LEN=SCCLEN3), ALLOCATABLE, PUBLIC :: CHRT02( : )

!.........  FIPS code=0, SCC=all  

        INTEGER               , ALLOCATABLE, PUBLIC :: MPRT03( :,: )
        INTEGER               , ALLOCATABLE, PUBLIC :: WPRT03( :,: )
        INTEGER               , ALLOCATABLE, PUBLIC :: DPRT03( :,: )
        CHARACTER(LEN=SCCLEN3), ALLOCATABLE, PUBLIC :: CHRT03( : )

!.........  FIPS code=state default, SCC=0  

        INTEGER               , ALLOCATABLE, PUBLIC :: MPRT04( : )
        INTEGER               , ALLOCATABLE, PUBLIC :: WPRT04( : )
        INTEGER               , ALLOCATABLE, PUBLIC :: DPRT04( : )
        CHARACTER(LEN=STALEN3), ALLOCATABLE, PUBLIC :: CHRT04( : )

!.........  FIPS code=state default, SCC=left

        INTEGER               , ALLOCATABLE, PUBLIC :: MPRT05( :,: )
        INTEGER               , ALLOCATABLE, PUBLIC :: WPRT05( :,: )
        INTEGER               , ALLOCATABLE, PUBLIC :: DPRT05( :,: )
        CHARACTER(LEN=STSLEN3), ALLOCATABLE, PUBLIC :: CHRT05( : )

!.........  FIPS code=state, SCC=all
        INTEGER               , ALLOCATABLE, PUBLIC :: MPRT06( :,: )
        INTEGER               , ALLOCATABLE, PUBLIC :: WPRT06( :,: )
        INTEGER               , ALLOCATABLE, PUBLIC :: DPRT06( :,: )
        CHARACTER(LEN=STSLEN3), ALLOCATABLE, PUBLIC :: CHRT06( : )

!.........  FIPS code=all, SCC=0
        INTEGER               , ALLOCATABLE, PUBLIC :: MPRT07( : )
        INTEGER               , ALLOCATABLE, PUBLIC :: WPRT07( : )
        INTEGER               , ALLOCATABLE, PUBLIC :: DPRT07( : )
        CHARACTER(LEN=FIPLEN3), ALLOCATABLE, PUBLIC :: CHRT07( : )
        
!.........  FIPS code=all, SCC=left
        INTEGER               , ALLOCATABLE, PUBLIC :: MPRT08( :,: )
        INTEGER               , ALLOCATABLE, PUBLIC :: WPRT08( :,: )
        INTEGER               , ALLOCATABLE, PUBLIC :: DPRT08( :,: )
        CHARACTER(LEN=FPSLEN3), ALLOCATABLE, PUBLIC :: CHRT08( : )

!.........  FIPS code=all, SCC=all
        INTEGER               , ALLOCATABLE, PUBLIC :: MPRT09( :,: )
        INTEGER               , ALLOCATABLE, PUBLIC :: WPRT09( :,: )
        INTEGER               , ALLOCATABLE, PUBLIC :: DPRT09( :,: )
        CHARACTER(LEN=FPSLEN3), ALLOCATABLE, PUBLIC :: CHRT09( : )

!.........  PLANT=non-blank, SCC=0
        INTEGER               , ALLOCATABLE, PUBLIC :: MPRT10( :,: )
        INTEGER               , ALLOCATABLE, PUBLIC :: WPRT10( :,: )
        INTEGER               , ALLOCATABLE, PUBLIC :: DPRT10( :,: )
        CHARACTER(LEN=FPLLEN3), ALLOCATABLE, PUBLIC :: CHRT10( : )

!.........  Plant=non-blank, SCC=all
        INTEGER               , ALLOCATABLE, PUBLIC :: MPRT11( :,: )
        INTEGER               , ALLOCATABLE, PUBLIC :: WPRT11( :,: )
        INTEGER               , ALLOCATABLE, PUBLIC :: DPRT11( :,: )
        CHARACTER(LEN=SS0LEN3), ALLOCATABLE, PUBLIC :: CHRT11( : )

!.........  CHAR1=non-blank, SCC=all
        INTEGER               , ALLOCATABLE, PUBLIC :: MPRT12( :,: )
        INTEGER               , ALLOCATABLE, PUBLIC :: WPRT12( :,: )
        INTEGER               , ALLOCATABLE, PUBLIC :: DPRT12( :,: )
        CHARACTER(LEN=SS1LEN3), ALLOCATABLE, PUBLIC :: CHRT12( : )

!.........  CHAR2=non-blank, SCC=all
        INTEGER               , ALLOCATABLE, PUBLIC :: MPRT13( :,: )
        INTEGER               , ALLOCATABLE, PUBLIC :: WPRT13( :,: )
        INTEGER               , ALLOCATABLE, PUBLIC :: DPRT13( :,: )
        CHARACTER(LEN=SS2LEN3), ALLOCATABLE, PUBLIC :: CHRT13( : )

!.........  CHAR3=non-blank, SCC=all
        INTEGER               , ALLOCATABLE, PUBLIC :: MPRT14( :,: )
        INTEGER               , ALLOCATABLE, PUBLIC :: WPRT14( :,: )
        INTEGER               , ALLOCATABLE, PUBLIC :: DPRT14( :,: )
        CHARACTER(LEN=SS3LEN3), ALLOCATABLE, PUBLIC :: CHRT14( : )

!.........  CHAR4=non-blank, SCC=all
        INTEGER               , ALLOCATABLE, PUBLIC :: MPRT15( :,: )
        INTEGER               , ALLOCATABLE, PUBLIC :: WPRT15( :,: )
        INTEGER               , ALLOCATABLE, PUBLIC :: DPRT15( :,: )
        CHARACTER(LEN=SS4LEN3), ALLOCATABLE, PUBLIC :: CHRT15( : )

!.........  CHAR5=non-blank, SCC=all
        INTEGER               , ALLOCATABLE, PUBLIC :: MPRT16( :,: )
        INTEGER               , ALLOCATABLE, PUBLIC :: WPRT16( :,: )
        INTEGER               , ALLOCATABLE, PUBLIC :: DPRT16( :,: )
        CHARACTER(LEN=SS5LEN3), ALLOCATABLE, PUBLIC :: CHRT16( : )

!.........  Unsorted, unprocessed cross-reference arrays
        INTEGER, ALLOCATABLE, PUBLIC:: INDXTA( : )  !  sorting index
        INTEGER, ALLOCATABLE, PUBLIC:: IFIPTA( : )  !  co/st/cy FIPS codes
        INTEGER, ALLOCATABLE, PUBLIC:: ISPTA ( : )  !  pollutant codes  
        INTEGER, ALLOCATABLE, PUBLIC:: MPRNA ( : )  !  monthly profile codes
        INTEGER, ALLOCATABLE, PUBLIC:: WPRNA ( : )  !  weekly profile codes
        INTEGER, ALLOCATABLE, PUBLIC:: DPRNA ( : )  !  diurnal profile codes

        CHARACTER(LEN=SPNLEN3), ALLOCATABLE, PUBLIC:: CSPRNA( : ) ! spec prof #
        CHARACTER(LEN=SCCLEN3), ALLOCATABLE, PUBLIC:: CSCCTA( : ) ! SCC

        CHARACTER(LEN=SS5LEN3+POLLEN3), ALLOCATABLE, PUBLIC:: CSRCTA(:)
                                             ! source chars // SCC // pollutant

        END MODULE MODXREF
