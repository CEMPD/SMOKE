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
! COPYRIGHT (C) 1999, MCNC--North Carolina Supercomputing Center
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

!.........  Per-source arrays with temporal profile indices to tables
        INTEGER, ALLOCATABLE :: MDEX ( :,: ) ! monthly       profile subscript
        INTEGER, ALLOCATABLE :: WDEX ( :,: ) ! weeky         profile subscript
        INTEGER, ALLOCATABLE :: DDEX ( :,: ) ! wkday diurnal profile subscript
        INTEGER, ALLOCATABLE :: EDEX ( :,: ) ! wkend-diurnal profile subscript
 
!.........  Per-source arrays with control packet indices to tables.  
        INTEGER, ALLOCATABLE :: CTGIDX( :,: ) ! index for CTG data tables
        INTEGER, ALLOCATABLE :: CTLIDX( :,: ) ! index for CONTROL data tables
        INTEGER, ALLOCATABLE :: ALWIDX( :,: ) ! index for ALLOWABLE data tables
        INTEGER, ALLOCATABLE :: ADDIDX( :,: ) ! index for ADD data tables

!.........  Per-source arrays for position in gridding surrogates table
        INTEGER, ALLOCATABLE :: SRGIDPOS( : ) ! surrogate position
        INTEGER, ALLOCATABLE :: SGFIPPOS( : ) ! cy/st/co code position

!.........  Per-source arrays with index for assigning emission factors (by 
!           pollutant). Index goes to IPSIA.
        INTEGER, ALLOCATABLE :: EFSIDX( :,: ) 

!.........  Number to add to monthly profile to indicate pollutant-specific

        INTEGER, PARAMETER :: ADDPS = 900000 ! multiple of 9

!.........  Count of number of entries in each (non-default) table

        INTEGER, PARAMETER, PUBLIC :: NXTYPES = 16
        INTEGER, PUBLIC            :: TXCNT( NXTYPES )

!.........  Sorted groups of cross-references
!.........  For temporal entries, integer values are the temporal profile codes
!.........  For speciation, CSPT* values are the speciation profile codes
!.........  For control entries, integer values are the stored position in table
!.........  For emission factors, 

!.........  Default FIPS code=0, SCC=0

        INTEGER               , ALLOCATABLE, PUBLIC :: MPRT01( : ) ! mon-of-yr
        INTEGER               , ALLOCATABLE, PUBLIC :: WPRT01( : ) ! day-of-week
        INTEGER               , ALLOCATABLE, PUBLIC :: DPRT01( : ) ! diurnal
        INTEGER               , ALLOCATABLE, PUBLIC :: ICTL01( : ) ! for cntrls
        INTEGER               , ALLOCATABLE, PUBLIC :: IEFS01( : ) ! for EFs
        INTEGER               ,              PUBLIC :: ISRG01      ! for surgts
        CHARACTER(LEN=SPNLEN3), ALLOCATABLE, PUBLIC :: CSPT01( : ) ! spec prof

!.........  FIPS code=0, SCC=left  (dim: n<entries> * npol)

        INTEGER               , ALLOCATABLE, PUBLIC :: MPRT02( :,: )
        INTEGER               , ALLOCATABLE, PUBLIC :: WPRT02( :,: )
        INTEGER               , ALLOCATABLE, PUBLIC :: DPRT02( :,: )
        INTEGER               , ALLOCATABLE, PUBLIC :: ICTL02( :,: )
        INTEGER               , ALLOCATABLE, PUBLIC :: IEFS02( :,: )
        INTEGER               , ALLOCATABLE, PUBLIC :: ISRG02( : )
        CHARACTER(LEN=SPNLEN3), ALLOCATABLE, PUBLIC :: CSPT02( :,: )
        CHARACTER(LEN=SCCLEN3), ALLOCATABLE, PUBLIC :: CHRT02( : )

!.........  FIPS code=0, SCC=all  

        INTEGER               , ALLOCATABLE, PUBLIC :: MPRT03( :,: )
        INTEGER               , ALLOCATABLE, PUBLIC :: WPRT03( :,: )
        INTEGER               , ALLOCATABLE, PUBLIC :: DPRT03( :,: )
        INTEGER               , ALLOCATABLE, PUBLIC :: ICTL03( :,: )
        INTEGER               , ALLOCATABLE, PUBLIC :: IEFS03( :,: )
        INTEGER               , ALLOCATABLE, PUBLIC :: ISRG03( : )
        CHARACTER(LEN=SPNLEN3), ALLOCATABLE, PUBLIC :: CSPT03( :,: )
        CHARACTER(LEN=SCCLEN3), ALLOCATABLE, PUBLIC :: CHRT03( : )

!.........  FIPS code=state default, SCC=0  

        INTEGER               , ALLOCATABLE, PUBLIC :: MPRT04( : )
        INTEGER               , ALLOCATABLE, PUBLIC :: WPRT04( : )
        INTEGER               , ALLOCATABLE, PUBLIC :: DPRT04( : )
        INTEGER               , ALLOCATABLE, PUBLIC :: ICTL04( :,: )
        INTEGER               , ALLOCATABLE, PUBLIC :: IEFS04( :,: )
        INTEGER               , ALLOCATABLE, PUBLIC :: ISRG04( : )
        CHARACTER(LEN=SPNLEN3), ALLOCATABLE, PUBLIC :: CSPT04( :,: )
        CHARACTER(LEN=STALEN3), ALLOCATABLE, PUBLIC :: CHRT04( : )

!.........  FIPS code=state default, SCC=left

        INTEGER               , ALLOCATABLE, PUBLIC :: MPRT05( :,: )
        INTEGER               , ALLOCATABLE, PUBLIC :: WPRT05( :,: )
        INTEGER               , ALLOCATABLE, PUBLIC :: DPRT05( :,: )
        INTEGER               , ALLOCATABLE, PUBLIC :: ICTL05( :,: )
        INTEGER               , ALLOCATABLE, PUBLIC :: IEFS05( :,: )
        INTEGER               , ALLOCATABLE, PUBLIC :: ISRG05( : )
        CHARACTER(LEN=SPNLEN3), ALLOCATABLE, PUBLIC :: CSPT05( :,: )
        CHARACTER(LEN=STSLEN3), ALLOCATABLE, PUBLIC :: CHRT05( : )

!.........  FIPS code=state, SCC=all
        INTEGER               , ALLOCATABLE, PUBLIC :: MPRT06( :,: )
        INTEGER               , ALLOCATABLE, PUBLIC :: WPRT06( :,: )
        INTEGER               , ALLOCATABLE, PUBLIC :: DPRT06( :,: )
        INTEGER               , ALLOCATABLE, PUBLIC :: ICTL06( :,: )
        INTEGER               , ALLOCATABLE, PUBLIC :: IEFS06( :,: )
        INTEGER               , ALLOCATABLE, PUBLIC :: ISRG06( : )
        CHARACTER(LEN=SPNLEN3), ALLOCATABLE, PUBLIC :: CSPT06( :,: )
        CHARACTER(LEN=STSLEN3), ALLOCATABLE, PUBLIC :: CHRT06( : )

!.........  FIPS code=all, SCC=0
        INTEGER               , ALLOCATABLE, PUBLIC :: MPRT07( : )
        INTEGER               , ALLOCATABLE, PUBLIC :: WPRT07( : )
        INTEGER               , ALLOCATABLE, PUBLIC :: DPRT07( : )
        INTEGER               , ALLOCATABLE, PUBLIC :: ICTL07( :,: )
        INTEGER               , ALLOCATABLE, PUBLIC :: IEFS07( :,: )
        INTEGER               , ALLOCATABLE, PUBLIC :: ISRG07( : )
        CHARACTER(LEN=SPNLEN3), ALLOCATABLE, PUBLIC :: CSPT07( :,: )
        CHARACTER(LEN=FIPLEN3), ALLOCATABLE, PUBLIC :: CHRT07( : )
        
!.........  FIPS code=all, SCC=left
        INTEGER               , ALLOCATABLE, PUBLIC :: MPRT08( :,: )
        INTEGER               , ALLOCATABLE, PUBLIC :: WPRT08( :,: )
        INTEGER               , ALLOCATABLE, PUBLIC :: DPRT08( :,: )
        INTEGER               , ALLOCATABLE, PUBLIC :: ICTL08( :,: )
        INTEGER               , ALLOCATABLE, PUBLIC :: IEFS08( :,: )
        INTEGER               , ALLOCATABLE, PUBLIC :: ISRG08( : )
        CHARACTER(LEN=SPNLEN3), ALLOCATABLE, PUBLIC :: CSPT08( :,: )
        CHARACTER(LEN=FPSLEN3), ALLOCATABLE, PUBLIC :: CHRT08( : )

!.........  FIPS code=all, SCC=all
        INTEGER               , ALLOCATABLE, PUBLIC :: MPRT09( :,: )
        INTEGER               , ALLOCATABLE, PUBLIC :: WPRT09( :,: )
        INTEGER               , ALLOCATABLE, PUBLIC :: DPRT09( :,: )
        INTEGER               , ALLOCATABLE, PUBLIC :: ICTL09( :,: )
        INTEGER               , ALLOCATABLE, PUBLIC :: IEFS09( :,: )
        INTEGER               , ALLOCATABLE, PUBLIC :: ISRG09( : )
        CHARACTER(LEN=SPNLEN3), ALLOCATABLE, PUBLIC :: CSPT09( :,: )
        CHARACTER(LEN=FPSLEN3), ALLOCATABLE, PUBLIC :: CHRT09( : )

!.........  PLANT=non-blank, SCC=0
        INTEGER               , ALLOCATABLE, PUBLIC :: MPRT10( :,: )
        INTEGER               , ALLOCATABLE, PUBLIC :: WPRT10( :,: )
        INTEGER               , ALLOCATABLE, PUBLIC :: DPRT10( :,: )
        INTEGER               , ALLOCATABLE, PUBLIC :: ICTL10( :,: )
        INTEGER               , ALLOCATABLE, PUBLIC :: IEFS10( :,: )
        CHARACTER(LEN=SPNLEN3), ALLOCATABLE, PUBLIC :: CSPT10( :,: )
        CHARACTER(LEN=FPLLEN3), ALLOCATABLE, PUBLIC :: CHRT10( : )

!.........  Plant=non-blank, SCC=all
        INTEGER               , ALLOCATABLE, PUBLIC :: MPRT11( :,: )
        INTEGER               , ALLOCATABLE, PUBLIC :: WPRT11( :,: )
        INTEGER               , ALLOCATABLE, PUBLIC :: DPRT11( :,: )
        INTEGER               , ALLOCATABLE, PUBLIC :: ICTL11( :,: )
        INTEGER               , ALLOCATABLE, PUBLIC :: IEFS11( :,: )
        CHARACTER(LEN=SPNLEN3), ALLOCATABLE, PUBLIC :: CSPT11( :,: )
        CHARACTER(LEN=SS0LEN3), ALLOCATABLE, PUBLIC :: CHRT11( : )

!.........  CHAR1=non-blank, SCC=all
        INTEGER               , ALLOCATABLE, PUBLIC :: MPRT12( :,: )
        INTEGER               , ALLOCATABLE, PUBLIC :: WPRT12( :,: )
        INTEGER               , ALLOCATABLE, PUBLIC :: DPRT12( :,: )
        INTEGER               , ALLOCATABLE, PUBLIC :: ICTL12( :,: )
        INTEGER               , ALLOCATABLE, PUBLIC :: IEFS12( :,: )
        CHARACTER(LEN=SPNLEN3), ALLOCATABLE, PUBLIC :: CSPT12( :,: )
        CHARACTER(LEN=SS1LEN3), ALLOCATABLE, PUBLIC :: CHRT12( : )

!.........  CHAR2=non-blank, SCC=all
        INTEGER               , ALLOCATABLE, PUBLIC :: MPRT13( :,: )
        INTEGER               , ALLOCATABLE, PUBLIC :: WPRT13( :,: )
        INTEGER               , ALLOCATABLE, PUBLIC :: DPRT13( :,: )
        INTEGER               , ALLOCATABLE, PUBLIC :: ICTL13( :,: )
        INTEGER               , ALLOCATABLE, PUBLIC :: IEFS13( :,: )
        CHARACTER(LEN=SPNLEN3), ALLOCATABLE, PUBLIC :: CSPT13( :,: )
        CHARACTER(LEN=SS2LEN3), ALLOCATABLE, PUBLIC :: CHRT13( : )

!.........  CHAR3=non-blank, SCC=all
        INTEGER               , ALLOCATABLE, PUBLIC :: MPRT14( :,: )
        INTEGER               , ALLOCATABLE, PUBLIC :: WPRT14( :,: )
        INTEGER               , ALLOCATABLE, PUBLIC :: DPRT14( :,: )
        INTEGER               , ALLOCATABLE, PUBLIC :: ICTL14( :,: )
        INTEGER               , ALLOCATABLE, PUBLIC :: IEFS14( :,: )
        CHARACTER(LEN=SPNLEN3), ALLOCATABLE, PUBLIC :: CSPT14( :,: )
        CHARACTER(LEN=SS3LEN3), ALLOCATABLE, PUBLIC :: CHRT14( : )

!.........  CHAR4=non-blank, SCC=all
        INTEGER               , ALLOCATABLE, PUBLIC :: MPRT15( :,: )
        INTEGER               , ALLOCATABLE, PUBLIC :: WPRT15( :,: )
        INTEGER               , ALLOCATABLE, PUBLIC :: DPRT15( :,: )
        INTEGER               , ALLOCATABLE, PUBLIC :: ICTL15( :,: )
        INTEGER               , ALLOCATABLE, PUBLIC :: IEFS15( :,: )
        CHARACTER(LEN=SPNLEN3), ALLOCATABLE, PUBLIC :: CSPT15( :,: )
        CHARACTER(LEN=SS4LEN3), ALLOCATABLE, PUBLIC :: CHRT15( : )

!.........  CHAR5=non-blank, SCC=all
        INTEGER               , ALLOCATABLE, PUBLIC :: MPRT16( :,: )
        INTEGER               , ALLOCATABLE, PUBLIC :: WPRT16( :,: )
        INTEGER               , ALLOCATABLE, PUBLIC :: DPRT16( :,: )
        INTEGER               , ALLOCATABLE, PUBLIC :: ICTL16( :,: )
        INTEGER               , ALLOCATABLE, PUBLIC :: IEFS16( :,: )
        CHARACTER(LEN=SPNLEN3), ALLOCATABLE, PUBLIC :: CSPT16( :,: )
        CHARACTER(LEN=SS5LEN3), ALLOCATABLE, PUBLIC :: CHRT16( : )

!.........  Unsorted, unprocessed cross-reference arrays
        INTEGER, ALLOCATABLE, PUBLIC:: INDXTA ( : ) !  sorting index
        INTEGER, ALLOCATABLE, PUBLIC:: IFIPTA ( : ) !  co/st/cy FIPS codes
        INTEGER, ALLOCATABLE, PUBLIC:: ISPTA  ( : ) !  pollutant index to EINAM
        INTEGER, ALLOCATABLE, PUBLIC:: ISRGCDA( : ) !  spatial surrogate codes
        INTEGER, ALLOCATABLE, PUBLIC:: MPRNA  ( : ) !  monthly profile codes
        INTEGER, ALLOCATABLE, PUBLIC:: WPRNA  ( : ) !  weekly profile codes
        INTEGER, ALLOCATABLE, PUBLIC:: DPRNA  ( : ) !  diurnal profile codes
        INTEGER, ALLOCATABLE, PUBLIC:: IPSIA( :,: ) !  24 hours code for EFs

        CHARACTER(LEN=SPNLEN3), ALLOCATABLE, PUBLIC:: CSPRNA( : ) ! spec prof #
        CHARACTER(LEN=SCCLEN3), ALLOCATABLE, PUBLIC:: CSCCTA( : ) ! SCC

        CHARACTER(LEN=SS5LEN3+POLLEN3), ALLOCATABLE, PUBLIC:: CSRCTA(:)
                                             ! source chars // SCC // pollutant

        END MODULE MODXREF
