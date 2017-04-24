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

!.........  Turn on/off xref duplicate-checking in "lib/xreftbl.f":
        LOGICAL, SAVE :: XDUPCHK = .TRUE.

!.........  Per-source arrays with temporal profile indices to tables
        INTEGER, ALLOCATABLE :: MDEX ( :,: ) ! monthly       profile subscript
        INTEGER, ALLOCATABLE :: WDEX ( :,: ) ! weeky         profile subscript
        INTEGER, ALLOCATABLE :: DDEX ( :,: ) ! wkday diurnal profile subscript
        INTEGER, ALLOCATABLE :: EDEX ( :,: ) ! wkend-diurnal profile subscript
 
!.........  Per-source arrays with control packet indices to tables.  
        INTEGER, ALLOCATABLE :: ASGNINDX( : ) ! index for projection/controls

!.........  Per-source arrays for position in gridding surrogates table
        INTEGER, ALLOCATABLE, PUBLIC :: ASRGID   ( : ) ! tmp surrgate codes
        INTEGER, ALLOCATABLE, PUBLIC :: SRGFIPIDX( : ) ! tmp surrgate codes

!.........  Per-source arrays with index for assigning emission factors (by 
!           pollutant). Index goes to IPSIA.
        INTEGER, ALLOCATABLE :: EFSIDX( :,: ) 

!.........  Per-source arrays with index and number of entries for
!           assigning multiple point locations to area sources
        INTEGER, ALLOCATABLE :: AR2PTTBL( : )  ! index to table number
        INTEGER, ALLOCATABLE :: AR2PTIDX( : )  ! index to position in table
        INTEGER, ALLOCATABLE :: AR2PTCNT( : )  ! count of entries in table

!.........  Per-source array to indicate whether source receives 
!           NONHAPVOC or NONHAPTOG calculation
        LOGICAL, ALLOCATABLE :: LNONHAP( : )   ! true: NONHAP is computed
        LOGICAL              :: NHAP_EXCL      ! true: Exclude NONHAP computation
        CHARACTER(10)        :: PROC_HAPS      ! how-to process HAPs[all,none,partial] 

!.........  Per-source array with speed profile ID (-1 indicates use inventory speed)
        INTEGER, ALLOCATABLE :: SPDPROFID( : )  ! speed profile code number

!.........  Number to add to monthly profile to indicate pollutant-specific

        INTEGER, PARAMETER :: ADDPS = 900000 ! multiple of 9

!.........  Count of number of entries in each (non-default) table
!           NOTE- Added 9 (from 16-25) for special SCC-level matching
!           NOTE- Added 6 (from 26-31) for special SIC matching
!           NOTE- Added 6 (from 32-38) for special MACT matching
        INTEGER, PARAMETER, PUBLIC :: NXTYPES = 38
        INTEGER, PUBLIC            :: TXCNT( NXTYPES )

!.........  Sorted groups of cross-references
!.........  For temporal entries, integer values are the temporal profile codes
!.........  For speciation, CSPT* values are the speciation profile codes
!.........  For control entries, integer values are the stored position in table
!.........  For emission factors, 

!.........  Default FIPS code=0, SCC=0
        INTEGER           , ALLOCATABLE, PUBLIC :: MPRT01( : ) ! mon-of-yr
        INTEGER           , ALLOCATABLE, PUBLIC :: WPRT01( : ) ! day-of-week
        INTEGER           , ALLOCATABLE, PUBLIC :: DPRT01( : ) ! diurnal
        INTEGER           , ALLOCATABLE, PUBLIC :: ICTL01( : ) ! for cntrls
        INTEGER           , ALLOCATABLE, PUBLIC :: IEFS01( : ) ! for EFs
        INTEGER           ,              PUBLIC :: ISRG01      ! for surgts
        INTEGER           ,              PUBLIC :: IMVS01      ! for VMT mix
        INTEGER           ,              PUBLIC :: ISPD01      ! for speeds
        CHARACTER(SPNLEN3), ALLOCATABLE, PUBLIC :: CSPT01( : ) ! spec prof

!.........  FIPS code=0, SCC=left  (dim: n<entries> * npol)

        INTEGER           , ALLOCATABLE, PUBLIC :: MPRT02( :,: )
        INTEGER           , ALLOCATABLE, PUBLIC :: WPRT02( :,: )
        INTEGER           , ALLOCATABLE, PUBLIC :: DPRT02( :,: )
        INTEGER           , ALLOCATABLE, PUBLIC :: ICTL02( :,: )
        INTEGER           , ALLOCATABLE, PUBLIC :: IEFS02( :,: )
        INTEGER           , ALLOCATABLE, PUBLIC :: ISRG02( : )
        INTEGER           , ALLOCATABLE, PUBLIC :: IMVS02( : )
        INTEGER           , ALLOCATABLE, PUBLIC :: ISPD02( : )
        CHARACTER(SPNLEN3), ALLOCATABLE, PUBLIC :: CSPT02( :,: )
        CHARACTER(SCCLEN3), ALLOCATABLE, PUBLIC :: CHRT02( : )

!.........  FIPS code=0, SCC=all (level 4)

        INTEGER           , ALLOCATABLE, PUBLIC :: MPRT03( :,: )
        INTEGER           , ALLOCATABLE, PUBLIC :: WPRT03( :,: )
        INTEGER           , ALLOCATABLE, PUBLIC :: DPRT03( :,: )
        INTEGER           , ALLOCATABLE, PUBLIC :: ICTL03( :,: )
        INTEGER           , ALLOCATABLE, PUBLIC :: IEFS03( :,: )
        INTEGER           , ALLOCATABLE, PUBLIC :: ISRG03( : )
        INTEGER           , ALLOCATABLE, PUBLIC :: IMVS03( : )
        INTEGER           , ALLOCATABLE, PUBLIC :: ISPD03( : )
        CHARACTER(SPNLEN3), ALLOCATABLE, PUBLIC :: CSPT03( :,: )
        CHARACTER(SCCLEN3), ALLOCATABLE, PUBLIC :: CHRT03( : )

!.........  FIPS code=state default, SCC=0  

        INTEGER           , ALLOCATABLE, PUBLIC :: MPRT04( :,: )
        INTEGER           , ALLOCATABLE, PUBLIC :: WPRT04( :,: )
        INTEGER           , ALLOCATABLE, PUBLIC :: DPRT04( :,: )
        INTEGER           , ALLOCATABLE, PUBLIC :: ICTL04( :,: )
        INTEGER           , ALLOCATABLE, PUBLIC :: IEFS04( :,: )
        INTEGER           , ALLOCATABLE, PUBLIC :: ISRG04( : )
        INTEGER           , ALLOCATABLE, PUBLIC :: IMVS04( : )
        INTEGER           , ALLOCATABLE, PUBLIC :: ISPD04( : )
        CHARACTER(SPNLEN3), ALLOCATABLE, PUBLIC :: CSPT04( :,: )
        CHARACTER(STALEN3), ALLOCATABLE, PUBLIC :: CHRT04( : )

!.........  FIPS code=state default, SCC=left

        INTEGER           , ALLOCATABLE, PUBLIC :: MPRT05( :,: )
        INTEGER           , ALLOCATABLE, PUBLIC :: WPRT05( :,: )
        INTEGER           , ALLOCATABLE, PUBLIC :: DPRT05( :,: )
        INTEGER           , ALLOCATABLE, PUBLIC :: ICTL05( :,: )
        INTEGER           , ALLOCATABLE, PUBLIC :: IEFS05( :,: )
        INTEGER           , ALLOCATABLE, PUBLIC :: ISRG05( : )
        INTEGER           , ALLOCATABLE, PUBLIC :: IMVS05( : )
        INTEGER           , ALLOCATABLE, PUBLIC :: ISPD05( : )
        CHARACTER(SPNLEN3), ALLOCATABLE, PUBLIC :: CSPT05( :,: )
        CHARACTER(STSLEN3), ALLOCATABLE, PUBLIC :: CHRT05( : )

!.........  FIPS code=state, SCC=all (level 4)
        INTEGER           , ALLOCATABLE, PUBLIC :: MPRT06( :,: )
        INTEGER           , ALLOCATABLE, PUBLIC :: WPRT06( :,: )
        INTEGER           , ALLOCATABLE, PUBLIC :: DPRT06( :,: )
        INTEGER           , ALLOCATABLE, PUBLIC :: ICTL06( :,: )
        INTEGER           , ALLOCATABLE, PUBLIC :: IEFS06( :,: )
        INTEGER           , ALLOCATABLE, PUBLIC :: ISRG06( : )
        INTEGER           , ALLOCATABLE, PUBLIC :: IMVS06( : )
        INTEGER           , ALLOCATABLE, PUBLIC :: ISPD06( : )
        CHARACTER(SPNLEN3), ALLOCATABLE, PUBLIC :: CSPT06( :,: )
        CHARACTER(STSLEN3), ALLOCATABLE, PUBLIC :: CHRT06( : )

!.........  FIPS code=all, SCC=0
        INTEGER           , ALLOCATABLE, PUBLIC :: MPRT07( :,: )
        INTEGER           , ALLOCATABLE, PUBLIC :: WPRT07( :,: )
        INTEGER           , ALLOCATABLE, PUBLIC :: DPRT07( :,: )
        INTEGER           , ALLOCATABLE, PUBLIC :: ICTL07( :,: )
        INTEGER           , ALLOCATABLE, PUBLIC :: IEFS07( :,: )
        INTEGER           , ALLOCATABLE, PUBLIC :: ISRG07( : )
        INTEGER           , ALLOCATABLE, PUBLIC :: IMVS07( : )
        INTEGER           , ALLOCATABLE, PUBLIC :: ISPD07( : )
        CHARACTER(SPNLEN3), ALLOCATABLE, PUBLIC :: CSPT07( :,: )
        CHARACTER(FIPLEN3), ALLOCATABLE, PUBLIC :: CHRT07( : )

!.........  FIPS code=all, SCC=left
        INTEGER           , ALLOCATABLE, PUBLIC :: MPRT08( :,: )
        INTEGER           , ALLOCATABLE, PUBLIC :: WPRT08( :,: )
        INTEGER           , ALLOCATABLE, PUBLIC :: DPRT08( :,: )
        INTEGER           , ALLOCATABLE, PUBLIC :: ICTL08( :,: )
        INTEGER           , ALLOCATABLE, PUBLIC :: IEFS08( :,: )
        INTEGER           , ALLOCATABLE, PUBLIC :: ISRG08( : )
        INTEGER           , ALLOCATABLE, PUBLIC :: IMVS08( : )
        INTEGER           , ALLOCATABLE, PUBLIC :: ISPD08( : )
        INTEGER           , ALLOCATABLE, PUBLIC :: ARPT08( :,: ) ! tbl num, row, & count
        CHARACTER(SPNLEN3), ALLOCATABLE, PUBLIC :: CSPT08( :,: )
        CHARACTER(FPSLEN3), ALLOCATABLE, PUBLIC :: CHRT08( : )

!.........  FIPS code=all, SCC=all (level 4)
        INTEGER           , ALLOCATABLE, PUBLIC :: MPRT09( :,: )
        INTEGER           , ALLOCATABLE, PUBLIC :: WPRT09( :,: )
        INTEGER           , ALLOCATABLE, PUBLIC :: DPRT09( :,: )
        INTEGER           , ALLOCATABLE, PUBLIC :: ICTL09( :,: )
        INTEGER           , ALLOCATABLE, PUBLIC :: IEFS09( :,: )
        INTEGER           , ALLOCATABLE, PUBLIC :: ISRG09( : )
        INTEGER           , ALLOCATABLE, PUBLIC :: IMVS09( : )
        INTEGER           , ALLOCATABLE, PUBLIC :: ISPD09( : )
        INTEGER           , ALLOCATABLE, PUBLIC :: ARPT09( :,: ) ! tbl num, row, & count
        CHARACTER(SPNLEN3), ALLOCATABLE, PUBLIC :: CSPT09( :,: )
        CHARACTER(FPSLEN3), ALLOCATABLE, PUBLIC :: CHRT09( : )

!.........  PLANT=non-blank, SCC=0
        INTEGER           , ALLOCATABLE, PUBLIC :: MPRT10( :,: )
        INTEGER           , ALLOCATABLE, PUBLIC :: WPRT10( :,: )
        INTEGER           , ALLOCATABLE, PUBLIC :: DPRT10( :,: )
        INTEGER           , ALLOCATABLE, PUBLIC :: ICTL10( :,: )
        INTEGER           , ALLOCATABLE, PUBLIC :: IEFS10( :,: )
        INTEGER           , ALLOCATABLE, PUBLIC :: IMVS10( : )
        INTEGER           , ALLOCATABLE, PUBLIC :: ISPD10( : )
        CHARACTER(SPNLEN3), ALLOCATABLE, PUBLIC :: CSPT10( :,: )
        CHARACTER(FPLLEN3), ALLOCATABLE, PUBLIC :: CHRT10( : )

!.........  Plant=non-blank, SCC=all
        INTEGER           , ALLOCATABLE, PUBLIC :: MPRT11( :,: )
        INTEGER           , ALLOCATABLE, PUBLIC :: WPRT11( :,: )
        INTEGER           , ALLOCATABLE, PUBLIC :: DPRT11( :,: )
        INTEGER           , ALLOCATABLE, PUBLIC :: ICTL11( :,: )
        INTEGER           , ALLOCATABLE, PUBLIC :: IEFS11( :,: )
        INTEGER           , ALLOCATABLE, PUBLIC :: IMVS11( : )
        INTEGER           , ALLOCATABLE, PUBLIC :: ISPD11( : )
        CHARACTER(SPNLEN3), ALLOCATABLE, PUBLIC :: CSPT11( :,: )
        CHARACTER(SS0LEN3), ALLOCATABLE, PUBLIC :: CHRT11( : )

!.........  CHAR1=non-blank, SCC=all or zero
        INTEGER           , ALLOCATABLE, PUBLIC :: MPRT12( :,: )
        INTEGER           , ALLOCATABLE, PUBLIC :: WPRT12( :,: )
        INTEGER           , ALLOCATABLE, PUBLIC :: DPRT12( :,: )
        INTEGER           , ALLOCATABLE, PUBLIC :: ICTL12( :,: )
        INTEGER           , ALLOCATABLE, PUBLIC :: IEFS12( :,: )
        INTEGER           , ALLOCATABLE, PUBLIC :: IMVS12( : )
        INTEGER           , ALLOCATABLE, PUBLIC :: ISPD12( : )
        CHARACTER(SPNLEN3), ALLOCATABLE, PUBLIC :: CSPT12( :,: )
        CHARACTER(SS1LEN3), ALLOCATABLE, PUBLIC :: CHRT12( : )

!.........  CHAR2=non-blank, SCC=all or zero
        INTEGER           , ALLOCATABLE, PUBLIC :: MPRT13( :,: )
        INTEGER           , ALLOCATABLE, PUBLIC :: WPRT13( :,: )
        INTEGER           , ALLOCATABLE, PUBLIC :: DPRT13( :,: )
        INTEGER           , ALLOCATABLE, PUBLIC :: ICTL13( :,: )
        INTEGER           , ALLOCATABLE, PUBLIC :: IEFS13( :,: )
        CHARACTER(SPNLEN3), ALLOCATABLE, PUBLIC :: CSPT13( :,: )
        CHARACTER(SS2LEN3), ALLOCATABLE, PUBLIC :: CHRT13( : )

!.........  CHAR3=non-blank, SCC=all or zero
        INTEGER           , ALLOCATABLE, PUBLIC :: MPRT14( :,: )
        INTEGER           , ALLOCATABLE, PUBLIC :: WPRT14( :,: )
        INTEGER           , ALLOCATABLE, PUBLIC :: DPRT14( :,: )
        INTEGER           , ALLOCATABLE, PUBLIC :: ICTL14( :,: )
        INTEGER           , ALLOCATABLE, PUBLIC :: IEFS14( :,: )
        CHARACTER(SPNLEN3), ALLOCATABLE, PUBLIC :: CSPT14( :,: )
        CHARACTER(SS3LEN3), ALLOCATABLE, PUBLIC :: CHRT14( : )

!.........  CHAR4=non-blank, SCC=all or zero
        INTEGER           , ALLOCATABLE, PUBLIC :: MPRT15( :,: )
        INTEGER           , ALLOCATABLE, PUBLIC :: WPRT15( :,: )
        INTEGER           , ALLOCATABLE, PUBLIC :: DPRT15( :,: )
        INTEGER           , ALLOCATABLE, PUBLIC :: ICTL15( :,: )
        INTEGER           , ALLOCATABLE, PUBLIC :: IEFS15( :,: )
        CHARACTER(SPNLEN3), ALLOCATABLE, PUBLIC :: CSPT15( :,: )
        CHARACTER(SS4LEN3), ALLOCATABLE, PUBLIC :: CHRT15( : )

!.........  CHAR5=non-blank, SCC=all or zero
        INTEGER           , ALLOCATABLE, PUBLIC :: MPRT16( :,: )
        INTEGER           , ALLOCATABLE, PUBLIC :: WPRT16( :,: )
        INTEGER           , ALLOCATABLE, PUBLIC :: DPRT16( :,: )
        INTEGER           , ALLOCATABLE, PUBLIC :: ICTL16( :,: )
        INTEGER           , ALLOCATABLE, PUBLIC :: IEFS16( :,: )
        CHARACTER(SPNLEN3), ALLOCATABLE, PUBLIC :: CSPT16( :,: )
        CHARACTER(SS5LEN3), ALLOCATABLE, PUBLIC :: CHRT16( : )

!.........  Additional groups for special SCC handling. Note that these are
!              put at the end since they were added later and added for
!              Cntlmat only at first.  Did not want to mess up the numbering
!              scheme for TXCNT
!.........  FIPS code = 0, SCC = level 1 (Type 17)
        INTEGER           , ALLOCATABLE, PUBLIC :: ICTL02A( :,: )
        CHARACTER(SCCLEN3), ALLOCATABLE, PUBLIC :: CHRT02A( : )

!.........  FIPS code = 0, SCC = level 2 (Type 18)
        INTEGER           , ALLOCATABLE, PUBLIC :: ICTL02B( :,: )
        CHARACTER(SCCLEN3), ALLOCATABLE, PUBLIC :: CHRT02B( : )

!.........  FIPS code = 0, SCC = level 3 (Type 19)
        INTEGER           , ALLOCATABLE, PUBLIC :: ICTL02C( :,: )
        CHARACTER(SCCLEN3), ALLOCATABLE, PUBLIC :: CHRT02C( : )

!.........  FIPS code=state default, SCC=level 1 (Type 20)
        INTEGER           , ALLOCATABLE, PUBLIC :: ICTL05A( :,: )
        CHARACTER(STSLEN3), ALLOCATABLE, PUBLIC :: CHRT05A( : )

!.........  FIPS code=state default, SCC=level 2 (Type 21)
        INTEGER           , ALLOCATABLE, PUBLIC :: ICTL05B( :,: )
        CHARACTER(STSLEN3), ALLOCATABLE, PUBLIC :: CHRT05B( : )

!.........  FIPS code=state default, SCC=level 3 (Type 22)
        INTEGER           , ALLOCATABLE, PUBLIC :: ICTL05C( :,: )
        CHARACTER(STSLEN3), ALLOCATABLE, PUBLIC :: CHRT05C( : )

!.........  FIPS code=all, SCC=level 1 (Type 23)
        INTEGER           , ALLOCATABLE, PUBLIC :: ICTL08A( :,: )
        CHARACTER(FPSLEN3), ALLOCATABLE, PUBLIC :: CHRT08A( : )

!.........  FIPS code=all, SCC=level 2 (Type 24)
        INTEGER           , ALLOCATABLE, PUBLIC :: ICTL08B( :,: )
        CHARACTER(FPSLEN3), ALLOCATABLE, PUBLIC :: CHRT08B( : )

!.........  FIPS code=all, SCC=level 3 (Type 25)
        INTEGER           , ALLOCATABLE, PUBLIC :: ICTL08C( :,: )
        CHARACTER(FPSLEN3), ALLOCATABLE, PUBLIC :: CHRT08C( : )

!.........  Additional groups for special SIC handling. Note that these are
!              put at the end since they were added later and added for
!              Cntlmat only at first.  Did not want to mess up the numbering
!              scheme for TXCNT
!.........  FIPS code = 0, SIC = 2-digit (Type 26)
        INTEGER           , ALLOCATABLE, PUBLIC :: ICTL26( :,: )
        CHARACTER(SPNLEN3), ALLOCATABLE, PUBLIC :: CSPT26( :,: )
        CHARACTER(SICLEN3), ALLOCATABLE, PUBLIC :: CHRT26( : )

!.........  FIPS code = 0, SIC = all (Type 27)
        INTEGER           , ALLOCATABLE, PUBLIC :: ICTL27( :,: )
        CHARACTER(SPNLEN3), ALLOCATABLE, PUBLIC :: CSPT27( :,: )
        CHARACTER(SICLEN3), ALLOCATABLE, PUBLIC :: CHRT27( : )

!.........  FIPS code = state, SIC = 2-digit (Type 28)
        INTEGER           , ALLOCATABLE, PUBLIC :: ICTL28( :,: )
        CHARACTER(SPNLEN3), ALLOCATABLE, PUBLIC :: CSPT28( :,: )
        CHARACTER(STILEN3), ALLOCATABLE, PUBLIC :: CHRT28( : )

!.........  FIPS code = state, SIC = all (Type 29)
        INTEGER           , ALLOCATABLE, PUBLIC :: ICTL29( :,: )
        CHARACTER(SPNLEN3), ALLOCATABLE, PUBLIC :: CSPT29( :,: )
        CHARACTER(STILEN3), ALLOCATABLE, PUBLIC :: CHRT29( : )

!.........  FIPS code = all, SIC = 2-digit (Type 30)
        INTEGER           , ALLOCATABLE, PUBLIC :: ICTL30( :,: )
        CHARACTER(SPNLEN3), ALLOCATABLE, PUBLIC :: CSPT30( :,: )
        CHARACTER(FPILEN3), ALLOCATABLE, PUBLIC :: CHRT30( : )

!.........  FIPS code = all, SIC = all (Type 31)
        INTEGER           , ALLOCATABLE, PUBLIC :: ICTL31( :,: )
        CHARACTER(SPNLEN3), ALLOCATABLE, PUBLIC :: CSPT31( :,: )
        CHARACTER(FPILEN3), ALLOCATABLE, PUBLIC :: CHRT31( : )

!.........  Additional groups for MACT processing

!.........  FIPS code = 0, SCC = 0, MACT = all (Type 32)
        INTEGER           , ALLOCATABLE, PUBLIC :: ICTL32( :,: )
        CHARACTER(SPNLEN3), ALLOCATABLE, PUBLIC :: CSPT32( :,: )
        CHARACTER(MACLEN3), ALLOCATABLE, PUBLIC :: CHRT32( : )
        
!.........  FIPS code = 0, SCC = all, MACT = all (Type 33)
        INTEGER           , ALLOCATABLE, PUBLIC :: ICTL33( :,: )
        CHARACTER(SPNLEN3), ALLOCATABLE, PUBLIC :: CSPT33( :,: )
        CHARACTER(MSCLEN3), ALLOCATABLE, PUBLIC :: CHRT33( : )
        
!.........  FIPS code = state, SCC = 0, MACT = all (Type 34)
        INTEGER           , ALLOCATABLE, PUBLIC :: ICTL34( :,: )
        CHARACTER(SPNLEN3), ALLOCATABLE, PUBLIC :: CSPT34( :,: )
        CHARACTER(MSTLEN3), ALLOCATABLE, PUBLIC :: CHRT34( : )
        
!.........  FIPS code = state, SCC = all, MACT = all (Type 35)
        INTEGER           , ALLOCATABLE, PUBLIC :: ICTL35( :,: )
        CHARACTER(SPNLEN3), ALLOCATABLE, PUBLIC :: CSPT35( :,: )
        CHARACTER(MSSLEN3), ALLOCATABLE, PUBLIC :: CHRT35( : )
        
!.........  FIPS code = all, SCC = 0, MACT = all (Type 36)
        INTEGER           , ALLOCATABLE, PUBLIC :: ICTL36( :,: )
        CHARACTER(SPNLEN3), ALLOCATABLE, PUBLIC :: CSPT36( :,: )
        CHARACTER(MFPLEN3), ALLOCATABLE, PUBLIC :: CHRT36( : )
        
!.........  FIPS code = all, SCC = all, MACT = all (Type 37)
        INTEGER           , ALLOCATABLE, PUBLIC :: ICTL37( :,: )
        CHARACTER(SPNLEN3), ALLOCATABLE, PUBLIC :: CSPT37( :,: )
        CHARACTER(MFSLEN3), ALLOCATABLE, PUBLIC :: CHRT37( : )

!.........  FIPS code = all, Plant=non-blank, SCC=0, MACT = all (Type 38)
        INTEGER           , ALLOCATABLE, PUBLIC :: ICTL38( :,: )
        CHARACTER(SPNLEN3), ALLOCATABLE, PUBLIC :: CSPT38( :,: )
        CHARACTER(FPMLEN3), ALLOCATABLE, PUBLIC :: CHRT38( : )

!.........  Unsorted, unprocessed cross-reference arrays
        INTEGER, ALLOCATABLE, PUBLIC:: INDXTA ( : ) !  sorting index
        INTEGER, ALLOCATABLE, PUBLIC:: ISPTA  ( : ) !  pollutant index to EINAM
        INTEGER, ALLOCATABLE, PUBLIC:: ISRGCDA( : ) !  spatial surrogate codes
        INTEGER, ALLOCATABLE, PUBLIC:: SRGIDAA( : ) !  assigned spatial surrogate codes
        INTEGER, ALLOCATABLE, PUBLIC:: MPRNA  ( : ) !  monthly profile codes
        INTEGER, ALLOCATABLE, PUBLIC:: WPRNA  ( : ) !  weekly profile codes
        INTEGER, ALLOCATABLE, PUBLIC:: DPRNA  ( : ) !  diurnal profile codes
        INTEGER, ALLOCATABLE, PUBLIC:: IPSIA( :,: ) !  24 hours code for EFs
        INTEGER, ALLOCATABLE, PUBLIC:: IARPTA( :,: )!  tbl num, row, & count
        INTEGER, ALLOCATABLE, PUBLIC:: ISPDCDA( : ) !  speed profile codes

        CHARACTER(FIPLEN3), ALLOCATABLE, PUBLIC:: CFIPTA( : ) ! FIPS codes
        CHARACTER(SPNLEN3), ALLOCATABLE, PUBLIC:: CSPRNA( : ) ! spec prof #
        CHARACTER(SCCLEN3), ALLOCATABLE, PUBLIC:: CSCCTA( : ) ! SCC
        CHARACTER(MACLEN3), ALLOCATABLE, PUBLIC:: CMACTA( : ) ! MACT
        CHARACTER(SICLEN3), ALLOCATABLE, PUBLIC:: CISICA( : ) ! SIC
        CHARACTER(TAGLEN3), ALLOCATABLE, PUBLIC:: CTAGNA( : ) ! Tag labels
        CHARACTER(IOVLEN3+TAGLEN3+1),
     &                      ALLOCATABLE, PUBLIC:: CSPCTAGNA( : ) ! Species name - Tag labels

        CHARACTER(SSMLEN3+POLLEN3), ALLOCATABLE, PUBLIC:: CSRCTA(:)
                                          ! source chars // SCC // MACT // SIC // pollutant

        END MODULE MODXREF
