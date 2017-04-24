
        MODULE MODTMPRL

!***********************************************************************
!  Module body starts at line
!
!  DESCRIPTION:
!     This module contains temporal allocation information.
!
!  PRECONDITIONS REQUIRED:
!
!  SUBROUTINES AND FUNCTIONS CALLED:
!
!  REVISION HISTORY:
!     Created 1/99 by M. Houyoux for temporal profiles only
!     Changed 10/2000 by MRH to include additional temporal info.
!     Modified 2012 by BH Baek for new Gentpro met-based temporal profiles
!     Modified 07/2014 by C.Coats for new GENTPRO CSV profiles and cross-references
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

C.........  Hourly-emissions file information
        INTEGER, PUBLIC :: NTPDAT   ! No. data values in hourly emissions file

        LOGICAL, PUBLIC :: LTFLAG = .FALSE.    ! flag for output hourly emissions in local time

        CHARACTER(IOVLEN3), ALLOCATABLE, PUBLIC :: TPNAME( : )  ! data names
        CHARACTER(IOULEN3), ALLOCATABLE, PUBLIC :: TPUNIT( : )  ! data units
        CHARACTER(IODLEN3), ALLOCATABLE, PUBLIC :: TPDESC( : )  ! data descs

C.........  Holiday dates arrays

        INTEGER, PUBLIC :: NHOLIDAY  ! number of holidays

        INTEGER, ALLOCATABLE, PUBLIC :: HOLREGN ( : ) ! Region code of holiday
        INTEGER, ALLOCATABLE, PUBLIC :: HOLJDATE( : ) ! Julian date of holidays
        INTEGER, ALLOCATABLE, PUBLIC :: HOLALTDY( : ) ! alternative day of week

C.........  Sorted temporal profiles

        INTEGER,        PUBLIC :: NMETPROF                ! number of Met-based temp profiles
        LOGICAL,        PUBLIC :: METPROFLAG = .FALSE.    ! flag for Met-based temp profiles
        CHARACTER( 8 ), PUBLIC :: METPROTYPE = ' '        ! type of Met-based temp profiles
        CHARACTER( 8 ), PUBLIC :: HOUR_TPROF = ' '        ! hour of day, hour of month, hour of year

        INTEGER, PUBLIC :: NMON   ! number of month-of-year profiles
        INTEGER, PUBLIC :: NWEK   ! number of day-of-week   profiles
        INTEGER, PUBLIC :: NHRL   ! number of hour-of-day   profiles
        INTEGER, PUBLIC :: NDOM   ! number of day-of-month  profiles

        LOGICAL, PUBLIC :: MONFLAG( 12 )       !  true iff this month   active for SDATE:EDATE
        LOGICAL, PUBLIC :: DAYFLAG(  7 )       !  true iff this  day    active for SDATE:EDATE (as for I/O API WKDAY())
        LOGICAL, PUBLIC :: WKDFLAG             !  true iff some weekday active for SDATE:EDATE (as for I/O API WKDAY())
        LOGICAL, PUBLIC :: WKEFLAG             !  true iff some weekend active for SDATE:EDATE (as for I/O API WKDAY())

        INTEGER, ALLOCATABLE, PUBLIC :: IPOL2D( :,: )   !  (NGSZ,NGRP) EANAM-subscripts per group, species

        LOGICAL, ALLOCATABLE, PUBLIC :: POLREFFLAG( : ) ! (0:NIPPA) -- does pollutant exist in XREF file
        LOGICAL, ALLOCATABLE, PUBLIC :: METREFFLAG( : ) ! (0:NIPPA) -- does pollutant exist in XREF file

C.............  Sorted profile data structures filtered from the input files

        CHARACTER(16), ALLOCATABLE :: MTHIDP( : )       !  month-of-year profile IDs for TPRO_MONTHLY
        CHARACTER(16), ALLOCATABLE :: WEKIDP( : )       !  day-of-week   profile IDs for TPRO_WEEKLY
        CHARACTER(16), ALLOCATABLE :: DOMIDP( : )       !  day-of-month  profile IDs for TPRO_DAILY
        CHARACTER(16), ALLOCATABLE :: HRLIDP( : )       !  hour-of-day   profile IDs for TPRO_HOURLY

        INTEGER, ALLOCATABLE, PUBLIC :: MONREF( : )   ! (NMON) month-of-year codes
        INTEGER, ALLOCATABLE, PUBLIC :: WEKREF( : )   ! (NWEK) day-of-week   codes
        INTEGER, ALLOCATABLE, PUBLIC :: HRLREF( : )   ! (NHRL) hour-of-day   codes
        INTEGER, ALLOCATABLE, PUBLIC :: DAYREF( : )   ! (NDOM) day-of-month  codes

        INTEGER, PUBLIC :: MTHCOUNT             !  counts for each XREF category:  month-of-year
        INTEGER, PUBLIC :: WEKCOUNT             !  day-of-week
        INTEGER, PUBLIC :: DOMCOUNT             !  day-of-month
        INTEGER, PUBLIC :: MONCOUNT             !  Monday hour-of-day
        INTEGER, PUBLIC :: TUECOUNT
        INTEGER, PUBLIC :: WEDCOUNT
        INTEGER, PUBLIC :: THUCOUNT
        INTEGER, PUBLIC :: FRICOUNT
        INTEGER, PUBLIC :: SATCOUNT
        INTEGER, PUBLIC :: SUNCOUNT
        INTEGER, PUBLIC :: METCOUNT

        INTEGER, ALLOCATABLE, PUBLIC :: MTHPDEX( : )        ! (MTHCOUNT)  profile-subscripts per XREF category
        INTEGER, ALLOCATABLE, PUBLIC :: WEKPDEX( : )        ! (WEKCOUNT)
        INTEGER, ALLOCATABLE, PUBLIC :: DOMPDEX( : )        ! etc...
        INTEGER, ALLOCATABLE, PUBLIC :: MONPDEX( : )
        INTEGER, ALLOCATABLE, PUBLIC :: TUEPDEX( : )
        INTEGER, ALLOCATABLE, PUBLIC :: WEDPDEX( : )
        INTEGER, ALLOCATABLE, PUBLIC :: THUPDEX( : )
        INTEGER, ALLOCATABLE, PUBLIC :: FRIPDEX( : )
        INTEGER, ALLOCATABLE, PUBLIC :: SATPDEX( : )
        INTEGER, ALLOCATABLE, PUBLIC :: SUNPDEX( : )

        CHARACTER(ALLLEN3), ALLOCATABLE, PUBLIC :: MTHKEYS( : )     !  (MTHCOUNT) sorted cross-reference keys  per XREF category
        CHARACTER(ALLLEN3), ALLOCATABLE, PUBLIC :: WEKKEYS( : )
        CHARACTER(ALLLEN3), ALLOCATABLE, PUBLIC :: DOMKEYS( : )
        CHARACTER(ALLLEN3), ALLOCATABLE, PUBLIC :: MONKEYS( : )
        CHARACTER(ALLLEN3), ALLOCATABLE, PUBLIC :: TUEKEYS( : )
        CHARACTER(ALLLEN3), ALLOCATABLE, PUBLIC :: WEDKEYS( : )
        CHARACTER(ALLLEN3), ALLOCATABLE, PUBLIC :: THUKEYS( : )
        CHARACTER(ALLLEN3), ALLOCATABLE, PUBLIC :: FRIKEYS( : )
        CHARACTER(ALLLEN3), ALLOCATABLE, PUBLIC :: SATKEYS( : )
        CHARACTER(ALLLEN3), ALLOCATABLE, PUBLIC :: SUNKEYS( : )
        CHARACTER(ALLLEN3), ALLOCATABLE, PUBLIC :: METKEYS( : )

        INTEGER, ALLOCATABLE, PUBLIC :: METPROF( :,: )          !  (NSRC,NIPPA) METFACS-subscripts per source, species
        INTEGER, ALLOCATABLE, PUBLIC :: MTHPROF( :,: )          !  (NSRC,NIPPA) MONFAC-subscripts per source, species
        INTEGER, ALLOCATABLE, PUBLIC :: WEKPROF( :,: )          !  (NSRC,NIPPA) WEKFAC-subscripts, or 0
        INTEGER, ALLOCATABLE, PUBLIC :: DOMPROF( :,: )          !  (NSRC,NIPPA) DOMFAC-subscripts, or 0
        INTEGER, ALLOCATABLE, PUBLIC :: HRLPROF( :,:,: )        !  (NSRC,7,NIPPA) HRLFAC-subscripts per source, species

        INTEGER, ALLOCATABLE, PUBLIC :: STDATE( : )     ! Start date
        INTEGER, ALLOCATABLE, PUBLIC :: STTIME( : )     ! Start time
        INTEGER, ALLOCATABLE, PUBLIC :: RUNLEN( : )     ! Run hours
        INTEGER, ALLOCATABLE, PUBLIC :: ITDATE( : )     ! Start julian date

        REAL   , ALLOCATABLE, PUBLIC :: MONFAC_ORG( :,: )   ! (12,NMON) Org monthly factors before monthly adjustment
        REAL   , ALLOCATABLE, PUBLIC :: MONFAC( :,: )       ! (12,NMON)    month-of-year factors
        REAL   , ALLOCATABLE, PUBLIC :: WEKFAC( :,: )       ! ( 7,NWEK)    day-of-week facs (week-norm)
        REAL   , ALLOCATABLE, PUBLIC :: XWKFAC( :,: )       ! ( 7,NWEK)    day-of-week facs (wkday-norm)
        REAL   , ALLOCATABLE, PUBLIC :: DOMFAC( :,:,: )     ! (31,12,NDOM) day-of-month facs
        REAL   , ALLOCATABLE, PUBLIC :: HRLFAC( :,: )       ! (24,NHRL)    hour-of-day factors
        REAL   , ALLOCATABLE, PUBLIC :: METFACS( :,:,: )    ! Met-based temporal profiles

        END MODULE MODTMPRL
