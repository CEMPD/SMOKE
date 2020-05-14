
        MODULE MODSOURC

!***********************************************************************
!  Module body starts at line 41
!
!  DESCRIPTION:
!     This module contains the public allocatable arrays for the source
!     characteristics (both sorted and unsorted)
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

        INCLUDE 'EMPRVT3.EXT'   !  emissions private parameters

!.........  Flag to indicate whether source is intergrated or not
        LOGICAL, SAVE :: INTGRFLAG = .FALSE.

!.........  Sorted list of point sources for SMOKE inventory file
        INTEGER, ALLOCATABLE, PUBLIC:: IRCLAS( : )  !  road class number
        INTEGER, ALLOCATABLE, PUBLIC:: IVTYPE( : )  !  vehicle type code
        INTEGER, ALLOCATABLE, PUBLIC:: CELLID( : )  !  Cell ID
        INTEGER, POINTER,     PUBLIC:: IPOSCOD( : ) !  positn of pol in INVPCOD
        INTEGER, ALLOCATABLE, PUBLIC:: TZONES( : )  !  time zones
        INTEGER, POINTER,     PUBLIC:: TPFLAG( : )  !  temporal profile types
        INTEGER, POINTER,     PUBLIC:: INVYR ( : )  !  inv year for this record
        INTEGER, POINTER,     PUBLIC:: NPCNT ( : )  !  No. of pols per raw rec
        INTEGER, ALLOCATABLE, PUBLIC:: FLTRDAYL( : )!  daylight time filter
        INTEGER, ALLOCATABLE, PUBLIC:: SRGID ( :,: )!  primary & fallbk surg ID
        INTEGER, ALLOCATABLE, PUBLIC:: IDIU  ( : )  !  Hr prof code per source
        INTEGER, ALLOCATABLE, PUBLIC:: IWEK  ( : )  !  Wk prof code per source
        INTEGER, ALLOCATABLE, PUBLIC:: IMON  ( : )  !  Mn prof code per source

        REAL*8 , ALLOCATABLE, PUBLIC:: XLOCA ( : )  !  lon X-location
        REAL*8 , ALLOCATABLE, PUBLIC:: YLOCA ( : )  !  lat Y-location
        REAL*8 , ALLOCATABLE, PUBLIC:: XLOC1 ( : )  !  lon X-location link start
        REAL*8 , ALLOCATABLE, PUBLIC:: YLOC1 ( : )  !  lat Y-location link start
        REAL*8 , ALLOCATABLE, PUBLIC:: XLOC2 ( : )  !  lon X-location link end
        REAL*8 , ALLOCATABLE, PUBLIC:: YLOC2 ( : )  !  lat Y-location link end
        REAL   , ALLOCATABLE, PUBLIC:: SPEED ( : )  !  speed
        REAL   , ALLOCATABLE, PUBLIC:: STKHT ( : )  !  stack height   (m)
        REAL   , ALLOCATABLE, PUBLIC:: STKDM ( : )  !  stack diameter (m)
        REAL   , ALLOCATABLE, PUBLIC:: STKTK ( : )  !  exhaust temp   (deg K)
        REAL   , ALLOCATABLE, PUBLIC:: STKVE ( : )  !  exhaust veloc  (m/s)
        REAL   , ALLOCATABLE, PUBLIC:: VPOP  ( : )  !  vehicle population
        REAL   , ALLOCATABLE, PUBLIC:: VMT   ( : )  !  vehicle miles traveled (miles/day)
        REAL   , ALLOCATABLE, PUBLIC:: HEATCONTENT ( : )  !  heatcontent constant (m/s)

        REAL   , POINTER,     PUBLIC:: POLVAL( :,: )!  pol-spec values by pol

        CHARACTER(FIPLEN3), POINTER,     PUBLIC:: CIFIP  ( : ) ! FIPS code
        CHARACTER(SCCLEN3), POINTER,     PUBLIC:: CSCC   ( : ) ! SCC
        CHARACTER(EXTLEN3), POINTER,     PUBLIC:: CEXTORL( : ) ! Additional Extended ORL vars
        CHARACTER(NEILEN3), ALLOCATABLE, PUBLIC:: CNEIUID( : ) ! NEI Unique ID
        CHARACTER(ORSLEN3), ALLOCATABLE, PUBLIC:: CORIS  ( : ) ! DOE plant ID
        CHARACTER(BLRLEN3), ALLOCATABLE, PUBLIC:: CBLRID ( : ) ! boiler ID
        CHARACTER(LNKLEN3), ALLOCATABLE, PUBLIC:: CLINK  ( : ) ! link
        CHARACTER(DSCLEN3), ALLOCATABLE, PUBLIC:: CPDESC ( : ) ! plant desc
        CHARACTER(ALLLEN3), POINTER,     PUBLIC:: CSOURC ( : ) ! concat src
        CHARACTER(VTPLEN3), ALLOCATABLE, PUBLIC:: CVTYPE ( : ) ! vehicle type
        CHARACTER(INTLEN3), POINTER,     PUBLIC:: CINTGR ( : ) ! integrate status
        CHARACTER(ERPLEN3), ALLOCATABLE, PUBLIC:: CERPTYP( : ) ! emission release point type
        CHARACTER(MACLEN3), POINTER,     PUBLIC:: CMACT  ( : ) ! MACT code
        CHARACTER(NAILEN3), POINTER,     PUBLIC:: CNAICS ( : ) ! NAICS code
        CHARACTER(STPLEN3), POINTER,     PUBLIC:: CSRCTYP( : ) ! source type code code
        CHARACTER(SICLEN3), POINTER,     PUBLIC:: CISIC  ( : ) ! SIC
        CHARACTER(SHPLEN3), POINTER,     PUBLIC:: CSHAPE ( : ) ! area-source SHAPE_ID

        !!  <source,pollutant> :: speciation-profiles/fractions matrix, read by RDSSUP()
        INTEGER           ,              PUBLIC:: NSPFRC       ! total # of profiles&fractions
        INTEGER           , ALLOCATABLE, PUBLIC:: SPPNLO( : )  ! spec prof counts    (nsrc)
        INTEGER           , ALLOCATABLE, PUBLIC:: SPPNHI( : )  ! spec prof counts    (nsrc)
        CHARACTER(SPNLEN3), ALLOCATABLE, PUBLIC:: SPPROF( : )  ! spec prof codes     (nspfrc)
        REAL              , ALLOCATABLE, PUBLIC:: SPFRAC( : )  ! spec prof fractions (nspfrc)

        CHARACTER(TMPLEN3), ALLOCATABLE, PUBLIC:: CMON   ( : ) ! monthly profile code
        CHARACTER(TMPLEN3), ALLOCATABLE, PUBLIC:: CWEK   ( : ) ! weekly profile code
        CHARACTER(TMPLEN3), ALLOCATABLE, PUBLIC:: CDOM   ( : ) ! day of month profile code
        CHARACTER(TMPLEN3), ALLOCATABLE, PUBLIC:: CMND   ( : ) ! Monday profile code
        CHARACTER(TMPLEN3), ALLOCATABLE, PUBLIC:: CTUE   ( : ) ! Tuesday profile code
        CHARACTER(TMPLEN3), ALLOCATABLE, PUBLIC:: CWED   ( : ) ! Wednesday profile code
        CHARACTER(TMPLEN3), ALLOCATABLE, PUBLIC:: CTHU   ( : ) ! Thursday profile code
        CHARACTER(TMPLEN3), ALLOCATABLE, PUBLIC:: CFRI   ( : ) ! Friday profile code
        CHARACTER(TMPLEN3), ALLOCATABLE, PUBLIC:: CSAT   ( : ) ! Saturday profile code
        CHARACTER(TMPLEN3), ALLOCATABLE, PUBLIC:: CSUN   ( : ) ! sunday profile code
        CHARACTER(TMPLEN3), ALLOCATABLE, PUBLIC:: CMET   ( : ) ! sunday profile code

!.........  Unsorted list of point sources for SMOKE inventory file
        INTEGER, POINTER,     PUBLIC:: INDEXA( : ) !  subscript table for SORTIC
        INTEGER, ALLOCATABLE, PUBLIC:: IRCLASA( : )!  road class number
        INTEGER, ALLOCATABLE, PUBLIC:: IVTYPEA( : )!  vehicle type code
        INTEGER, POINTER,     PUBLIC:: IPOSCODA(:) !  positn of pol in INVPCOD
        INTEGER, POINTER,     PUBLIC:: ICASCODA(:) !  positn of CAS num. in UNIQCAS
        INTEGER, POINTER,     PUBLIC:: TPFLGA( : ) !  temporal resolution code
        INTEGER, POINTER,     PUBLIC:: INVYRA( : ) !  inventory year
        INTEGER, ALLOCATABLE, PUBLIC:: IDIUA ( : ) !  Hrly prof code per source
        INTEGER, ALLOCATABLE, PUBLIC:: IWEKA ( : ) !  Wkly prof code per source
        INTEGER, POINTER,     PUBLIC:: INRECA( : ) !  Input record per src x pol
        INTEGER, POINTER,     PUBLIC:: SRCIDA( : ) !  Source ID

        REAL*8 , POINTER,     PUBLIC:: XLOCAA( : ) !  UTM X-location (m)
        REAL*8 , POINTER,     PUBLIC:: YLOCAA( : ) !  UTM Y-location (m)
        REAL*8 , ALLOCATABLE, PUBLIC:: XLOC1A( : ) !  lon X-location link start
        REAL*8 , ALLOCATABLE, PUBLIC:: YLOC1A( : ) !  lat Y-location link start
        REAL*8 , ALLOCATABLE, PUBLIC:: XLOC2A( : ) !  lon X-location link end
        REAL*8 , ALLOCATABLE, PUBLIC:: YLOC2A( : ) !  lat Y-location link end
        REAL   , ALLOCATABLE, PUBLIC:: SPEEDA( : ) !  speed
        REAL   , ALLOCATABLE, PUBLIC:: STKHTA( : ) !  stack height   (m)
        REAL   , ALLOCATABLE, PUBLIC:: STKDMA( : ) !  stack diameter (m)
        REAL   , ALLOCATABLE, PUBLIC:: STKTKA( : ) !  exhaust temperature (deg K)
        REAL   , ALLOCATABLE, PUBLIC:: STKVEA( : ) !  exhaust velocity    (m/s)
        REAL   , POINTER,     PUBLIC:: POLVLA( :,: )! emis-spec values. See BLDENAMS.
        REAL   , ALLOCATABLE, PUBLIC:: VMTA  ( : ) !  vehicle miles traveled

        REAL   , ALLOCATABLE, PUBLIC:: FUGHGT( : ) !  fugitive emissions height
        REAL   , ALLOCATABLE, PUBLIC:: FUGWID( : ) !  fugitive emissions width (YDIM)
        REAL   , ALLOCATABLE, PUBLIC:: FUGLEN( : ) !  fugitive emissions length (XDIM)
        REAL   , ALLOCATABLE, PUBLIC:: FUGANG( : ) !  fugitive emissions angle

        CHARACTER(SCCLEN3), POINTER,     PUBLIC:: CSCCA  ( : ) ! SCC
        CHARACTER(ORSLEN3), ALLOCATABLE, PUBLIC:: CORISA ( : ) ! DOE plant ID
        CHARACTER(BLRLEN3), ALLOCATABLE, PUBLIC:: CBLRIDA( : ) ! boiler ID
        CHARACTER(LNKLEN3), ALLOCATABLE, PUBLIC:: CLINKA ( : ) ! link
        CHARACTER(DSCLEN3), ALLOCATABLE, PUBLIC:: CPDESCA( : ) ! plant desc
        CHARACTER(ALLCAS3), POINTER,     PUBLIC:: CSOURCA( : ) ! concat src
        CHARACTER(VTPLEN3), ALLOCATABLE, PUBLIC:: CVTYPEA( : ) ! vehicle type

!.........  MEDS-gridded inventory related arrays
        INTEGER,              PUBLIC:: NMEDGRD
        INTEGER,              PUBLIC:: NMEDGAI

        CHARACTER(CHRLEN3), ALLOCATABLE, PUBLIC:: CMEDGRD( :,: )   ! MEDS grid row/col coord
        CHARACTER(FIPLEN3), ALLOCATABLE, PUBLIC:: COABDST( :,: )   ! MEDS GAI CO-ABS-DIST code

!.........  Unsorted list of file numbers and records by source
        INTEGER, PUBLIC :: NSTRECS                      ! size of SRCSBYREC
        INTEGER, ALLOCATABLE, PUBLIC:: SRCSBYREC( :,: ) ! file number, record number, and
                                                        ! src number for each inventory record
        INTEGER, ALLOCATABLE, PUBLIC:: RECIDX( : )      ! index for SRCSBYREC

        END MODULE MODSOURC
