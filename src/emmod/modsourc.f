
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

        INCLUDE 'EMPRVT3.EXT'   !  emissions private parameters

!.........  Sorted list of point sources for SMOKE inventory file
        INTEGER, POINTER,     PUBLIC:: IFIP  ( : )  !  source FIPS (county) ID
        INTEGER, ALLOCATABLE, PUBLIC:: ISIC  ( : )  !  source SIC
        INTEGER, ALLOCATABLE, PUBLIC:: IRCLAS( : )  !  road class number
        INTEGER, ALLOCATABLE, PUBLIC:: IVTYPE( : )  !  vehicle type code
        INTEGER, ALLOCATABLE, PUBLIC:: CELLID( : )  !  Cell ID
        INTEGER, POINTER,     PUBLIC:: IPOSCOD( : ) !  positn of pol in INVPCOD
        INTEGER, ALLOCATABLE, PUBLIC:: TZONES( : )  !  time zones
        INTEGER, POINTER,     PUBLIC:: TPFLAG( : )  !  temporal profile types
        INTEGER, POINTER,     PUBLIC:: INVYR ( : )  !  inv year for this record
        INTEGER, ALLOCATABLE, PUBLIC:: IDIU  ( : )  !  Hr prof code per source
        INTEGER, ALLOCATABLE, PUBLIC:: IWEK  ( : )  !  Wk prof code per source
        INTEGER, ALLOCATABLE, PUBLIC:: IMON  ( : )  !  Mn prof code per source
        INTEGER, POINTER,     PUBLIC:: NPCNT ( : )  !  No. of pols per raw rec
        INTEGER, ALLOCATABLE, PUBLIC:: FLTRDAYL( : )!  daylight time filter
        INTEGER, ALLOCATABLE, PUBLIC:: SRGID ( :,: )!  primary & fallbk surg ID

        REAL   , ALLOCATABLE, PUBLIC:: XLOCA ( : )  !  lon X-location 
        REAL   , ALLOCATABLE, PUBLIC:: YLOCA ( : )  !  lat Y-location 
        REAL   , ALLOCATABLE, PUBLIC:: XLOC1 ( : )  !  lon X-location link start 
        REAL   , ALLOCATABLE, PUBLIC:: YLOC1 ( : )  !  lat Y-location link start
        REAL   , ALLOCATABLE, PUBLIC:: XLOC2 ( : )  !  lon X-location link end 
        REAL   , ALLOCATABLE, PUBLIC:: YLOC2 ( : )  !  lat Y-location link end
        REAL   , ALLOCATABLE, PUBLIC:: SPEED ( : )  !  speed
        REAL   , ALLOCATABLE, PUBLIC:: STKHT ( : )  !  stack height   (m)
        REAL   , ALLOCATABLE, PUBLIC:: STKDM ( : )  !  stack diameter (m)
        REAL   , ALLOCATABLE, PUBLIC:: STKTK ( : )  !  exhaust temp   (deg K)
        REAL   , ALLOCATABLE, PUBLIC:: STKVE ( : )  !  exhaust veloc  (m/s)
        REAL   , ALLOCATABLE, PUBLIC:: VMT   ( : )  !  vehicle miles traveled (miles/day)

        REAL   , POINTER,     PUBLIC:: POLVAL( :,: )!  pol-spec values by pol

        CHARACTER(LEN=SCCLEN3), POINTER,     PUBLIC:: CSCC  ( : ) ! SCC
        CHARACTER(LEN=ORSLEN3), ALLOCATABLE, PUBLIC:: CORIS ( : ) ! DOE plant ID
        CHARACTER(LEN=BLRLEN3), ALLOCATABLE, PUBLIC:: CBLRID( : ) ! boiler ID
        CHARACTER(LEN=LNKLEN3), ALLOCATABLE, PUBLIC:: CLINK ( : ) ! link
        CHARACTER(LEN=DSCLEN3), ALLOCATABLE, PUBLIC:: CPDESC( : ) ! plant desc
        CHARACTER(LEN=ALLLEN3), POINTER,     PUBLIC:: CSOURC( : ) ! concat src
        CHARACTER(LEN=VTPLEN3), ALLOCATABLE, PUBLIC:: CVTYPE( : ) ! vehicle type
        CHARACTER(LEN=SPNLEN3), ALLOCATABLE, PUBLIC:: SPPROF( :,: ) ! spec prof

!.........  Unsorted list of point sources for SMOKE inventory file
        INTEGER, POINTER,     PUBLIC:: INDEXA( : ) !  subscript table for SORTIC
        INTEGER, POINTER,     PUBLIC:: IFIPA ( : ) !  raw state/county FIPS code
        INTEGER, ALLOCATABLE, PUBLIC:: ISICA ( : ) !  raw SIC
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

        REAL   , POINTER,     PUBLIC:: XLOCAA( : ) !  UTM X-location (m)
        REAL   , POINTER,     PUBLIC:: YLOCAA( : ) !  UTM Y-location (m)
        REAL   , ALLOCATABLE, PUBLIC:: XLOC1A( : ) !  lon X-location link start 
        REAL   , ALLOCATABLE, PUBLIC:: YLOC1A( : ) !  lat Y-location link start
        REAL   , ALLOCATABLE, PUBLIC:: XLOC2A( : ) !  lon X-location link end 
        REAL   , ALLOCATABLE, PUBLIC:: YLOC2A( : ) !  lat Y-location link end
        REAL   , ALLOCATABLE, PUBLIC:: SPEEDA( : ) !  speed
        REAL   , ALLOCATABLE, PUBLIC:: STKHTA( : ) !  stack height   (m)
        REAL   , ALLOCATABLE, PUBLIC:: STKDMA( : ) !  stack diameter (m)
        REAL   , ALLOCATABLE, PUBLIC:: STKTKA( : ) !  exhaust temperature (deg K)
        REAL   , ALLOCATABLE, PUBLIC:: STKVEA( : ) !  exhaust velocity    (m/s)
        REAL   , POINTER,     PUBLIC:: POLVLA( :,: )! emis-spec values. See BLDENAMS.
        REAL   , ALLOCATABLE, PUBLIC:: VMTA  ( : ) !  vehicle miles traveled

        CHARACTER(LEN=SCCLEN3), POINTER,     PUBLIC:: CSCCA  ( : ) ! SCC
        CHARACTER(LEN=ORSLEN3), ALLOCATABLE, PUBLIC:: CORISA ( : ) ! DOE plant ID
        CHARACTER(LEN=BLRLEN3), ALLOCATABLE, PUBLIC:: CBLRIDA( : ) ! boiler ID
        CHARACTER(LEN=LNKLEN3), ALLOCATABLE, PUBLIC:: CLINKA ( : ) ! link
        CHARACTER(LEN=DSCLEN3), ALLOCATABLE, PUBLIC:: CPDESCA( : ) ! plant desc
        CHARACTER(LEN=ALLCAS3), POINTER,     PUBLIC:: CSOURCA( : ) ! concat src
        CHARACTER(LEN=VTPLEN3), ALLOCATABLE, PUBLIC:: CVTYPEA( : ) ! vehicle type

!.........  Unsorted list of file numbers and records by source
        INTEGER, PUBLIC :: NSTRECS                      ! size of SRCSBYREC
        INTEGER, ALLOCATABLE, PUBLIC:: SRCSBYREC( :,: ) ! file number, record number, and
                                                        ! src number for each inventory record
        INTEGER, ALLOCATABLE, PUBLIC:: RECIDX( : )      ! index for SRCSBYREC

        END MODULE MODSOURC
