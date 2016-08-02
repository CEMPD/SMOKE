
        MODULE MODLISTS

!***********************************************************************
!  Module body starts at line
!
!  DESCRIPTION:
!     This module contains the public allocatable arrays for various lists
!     from the inventory, such as for SCCs, SICs, inventory pollutants, 
!     and other criteria for selecting the parts of cross-reference files that
!     are relevant.
!
!  PRECONDITIONS REQUIRED:
!
!  SUBROUTINES AND FUNCTIONS CALLED:
!
!  REVISION HISTORY:
!     Created 3/99 by M. Houyoux
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

!.........  Sizes of the arrays...

        INTEGER, PUBLIC :: NINVSIC   ! no. unique SICs in inventory
        INTEGER, PUBLIC :: NINVSIC2  ! no. unique 2-digit SICs in inventory
        INTEGER, PUBLIC :: NINVSCC   ! no. unique SCCs in inventory
        INTEGER, PUBLIC :: NINVSCL   ! no. unique left SCCs in inventory
        INTEGER, PUBLIC :: NINVIFIP  ! no. unique country/state/county codes
        INTEGER, PUBLIC :: NINVVTYP  ! no. unique vehicle types in inventory
        INTEGER, PUBLIC :: NINVMACT  ! no. unique MACTs in inventory
        INTEGER, PUBLIC :: NINVNAICS ! no. unique NAICS in inventory
        INTEGER, PUBLIC :: NINVORIS  ! no. unique ORIS 
        INTEGER, PUBLIC :: NORISBLR  ! no. unique ORIS // boilers
        INTEGER, PUBLIC :: NOBLRSRC  ! no. sources matching ORIS/boiler combinations

!.........  Controllers
        LOGICAL, PUBLIC :: ORISFLAG  ! true: create ORIS-based arrays

!.........  Information about the inventory files
        INTEGER,        ALLOCATABLE :: FILFMT( : )  ! format of inventory file(s)
        CHARACTER(300), ALLOCATABLE :: LSTSTR( : )  ! contents of list-fmt inventory file

!.........  Inventory format flags
        LOGICAL, PUBLIC :: FIREFF10  ! true: FF10 FIRE format inventory
        LOGICAL, PUBLIC :: FIREFLAG  ! true: ORL FIRE format inventory
        LOGICAL, PUBLIC :: FF10FLAG  ! true: FF10 format inventories
        LOGICAL, PUBLIC :: MEDSFLAG  ! true: Pregridded FF10 format inventories
        LOGICAL, PUBLIC :: APIFLAG   ! true: Pregridded IOAPI format inventories
        LOGICAL, PUBLIC :: NCDFLAG   ! true: Pregridded raw NetCDF format inventories

!.........  Unique lists of source characteristics and associated arrays...

!.........  SIC arrays dimensioned by NINVSIC
        CHARACTER(SICLEN3), ALLOCATABLE, PUBLIC :: INVSIC ( : ) ! SICs
        CHARACTER(SICLEN3), ALLOCATABLE, PUBLIC :: INVSIC2( : ) ! 2-digit SICs
        CHARACTER(SDSLEN3), ALLOCATABLE, PUBLIC :: SICDESC( : ) ! descrptn

!.........  SCC arrays dimensioned by NINVSCC
        INTEGER, ALLOCATABLE, PUBLIC :: SCCDLEV( :,: ) ! SCC desc length for levels
        CHARACTER(SCCLEN3), ALLOCATABLE, PUBLIC :: INVSCC( : ) ! SCCs
        CHARACTER(SCCLEN3), ALLOCATABLE, PUBLIC :: INVSCL( : ) ! left SCCs
        CHARACTER(SDSLEN3), ALLOCATABLE, PUBLIC :: SCCDESC( : ) ! descrptn

!.........  Country/state/county codes dimensioned by NINVIFIP
        CHARACTER(FIPLEN3), ALLOCATABLE, PUBLIC :: INVCFIP( : )

!.........  Vehicle types dimensioned by NINVVTYP
        CHARACTER(VTPLEN3), ALLOCATABLE, PUBLIC :: INVVTYP( : )

!.........  MACT codes dimensioned by NINVMACT
        CHARACTER(MACLEN3), ALLOCATABLE, PUBLIC :: INVMACT( : ) ! MACT codes
        CHARACTER(SDSLEN3), ALLOCATABLE, PUBLIC :: MACTDESC( : ) ! descrptn

!.........  NAICS codes dimensioned by NINVNAICS
        CHARACTER(NAILEN3), ALLOCATABLE, PUBLIC :: INVNAICS( : ) ! NAICS codes
        CHARACTER(SDSLEN3), ALLOCATABLE, PUBLIC :: NAICSDESC( : ) ! descrptn

!.........  ORIS arrays
        INTEGER,            ALLOCATABLE, PUBLIC :: OBSRCBG( : ) ! 1st source per ORIS/boiler
        INTEGER,            ALLOCATABLE, PUBLIC :: OBSRCNT( : ) ! source count per ORIS/boiler
        INTEGER,            ALLOCATABLE, PUBLIC :: OBSRCNM( : ) ! src numbers matching ORIS/boiler
        LOGICAL,            ALLOCATABLE, PUBLIC :: IORSMTCH( : ) ! true: inventory ORIS matched to CEM
        CHARACTER(ORSLEN3), ALLOCATABLE, PUBLIC :: INVORIS( : ) ! unique ORIS
        CHARACTER(FIPLEN3), ALLOCATABLE, PUBLIC :: INVORFP( : ) ! FIPS for ORIS in inventory
        CHARACTER(DSCLEN3), ALLOCATABLE, PUBLIC :: INVODSC( : ) ! plant description from inventory
        CHARACTER(OBRLEN3), ALLOCATABLE, PUBLIC :: ORISBLR( : ) ! ORIS // boiler

!.........  For valid pollutants and activities...
        INTEGER, PUBLIC :: MXIDAT = 0   ! Max no of inv pols & acvtys
        INTEGER, PUBLIC :: NINVTBL = 0  ! Number of entries in inventory table
        INTEGER, PUBLIC :: NUNIQCAS = 0 ! Number of unique CAS codes
        INTEGER, PUBLIC :: NINVKEEP = 0 ! Number of kept INVTABLE entries (KEEP=Y)
        INTEGER, PUBLIC :: NINVDROP = 0 ! Number of dropped INVTABLE entries (KEEP=N)

C.........  Full list of inventory pollutants/activities (in output order)
C.........  Dimensioned by MXIDAT
        INTEGER, ALLOCATABLE, PUBLIC :: INVDCOD( : ) ! 5-digit SPECIATE4 ID (if any)
        INTEGER, ALLOCATABLE, PUBLIC :: INVSTAT( : ) ! Status (<0 activity; >0 pol)

        REAL   , ALLOCATABLE, PUBLIC :: INVDCNV( : ) ! local conversion factor

        CHARACTER,          ALLOCATABLE, PUBLIC :: INVDVTS( : ) ! V=VOC, T=TOG, N=not
        CHARACTER(IOVLEN3), ALLOCATABLE, PUBLIC :: INVDNAM( : ) ! data name 
        CHARACTER(IOULEN3), ALLOCATABLE, PUBLIC :: INVDUNT( : ) ! units for SMOKE intmdt inventory
        CHARACTER(DDSLEN3), ALLOCATABLE, PUBLIC :: INVDDSC( : ) ! inventory data description

C.........  Inventory table arrays - unsorted raw data, dimensioned by NINVTBL
        INTEGER, ALLOCATABLE, PUBLIC :: ITIDXA ( : ) ! Sorting index 1
        INTEGER, ALLOCATABLE, PUBLIC :: ITIDXA2( : ) ! Sorting index 2
        INTEGER, ALLOCATABLE, PUBLIC :: ITLINNO( : ) ! Line number of input file for record
        INTEGER, ALLOCATABLE, PUBLIC :: ITCODA ( : ) ! 5-digit SPECIATE4 ID (if any)
        INTEGER, ALLOCATABLE, PUBLIC :: ITNTIA ( : ) ! NTI HAP number
        INTEGER, ALLOCATABLE, PUBLIC :: ITREAA ( : ) ! Reactivity group
        INTEGER, ALLOCATABLE, PUBLIC :: ITSTATA( : ) ! Status (<0 activity; >0 pol)

        LOGICAL, ALLOCATABLE, PUBLIC :: ITKEEPA( : ) ! true: keep record data
        LOGICAL, ALLOCATABLE, PUBLIC :: ITMSPC ( : ) ! true: pollutant is a model species
        LOGICAL, ALLOCATABLE, PUBLIC :: ITEXPL ( : ) ! true: pollutant is explicit in the mechanism

        REAL   , ALLOCATABLE, PUBLIC :: ITFACA ( : ) ! fraction of CAS emissions in pollutant

        CHARACTER,          ALLOCATABLE, PUBLIC :: ITVTSA ( : ) ! V=VOC, T=TOG, N=not
        CHARACTER(IOVLEN3), ALLOCATABLE, PUBLIC :: ITNAMA ( : ) ! data name 
        CHARACTER(IOULEN3), ALLOCATABLE, PUBLIC :: ITUNTA ( : ) ! units for SMOKE intmdt inventory
        CHARACTER(CASLEN3), ALLOCATABLE, PUBLIC :: ITCASA ( : ) ! CAS code (left justified)
        CHARACTER(DDSLEN3), ALLOCATABLE, PUBLIC :: ITDSCA ( : ) ! inventory data description
        CHARACTER(DDSLEN3), ALLOCATABLE, PUBLIC :: ITCASDSCA( : ) ! CAS code description
        CHARACTER(CDTLEN3), ALLOCATABLE, PUBLIC :: ITCASDNMA( : ) ! CAS code // data name

C.........  CAS codes in sorted order WITH DUPLICATES and reference index back to 
C           input order - dimensioned by NINVTBL
        INTEGER,            ALLOCATABLE, PUBLIC :: SCASIDX( : ) ! index to ITNAMA
        CHARACTER(CASLEN3), ALLOCATABLE, PUBLIC :: SORTCAS( : ) ! CAS code (left justified)
 
C.........  CAS codes in sorted order WITHOUT DUPLICATES and count of pollutants
C           for each CAS code (0 pollutants indicates that no pollutants for that
C           CAS are kept) - dimensioned by NUNIQCAS
        INTEGER,            ALLOCATABLE, PUBLIC :: UCASIDX ( : ) ! index to first entry in SORTCAS
        INTEGER,            ALLOCATABLE, PUBLIC :: UCASNPOL( : ) ! pol count per CAS code
        INTEGER,            ALLOCATABLE, PUBLIC :: UCASNKEP( : ) ! kept pol count per CAS code
        CHARACTER(CASLEN3), ALLOCATABLE, PUBLIC :: UNIQCAS ( : ) ! CAS code (left justified)

C.........  SPECIATE4 IDs in sorted order and reference index back to input order
C           Dimensioned by MXIDAT
        INTEGER, ALLOCATABLE, PUBLIC :: IDXCOD ( : ) ! index to INVDNAM
        INTEGER, ALLOCATABLE, PUBLIC :: SORTCOD( : ) ! SPECIATE4 ID

C.........  Arrays for per CAS reports - dimensioned by NUNIQCAS
C           The array order is the same as UNIQCAS
        INTEGER, ALLOCATABLE, PUBLIC :: RECSBYCAS( : ) ! no. records per CAS number
        REAL   , ALLOCATABlE, PUBLIC :: EMISBYCAS( : ) ! total emissions per CAS number

C.........  Arrays for CAS/pollutant combo reports - dimensioned by NINVTBL
C           The array order is the same as SORTCAS
        REAL   , ALLOCATABLE, PUBLIC :: EMISBYPOL( : ) ! total emissions per CAS/pollutant combo

        END MODULE MODLISTS
