        MODULE MODCNTRL

!***********************************************************************
!  Module body starts at line
!
!  DESCRIPTION:
!     This module contains the public allocatable arrays for control factors
!     and matrices
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

!.........  CONTROL PACKET DATA ARRAYS...

!.........  CTG-packet-specific data tables
        REAL   , ALLOCATABLE, PUBLIC:: CUTCTG ( : )  ! CTG Emissions cutoff
        REAL   , ALLOCATABLE, PUBLIC:: FACCTG ( : )  ! CTG control factor
        REAL   , ALLOCATABLE, PUBLIC:: FACMACT( : )  ! MACT control factor
        REAL   , ALLOCATABLE, PUBLIC:: FACRACT( : )  ! RACT control factor
        CHARACTER(120), ALLOCATABLE, PUBLIC:: CTGCOMT( : ) ! packet comment

!.........  CONTROL-packet-specific data tables
        INTEGER, ALLOCATABLE, PUBLIC:: ICTLEQUP( : )  ! prim. ctl equipment code
        CHARACTER(SICLEN3), ALLOCATABLE, PUBLIC:: CCTLSIC ( : ) ! SIC code
        REAL   , ALLOCATABLE, PUBLIC:: FACCEFF ( : )  ! control effcncy (0-100)
        REAL   , ALLOCATABLE, PUBLIC:: FACREFF ( : )  ! rule effectvnss (0-100)
        REAL   , ALLOCATABLE, PUBLIC:: FACRLPN ( : )  ! rule penetrtion (0-100)
        LOGICAL, ALLOCATABLE, PUBLIC:: CTLRPLC ( : )  ! replacement flag
                                                      !   true: replace, false: add
        CHARACTER(120), ALLOCATABLE, PUBLIC:: CTLCOMT( : ) ! packet comment

!.........  ALLOWABLE-packet-specific data tables
        CHARACTER(SICLEN3), ALLOCATABLE, PUBLIC:: CALWSIC ( : ) ! SIC code
        REAL   , ALLOCATABLE, PUBLIC:: FACALW  ( : )  ! allowable control factor
        REAL   , ALLOCATABLE, PUBLIC:: EMCAPALW( : )  ! emissions cap
        REAL   , ALLOCATABLE, PUBLIC:: EMREPALW( : )  ! replacement emissions
        CHARACTER(120), ALLOCATABLE, PUBLIC:: ALWCOMT( : ) ! packet comment

!.........  REACTIVITY-packet-specific data tables
        CHARACTER(SICLEN3), ALLOCATABLE, PUBLIC:: CREASIC ( : ) ! SIC code
        REAL   , ALLOCATABLE, PUBLIC:: EMREPREA( : )  ! replacement emissions
        REAL   , ALLOCATABLE, PUBLIC:: PRJFCREA( : )  ! projection factor
        REAL   , ALLOCATABLE, PUBLIC:: MKTPNREA( : )  ! market pen rt [frac/yr]
        CHARACTER(SCCLEN3), ALLOCATABLE, PUBLIC:: CSCCREA( : ) ! New SCC
        CHARACTER(SPNLEN3), ALLOCATABLE, PUBLIC:: CSPFREA( : ) ! SPROF
        CHARACTER(120), ALLOCATABLE, PUBLIC:: REACOMT( : ) ! packet comment

!.........  PROJECT PTS/AMS-packet-specific data tables
        CHARACTER(SICLEN3), ALLOCATABLE, PUBLIC:: CPRJSIC ( : ) ! SIC code
        REAL   , ALLOCATABLE, PUBLIC:: PRJFC   ( : )  ! projection factor
        CHARACTER(120), ALLOCATABLE, PUBLIC:: PRJCOMT( : ) ! packet comment

!.........  MACT-packet-specific data tables
        CHARACTER(STPLEN3), ALLOCATABLE, PUBLIC:: CMACSRCTYP( : )  ! source code type
        REAL   , ALLOCATABLE, PUBLIC:: MACEXEFF( : )  ! CE for existing sources (0-100)
        REAL   , ALLOCATABLE, PUBLIC:: MACNWEFF( : )  ! CE for new sources (0-100)
        REAL   , ALLOCATABLE, PUBLIC:: MACNWFRC( : )  ! fraction of new sources (0-100)
        CHARACTER(120), ALLOCATABLE, PUBLIC:: MACCOMT( : ) ! packet comment

!..............................................................................
!.........  Cntlmat work arrays for computing matrices
        REAL   , ALLOCATABLE :: BACKOUT ( : )   ! factor used to account for pol
                                                ! specific control info that is
                                                ! already in the inventory
        REAL   , ALLOCATABLE :: DATVAL  ( :,: ) ! emissions and control settings
        REAL   , ALLOCATABLE :: FACTOR  ( : )   ! multiplicative controls

!..............................................................................
!.........  Cntlmat variables for reporting ...
        INTEGER, PARAMETER :: NCNTLRPT = 5    ! no. of Cntlmat reports
        INTEGER            :: RPTDEV( NCNTLRPT ) = ( / 0, 0, 0, 0, 0 / )

        INTEGER, ALLOCATABLE :: GRPINDX ( : )   ! index from sources to groups
        INTEGER, ALLOCATABLE :: GRPSTIDX( : )   ! sorting index

        REAL   , ALLOCATABLE :: GRPINEM ( :,: ) ! initial emissions
        REAL   , ALLOCATABLE :: GRPOUTEM( :,: ) ! controlled emissions

        LOGICAL, ALLOCATABLE :: GRPFLAG ( : )   ! true: group controlled

        CHARACTER(STALEN3+SCCLEN3), ALLOCATABLE :: GRPCHAR( : ) ! group chars


!..............................................................................
!.........  CONTROL MATRICES...

!.........  Output pollutants (variable names) for control matrices
        INTEGER, PUBLIC :: NVCMULT = 0  ! number of multultiplicative variables
        INTEGER, PUBLIC :: NVPROJ  = 0  ! number of projection variables

        LOGICAL, PUBLIC :: POLSFLAG = .FALSE. ! true: proj data-spec assignment
        LOGICAL, ALLOCATABLE, PUBLIC :: PCTLFLAG( :,: ) ! control flags

        CHARACTER(IOVLEN3), ALLOCATABLE, PUBLIC :: PNAMMULT( : ) ! mult
        CHARACTER(IOVLEN3), ALLOCATABLE, PUBLIC :: PNAMPROJ( : ) ! projectn

!.........  Generic projection matrix (used in Smkreport, structure
!           plans for future 2-d structure)
        REAL   , ALLOCATABLE, PUBLIC :: PRMAT  ( :,: )

!.........  Multiplicative control matrices, dim n*src, mxpolpgp
        REAL   , ALLOCATABLE, PUBLIC :: ACUMATX( :,: ) ! area
        REAL   , ALLOCATABLE, PUBLIC :: MCUMATX( :,: ) ! mobile
        REAL   , ALLOCATABLE, PUBLIC :: PCUMATX( :,: ) ! point

!.........  Reactivity control matrices...
!.........  Indices for source IDs with reactivity controls
        INTEGER, ALLOCATABLE, PUBLIC :: ACRIDX( : ) ! area, dim: nasreac
        INTEGER, ALLOCATABLE, PUBLIC :: MCRIDX( : ) ! mobile, dim: nmsreac
        INTEGER, ALLOCATABLE, PUBLIC :: PCRIDX( : ) ! point, dim: npsreac

!.........  Reactivity-based base-year inventory emissions 
        REAL   , ALLOCATABLE, PUBLIC :: ACRREPEM( : ) ! area
        REAL   , ALLOCATABLE, PUBLIC :: MCRREPEM( : ) ! mobile
        REAL   , ALLOCATABLE, PUBLIC :: PCRREPEM( : ) ! point

!......... Base-year inventory emissions for reactivity sources
        REAL   , ALLOCATABLE, PUBLIC :: ACRBASEM( : ) ! area
        REAL   , ALLOCATABLE, PUBLIC :: MCRBASEM( : ) ! mobile
        REAL   , ALLOCATABLE, PUBLIC :: PCRBASEM( : ) ! point

!.........  Reactivity-based projection factors to future year
        REAL   , ALLOCATABLE, PUBLIC :: ACRPRJFC( : ) ! area
        REAL   , ALLOCATABLE, PUBLIC :: MCRPRJFC( : ) ! mobile
        REAL   , ALLOCATABLE, PUBLIC :: PCRPRJFC( : ) ! point

!.........  Market penetration to future year [fraction]
        REAL   , ALLOCATABLE, PUBLIC :: ACRMKTPN( : ) ! area
        REAL   , ALLOCATABLE, PUBLIC :: MCRMKTPN( : ) ! mobile
        REAL   , ALLOCATABLE, PUBLIC :: PCRMKTPN( : ) ! point

!.........  Reactivity-based source category codes (SCCs) 
        CHARACTER(SCCLEN3), ALLOCATABLE, PUBLIC :: ACRCSCC( : ) ! area
        CHARACTER(SCCLEN3), ALLOCATABLE, PUBLIC :: MCRCSCC( : ) ! mobile
        CHARACTER(SCCLEN3), ALLOCATABLE, PUBLIC :: PCRCSCC( : ) ! point

!.........  Reactivity-based speciation profile codes 
        CHARACTER(SPNLEN3), ALLOCATABLE, PUBLIC :: ACRSPROF( : ) ! area
        CHARACTER(SPNLEN3), ALLOCATABLE, PUBLIC :: MCRSPROF( : ) ! mobile
        CHARACTER(SPNLEN3), ALLOCATABLE, PUBLIC :: PCRSPROF( : ) ! point

!.........  Reactivity-based speciation factors (mole-based or mass-based )
        REAL   , ALLOCATABLE, PUBLIC :: ACRFAC( :,: ) ! area: ansreac, ansmatv
        REAL   , ALLOCATABLE, PUBLIC :: MCRFAC( :,: ) ! mobile: mnsreac, mnsmatv
        REAL   , ALLOCATABLE, PUBLIC :: PCRFAC( :,: ) ! point: pnsreac, pnsmatv
        REAL   , ALLOCATABLE, PUBLIC :: RMTXMASS( :,: ) ! general mass-based
        REAL   , ALLOCATABLE, PUBLIC :: RMTXMOLE( :,: ) ! general mole-based

        END MODULE MODCNTRL
