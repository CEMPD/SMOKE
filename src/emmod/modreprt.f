
        MODULE MODREPRT

!***********************************************************************
!  Module body starts at line
!
!  DESCRIPTION:
!     This module contains the public data that are used for Smkreport
!
!  PRECONDITIONS REQUIRED:
!
!  SUBROUTINES AND FUNCTIONS CALLED:
!
!  REVISION HISTORY:
!     Created 7/2000 by M. Houyoux
!
!***************************************************************************
!
! Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
!                System
! File: @(#)$Id$
!
! COPYRIGHT (C) 2000, MCNC--North Carolina Supercomputing Center
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

        INCLUDE 'EMPRVT3.EXT'

!.........  Packet-specific parameters
        INTEGER, PARAMETER, PUBLIC :: NALLPCKT = 7
        INTEGER, PARAMETER, PUBLIC :: RPKTLEN  = 21
        INTEGER, PARAMETER, PUBLIC :: RPT_IDX  = 1  ! pos. report pkt in master
        INTEGER, PARAMETER, PUBLIC :: TIM_IDX  = 2  ! pos. report time in mstr
        INTEGER, PARAMETER, PUBLIC :: O3S_IDX  = 3  ! pos. O3 season in master
        INTEGER, PARAMETER, PUBLIC :: FIL_IDX  = 4  ! pos. file pkt in master
        INTEGER, PARAMETER, PUBLIC :: DEL_IDX  = 5  ! pos. delimiter pkt in mstr
        INTEGER, PARAMETER, PUBLIC :: REG_IDX  = 6  ! pos. region pkt in master
        INTEGER, PARAMETER, PUBLIC :: SBG_IDX  = 7  ! pos. subgrid pkt in master
        INTEGER, PARAMETER, PUBLIC :: ELEVOUT3 = 1  ! code for elevated only
        INTEGER, PARAMETER, PUBLIC :: PINGOUT3 = 2  ! code for PinG only
        INTEGER, PARAMETER, PUBLIC :: NOELOUT3 = 3  ! code for low-level only
        INTEGER, PARAMETER, PUBLIC :: LENLAB3  = 200! length of group labels
        INTEGER, PARAMETER, PUBLIC :: LENELV3  = 1  ! length of elev status
        INTEGER, PARAMETER, PUBLIC :: LENTTL3  = 300! length of titles
        INTEGER, PARAMETER, PUBLIC :: QAFMTL3  = 2000 ! lngth of format statmnt

        CHARACTER(LEN=RPKTLEN), PARAMETER, PUBLIC :: 
     &                          ALLPCKTS( NALLPCKT ) = 
     &                                  ( / '/CREATE REPORT/      ',
     &                                      '/REPORT TIME/        ',
     &                                      '/O3SEASON/           ',
     &                                      '/NEWFILE/            ',
     &                                      '/DELIMITER/          ',
     &                                      '/DEFINE GROUP REGION/',
     &                                      '/DEFINE SUBGRID/     '  / )

!.........  Define types needed for module
        TYPE :: EACHRPT

            SEQUENCE

            INTEGER       :: BEGSUMHR      ! start hour for summing
            INTEGER       :: ELEVSTAT      ! elev/no-elev/all status
            INTEGER       :: OUTTIME       ! end hour for summing
            INTEGER       :: NUMDATA       ! number of data values
            INTEGER       :: NUMTITLE      ! number of titles
            INTEGER       :: RENDLIN       ! rpt packet end line
            INTEGER       :: RSTARTLIN     ! rpt packet start line
            INTEGER       :: SCCRES        ! SCC resolution

            LOGICAL       :: BYCELL        ! true: by cell
            LOGICAL       :: BYCNRY        ! true: by country code
            LOGICAL       :: BYCNTY        ! true: by county code
            LOGICAL       :: BYCONAM       ! true: by country name
            LOGICAL       :: BYCYNAM       ! true: by county name
            LOGICAL       :: BYDATE        ! true: by date
            LOGICAL       :: BYELEV        ! true: by elev status
            LOGICAL       :: BYHOUR        ! true: by hour
            LOGICAL       :: BYSCC         ! true: by SCC 
            LOGICAL       :: BYSRC         ! true: by source 
            LOGICAL       :: BYSTAT        ! true: by state code
            LOGICAL       :: BYSTNAM       ! true: by state name
            LOGICAL       :: BYRCL         ! true: by road class (mb)
            LOGICAL       :: LAYFRAC       ! true: use PLAY file
            LOGICAL       :: NORMCELL      ! true: output SCC name
            LOGICAL       :: O3SEASON      ! true: use O3 seas data
            LOGICAL       :: SCCNAM        ! true: output SCC name
            LOGICAL       :: SRCNAM        ! true: output facility nm
            LOGICAL       :: STKPARM       ! true: output stack parms
            LOGICAL       :: USEGMAT       ! true: use gridding
            LOGICAL       :: USEHOUR       ! true: use hourly data
            LOGICAL       :: USESLMAT      ! true: use mole spec
            LOGICAL       :: USESSMAT      ! true: use mass spec

            CHARACTER*1            :: DELIM         ! output delimeter
            CHARACTER*20           :: DATAFMT       ! data format
            CHARACTER(LEN=LENLAB3) :: REGNNAM       ! region names
            CHARACTER(LEN=LENLAB3) :: SUBGNAM       ! subgrid names
            CHARACTER*300          :: OFILENAM      ! output names

        END TYPE

!.........  Input file characteristics not available in MODINFO
        INTEGER, PUBLIC :: NCOLS  = 0       ! no. grid columns
        INTEGER, PUBLIC :: NGRID  = 0       ! no. grid cells
        INTEGER, PUBLIC :: NMAJOR = 0       ! no. major sources
        INTEGER, PUBLIC :: NMATX  = 0       ! size of gridding matrix
        INTEGER, PUBLIC :: NPING  = 0       ! no. PinG sources
        INTEGER, PUBLIC :: NROWS  = 0       ! no. grid rows
        INTEGER, PUBLIC :: NSTEPS = 0       ! no. time steps
        INTEGER, PUBLIC :: SDATE  = 0       ! Julian start date of run
        INTEGER, PUBLIC :: STIME  = 0       ! start time of run (HHMMSS)
        INTEGER, PUBLIC :: TSTEP  = 0       ! time step (HHMMSS)
        INTEGER, PUBLIC :: TZONE  = 0       ! time zone of hourly data (0-23)

!.........  Controls for whether input files are needed
!.........  These variables are set once, and never reset
        LOGICAL, PUBLIC :: CUFLAG = .FALSE. ! true: read in multipl. control matrix
        LOGICAL, PUBLIC :: CAFLAG = .FALSE. ! true: read in additive control matrix
        LOGICAL, PUBLIC :: CRFLAG = .FALSE. ! true: read in reactivity control matrix
        LOGICAL, PUBLIC :: GFLAG  = .FALSE. ! true: read in grd matrix and G_GRIDPATH
        LOGICAL, PUBLIC :: LFLAG  = .FALSE. ! true: read in layer fracs file
        LOGICAL, PUBLIC :: NFLAG  = .FALSE. ! true: read in SCC names file
        LOGICAL, PUBLIC :: SLFLAG = .FALSE. ! true: read in mole speciation matrix
        LOGICAL, PUBLIC :: SSFLAG = .FALSE. ! true: read in mass speciation matrix
        LOGICAL, PUBLIC :: TFLAG  = .FALSE. ! true: read in hourly emissions
        LOGICAL, PUBLIC :: VFLAG  = .FALSE. ! true: read in elevated source file PELV
        LOGICAL, PUBLIC :: YFLAG  = .FALSE. ! true: read in country, state, county info

!.........  REPCONFIG file characteristics
        INTEGER, PUBLIC :: MXGRPREC = 0   ! max no. of raw recs for any group
        INTEGER, PUBLIC :: MXINDAT  = 0   ! max no. data vars listed per rep
        INTEGER, PUBLIC :: MXOUTDAT = 0   ! max no. data vars for output
        INTEGER, PUBLIC :: MXTITLE  = 0   ! max number of titles per report
        INTEGER, PUBLIC :: NFILE    = 0   ! no. of output filess
        INTEGER, PUBLIC :: NLINE_RC = 0   ! no. of lines in the file
        INTEGER, PUBLIC :: NREPORT  = 0   ! no. of reports
        INTEGER, PUBLIC :: NREGRAW  = 0   ! no. raw file region groups
        INTEGER, PUBLIC :: NSBGRAW  = 0   ! no. raw file subgrids 
	INTEGER, PUBLIC :: MINC     = 0   ! minimum output source chars

        LOGICAL, PUBLIC :: RC_ERROR = .FALSE.  ! true: error found reading file

!.........  Group dimensions of group characteristic arrays
        INTEGER, PUBLIC :: MXREGREC = 0   ! max no. records in full region grp
        INTEGER, PUBLIC :: MXSUBREC = 0   ! max no. records in full subgrids
        INTEGER, PUBLIC :: NREGNGRP = 0   ! no. region groups
        INTEGER, PUBLIC :: NSUBGRID = 0   ! no. subgrids 

!.........  Group characteristics arrays
        INTEGER, ALLOCATABLE, PUBLIC :: VALIDCEL( :,: )   ! valid cell numbers
        INTEGER, ALLOCATABLE, PUBLIC :: VALIDRGN( :,: )   ! valid region numbers

!.........  Group label arrays
        CHARACTER(LEN=LENLAB3), ALLOCATABLE, PUBLIC :: REGNNAM( : ) ! region group names
        CHARACTER(LEN=LENLAB3), ALLOCATABLE, PUBLIC :: SUBGNAM( : ) ! subgrid names

!.........  Report characteristics arrays, dimenioned by NREPORT

        TYPE( EACHRPT ), ALLOCATABLE, PUBLIC :: ALLRPT( : ) ! integer recs

        LOGICAL        , ALLOCATABLE, PUBLIC :: ALLOUTHR( :,: ) ! true: write

        CHARACTER(LEN=LENTTL3), ALLOCATABLE, PUBLIC :: TITLES ( :,: ) ! report titles
        CHARACTER(LEN=IOVLEN3), ALLOCATABLE, PUBLIC :: INDNAM ( :,: ) ! var nams
        CHARACTER(LEN=IOVLEN3), ALLOCATABLE, PUBLIC :: OUTDNAM( :,: ) ! var nams
        CHARACTER(LEN=IOULEN3), ALLOCATABLE, PUBLIC :: ALLUNIT( :,: ) ! units

!.........  Temporary output file-specific settings
!.........  All widths include leading blanks and trailing commas
        INTEGER      , PUBLIC :: CELLWIDTH=0 ! width of cell output columns
        INTEGER      , PUBLIC :: CHARWIDTH=0 ! width of source char output cols
        INTEGER      , PUBLIC :: COWIDTH  =0 ! width of country name column
        INTEGER      , PUBLIC :: CYWIDTH  =0 ! width of county name column
        INTEGER      , PUBLIC :: DATEWIDTH=0 ! width of date column
        INTEGER      , PUBLIC :: ELEVWIDTH=0 ! width of elevated srcs flag col
        INTEGER      , PUBLIC :: HOURWIDTH=0 ! width of hour column
        INTEGER      , PUBLIC :: PDSCWIDTH=0 ! width of plant description col
        INTEGER      , PUBLIC :: REGNWIDTH=0 ! width of region column
        INTEGER      , PUBLIC :: SCCWIDTH =0 ! width of SCC
        INTEGER      , PUBLIC :: SDSCWIDTH=0 ! width of SCC description column
        INTEGER      , PUBLIC :: SRCWIDTH =0 ! width of source IDs column
        INTEGER      , PUBLIC :: STWIDTH  =0 ! width of state name column
        INTEGER      , PUBLIC :: STKPWIDTH=0 ! width of stack parameters columns

        CHARACTER*50 , PUBLIC :: CELLFMT     ! format string for cell columns
        CHARACTER*50 , PUBLIC :: CELLFMT     ! format string for cell columns
        CHARACTER*50 , PUBLIC :: HOURFMT     ! format string for hour column
        CHARACTER*50 , PUBLIC :: REGNFMT     ! format string for region column
        CHARACTER*50 , PUBLIC :: SRCFMT      ! format string for source IDs
        CHARACTER*100, PUBLIC :: STKPFMT     ! format string for stack params
        CHARACTER*200, PUBLIC :: CHARFMT     ! format string for source chars
        CHARACTER*300, PUBLIC :: FIL_ONAME   ! output file, physical or logical

!.........  Temporary packet-specific settings
        INTEGER, PUBLIC :: PKT_IDX  = 0       ! index to ALLPCKTS for current
        INTEGER, PUBLIC :: PKTEND   = 0       ! ending line number of packet
        INTEGER, PUBLIC :: PKTSTART = 0       ! starting line number of packet

        LOGICAL, PUBLIC :: INGROUP  = .FALSE. ! true: currently in a group defn
        LOGICAL, PUBLIC :: INPACKET = .FALSE. ! true: currently in a packet
        LOGICAL, PUBLIC :: INREPORT = .FALSE. ! true: currently in a report defn

        INTEGER, PUBLIC :: PKTCOUNT ( NALLPCKT ) ! no. of packets of given type
        LOGICAL, PUBLIC :: PKTSTATUS( NALLPCKT ) ! currently active packet

        CHARACTER(LEN=RPKTLEN), PUBLIC :: PCKTNAM

!.........  Temporary group-specific settings
        INTEGER, PUBLIC :: GRPNRECS = 0          ! no. records in a group

        LOGICAL, PUBLIC :: GRP_INCLSTAT = .TRUE. ! true=include; false=exclude

        CHARACTER(LEN=LENLAB3), PUBLIC :: GRP_LABEL = ' ' ! generic group label

!.........  Temporary report-specific settings
        TYPE( EACHRPT ), PUBLIC :: RPT_

        LOGICAL, PUBLIC :: LSUBGRID      ! true: select with a subgrid
        LOGICAL, PUBLIC :: LREGION       ! true: select with a region group

        CHARACTER(LEN=LENTTL3), PUBLIC :: TITLE    ! current line title

        ! Input units for each output data column
        CHARACTER(LEN=IOULEN3), ALLOCATABLE, PUBLIC :: INUNIT( : ) 

        ! Conversion factors for each output data column
        REAL, ALLOCATABLE, PUBLIC :: UCNVFAC( : )

!.........  Temporary line-specific settings
        LOGICAL, PUBLIC :: LIN_SUBDATA      ! true: line is SELECT DATA
        LOGICAL, PUBLIC :: LIN_TITLE        ! true: line is TITLE
        LOGICAL, PUBLIC :: LIN_UNIT         ! true: line is UNITS

        END MODULE MODREPRT
