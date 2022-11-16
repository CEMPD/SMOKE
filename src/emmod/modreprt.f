
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
!     Revised 7/2003 by A. Holland
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

!.........  Packet-specific parameters
        INTEGER, PARAMETER, PUBLIC :: NALLPCKT = 11
        INTEGER, PARAMETER, PUBLIC :: RPKTLEN  = 21
        INTEGER, PARAMETER, PUBLIC :: NCRTSYBL = 13
        INTEGER, PARAMETER, PUBLIC :: RPT_IDX  = 1  ! pos. report pkt in master
        INTEGER, PARAMETER, PUBLIC :: TIM_IDX  = 2  ! pos. report time in mstr
        INTEGER, PARAMETER, PUBLIC :: ADY_IDX  = 3  ! pos. average day in master
        INTEGER, PARAMETER, PUBLIC :: FIL_IDX  = 4  ! pos. file pkt in master
        INTEGER, PARAMETER, PUBLIC :: DEL_IDX  = 5  ! pos. delimiter pkt in mstr
        INTEGER, PARAMETER, PUBLIC :: REG_IDX  = 6  ! pos. region pkt in master
        INTEGER, PARAMETER, PUBLIC :: SBG_IDX  = 7  ! pos. subgrid pkt in master
        INTEGER, PARAMETER, PUBLIC :: ELG_IDX  = 8  ! pos. elevated grp in mstr
        INTEGER, PARAMETER, PUBLIC :: PNG_IDX  = 9  ! pos. ping in master
        INTEGER, PARAMETER, PUBLIC :: ELV_IDX  = 10 ! pos. elevated in master
        INTEGER, PARAMETER, PUBLIC :: LAB_IDX  = 11 ! pos. label in master
        INTEGER, PARAMETER, PUBLIC :: NOELOUT3 = 1  ! code for low-level only
        INTEGER, PARAMETER, PUBLIC :: ELEVOUT3 = 2  ! code for elevated only
        INTEGER, PARAMETER, PUBLIC :: PINGOUT3 = 3  ! code for PinG only
        INTEGER, PARAMETER, PUBLIC :: LENLAB3  = 200! length of group labels
        INTEGER, PARAMETER, PUBLIC :: LENELV3  = 1  ! length of elev status
        INTEGER, PARAMETER, PUBLIC :: LENTTL3  = 300! length of titles
        INTEGER, PARAMETER, PUBLIC :: QAFMTL3  = 3500 ! lngth of format statmnt

        CHARACTER(4), PARAMETER, PUBLIC :: CHKPFX = 'chk_'

        CHARACTER(RPKTLEN), PARAMETER, PUBLIC :: 
     &                          ALLPCKTS( NALLPCKT ) = 
     &                                  ( / '/CREATE REPORT/      ',
     &                                      '/REPORT TIME/        ',
     &                                      '/AVEDAY/             ',
     &                                      '/NEWFILE/            ',
     &                                      '/DELIMITER/          ',
     &                                      '/DEFINE GROUP REGION/',
     &                                      '/DEFINE SUBGRID/     ',
     &                                      '/SPECIFY ELEV GROUPS/',
     &                                      '/SPECIFY PING/       ',
     &                                      '/SPECIFY ELEV/       ',
     &                                      '/SET LABEL/          ' / )
        
        CHARACTER(4), PARAMETER, PUBLIC :: CRTSYBL( NCRTSYBL ) = 
     &                                  ( / '=   ',
     &                                      '==  ',
     &                                      '=>  ',
     &                                      '>=  ',
     &                                      '=<  ',
     &                                      '<=  ',
     &                                      '<   ',
     &                                      '>   ',
     &                                      '+/- ',
     &                                      '-/+ ',
     &                                      'TOP ',
     &                                      'IS  ',
     &                                      '%   '  / )

!.........  Define types needed for module
        TYPE :: EACHRPT
            SEQUENCE
            INTEGER       :: BEGSUMHR      ! start hour for summing
            INTEGER       :: ELEVSTAT      ! elev/no-elev/all status
            INTEGER       :: OUTTIME       ! end hour for summing
            INTEGER       :: NUMDATA       ! number of data values
            INTEGER       :: NUMFILES      ! number of files per report
            INTEGER       :: NUMSECT       ! number of sections per report
            INTEGER       :: NUMTITLE      ! number of titles
            INTEGER       :: RENDLIN       ! rpt packet end line
            INTEGER       :: RPTMODE       ! rpt mode
            INTEGER       :: RPTNVAR       ! no. of variables per rpt or section
            INTEGER       :: RSTARTLIN     ! rpt packet start line
            INTEGER       :: SCCRES        ! SCC resolution
            INTEGER       :: SRGRES        ! surrogate resolution (1st or 2nd)
            LOGICAL       :: BYCELL        ! true: by cell
            LOGICAL       :: BYCNRY        ! true: by country code (geocode 2)
            LOGICAL       :: BYCNTY        ! true: by county code (geocode 4)
            LOGICAL       :: BYCONAM       ! true: by country name
            LOGICAL       :: BYCYNAM       ! true: by county name
            LOGICAL       :: BYDATE        ! true: by date
            LOGICAL       :: BYELEV        ! true: by elev status
            LOGICAL       :: BYERPTYP      ! true: by emissions release point type
            LOGICAL       :: ELVSTKGRP     ! true: stack gourp ID by elev status
            LOGICAL       :: BYGEO1        ! true: by geocode 1
            LOGICAL       :: BYGEO1NAM     ! true: by geocode 1 name
            LOGICAL       :: BYHOUR        ! true: by hour
            LOGICAL       :: BYLAYER       ! true: by layer
            LOGICAL       :: BYLATLON      ! true: by lat-lon coordinates
            LOGICAL       :: BYMON         ! true: by monthly temporal code
            LOGICAL       :: BYWEK         ! true: by weekly temporal code
            LOGICAL       :: BYDOM         ! true: by day of month temporal code
            LOGICAL       :: BYMND         ! true: by Monday diurnal temporal code
            LOGICAL       :: BYTUE         ! true: by Tuesday diurnal temporal code
            LOGICAL       :: BYWED         ! true: by Wednesday diurnal temporal code
            LOGICAL       :: BYTHU         ! true: by Thursday diurnal temporal code
            LOGICAL       :: BYFRI         ! true: by Friday diurnal temporal code
            LOGICAL       :: BYSAT         ! true: by Saturday diurnal temporal code
            LOGICAL       :: BYSUN         ! true: by Sunday diurnal temporal code
            LOGICAL       :: BYMET         ! true: by Met-based hourly temporal code
            LOGICAL       :: BYPLANT       ! true: by plant 
            LOGICAL       :: BYFACILITY    ! true: by Facility
            LOGICAL       :: BYUNIT        ! true: by Unit ID
            LOGICAL       :: BYSCC         ! true: by SCC 
            LOGICAL       :: BYSIC         ! true: by SIC 
            LOGICAL       :: BYINTGR       ! true: by INTEGRATE 
            LOGICAL       :: BYMACT        ! true: by MACT
            LOGICAL       :: BYNAICS       ! true: by NAICS
            LOGICAL       :: BYORIS        ! true: by ORIS 
            LOGICAL       :: BYBOILER      ! true: by boiler
            LOGICAL       :: BYSRCTYP      ! true: by source type
            LOGICAL       :: BYSPC         ! true: by speciation codes 
            LOGICAL       :: BYSRC         ! true: by source 
            LOGICAL       :: BYSTACK       ! true: by stack
            LOGICAL       :: BYSTKPARM     ! true: by stack and fugutive params
            LOGICAL       :: BYSTAT        ! true: by state code (geocode 3)
            LOGICAL       :: BYSTNAM       ! true: by state name
            LOGICAL       :: BYSRG         ! true: by surrogate codes
            LOGICAL       :: BYRCL         ! true: by road class (mb)
            LOGICAL       :: CARB          ! true: output CARB summary QA report
            LOGICAL       :: CHKPROJ       ! true: check projctns vs. rpt
            LOGICAL       :: CHKCNTL       ! true: check controls vs. rpt
            LOGICAL       :: LATLON        ! true: output stack coordinates
            LOGICAL       :: GRDCOR        ! true: output grid coordinates 
            LOGICAL       :: GRDPNT        ! true: output grid corner coordinates
            LOGICAL       :: LAYFRAC       ! true: use PLAY file
            LOGICAL       :: NORMCELL      ! true: normalize by cell area
            LOGICAL       :: NORMPOP       ! true: normalize by county pop
            LOGICAL       :: AVEDAY        ! true: use average day data
            LOGICAL       :: SCCNAM        ! true: output SCC name
            LOGICAL       :: SICNAM        ! true: output SIC name
            LOGICAL       :: GSPRONAM      ! true: output GSPRO name
            LOGICAL       :: MACTNAM       ! true: output MACT name
            LOGICAL       :: NAICSNAM      ! true: output NAICS name
            LOGICAL       :: ORISNAM       ! true: output ORIS name  
            LOGICAL       :: SRCNAM        ! true: output facility nm
            LOGICAL       :: STKPARM       ! true: output stack parms
            LOGICAL       :: FUGPARM       ! true: output fugitive parms
            LOGICAL       :: USEASCELEV    ! true: use ascii elevation file
            LOGICAL       :: USECRMAT      ! true: use reactivity controls
            LOGICAL       :: USECUMAT      ! true: use multiplicative controls
            LOGICAL       :: USEGMAT       ! true: use gridding
            LOGICAL       :: USEHOUR       ! true: use hourly data
            LOGICAL       :: USELABEL      ! true: use user-defined label
            LOGICAL       :: USEPRMAT      ! true: use projection matrix
            LOGICAL       :: USESLMAT      ! true: use mole spec
            LOGICAL       :: USESSMAT      ! true: use mass spec
            LOGICAL       :: SRCMAP        ! true: output src mapping

            CHARACTER          :: DELIM         ! output delimeter
            CHARACTER(20)      :: DATAFMT       ! data format
            CHARACTER(LENLAB3) :: LABEL         ! user-defined label
            CHARACTER(LENLAB3) :: REGNNAM       ! region names
            CHARACTER(IOVLEN3) :: SPCPOL        ! pollutant for BYSPC
            CHARACTER(LENLAB3) :: SUBGNAM       ! subgrid names
            CHARACTER(300)     :: OFILENAM      ! output names
        END TYPE

!.........  Input file characteristics not available in MODINFO
        INTEGER, PUBLIC :: EMLAYS = 1       ! no. emissions layers
        INTEGER, PUBLIC :: NMAJOR = 0       ! no. major sources
        INTEGER, PUBLIC :: NMATX  = 1       ! size of gridding matrix
        INTEGER, PUBLIC :: NPING  = 0       ! no. PinG sources
        INTEGER, PUBLIC :: NSTEPS = 0       ! no. time steps in data file
        INTEGER, PUBLIC :: PYEAR  = 0       ! projected inventory year
        INTEGER, PUBLIC :: PRBYR  = 0       ! proj matrix base yr
        INTEGER, PUBLIC :: PRPYR  = 0       ! proj matrix projected yr
        INTEGER, PUBLIC :: SDATE  = 0       ! Julian start date of run
        INTEGER, PUBLIC :: STIME  = 0       ! start time of run (HHMMSS)
        INTEGER, PUBLIC :: TSTEP  = 0       ! time step (HHMMSS)
        INTEGER, PUBLIC :: TZONE  = 0       ! time zone of hourly data (0-23)

        INTEGER, PUBLIC, ALLOCATABLE :: STKX( : )   ! x cell no. of stack
        INTEGER, PUBLIC, ALLOCATABLE :: STKY( : )   ! y cell no. of stack

!.........  Controls for whether input files are needed
!.........  These variables are set once, and never reset
        LOGICAL, PUBLIC :: AFLAG  = .FALSE. ! true: read in ASCII elevated file
        LOGICAL, PUBLIC :: CUFLAG = .FALSE. ! true: read in multipl. control matrix
        LOGICAL, PUBLIC :: CURPTFLG = .FALSE. ! true: read mult. cntl report
        LOGICAL, PUBLIC :: CRFLAG = .FALSE. ! true: read in reactivity control matrix
        LOGICAL, PUBLIC :: GFLAG  = .FALSE. ! true: read in grd matrix and G_GRIDPATH
        LOGICAL, PUBLIC :: GSFLAG = .FALSE. ! true: read gridding supplementary file
        LOGICAL, PUBLIC :: LFLAG  = .FALSE. ! true: read in layer fracs file
        LOGICAL, PUBLIC :: NFLAG  = .FALSE. ! true: read in SCC names file
        LOGICAL, PUBLIC :: SDFLAG = .FALSE. ! true: read in GSPRO names file
        LOGICAL, PUBLIC :: NIFLAG = .FALSE. ! true: read in SIC names file
        LOGICAL, PUBLIC :: NMFLAG = .FALSE. ! true: read in MACT names file
        LOGICAL, PUBLIC :: NNFLAG = .FALSE. ! true: read in NAICS names file
        LOGICAL, PUBLIC :: NOFLAG = .FALSE. ! true: read in ORIS names file
        LOGICAL, PUBLIC :: SLFLAG = .FALSE. ! true: read in mole speciation matrix
        LOGICAL, PUBLIC :: SSFLAG = .FALSE. ! true: read in mass speciation matrix
        LOGICAL, PUBLIC :: PRFLAG = .FALSE. ! true: read projection matrix
        LOGICAL, PUBLIC :: PRRPTFLG = .FALSE. ! true: read projectn report
        LOGICAL, PUBLIC :: PSFLAG = .FALSE. ! true: read spec supplementary
        LOGICAL, PUBLIC :: TFLAG  = .FALSE. ! true: read in hourly emissions
        LOGICAL, PUBLIC :: TSFLAG = .FALSE. ! true: read tmprl supplementary
        LOGICAL, PUBLIC :: VFLAG  = .FALSE. ! true: read in elevated source file PELV
        LOGICAL, PUBLIC :: YFLAG  = .FALSE. ! true: read in country, state, county info
        LOGICAL, PUBLIC :: DLFLAG = .FALSE. ! true: writing daily emission by layer
        LOGICAL, PUBLIC :: NFDFLAG= .FALSE. ! true: read in NFDRSCODE header
        LOGICAL, PUBLIC :: MATFLAG= .FALSE. ! true: read in MATBURNED header
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
        INTEGER, PUBLIC :: NSPCPOL  = 0   ! no. pollutants specified for BYSPC
        INTEGER, PUBLIC :: MINC     = 0   ! minimum output source chars
        INTEGER, PUBLIC :: MXRPTNVAR= 30  ! max. number of variables per report

        LOGICAL, PUBLIC :: RC_ERROR = .FALSE.  ! true: error found reading file
        LOGICAL, PUBLIC :: DATAMISS = .FALSE.  ! true: no SELECT DATA instrs
        LOGICAL, PUBLIC :: POFLAG   = .FALSE.  ! one or more rpts use pop data

!.........  Group dimensions of group characteristic arrays
        INTEGER, PUBLIC :: MXREGREC = 0   ! max no. records in full region grp
        INTEGER, PUBLIC :: MXSUBREC = 0   ! max no. records in full subgrids
        INTEGER, PUBLIC :: NREGNGRP = 0   ! no. region groups
        INTEGER, PUBLIC :: NSUBGRID = 0   ! no. subgrids
        
!.........  Group characteristics arrays
        INTEGER, ALLOCATABLE, PUBLIC :: NREGREC ( : )     ! no. recs per region grp
c        INTEGER, ALLOCATABLE, PUBLIC :: NSUBREC ( : )     ! no. recs per subgrid
        INTEGER, ALLOCATABLE, PUBLIC :: VALIDCEL( :,: )   ! valid cell numbers
        CHARACTER(FIPLEN3), ALLOCATABLE, PUBLIC :: EXCLDRGN( :,: )   ! excluded region codes

!.........  Group label arrays
        CHARACTER(LENLAB3), ALLOCATABLE, PUBLIC :: REGNNAM( : ) ! region group names
        CHARACTER(LENLAB3), ALLOCATABLE, PUBLIC :: SUBGNAM( : ) ! subgrid names

!.........  Allocatable arrays available across all reports
        INTEGER, ALLOCATABLE, PUBLIC :: LOC_BEGP( : )  ! actual src char string starts
        INTEGER, ALLOCATABLE, PUBLIC :: LOC_ENDP( : )  ! actual src char string ends

        LOGICAL, ALLOCATABLE, PUBLIC :: LSPCPOL ( : )  ! true: spc pol in *SSUP file
        CHARACTER(IOVLEN3), ALLOCATABLE, PUBLIC :: SPCPOL( : ) ! pols for BYSPC
        CHARACTER(IOVLEN3), ALLOCATABLE, PUBLIC :: ASCNAM( : ) ! pols from ASCII elevated file

!.........  Report characteristics arrays, dimensioned by NREPORT

        TYPE( EACHRPT ), ALLOCATABLE, PUBLIC :: ALLRPT( : )     ! integer recs

        LOGICAL        , ALLOCATABLE, PUBLIC :: ALLOUTHR( :,: ) ! true: write

        CHARACTER(LENTTL3), ALLOCATABLE, PUBLIC :: TITLES ( :,: ) ! report titles
        CHARACTER(IOVLEN3), ALLOCATABLE, PUBLIC :: INDNAM ( :,: ) ! var nams
        CHARACTER(IOVLEN3), ALLOCATABLE, PUBLIC :: OUTDNAM( :,: ) ! var nams
        CHARACTER(IOVLEN3), ALLOCATABLE, PUBLIC :: RDNAMES( :,: ) ! for reads
        CHARACTER(IODLEN3), ALLOCATABLE, PUBLIC :: ALLUSET( :,: ) ! units

!.........  Temporary output file-specific settings
!.........  All widths include leading blanks and trailing commas
        INTEGER      , PUBLIC :: CELLWIDTH =0 ! width of cell output columns
        INTEGER      , PUBLIC :: CHARWIDTH =0 ! width of source char output cols
        INTEGER      , PUBLIC :: COWIDTH   =0 ! width of country name column
        INTEGER      , PUBLIC :: CYWIDTH   =0 ! width of county name column
        INTEGER      , PUBLIC :: DATEWIDTH =0 ! width of date column
        INTEGER      , PUBLIC :: ELEVWIDTH =0 ! width of elevated srcs flag col
        INTEGER      , PUBLIC :: ERTYPWIDTH=0 ! width of emissions release point type col
        INTEGER      , PUBLIC :: GEO1WIDTH =0 ! width of geo level 1 name column
        INTEGER      , PUBLIC :: HOURWIDTH =0 ! width of hour column
        INTEGER      , PUBLIC :: LTLNWIDTH =0 ! width of lat/lon columns
        INTEGER      , PUBLIC :: LAMBWIDTH =0 ! width of lambert coord columns
        INTEGER      , PUBLIC :: LLGRDWIDTH=0 ! width of lat/lon grid coords columns
        INTEGER      , PUBLIC :: LABELWIDTH=0 ! width of user-defined label
        INTEGER      , PUBLIC :: LAYRWIDTH =0 ! width of layer number label
        INTEGER      , PUBLIC :: PDSCWIDTH =0 ! width of plant description col
        INTEGER      , PUBLIC :: REGNWIDTH =0 ! width of region column
        INTEGER      , PUBLIC :: SCCWIDTH  =0 ! width of SCC
        INTEGER      , PUBLIC :: SDSCWIDTH =0 ! width of SCC description column
        INTEGER      , PUBLIC :: SICWIDTH  =0 ! width of SIC
        INTEGER      , PUBLIC :: SIDSWIDTH =0 ! width of SIC description column
        INTEGER      , PUBLIC :: INTGRWIDTH=0 ! width of INT_STAT (INTEGRATE) 
        INTEGER      , PUBLIC :: MACTWIDTH =0 ! width of MACT
        INTEGER      , PUBLIC :: MACDSWIDTH=0 ! width of MACT description column
        INTEGER      , PUBLIC :: NAIWIDTH  =0 ! width of NAICS
        INTEGER      , PUBLIC :: NAIDSWIDTH=0 ! width of NAICS description column
        INTEGER      , PUBLIC :: ORSWIDTH  =0 ! width of ORIS
        INTEGER      , PUBLIC :: ORSDSWIDTH=0 ! width of ORIS description column
        INTEGER      , PUBLIC :: STYPWIDTH =0 ! width of source type code
        INTEGER      , PUBLIC :: STKGWIDTH =0 ! width of stack group ID code
        INTEGER      , PUBLIC :: SPCWIDTH  =0 ! width of speciation profile label
        INTEGER      , PUBLIC :: SPDSWIDTH =0 ! width of speciation profile description
        INTEGER      , PUBLIC :: SRCWIDTH  =0 ! width of source IDs column
        INTEGER      , PUBLIC :: SRG1WIDTH =0 ! width of primary surg column
        INTEGER      , PUBLIC :: SRG2WIDTH =0 ! width of fallback surg column
        INTEGER      , PUBLIC :: STWIDTH   =0 ! width of state name column
        INTEGER      , PUBLIC :: STKPWIDTH =0 ! width of stack parameters columns
        INTEGER      , PUBLIC :: FUGPWIDTH =0 ! width of fugitive parameters columns
        INTEGER      , PUBLIC :: UNITWIDTH =0 ! width of unit column
        INTEGER      , PUBLIC :: VARWIDTH  =0 ! width of variable column
        INTEGER      , PUBLIC :: MONWIDTH  =0 ! width of monthly profile label
        INTEGER      , PUBLIC :: WEKWIDTH  =0 ! width of weekly profile label
        INTEGER      , PUBLIC :: DOMWIDTH  =0 ! width of day of month profile label
        INTEGER      , PUBLIC :: MNDWIDTH  =0 ! width of monday profile label
        INTEGER      , PUBLIC :: TUEWIDTH  =0 ! width of tuesday profile label
        INTEGER      , PUBLIC :: WEDWIDTH  =0 ! width of wednesday profile label
        INTEGER      , PUBLIC :: THUWIDTH  =0 ! width of thursday profile label
        INTEGER      , PUBLIC :: FRIWIDTH  =0 ! width of friday profile label
        INTEGER      , PUBLIC :: SATWIDTH  =0 ! width of saturday profile label
        INTEGER      , PUBLIC :: SUNWIDTH  =0 ! width of sunday profile label
        INTEGER      , PUBLIC :: METWIDTH  =0 ! width of met-based profile label
        INTEGER      , PUBLIC :: UNITIDWIDTH =0 ! width of unit ID column

        CHARACTER(50),  PUBLIC :: CELLFMT     ! format string for cell columns
        CHARACTER(50),  PUBLIC :: DATEFMT     ! format string for date column
        CHARACTER(50),  PUBLIC :: HOURFMT     ! format string for hour column
        CHARACTER(50),  PUBLIC :: LTLNFMT     ! format string for lat/lons
        CHARACTER(70),  PUBLIC :: LAMBFMT     ! format string for lambert coord 
        CHARACTER(120), PUBLIC :: LLGRDFMT    ! format string for lat/lon grid coords
        CHARACTER(50),  PUBLIC :: LAYRFMT     ! format string for layer column
        CHARACTER(50),  PUBLIC :: REGNFMT     ! format string for region column
        CHARACTER(50),  PUBLIC :: STKGFMT     ! format string for stack group IDs 
        CHARACTER(50),  PUBLIC :: SRCFMT      ! format string for source IDs
        CHARACTER(50),  PUBLIC :: SRG1FMT     ! format string for primary surg
        CHARACTER(50),  PUBLIC :: SRG2FMT     ! format string for fallback surg
        CHARACTER(100), PUBLIC :: STKPFMT     ! format string for stack params
        CHARACTER(100), PUBLIC :: FUGPFMT     ! format string for fugitive params
        CHARACTER(200), PUBLIC :: CHARFMT     ! format string for source chars
        CHARACTER(300), PUBLIC :: FIL_ONAME   ! output file, physical or logical

!.........  Temporary packet-specific settings
        INTEGER, PUBLIC :: PKT_IDX  = 0       ! index to ALLPCKTS for current
        INTEGER, PUBLIC :: PKTEND   = 0       ! ending line number of packet
        INTEGER, PUBLIC :: PKTSTART = 0       ! starting line number of packet

        LOGICAL, PUBLIC :: INGROUP  = .FALSE. ! true: currently in a group defn
        LOGICAL, PUBLIC :: INPACKET = .FALSE. ! true: currently in a packet
        LOGICAL, PUBLIC :: INREPORT = .FALSE. ! true: currently in a report defn
        LOGICAL, PUBLIC :: INSPCIFY = .FALSE. ! true: currently in a group defn

        INTEGER, PUBLIC :: PKTCOUNT ( NALLPCKT ) ! no. of packets of given type
        LOGICAL, PUBLIC :: PKTSTATUS( NALLPCKT ) ! currently active packet

        CHARACTER(RPKTLEN), PUBLIC :: PCKTNAM

!.........  Temporary group-specific settings
        INTEGER, PUBLIC :: GRPNRECS = 0          ! no. records in a group

        LOGICAL, PUBLIC :: GRP_INCLSTAT = .TRUE. ! true=include; false=exclude

        CHARACTER(LENLAB3), PUBLIC :: GRP_LABEL = ' ' ! generic group label

!.........  Temporary specify-specific settings
        INTEGER, PUBLIC :: SPCF_NAND = 0
        INTEGER, PUBLIC :: SPCF_NOR  = 0

!.........  Temporary report-specific settings
        TYPE( EACHRPT ), PUBLIC :: RPT_

        INTEGER, PUBLIC :: ASCDATA       ! no. of data from ASCII elevated file
        INTEGER, PUBLIC :: ASCREC        ! line no. of ASCII elevated file
        INTEGER, PUBLIC :: EDATE         ! Julian ending date
        INTEGER, PUBLIC :: ETIME         ! ending time (HHMMSS)
        INTEGER, PUBLIC :: RPTNSTEP      ! no. of time steps for current report

        LOGICAL, PUBLIC :: LSUBGRID      ! true: select with a subgrid
        LOGICAL, PUBLIC :: LREGION       ! true: select with a region group

        CHARACTER(IODLEN3), PUBLIC :: UNITSET  ! current line units
        CHARACTER(LENTTL3), PUBLIC :: TITLE    ! current line title

        ! Output units for each output data column
        CHARACTER(IODLEN3), ALLOCATABLE, PUBLIC :: OUTUNIT( : ) 

        ! Conversion factors for each output data column
        REAL, ALLOCATABLE, PUBLIC :: UCNVFAC( : )

!.........  Temporary line-specific settings
        LOGICAL, PUBLIC :: LIN_DEFGRP       ! true: line is DEFINE GROUP
        LOGICAL, PUBLIC :: LIN_GROUP        ! true: line is group entry
        LOGICAL, PUBLIC :: LIN_SPCIFY       ! true: line is specification entry
        LOGICAL, PUBLIC :: LIN_SUBDATA      ! true: line is SELECT DATA
        LOGICAL, PUBLIC :: LIN_SUBGRID      ! true: line is SELECT SUBGRID
        LOGICAL, PUBLIC :: LIN_SUBREGN      ! true: line is SELECT REGION
        LOGICAL, PUBLIC :: LIN_TITLE        ! true: line is TITLE
        LOGICAL, PUBLIC :: LIN_UNIT         ! true: line is UNITS

        END MODULE MODREPRT
