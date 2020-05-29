
        SUBROUTINE WRREPHDR( FDEV, RCNT, FILENUM, LH, OUTFMT )

C***********************************************************************
C  subroutine body starts at line 
C
C  DESCRIPTION:
C     This subroutine writes the header lines for a report based on the 
C     settings of the report-specific flags.  The header lines include lines 
C     for identifying the type of report as well as the column headers, 
C     which are set based on the labels for the generic output data columns.
C     It also determines the column widths
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C     Created 7/2000 by M Houyoux
C     Revised 7/2003 by A. Holland
C
C***********************************************************************
C  
C Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
C                System
C File: @(#)$Id$
C  
C COPYRIGHT (C) 2004, Environmental Modeling for Policy Development
C All Rights Reserved
C 
C Carolina Environmental Program
C University of North Carolina at Chapel Hill
C 137 E. Franklin St., CB# 6116
C Chapel Hill, NC 27599-6116
C 
C smoke@unc.edu
C  
C Pathname: $Source$
C Last updated: $Date$ 
C  
C***********************************************************************

C.........  MODULES for public variables
C...........   This module is the inventory arrays
        USE MODSOURC, ONLY: STKHT, STKDM, STKTK, STKVE, CPDESC, FUGHGT,
     &                      FUGWID, FUGLEN, FUGANG

C.........  This module contains the lists of unique source characteristics
        USE MODLISTS, ONLY: NINVSCC, SCCDESC, SCCDLEV, NINVSIC, SICDESC,
     &                      NINVMACT, MACTDESC, NINVNAICS, NAICSDESC

C.........  This module contains Smkreport-specific settings
        USE MODREPRT, ONLY: QAFMTL3, AFLAG, OUTDNAM, RPT_, LREGION,
     &                      PDSCWIDTH, VARWIDTH, DATEFMT, DATEWIDTH,
     &                      HOURFMT, HOURWIDTH, LAYRFMT, LAYRWIDTH,
     &                      CELLFMT, CELLWIDTH, SRCFMT, SRCWIDTH,
     &                      REGNFMT, REGNWIDTH, COWIDTH, STWIDTH,
     &                      CYWIDTH, SCCWIDTH, SRG1FMT, SRG1WIDTH,
     &                      SRG2FMT, SRG2WIDTH, MONWIDTH, WEKWIDTH,
     &                      DOMWIDTH, MNDWIDTH, TUEWIDTH, WEDWIDTH,
     &                      THUWIDTH, FRIWIDTH, SATWIDTH, SUNWIDTH,
     &                      METWIDTH,
     &                      CHARFMT, CHARWIDTH, STKPFMT, STKPWIDTH,
     &                      SPCWIDTH, ELEVWIDTH, SDSCWIDTH, UNITWIDTH,
     &                      MINC, LENELV3, SDATE, STIME, EDATE, ETIME,
     &                      PYEAR, PRBYR, PRPYR, OUTUNIT, TITLES,
     &                      ALLRPT, LOC_BEGP, LOC_ENDP,
     &                      SICWIDTH, SIDSWIDTH, MACTWIDTH, MACDSWIDTH,
     &                      NAIWIDTH, NAIDSWIDTH, STYPWIDTH, SPDSWIDTH,
     &                      LTLNFMT, LTLNWIDTH, LABELWIDTH, DLFLAG,
     &                      NFDFLAG, MATFLAG, ORSWIDTH, ORSDSWIDTH,
     &                      STKGWIDTH, STKGFMT, INTGRWIDTH, GEO1WIDTH,
     &                      ERTYPWIDTH, FUGPFMT, FUGPWIDTH, LAMBWIDTH,
     &                      LAMBFMT, LLGRDFMT, LLGRDWIDTH, POINTWIDTH

C.........  This module contains report arrays for each output bin
        USE MODREPBN, ONLY: NOUTBINS, BINX, BINY, BINSMKID, BINREGN,
     &                      BINSRGID1, BINSRGID2, BINMONID, BINWEKID,
     &                      BINDOMID, BINMNDID, BINTUEID, BINWEDID,
     &                      BINTHUID, BINFRIID, BINSATID, BINSUNID,
     &                      BINMETID, BINRCL, BINDATA, BINSNMIDX,
     &                      BINCYIDX, BINSTIDX, BINCOIDX, BINSPCID,
     &                      BINPLANT, BINSIC, BINSICIDX, BINMACT, 
     &                      BINMACIDX, BINNAICS, BINNAIIDX, BINSRCTYP,
     &                      BINORIS, BINORSIDX, BINSTKGRP, BININTGR,
     &                      BINGEO1IDX, BINUNIT, BINSPCIDX

C.........  This module contains the arrays for state and county summaries
        USE MODSTCY, ONLY: NCOUNTRY, NSTATE, NCOUNTY, STCYPOPYR,
     &                     CTRYNAM, STATNAM, CNTYNAM, ORISDSC, NORIS,
     &                     NGEOLEV1, GEOLEV1NAM

C.........  This module contains the global variables for the 3-d grid
        USE MODGRID, ONLY: GRDNM

CC...........  This module contains the information about the source category
        USE MODINFO, ONLY: NCHARS, CATEGORY, CATDESC, BYEAR, INVPIDX,
     &                     EANAM, ATTRUNIT

C.........  This module contsin the speciation variables
        USE MODSPRO,  ONLY : NSPROF, SPROFN, SPCDESC


        IMPLICIT NONE

C...........   INCLUDES
        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters

C...........  EXTERNAL FUNCTIONS and their descriptions:
        CHARACTER(2)    CRLF
        INTEGER         STR2INT
        CHARACTER(14)   MMDDYY
        INTEGER         WKDAY
        LOGICAL         USEEXPGEO

        EXTERNAL   CRLF, STR2INT, MMDDYY, WKDAY, USEEXPGEO

C...........   SUBROUTINE ARGUMENTS
        INTEGER     , INTENT (IN) :: FDEV       ! output file unit number
        INTEGER     , INTENT (IN) :: RCNT       ! report count
        INTEGER     , INTENT (IN) :: FILENUM    ! file number
        INTEGER     , INTENT(OUT) :: LH         ! header width
        CHARACTER(QAFMTL3),
     &                INTENT(OUT) :: OUTFMT     ! output record format

C...........   Local parameters
        INTEGER, PARAMETER :: OLINELEN = 10000
        INTEGER, PARAMETER :: IHDRDATE = 1
        INTEGER, PARAMETER :: IHDRHOUR = 2
        INTEGER, PARAMETER :: IHDRLAYR = 3
        INTEGER, PARAMETER :: IHDRCOL  = 4
        INTEGER, PARAMETER :: IHDRROW  = 5
        INTEGER, PARAMETER :: IHDRSRC  = 6
        INTEGER, PARAMETER :: IHDRREGN = 7
        INTEGER, PARAMETER :: IHDRCNRY = 8
        INTEGER, PARAMETER :: IHDRSTAT = 9
        INTEGER, PARAMETER :: IHDRCNTY = 10
        INTEGER, PARAMETER :: IHDRSCC  = 11
        INTEGER, PARAMETER :: IHDRSIC  = 12
        INTEGER, PARAMETER :: IHDRMACT = 13
        INTEGER, PARAMETER :: IHDRNAI  = 14
        INTEGER, PARAMETER :: IHDRSTYP = 15
        INTEGER, PARAMETER :: IHDRSRG1 = 16
        INTEGER, PARAMETER :: IHDRSRG2 = 17
        INTEGER, PARAMETER :: IHDRMON  = 18
        INTEGER, PARAMETER :: IHDRWEK  = 19
        INTEGER, PARAMETER :: IHDRDOM  = 20
        INTEGER, PARAMETER :: IHDRMND  = 21
        INTEGER, PARAMETER :: IHDRTUE  = 22
        INTEGER, PARAMETER :: IHDRWED  = 23
        INTEGER, PARAMETER :: IHDRTHU  = 24
        INTEGER, PARAMETER :: IHDRFRI  = 25
        INTEGER, PARAMETER :: IHDRSAT  = 26
        INTEGER, PARAMETER :: IHDRSUN  = 27
        INTEGER, PARAMETER :: IHDRMET  = 28
        INTEGER, PARAMETER :: IHDRSPC  = 29
        INTEGER, PARAMETER :: IHDRHT   = 30
        INTEGER, PARAMETER :: IHDRDM   = 31
        INTEGER, PARAMETER :: IHDRTK   = 32
        INTEGER, PARAMETER :: IHDRVE   = 33
        INTEGER, PARAMETER :: IHDRLAT  = 34
        INTEGER, PARAMETER :: IHDRLON  = 35
        INTEGER, PARAMETER :: IHDRELEV = 36
        INTEGER, PARAMETER :: IHDRSTKG = 37
        INTEGER, PARAMETER :: IHDRPNAM = 38
        INTEGER, PARAMETER :: IHDRSNAM = 39
        INTEGER, PARAMETER :: IHDRINAM = 40    ! SIC name
        INTEGER, PARAMETER :: IHDRMNAM = 41
        INTEGER, PARAMETER :: IHDRNNAM = 42
        INTEGER, PARAMETER :: IHDRVAR  = 43
        INTEGER, PARAMETER :: IHDRDATA = 44
        INTEGER, PARAMETER :: IHDRUNIT = 45
        INTEGER, PARAMETER :: IHDRLABL = 46
        INTEGER, PARAMETER :: IHDRNFDRS= 47
        INTEGER, PARAMETER :: IHDRMATBN= 48
        INTEGER, PARAMETER :: IHDRORIS = 49
        INTEGER, PARAMETER :: IHDRORNM = 50
        INTEGER, PARAMETER :: IHDRINTGR= 51
        INTEGER, PARAMETER :: IHDRGEO1 = 52
        INTEGER, PARAMETER :: IHDRGEO2 = 53
        INTEGER, PARAMETER :: IHDRGEO3 = 54
        INTEGER, PARAMETER :: IHDRGEO4 = 55
        INTEGER, PARAMETER :: IHDRCARB = 56
        INTEGER, PARAMETER :: IHDRCADT = 57
        INTEGER, PARAMETER :: IHDRERTYP= 58
        INTEGER, PARAMETER :: IHDPOINT = 59
        INTEGER, PARAMETER :: IHDRFUGHT= 60
        INTEGER, PARAMETER :: IHDRFUGWD= 61
        INTEGER, PARAMETER :: IHDRFUGLN= 62
        INTEGER, PARAMETER :: IHDRFUGAN= 63
        INTEGER, PARAMETER :: IHDRLAMBX= 64
        INTEGER, PARAMETER :: IHDRLAMBY= 65
        INTEGER, PARAMETER :: IHDRUTMX = 66
        INTEGER, PARAMETER :: IHDRUTMY = 67
        INTEGER, PARAMETER :: IHDRUTMZ = 68
        INTEGER, PARAMETER :: IHDRSWLAT= 69
        INTEGER, PARAMETER :: IHDRSWLON= 70
        INTEGER, PARAMETER :: IHDRNWLAT= 71
        INTEGER, PARAMETER :: IHDRNWLON= 72
        INTEGER, PARAMETER :: IHDRNELAT= 73
        INTEGER, PARAMETER :: IHDRNELON= 74
        INTEGER, PARAMETER :: IHDRSELAT= 75
        INTEGER, PARAMETER :: IHDRSELON= 76
        INTEGER, PARAMETER :: IHDRFNAM = 77   ! GSPRO name
        INTEGER, PARAMETER :: NHEADER  = 77

        CHARACTER(12), PARAMETER :: MISSNAME = 'Missing Name'

        CHARACTER(17), PARAMETER :: HEADERS( NHEADER ) = 
     &                          ( / 'Date             ',
     &                              'Hour             ',
     &                              'Layer            ',
     &                              'X cell           ',
     &                              'Y cell           ',
     &                              'Source ID        ',
     &                              'Region           ',
     &                              'Country          ',
     &                              'State            ',
     &                              'County           ',
     &                              'SCC              ',
     &                              'SIC              ',
     &                              'MACT             ',
     &                              'NAICS            ',
     &                              'Source type      ',
     &                              'Primary Srg      ',
     &                              'Fallbk Srg       ',
     &                              'Monthly Prf      ',
     &                              'Weekly  Prf      ',
     &                              'Day-Month Prf    ',
     &                              'Mon Diu Prf      ',
     &                              'Tue Diu Prf      ',
     &                              'Wed Diu Prf      ',
     &                              'Thu Diu Prf      ',
     &                              'Fri Diu Prf      ',
     &                              'Sat Diu Prf      ',
     &                              'Sun Diu Prf      ',
     &                              'Met-Diu Prf      ',
     &                              'Speciation Prf   ',
     &                              'Stk Ht           ',
     &                              'Stk Dm           ',
     &                              'Stk Tmp          ',
     &                              'Stk Vel          ',
     &                              'Latitude         ',
     &                              'Longitude        ',
     &                              'Elevstat         ',
     &                              'Stack Groups     ',
     &                              'Plt Name         ',
     &                              'SCC Description  ',
     &                              'SIC Description  ',
     &                              'MACT Description ',
     &                              'NAICS Description',
     &                              'Variable         ',
     &                              'Data value       ',
     &                              'Units            ',
     &                              'Label            ',
     &                              'NFDRS            ',
     &                              'MATBURNED        ',
     &                              'ORIS             ',
     &                              'ORIS Description ',
     &                              'INT_STAT         ',
     &                              'Geo Regn Level 1 ',
     &                              'Geo Regn Level 2 ',
     &                              'Geo Regn Level 3 ',
     &                              'Geo Regn Level 4 ',
     &                              'T,Yr,Mon,Jday,Dow',
     &                              'Basin,Dist,Cnty  ',
     &                              'Point Unit Type  ',
     &                              'Emis Release Type',
     &                              'Fug Ht           ',
     &                              'Fug Wdt          ',
     &                              'Fug Len          ',
     &                              'Fug Ang          ',
     &                              'Lambert-X        ',
     &                              'Lambert-Y        ',
     &                              'UTM_X            ',
     &                              'UTM_Y            ',
     &                              'UTM Zone         ',
     &                              'SW Latitude      ',
     &                              'SW Longitude     ',
     &                              'NW Latitude      ',
     &                              'NW Longitude     ',
     &                              'NE Latitude      ',
     &                              'NE Longitude     ',
     &                              'SE Latitude      ',
     &                              'SE Longitude     ',
     &                              'GSPRO Description' / )

C...........   Local variables that depend on module variables
        LOGICAL    LGEO1USE ( NGEOLEV1 )
        LOGICAL    LCTRYUSE ( NCOUNTRY )
        LOGICAL    LSTATUSE ( NSTATE )
        LOGICAL    LCNTYUSE ( NCOUNTY )
        LOGICAL    LSCCUSE  ( NINVSCC )
        LOGICAL    LSICUSE  ( NINVSIC )
        LOGICAL    LSPCUSE  ( NSPROF )
        LOGICAL    LMACTUSE ( NINVMACT )
        LOGICAL    LNAICSUSE( NINVNAICS )
        LOGICAL    LORISUSE ( NORIS )

        CHARACTER(10) CHRHDRS( NCHARS )  ! Source characteristics headers

C...........   Other local arrays
        INTEGER       PWIDTH( 8 )

C...........   Other local variables
        INTEGER     I, J, K, K1, K2, L, L1, L2, S, V, IOS

        INTEGER     IHDRIDX         ! tmp header index
        INTEGER     LN              ! length of single units entry
        INTEGER     LU              ! cumulative width of units header
        INTEGER     LV              ! width of delimiter
        INTEGER     NC              ! tmp no. src chars
        INTEGER     NDECI           ! no decimal place of data format
        INTEGER     NLEFT           ! value of left part of data format
        INTEGER     NWIDTH          ! tmp with
        INTEGER     W1, W2          ! tmp widths
        INTEGER     STIDX           ! starting index of loop
        INTEGER     EDIDX           ! ending index of loop
        INTEGER, SAVE:: PDEV  = 0   ! previous output file unit number

        REAL        VAL             ! tmp data value
        REAL        PREVAL          ! tmp previous data value

        LOGICAL  :: GEO1MISS              ! true: missing geo level 1 name
        LOGICAL  :: CNRYMISS              ! true: >=1 missing country name
        LOGICAL  :: CNTYMISS              ! true: >=1 missing county name
        LOGICAL  :: DATFLOAT              ! true: use float output format
        LOGICAL  :: ORISMISS              ! true: >=1 missing ORIS name
        LOGICAL  :: STATMISS              ! true: >=1 missing state name
        LOGICAL  :: SCCMISS               ! true: >=1 missing SCC name
        LOGICAL  :: SICMISS               ! true: >=1 missing SIC name
        LOGICAL  :: SPCMISS               ! true: >=1 missing GSPRO name
        LOGICAL  :: MACTMISS              ! true: >=1 missing MACT name
        LOGICAL  :: NAICSMISS             ! true: >=1 missing NAICS name
        LOGICAL  :: DYSTAT                ! true: write average day header

        LOGICAL  :: FIRSTIME = .TRUE.     ! true: first time routine called

        CHARACTER(50)  :: BUFFER      ! write buffer
        CHARACTER(50)  :: LINFMT      ! header line of '-'
        CHARACTER(300) :: MESG        ! message buffer
        CHARACTER(IODLEN3)  :: TMPUNIT    ! tmp units buffer
        CHARACTER(OLINELEN) :: HDRBUF     ! labels line buffer
        CHARACTER(OLINELEN) :: UNTBUF     ! units line buffer
        CHARACTER(QAFMTL3)  :: TMPFMT ! temporary format for Linux PG compiler

        CHARACTER(16) :: PROGNAME = 'WRREPHDR' ! program name

C***********************************************************************
C   begin body of subroutine WRREPHDR
        
C.........  Initialize output subroutine arguments
        LH     = 0
        OUTFMT = ' '

C.........  Initialize local variables for current report
        GEO1MISS = .FALSE.
        CNRYMISS = .FALSE.
        STATMISS = .FALSE.
        CNTYMISS = .FALSE.
        SCCMISS  = .FALSE.
        SICMISS  = .FALSE.
        SPCMISS  = .FALSE.
        MACTMISS = .FALSE.
        NAICSMISS= .FALSE.
        ORISMISS = .FALSE.

        LGEO1USE = .FALSE.    ! array
        LCTRYUSE = .FALSE.    ! array
        LSTATUSE = .FALSE.    ! array
        LCNTYUSE = .FALSE.    ! array
        LSCCUSE  = .FALSE.    ! array
        LSICUSE  = .FALSE.    ! array
        LSPCUSE  = .FALSE.    ! array
        LMACTUSE = .FALSE.    ! array
        LNAICSUSE= .FALSE.    ! array
        LORISUSE = .FALSE.    ! array
        
        PWIDTH   = 0          ! array
        LU       = 0

        IF( AFLAG ) THEN
            DEALLOCATE( OUTDNAM )
            ALLOCATE( OUTDNAM( RPT_%NUMDATA, RCNT ), STAT=IOS )
            CALL CHECKMEM( IOS, 'OUTDNAM', PROGNAME )
            OUTDNAM  = ''       ! array

            DO I = 1, RPT_%NUMDATA
                IF( AFLAG ) OUTDNAM( I, RCNT ) = EANAM( I )
            END DO
        END IF

C.........  Initialize report-specific settings
        RPT_ = ALLRPT( RCNT )  ! many-values

        LREGION = ( RPT_%BYGEO1 .OR. RPT_%BYCNRY .OR. 
     &              RPT_%BYSTAT .OR. RPT_%BYCNTY )

C.........  Define source-category specific header
C.........  NOTE that (1) will not be used and none will be for area sources
        CHRHDRS( 1 ) = HEADERS( IHDRREGN )
        SELECT CASE( CATEGORY )
        CASE( 'AREA' )
            CHRHDRS( 2 ) = HEADERS( IHDRSCC )

        CASE( 'MOBILE' )
            CHRHDRS( 2 ) = 'Road '
            CHRHDRS( 3 ) = 'Link'
            CHRHDRS( 4 ) = 'Veh Type'
            CHRHDRS( 5 ) = 'SCC'

        CASE( 'POINT' )
            CHRHDRS( 2 ) = 'Plant ID'
            IF ( NCHARS .GE. 3 ) THEN
                IF( .NOT. AFLAG ) THEN
                    CHRHDRS( 3 ) = 'Char 1'
                ELSE
                    CHRHDRS( 3 ) = 'Stack ID'
                END IF
            END IF
            IF ( NCHARS .GE. 4 ) CHRHDRS( 4 ) = 'Char 2'
            IF ( NCHARS .GE. 5 ) CHRHDRS( 5 ) = 'Char 3'
            IF ( NCHARS .GE. 6 ) CHRHDRS( 6 ) = 'Char 4'
            IF ( NCHARS .GE. 7 ) CHRHDRS( 7 ) = 'Char 5'

        END SELECT

C............................................................................
C.........  Pre-process output bins to determine the width of the stack 
C           parameter and variable-length string columns.
C.........  For country, state, county, SCC names, and SIC names only flag 
C           which ones are being used by the selected sources.
C............................................................................
        PDSCWIDTH = 1
        DO I = 1, NOUTBINS

C.............  Include geo code level 1 name in string
            IF( RPT_%BYGEO1NAM ) THEN
                J = BINGEO1IDX( I )
                IF( J .GT. 0 ) LGEO1USE( J ) = .TRUE.
                IF( J .LE. 0 ) GEO1MISS = .TRUE.
            END IF

C.............  Include country name in string
            IF( RPT_%BYCONAM ) THEN
                J = BINCOIDX( I )
                IF( J .GT. 0 ) LCTRYUSE( J ) = .TRUE.
                IF( J .LE. 0 ) CNRYMISS = .TRUE.
            END IF

C.............  Include state name in string
            IF( RPT_%BYSTNAM ) THEN
                J = BINSTIDX( I )
                IF( J .GT. 0 ) LSTATUSE( J ) = .TRUE.
                IF( J .LE. 0 ) STATMISS = .TRUE.
            END IF

C.............  Include county name in string
            IF( RPT_%BYCYNAM ) THEN
                J = BINCYIDX( I )
                IF( J .GT. 0 ) LCNTYUSE( J ) = .TRUE.
                IF( J .LE. 0 ) CNTYMISS = .TRUE.
            END IF

C.............  Include stack parameters
            IF( RPT_%STKPARM ) THEN
                S = BINSMKID( I )
 
                BUFFER = ' '
                WRITE( BUFFER, '(F30.0)' ) STKHT( S )
                BUFFER = ADJUSTL( BUFFER )
                PWIDTH( 1 ) = MAX( PWIDTH( 1 ), LEN_TRIM( BUFFER ) )

                BUFFER = ' '
                WRITE( BUFFER, '(F30.0)' ) STKDM( S )
                BUFFER = ADJUSTL( BUFFER )
                PWIDTH( 2 ) = MAX( PWIDTH( 2 ), LEN_TRIM( BUFFER ) )

                BUFFER = ' '
                WRITE( BUFFER, '(F30.0)' ) STKTK( S )
                BUFFER = ADJUSTL( BUFFER )
                PWIDTH( 3 ) = MAX( PWIDTH( 3 ), LEN_TRIM( BUFFER ) )

                BUFFER = ' '
                WRITE( BUFFER, '(F30.0)' ) STKVE( S )
                BUFFER = ADJUSTL( BUFFER )
                PWIDTH( 4 ) = MAX( PWIDTH( 4 ), LEN_TRIM( BUFFER ) )

            END IF

C.............  Include fugitive parameters
            IF( RPT_%FUGPARM ) THEN
                S = BINSMKID( I )

                BUFFER = ' '
                WRITE( BUFFER, '(F30.0)' ) FUGHGT( S )
                BUFFER = ADJUSTL( BUFFER )
                PWIDTH( 5 ) = MAX( PWIDTH( 5 ), LEN_TRIM( BUFFER ) )

                BUFFER = ' '
                WRITE( BUFFER, '(F30.0)' ) FUGWID( S )
                BUFFER = ADJUSTL( BUFFER )
                PWIDTH( 6 ) = MAX( PWIDTH( 6 ), LEN_TRIM( BUFFER ) )

                BUFFER = ' '
                WRITE( BUFFER, '(F30.0)' ) FUGLEN( S )
                BUFFER = ADJUSTL( BUFFER )
                PWIDTH( 7 ) = MAX( PWIDTH( 7 ), LEN_TRIM( BUFFER ) )

                BUFFER = ' '
                WRITE( BUFFER, '(F30.0)' ) FUGANG( S )
                BUFFER = ADJUSTL( BUFFER )
                PWIDTH( 8 ) = MAX( PWIDTH( 8 ), LEN_TRIM( BUFFER ) )

            END IF

C.............  Include plant description (for point sources)
            IF( RPT_%SRCNAM ) THEN
                S = BINSMKID( I )
                PDSCWIDTH = MAX( PDSCWIDTH, LEN_TRIM( CPDESC( S ) ) )
            END IF

C.............  Include SCC description
            IF( RPT_%SCCNAM ) THEN
                J = BINSNMIDX( I ) 
                IF( J .GT. 0 ) LSCCUSE( J ) = .TRUE.
            END IF

C.............  Include SIC description
            IF( RPT_%SICNAM ) THEN
                J = BINSICIDX( I ) 
                IF( J .GT. 0 ) LSICUSE( J ) = .TRUE.
            END IF

C.............  Include GSPRO description
            IF( RPT_%SPCNAM ) THEN
                J = BINSPCIDX( I )
                IF( J .GT. 0 ) LSPCUSE( J ) = .TRUE.
            END IF
 
C.............  Include MACT description
            IF( RPT_%MACTNAM ) THEN
                J = BINMACIDX( I ) 
                IF( J .GT. 0 ) LMACTUSE( J ) = .TRUE.
            END IF

C.............  Include NAICS description
            IF( RPT_%NAICSNAM ) THEN
                J = BINNAIIDX( I ) 
                IF( J .GT. 0 ) LNAICSUSE( J ) = .TRUE.
            END IF

C.............  Include ORIS description
            IF( RPT_%ORISNAM ) THEN
                J = BINORSIDX( I ) 
                IF( J .GT. 0 ) LORISUSE( J ) = .TRUE.
            END IF

       END DO  ! End loop through bins

C............................................................................
C.........  Set the widths of each output column, while including the
C           width of the column header.
C.........  Build the formats for the data in each column
C.........  Build the header as we go along
C............................................................................

C.........  The extra length for each variable is 1 space and 1 delimiter width
        LV = LEN_TRIM( RPT_%DELIM ) + 1

C.........  User-defined label
        IF( RPT_%USELABEL ) THEN
            J = MAX( LEN_TRIM( RPT_%LABEL ), 
     &               LEN_TRIM( HEADERS(IHDRLABL) ) )
            LABELWIDTH = J + LV

            CALL ADD_TO_HEADER( J, HEADERS(IHDRLABL), LH, HDRBUF )
            CALL ADD_TO_HEADER( J, ' ', LU, UNTBUF )

        END IF

C.........  Date column
        IF( RPT_%CARB ) THEN
            J = 17  ! header width is AGTYPE,YEAR,MM,JDAY,DOW
            WRITE( DATEFMT, 94220 ) RPT_%DELIM  ! leading zeros
            DATEWIDTH = J + LV

            CALL ADD_TO_HEADER( J, HEADERS(IHDRCARB), LH, HDRBUF )
            CALL ADD_TO_HEADER( J, ' ', LU, UNTBUF )

        END IF

C.........  Date column
        IF( RPT_%BYDATE .AND. .NOT. RPT_%CARB ) THEN
            J = 10  ! header width is MM/DD/YYYY
            WRITE( DATEFMT, 94620 ) RPT_%DELIM  ! leading zeros
            DATEWIDTH = J + LV

            CALL ADD_TO_HEADER( J, HEADERS(IHDRDATE), LH, HDRBUF )
            CALL ADD_TO_HEADER( J, ' ', LU, UNTBUF )

        END IF

C.........  Hour column
        IF( RPT_%BYHOUR ) THEN
            IF( .NOT. DLFLAG ) THEN
                J = LEN_TRIM( HEADERS( IHDRHOUR ) )  ! header width
                WRITE( HOURFMT, 94630 ) J, 2, RPT_%DELIM  ! leading zeros
                J = MAX( 2, J )
                HOURWIDTH = J + LV

                CALL ADD_TO_HEADER( J, HEADERS(IHDRHOUR), LH, HDRBUF )
                CALL ADD_TO_HEADER( J, ' ', LU, UNTBUF )

            END IF
        ELSE IF( .NOT. RPT_%BYHOUR .AND. RPT_%CARB ) THEN
            IF( .NOT. DLFLAG ) THEN
                J = LEN_TRIM( HEADERS( IHDRHOUR ) )  ! header width
                WRITE( HOURFMT, 94630 ) J, 2, RPT_%DELIM  ! leading zeros
                J = MAX( 2, J )
                HOURWIDTH = J + LV

                CALL ADD_TO_HEADER( J, HEADERS(IHDRHOUR), LH, HDRBUF )
                CALL ADD_TO_HEADER( J, ' ', LU, UNTBUF )

            END IF
        END IF

C.........  Layer column
        IF( RPT_%BYLAYER ) THEN
            J = LEN_TRIM( HEADERS( IHDRLAYR ) )  ! header width
            WRITE( LAYRFMT, 94630 ) J, 2, RPT_%DELIM  ! leading zeros
            J = MAX( 2, J )
            LAYRWIDTH = J + LV

            CALL ADD_TO_HEADER( J, HEADERS(IHDRLAYR), LH, HDRBUF )
            CALL ADD_TO_HEADER( J, ' ', LU, UNTBUF )

        END IF

C.........  Cell columns
        IF( RPT_%BYCELL ) THEN

C.............  X-cell
            J = LEN_TRIM( HEADERS( IHDRCOL ) )
            W1 = INTEGER_COL_WIDTH( NOUTBINS, BINX )
            W1 = MAX( W1, J )
            CALL ADD_TO_HEADER( W1, HEADERS(IHDRCOL), LH, HDRBUF )
            CALL ADD_TO_HEADER( W1, ' ', LU, UNTBUF )

C.............  Y-cell
            J = LEN_TRIM( HEADERS( IHDRROW ) )
            W2 = INTEGER_COL_WIDTH( NOUTBINS, BINY )
            W2 = MAX( W2, J )

            CALL ADD_TO_HEADER( W2, HEADERS(IHDRROW), LH, HDRBUF )
            CALL ADD_TO_HEADER( W2, ' ', LU, UNTBUF )

C.............  Write format to include both x-cell and y-cell
            WRITE( CELLFMT, 94635 ) W1, RPT_%DELIM, W2, RPT_%DELIM
            CELLWIDTH = W1 + W2 + 2*LV

        ELSE IF( .NOT. RPT_%BYCELL .AND. RPT_%CARB ) THEN

C.............  X-cell
            J = LEN_TRIM( HEADERS( IHDRCOL ) )
            W1 = J
            CALL ADD_TO_HEADER( W1, HEADERS(IHDRCOL), LH, HDRBUF )
            CALL ADD_TO_HEADER( W1, ' ', LU, UNTBUF )

C.............  Y-cell
            J = LEN_TRIM( HEADERS( IHDRROW ) )
            W2 = J
            CALL ADD_TO_HEADER( W2, HEADERS(IHDRROW), LH, HDRBUF )
            CALL ADD_TO_HEADER( W2, ' ', LU, UNTBUF )

C.............  Write format to include both x-cell and y-cell
            WRITE( CELLFMT, 94635 ) W1, RPT_%DELIM, W2, RPT_%DELIM
            CELLWIDTH = W1 + W2 + 2*LV

        END IF

C.........  Source ID column
        IF( RPT_%BYSRC ) THEN

            J = LEN_TRIM( HEADERS( IHDRSRC ) )
            W1 = INTEGER_COL_WIDTH( NOUTBINS, BINSMKID )
            W1 = MAX( W1, J )

            CALL ADD_TO_HEADER( W1, HEADERS(IHDRSRC), LH, HDRBUF )
            CALL ADD_TO_HEADER( W1, ' ', LU, UNTBUF )

            WRITE( SRCFMT, 94625 ) W1, RPT_%DELIM
            SRCWIDTH = W1 + LV

        END IF

C.........  Region code column
        IF( LREGION ) THEN
            IF( RPT_%CARB ) THEN
                J  = LEN_TRIM( HEADERS( IHDRCADT ) )
            ELSE
                J  = LEN_TRIM( HEADERS( IHDRREGN ) )
            ENDIF
            W1 = 0
            DO I = 1, NOUTBINS
                W1 = MAX( W1, LEN_TRIM( BINREGN( I ) ) )
            END DO
            W1  = MAX( W1, J )

            IF( RPT_%CARB ) THEN
                CALL ADD_TO_HEADER( W1, HEADERS(IHDRCADT), LH, HDRBUF)
            ELSE
                CALL ADD_TO_HEADER( W1, HEADERS(IHDRREGN), LH, HDRBUF)
            ENDIF
            CALL ADD_TO_HEADER( W1, ' ', LU, UNTBUF )

            WRITE( REGNFMT, 94630 ) W1, FIPLEN3, RPT_%DELIM     ! leading zeros
            REGNWIDTH = W1 + LV

        ELSE IF( .NOT. LREGION .AND. RPT_%CARB ) THEN
            J  = LEN_TRIM( HEADERS( IHDRCADT ) )
            W1 = J

            CALL ADD_TO_HEADER( W1, HEADERS(IHDRCADT), LH, HDRBUF)
            CALL ADD_TO_HEADER( W1, ' ', LU, UNTBUF )

            WRITE( REGNFMT, 94630 ) W1, FIPLEN3, RPT_%DELIM     ! leading zeros
            REGNWIDTH = W1 + LV

        END IF

C.........  Set widths and build formats for country, state, and county names.
C           These are done on loops of unique lists of these names
C           so that the LEN_TRIMs can be done on the shortest possible list
C           of entries instead of on all entries in the bins list.

C.........  Geo code level 1 names
        IF( RPT_%BYGEO1NAM ) THEN

C.............  For regions in the inventory, get max name width
            NWIDTH = 0
            DO I = 1, NGEOLEV1
                IF( LGEO1USE( I ) ) THEN
                    NWIDTH = MAX( NWIDTH, LEN_TRIM( GEOLEV1NAM( I ) ) )
                END IF
            END DO

C.............  If any missing region names, check widths
            IF( GEO1MISS ) NWIDTH = MAX( NWIDTH, LEN_TRIM( MISSNAME ) )

C.............  Set geo code name column width 
            J = LEN_TRIM( HEADERS( IHDRGEO1 ) )
            J = MAX( NWIDTH, J )

            CALL ADD_TO_HEADER( J, HEADERS(IHDRGEO1), LH, HDRBUF )
            CALL ADD_TO_HEADER( J, ' ', LU, UNTBUF )

            GEO1WIDTH = J + LV

        END IF

C.........  Country names
        IF( RPT_%BYCONAM ) THEN

C.............  For countries in the inventory, get max name width
            NWIDTH = 0
            DO I = 1, NCOUNTRY
                IF( LCTRYUSE( I ) ) THEN
                    NWIDTH = MAX( NWIDTH, LEN_TRIM( CTRYNAM( I ) ) )
                END IF
            END DO

C.............  If any missing country names, check widths
            IF( CNRYMISS ) NWIDTH = MAX( NWIDTH, LEN_TRIM( MISSNAME ) )

C.............  Set country name column width 
            IHDRIDX = IHDRCNRY
            IF( USEEXPGEO() ) IHDRIDX = IHDRGEO2
            J = LEN_TRIM( HEADERS( IHDRIDX ) )
            J = MAX( NWIDTH, J )

            CALL ADD_TO_HEADER( J, HEADERS(IHDRIDX), LH, HDRBUF )
            CALL ADD_TO_HEADER( J, ' ', LU, UNTBUF )

            COWIDTH = J + LV

        END IF

C.........  State names
        IF( RPT_%BYSTNAM ) THEN

C.............  For states in the inventory, get max name width
            NWIDTH = 0
            DO I = 1, NSTATE
                IF( LSTATUSE( I ) ) THEN
                    NWIDTH = MAX( NWIDTH, LEN_TRIM( STATNAM( I ) ) )
                END IF
            END DO

C.............  If any missing state names, check widths
            IF( STATMISS ) NWIDTH = MAX( NWIDTH, LEN_TRIM( MISSNAME ) )

C.............  Set state name column width 
            IHDRIDX = IHDRSTAT
            IF( USEEXPGEO() ) IHDRIDX = IHDRGEO3
            J = LEN_TRIM( HEADERS( IHDRIDX ) )
            J = MAX( NWIDTH, J )

            CALL ADD_TO_HEADER( J, HEADERS(IHDRIDX), LH, HDRBUF )
            CALL ADD_TO_HEADER( J, ' ', LU, UNTBUF )

            STWIDTH = J + LV

        END IF

C.........  County names
        IF( RPT_%BYCYNAM ) THEN

C.............  For counties in the inventory, get max name width
            NWIDTH = 0
            DO I = 1, NCOUNTY
                IF( LCNTYUSE( I ) ) THEN
                    NWIDTH = MAX( NWIDTH, LEN_TRIM( CNTYNAM( I ) ) )
                END IF
            END DO

C.............  If any missing county names, check widths
            IF( CNTYMISS ) NWIDTH = MAX( NWIDTH, LEN_TRIM( MISSNAME ) )

C.............  Set county name column width 
            IHDRIDX = IHDRCNTY
            IF( USEEXPGEO() ) IHDRIDX = IHDRGEO4
            J = LEN_TRIM( HEADERS( IHDRIDX ) )
            J = MAX( NWIDTH, J )

            CALL ADD_TO_HEADER( J, HEADERS(IHDRIDX), LH, HDRBUF )
            CALL ADD_TO_HEADER( J, ' ', LU, UNTBUF )

            CYWIDTH = J + LV

        END IF

C.........  SCC column
        IF( RPT_%BYSCC ) THEN
            J = LEN_TRIM( HEADERS( IHDRSCC ) )
            IF( RPT_%SCCRES < NSCCLV3-1 ) J = J + 7  ! Plus " Tier #"
            J = MAX( SCCLEN3, J )
   
            IF( RPT_%SCCRES < NSCCLV3-1 ) THEN
                WRITE( BUFFER, '(A,I1)' ) TRIM( HEADERS(IHDRSCC) ) // 
     &                                    ' Tier ', RPT_%SCCRES
            ELSE
                BUFFER = TRIM( HEADERS(IHDRSCC) )
            END IF

            CALL ADD_TO_HEADER( J, BUFFER, LH, HDRBUF )
            CALL ADD_TO_HEADER( J, ' ', LU, UNTBUF )

            SCCWIDTH = J + LV
        END IF

C.........  SIC column
        IF( RPT_%BYSIC ) THEN
            IF( MATFLAG ) THEN
                J = LEN_TRIM( HEADERS( IHDRMATBN ) )
                J = MAX( SICLEN3, J )
                
                CALL ADD_TO_HEADER( J, HEADERS(IHDRMATBN), LH, HDRBUF )
                CALL ADD_TO_HEADER( J, ' ', LU, UNTBUF )

                SICWIDTH = J + LV
            ELSE
                J = LEN_TRIM( HEADERS( IHDRSIC ) )
                J = MAX( SICLEN3, J )
                
                CALL ADD_TO_HEADER( J, HEADERS(IHDRSIC), LH, HDRBUF )
                CALL ADD_TO_HEADER( J, ' ', LU, UNTBUF )

                SICWIDTH = J + LV
            END IF
        END IF

C.........  INTEGRATE column
        IF( RPT_%BYINTGR ) THEN
            J = LEN_TRIM( HEADERS( IHDRINTGR ) )
            J = MAX( INTLEN3, J )
            CALL ADD_TO_HEADER( J, HEADERS(IHDRINTGR), LH, HDRBUF )
            CALL ADD_TO_HEADER( J, ' ', LU, UNTBUF )
            INTGRWIDTH = J + LV
        END IF

C.........  MACT column
        IF( RPT_%BYMACT ) THEN
            IF( NFDFLAG ) THEN
                J = LEN_TRIM( HEADERS( IHDRNFDRS ) )
                J = MAX( MACLEN3, J )
                CALL ADD_TO_HEADER( J, HEADERS(IHDRNFDRS), LH, HDRBUF )
                CALL ADD_TO_HEADER( J, ' ', LU, UNTBUF )
                MACTWIDTH = J + LV
            ELSE
                J = LEN_TRIM( HEADERS( IHDRMACT ) )
                J = MAX( MACLEN3, J )
                CALL ADD_TO_HEADER( J, HEADERS(IHDRMACT), LH, HDRBUF )
                CALL ADD_TO_HEADER( J, ' ', LU, UNTBUF )
                MACTWIDTH = J + LV
            END IF
        END IF

C.........  NAICS column
        IF( RPT_%BYNAICS ) THEN
            J = LEN_TRIM( HEADERS( IHDRNAI ) )
            J = MAX( NAILEN3, J )
    
            CALL ADD_TO_HEADER( J, HEADERS(IHDRNAI), LH, HDRBUF )
            CALL ADD_TO_HEADER( J, ' ', LU, UNTBUF )

            NAIWIDTH = J + LV
        END IF

C.........  SRCTYP column
        IF( RPT_%BYSRCTYP ) THEN
            J = LEN_TRIM( HEADERS( IHDRSTYP ) )
            J = MAX( STPLEN3, J )
    
            CALL ADD_TO_HEADER( J, HEADERS(IHDRSTYP), LH, HDRBUF )
            CALL ADD_TO_HEADER( J, ' ', LU, UNTBUF )

            STYPWIDTH = J + LV
        END IF

C.........  Primary surrogates column
        IF( RPT_%SRGRES .EQ. 1 ) THEN
            J = LEN_TRIM( HEADERS( IHDRSRG1 ) )
            W1 = INTEGER_COL_WIDTH( NOUTBINS, BINSRGID1 )
            W1  = MAX( W1, J )
            CALL ADD_TO_HEADER( W1, HEADERS(IHDRSRG1), LH, HDRBUF )
            CALL ADD_TO_HEADER( W1, ' ', LU, UNTBUF )

            WRITE( SRG1FMT, 94650 ) W1, RPT_%DELIM 
            SRG1WIDTH = W1 + LV
        END IF

C.........  Fallback surrogates column
        IF( RPT_%SRGRES .GE. 1 ) THEN
            J = LEN_TRIM( HEADERS( IHDRSRG2 ) )
            W1 = INTEGER_COL_WIDTH( NOUTBINS, BINSRGID2 )
            W1  = MAX( W1, J )
            CALL ADD_TO_HEADER( W1, HEADERS(IHDRSRG2), LH, HDRBUF )
            CALL ADD_TO_HEADER( W1, ' ', LU, UNTBUF )

            WRITE( SRG2FMT, 94650 ) W1, RPT_%DELIM 
            SRG2WIDTH = W1 + LV
        END IF

C.........  Temporal profiles columns
        IF( RPT_%BYMON ) THEN          ! Monthly

            NWIDTH = 0
            DO I = 1, NOUTBINS
                NWIDTH = MAX( NWIDTH, LEN_TRIM( BINMONID( I ) ) )
            END DO

C.............  Set speciation profiles column width 
            J = LEN_TRIM( HEADERS( IHDRMON ) )
            J = MAX( NWIDTH, J )

            CALL ADD_TO_HEADER( J, HEADERS(IHDRMON), LH, HDRBUF )
            CALL ADD_TO_HEADER( J, ' ', LU, UNTBUF )

            MONWIDTH = J + LV

        END IF

C.........  Temporal profiles columns
        IF( RPT_%BYWEK ) THEN          ! Weekly

            NWIDTH = 0
            DO I = 1, NOUTBINS
                NWIDTH = MAX( NWIDTH, LEN_TRIM( BINWEKID( I ) ) )
            END DO

C.............  Set speciation profiles column width 
            J = LEN_TRIM( HEADERS( IHDRWEK ) )
            J = MAX( NWIDTH, J )

            CALL ADD_TO_HEADER( J, HEADERS(IHDRWEK), LH, HDRBUF )
            CALL ADD_TO_HEADER( J, ' ', LU, UNTBUF )

            WEKWIDTH = J + LV

        END IF

C.........  Temporal profiles columns
        IF( RPT_%BYDOM ) THEN          ! Day of Month 

            NWIDTH = 0
            DO I = 1, NOUTBINS
                NWIDTH = MAX( NWIDTH, LEN_TRIM( BINDOMID( I ) ) )
            END DO

C.............  Set speciation profiles column width 
            J = LEN_TRIM( HEADERS( IHDRDOM ) )
            J = MAX( NWIDTH, J )

            CALL ADD_TO_HEADER( J, HEADERS(IHDRDOM), LH, HDRBUF )
            CALL ADD_TO_HEADER( J, ' ', LU, UNTBUF )

            DOMWIDTH = J + LV

        END IF

C.........  Temporal profiles columns
        IF( RPT_%BYMND ) THEN          ! Monday

            NWIDTH = 0
            DO I = 1, NOUTBINS
                NWIDTH = MAX( NWIDTH, LEN_TRIM( BINMNDID( I ) ) )
            END DO

C.............  Set speciation profiles column width 
            J = LEN_TRIM( HEADERS( IHDRMND ) )
            J = MAX( NWIDTH, J )

            CALL ADD_TO_HEADER( J, HEADERS(IHDRMND), LH, HDRBUF )
            CALL ADD_TO_HEADER( J, ' ', LU, UNTBUF )

            MNDWIDTH = J + LV

        END IF

C.........  Temporal profiles columns
        IF( RPT_%BYTUE ) THEN          ! Tuesday

            NWIDTH = 0
            DO I = 1, NOUTBINS
                NWIDTH = MAX( NWIDTH, LEN_TRIM( BINTUEID( I ) ) )
            END DO

C.............  Set speciation profiles column width 
            J = LEN_TRIM( HEADERS( IHDRTUE ) )
            J = MAX( NWIDTH, J )

            CALL ADD_TO_HEADER( J, HEADERS(IHDRTUE), LH, HDRBUF )
            CALL ADD_TO_HEADER( J, ' ', LU, UNTBUF )

            TUEWIDTH = J + LV

        END IF

C.........  Temporal profiles columns
        IF( RPT_%BYWED ) THEN          ! Wednesday

            NWIDTH = 0
            DO I = 1, NOUTBINS
                NWIDTH = MAX( NWIDTH, LEN_TRIM( BINWEDID( I ) ) )
            END DO

C.............  Set speciation profiles column width 
            J = LEN_TRIM( HEADERS( IHDRWED ) )
            J = MAX( NWIDTH, J )

            CALL ADD_TO_HEADER( J, HEADERS(IHDRWED), LH, HDRBUF )
            CALL ADD_TO_HEADER( J, ' ', LU, UNTBUF )

            WEDWIDTH = J + LV

        END IF

C.........  Temporal profiles columns
        IF( RPT_%BYTHU ) THEN          ! Thursday

            NWIDTH = 0
            DO I = 1, NOUTBINS
                NWIDTH = MAX( NWIDTH, LEN_TRIM( BINTHUID( I ) ) )
            END DO

C.............  Set speciation profiles column width 
            J = LEN_TRIM( HEADERS( IHDRTHU ) )
            J = MAX( NWIDTH, J )

            CALL ADD_TO_HEADER( J, HEADERS(IHDRTHU), LH, HDRBUF )
            CALL ADD_TO_HEADER( J, ' ', LU, UNTBUF )

            THUWIDTH = J + LV

        END IF

C.........  Temporal profiles columns
        IF( RPT_%BYFRI ) THEN          ! Friday

            NWIDTH = 0
            DO I = 1, NOUTBINS
                NWIDTH = MAX( NWIDTH, LEN_TRIM( BINFRIID( I ) ) )
            END DO

C.............  Set speciation profiles column width 
            J = LEN_TRIM( HEADERS( IHDRFRI ) )
            J = MAX( NWIDTH, J )

            CALL ADD_TO_HEADER( J, HEADERS(IHDRFRI), LH, HDRBUF )
            CALL ADD_TO_HEADER( J, ' ', LU, UNTBUF )

            FRIWIDTH = J + LV

        END IF

C.........  Temporal profiles columns
        IF( RPT_%BYSAT ) THEN          ! Saturday

            NWIDTH = 0
            DO I = 1, NOUTBINS
                NWIDTH = MAX( NWIDTH, LEN_TRIM( BINSATID( I ) ) )
            END DO

C.............  Set speciation profiles column width 
            J = LEN_TRIM( HEADERS( IHDRSAT ) )
            J = MAX( NWIDTH, J )

            CALL ADD_TO_HEADER( J, HEADERS(IHDRSAT), LH, HDRBUF )
            CALL ADD_TO_HEADER( J, ' ', LU, UNTBUF )

            SATWIDTH = J + LV

        END IF

C.........  Temporal profiles columns
        IF( RPT_%BYSUN ) THEN          ! Sunday

            NWIDTH = 0
            DO I = 1, NOUTBINS
                NWIDTH = MAX( NWIDTH, LEN_TRIM( BINSUNID( I ) ) )
            END DO

C.............  Set speciation profiles column width 
            J = LEN_TRIM( HEADERS( IHDRSUN ) )
            J = MAX( NWIDTH, J )

            CALL ADD_TO_HEADER( J, HEADERS(IHDRSUN), LH, HDRBUF )
            CALL ADD_TO_HEADER( J, ' ', LU, UNTBUF )

            SUNWIDTH = J + LV

        END IF

C.........  Temporal profiles columns
        IF( RPT_%BYMET ) THEN          ! Met-based

            NWIDTH = 0
            DO I = 1, NOUTBINS
                NWIDTH = MAX( NWIDTH, LEN_TRIM( BINMETID( I ) ) )
            END DO

C.............  Set speciation profiles column width 
            J = LEN_TRIM( HEADERS( IHDRMET ) )
            J = MAX( NWIDTH, J )

            CALL ADD_TO_HEADER( J, HEADERS(IHDRMET), LH, HDRBUF )
            CALL ADD_TO_HEADER( J, ' ', LU, UNTBUF )

            METWIDTH = J + LV

        END IF

C.........  Speciation profile column
        IF( RPT_%BYSPC ) THEN

            NWIDTH = 0
            DO I = 1, NOUTBINS
                NWIDTH = MAX( NWIDTH, LEN_TRIM( BINSPCID( I ) ) )
            END DO

C.............  Set speciation profiles column width 
            J = LEN_TRIM( HEADERS( IHDRSPC ) )
            J = MAX( NWIDTH, J )

            CALL ADD_TO_HEADER( J, HEADERS(IHDRSPC), LH, HDRBUF )
            CALL ADD_TO_HEADER( J, ' ', LU, UNTBUF )

            SPCWIDTH = J + LV

        END IF

C.........  Road class.  By roadclass can only be true if by source is not
C           being used.
        IF( RPT_%BYRCL ) THEN
            J  = LEN_TRIM( CHRHDRS( 2 ) )
            W1 = INTEGER_COL_WIDTH( NOUTBINS, BINRCL )
            W1  = MAX( W1, J )

            CALL ADD_TO_HEADER( W1, CHRHDRS( 2 ), LH, HDRBUF )
            CALL ADD_TO_HEADER( W1, ' ', LU, UNTBUF )

            WRITE( CHARFMT, 94645 ) W1, RPT_%DELIM 
            CHARWIDTH = W1 + LV

        END IF

C.........  Source characteristics. NOTE - the source characteristics have
C           already been rearranged and their widths reset based on the
C           inventory.  The SCC has been removed if its one of the source
C           characteristics, and NCHARS reset accordingly.
        IF( RPT_%BYSRC ) THEN        

            CHARWIDTH = 0
            CHARFMT = '('

            DO K = MINC, NCHARS

C.................  Build source characteristics output format for WRREPOUT
                TMPFMT = CHARFMT
                L  = LEN_TRIM( TMPFMT )
                J  = LEN_TRIM( CHRHDRS( K ) )
                W1 = MAX( LOC_ENDP( K ) - LOC_BEGP( K ) + 1, J )
                WRITE( CHARFMT, '(A,I2.2,A)' ) TMPFMT( 1:L )// 
     &                 '1X,A', W1, ',"'//RPT_%DELIM//'",'

                CALL ADD_TO_HEADER( W1, CHRHDRS( K ), LH, HDRBUF )
                CALL ADD_TO_HEADER( W1, ' ', LU, UNTBUF )

                CHARWIDTH = CHARWIDTH + W1 + LV

            END DO

            TMPFMT = CHARFMT
            L = LEN_TRIM( TMPFMT ) - 1   ! (minus 1 to remove trailing comma)
            CHARFMT = TMPFMT( 1:L ) // ')'

        END IF

C.........  Plant ID
        IF( RPT_%BYPLANT ) THEN
            NWIDTH = 0
            DO I = 1, NOUTBINS
                NWIDTH = MAX( NWIDTH, LEN_TRIM( BINPLANT( I ) ) )
            END DO

            J  = LEN_TRIM( CHRHDRS( 2 ) )
            W1 = MAX( NWIDTH, J )

            CALL ADD_TO_HEADER( W1, CHRHDRS( 2 ), LH, HDRBUF )
            CALL ADD_TO_HEADER( W1, ' ', LU, UNTBUF )

            WRITE( CHARFMT, 94645 ) W1, RPT_%DELIM
            CHARWIDTH = W1 + LV

        END IF

C.........  Point Unit ID
        IF( RPT_%BYUNIT ) THEN
            NWIDTH = 0
            DO I = 1, NOUTBINS
                NWIDTH = MAX( NWIDTH, LEN_TRIM( BINUNIT( I ) ) )
            END DO

            J  = LEN_TRIM( CHRHDRS( 3 ) )
            W1 = MAX( NWIDTH, J )

            CALL ADD_TO_HEADER( W1, CHRHDRS( 3 ), LH, HDRBUF )
            CALL ADD_TO_HEADER( W1, ' ', LU, UNTBUF )

            WRITE( CHARFMT, 94645 ) W1, RPT_%DELIM
            CHARWIDTH = W1 + LV

        END IF

C.........  ORIS ID
        IF( RPT_%BYORIS ) THEN
            NWIDTH = 0
            DO I = 1, NOUTBINS
                NWIDTH = MAX( NWIDTH, LEN_TRIM( BINORIS( I ) ) )
            END DO

            J = LEN_TRIM( HEADERS( IHDRORIS ) )
            J = MAX( NWIDTH, J )

            CALL ADD_TO_HEADER( J, HEADERS( IHDRORIS ), LH, HDRBUF )
            CALL ADD_TO_HEADER( J, ' ', LU, UNTBUF )

            ORSWIDTH = J + LV

        END IF

C.........  Stack parameters.  +3 for decimal and 2 significant figures
        IF( RPT_%STKPARM ) THEN

            J = LEN_TRIM( HEADERS( IHDRHT ) )
            PWIDTH( 1 ) = 10
            CALL ADD_TO_HEADER( PWIDTH( 1 ), HEADERS( IHDRHT ), 
     &                          LH, HDRBUF )
            CALL ADD_TO_HEADER( PWIDTH( 1 ), ATTRUNIT( 6 ), LU, UNTBUF )

            J = LEN_TRIM( HEADERS( IHDRDM ) )
            PWIDTH( 2 ) = 10
            CALL ADD_TO_HEADER( PWIDTH( 2 ), HEADERS( IHDRDM ), 
     &                          LH, HDRBUF )
            CALL ADD_TO_HEADER( PWIDTH( 2 ), ATTRUNIT( 7 ), LU, UNTBUF )

            J = LEN_TRIM( HEADERS( IHDRTK ) )
            PWIDTH( 3 ) = 10
            CALL ADD_TO_HEADER( PWIDTH( 3 ), HEADERS( IHDRTK ), 
     &                          LH, HDRBUF )
            CALL ADD_TO_HEADER( PWIDTH( 3 ), ATTRUNIT( 8 ), LU, UNTBUF )

            J = LEN_TRIM( HEADERS( IHDRVE ) )
            PWIDTH( 4 ) = 10
            CALL ADD_TO_HEADER( PWIDTH( 4 ), HEADERS( IHDRVE ), 
     &                          LH, HDRBUF )
            CALL ADD_TO_HEADER( PWIDTH( 4 ), ATTRUNIT( 9 ), LU, UNTBUF )

            WRITE( STKPFMT, 94640 ) PWIDTH( 1 ), RPT_%DELIM,
     &                              PWIDTH( 2 ), RPT_%DELIM,
     &                              PWIDTH( 3 ), RPT_%DELIM,
     &                              PWIDTH( 4 ), RPT_%DELIM

            STKPWIDTH = SUM( PWIDTH( 1:4 ) ) + 4*LV

        END IF

C.........  Fugitive parameters.  +3 for decimal and 2 significant figures
        IF( RPT_%FUGPARM ) THEN

            J = LEN_TRIM( HEADERS( IHDRFUGHT ) )
            PWIDTH( 5 ) = 10 
            CALL ADD_TO_HEADER( PWIDTH( 5 ), HEADERS( IHDRFUGHT ),
     &                          LH, HDRBUF )
            CALL ADD_TO_HEADER( PWIDTH( 5 ), ATTRUNIT( 6 ), LU, UNTBUF )

            J = LEN_TRIM( HEADERS( IHDRFUGWD ) )
            PWIDTH( 6 ) = 10 
            CALL ADD_TO_HEADER( PWIDTH( 6 ), HEADERS( IHDRFUGWD ),
     &                          LH, HDRBUF )
            CALL ADD_TO_HEADER( PWIDTH( 6 ), ATTRUNIT( 7 ), LU, UNTBUF )

            J = LEN_TRIM( HEADERS( IHDRFUGLN ) )
            PWIDTH( 7 ) = 10
            CALL ADD_TO_HEADER( PWIDTH( 7 ), HEADERS( IHDRFUGLN ),
     &                          LH, HDRBUF )
            CALL ADD_TO_HEADER( PWIDTH( 7 ), ATTRUNIT( 8 ), LU, UNTBUF )

            J = LEN_TRIM( HEADERS( IHDRFUGAN ) )
            PWIDTH( 8 ) = 10
            CALL ADD_TO_HEADER( PWIDTH( 8 ), HEADERS( IHDRFUGAN ),
     &                          LH, HDRBUF )
            CALL ADD_TO_HEADER( PWIDTH( 8 ), ATTRUNIT( 9 ), LU, UNTBUF )

            WRITE( FUGPFMT, 94640 ) PWIDTH( 5 ), RPT_%DELIM,
     &                              PWIDTH( 6 ), RPT_%DELIM,
     &                              PWIDTH( 7 ), RPT_%DELIM,
     &                              PWIDTH( 8 ), RPT_%DELIM

            FUGPWIDTH = SUM( PWIDTH( 5:8 ) ) + 4*LV

        END IF

C.........  Point-source latitude and longitude
        IF( RPT_%LATLON ) THEN
        
            J = LEN_TRIM( HEADERS( IHDRLAT ) )
            PWIDTH( 1 ) = 13
            CALL ADD_TO_HEADER( PWIDTH( 1 ), HEADERS( IHDRLAT ),
     &                          LH, HDRBUF )
            CALL ADD_TO_HEADER( PWIDTH( 1 ), '    ', LU, UNTBUF )
            
            J = LEN_TRIM( HEADERS( IHDRLON ) )
            PWIDTH( 2 ) = 13
            CALL ADD_TO_HEADER( PWIDTH( 2 ), HEADERS( IHDRLON ),
     &                          LH, HDRBUF )
            CALL ADD_TO_HEADER( PWIDTH( 2 ), '     ', LU, UNTBUF )
            
            WRITE( LTLNFMT, 94642 ) PWIDTH( 1 ), RPT_%DELIM,
     &                              PWIDTH( 2 ), RPT_%DELIM
     
            LTLNWIDTH = SUM( PWIDTH( 1:2 ) ) + 2*LV
            
        END IF

C.........  Point-source grid lamber-x&y, and grid utm_x,y,zone 
        IF( RPT_%GRDCOR ) THEN

            J = LEN_TRIM( HEADERS( IHDRLAMBX ) )
            PWIDTH( 1 ) = 13
            CALL ADD_TO_HEADER( PWIDTH( 1 ), HEADERS( IHDRLAMBX ),
     &                          LH, HDRBUF )
            CALL ADD_TO_HEADER( PWIDTH( 1 ), '    ', LU, UNTBUF )

            J = LEN_TRIM( HEADERS( IHDRLAMBY ) )
            PWIDTH( 2 ) = 13
            CALL ADD_TO_HEADER( PWIDTH( 2 ), HEADERS( IHDRLAMBY ),
     &                          LH, HDRBUF )
            CALL ADD_TO_HEADER( PWIDTH( 2 ), '     ', LU, UNTBUF )

            J = LEN_TRIM( HEADERS( IHDRUTMX ) )
            PWIDTH( 3 ) = 13
            CALL ADD_TO_HEADER( PWIDTH( 3 ), HEADERS( IHDRUTMX ),
     &                          LH, HDRBUF )
            CALL ADD_TO_HEADER( PWIDTH( 3 ), '    ', LU, UNTBUF )

            J = LEN_TRIM( HEADERS( IHDRUTMY ) )
            PWIDTH( 4 ) = 13
            CALL ADD_TO_HEADER( PWIDTH( 4 ), HEADERS( IHDRUTMY ),
     &                          LH, HDRBUF )
            CALL ADD_TO_HEADER( PWIDTH( 4 ), '     ', LU, UNTBUF )

            J = LEN_TRIM( HEADERS( IHDRUTMZ ) )
            PWIDTH( 5 ) = 13
            CALL ADD_TO_HEADER( PWIDTH( 5 ), HEADERS( IHDRUTMZ ),
     &                          LH, HDRBUF )
            CALL ADD_TO_HEADER( PWIDTH( 5 ), '     ', LU, UNTBUF )

            WRITE( LAMBFMT, 94643 ) PWIDTH( 1 ), RPT_%DELIM,
     &                              PWIDTH( 2 ), RPT_%DELIM,
     &                              PWIDTH( 3 ), RPT_%DELIM,
     &                              PWIDTH( 4 ), RPT_%DELIM,
     &                              PWIDTH( 5 ), RPT_%DELIM

            LAMBWIDTH = SUM( PWIDTH( 1:5 ) ) + 5*LV

        END IF

C.........  Grid cell corner coordinates
        IF( RPT_%GRDPNT ) THEN
        
            J = LEN_TRIM( HEADERS( IHDRSWLAT ) )
            PWIDTH( 1 ) = 13
            CALL ADD_TO_HEADER( PWIDTH( 1 ), HEADERS( IHDRSWLAT ),
     &                          LH, HDRBUF )
            CALL ADD_TO_HEADER( PWIDTH( 1 ), '    ', LU, UNTBUF )
        
            J = LEN_TRIM( HEADERS( IHDRSWLON ) )
            PWIDTH( 2 ) = 13
            CALL ADD_TO_HEADER( PWIDTH( 2 ), HEADERS( IHDRSWLON ),
     &                          LH, HDRBUF )
            CALL ADD_TO_HEADER( PWIDTH( 2 ), '    ', LU, UNTBUF )
        
            J = LEN_TRIM( HEADERS( IHDRNWLAT ) )
            PWIDTH( 3 ) = 13
            CALL ADD_TO_HEADER( PWIDTH( 3 ), HEADERS( IHDRNWLAT ),
     &                          LH, HDRBUF )
            CALL ADD_TO_HEADER( PWIDTH( 3 ), '    ', LU, UNTBUF )
        
            J = LEN_TRIM( HEADERS( IHDRNWLON ) )
            PWIDTH( 4 ) = 13
            CALL ADD_TO_HEADER( PWIDTH( 4 ), HEADERS( IHDRNWLON ),
     &                          LH, HDRBUF )
            CALL ADD_TO_HEADER( PWIDTH( 4 ), '    ', LU, UNTBUF )
        
            J = LEN_TRIM( HEADERS( IHDRNELAT ) )
            PWIDTH( 5 ) = 13
            CALL ADD_TO_HEADER( PWIDTH( 5 ), HEADERS( IHDRNELAT ),
     &                          LH, HDRBUF )
            CALL ADD_TO_HEADER( PWIDTH( 5 ), '    ', LU, UNTBUF )
        
            J = LEN_TRIM( HEADERS( IHDRNELON ) )
            PWIDTH( 6 ) = 13
            CALL ADD_TO_HEADER( PWIDTH( 6 ), HEADERS( IHDRNELON ),
     &                          LH, HDRBUF )
            CALL ADD_TO_HEADER( PWIDTH( 6 ), '    ', LU, UNTBUF )
        
            J = LEN_TRIM( HEADERS( IHDRSELAT ) )
            PWIDTH( 7 ) = 13
            CALL ADD_TO_HEADER( PWIDTH( 7 ), HEADERS( IHDRSELAT ),
     &                          LH, HDRBUF )
            CALL ADD_TO_HEADER( PWIDTH( 7 ), '    ', LU, UNTBUF )
        
            J = LEN_TRIM( HEADERS( IHDRSELON ) )
            PWIDTH( 8 ) = 13
            CALL ADD_TO_HEADER( PWIDTH( 8 ), HEADERS( IHDRSELON ),
     &                          LH, HDRBUF )
            CALL ADD_TO_HEADER( PWIDTH( 8 ), '    ', LU, UNTBUF )
            
            WRITE( LLGRDFMT, 94644 ) PWIDTH( 1 ), RPT_%DELIM,
     &                               PWIDTH( 2 ), RPT_%DELIM,
     &                               PWIDTH( 3 ), RPT_%DELIM,
     &                               PWIDTH( 4 ), RPT_%DELIM,
     &                               PWIDTH( 5 ), RPT_%DELIM,
     &                               PWIDTH( 6 ), RPT_%DELIM,
     &                               PWIDTH( 7 ), RPT_%DELIM,
     &                               PWIDTH( 8 ), RPT_%DELIM
        
            LLGRDWIDTH = SUM( PWIDTH( 1:8 ) ) + 8*LV
        
        END IF

C.........  Elevated flag column
        IF( RPT_%BYELEV ) THEN
            J = LEN_TRIM( HEADERS( IHDRELEV ) )
            J = MAX( LENELV3, J )

            CALL ADD_TO_HEADER( J, HEADERS(IHDRELEV), LH, HDRBUF )
            CALL ADD_TO_HEADER( J, ' ', LU, UNTBUF )

            ELEVWIDTH = J + LV

        END IF

C.........  Stack group IDs when BY ELEVSTAT (RPT_%BYELEV)
        IF( RPT_%ELVSTKGRP ) THEN
            J = LEN_TRIM( HEADERS( IHDRSTKG ) )
            W1 = INTEGER_COL_WIDTH( NOUTBINS, BINSTKGRP )
            W1 = MAX( W1, J )  
            CALL ADD_TO_HEADER( W1, HEADERS(IHDRSTKG), LH, HDRBUF )
            CALL ADD_TO_HEADER( W1, ' ', LU, UNTBUF )

            WRITE( STKGFMT, 94625 ) W1, RPT_%DELIM
            STKGWIDTH = W1 + LV
        END IF

C.........  Point unit type column
        IF( RPT_%BYUNIT ) THEN
            J = LEN_TRIM( HEADERS( IHDPOINT ) )
            J = MAX( CHRLEN3, J )

            CALL ADD_TO_HEADER( J, HEADERS(IHDPOINT), LH, HDRBUF )
            CALL ADD_TO_HEADER( J, ' ', LU, UNTBUF )

            POINTWIDTH = J + LV
        END IF

C.........  Emissions release point type column
        IF( RPT_%BYERPTYP ) THEN
            J = LEN_TRIM( HEADERS( IHDRERTYP ) )
            J = MAX( ERPLEN3, J )
    
            CALL ADD_TO_HEADER( J, HEADERS(IHDRERTYP), LH, HDRBUF )
            CALL ADD_TO_HEADER( J, ' ', LU, UNTBUF )

            ERTYPWIDTH = J + LV
        END IF

C.........  Plant descriptions
        IF( RPT_%SRCNAM ) THEN
            J = LEN_TRIM( HEADERS( IHDRPNAM ) )
            J = MAX( PDSCWIDTH, J )

            CALL ADD_TO_HEADER( J, HEADERS(IHDRPNAM), LH, HDRBUF )
            CALL ADD_TO_HEADER( J, ' ', LU, UNTBUF )

            PDSCWIDTH = J + LV
        END IF

C.........  ORIS descriptions
        IF( RPT_%ORISNAM ) THEN
            NWIDTH = 0
            DO I = 1, NORIS
                IF( LORISUSE( I ) ) THEN
                    L = LEN( ORISDSC( I ) )
                    NWIDTH = MAX( NWIDTH, L )
                    IF ( L .EQ. 0 ) ORISMISS = .TRUE.
                END IF
            END DO

C.............  If any missing ORIS names, check widths
            IF( ORISMISS ) NWIDTH = MAX( NWIDTH, LEN_TRIM( MISSNAME ) )

C.............  Set ORIS name column width 
            J = LEN_TRIM( HEADERS( IHDRORNM ) )
            J = MAX( NWIDTH, J ) + 2 ! two for quotes

            CALL ADD_TO_HEADER( J, HEADERS( IHDRORNM ), LH, HDRBUF )
            CALL ADD_TO_HEADER( J, ' ', LU, UNTBUF )

            ORSDSWIDTH = J + LV - 2       ! quotes in count for header print
        END IF

C.........  SCC names
        IF( RPT_%SCCNAM ) THEN

C.............  For SCC descriptions in the inventory, get max name 
C               width
            NWIDTH = 0
            DO I = 1, NINVSCC
                IF( LSCCUSE( I ) ) THEN
                    J = RPT_%SCCRES
                    L = SCCDLEV( I,J )
                    NWIDTH = MAX( NWIDTH, LEN( SCCDESC( I )( 1:L ) ) )
                    IF ( NWIDTH .EQ. 0 ) SCCMISS = .TRUE.
                END IF
            END DO

C.............  If any missing SCC names, check widths
            IF( SCCMISS ) NWIDTH = MAX( NWIDTH, LEN_TRIM( MISSNAME ) )

C.............  Set SCC name column width 
            J = LEN_TRIM( HEADERS( IHDRSNAM ) )
            J = MAX( NWIDTH, J ) + 2     ! two for quotes

            CALL ADD_TO_HEADER( J, HEADERS(IHDRSNAM), LH, HDRBUF )
            CALL ADD_TO_HEADER( J, ' ', LU, UNTBUF )

            SDSCWIDTH = J + LV - 2       ! quotes in count for header print

        END IF

C.........  SIC names
        IF( RPT_%SICNAM ) THEN

C.............  For SIC descriptions in the inventory, get max name 
C               width
            NWIDTH = 0
            DO I = 1, NINVSIC
                IF( LSICUSE( I ) ) THEN
                    NWIDTH = MAX( NWIDTH, LEN_TRIM( SICDESC( I ) ) )
                    IF ( NWIDTH .EQ. 0 ) SICMISS = .TRUE.
                END IF
            END DO

C.............  If any missing SIC names, check widths
            IF( SICMISS ) NWIDTH = MAX( NWIDTH, LEN_TRIM( MISSNAME ) )

C.............  Set SIC name column width 
            J = LEN_TRIM( HEADERS( IHDRINAM ) )
            J = MAX( NWIDTH, J ) + 2     ! two for quotes

            CALL ADD_TO_HEADER( J, HEADERS(IHDRINAM), LH, HDRBUF )
            CALL ADD_TO_HEADER( J, ' ', LU, UNTBUF )

            SIDSWIDTH = J + LV - 2       ! quotes in count for header print

        END IF

C.........  GSPRO names
        IF( RPT_%SPCNAM ) THEN

C.............  For GSPRO descriptions in the inventory, get max name
C               width
            NWIDTH = 0
            DO I = 1, NSPROF
                IF( LSPCUSE( I ) ) THEN
                    NWIDTH = MAX( NWIDTH, LEN_TRIM( SPCDESC( I ) ) )
                    IF ( NWIDTH .EQ. 0 ) SPCMISS = .TRUE.
                END IF
            END DO

C.............  If any missing GSPRO desc names, check widths
            IF( SPCMISS ) NWIDTH = MAX( NWIDTH, LEN_TRIM( MISSNAME ) )

C.............  Set GSPRO desc name column width
            J = LEN_TRIM( HEADERS( IHDRFNAM ) )
            J = MAX( NWIDTH, J ) + 2     ! two for quotes

            CALL ADD_TO_HEADER( J, HEADERS(IHDRFNAM), LH, HDRBUF )
            CALL ADD_TO_HEADER( J, ' ', LU, UNTBUF )

            SPDSWIDTH = J + LV - 2       ! quotes in count for header print

        END IF

C.........  MACT names
        IF( RPT_%MACTNAM ) THEN

C.............  For MACT descriptions in the inventory, get max name 
C               width
            NWIDTH = 0
            DO I = 1, NINVMACT
                IF( LMACTUSE( I ) ) THEN
                    NWIDTH = MAX( NWIDTH, LEN_TRIM( MACTDESC( I ) ) )
                    IF ( NWIDTH .EQ. 0 ) MACTMISS = .TRUE.
                END IF
            END DO

C.............  If any missing MACT names, check widths
            IF( MACTMISS ) NWIDTH = MAX( NWIDTH, LEN_TRIM( MISSNAME ) )

C.............  Set MACT name column width 
            J = LEN_TRIM( HEADERS( IHDRMNAM ) )
            J = MAX( NWIDTH, J ) + 2     ! two for quotes

            CALL ADD_TO_HEADER( J, HEADERS(IHDRMNAM), LH, HDRBUF )
            CALL ADD_TO_HEADER( J, ' ', LU, UNTBUF )

            MACDSWIDTH = J + LV - 2       ! quotes in count for header print

        END IF

C.........  NAICS names
        IF( RPT_%NAICSNAM ) THEN

C.............  For NAICS descriptions in the inventory, get max name 
C               width
            NWIDTH = 0
            DO I = 1, NINVNAICS
                IF( LNAICSUSE( I ) ) THEN
                    L = LEN_TRIM( NAICSDESC( I ) )
                    NWIDTH = MAX( NWIDTH, L )
                    IF ( L .EQ. 0 ) NAICSMISS = .TRUE.
                END IF
            END DO

C.............  If any missing NAICS names, check widths
            IF( NAICSMISS ) NWIDTH = MAX( NWIDTH, LEN_TRIM( MISSNAME ) )

C.............  Set NAICS name column width 
            J = LEN_TRIM( HEADERS( IHDRNNAM ) )
            J = MAX( NWIDTH, J ) + 2     ! two for quotes

            CALL ADD_TO_HEADER( J, HEADERS(IHDRNNAM), LH, HDRBUF )
            CALL ADD_TO_HEADER( J, ' ', LU, UNTBUF )

            NAIDSWIDTH = J + LV - 2       ! quotes in count for header print

        END IF

C.........  Variable column
        IF( RPT_%RPTMODE .EQ. 3 ) THEN
            NWIDTH = 0
            DO I = 1, RPT_%NUMDATA
                NWIDTH = MAX( NWIDTH, LEN_TRIM( OUTDNAM( I, RCNT ) ) )
            END DO

            J = LEN_TRIM( HEADERS( IHDRVAR ) )
            J = MAX( NWIDTH, J )

            CALL ADD_TO_HEADER( J, HEADERS(IHDRVAR), LH, HDRBUF )

            VARWIDTH = J + LV

C.............  Units
            NWIDTH = 0
            DO I = 1, RPT_%NUMDATA
                NWIDTH = MAX( NWIDTH, LEN_TRIM( OUTUNIT( I ) ) )
            END DO

            J = LEN_TRIM( HEADERS( IHDRUNIT ) )
            J = MAX( NWIDTH, J )

            CALL ADD_TO_HEADER( J, HEADERS(IHDRUNIT), LH, HDRBUF )

            UNITWIDTH = J + LV

        END IF

C.........  Determine the format type requested (if any) - either float or
C           scientific. Also determine the number of decimal places 
C           requested.
C.........  The data format will already have been QAed so not need to worry
C           about that here.
        J = INDEX( RPT_%DATAFMT, '.' )
        L = LEN_TRIM( RPT_%DATAFMT )
        NLEFT = STR2INT( RPT_%DATAFMT(   2:J-1 ) )
        NDECI = STR2INT( RPT_%DATAFMT( J+1:L   ) )

        J = INDEX( RPT_%DATAFMT, 'F' )
        DATFLOAT = ( J .GT. 0 )

C.........  Data values. Get width for columns that use the "F" format instead
C           of the "E" format.  The code will not permit the user to specify
C           a width that is too small for the value requested.
        IF( RPT_%NUMDATA .GT. 0 ) THEN

            OUTFMT = '(A,1X,'            

            IF( RPT_%NUMFILES .EQ. 1 ) THEN

                IF( RPT_%RPTMODE .EQ. 2 ) THEN

                    IF( FIRSTIME ) THEN
                        STIDX = 1
                        EDIDX = RPT_%RPTNVAR
                        FIRSTIME = .FALSE.

                    ELSE
                        IF( EDIDX + RPT_%RPTNVAR .GT.
     &                                 RPT_%NUMDATA ) THEN

                            STIDX = EDIDX + 1
                            EDIDX = RPT_%NUMDATA

                        ELSE

                            STIDX = EDIDX + 1
                            EDIDX = EDIDX + RPT_%RPTNVAR

                        END IF

                    END IF

                ELSE IF( RPT_%RPTMODE .EQ. 1 ) THEN

                    STIDX = 1
                    EDIDX = RPT_%NUMDATA

                ELSE IF( RPT_%RPTMODE .EQ. 3 ) THEN

                    STIDX = 1
                    EDIDX = 1

                ELSE

                    STIDX = 1
                    EDIDX = RPT_%NUMDATA
        
                END IF

            ELSE

                IF( FIRSTIME ) THEN
                    STIDX = 1
                    EDIDX = RPT_%RPTNVAR
                    FIRSTIME = .FALSE.

                ELSE
                    IF( EDIDX + RPT_%RPTNVAR .GT. 
     &                                 RPT_%NUMDATA ) THEN
                
                        STIDX = EDIDX + 1
                        EDIDX = RPT_%NUMDATA

                    ELSE

                        STIDX = EDIDX + 1
                        EDIDX = EDIDX + RPT_%RPTNVAR
 
                    END IF

                END IF

            END IF

            DO J = STIDX, EDIDX

C.................  Build temporary units fields and get final width
                L = LEN_TRIM( OUTUNIT( J ) )
                TMPUNIT = '[' // OUTUNIT( J )( 1:L ) // ']'
                LN = LEN_TRIM( TMPUNIT )

C.................  If float format
                IF ( DATFLOAT ) THEN

C.....................  For database output get max data value for all columns
                    IF( RPT_%RPTMODE .EQ. 3 ) THEN

                        PREVAL = 0.0
                        DO I = 1, RPT_%NUMDATA

                            VAL = MAXVAL( BINDATA( :,J ) )
                            IF( VAL .LE. PREVAL ) CYCLE
                            PREVAL = VAL

                        END DO

                        VAL = PREVAL

                    ELSE

C.....................  Get maximum data value for this column
                        VAL = MAXVAL( BINDATA( :,J ) )

                    END IF

                    BUFFER = ' '
                    WRITE( BUFFER, '(F30.0)' ) VAL
                    BUFFER = ADJUSTL( BUFFER )

C.....................  Store the minimum width of the left part of the format. 
                    W1 = LEN_TRIM( BUFFER )

C.....................  Increase the width to include the decimal places
                    W1 = W1 + NDECI + 1           ! +1 for decimal point

C.....................  Set the left part of the format.  Compare needed width 
C                       with requested width and width of the column header
C                       and units header
                    IF( RPT_%RPTMODE .EQ. 3 ) THEN
                        L2 = LEN_TRIM( HEADERS( IHDRDATA ) )
                        W1 = MAX( NLEFT, W1, L2 )
                    ELSE
                        L2 = LEN_TRIM( OUTDNAM( J,RCNT ) )
                        W1  = MAX( NLEFT, W1, L2, LN )
                    END IF

C.....................  Build the array of output formats for the data in 
C                       current report
                    TMPFMT = OUTFMT 
                    L2 = LEN_TRIM( TMPFMT )
                    WRITE( OUTFMT, '(A,I2.2,A,I2.2)' ) 
     &                     TMPFMT( 1:L2 ) // 'F', W1, '.', NDECI

C.................  If exponential output format
                ELSE

                    W1 = 0 ! added by GAP 1/17/07

C.....................  Set the left part of the format.  Compare needed width 
C                       with requested width and width of the column header
C                       and units header
                    IF( RPT_%RPTMODE .EQ. 3 ) THEN
                        L2 = LEN_TRIM( HEADERS( IHDRDATA ) )
                        W1 = MAX( NLEFT, W1, L2 )
                    ELSE
                        L2 = LEN_TRIM( OUTDNAM( J,RCNT ) )
                        W1  = MAX( NLEFT, W1, L2, LN )
                    END IF

                    L1 = LEN_TRIM( RPT_%DATAFMT )
                    TMPFMT = OUTFMT
                    L2 = LEN_TRIM( TMPFMT )
                    WRITE( OUTFMT, '(A,I2.2,A,I2.2)' ) 
     &                     TMPFMT( 1:L2 ) // 'E', W1, '.', NDECI

                END IF

C.................  Add delimeter to output formats except for last value
                TMPFMT = OUTFMT
                L1 = LEN_TRIM( TMPFMT )
                IF( J .NE. EDIDX ) THEN
                    IF( L1 .LT. QAFMTL3-8 ) THEN
                        OUTFMT = TMPFMT( 1:L1 ) // ',"' // 
     &                           RPT_%DELIM // '",1X,'
                    ELSE
                        GO TO 988
                    END IF

C.................  Otherwise make sure there is no comma on the end and
C                   add the ending parenthese
                ELSE
                    IF( L1 .LT. QAFMTL3-1 ) THEN      
                        IF( TMPFMT( L1:L1 ) .EQ. ',' ) L1 = L1 - 1
                        IF( RPT_%RPTMODE .EQ. 3 ) THEN
C                            OUTFMT = TMPFMT( 1:L1 ) // ',A,1X,A)'
                            OUTFMT = TMPFMT( 1:L1 ) // ')'
                        ELSE
                            OUTFMT = TMPFMT( 1:L1 ) // ')'
                        END IF

                    ELSE
                        GO TO 988
                    END IF

                END IF
             
C.................  Add next entry to header buffers
                IF( RPT_%RPTMODE .EQ. 3 ) THEN
                    CALL ADD_TO_HEADER( W1, HEADERS( IHDRDATA ),
     &                                  LH, HDRBUF )

                ELSE
                    CALL ADD_TO_HEADER( W1, OUTDNAM( J,RCNT ), 
     &                                  LH, HDRBUF )

                END IF

C.................  Add next entry to units line buffer
                CALL ADD_TO_HEADER( W1, TMPUNIT, LU, UNTBUF )

            END DO

        END IF     ! End if any data to output or not

C............................................................................
C.........  Write out the header to the report
C............................................................................

C.........  Write line to separate reports from each other and from metadata
        IF( PDEV == FDEV ) WRITE( FDEV, '(/,A,/)' ) REPEAT( '#', LH )

C.........  If multifile report, write out number of current file
        IF( RPT_%NUMFILES .GT. 1 ) THEN
            WRITE( MESG,94020 ) '# File', FILENUM, 'of', RPT_%NUMFILES
            WRITE( FDEV,93000 ) TRIM( MESG )

        END IF

C.........  User Titles  ....................................................

C.........  Loop through user-defined titles for current report, and write
C           to the report verbatim.
        DO I = 1, RPT_%NUMTITLE

            L2 = LEN_TRIM( TITLES( I,RCNT ) )
            WRITE( FDEV,93000 ) '# ' // TITLES( I,RCNT )( 1:L2 )

        END DO

C.........  Automatic Titles  ...............................................

C.........  Source category processed
        WRITE( FDEV,93000 ) '# Processed as ' // TRIM( CATDESC ) // 
     &                      ' sources'

C.........  The year of the inventory
        WRITE( MESG,94010 ) '# Base inventory year', BYEAR
        WRITE( FDEV,93000 ) TRIM( MESG )

        IF( PYEAR .NE. 0 ) THEN 
            WRITE( MESG,94010 ) '# Projected inventory year', PYEAR
            WRITE( FDEV,93000 ) TRIM( MESG )
        END IF

C.........  Whether projection factors were applied and for what year
        IF( RPT_%USEPRMAT ) THEN
            WRITE( MESG,94010 ) '# Projection factors applied to ' //
     &             'inventory for converting from', PRBYR, 'to', PRPYR
            WRITE( FDEV,93000 ) TRIM( MESG )
        END IF

C.........  Whether multiplicative control factors were applied
        IF( RPT_%USECUMAT ) THEN
            WRITE( FDEV,93000 ) '# Multiplicative control factors ' //
     &             'applied'
        END IF

C.........  Whether a gridding matrix was applied and the grid name
        IF( RPT_%USEGMAT .OR. AFLAG ) THEN
            WRITE( FDEV,93000 ) '# Gridding matrix applied for grid' // 
     &                          TRIM( GRDNM )
        ELSE
            WRITE( FDEV,93000 ) '# No gridding matrix applied'
        END IF

C.........  Whether a speciation matrix was applied and mole- or mass-based
        IF( RPT_%USESLMAT .OR. AFLAG ) THEN
            WRITE( FDEV,93000 ) '# Molar speciation matrix applied'

        ELSE IF( RPT_%USESSMAT ) THEN
            WRITE( FDEV,93000 ) '# Mass speciation matrix applied'

        ELSE
            WRITE( FDEV,93000 ) '# No speciation matrix applied'

        END IF

C.........  What pollutant was used for speciation profiles
        IF( RPT_%BYSPC ) THEN
            L = LEN_TRIM( RPT_%SPCPOL )
            WRITE( FDEV,93000 )'# Speciation profiles for pollutant "'//
     &                          RPT_%SPCPOL( 1:L ) // '"'
        END IF

C.........  Whether hourly data or inventory data were input 
C.........  For hourly data, the time period processed
        IF( RPT_%USEHOUR .OR. AFLAG ) THEN

            K1 = WKDAY( SDATE )
            K2 = WKDAY( EDATE )
            L1 = LEN_TRIM( DAYS( K1 ) )
            L2 = LEN_TRIM( DAYS( K2 ) )

            WRITE( FDEV,93010 ) 
     &            '# Temporal factors applied for episode from'
            WRITE( FDEV,93010 ) '# ' // BLANK5 // 
     &             DAYS( K1 )( 1:L1 ) // ' ' // MMDDYY( SDATE ) //
     &             ' at', STIME, 'to'

            WRITE( FDEV,93010 ) '# ' // BLANK5 // 
     &             DAYS( K2 )( 1:L2 ) // ' '// MMDDYY( EDATE ) //
     &             ' at', ETIME

C.............  Compare average day setting in configuration file with what
C               is actually available in the hourly emissions file.  Give
C               messages and titles accordingly.
            DYSTAT = .FALSE.
            IF( INVPIDX .EQ. 1 ) DYSTAT = .TRUE.

        ELSE
            WRITE( FDEV,93000 ) '# No temporal factors applied'

            DYSTAT = .FALSE.
            IF( RPT_%AVEDAY ) DYSTAT = .TRUE.

        END IF

C.........  Write average day status
        IF( DYSTAT ) THEN
            WRITE( FDEV,93000 ) '# Average day data basis in report'
        ELSE
            WRITE( FDEV,93000 ) '# Annual total data basis in report'
        END IF

C.........  Write normalization status
        IF( RPT_%NORMCELL ) THEN
            L = LEN_TRIM( GRDNM )
            WRITE( FDEV, 93000 ) '# Data divided by grid cell areas '//
     &                           'based on grid ' // GRDNM( 1:L )
        END IF

        IF( RPT_%NORMPOP ) THEN
            WRITE( FDEV, 93020 ) '# Data divided by year ', STCYPOPYR,
     &                           ' population'
        END IF

C.........  The name of the group used to select the data
        IF( RPT_%REGNNAM .NE. ' ' ) THEN

            L = LEN_TRIM( RPT_%REGNNAM )
            WRITE( FDEV,93000 ) '# Region group "'//RPT_%REGNNAM( 1:L )
     &                          // '" applied'

        END IF

C.........  The name of the subgrid used to select the data
        IF( RPT_%SUBGNAM .NE. ' ' ) THEN

            L = LEN_TRIM( RPT_%SUBGNAM )
            WRITE( FDEV,93000 ) '# Subgrid "' // RPT_%SUBGNAM( 1:L )
     &                          // '" applied'

        END IF

C.........  Column headers  .................................................

        IF( RPT_%RPTMODE .NE. 3 ) THEN
C.........  Remove leading spaces from column units
            L = LEN_TRIM( UNTBUF )
            UNTBUF = UNTBUF( 2:L )
            L = L - 1

C.........   Write data output units
            WRITE( FDEV, 93000 ) UNTBUF( 1:L )

        END IF

C.........  Remove leading spaces from column headers
        L = LEN_TRIM( HDRBUF )
        HDRBUF = HDRBUF( 2:L )
        L = L - 1

C.........  Write column headers
        WRITE( FDEV, 93000 ) HDRBUF( 1:L )

C.........  Store previous output file unit number
        PDEV = FDEV

C.........  Successful completion of routine
        RETURN

C.........  Unsuccessful completion of routine
988     WRITE( MESG,94010 ) 'INTERNAL ERROR: Allowable length ' //
     &         'of format statement (', QAFMTL3, ') exceeded' //
     &         CRLF() // BLANK10 // 'at output data field "'//
     &         OUTDNAM( J,RCNT )( 1:LEN_TRIM( OUTDNAM( J,RCNT ) ) ) //
     &         '". Must rerun with fewer outputs or change' // CRLF() //
     &         BLANK10 // 'value of QAFMTL3 in modreprt.f and ' //
     &         'recompile SMOKE library and Smkreport.'
       CALL M3MSG2( MESG )

       CALL M3EXIT( PROGNAME, 0, 0, ' ', 2 )

C******************  FORMAT  STATEMENTS   ******************************

C...........   Formatted file I/O formats............ 93xxx

93000   FORMAT( A )

93010   FORMAT( 10( A, :, 1X, I6.6, :, 1X ) )

93020   FORMAT( A, I4.4, A )

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I10, :, 1X ) )

94020   FORMAT( 10( A, :, I3, :, 1X ) )

94220   FORMAT( '(1X,A,",",I4.4,",",I2.2,",",I3.3,",",I1,"', A, '")' )

94620   FORMAT( '(1X,I2.2,"/",I2.2,"/",I4.4,"', A, '")' )

94625   FORMAT( '(1X,I', I2.2, ',"', A, '")' )

94630   FORMAT( '(1X,I', I2.2, '.', I2.2, ',"', A, '")' )

94635   FORMAT( '(1X,', 'I',I2.2, ',"',A,'", I',I2.2, ',"',A,'")' )

94640   FORMAT( '(', 3('1X,F', I2.2, '.5,"', A, '",'), 
     &          '1X,F', I2.2, '.5,"', A, '")' )

94642   FORMAT( '(1X,F',I2.2,'.8,"', A,'",1X,F',I2.2,'.8,"',A,'")' )  ! lat/lons

94643   FORMAT( '(', 4('1X,F', I2.2, '.2,"', A, '",'), 
     &          '1X,F', I2.2, '.2,"', A, '")' )

94644   FORMAT( '(', 8('1X,F', I2.2, '.8,"', A, '",'), ')' )

94645   FORMAT( '(I', I1, ',"', A, '")' )

94650   FORMAT( '(I', I3.3, ',"', A, '")' )

C******************  INTERNAL SUBPROGRAMS  *****************************
 
        CONTAINS
 
C.............  This internal subprogram builds the report header
            SUBROUTINE ADD_TO_HEADER( LCOL, LABEL, LHDR, HDRBUF )

C.............  Subprogram arguments
            INTEGER     , INTENT (IN)     :: LCOL   ! width of current column
            CHARACTER(*), INTENT (IN)     :: LABEL  ! column label
            INTEGER     , INTENT (IN OUT) :: LHDR   ! header length
            CHARACTER(*), INTENT (IN OUT) :: HDRBUF ! header

C----------------------------------------------------------------------

C.............  If this is the firstime for this report
            IF( LHDR .EQ. 0 ) THEN

C.................  Initialize header and its length
                HDRBUF = ' ' // '#' // LABEL
                LHDR   = LCOL + LV

C.............  If not a new report...
            ELSE

                HDRBUF = HDRBUF( 1:LHDR ) // RPT_%DELIM // ' ' // LABEL
                LHDR = LHDR + LCOL + LV     ! space included in LV

            END IF
 
            END SUBROUTINE ADD_TO_HEADER

C----------------------------------------------------------------------
C----------------------------------------------------------------------

C.............  This internal subprogram finds the width of the largest
C               integer in an array
            INTEGER FUNCTION INTEGER_COL_WIDTH( NVAL, IARRAY )

C.............  Subprogram arguments
            INTEGER, INTENT (IN) :: NVAL             ! size of array
            INTEGER, INTENT (IN) :: IARRAY( NVAL )   ! integer array

C.............  Local subprogram variables
            INTEGER          M            ! tmp integer value

            CHARACTER(16)    NUMBUF       ! tmp number string

C----------------------------------------------------------------------

C.............  Find maximum integer value in list
            M = MAXVAL( IARRAY )

C.............  Write integer to character string
            WRITE( NUMBUF, '(I16)' ) M

C.............  Find its width
            NUMBUF = ADJUSTL( NUMBUF )
            INTEGER_COL_WIDTH = LEN_TRIM( NUMBUF )
 
            END FUNCTION INTEGER_COL_WIDTH

        END SUBROUTINE WRREPHDR
