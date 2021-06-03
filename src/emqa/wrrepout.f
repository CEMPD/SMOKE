
        SUBROUTINE WRREPOUT( FDEV, RCNT, NDATA, JDATE, JTIME, 
     &                       LAYER, DELIM, OUTFMT, ZEROFLAG, EFLAG )

C***********************************************************************
C  subroutine body starts at line 
C
C  DESCRIPTION:
C    The WRREPOUT routine outputs the report information, depending on the
C    specifications of the output reports and the contents of the bins.C    
C
C  PRECONDITIONS REQUIRED:
C    The output file FDEV is opened
C    The column widths and internal write formats have been created
C    The output delimeter has been specified,
C    The report count RCNT is set
C    The bins have been populated with characteristics and data values
C    The number of data values NDATA has been set
C    The output date string has been set, if needed
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
        USE MODSOURC, ONLY: CPDESC, CSOURC, STKHT, STKDM, STKTK, STKVE,
     &                      XLOCA, YLOCA, FUGHGT, FUGWID, FUGLEN,
     &                      FUGANG, NGSPRO, GSPROID, GSPRDESC

C.........  This module contains the global variables for the 3-d grid
        USE MODGRID, ONLY:  NCOLS, NROWS, XORIG, YORIG, XOFF, YOFF,
     &                      GDTYP, XCELL, YCELL, XCENT, YCENT,
     &                      P_ALP, P_BET, P_GAM, OFFLAG, GRDNM,
     &                      XOFF_A, YOFF_A, XDIFF, YDIFF, NGRID

C.........  This module contains the lists of unique source characteristics
        USE MODLISTS, ONLY: SCCDESC, SCCDLEV, SICDESC, MACTDESC, 
     &                      NAICSDESC

C.........  This module contains Smkreport-specific settings
        USE MODREPRT, ONLY: RPT_, LREGION, VARWIDTH,
     &                      DATEFMT, DATEWIDTH, HOURFMT, HOURWIDTH,
     &                      LAYRFMT, LAYRWIDTH, CELLFMT, CELLWIDTH,
     &                      SRCFMT, SRCWIDTH, REGNFMT, REGNWIDTH,
     &                      CYWIDTH, STWIDTH, COWIDTH, SCCWIDTH,
     &                      SRG1FMT, SRG1WIDTH, SRG2FMT, SRG2WIDTH,
     &                      MONWIDTH, WEKWIDTH, DOMWIDTH, MNDWIDTH,
     &                      TUEWIDTH, WEDWIDTH, THUWIDTH, FRIWIDTH,
     &                      SATWIDTH, SUNWIDTH, METWIDTH,
     &                      CHARFMT, CHARWIDTH, SPDSWIDTH,
     &                      STKPFMT, STKPWIDTH, ELEVWIDTH,
     &                      PDSCWIDTH, SDSCWIDTH, SPCWIDTH, MINC,
     &                      LOC_BEGP, LOC_ENDP, OUTDNAM, OUTUNIT,
     &                      ALLRPT, SICWIDTH, SIDSWIDTH, UNITWIDTH,
     &                      MACTWIDTH, MACDSWIDTH, NAIWIDTH, LABELWIDTH,
     &                      NAIDSWIDTH, STYPWIDTH, LTLNFMT,
     &                      LTLNWIDTH, DLFLAG, ORSWIDTH, ORSDSWIDTH,
     &                      STKGWIDTH, STKGFMT, INTGRWIDTH, GEO1WIDTH,
     &                      ERTYPWIDTH, FUGPFMT, FUGPWIDTH, LAMBWIDTH,
     &                      LAMBFMT, LLGRDFMT, LLGRDWIDTH

C.........  This module contains report arrays for each output bin
        USE MODREPBN, ONLY: NOUTBINS, BINDATA, BINSCC, BINPLANT,
     &                      BINX, BINY, BINSMKID, BINREGN, 
     &                      BINCOIDX, BINSTIDX, BINCYIDX,
     &                      BINMONID, BINWEKID, BINDOMID, BINMNDID,
     &                      BINTUEID, BINWEDID, BINTHUID, BINFRIID,
     &                      BINSATID, BINSUNID, BINMETID, BINSPCIDX,
     &                      BINSRGID1, BINSRGID2, BINSPCID, BINRCL,
     &                      BINELEV, BINSNMIDX, BINBAD, BINSIC, 
     &                      BINSICIDX, BINMACT, BINMACIDX, BINNAICS,
     &                      BINNAIIDX, BINSRCTYP, BINORIS, BINORSIDX,
     &                      BINORIS, BINORSIDX, BINSTKGRP, BININTGR,
     &                      BINGEO1IDX, BINERPTYP, BINFACILITY

C.........  This module contains the arrays for state and county summaries
        USE MODSTCY, ONLY: CTRYNAM, STATNAM, CNTYNAM, NORIS, ORISDSC,
     &                     GEOLEV1NAM

C.........  This module contains the information about the source category
        USE MODINFO, ONLY: MXCHRS, NCHARS, BYEAR

C.........  MODULES for I/O API INTERFACEs, geo-transform codes:
        USE M3UTILIO 
        USE MODGCTP

        IMPLICIT NONE

C...........   INCLUDES
        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters

C...........   SUBROUTINE ARGUMENTS
        INTEGER     , INTENT(IN   ) :: FDEV
        INTEGER     , INTENT(IN   ) :: RCNT
        INTEGER     , INTENT(IN   ) :: NDATA
        INTEGER     , INTENT(IN   ) :: JDATE
        INTEGER     , INTENT(IN   ) :: JTIME
        INTEGER     , INTENT(IN   ) :: LAYER    ! layer number for output
        CHARACTER(*), INTENT(IN   ) :: DELIM
        CHARACTER(*), INTENT(IN   ) :: OUTFMT
        LOGICAL     , INTENT(IN   ) :: ZEROFLAG
        LOGICAL     , INTENT(INOUT) :: EFLAG

C...........   Local parameters
        INTEGER, PARAMETER :: STRLEN = 10000   ! Maximum info string length

C...........   Arrays for source characteristics output formatting
        CHARACTER(300) CHARS ( MXCHRS ) !  source fields for output

        LOGICAL, ALLOCATABLE, SAVE :: LF ( : ) ! true if column should be output

C...........   Other local variables
        INTEGER     I, J, K, L, L1, N, S, V           ! counters and indices

        INTEGER     DAY, JDAY, WDAY           ! day of month, julian day, Weekday
        INTEGER     IOS                       ! i/o status
        INTEGER     LE                        ! output string accum. length
        INTEGER     LV                        ! width of delimiter
        INTEGER     LX                        ! extra space for first column
        INTEGER     MONTH                     ! month number
        INTEGER     MXLE                      ! max value of LE
        INTEGER     NC                        ! no. source char fields
        INTEGER     OUTHOUR                   ! output hour
        INTEGER     YEAR                      ! 4-digit year
        INTEGER     STIDX                     ! starting index of loop
        INTEGER     EDIDX                     ! ending index of loop

        REAL*8      LAMBX, LAMBY, UTM_X, UTM_Y, UTM_Z ! grid lambert/utm coord for point sources
        REAL*8      DLAT, DLON, DXORIG, DYORIG    ! grid lambert/utm coord for cell center 
        REAL*8      X0(1,1), Y0(1,1)              ! grid lambert/utm coord for cell center 
        
        REAL*8      SWLAT, SWLON, NWLAT, NWLON, NELAT, NELON, SELAT, SELON ! lat-lon coords of grid cell

        INTEGER, SAVE :: PRCNT = 0

        LOGICAL, SAVE :: FIRSTIME  = .TRUE.   ! true: first time routine called

        REAL        ECHECK                    ! tmp sum of emissions in a bin
        
        CHARACTER(1)        AGT               ! tmp aggregation type for CARB QADEF report
        CHARACTER(17)       OUTCARB           !  output date string for CARB QADEF report
        CHARACTER(12)       OUTDATE           !  output date string
        CHARACTER(100)   :: BADRGNM = 'Name unknown'
        CHARACTER(200)      BUFFER            !  string building buffer
        CHARACTER(300)      MESG              !  message buffer
        CHARACTER(STRLEN)   STRING            !  output string
        CHARACTER(SCCLEN3)  TSCC              ! tmp SCC string
        CHARACTER(FIPLEN3)  TFIPS             ! tmp FIPS string

        CHARACTER(16) :: PROGNAME = 'WRREPOUT' ! program name

C***********************************************************************
C   begin body of subroutine WRREPOUT

C.........  Create hour for output
        OUTHOUR = JTIME / 10000 + 1

C.............  Width of delimeter
        LV = LEN_TRIM( DELIM )

C.........  When a new report is starting...
        IF( RCNT .NE. PRCNT ) THEN

C.............  Transfer array info to scalar info for this report
            RPT_ = ALLRPT( RCNT )  ! multi-value 

            LREGION = ( RPT_%BYGEO1 .OR. RPT_%BYCNRY .OR. RPT_%BYSTAT .OR. RPT_%BYCNTY )

C.............  Allocate memory for LF if not available already
            IF( .NOT. ALLOCATED( LF ) ) THEN
                ALLOCATE( LF( MXCHRS ), STAT=IOS )
                CALL CHECKMEM( IOS, 'LF', PROGNAME )
            END IF

C.............  Initialize output status of source characteristics
            LF = .FALSE.    ! array

C.............  Update logical source-characteristics fields
C.............  In future, there can be different cases here for "BY STACK", for
C               example
            IF( RPT_%BYSRC ) THEN
                LF( 1:NCHARS ) = .TRUE.
            END IF

        END IF

C.........  Loop through variables for the database format
        DO V = 1, NDATA

C.........  Loop through entries for all bins for current date and hour
            DO I = 1, NOUTBINS

C.............  Check for zero emissions if flag is not set
                IF( .NOT. ZEROFLAG ) THEN
                    ECHECK = SUM( BINDATA( I,1:NDATA ) )
                    IF ( RPT_%BYSPC ) ECHECK = BINDATA( I,1 )
                    IF ( ECHECK .EQ. 0. ) CYCLE
                END IF

C.............  Build tmp string based on date, hour, and other columns that
c               are included in the output file.  Whether these are included
c               is determined by the report settings.
                MXLE   = 1
                STRING = ' '
                LE     = 1
                LX     = 1

C..............  Include user-defined label in string
                IF( RPT_%USELABEL ) THEN

                    L = LABELWIDTH
                    L1 = L - LV 
                    STRING = STRING( 1:LE )// RPT_%LABEL( 1:L1 )// DELIM
                    MXLE = MXLE + L + LX + LV
                    LE = MIN( MXLE, STRLEN )
                    LX = 0

                END IF

C.............  Include date in string for CARB QADEF report
                IF( RPT_%CARB ) THEN

C.................  Temporal resolution (A:Annual, D:Daily, H:Hourly)
                    IF( RPT_%BYDATE .AND. .NOT. RPT_%BYHOUR ) THEN
                        AGT = 'D'
                    ELSE IF( RPT_%BYHOUR ) THEN
                        AGT = 'H'
                    ELSE
                        AGT = 'A'
                    ENDIF


C.................  Compute year, month, julidan day, day of week (1,,,,,,7)
                    YEAR = JDATE / 1000
                    IF( AGT == 'A' ) THEN
                        MONTH = 0
                        WDAY = 0
                        JDAY = 0
                    ELSE
                        CALL DAYMON( JDATE, MONTH, DAY )
                        WDAY = WKDAY( JDATE )
                        JDAY = JDATE - ( YEAR * 1000 )
                    ENDIF

C.................  Add date field to header
                    OUTCARB = ' '
                    WRITE( OUTCARB, DATEFMT ) AGT, YEAR, MONTH, JDAY, WDAY

                    STRING = STRING( 1:LE ) // OUTCARB
                    MXLE = MXLE + DATEWIDTH + LX
                    LE = MIN( MXLE, STRLEN )
                    LX = 0

                END IF

C.............  Include date in string
                IF( RPT_%BYDATE .AND. .NOT. RPT_%CARB ) THEN

C.................  Get month and day from Julian date
                    CALL DAYMON( JDATE, MONTH, DAY )

C.................  Compute year
                    YEAR = JDATE / 1000

C.................  Add date field to header
                    OUTDATE = ' '
                    WRITE( OUTDATE, DATEFMT ) MONTH, DAY, YEAR

                    STRING = STRING( 1:LE ) // OUTDATE
                    MXLE = MXLE + DATEWIDTH + LX
                    LE = MIN( MXLE, STRLEN )
                    LX = 0

                END IF

C.............  Include hour in string
                IF( RPT_%BYHOUR ) THEN
                    IF( .NOT. DLFLAG ) THEN
                        BUFFER = ' '
                        WRITE( BUFFER, HOURFMT ) OUTHOUR  ! Integer
                        STRING = STRING( 1:LE ) // BUFFER
                        MXLE = MXLE + HOURWIDTH + LX
                        LE = MIN( MXLE, STRLEN )
                        LX = 0
                    END IF

                ELSE IF( .NOT. RPT_%BYHOUR .AND. RPT_%CARB ) THEN
                    IF( .NOT. DLFLAG ) THEN
                        BUFFER = ' '
                        OUTHOUR = 0
                        WRITE( BUFFER, HOURFMT ) OUTHOUR  ! Integer
                        STRING = STRING( 1:LE ) // BUFFER
                        MXLE = MXLE + HOURWIDTH + LX
                        LE = MIN( MXLE, STRLEN )
                        LX = 0
                    END IF

                END IF

C.............  Include layer in string
                IF( RPT_%BYLAYER ) THEN
                    BUFFER = ' '
                    WRITE( BUFFER, LAYRFMT ) LAYER    ! Integer
                    STRING = STRING( 1:LE ) // BUFFER
                    MXLE = MXLE + LAYRWIDTH + LX
                    LE = MIN( MXLE, STRLEN )
                    LX = 0
                END IF

C.............  Include cell numbers in string
                IF( RPT_%BYCELL ) THEN
                    BUFFER = ' '
                    WRITE( BUFFER, CELLFMT ) BINX( I ), BINY( I )  ! Integers
                    STRING = STRING( 1:LE ) // BUFFER
                    MXLE = MXLE + CELLWIDTH + LX
                    LE = MIN( MXLE, STRLEN )
                    LX = 0
                ELSE IF( .NOT. RPT_%BYCELL .AND. RPT_%CARB ) THEN
                    BUFFER = ' '
                    WRITE( BUFFER, CELLFMT ) 0,0 ! Integers
                    STRING = STRING( 1:LE ) // BUFFER
                    MXLE = MXLE + CELLWIDTH + LX
                    LE = MIN( MXLE, STRLEN )
                    LX = 0
                END IF

C.............  Include source number in string
                IF( RPT_%BYSRC ) THEN
                    BUFFER = ' '
                    WRITE( BUFFER, SRCFMT ) BINSMKID( I )  ! Integer
                    STRING = STRING( 1:LE ) // BUFFER
                    MXLE = MXLE + SRCWIDTH + LX
                    LE = MIN( MXLE, STRLEN )
                    LX = 0
                END IF

C.............  Include country/state/county code in string
                IF( LREGION ) THEN

                    IF( RPT_%CARB ) THEN
                        TFIPS = BINREGN(I)(1:3)//','//BINREGN(I)(10:12)//','//BINREGN(I)(7:9)
                    ELSE
                        TFIPS = BINREGN( I )
                    ENDIF

                    STRING = STRING( 1:LE ) // TFIPS // DELIM
                    MXLE = MXLE + REGNWIDTH + LX
                    LE = MIN( MXLE, STRLEN )
                    LX = 0
                ELSE IF( .NOT. LREGION .AND. RPT_%CARB ) THEN
                    STRING = STRING( 1:LE ) // '      ' // DELIM
                    MXLE = MXLE + REGNWIDTH + LX
                    LE = MIN( MXLE, STRLEN )
                    LX = 0
                END IF

C.............  Include region level 1 name in string
                IF( RPT_%BYGEO1NAM ) THEN
                    J = BINGEO1IDX( I )
                    L = GEO1WIDTH
                    L1 = L - LV - 1                        ! 1 for space
                    IF( J .LE. 0 ) THEN
                        STRING = STRING( 1:LE ) // 
     &                           BADRGNM( 1:L1 ) // DELIM
                    ELSE
                        STRING = STRING( 1:LE ) // 
     &                           GEOLEV1NAM( J )( 1:L1 ) // DELIM
                    END IF

                    MXLE = MXLE + L
                    LE = MIN( MXLE, STRLEN )
                END IF

C.............  Include country name in string
                IF( RPT_%BYCONAM ) THEN
                    J = BINCOIDX( I )
                    L = COWIDTH
                    L1 = L - LV - 1                        ! 1 for space
                    IF( J .LE. 0 ) THEN
                        STRING = STRING( 1:LE ) // 
     &                           BADRGNM( 1:L1 ) // DELIM
                    ELSE
                        STRING = STRING( 1:LE ) // 
     &                           CTRYNAM( J )( 1:L1 ) // DELIM
                    END IF

                    MXLE = MXLE + L
                    LE = MIN( MXLE, STRLEN )
                END IF

C.............  Include state name in string
                IF( RPT_%BYSTNAM ) THEN
                    J = BINSTIDX( I )
                    L = STWIDTH
                    L1 = L - LV - 1                        ! 1 for space
                    IF( J .LE. 0 ) THEN
                        STRING = STRING( 1:LE ) // 
     &                           BADRGNM( 1:L1 ) // DELIM
                    ELSE
                        STRING = STRING( 1:LE ) // 
     &                           STATNAM( J )( 1:L1 ) // DELIM
                    END IF
                    MXLE = MXLE + L
                    LE = MIN( MXLE, STRLEN )
                END IF

C.............  Include county name in string
                IF( RPT_%BYCYNAM ) THEN
                    J = BINCYIDX( I )
                    L = CYWIDTH
                    L1 = L - LV - 1                        ! 1 for space
                    IF( J .LE. 0 ) THEN
                        STRING = STRING( 1:LE ) // 
     &                           BADRGNM( 1:L1 ) // DELIM
                    ELSE
                        STRING = STRING( 1:LE ) // 
     &                           CNTYNAM( J )( 1:L1 ) // DELIM
                    END IF
                    MXLE = MXLE + L
                    LE = MIN( MXLE, STRLEN )
                END IF

C.............  Include SCC code in string
                IF( RPT_%BYSCC ) THEN
                    L = SCCWIDTH
                    L1 = L - LV - 1                        ! 1 for space
                    TSCC = BINSCC( I )
c                    IF( TSCC(1:2) .EQ. '00' ) TSCC='  '//TSCC(3:SCCLEN3)
                    STRING = STRING( 1:LE ) // 
     &                       TSCC( 1:L1 ) // DELIM
                    MXLE = MXLE + L + LX
                    LE = MIN( MXLE, STRLEN )
                    LX = 0
                END IF

C.............  Include SIC code
                IF( RPT_%BYSIC ) THEN
                    L = SICWIDTH
                    L1 = L - LV - 1
                    STRING = STRING( 1:LE ) // 
     &                       BINSIC( I )( 1:MIN(L1,SICLEN3) ) // DELIM
                    MXLE = MXLE + L + LX
                    LE = MIN( MXLE, STRLEN )
                    LX = 0
                END IF

C.............  Include INTEGRATE code in string
                IF( RPT_%BYINTGR ) THEN
                    L = INTGRWIDTH
                    L1 = L - LV - 1                        ! 1 for space
                    STRING = STRING( 1:LE ) // 
     &                       BININTGR( I )( 1:MIN(L1,INTLEN3) ) // DELIM
                    MXLE = MXLE + L + LX
                    LE = MIN( MXLE, STRLEN )
                    LX = 0
                END IF

C.............  Include MACT code in string
                IF( RPT_%BYMACT ) THEN
                    L = MACTWIDTH
                    L1 = L - LV - 1                        ! 1 for space
                    STRING = STRING( 1:LE ) // 
     &                       BINMACT( I )( 1:MIN(L1,MACLEN3) ) // DELIM
                    MXLE = MXLE + L + LX
                    LE = MIN( MXLE, STRLEN )
                    LX = 0
                END IF

C.............  Include NAICS code in string
                IF( RPT_%BYNAICS ) THEN
                    L = NAIWIDTH
                    L1 = L - LV - 1                        ! 1 for space
                    STRING = STRING( 1:LE ) // 
     &                       BINNAICS( I )( 1:MIN(L1,NAILEN3) ) // DELIM
                    MXLE = MXLE + L + LX
                    LE = MIN( MXLE, STRLEN )
                    LX = 0
                END IF

C.............  Include SRCTYP code in string
                IF( RPT_%BYSRCTYP ) THEN
                    L = STYPWIDTH
                    L1 = L - LV - 1                        ! 1 for space
                    STRING = STRING( 1:LE+8 ) // 
     &                    BINSRCTYP( I )( 1:L1-8 ) // DELIM
                    MXLE = MXLE + L + LX
                    LE = MIN( MXLE, STRLEN )
                    LX = 0
                END IF

C.............  Include primary surrogate code
                IF( RPT_%SRGRES .EQ. 1 ) THEN
                    BUFFER = ' '
                    WRITE( BUFFER, SRG1FMT ) BINSRGID1( I )  ! Integer
                    STRING = STRING( 1:LE ) // BUFFER
                    MXLE = MXLE + SRG1WIDTH + LX
                    LE = MIN( MXLE, STRLEN )
                    LX = 0
                END IF

C.............  Include fallback surrogate code
                IF( RPT_%SRGRES .GE. 1 ) THEN
                    BUFFER = ' '
                    WRITE( BUFFER, SRG2FMT ) BINSRGID2( I )  ! Integer
                    STRING = STRING( 1:LE ) // BUFFER
                    MXLE = MXLE + SRG2WIDTH + LX
                    LE = MIN( MXLE, STRLEN )
                    LX = 0
                END IF

C.............  Include monthly temporal profile
                IF( RPT_%BYMON ) THEN
                    L = MONWIDTH
                    L1 = L - LV - 1 - 1                  ! 1 for space                
                    STRING = STRING( 1:LE ) // ' ' //
     &                       BINMONID( I )( 1:L1 )// DELIM
                    MXLE = MXLE + L + LX
                    LE = MIN( MXLE, STRLEN )
                    LX = 0
                END IF

C.............  Include weekly temporal profile
                IF( RPT_%BYWEK ) THEN
                    L = WEKWIDTH
                    L1 = L - LV - 1 - 1                  ! 1 for space                
                    STRING = STRING( 1:LE ) // ' ' //
     &                       BINWEKID( I )( 1:L1 )// DELIM
                    MXLE = MXLE + L + LX
                    LE = MIN( MXLE, STRLEN )
                    LX = 0
                END IF

C.............  Include day of month temporal profile
                IF( RPT_%BYDOM ) THEN
                    L = DOMWIDTH
                    L1 = L - LV - 1 - 1                  ! 1 for space                
                    STRING = STRING( 1:LE ) // ' ' //
     &                       BINDOMID( I )( 1:L1 )// DELIM
                    MXLE = MXLE + L + LX
                    LE = MIN( MXLE, STRLEN )
                    LX = 0
                END IF

C.............  Include Monday diurnal temporal profile
                IF( RPT_%BYMND ) THEN
                    L = MNDWIDTH
                    L1 = L - LV - 1 - 1                  ! 1 for space                
                    STRING = STRING( 1:LE ) // ' ' //
     &                       BINMNDID( I )( 1:L1 )// DELIM
                    MXLE = MXLE + L + LX
                    LE = MIN( MXLE, STRLEN )
                    LX = 0
                END IF

C.............  Include Tuesday diurnal temporal profile
                IF( RPT_%BYTUE ) THEN
                    L = TUEWIDTH
                    L1 = L - LV - 1 - 1                  ! 1 for space                
                    STRING = STRING( 1:LE ) // ' ' //
     &                       BINTUEID( I )( 1:L1 )// DELIM
                    MXLE = MXLE + L + LX
                    LE = MIN( MXLE, STRLEN )
                    LX = 0
                END IF

C.............  Include Wednesday diurnal temporal profile
                IF( RPT_%BYWED ) THEN
                    L = WEDWIDTH
                    L1 = L - LV - 1 - 1                  ! 1 for space                
                    STRING = STRING( 1:LE ) // ' ' //
     &                       BINWEDID( I )( 1:L1 )// DELIM
                    MXLE = MXLE + L + LX
                    LE = MIN( MXLE, STRLEN )
                    LX = 0
                END IF

C.............  Include Thursday diurnal temporal profile
                IF( RPT_%BYTHU ) THEN
                    L = THUWIDTH
                    L1 = L - LV - 1 - 1                  ! 1 for space                
                    STRING = STRING( 1:LE ) // ' ' //
     &                       BINTHUID( I )( 1:L1 )// DELIM
                    MXLE = MXLE + L + LX
                    LE = MIN( MXLE, STRLEN )
                    LX = 0
                END IF

C.............  Include Friday diurnal temporal profile
                IF( RPT_%BYFRI ) THEN
                    L = FRIWIDTH
                    L1 = L - LV - 1 - 1                  ! 1 for space                
                    STRING = STRING( 1:LE ) // ' ' //
     &                       BINFRIID( I )( 1:L1 )// DELIM
                    MXLE = MXLE + L + LX
                    LE = MIN( MXLE, STRLEN )
                    LX = 0
                END IF

C.............  Include Saturday diurnal temporal profile
                IF( RPT_%BYSAT ) THEN
                    L = SATWIDTH
                    L1 = L - LV - 1 - 1                  ! 1 for space                
                    STRING = STRING( 1:LE ) // ' ' //
     &                       BINSATID( I )( 1:L1 )// DELIM
                    MXLE = MXLE + L + LX
                    LE = MIN( MXLE, STRLEN )
                    LX = 0
                END IF

C.............  Include Sunday diurnal temporal profile
                IF( RPT_%BYSUN ) THEN
                    L = SUNWIDTH
                    L1 = L - LV - 1 - 1                  ! 1 for space                
                    STRING = STRING( 1:LE ) // ' ' //
     &                       BINSUNID( I )( 1:L1 )// DELIM
                    MXLE = MXLE + L + LX
                    LE = MIN( MXLE, STRLEN )
                    LX = 0
                END IF

C.............  Include Met-based diurnal temporal profile
                IF( RPT_%BYMET ) THEN
                    L = METWIDTH
                    L1 = L - LV - 1 - 1                  ! 1 for space                
                    STRING = STRING( 1:LE ) // ' ' //
     &                       BINMETID( I )( 1:L1 )// DELIM
                    MXLE = MXLE + L + LX
                    LE = MIN( MXLE, STRLEN )
                    LX = 0
                END IF

C.............  Include speciation profile
                IF( RPT_%BYSPC ) THEN
                    L = SPCWIDTH
                    L1 = L - LV - 1 - SPNLEN3                  ! 1 for space                
                    STRING = STRING( 1:LE ) //
     &                       BINSPCID( I )// BLANK16( 1:L1 )// DELIM
                    MXLE = MXLE + L + LX
                    LE = MIN( MXLE, STRLEN )
                    LX = 0
                END IF

C.............  Include plant ID
                IF( RPT_%BYPLANT ) THEN
                    L = CHARWIDTH
                    L1 = L - LV - 1                        ! 1 for space
                    STRING = STRING( 1:LE ) //
     &                       BINPLANT( I )( 1:L1 ) // DELIM
                    MXLE = MXLE + L + LX
                    LE = MIN( MXLE, STRLEN )
                    LX = 0
                END IF

C.............  Include Facility ID
                IF( RPT_%BYFACILITY ) THEN
                    L = CHARWIDTH
                    L1 = L - LV - 1                        ! 1 for space
                    STRING = STRING( 1:LE ) //
     &                       BINFACILITY( I )( 1:L1 ) // DELIM
                    MXLE = MXLE + L + LX
                    LE = MIN( MXLE, STRLEN )
                    LX = 0
                END IF


C.............  Include road class code
                IF( RPT_%BYRCL ) THEN

C.................  Write characteristics
                    BUFFER = ' '
                    WRITE( BUFFER, CHARFMT ) BINRCL( I )
                    STRING = STRING( 1:LE ) // BUFFER
                    MXLE = MXLE + CHARWIDTH + LX
                    LE = MIN( MXLE, STRLEN )
                    LX = 0
                END IF

C.............  Include source characteristics
                IF( RPT_%BYSRC ) THEN
                    S = BINSMKID( I ) 

C.................  Disaggregate source characteristics
                    CALL PARSCSRC( CSOURC( S ), MXCHRS, LOC_BEGP,
     &                             LOC_ENDP, LF, NC, CHARS )

C.................  Write characteristics
                    BUFFER = ' '
                    IF( MINC < NC ) THEN    ! only for mobile, point 
                        WRITE( BUFFER, CHARFMT )
     &                                ( CHARS( K ), K = MINC, NC )
                    END IF
                    STRING = STRING( 1:LE ) // BUFFER
                    MXLE = MXLE + CHARWIDTH + LX
                    LE = MIN( MXLE, STRLEN )
                    LX = 0
                END IF

C.............  Include ORIS code in string
                IF( RPT_%BYORIS ) THEN
                    L = ORSWIDTH
                    L1 = L - LV - 1                        ! 1 for space
                    STRING = STRING( 1:LE ) // 
     &                       BINORIS( I )( 1:MIN(L1,ORSLEN3) ) // DELIM
                    MXLE = MXLE + L + LX
                    LE = MIN( MXLE, STRLEN )
                    LX = 0
                END IF

C.............  Include stack parameters
                IF( RPT_%STKPARM ) THEN
                    S = BINSMKID( I ) 
                    BUFFER = ' '
                    WRITE( BUFFER, STKPFMT ) STKHT( S ), STKDM( S ),
     &                                       STKTK( S ), STKVE( S )
                    STRING = STRING( 1:LE ) // BUFFER
                    MXLE = MXLE + STKPWIDTH
                    LE = MIN( MXLE, STRLEN )
                END IF

C.............  Include fugitive parameters
                IF( RPT_%FUGPARM ) THEN
                    S = BINSMKID( I )
                    BUFFER = ' '
                    WRITE( BUFFER, FUGPFMT ) FUGHGT( S ), FUGWID( S ),
     &                                       FUGLEN( S ), FUGANG( S )
                    STRING = STRING( 1:LE ) // BUFFER
                    MXLE = MXLE + FUGPWIDTH
                    LE = MIN( MXLE, STRLEN )
                END IF

C.............  Include lat/lons for point sources
                IF( RPT_%LATLON ) THEN
                    S = BINSMKID( I )
                    BUFFER = ' '
                    WRITE( BUFFER, LTLNFMT ) YLOCA( S ), XLOCA( S )
                    STRING = STRING( 1:LE ) // BUFFER
                    MXLE = MXLE + LTLNWIDTH
                    LE = MIN( MXLE, STRLEN )
                END IF

C.............  Include grid lambert/utm coordinates for point source
                IF( RPT_%GRDCOR ) THEN

                    S = BINSMKID( I )
                    LAMBX = 0.0D0
                    LAMBY = 0.0D0
                    UTM_X = 0.0D0
                    UTM_Y = 0.0D0
                    UTM_Z = 0.0D0
                    DLON  = XLOCA( S )
                    DLAT  = YLOCA( S )

C.....................  Compute lamber x&Y and utm x&y coordinates for point source
                    CALL XY2XY( GDTYP  , P_ALP, P_BET, P_GAM, XCENT, YCENT,
     &                          LATGRD3, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0,
     &                          DLON, DLAT, LAMBX, LAMBY )

C.....................  Compute the LL for grid cell center UTM zone
                    DXORIG = XORIG + XCELL * ( BINX( I ) - 1 )
                    DYORIG = YORIG + YCELL * ( BINY( I ) - 1 )
                    CALL GRID2XY( LATGRD3, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0,
     &                            GDTYP, P_ALP, P_BET, P_GAM, XCENT, YCENT,
     &                            1, 1, DXORIG, DYORIG, XCELL, YCELL,
     &                            X0, Y0 )
                    UTM_Z = DBLE( FLOOR( ( ( X0(1,1) + 180.0 ) / 6.0 ) + 1.0 ) )

                    CALL XY2XY( UTMGRD3, UTM_Z, 0.0D0, 0.0D0, 0.0D0, 0.0D0,
     &                          LATGRD3, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0,
     &                          DLON, DLAT, UTM_X, UTM_Y )

                    BUFFER = ' '
                    WRITE( BUFFER, LAMBFMT ) LAMBX, LAMBY, UTM_X, UTM_Y, UTM_Z 
                    STRING = STRING( 1:LE ) // BUFFER
                    MXLE = MXLE + LAMBWIDTH
                    LE = MIN( MXLE, STRLEN )

                END IF

C.............  Include grid cell corner coordinates
                IF( RPT_%GRDPNT ) THEN
                
                    DXORIG = XORIG + XCELL * ( BINX( I ) - 1 )  ! SW corner
                    DYORIG = YORIG + YCELL * ( BINY( I ) - 1 )
                    
                    CALL XY2XY( LATGRD3, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0,
     &                          GDTYP  , P_ALP, P_BET, P_GAM, XCENT, YCENT,
     &                          DXORIG, DYORIG, SWLON, SWLAT )
                    
                    DYORIG = DYORIG + YCELL  ! NW corner
                    
                    CALL XY2XY( LATGRD3, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0,
     &                          GDTYP  , P_ALP, P_BET, P_GAM, XCENT, YCENT,
     &                          DXORIG, DYORIG, NWLON, NWLAT )
                    
                    DXORIG = DXORIG + XCELL  ! NE corner
                    
                    CALL XY2XY( LATGRD3, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0,
     &                          GDTYP  , P_ALP, P_BET, P_GAM, XCENT, YCENT,
     &                          DXORIG, DYORIG, NELON, NELAT )
                    
                    DYORIG = DYORIG - YCELL  ! SE corner
                    
                    CALL XY2XY( LATGRD3, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0,
     &                          GDTYP  , P_ALP, P_BET, P_GAM, XCENT, YCENT,
     &                          DXORIG, DYORIG, SELON, SELAT )
                    
                    BUFFER = ' '
                    WRITE( BUFFER, LLGRDFMT ) SWLAT, SWLON, NWLAT, NWLON, NELAT, NELON, SELAT, SELON
                    STRING = STRING( 1:LE ) // BUFFER
                    MXLE = MXLE + LLGRDWIDTH
                    LE = MIN( MXLE, STRLEN )
                
                END IF

C.............  Include elevated sources flag
                IF( RPT_%BYELEV ) THEN
                    L = ELEVWIDTH
                    L1 = L - LV  - 2      ! 1 for space, minus 1 for how used
                    STRING = STRING( 1:LE ) // 
     &                       BLANK16( 1:L1 ) // BINELEV( I ) // DELIM
                    MXLE = MXLE + L + LX
                    LE = MIN( MXLE, STRLEN )
                    LX = 0
                END IF

C.............  Include stack group IDs when elevated sources flag = Y
                IF( RPT_%ELVSTKGRP ) THEN
                    BUFFER = ' '
                    WRITE( BUFFER, STKGFMT ) BINSTKGRP( I )    ! Integer
                    STRING = STRING( 1:LE ) // BUFFER
                    MXLE = MXLE + STKGWIDTH + LX
                    LE = MIN( MXLE, STRLEN )
                    LX = 0

                END IF

C.............  Include emissions release point type
                IF( RPT_%BYERPTYP ) THEN
                    L = ERTYPWIDTH
                    L1 = L - LV - 1                        ! 1 for space
                    STRING = STRING( 1:LE+15 ) //
     &                       BINERPTYP( I )( 1:L1-15 ) // DELIM
                    MXLE = MXLE + L + LX
                    LE = MIN( MXLE, STRLEN )
                    LX = 0
                END IF

C.............  Include plant description (for point sources)
                IF( RPT_%SRCNAM ) THEN
                    S = BINSMKID( I )
                    L = PDSCWIDTH
                    L1 = L - LV - 1                        ! 1 for space
                    STRING = STRING( 1:LE ) // 
     &                       CPDESC( S )( 1:L1 ) // DELIM
                    MXLE = MXLE + L
                    LE = MIN( MXLE, STRLEN )
                END IF

C.............  Include ORIS description
C.............  This is knowingly including extra blanks before final quote
                IF( RPT_%ORISNAM ) THEN
                    J = BINORSIDX( I ) 
                    L = ORSDSWIDTH
                    L1 = L - LV - 1                        ! 1 for space
                    IF ( J > 0 ) THEN
                        BUFFER = '"'//ORISDSC( J )(1:L1)//'"'
                    ELSE
                        BUFFER = ' '  ! Leave field blank without quotes
                    ENDIF

                    STRING = STRING( 1:LE ) // BUFFER(1:L1+2) // DELIM

                    MXLE = MXLE + L + 2
                    LE = MIN( MXLE, STRLEN )
                END IF

C.............  Include SCC description
C.............  This is knowingly including extra blanks before final quote
                IF( RPT_%SCCNAM ) THEN
                    J = BINSNMIDX( I ) 
                    L = SDSCWIDTH
                    L1 = L - LV - 1                        ! 1 for space
                    K = SCCDLEV( J, RPT_%SCCRES )
                    STRING = STRING( 1:LE ) // 
     &                       SCCDESC( J )( 1:K ) //
     &                       REPEAT( " ", L1-K ) // DELIM
                    MXLE = MXLE + L + 2
                    LE = MIN( MXLE, STRLEN )

C.....................  Write warning msg when the description is unavailable
                    N = INDEX( SCCDESC( J ), 'Description unavailable' )
                    MESG = 'WARNING: Description of SCC ' // 
     &                      TRIM(BINSCC( I )) // ' is not available'
                    IF ( N .GT. 0 ) CALL M3MESG( MESG )

                END IF

C.............  Include SIC description
C.............  This is knowingly including extra blanks before final quote
                IF( RPT_%SICNAM ) THEN
                    J = BINSICIDX( I ) 
                    L = SIDSWIDTH
                    L1 = L - LV - 1                        ! 1 for space
                    STRING = STRING( 1:LE ) // 
     &                       '"'// SICDESC( J )( 1:L1 )// '"' // DELIM
                    MXLE = MXLE + L + 2
                    LE = MIN( MXLE, STRLEN )

C.....................  Write warning msg when the description is unavailable
                    N = INDEX( SICDESC( J ), 'Description unavailable' )
                    MESG = 'WARNING: Description of SIC ' // 
     &                      BINSIC( I ) // ' is not available'
                    IF ( N .GT. 0 ) CALL M3MESG( MESG )

                END IF

C.............  Include GSPRO description
C.............  This is knowingly including extra blanks before final quote
                IF( RPT_%GSPRONAM ) THEN
                    J = BINSPCIDX( I )
                    L = SPDSWIDTH
                    L1 = L - LV - 1                        ! 1 for space
                    STRING = STRING( 1:LE ) //
     &                       '"'// GSPRDESC( J )( 1:L1 )// '"' // DELIM
                    MXLE = MXLE + L + 2
                    LE = MIN( MXLE, STRLEN )

C.....................  Write warning msg when the description is unavailable
                    N = INDEX( GSPRDESC( J ), 'Description unavailable' )
                    MESG = 'WARNING: Description of GSPRO ' //
     &                      BINSPCID( I ) // ' is not available'
                    IF ( N .GT. 0 ) CALL M3MESG( MESG )

                END IF

C.............  Include MACT description
C.............  This is knowingly including extra blanks before final quote
                IF( RPT_%MACTNAM ) THEN
                    J = BINMACIDX( I ) 
                    L = MACDSWIDTH
                    L1 = L - LV - 1                        ! 1 for space
                    STRING = STRING( 1:LE ) // 
     &                       '"'// MACTDESC( J )( 1:L1 )// '"' // DELIM
                    MXLE = MXLE + L + 2
                    LE = MIN( MXLE, STRLEN )

C.....................  Write warning msg when the description is unavailable
                    N = INDEX( MACTDESC( J ), 'Description unavailable')
                    MESG = 'WARNING: Description of MACT ' // 
     &                      TRIM(BINMACT( I )) // ' is not available'
                    IF ( N .GT. 0 ) CALL M3MESG( MESG )

                END IF

C.............  Include NAICS description
C.............  This is knowingly including extra blanks before final quote
                IF( RPT_%NAICSNAM ) THEN
                    J = BINNAIIDX( I ) 
                    L = NAIDSWIDTH
                    L1 = L - LV - 1                        ! 1 for space
                    STRING = STRING( 1:LE ) // 
     &                       '"'// NAICSDESC( J )( 1:L1 )// '"' // DELIM
                    MXLE = MXLE + L + 2
                    LE = MIN( MXLE, STRLEN )

C.....................  Write warning msg when the description is unavailable
                    N = INDEX( NAICSDESC( J ),'Description unavailable')
                    MESG = 'WARNING: Description of NAICS ' // 
     &                      TRIM(BINNAICS( I )) // ' is not available'
                    IF ( N .GT. 0 ) CALL M3MESG( MESG )

                END IF

C.............  Include variable in string
                IF( RPT_%RPTMODE .EQ. 3 ) THEN

                    L = VARWIDTH
                    L1 = L - LV
                    STRING = STRING( 1:LE ) //
     &                       OUTDNAM( V, RCNT )( 1:L1 ) // DELIM
                    MXLE = MXLE + L
                    LE = MIN( MXLE, STRLEN )

                    L = UNITWIDTH
                    L1 = L - LV
                    STRING = STRING( 1:LE ) //
     &                       OUTUNIT( V )( 1:L1 ) // DELIM
                    MXLE = MXLE + L
                    LE = MIN( MXLE, STRLEN )

                END IF

C.............  Remove leading spaces and get new length
                STRING = STRING( 2:LE )
                LE = LE - 1

C.............  Output error message of string is getting shortened
                IF( MXLE .GT. STRLEN ) THEN
                    EFLAG = .TRUE.
                    WRITE( MESG,94010 ) 'INTERNAL ERROR: Output ' //
     &                     'string getting truncated in report', RCNT
                    CALL M3MSG2( MESG )
                END IF

C.............  If current bin has bad values, update those now.
                IF( BINBAD( I ) .GT. 0 ) BINDATA( I,1:NDATA ) = -9.

C.............  Write out this record
                IF( RPT_%NUMFILES .EQ. 1 ) THEN

                    IF( RPT_%RPTMODE .EQ. 2 ) THEN

                        IF( FIRSTIME ) THEN
                            STIDX = 1
                            EDIDX = RPT_%RPTNVAR
                            FIRSTIME = .FALSE.

                        ELSE
                            IF( EDIDX + RPT_%RPTNVAR .GT.
     &                                     NDATA ) THEN

                                STIDX = EDIDX + 1
                                EDIDX = NDATA

                            ELSE
                                STIDX = EDIDX + 1
                                EDIDX = EDIDX + RPT_%RPTNVAR

                            END IF

                        END IF

                    ELSE IF( RPT_%RPTMODE .EQ. 1 ) THEN

                        STIDX = 1
                        EDIDX = NDATA

                    ELSE IF( RPT_%RPTMODE .EQ. 3 ) THEN

                        STIDX = 1
                        EDIDX = 1

                    ELSE

                        STIDX = 1
                        EDIDX = NDATA

                    END IF

                ELSE
                    IF( FIRSTIME ) THEN
                        STIDX = 1
                        EDIDX = RPT_%RPTNVAR
                        FIRSTIME = .FALSE.

                    ELSE
                        IF( EDIDX + RPT_%RPTNVAR .GT.
     &                                     NDATA ) THEN

                            STIDX = EDIDX + 1
                            EDIDX = NDATA

                        ELSE
                            STIDX = EDIDX + 1
                            EDIDX = EDIDX + RPT_%RPTNVAR

                        END IF

                    END IF

                END IF

                IF( RPT_%RPTMODE .NE. 3 ) THEN

                    WRITE( FDEV, OUTFMT ) STRING( 1:LE ), 
     &                              ( BINDATA( I,J ), J=STIDX, EDIDX )

                ELSE

                    WRITE( FDEV, OUTFMT ) STRING( 1:LE ), BINDATA(I,V)

                END IF

            END DO  ! End loop through bins

            IF( RPT_%RPTMODE .NE. 3 ) GOTO 777

        END DO   ! End loop through variables

C.........  Save report number for next time routine is called
777     PRCNT = RCNT

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Formatted file I/O formats............ 93xxx

93000   FORMAT( A )

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I10, :, 1X ) )

94180   FORMAT( I2.2, '/', I2.2, '/', I4.4 )

        END SUBROUTINE WRREPOUT







