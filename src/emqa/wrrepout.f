
        SUBROUTINE WRREPOUT( FDEV, RCNT, NDATA, JDATE, JTIME, 
     &                       LAYER, DELIM, OUTFMT, EFLAG )

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
C
C***********************************************************************
C  
C Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
C                System
C File: @(#)$Id$
C  
C COPYRIGHT (C) 2000, MCNC--North Carolina Supercomputing Center
C All Rights Reserved
C  
C See file COPYRIGHT for conditions of use.
C  
C Environmental Programs Group
C MCNC--North Carolina Supercomputing Center
C P.O. Box 12889
C Research Triangle Park, NC  27709-2889
C  
C env_progs@mcnc.org
C  
C Pathname: $Source$
C Last updated: $Date$ 
C  
C***********************************************************************

C.........  MODULES for public variables
C...........   This module is the inventory arrays
        USE MODSOURC

C.........  This module contains the lists of unique source characteristics
        USE MODLISTS

C.........  This module contains Smkreport-specific settings
        USE MODREPRT

C.........  This module contains report arrays for each output bin
        USE MODREPBN

C.........  This module contains the arrays for state and county summaries
        USE MODSTCY

C.........  This module contains the information about the source category
        USE MODINFO

        IMPLICIT NONE

C...........   INCLUDES
        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters

C...........  EXTERNAL FUNCTIONS and their descriptions:
        CHARACTER*2   CRLF
        EXTERNAL   CRLF

C...........   SUBROUTINE ARGUMENTS
        INTEGER     , INTENT (IN) :: FDEV
        INTEGER     , INTENT (IN) :: RCNT
        INTEGER     , INTENT (IN) :: NDATA
        INTEGER     , INTENT (IN) :: JDATE
        INTEGER     , INTENT (IN) :: JTIME
        INTEGER     , INTENT (IN) :: LAYER    ! layer number for output
        CHARACTER(*), INTENT (IN) :: DELIM
        CHARACTER(*), INTENT (IN) :: OUTFMT
        LOGICAL     , INTENT(OUT) :: EFLAG

C...........   Local parameters
        INTEGER, PARAMETER :: STRLEN = 2000   ! Maximum info string length

C...........   Arrays for source characteristics output formatting
        CHARACTER*300 CHARS ( MXCHRS ) !  source fields for output

        LOGICAL, ALLOCATABLE, SAVE :: LF ( : ) ! true if column should be output

C...........   Other local variables
        INTEGER     I, J, K, L, L1, S              ! counters and indices

        INTEGER     DAY                       ! day of month
        INTEGER     IOS                       ! i/o status
        INTEGER     LE                        ! output string accum. length
        INTEGER     LV                        ! width of delimiter
        INTEGER     LX                        ! extra space for first column
        INTEGER     MONTH                     ! month number
        INTEGER     MXLE                      ! max value of LE
        INTEGER     NC                        ! no. source char fields
        INTEGER     OUTHOUR                   ! output hour
        INTEGER     YEAR                      ! 4-digit year

        INTEGER, SAVE :: PRCNT = 0

        LOGICAL, SAVE :: FIRSTIME  = .TRUE.   ! true: first time routine called

        REAL        ECHECK                    ! tmp sum of emissions in a bin

        CHARACTER*12            OUTDATE           !  output date string
        CHARACTER*100        :: BADRGNM = 'Name unknown'
        CHARACTER*100           BUFFER            !  string building buffer
        CHARACTER*300           MESG              !  message buffer
        CHARACTER(LEN=STRLEN)   STRING            !  output string

        CHARACTER*16 :: PROGNAME = 'WRREPOUT' ! program name

C***********************************************************************
C   begin body of subroutine WRREPOUT

C.........  Create hour for output
        OUTHOUR = JTIME / 10000 + 1

C.........  When a new report is starting...
        IF( RCNT .NE. PRCNT ) THEN

C.............  Transfer array info to scalar info for this report
            RPT_ = ALLRPT( RCNT )  ! multi-value 

            LREGION = ( RPT_%BYCNRY .OR. RPT_%BYSTAT .OR. RPT_%BYCNTY )

C.............  Allocate memory for LF if not available already
            IF( .NOT. ALLOCATED( LF ) ) THEN
        	ALLOCATE( LF( MXCHRS ), STAT=IOS )
        	CALL CHECKMEM( IOS, 'LF', PROGNAME )
            END IF

C.............  Initialize output status of source characteristics
            LF = .FALSE.    ! array

C.............  Width of delimeter
            LV = LEN_TRIM( DELIM )

C.............  Update logical source-characteristics fields
C.............  In future, there can be different cases here for "BY STACK", for
C               example
            IF( RPT_%BYSRC ) THEN
                LF( 1:NCHARS ) = .TRUE.
            END IF

        END IF

C.........  Loop through entries for all bins for current date and hour
        DO I = 1, NOUTBINS

C.............  Check if emissions are zero and skip record if they are.
            ECHECK = SUM( BINDATA( I,1:NDATA ) )
            IF ( ECHECK .EQ. 0. ) CYCLE

C.............  Build tmp string based on date, hour, and other columns that
c               are included in the output file.  Whether these are included
c               is determined by the report settings.
            MXLE   = 1
            STRING = ' '
            LE     = 1
            LX     = 1

C.............  Include date in string
            IF( RPT_%BYDATE ) THEN

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
                BUFFER = ' '
                WRITE( BUFFER, HOURFMT ) OUTHOUR  ! Integer
                STRING = STRING( 1:LE ) // BUFFER
                MXLE = MXLE + HOURWIDTH + LX
                LE = MIN( MXLE, STRLEN )
                LX = 0
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
                BUFFER = ' '
                WRITE( BUFFER, REGNFMT ) BINREGN( I )  ! Integer
                STRING = STRING( 1:LE ) // BUFFER
                MXLE = MXLE + REGNWIDTH + LX
                LE = MIN( MXLE, STRLEN )
                LX = 0
            END IF


C.............  Include country name in string
            IF( RPT_%BYCONAM ) THEN
                J = BINCOIDX( I )
                L = COWIDTH
                L1 = L - LV - 1                        ! 1 for space
                IF( J .LE. 0 ) THEN
                    STRING = STRING( 1:LE ) // 
     &                       BADRGNM( 1:L1 ) // DELIM
                ELSE
                    STRING = STRING( 1:LE ) // 
     &                       CTRYNAM( J )( 1:L1 ) // DELIM
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
     &                       BADRGNM( 1:L1 ) // DELIM
                ELSE
                    STRING = STRING( 1:LE ) // 
     &                       STATNAM( J )( 1:L1 ) // DELIM
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
     &                       BADRGNM( 1:L1 ) // DELIM
                ELSE
                    STRING = STRING( 1:LE ) // 
     &                       CNTYNAM( J )( 1:L1 ) // DELIM
                END IF
                MXLE = MXLE + L
                LE = MIN( MXLE, STRLEN )
            END IF

C.............  Include SCC code in string
            IF( RPT_%BYSCC ) THEN
                L = SCCWIDTH
                L1 = L - LV - 1                        ! 1 for space
                STRING = STRING( 1:LE ) // 
     &                   ' ' // BINSCC( I )( 1:L1 ) // DELIM
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
                BUFFER = ' '
                WRITE( BUFFER, MONFMT ) BINMONID( I )  ! Integer
                STRING = STRING( 1:LE ) // BUFFER
                MXLE = MXLE + MONWIDTH + LX
                LE = MIN( MXLE, STRLEN )
                LX = 0
            END IF

C.............  Include weekly temporal profile
            IF( RPT_%BYWEK ) THEN
                BUFFER = ' '
                WRITE( BUFFER, WEKFMT ) BINWEKID( I )  ! Integer
                STRING = STRING( 1:LE ) // BUFFER
                MXLE = MXLE + WEKWIDTH + LX
                LE = MIN( MXLE, STRLEN )
                LX = 0
            END IF

C.............  Include diurnal temporal profile
            IF( RPT_%BYDIU ) THEN
                BUFFER = ' '
                WRITE( BUFFER, DIUFMT ) BINDIUID( I )  ! Integer
                STRING = STRING( 1:LE ) // BUFFER
                MXLE = MXLE + DIUWIDTH + LX
                LE = MIN( MXLE, STRLEN )
                LX = 0
            END IF

C.............  Include speciation profile
            IF( RPT_%BYSPC ) THEN
                L = SPCWIDTH
                L1 = L - LV - 1 - SPNLEN3                  ! 1 for space                
                STRING = STRING( 1:LE ) // ' ' //
     &                   BINSPCID( I )// BLANK16( 1:L1 )// DELIM
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
                CALL PARSCSRC( CSOURC( S ), MXCHRS, LOC_BEGP, LOC_ENDP, 
     &                         LF, NC, CHARS )

C.................  Write characteristics
                BUFFER = ' '
                WRITE( BUFFER, CHARFMT ) ( CHARS( K ), K = MINC, NC )
                STRING = STRING( 1:LE ) // BUFFER
                MXLE = MXLE + CHARWIDTH + LX
                LE = MIN( MXLE, STRLEN )
                LX = 0
            END IF

C.............  Include stack parameters
            IF( RPT_%STKPARM ) THEN
                S = BINSMKID( I ) 
                BUFFER = ' '
                WRITE( BUFFER, STKPFMT ) STKHT( S ), STKDM( S ),
     &                                   STKTK( S ), STKVE( S )
                STRING = STRING( 1:LE ) // BUFFER
                MXLE = MXLE + STKPWIDTH
                LE = MIN( MXLE, STRLEN )
            END IF

C.............  Include elevated sources flag
            IF( RPT_%BYELEV ) THEN
                L = ELEVWIDTH
                L1 = L - LV  - 2      ! 1 for space, minus 1 for how used
                STRING = STRING( 1:LE ) // 
     &                   BLANK16( 1:L1 ) // BINELEV( I ) // DELIM
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
     &                   ' ' // CPDESC( S )( 1:L1 ) // DELIM
                MXLE = MXLE + L
                LE = MIN( MXLE, STRLEN )
            END IF

C.............  Include SCC description
C.............  This is knowingly including extra blanks before final quote
            IF( RPT_%SCCNAM ) THEN
                J = BINSNMIDX( I ) 
                L = SDSCWIDTH
                L1 = L - LV - 1                        ! 1 for space
                STRING = STRING( 1:LE ) // 
     &                   ' "'// SCCDESC( J )( 1:L1 )// '"' // DELIM
                MXLE = MXLE + L
                LE = MIN( MXLE, STRLEN )
            END IF

C.............  Remove leading spaces and get new length
            STRING = STRING( 2:LE )
            LE = LE - 1

C.............  Output error message of string is getting shortened
            IF( MXLE .GT. STRLEN ) THEN
                EFLAG = .TRUE.
                WRITE( MESG,94010 ) 'INTERNAL ERROR: Output string ' //
     &                 'getting truncated in report', RCNT
                CALL M3MSG2( MESG )
            END IF

C.............  If current bin has bad values, update those now.
            IF( BINBAD( I ) .GT. 0 ) BINDATA( I,1:NDATA ) = -9.

C.............  Write out this record
            WRITE( FDEV, OUTFMT ) STRING( 1:LE ), 
     &                          ( BINDATA( I,J ), J=1, NDATA )

        END DO  ! End loop through bins

C.........  Save report number for next time routine is called
        PRCNT = RCNT

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Formatted file I/O formats............ 93xxx

93000   FORMAT( A )

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I10, :, 1X ) )

94180   FORMAT( I2.2, '/', I2.2, '/', I4.4 )

        END SUBROUTINE WRREPOUT







