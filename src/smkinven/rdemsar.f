
        SUBROUTINE RDEMSAR( EDEV, INY, NRAWIN, WKSET, 
     &                      NRAWOUT, IOS, IREC, ERFILDSC,
     &                      EFLAG, NDROP, EDROP )

C***********************************************************************
C  subroutine body starts at line 120
C
C  DESCRIPTION:
C      This subroutine reads the EMS-95 format for the area source formatted
C      files. It can be called multiple times for multiple files.
C
C  PRECONDITIONS REQUIRED:
C      Files must be opened and their unit numbers stored in EDEV() in the
C      order listed in the description.
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C      Subroutines: I/O API subroutines, BLDCSRC, CHECKMEM 
C      Functions: I/O API functions, YR2DAY
C
C  REVISION  HISTORY:
C      Copied from rdemspt.f by M. Houyoux (11/99)
C
C****************************************************************************
C
C Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
C                System
C File: @(#)$Id$
C
C COPYRIGHT (C) 2002, MCNC Environmental Modeling Center
C All Rights Reserved
C
C See file COPYRIGHT for conditions of use.
C
C Environmental Modeling Center
C MCNC
C P.O. Box 12889
C Research Triangle Park, NC  27709-2889
C
C smoke@emc.mcnc.org
C
C Pathname: $Source$
C Last updated: $Date$ 
C
C***************************************************************************

C...........   MODULES for public variables
C...........   This module is the inventory arrays
        USE MODSOURC

C.........  This module contains the lists of unique inventory information
        USE MODLISTS

C.........  This module contains the information about the source category
        USE MODINFO

        IMPLICIT NONE

C...........   INCLUDES

        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
        INCLUDE 'CONST3.EXT'    !  physical constants
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
        INCLUDE 'FDESC3.EXT'    !  I/O API file description data structures.

C...........   EXTERNAL FUNCTIONS and their descriptions:
        CHARACTER*2            CRLF
        INTEGER                INDEX1
        INTEGER                STR2INT
        REAL                   STR2REAL
        REAL                   YR2DAY  

        EXTERNAL CRLF, INDEX1, STR2INT, STR2REAL, YR2DAY

C...........   SUBROUTINE ARGUMENTS
        INTEGER     , INTENT (IN) :: EDEV      !  unit no. for input file
        INTEGER     , INTENT (IN) :: INY       !  inv year for this set of files
        INTEGER     , INTENT (IN) :: NRAWIN    !  total raw record-count
        INTEGER     , INTENT (IN) :: WKSET     !  weekly profile interpretation
        INTEGER     , INTENT(OUT) :: NRAWOUT   ! valid raw record-count
        INTEGER     , INTENT(OUT) :: IOS       ! I/O status
        INTEGER     , INTENT(OUT) :: IREC      ! line number
        CHARACTER(*), INTENT(OUT) :: ERFILDSC  ! file desc of file in error
        LOGICAL     , INTENT(OUT) :: EFLAG     ! error flag 
        INTEGER    , INTENT(INOUT):: NDROP     ! number of records dropped
        REAL       , INTENT(INOUT):: EDROP( MXIDAT ) ! emis dropped per pol

C...........   Local parameters
        INTEGER, PARAMETER :: MXDATFIL = 60  ! arbitrary max no. data variables

C...........   Other local variables

        REAL             CEFF    !  Temporary control effectiveness
        REAL             DAY2YR  !  Local, leap-year-able, DAY to YEAR factor
        REAL             EMIS    !  Temporary emission value
        REAL             REFF    !  Temporary rule effectiveness
        REAL             RPEN    !  Temporary rule penetration

        INTEGER          I, J, L                ! counters and indices

        INTEGER          COD         !  temporary pollutant code number
        INTEGER          ES          !  counter for emission file
        INTEGER          FIP         !  temporary fip, scc, sic
        INTEGER       :: ICC     = 0 !  country code, def = 0
        INTEGER       :: NPOA    = 0 !  number of input data variables
        INTEGER, SAVE :: NSRCSAV = 0 ! cumulative source count
        INTEGER          TPF         !  temporary temporal ID
        INTEGER          YR4         !  tmp year from current file

        CHARACTER*2            TMPAA !  tmp time period code
        CHARACTER*10 , SAVE :: FIPFMT! formt to write co/st/cy to string
        CHARACTER*300          LINE  !  Input line from POINT file
        CHARACTER*300          MESG  !  Text for M3EXIT()
        CHARACTER(LEN=IOVLEN3) CPOL  !  Temporary pollutant code
        CHARACTER(LEN=FIPLEN3) CFIP  !  Character FIP code
        CHARACTER(LEN=POLLEN3) CCOD  !  Character pollutant index to INVDNAM
        CHARACTER(LEN=SCCLEN3) TSCC  !  Temporary character SCC

        CHARACTER*16 :: PROGNAME = 'RDEMSAR' ! Program name

C***********************************************************************
C   begin body of subroutine RDEMSAR

C.........  Set internal formats       
        WRITE( FIPFMT, '("(I",I2.2,".",I2.2,")")' ) FIPLEN3, FIPLEN3

C........  Set rule effectiveness and rule penetration to inventory defaults
        REFF = 100.0
        RPEN = 100.0

C........................................................................
C.............  Head of the input file read loop  .......................
C........................................................................

        ES   = NSRCSAV
        IREC = 0
        ERFILDSC = 'emission'
        DO 

C.............  Read a line of emission.pt file and check input status

            READ( EDEV, 93000, END=199, IOSTAT=IOS ) LINE
            IREC = IREC + 1

C.............  Check read i/o status
            IF ( IOS .NE. 0 ) THEN

                EFLAG = .TRUE.
                WRITE( MESG, 94010 )
     &              'I/O error', IOS, 
     &              'reading inventory file at line', IREC
                CALL M3MESG( MESG )
                CYCLE

            END IF

            IF ( LINE .EQ. ' ' ) CYCLE      ! skip if line is blank

C.............  Scan for header lines and check to ensure all are set 
C               properly (no header fields are required)
            CALL GETHDR( MXDATFIL, .FALSE., .FALSE., .FALSE., 
     &                   LINE, ICC, YR4, NPOA, IOS )

C.............  Set error flag for error from header reading routine
            IF( IOS .EQ. 4 ) THEN
                WRITE( MESG,94010 ) 
     &                 'Maximum allowed data variables ' //
     &                 '(MXDATFIL=', MXDATFIL, CRLF() // BLANK10 //
     &                 ') exceeded in input file'
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

            ELSE IF( IOS .GT. 0 ) THEN
                EFLAG = .TRUE.

            END IF

C.............  Resolve any differences between year from calling program
C               and year in file header
            IF ( IOS .EQ. 0   .AND.
     &           INY .GT. 0   .AND.
     &           YR4 .GT. 0   .AND.
     &           INY .NE. YR4       ) THEN
                WRITE( MESG,94010 ) 'NOTE: Using year', INY,
     &                 'from list file, and not year', YR4,
     &                 'from inventory file.'
                CALL M3MSG2( MESG )
                YR4 = INY
            ELSE IF ( INY .GT. 0 .AND.
     &                YR4 .LE. 1       ) THEN
                YR4 = INY
            END IF

C.............  If a header line was encountered, go to next line
            IF( IOS .GE. 0 ) CYCLE

C.............  Define day to year conversion factor and real type for integer 
C               missing value
            DAY2YR  = 1. / YR2DAY( YR4 )

C.............  Find pollutant name in master list to set index COD
C.............  NOTE- Pollutant names here and in INVDNAM are uppercase
            CPOL = ADJUSTL( LINE( 21:25 ) )
            CALL UPCASE( CPOL )
            COD  = INDEX1( CPOL, MXIDAT, INVDNAM )

            IF( COD .LE. 0 ) THEN
                L = LEN_TRIM( CPOL )
                WRITE( MESG,94010 )  'Source dropped: ' //
     &                 'pollutant name "' // CPOL( 1:L ) // 
     &                 '" in emission file at line', IREC,
     &                 CRLF() // BLANK5 // 
     &                 'is not in master pollutants list'
                CALL M3MESG( MESG )
                CYCLE      !  to head of loop
            END IF

C.............  Check and set emissions value

            EMIS = STR2REAL( LINE( 52:65 ) )
            IF ( EMIS .LT. 0.0 )  THEN
                WRITE( MESG,94010 )  'Source dropped: ' //
     &                 'bad emissions value "' // LINE( 52:65 ) //
     &                 '" in emission file at line', IREC
                CALL M3MESG( MESG )
                CYCLE
            END IF

C.............  Set country/state/county code
            FIP  = ICC  * 100000 +
     &             1000 * STR2INT( LINE( 1:2 ) ) +
     &                    STR2INT( LINE( 3:5 ) )
            WRITE( CFIP,FIPFMT ) FIP

C.............  Set source category code
            TSCC = ADJUSTL( LINE( 6:20 ) )

C.............  Set control efficiency
            CEFF = STR2REAL( LINE( 88 : 94 ) )

C.............  Check and set time period type (Year/day/hourly)
            TMPAA = LINE( 95:96 )
            CALL UPCASE( TMPAA )
            IF ( TMPAA .EQ. 'AA' ) THEN 

                TPF = MTPRFAC * WTPRFAC       !  use month, week profiles

            ELSE IF ( TMPAA .EQ. 'AD' ) THEN 

                TPF  = WKSET                !  use week profiles
                EMIS = DAY2YR * EMIS

            ELSE IF ( TMPAA .EQ. 'DS' ) THEN

                TPF = 1                     !  use only hourly profiles
                EMIS = DAY2YR * EMIS

            ELSE                            !  unrecognized type

                NDROP = NDROP + 1
                EDROP( COD ) = EDROP( COD ) + EMIS
                WRITE( MESG,94010 )  'Source dropped: unsupported ' //
     &                 'time period type "' // TMPAA //
     &                 '" in emission file at line', IREC
                CALL M3MESG( MESG )
                CYCLE          !  to head of MDEV-read loop

            END IF          !  tests on record type line( 57:58 )

C.............  Time to store data in unsorted lists if we've made it this far
            ES = ES + 1

            IF ( ES .LE. NRAWIN ) THEN

                IFIPA  ( ES ) = FIP
                TPFLGA ( ES ) = TPF
                INVYRA ( ES ) = YR4
                CSCCA  ( ES ) = TSCC
                POLVLA ( ES,NEM ) = EMIS
                POLVLA ( ES,NCE ) = CEFF
                POLVLA ( ES,NRE ) = REFF
                POLVLA ( ES,NRP ) = RPEN

                WRITE( CCOD,94125 ) COD
 
                CALL BLDCSRC( CFIP, TSCC, CHRBLNK3, CHRBLNK3, CHRBLNK3, 
     &                        CHRBLNK3, CHRBLNK3, CCOD, CSOURCA( ES ) )

            END IF      !  if ES in range

        ENDDO           !  to head of file read loop

199     CONTINUE        !  end of the file read loop

        CLOSE( EDEV )

        WRITE( MESG,94010 ) 
     &         'EMISSION FILE processed:'  // CRLF() // BLANK5 //
     &         '   This-file  EMS-95 SOURCE  record-count', ES-NSRCSAV,
     &         CRLF() // BLANK5 //
     &         '   Cumulative EMS-95 SOURCE  record-count', ES

        CALL M3MSG2( MESG )

        NSRCSAV = ES        !  cumulative emissions lines

C........................................................................

        IF( NSRCSAV .GT. NRAWIN ) THEN

            EFLAG = .TRUE.
            MESG = 'Memory allocation insufficient for EMS-95 inventory'
            CALL M3MSG2( MESG )

        ELSE
            NRAWOUT = NSRCSAV

        END IF		!  if overflow or if errors

999     CONTINUE

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Formatted file I/O formats............ 93xxx

93000   FORMAT( A )

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

94125   FORMAT( I5 )

94300   FORMAT( A, I2.2, A, I2.2, A )

        END SUBROUTINE RDEMSAR
