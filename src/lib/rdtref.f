
        SUBROUTINE RDTREF( FDEV, FFORMAT )

C***********************************************************************
C  subroutine body starts at line 
C
C  DESCRIPTION:
C     Reads the temporal cross-reference file for any source category.  It
C     allocates memory (locally) for reading the unsorted x-refs. It sorts the
C     x-refs for processing. It allocates memory for the appropriate x-ref 
C     tables and populates the tables (passed via modules).
C
C  PRECONDITIONS REQUIRED:
C     File unit FDEV already is opened... MORE
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C     Created 1/99 by M. Houyoux
C
C****************************************************************************/
C
C Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
C                System
C File: @(#)$Id$
C
C COPYRIGHT (C) 1999, MCNC--North Carolina Supercomputing Center
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
C***************************************************************************

C.........  MODULES for public variables
C...........   This module is for cross reference tables
        USE MODXREF

C.........  This module contains the lists of unique source characteristics
        USE MODLISTS

C.........  This module contains the information about the source category
        USE MODINFO

        IMPLICIT NONE

C...........   INCLUDES

        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters

C...........   EXTERNAL FUNCTIONS and their descriptions:
        LOGICAL         CHKINT
        CHARACTER*2     CRLF
        INTEGER         GETNLIST
        INTEGER         GETFLINE
        INTEGER         INDEX1
        INTEGER         STR2INT

        EXTERNAL  CHKINT, CRLF, GETNLIST, GETFLINE, INDEX1, STR2INT

C...........   SUBROUTINE ARGUMENTS
        INTEGER     , INTENT (IN) :: FDEV              ! x-ref file unit no.
        CHARACTER(*), INTENT(OUT) :: FFORMAT           ! format of input file
 
C...........   Local parameters
        INTEGER    , PARAMETER :: AREATYP  = 1
        INTEGER    , PARAMETER :: MOBILTYP = 2
        INTEGER    , PARAMETER :: POINTTYP = 3
        INTEGER    , PARAMETER :: MXTCOL   = 12

        CHARACTER*6, PARAMETER :: LOCCATS( 3 ) = 
     &                         ( / 'AREA  ', 'MOBILE', 'POINT ' / )

C...........   Sorted pollutant/emission type names
        INTEGER                   INDXP  ( NIPPA ) !  sorting index for pol/etyp
        CHARACTER(LEN=IOVLEN3) :: SRTINAM( NIPPA ) !  sorted pol/act names

C...........   Array of point source plant characeristics
        CHARACTER(LEN=CHRLEN3) CHARS( 5 )

C...........   Array for reading temporal x-ref fields
        CHARACTER*20            SEGMENT( MXTCOL )

C...........   Other local variables
        INTEGER         I, J, J1, J2, J3, K, L, N    !  counters and indices

        INTEGER         COD     !  temporary pollutant/emission type code
        INTEGER         FIP     !  temporary FIPS code
        INTEGER         IDIU    !  temporary diurnal profile code
        INTEGER      :: IDUM = 0!  tmp dummy integer
        INTEGER         IMON    !  temporary monthly profile code
        INTEGER         IOS     !  i/o status
        INTEGER         IWEK    !  temporary weekly profile code
        INTEGER         IREC    !  record counter
        INTEGER      :: JS = 0  !  position of SCC in source chars in x-ref file
        INTEGER         JSPC    !  tmp index to master pollutant/etype list
        INTEGER         LINTYPE !  temporary source category code
        INTEGER         LPCK    !  length of point definition packet
        INTEGER      :: NCP = 0 !  input point source header parm
        INTEGER         NFIELD  !  tmp number of fields in LINE
        INTEGER         NLINES  !  number of lines
        INTEGER         NREF    !  number of x-ref entries before filtering
        INTEGER         NXREF   !  number of valid x-ref entries
        INTEGER         RDT     !  temporary road class code
        INTEGER         VTYPE   !  temporary vehicle type number

        LOGICAL      :: EFLAG = .FALSE.   !  true: error occurred
        LOGICAL      :: HFLAG = .FALSE.   !  true: pt defn header encountered
        LOGICAL      :: PFLAG = .FALSE.   !  true: pol/act-spec entries skipped
        LOGICAL      :: SKIPREC = .FALSE. !  true: skip this x-ref entry

        CHARACTER*1            SCC1     !  1st character of SCC
        CHARACTER*5            CPOS     !  tmp sorted position of pol/act
        CHARACTER*300          LINE     !  line buffer
        CHARACTER*300          MESG     !  message buffer

        CHARACTER(LEN=SICLEN3) :: CDUM = '0' !  dummy character field for SIC
        CHARACTER(LEN=LNKLEN3) CLNK     !  temporary link code
        CHARACTER(LEN=ALLLEN3) CSRCALL  !  buffer for source char, incl pol/act
        CHARACTER(LEN=FIPLEN3) CFIP     !  buffer for CFIPS code
        CHARACTER(LEN=FIPLEN3) FIPZERO  !  buffer for zero FIPS code
        CHARACTER(LEN=SCCLEN3) TSCC     !  temporary SCC
        CHARACTER(LEN=SCCLEN3) SCCZERO  !  buffer for zero SCC
        CHARACTER(LEN=SCCLEN3) SCCTEST  !  buffer for tester SCC for FLTRXREF
        CHARACTER(LEN=PLTLEN3) PLT      !  tmp plant ID
        CHARACTER(LEN=SCCLEN3) PSCCL    !  left digits of TSCC of prev iteration
        CHARACTER(LEN=SCCLEN3) SCCL     !  left digits of TSCC
        CHARACTER(LEN=SCCLEN3) SCCR     !  right 5 digits of TSCC
        CHARACTER(LEN=SCCLEN3) SCRZERO  !  buffer for zero SCCR
        CHARACTER(LEN=IOVLEN3) CPOA     !  temporary pollutant/emission type
        CHARACTER(LEN=RWTLEN3) CRWT     !  roadway type no.
        CHARACTER(LEN=VIDLEN3) CVID     !  vehicle type ID no.

        CHARACTER*16 :: PROGNAME = 'RDTREF' ! program name

C***********************************************************************
C   begin body of subroutine RDTREF

C.........  Ensure that the CATEGORY is valid
        I = INDEX1( CATEGORY, 3, LOCCATS )

        IF( I .LE. 0 ) THEN
            L = LEN_TRIM( CATEGORY )
            MESG = 'INTERNAL ERROR: category "' // CATEGORY( 1:L ) // 
     &             '" is not valid in routine ' // PROGNAME
            CALL M3MSG2( MESG ) 
            CALL M3EXIT( PROGNAME, 0, 0, ' ', 2 ) 

        ENDIF

C.........  Set up zero strings for FIPS code of zero and SCC code of zero
        FIPZERO = REPEAT( '0', FIPLEN3 )
        SCCZERO = REPEAT( '0', SCCLEN3 )
        SCRZERO = REPEAT( '0', SCCLEN3 - LSCCEND )

C.........  Sort the actual list of pollutant/emission type names and store it
        DO I = 1, NIPPA
            INDXP( I ) = I
        END DO

        CALL SORTIC( NIPPA, INDXP, EANAM )

        DO I = 1, NIPPA
            J = INDXP( I )
            SRTINAM( I ) = EANAM( J )
        END DO

C.........  Write status message
        MESG = 'Reading temporal cross-reference file...'
        CALL M3MSG2( MESG )

C.........  Set up constants for loop...

C.........  Length of point definition packet, plus one
        LPCK = LEN_TRIM( PDEFPCKT ) + 1 

C.........  Get the number of lines in the file
        NLINES = GETFLINE( FDEV, 'Temporal cross reference file' )

C.........  Initialize character strings
        SEGMENT = ' '  ! array

C.........   First pass through file.  Determine format and count the number
C            of lines matching SCC list and pol/act list.  Do this so that 
C            we know how much memory to allocate for the unsorted, unprocessed 
C            arrays.
        IREC = 0
        NREF = 0
        DO I = 1, NLINES

            READ( FDEV, 93000, END=999, IOSTAT=IOS ) LINE
            IREC = IREC + 1

            IF ( IOS .NE. 0 ) THEN
                EFLAG = .TRUE.
                WRITE( MESG,94010 ) 
     &              'I/O error', IOS, 
     &              'reading TEMPORAL XREF file at line', IREC
                CALL M3MESG( MESG )
                CYCLE
            END IF

            L = LEN_TRIM( LINE )
            NFIELD = GETNLIST( L,LINE ) 

C.............  Skip blank lines
            IF( LINE .EQ. ' ' ) THEN
                CYCLE

C.............  Auto check for LIST or EPS2.0 formated temporal x-ref...
C.............  If header is found, read point source header information
            ELSE IF( INDEX( LINE, PDEFPCKT ) .GT. 0 ) THEN
                HFLAG = .TRUE.


                READ( LINE( LPCK:L ), * ) NCP, JS

C.................  Adjust for FIPS code and Plant ID, which are always there
                NCP = NCP + 2
                IF( JS .GT. 0 ) JS = JS + 2

C.................  Compare point source definition from header to inventory
                IF ( CATEGORY .EQ. 'POINT' ) CALL CHKPTDEF( NCP, JS )

                CYCLE

C.............  The source-formatted will have only 3 columns                
            ELSE IF( NFIELD .LE. 3 ) THEN

                FFORMAT = 'SOURCE' 
                NREF = NLINES     ! equals the number of sources
                EXIT              ! Exit from read loop

            ELSE
                FFORMAT = 'STANDARD'

                CALL PARSLINE( LINE, MXTCOL, SEGMENT )

C.................  Make sure SCC is set to SCCZERO if it is missing
                TSCC = ADJUSTL( SEGMENT( 1 ) )
                CALL FLTRNEG( TSCC )
                CALL PADZERO( TSCC )

C.................  Get SCC from source definition if it is defined already
                IF( JS .GT. 0 .AND. TSCC .EQ. SCCZERO ) THEN
                    TSCC = SEGMENT( 5 + JS )    ! from source definition
                END IF

                CPOA = SEGMENT( 5 )   ! pollutant/emission type name
                CFIP = SEGMENT( 6 )   ! country/state/county code

                IF( CATEGORY .NE. 'MOBILE' ) THEN
                    SCCTEST = TSCC
                ELSE
                    SCCTEST = SCCZERO
                END IF

C.................  Post-process x-ref information to scan for '-9', pad
C                   with zeros, compare SCC version master list, compare
C                   SIC version to master list, and compare pol/act name 
C                   with master list.
                CALL FLTRXREF( CFIP, CDUM, SCCTEST, CPOA, IDUM, 
     &                         IDUM, JSPC, PFLAG, SKIPREC  )

C.................  Skip lines that are not valid for this inven and src cat
                IF( SKIPREC ) CYCLE

C.................  Ensure that header is present for point sources and
                IF( CATEGORY .EQ. 'POINT' .AND. .NOT. HFLAG ) THEN
                    EFLAG = .TRUE.
                    HFLAG = .TRUE.  ! To turn off error message
                    MESG = 'ERROR: ' // PDEFPCKT( 1:LPCK ) // 
     &                     ' header is not present before first ' //
     &                     'point source line.'
                    CALL M3MSG2( MESG )

                ELSE
                    NREF = NREF + 1

                END IF

            END IF     ! End format of temporal x-ref file

        END DO         ! End first pass through file

        REWIND( FDEV )

C.........  Allocate memory for unsorted data used in all source categories and
C           input formats
        ALLOCATE( MPRNA( NREF ), STAT=IOS )
        CALL CHECKMEM( IOS, 'MPRNA', PROGNAME )
        ALLOCATE( WPRNA( NREF ), STAT=IOS )
        CALL CHECKMEM( IOS, 'WPRNA', PROGNAME )
        ALLOCATE( DPRNA( NREF ), STAT=IOS )
        CALL CHECKMEM( IOS, 'DPRNA', PROGNAME )

C.........  Special section for list-formatted cross-reference file (used only
C           with EMS-95 inputs).
        IF( FFORMAT .EQ. 'SOURCE' ) THEN

            DO I = 1, NLINES

                READ( FDEV, *, IOSTAT=IOS ) IMON, IWEK, IDIU

                IF( IOS .GT. 0 ) THEN
                    EFLAG = .TRUE.
                    WRITE( MESG,94010 ) 
     &                  'I/O error', IOS, 
     &                  'reading TEMPORAL XREF file at line', IREC
                    CALL M3MESG( MESG )
                    CYCLE
                END IF

                MPRNA( I ) = IMON
                WPRNA( I ) = IWEK
                DPRNA( I ) = IDIU

            END DO  !  End of loop for reading list-formatted xref file

        END IF

C.........  Check for errors from reading either format
        IF( EFLAG ) THEN
            MESG = 'Problem reading temporal cross-reference file.'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

C.........  Report format for XREF file
        ELSE
            MESG = 'NOTE: File read in as ' // FFORMAT // ' format.'
            CALL M3MSG2( MESG )

        END IF

C.........  Leave routine for list-formatted file because no grouping needed
        IF( FFORMAT .EQ. 'SOURCE' ) RETURN

C.........  FFORMAT = 'SOURCE' DOES NOT APPLY AFTER THIS POINT

C.........  Allocate memory for unsorted data used in all source categories
        ALLOCATE( INDXTA( NREF ), STAT=IOS )
        CALL CHECKMEM( IOS, 'INDXTA', PROGNAME )
        ALLOCATE( ISPTA( NREF ), STAT=IOS )
        CALL CHECKMEM( IOS, 'ISPTA', PROGNAME )
        ALLOCATE( CSCCTA( NREF ), STAT=IOS )
        CALL CHECKMEM( IOS, 'CSCCTA', PROGNAME )
        ALLOCATE( CSRCTA( NREF ), STAT=IOS )
        CALL CHECKMEM( IOS, 'CSRCTA', PROGNAME )

C.........  Initialize character strings
        CHARS   = ' '  ! array
        SEGMENT = ' '  ! array

C.........  Second pass through file: read lines and store unsorted data for
C           the source category of interest
        IREC = 0
        N    = 0
        DO I = 1, NLINES

            READ( FDEV, 93000, END=999, IOSTAT=IOS ) LINE
            IREC = IREC + 1

            IF ( IOS .NE. 0 ) THEN
                EFLAG = .TRUE.
                WRITE( MESG,94010 ) 
     &              'I/O error', IOS, 
     &              'reading temporal cross-reference file at line',IREC
                CALL M3MESG( MESG )
                CYCLE
            END IF

C.............  Skip blank lines
            IF( LINE .EQ. ' ' ) CYCLE

            J = INDEX( LINE, PDEFPCKT ) ! can be in middle of file
            L = LEN_TRIM( LINE )

C.............  Skip point source header information and filter the line again 
C               for storing info
            IF( J .LE. 0 ) THEN

C.................  Separate line from file into fields
                CALL PARSLINE( LINE, MXTCOL, SEGMENT )

C.................  Make sure SCC is set to SCCZERO if it is missing
                TSCC = ADJUSTL( SEGMENT( 1 ) )
                CALL FLTRNEG( TSCC )
                CALL PADZERO( TSCC )

C.................  Get SCC from source definition if it is defined already
                IF( JS .GT. 0 .AND. TSCC .EQ. SCCZERO ) THEN
                    TSCC = SEGMENT( 5 + JS )    ! from source definition
                END IF

                CPOA = SEGMENT( 5 )   ! pollutant/emission type name
                CFIP = SEGMENT( 6 )   ! country/state/county code

                IF( CATEGORY .NE. 'MOBILE' ) THEN
                    SCCTEST = TSCC
                ELSE
                    SCCTEST = SCCZERO
                END IF

C.................  Post-process x-ref information to scan for '-9', pad
C                   with zeros, compare SCC version master list, compare
C                   SIC version to master list, and compare pol/act name 
C                   with master list.
                CALL FLTRXREF( CFIP, CDUM, SCCTEST, CPOA, IDUM, 
     &                         IDUM, JSPC, PFLAG, SKIPREC    )

C.................  Skip lines that are not valid for this inven and src cat
                IF( SKIPREC ) CYCLE

C.................  Check for integers for temporal profile numbers
                IF( .NOT. CHKINT( SEGMENT( 2 ) ) .OR. 
     &              .NOT. CHKINT( SEGMENT( 3 ) ) .OR.
     &              .NOT. CHKINT( SEGMENT( 4 ) )      ) THEN
                    EFLAG = .TRUE.
                    WRITE( MESG,94010 )'ERROR: temporal profile '//
     &                     'code(s) at line', IREC, 'of temporal'//
     &                     CRLF() // BLANK10 // 
     &                     'cross-reference file are non-integer.'
                    CYCLE

C.................  Convert temporal profile numbers
                ELSE
                    IMON = STR2INT( SEGMENT( 2 ) )
                    IWEK = STR2INT( SEGMENT( 3 ) )
                    IDIU = STR2INT( SEGMENT( 4 ) )

C.....................  Check for bad cross-reference code
                    IF( IMON .LE. 0 .OR. IWEK .LE. 0 .OR.
     &                  IDIU .LE. 0 ) THEN
                        WRITE( MESG, 94010 ) 
     &                   'WARNING: Skipping bad cross-reference code '//
     &                   'at line. Values are:' // CRLF() // BLANK16 //
     &                   'Monthly=', IMON, 'Weekly=', IWEK, 'Diurnal=',
     &                   IDIU
                        CALL M3MESG( MESG )
                        CYCLE

                    END IF

                END IF
 
C.................  Write pol/act position to a character string
                WRITE( CPOS, '(I5)' ) JSPC  

                N = N + 1
                IF( N .GT. NREF ) CYCLE  ! Ensure no overflow

                CSRCALL = ' '
C.................  Store sorting criteria as right-justified in fields
C.................  For mobile sources, retrieve link from character field
C                   and extract road class and vehicle type from SCC
C.................  For point sources, retrieve plant + characteristics
                SELECT CASE( CATEGORY )

                CASE( 'AREA' )

                    CALL BLDCSRC( CFIP, CHRBLNK3, CHRBLNK3, CHRBLNK3,
     &                            CHRBLNK3, CHRBLNK3, CHRBLNK3,
     &                            POLBLNK3, CSRCALL )

                CASE( 'MOBILE' )

C.....................  Convert TSCC to internal value
                    CALL MBSCCADJ( IREC, TSCC, CRWT, CVID, TSCC, EFLAG )
    
                    CLNK = SEGMENT( 7 )
                    CALL FLTRNEG( CLNK )

C.....................  Reset road class to blank if no link.  Road class is
C                       now stored as part of SCC
                    IF( CLNK .EQ. ' ' ) CRWT = ' '

                    CALL BLDCSRC( CFIP, CRWT, CLNK, CHRBLNK3,
     &                            CHRBLNK3, CHRBLNK3, CHRBLNK3,
     &                            POLBLNK3, CSRCALL )

                CASE( 'POINT' )

C.....................  Store string source characteristics 
                    PLT = SEGMENT ( 7 )
                    CHARS( 1:5 ) = SEGMENT( 8:MXTCOL )

                    CALL BLDCSRC( CFIP, PLT, CHARS(1),
     &                            CHARS(2), CHARS(3), CHARS(4),
     &                            CHARS(5), POLBLNK3, CSRCALL )

                END SELECT

C.................  Store fields
                CSCCTA( N ) = TSCC
                CSRCTA( N ) = CSRCALL( 1:SRCLEN3 ) // TSCC // CPOS

C.................  Store case-indpendent fields
                INDXTA( N ) = N
                ISPTA ( N ) = JSPC ! Save index to original EANAM or zero
                MPRNA ( N ) = IMON
                WPRNA ( N ) = IWEK
                DPRNA ( N ) = IDIU

            END IF  !  This line matches source category of interest

        END DO      ! End of loop on I for reading in temporal x-ref file

C.........  Set actual number of cross-reference entries
        NXREF = N

C.........  Check if no cross-reference entries were found
        IF( NXREF .EQ. 0 ) THEN
            MESG = 'ERROR: No valid temporal cross-reference entries '//
     &             'were found'
            CALL M3MSG2( MESG )
            EFLAG = .TRUE.
        END IF

C.........  Check for errors reading cross-reference file, and abort
        IF( EFLAG ) THEN
            MESG = 'Problem reading temporal cross-reference file.'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

        CALL M3MSG2( 'Processing temporal cross-reference file...' )

C.........  Sort temporal cross-reference entries. Since CPOS was used in 
C           building CSRCTA, and CPOS will equal "0" when the x-ref entry is
C           not pol/act-specific, the non-pol/act-specific entries will
C           always appear first.  This is necessary for the table-generating
C           subroutines.

        CALL SORTIC( NXREF, INDXTA, CSRCTA )

        CALL XREFTBL( 'TEMPORAL', NXREF )

C.........  Deallocate temporary unsorted arrays
        DEALLOCATE( INDXTA, ISPTA, CSCCTA, CSRCTA, MPRNA, WPRNA, DPRNA )

C.........  Rewind file
        REWIND( FDEV )

        RETURN

C.........  Error message for reaching the end of file too soon
999     MESG = 'End of file reached unexpectedly. ' //
     &         'Check format of temporal' // CRLF() // BLANK5 //
     &         'cross reference file.'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

C******************  FORMAT  STATEMENTS   ******************************

C...........   Formatted file I/O formats............ 93xxx

93000   FORMAT( A )

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

        END SUBROUTINE RDTREF
