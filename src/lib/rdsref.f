
        SUBROUTINE RDSREF( FDEV, CATEGORY, NIPOL, NINVSCC, EINAM, 
     &                     INVSCC )

C***********************************************************************
C  subroutine body starts at line 
C
C  DESCRIPTION:
C     Reads the speciation cross-reference file for any source category.  It
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
C COPYRIGHT (C) 1998, MCNC--North Carolina Supercomputing Center
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

C...........   This module is for cross reference tables
        USE MODXREF

        IMPLICIT NONE

C...........   INCLUDES

        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters

C...........   EXTERNAL FUNCTIONS and their descriptions:
        CHARACTER*2     CRLF
        INTEGER         FIND1
        INTEGER         FINDC
        INTEGER         GETFLINE
        INTEGER         INDEX1
        INTEGER         STR2INT

        EXTERNAL  CRLF, FIND1, FINDC, GETFLINE, INDEX1, STR2INT

C...........   SUBROUTINE ARGUMENTS
        INTEGER     , INTENT (IN) :: FDEV              ! x-ref file unit no.
        CHARACTER(*), INTENT (IN) :: CATEGORY          ! source category name
        INTEGER     , INTENT (IN) :: NIPOL             ! no. of pols in inv
        INTEGER     , INTENT (IN) :: NINVSCC           ! no. of actl SCCs in inv
        CHARACTER(*), INTENT (IN) :: EINAM ( NIPOL )   ! pollutant names
        CHARACTER(*), INTENT (IN) :: INVSCC( NINVSCC ) ! list of SCCs in inven
 
C...........   Local parameters
        INTEGER    , PARAMETER :: AREATYP  = 1
        INTEGER    , PARAMETER :: MOBILTYP = 2
        INTEGER    , PARAMETER :: POINTTYP = 3
        INTEGER    , PARAMETER :: MXCHRS   = 7

        CHARACTER*6, PARAMETER :: LOCCATS( 3 ) = 
     &                         ( / 'AREA  ', 'MOBILE', 'POINT ' / )

C...........   Local arrays for reading speciation cross-referenc
        INTEGER, ALLOCATABLE :: LTYPE   ( : ) !  source category code per line

C...........   Left-portion of inventory SCCs
        INTEGER                :: NINVSCL
        CHARACTER(LEN=SCLLEN3) :: INVSCL( NINVSCC )

C...........   Sorted pollutant codes
        INTEGER                   INDXP  ( NIPOL ) !  sorting index for pols
        CHARACTER(LEN=IOVLEN3) :: SRTINAM( NIPOL ) !  sorted pollutant names

C...........   Array of source characeristics
        CHARACTER*100           CHARS( MXCHRS )
  
C...........   Arrays of column starts and ends (point sources)
        INTEGER :: PCS( MXCHRS ) = ( / 35, 42, 58, 74,  90, 106, 122 / )
        INTEGER :: PCE( MXCHRS ) = ( / 40, 56, 72, 88, 104, 120, 138 / )

C...........   Other local variables
        INTEGER         I, J, J1, J2, K, L, N    !  counters and indices

        INTEGER         COD     !  temporary pollutant code
        INTEGER         FIP     !  temporary FIPS code
        INTEGER         IDIU    !  temporary diurnal profile code     
        INTEGER         IMON    !  temporary monthly profile code     
        INTEGER         IOS     !  i/o status
        INTEGER         IWEK    !  temporary weekly profile code
        INTEGER         IREC    !  record counter
        INTEGER         JSCC    !  position of SCC in source characteristics
        INTEGER         JSPC    !  tmp index to master pollutant list
        INTEGER         LINTYPE !  temporary source category code
        INTEGER         LSA     !  ending position of left side of SCC
        INTEGER         LSB     !  starting position of right side of SCC
        INTEGER         LPCK    !  length of point definition packet
        INTEGER         NCHARS  !  number of source characteristics
        INTEGER      :: NCP = 0 !  input point source header parm
        INTEGER         NLINES  !  number of lines
        INTEGER         NREF    !  number of x-ref entries before filtering
        INTEGER         NXREF   !  number of valid x-ref entries
        INTEGER         RDT     !  temporary road class code
        INTEGER         THISTYP !  index in LOCCATS for CATEGORY
        INTEGER         VTYPE   !  temporary vehicle type number

        LOGICAL      :: EFLAG = .FALSE.   !  error flag
        LOGICAL      :: HFLAG = .FALSE.   !  header flag
        LOGICAL      :: PFLAG = .FALSE.   !  true when pol-spec entries skipped

        CHARACTER*1            SCC1       !  1st character of SCC
        CHARACTER*5            CCOD       !  temporary pol code or position

        CHARACTER*300          LINE     !  line buffer
        CHARACTER*300          MESG     !  message buffer
        CHARACTER(LEN=IOVLEN3) CBUF     !  temporary link code
        CHARACTER(LEN=ALLLEN3) CSRCALL  !  buffer for source char, incl pol
        CHARACTER(LEN=FIPLEN3) CFIP     !  buffer for CFIPS code
        CHARACTER(LEN=FIPLEN3) FIPZERO  !  buffer for zero FIPS code
        CHARACTER(LEN=SCCLEN3) TSCC     !  temporary SCC
        CHARACTER(LEN=SCCLEN3) SCCZERO  !  buffer for zero SCC
        CHARACTER(LEN=SCLLEN3) PSCCL    !  left digits of TSCC of prev iteration
        CHARACTER(LEN=SCLLEN3) SCCL5    !  left digits of TSCC
        CHARACTER(LEN=SCRLEN3) SCCR5    !  right 5 digits of TSCC
        CHARACTER(LEN=SCRLEN3) SCRZERO  !  buffer for zero SCCR5
        CHARACTER(LEN=SPNLEN3) SPCODE   !  tmp speciation profile code

        CHARACTER*16 :: PROGNAME = 'RDSREF' ! program name

C***********************************************************************
C   begin body of subroutine RDSREF

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
        SCRZERO = REPEAT( '0', SCRLEN3 )

C.........  Sort the actual list of pollutant names and store it
        DO I = 1, NIPOL
            INDXP( I ) = I
        ENDDO

        CALL SORTIC( NIPOL, INDXP, EINAM )

        DO I = 1, NIPOL
            J = INDXP( I )
            SRTINAM( I ) = EINAM( J )
        ENDDO

C.........  Starting position of right portion of SCC 
        LSB = SCCLEN3 - SCRLEN3 + 1
        LSA = LSB - 1

C.........  Create the list of the left portions of the SCCs and count
        J     = 0
        PSCCL = EMCMISS3
        DO I = 1, NINVSCC

            SCCL5 = INVSCC( I )( 1:LSA )
            IF( SCCL5 .NE. PSCCL ) THEN
                J = J + 1
                INVSCL( J ) = SCCL5
                PSCCL = SCCL5
            ENDIF

        ENDDO
        NINVSCL = J

C.........  Get the number of lines in the file
        NLINES = GETFLINE( FDEV, 'Speciation cross reference file' )

C.........  Allocate memory for storing the type of each line
        ALLOCATE( LTYPE( NLINES ), STAT=IOS )
        CALL CHECKMEM( IOS, 'LTYPE', PROGNAME )

C.........  Find index in list of categories for category of subroutine call.
C           Use it to  set LINTYPE when aren't sure if the record is for current
C           source category or not.
        THISTYP = INDEX1( CATEGORY, 3, LOCCATS )

C.........  First pass through file.  Determine format and count the number
C           of lines matching CATEGORY.  Do this so that we know how
C           much memory to allocate for the unsorted, unprocessed arrays.
        IREC = 0
        NREF = 0
        DO I = 1, NLINES

            READ( FDEV, 93000, END=999, IOSTAT=IOS ) LINE
            IREC = IREC + 1

            IF ( IOS .NE. 0 ) THEN
                EFLAG = .TRUE.
                WRITE( MESG,94010 ) 
     &              'I/O error', IOS, 
     &              'reading SPECIATION XREF file at line', IREC
                CALL M3MESG( MESG )
                CYCLE
            END IF

C.............  Skip blank lines
            IF( LINE .EQ. ' ' ) THEN
                CYCLE

C.............  Check for point source definition header
            ELSEIF( INDEX( LINE, PDEFPCKT ) .GT. 0 ) THEN
                HFLAG = .TRUE.
                CYCLE

            ELSE

                L    = LEN_TRIM( LINE )
                SCC1 = LINE( 1:1 )

C.................  Try to determine source category for this line
                IF( SCC1 .GT. '9' ) THEN  ! only mobile uses 'MV' in SCC
                    LINTYPE = MOBILTYP

                ELSEIF( L .GT. 40 ) THEN      ! only point has long lines
                    LINTYPE = POINTTYP

                ELSE                          ! either point or area
                    LINTYPE = THISTYP         ! ensures cycle IF below is false

                ENDIF

                LTYPE( I ) = LINTYPE

                IF( LOCCATS( LINTYPE ) .NE. CATEGORY ) CYCLE   ! to end of loop

                NREF = NREF + 1

C.................  Ensure that header is present for point sources and
                IF( LINTYPE .EQ. POINTTYP .AND. .NOT. HFLAG ) THEN
                    HFLAG = .TRUE.  ! To turn off error message
                    MESG = 'WARNING: ' // PDEFPCKT( 1:LPCK ) // 
     &                     ' header is not present before first ' //
     &                     'point source line.' // CRLF() // BLANK10 //
     &                     'Assuming inventory has IDA characteristics.'
                    CALL M3MSG2( MESG )
                ENDIF

            ENDIF

        ENDDO   ! End first pass through file

        REWIND( FDEV )

C.........  Allocate memory for unsorted data used in all source categories 
        ALLOCATE( CSPRNA( NREF ), STAT=IOS )
        CALL CHECKMEM( IOS, 'CSPRNA', PROGNAME )
        ALLOCATE( ISPTA( NREF ), STAT=IOS )
        CALL CHECKMEM( IOS, 'ISPTA', PROGNAME )
        ALLOCATE( INDXTA( NREF ), STAT=IOS )
        CALL CHECKMEM( IOS, 'INDXTA', PROGNAME )

C.........   Allocate memory for unsorted data used for a specific category
        SELECT CASE( CATEGORY )

        CASE( 'AREA' ) 

C NOTE: Insert these when the time comes

        CASE( 'MOBILE' )

C NOTE: Insert these when the time comes

        CASE( 'POINT' )
            ALLOCATE( CSCCTA( NREF ), STAT=IOS )
            CALL CHECKMEM( IOS, 'CSCCTA', PROGNAME )
            ALLOCATE( CSRCTA( NREF ), STAT=IOS )
            CALL CHECKMEM( IOS, 'CSRCTA', PROGNAME )

        END SELECT

C.........  Set up constants for loop.
C.........  Length of point definition packet, plus one
        LPCK = LEN_TRIM( PDEFPCKT ) + 1 

C.........  Format of point source fields - insert spaces
c        DO J = 1, 7
c            PSTARTIN( J ) = PTBEGL3( J ) +  J - 1	
c        ENDDO

C.........  Initialize character strings
        CHARS = ' ' ! array

C.........  Second pass through file: read lines and store unsorted data for
C           the source category of interest
        IREC   = 0
        N      = 0
        NCHARS = 3       ! IDA default
        JSCC   = 0       ! IDA default (although, not consistent w/ other files)
        DO I = 1, NLINES

            READ( FDEV, 93000, END=999, IOSTAT=IOS ) LINE
            IREC = IREC + 1

            IF ( IOS .NE. 0 ) THEN
                EFLAG = .TRUE.
                WRITE( MESG,94010 ) 
     &              'I/O error', IOS, 
     &              'reading SPECIATION XREF file at line', IREC
                CALL M3MESG( MESG )
                CYCLE
            END IF

C.............  Skip blank lines
            IF( LINE .EQ. ' ' ) CYCLE

            J = INDEX( LINE, PDEFPCKT ) ! can be in middle of file
            L = LEN_TRIM( LINE )

C.............  Read point source header information
            IF( J .GT. 0 ) THEN
                
                READ( LINE( LPCK:L ), * ) NCP, JSCC

C.................  Adjust for FIPS code and Plant ID
                NCHARS = NCP + 2
                IF( JSCC .GT. 0 ) JSCC = JSCC + 2   
                CYCLE

            ELSEIF( LOCCATS( LTYPE( I ) ) .EQ. CATEGORY ) THEN

C.................  Depending on source category, transfer line to temporary
C                   fields.  In cases where blanks are allowed, do not use
C                   STR2INT to prevent warning messages.
                SELECT CASE( CATEGORY )

                CASE( 'AREA' )    ! Same fields for area and mobile
                CASE( 'MOBILE' )

c                    FIP   = STR2INT( LINE(  1:6  ) )
c                    TSCC  =          LINE( 21:30 )
c                    CCOD  = ADJUSTL( LINE( 32:36 ) )

                CASE( 'POINT' )

C.....................  Right adjust SCC code and pad with zeros for blanks
                    TSCC = LINE( 1:10 )
                    CALL PADZERO( TSCC )

                    SPCODE = ADJUSTR( LINE( 12:16 ) ) !  profile number
                    CBUF   = ADJUSTL( LINE( 18:33 ) ) !  pollutant name

C.....................  Store string source characteristics 
                    DO J = 1, NCHARS
                        J1 = PCS( J )
                        J2 = PCE( J )
                        CHARS( J ) = LINE( J1:J2 )
                    ENDDO
                 
C.....................  Check for -9 in FIPS field and convert to zeros
                    CFIP = CHARS( 1 )
                    J = INDEX( CFIP, '-9' )
                    IF( J .GT. 0 .OR. CFIP .EQ. ' ' ) THEN
                        CFIP = FIPZERO
                    ENDIF
                    CALL PADZERO( CFIP )
                    CHARS( 1 ) = CFIP

                END SELECT

C................. Set left and right portions of SCC
                SCCL5 = TSCC(   1:LSA     )
                SCCR5 = TSCC( LSB:SCCLEN3 )

C................. Convert character SCC field to integer SCC number while
C                  allowing for case that SCC is blank.  If non-blank, compare
C                  with master SCC list for area and point sources.
                IF( TSCC .EQ. SCCZERO ) THEN
                    RDT   = 0
                    VTYPE = 0

                ELSEIF( CATEGORY .EQ. 'AREA' .OR. 
     &                  CATEGORY .EQ. 'POINT' ) THEN

                    IF( SCCR5 .EQ. SCRZERO ) THEN
                        J = FINDC( SCCL5, NINVSCL, INVSCL )
                    ELSE
                        J = FINDC( TSCC, NINVSCC, INVSCC )
                    ENDIF

                    IF( J .LE. 0 ) CYCLE  ! Skip entry if it doesn't apply

                ENDIF

C.................  Ensure that pollutant is in master list of pollutants or
C                   skip the pollutant-specific entry 
C.................  Filter the case where the pollutant code is not present
                IF( CBUF        .EQ. ' '  .OR. 
     &              CBUF( 1:2 ) .EQ. '-9' .OR.
     &              CBUF( 1:1 ) .EQ. '0'       ) THEN

                    EFLAG = .TRUE.
                    WRITE( MESG, 94010 ) 
     &                     'ERROR: Skipping cross-reference entry ' //
     &                     'at line', IREC, 
     &                     'because of missing pollutant.'
                    CALL M3MESG( MESG )
                    CYCLE

                ELSE
                    J = FINDC( CBUF, NIPOL, SRTINAM )

                    IF( J .LE. 0 ) THEN
                        PFLAG = .TRUE.  ! indicates skipped pol-spec entries
                        CYCLE
                    ELSE
                        JSPC = INDXP( J )
                    ENDIF

                ENDIF

C.................  Write pollutant position to character string
                WRITE( CCOD, '(I5.5)' ) JSPC  

C.................  Check for bad cross-reference code
                IF( SPCODE .EQ. ' ' ) THEN
                    WRITE( MESG, 94010 ) 
     &                'WARNING: Skipping blank profile code in cross-'//
     &                'reference file at line ', IREC
                    CALL M3MESG( MESG )
                    CYCLE

                ENDIF

                N = N + 1
                IF( N .GT. NREF ) CYCLE  ! Ensure no overflow

C.................  Store case-specific fields
                SELECT CASE( CATEGORY )

                CASE( 'AREA' )
                    IFIPTA( N ) = FIP
                    CSCCTA( N ) = TSCC

                CASE( 'MOBILE' )
C NOTE: Insert when the time comes

                CASE( 'POINT' )
 
                    CSCCTA( N ) = TSCC

C.....................  Store sorting criteria as right-justified in fields
                    CSRCALL = ' '
                    CALL BLDCSRC( CHARS(1), CHARS(2), CHARS(3),
     &                            CHARS(4), CHARS(5), CHARS(6),
     &                            CHARS(7), POLBLNK3, CSRCALL   )

                    CSRCTA( N ) = CSRCALL( 1:SRCLEN3 ) // TSCC // CCOD

                END SELECT

C.................  Store case-indpendent fields
                INDXTA( N ) = N
                ISPTA ( N ) = JSPC ! Save index to original PCODES or zero
                CSPRNA( N ) = SPCODE

            ENDIF  !  This line matches source category of interest

        ENDDO      ! End of loop on I for reading in speciation x-ref file

C.........  Reset number of cross-reference entries in case some were dropped
        NXREF = N

C.........  Write warning message for pollutants in cross-reference that are
C           not in master list
        IF( PFLAG ) THEN
            MESG = 'Pollutant-specific entries in the speciation ' //
     &             'cross-reference file have been skipped.'
            CALL M3WARN( PROGNAME, 0, 0, MESG )
        ENDIF

        IF( NXREF .EQ. 0 ) THEN
            EFLAG = .TRUE.
            MESG = 'No valid speciation cross-reference entries!'
            CALL M3MSG2( MESG )

        ELSEIF( NXREF .GT. NREF ) THEN
            EFLAG = .TRUE.
            WRITE( MESG,94010 ) 'INTERNAL ERROR: dimension for ' //
     &             'storing speciation cross-reference was', NREF,
     &             CRLF() // BLANK10 // 'but actually needed', NXREF
            CALL M3MSG2( MESG )

        ENDIF

C.......  Check for errors reading XREF file, and abort
        IF( EFLAG ) THEN
            MESG = 'Problem reading speciation cross-reference file.'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        ENDIF

        CALL M3MSG2( 'Processing speciation cross-reference file...' )

C.........  Sort speciation cross-reference entries. Since CCOD was used in 
C           building CSRCTA, and CCOD will equal "0" when the x-ref entry is
C           not pollutant-specific, the non-pollutant-specific entries will
C           always appear first.  This is necessary for the table-generating
C           subroutines.
        SELECT CASE( CATEGORY )

        CASE( 'AREA' ) 

        CASE( 'MOBILE' ) 

        CASE( 'POINT' ) 

            CALL SORTIC( NXREF, INDXTA, CSRCTA )

            CALL PXREFTBL( 'SPECIATION', NXREF, MXCHRS, NCHARS, JSCC,
     &                     NIPOL, EINAM )

C.............  Deallocate source-specific temporary unsorted arrays
            DEALLOCATE( CSCCTA, CSRCTA )

        END SELECT

C.........  Deallocate other temporary unsorted arrays
        DEALLOCATE( LTYPE, CSPRNA, INDXTA )

C.........  Rewind file
        REWIND( FDEV )

        RETURN

C.........  Error message for reaching the end of file too soon
999     MESG = 'End of file reached unexpectedly. ' //
     &         'Check format of speciation' // CRLF() // BLANK5 //
     &         'cross-reference file.'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

C******************  FORMAT  STATEMENTS   ******************************

C...........   Formatted file I/O formats............ 93xxx

93000   FORMAT( A )

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

        END SUBROUTINE RDSREF
