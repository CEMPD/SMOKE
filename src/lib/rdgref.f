
        SUBROUTINE RDGREF( FDEV )

C***********************************************************************
C  subroutine body starts at line 
C
C  DESCRIPTION:
C     Reads the gridding cross-reference file for area or mobile sources.  It
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
C     Created 4/99 by M. Houyoux
C     Modified 12/01 by Gabe Cano - deterministic mode
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
C.........  This module is for cross reference tables
        USE MODXREF

C.........  This module contains the lists of unique source characteristics
        USE MODLISTS

C.........  This module contains the information about the source category
        USE MODINFO

        IMPLICIT NONE

C...........   INCLUDES

        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters

C...........   EXTERNAL FUNCTIONS and their descriptions:
c        LOGICAL         BLKORCMT
        LOGICAL         CHKINT
        LOGICAL         CHKREAL
        CHARACTER*2     CRLF
        INTEGER         FIND1
        INTEGER         FINDC
        INTEGER         GETNLIST
        INTEGER         INDEX1
        INTEGER         STR2INT
        REAL            STR2REAL

c        EXTERNAL  BLKORCMT, CHKINT, CRLF, FIND1, FINDC, GETNLIST,
        EXTERNAL  CHKINT, CHKREAL, CRLF, FIND1, FINDC, GETNLIST,
     &            INDEX1, STR2INT, STR2REAL

C...........   SUBROUTINE ARGUMENTS
        INTEGER, INTENT (IN) :: FDEV   ! cross-reference file unit no.
 
C...........   Local parameters
        INTEGER    , PARAMETER :: AREATYP  = 1
        INTEGER    , PARAMETER :: BIOGTYP  = 2
        INTEGER    , PARAMETER :: MOBILTYP = 3

        CHARACTER*6, PARAMETER :: LOCCATS( 3 ) = 
     &                         ( / 'AREA  ', 'BIOG  ', 'MOBILE' / )

C...........   Array of input fields
        CHARACTER(LEN=SCCLEN3)  FIELDARR( 2 )
  
C...........   Other local variables
        INTEGER         C, I, J, J1, J2, L, N  !  counters and indices

        INTEGER         FIP     !  temporary FIPS code
        INTEGER         IDUM    !  dummy variable
        INTEGER         IOS     !  i/o status
        INTEGER         IREC    !  record counter
        INTEGER         ISRG    !  tmp surrogate ID
        INTEGER         JS      !  position of SCC in source chars in x-ref file
        INTEGER         LINTYPE !  temporary source category code
        INTEGER         LPCK    !  length of point definition packet
        INTEGER      :: NCP = 0 !  input point source header parm
        INTEGER         NLINES  !  maximum number of lines from input file
        INTEGER         NSPMAX  !  maximum number of SRC/PROB pairs from file
        INTEGER         NPAIR   !  number of SRC/PROB pairs from input file
        INTEGER         NCOLTOT !  total number of string fields in LINE
        INTEGER         NXREF   !  number of valid x-ref entries
        INTEGER         THISTYP !  index in LOCCATS for CATEGORY
        INTEGER         VTYPE   !  temporary vehicle type number

        LOGICAL      :: EFLAG = .FALSE.   !  true: error found
        LOGICAL      :: FIRSTIME = .TRUE. !  true: for loading NLINES and NSPMAX
        LOGICAL      :: LDUM  = .FALSE.   !  dummy
        LOGICAL         LOADED            !  true: after call to LOAD_FIELDARR
        LOGICAL      :: SKIPREC = .FALSE. !  true: record skipped in x-ref file

        REAL            RPROB   !  tmp surrogate probability

        CHARACTER*2            SCC2     !  1st & 2nd character of SCC
        CHARACTER*300          LINE     !  line buffer
        CHARACTER*300          MESG     !  message buffer
        CHARACTER(LEN=ALLLEN3) CSRCALL  !  buffer for source char, incl pol
        CHARACTER(LEN=FIPLEN3) CFIP     !  buffer for CFIPS code
        CHARACTER(LEN=RWTLEN3) CRWT     !  buffer for roadway type
        CHARACTER(LEN=SICLEN3) CDUM     !  dummy buffer for SIC code
        CHARACTER(LEN=SCCLEN3) CHKZERO  !  buffer to check for zero SCC
        CHARACTER(LEN=SCCLEN3) SCCZERO  !  zero SCC
        CHARACTER(LEN=SCCLEN3) TSCC     !  temporary SCC or roadway type
        CHARACTER(LEN=VIDLEN3) CVID     !  buffer for vehicle type ID

        CHARACTER*16 :: PROGNAME = 'RDGREF' ! program name

C***********************************************************************
C   begin body of subroutine RDGREF

C.........  Ensure that the CATEGORY is valid
C.........  Use THISTYP to  et LINTYPE when aren't sure if the record is for
C           current source category or not.
        THISTYP = INDEX1( CATEGORY, 3, LOCCATS )

        IF( THISTYP .LE. 0 ) THEN
            L = LEN_TRIM( CATEGORY )
            MESG = 'INTERNAL ERROR: category "' // CATEGORY( 1:L ) // 
     &             '" is not valid in routine ' // PROGNAME
            CALL M3MSG2( MESG ) 
            CALL M3EXIT( PROGNAME, 0, 0, ' ', 2 ) 

        ENDIF

C.........  Create the zero SCC
        SCCZERO = REPEAT( '0', SCCLEN3 )

C.........  Set up allocation constants
        NLINES = 0
        NSPMAX = 0

C.............  First pass: subfunction gets file line count, checks for 
C               the number of string field in LINE buffer, and gets the 
C               maximum number of
C.............  Second and suceeding passes through file: read lines and 
C               store unsorted data for the source category of interest
        DO 

C.............  Set up work constants for loop.
            IREC   = 0
            N      = 0
            CDUM   = ' ' 

            DO

                READ( FDEV, 93000, END=100, IOSTAT=IOS ) LINE
                IREC = IREC + 1

                IF ( IOS .NE. 0 ) THEN
                    EFLAG = .TRUE.
                    WRITE( MESG,94010 ) 
     &                  'I/O error', IOS, 
     &                  'reading GRIDDING X-REF file at line', IREC
                    CALL M3MESG( MESG )
                    CYCLE
                ENDIF

C.................  Skip blank lines
c                IF( BLKORCMT( LINE ) ) CYCLE
                IF( LINE .EQ. '!' .OR. LINE .EQ. ' ') CYCLE

C................. get number of column fields
                L = LEN_TRIM ( LINE )
                C = GETNLIST ( L , LINE )

                IF ( C .LT. 4 ) THEN

                    WRITE( MESG,94010 ) 
     &                  'WARNING: at line', IREC, 
     &                  'of GRIDDING X-REF file.  Expected ' //
     &                  'at least 4 data fields but found ', C
                    CALL M3MESG( MESG )
                    CYCLE

                ELSE IF ( MOD( C, 2 ) .GT. 0.0 ) THEN

                    WRITE( MESG,94010 ) 
     &                  'WARNING: at line', IREC, 
     &                  'of GRIDDING X-REF file.  Found ', C,
     &                  'data fields but expected an even number'
                    CALL M3MESG( MESG )
                    CYCLE

                ENDIF  ! end checks for column field count, C

C................. save actual number of column fields
                NCOLTOT = C

C................. count only surrogate and associated probability pairs
                NPAIR = INT( NCOLTOT / 2 ) - 1
                N = N + 1

                IF ( FIRSTIME ) THEN

                    IF ( NPAIR .GT. NSPMAX ) NSPMAX = NPAIR 
                    CYCLE

                ELSE

C.....................  Ensure no overflow
                    IF( N .GT. NLINES ) CYCLE  

C.....................  Prepare for first read from LINE buffer
                    LOADED = .FALSE.

C.................  Depending on source category, transfer line to temporary
C                   fields.  In cases where blanks are allowed, do not use
C                   STR2INT to prevent warning messages.
                    SELECT CASE( CATEGORY )
          
                    CASE( 'AREA' ) 
c                        CALL PARSLINE( LINE, 3 , FIELDARR )
                        CALL LOAD_FIELDARR

                        CFIP = FIELDARR( 1 )
                        TSCC = FIELDARR( 2 )

C.........................  Post-process x-ref information to scan for '-9', 
C                           pad with zeros, compare SCC version master list.
                        CALL FLTRXREF( CFIP, CDUM, TSCC, ' ', IDUM,
     &                                 IDUM, IDUM, LDUM, SKIPREC )

                    CASE( 'MOBILE' )
c                        CALL PARSLINE( LINE, 3, FIELDARR )
                        CALL LOAD_FIELDARR
          
                        L = LEN_TRIM( FIELDARR( 2 ) )

C.........................  Make sure SCC is full length
                        IF( L .NE. SCCLEN3 ) THEN
                            EFLAG = .TRUE.
                            WRITE( MESG,94010 ) 'ERROR: SCC value ' //
     &                             'is not', SCCLEN3, 'digits at line',
     &                             IREC
                            CALL M3MESG( MESG )

                        ENDIF
            
                        CFIP = FIELDARR( 1 )
                        TSCC = FIELDARR( 2 )
            
C.........................  Post-process x-ref information to scan for 
C                           '-9', pad with zeros.  Do not include SCC 
C                           in call below because right SCC will not work.
                        CALL FLTRXREF( CFIP, CDUM, SCCZERO, ' ', IDUM,
     &                                 IDUM, IDUM, LDUM, SKIPREC )
            
C.........................  Convert TSCC to internal value
                        CALL MBSCCADJ( IREC, TSCC, CRWT, CVID, 
     &                                 TSCC, EFLAG )
                    END SELECT

C.....................  Loop to iterate through all SRC/PROB pairs
                    DO J = 1, NPAIR
C.........................  Load SRC and PROB pair into FIELDARR( 1 )
C                           and  FIELDARR( 2 ), resp.
                        CALL LOAD_FIELDARR

C.........................  Make sure that the spatial surrogates code 
C                           is an integer
                        IF( .NOT. CHKINT( FIELDARR( 1 ) ) ) THEN
                            EFLAG = .TRUE.
                            WRITE( MESG,94010 ) 'ERROR: Spatial ' //
     &                             'surrogates code of pair ', J,
     &                             'is not an integer at line', IREC
                            CALL M3MESG( MESG )
                            EXIT

C.........................  Make sure that the spatial surrogate
C                           probability is a real number
                        ELSE IF( .NOT. CHKREAL( FIELDARR( 2 ) ) ) THEN
                            EFLAG = .TRUE.
                            WRITE( MESG,94010 ) 'ERROR: Spatial ' //
     &                             'surrogates code probability ' //
     &                             'of pair ', J,
     &                             'is not a real number at line', IREC
                            CALL M3MESG( MESG )
                            EXIT
                        ENDIF
          
C.........................  Convert surrogate code to an integer
                        ISRG = STR2INT( FIELDARR( 1 ) )

C.........................  Convert surrogate probability to a real
C                           number and make sure that the number is
C                           between 0.0 an 1.0
                        RPROB = STR2REAL( FIELDARR( 2 ) )
                        IF ( RPROB .LT. 0.0 .OR. RPROB .GT. 1.0 ) THEN

                            EFLAG = .TRUE.
                            WRITE( MESG,94010 ) 
     &                          'ERROR: at line', IREC, 
     &                          'of GRIDDING X-REF file, ' //
     &                          'probability value is not in (0,1)'
                            CALL M3MESG( MESG )
                            EXIT
                        ENDIF

                        ISRGCDA( N , J ) = ISRG
                        RSPROBA( N , J ) = RPROB

                    END DO  ! end J loop

C...................  Check for errors reading XREF file, and abort
                    IF( EFLAG ) THEN
                        MESG = 'Problem reading gridding cross-' //
     &                         'reference file.'
                        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                    ENDIF

c
C.....................  If this record is in error, or should be skipped because 
C                       it doesn't match any sources, go to next iteration
c                    IF( EFLAG .OR. SKIPREC ) CYCLE

           
C.....................  Store case-indpendent fields
                    INDXTA ( N ) = N
                    CSCCTA ( N ) = TSCC
                    INPAIRA ( N ) = NPAIR

C.....................  Store sorting criteria as right-justified in fields
                    CSRCALL = ' '
                    IF( CATEGORY .EQ. 'AREA' ) THEN
                        CALL BLDCSRC( CFIP, TSCC, CHRBLNK3, 
     &                                CHRBLNK3, CHRBLNK3, CHRBLNK3, 
     &                                CHRBLNK3, POLBLNK3, CSRCALL )
                    ELSE
                        CALL BLDCSRC( CFIP, RWTBLNK3, CHRBLNK3, 
     &                                CHRBLNK3, CHRBLNK3, CHRBLNK3, 
     &                                CHRBLNK3, POLBLNK3, CSRCALL )
                    ENDIF

                    CSRCTA( N ) = CSRCALL( 1:SC_ENDP( NCHARS ) ) // TSCC

                ENDIF  !  end of FIRSTIME IF block

            END DO  !  inner DO

 100        CONTINUE

C.............  Rewind file
            REWIND( FDEV )

            IF ( FIRSTIME ) THEN

                NLINES = N

C.................  Allocate memory for unsorted data used in all 
C                   source categories 
                ALLOCATE( INPAIRA( NLINES ), STAT=IOS )
                CALL CHECKMEM( IOS, 'INPAIRA', PROGNAME ) 
                ALLOCATE( CSCCTA( NLINES ), STAT=IOS )
                CALL CHECKMEM( IOS, 'CSCCTA', PROGNAME )
                ALLOCATE( CSRCTA( NLINES ), STAT=IOS )
                CALL CHECKMEM( IOS, 'CSRCTA', PROGNAME )
                ALLOCATE( INDXTA( NLINES ), STAT=IOS )
                CALL CHECKMEM( IOS, 'INDXTA', PROGNAME )
                ALLOCATE( ISRGCDA( NLINES, NSPMAX), STAT=IOS )
                CALL CHECKMEM( IOS, 'ISRGCDA', PROGNAME ) 
                ALLOCATE( RSPROBA( NLINES, NSPMAX), STAT=IOS )
                CALL CHECKMEM( IOS, 'RSPROBA', PROGNAME ) 

                ISRGCDA = IMISS3  !  array
                RSPROBA = AMISS3  !  array

                FIRSTIME = .FALSE.

            ELSE 

                EXIT  !  exit outer DO loop
 
            ENDIF

        END DO  ! outer DO


C.........  Reset number of cross-reference entries in case some were dropped
        NXREF = NLINES

C.........  Write errors for problems with input
        IF( NXREF .EQ. 0 ) THEN
            EFLAG = .TRUE.
            MESG = 'ERROR: No valid gridding cross-reference entries!'
            CALL M3MSG2( MESG )

        ELSE IF( NXREF .GT. NLINES ) THEN
            EFLAG = .TRUE.
            WRITE( MESG,94010 ) 'INTERNAL ERROR: dimension for ' //
     &             'storing gridding cross-reference was', NLINES,
     &             CRLF() // BLANK10 // 'but actually needed', NXREF
            CALL M3MSG2( MESG )

        ENDIF

C.......  Check for errors reading XREF file, and abort
        IF( EFLAG ) THEN
            MESG = 'Problem reading gridding cross-reference file.'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        ENDIF

        CALL M3MSG2( 'Processing gridding cross-reference file...' )

        CALL SORTIC( NXREF, INDXTA, CSRCTA )

C.........  Group cross-reference data into tables for different groups
        CALL XREFTBL( 'GRIDDING', NXREF )

C.........  Deallocate other temporary unsorted arrays
c        DEALLOCATE( ISRGCDA, CSRCTA, CSCCTA, INDXTA )
C......... Save ISRGCDA, RSPROBA, INPAIRA for later use 
        DEALLOCATE( CSRCTA, CSCCTA, INDXTA )

        RETURN


C.........  Error message for reaching the end of file too soon
999     MESG = 'End of file reached unexpectedly. ' //
     &         'Check format of gridding' // CRLF() // BLANK5 //
     &         'cross-reference file.'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

C******************  FORMAT  STATEMENTS   ******************************

C...........   Formatted file I/O formats............ 93xxx

93000   FORMAT( A )

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

C******************  INTERNAL SUBPROGRAMS  *****************************

        CONTAINS

C.............  This internal subprogram retrieves the input fields 
C               from LINE into LINEARR buffer.  Each successive call loads
C               FIELDARR with input field pairs from the surrogates file

            SUBROUTINE LOAD_FIELDARR

C.............  Local variables
            INTEGER          IOS2       ! I/O status
            INTEGER, SAVE :: PAIR       ! drives pairwise retrieval
            INTEGER          PGET1      ! low index to retrieve
            INTEGER          PGET2      ! high index to retrieve

C.............  Array of fields
            CHARACTER(LEN=SCCLEN3), ALLOCATABLE, SAVE :: LINEARR( : ) 

C----------------------------------------------------------------------

            IF ( .NOT. LOADED ) THEN

                PAIR = 0

                ALLOCATE( LINEARR( NCOLTOT ), STAT=IOS )
                CALL CHECKMEM( IOS, 'LINEARR', PROGNAME )

                CALL PARSLINE( LINE, NCOLTOT, LINEARR )
                LOADED = .TRUE.

            ENDIF

C.............  Get low and high indexes for LINEARR
            PGET1 = (PAIR * 2) + 1
            PGET2 = PGET1 + 1

            IF (PGET2 .GT. NCOLTOT) THEN

                FIELDARR( 1 ) = ' '
                FIELDARR( 2 ) = ' '

                EFLAG = .TRUE.
                WRITE( MESG,96010 ) 'INTERNAL ERROR: dimension for' //
     &                 ' storing gridding cross-reference was ', 
     &                 PGET1, ' and ', PGET2, ' but needed', NCOLTOT
                CALL M3MSG2( MESG )

                RETURN

            ENDIF

C.............  Retrieve specified field pair from LINEARR
            FIELDARR( 1 ) = LINEARR( PGET1 )  !  get SRC value
            FIELDARR( 2 ) = LINEARR( PGET2 )  !  get PROB value

C.............  Prepare for next data retrieval
            PAIR = PAIR + 1

C................. Remove LINEARR from memory
            IF (PGET2 .GE. NCOLTOT) DEALLOCATE( LINEARR )

            RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats............ 96xxx

96010   FORMAT( 10( A, :, I8, :, 1X ) )

            END SUBROUTINE LOAD_FIELDARR

        END SUBROUTINE RDGREF
