
        SUBROUTINE RDCODNAM( FDEV, MXDAT, NDAT, CODES, NAMES )

C***********************************************************************
C  subroutine body starts at line 92
C
C  DESCRIPTION:
C     Reads the master inventory pollutants list or activity list, concatenates
C     the names if they are longer than 16 characters, checks for duplicates,
C     renames if there are duplicates, and converts all names to uppercase. The 
C     subroutine returns the names in their original, unsorted order. If the
C     arrays are not empty on input, the subroutine appends the data from the
C     file to the input arrays.
C
C  PRECONDITIONS REQUIRED:
C     File unit FDEV already is opened
C     Memory allocated for CODES and NAMES
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C     Subroutines: Models-3 subroutines
C     Functions: Models-3 functions
C
C  REVISION  HISTORY:
C     Created 10/98 by M. Houyoux
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

        IMPLICIT NONE

C...........   INCLUDES

        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters

C...........   EXTERNAL FUNCTIONS and their descriptions:
        CHARACTER*2     CRLF
        INTEGER         LBLANK
        INTEGER         STR2INT
        INTEGER         TRIMLEN

        EXTERNAL CRLF, LBLANK, STR2INT, TRIMLEN

C...........   SUBROUTINE ARGUMENTS
        INTEGER     , INTENT    (IN) :: FDEV           ! input file unit no.
        INTEGER     , INTENT    (IN) :: MXDAT          ! max pols and/or actvties
        INTEGER     , INTENT(IN OUT) :: NDAT           ! no. cumulative recs
        INTEGER     , INTENT   (OUT) :: CODES( MXDAT ) ! pol or activity codes
        CHARACTER(*), INTENT   (OUT) :: NAMES( MXDAT ) ! pol or activity names

C...........   Unsorted pollutant records
        INTEGER       INDX1A( MXDAT )
        INTEGER       INDX2A( MXDAT )
        INTEGER       CODESA( MXDAT )
        CHARACTER(LEN=IOVLEN3) NAMESA( MXDAT )    !  pollutant names

C...........   Other local variables
        INTEGER         COD     !  tmp for pollutant code
        INTEGER         I, J    !  counters and indices
        INTEGER         IOS     !  i/o status
        INTEGER         IREC    !  record counter
        INTEGER         L, L2   !  length indices
        INTEGER         LCOD    !  previous pollutant code
        INTEGER         NDATSAV !  no. records in current file

        LOGICAL, SAVE :: FIRSTIME = .TRUE.

        CHARACTER*300   LNAM    !  previous pollutant name
        CHARACTER*300   PNAM    !  tmp for pollutant name (NOT correct length)
        CHARACTER*300   LINE    !  line buffer
        CHARACTER*300   MESG    !  message buffer

        CHARACTER*16 :: PROGNAME = 'RDCODNAM' ! program name

C***********************************************************************
C   begin body of subroutine RDCODNAM

C.........  Store initial value of counter
        NDATSAV = NDAT

C.........  Loop through input file...

        I    = 0
        IREC = 0
        DO

            READ( FDEV, 93000, END=22, IOSTAT=IOS ) LINE
            IREC = IREC + 1

            IF( IOS .GT. 0 ) THEN
                WRITE( MESG, 94010 ) 
     &                 'Error', IOS,  'reading names ' // 
     &                 'and codes file at line', IREC
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            ENDIF

C.............  Extract code and length of name from current line
            L2 = TRIMLEN( LINE )
            LINE = LINE( LBLANK( LINE ) + 1: L2 )
            COD = STR2INT( LINE( 1:5 ) )

            PNAM = ADJUSTL( LINE( 6:L2 ) )
            L = TRIMLEN( PNAM )

C.............  Truncate name to IOVLEN3 characters
            IF( L .GT. IOVLEN3 ) THEN

                PNAM = PNAM( 1:IOVLEN3 )

                WRITE( MESG,94010 )
     &                 'WARNING: Name too long for ' //
     &                 'I/O API variable names at line', IREC,
     &                 CRLF() // BLANK5 //
     &                 'in pollutant file. Truncating to "' //
     &                 PNAM( 1:IOVLEN3 ) // '"'
                CALL M3MESG( MESG )

            ELSEIF( COD .LE. 0 ) THEN
                WRITE( MESG,94010 ) 'WARNING: ' //
     &                 'Blank, alphabetic, or 0 data code at line',
     &                 IREC, CRLF() // BLANK10 //
     &                 'in code/names file. Skipping record.'
                CALL M3MESG( MESG )
                CYCLE

            ELSE                                        !  county-specific zone

                I = I + 1
                IF( I .LE. MXDAT ) THEN
                    INDX1A( I ) = I
                    INDX2A( I ) = I
                    CODESA( I ) = COD
                    NAMESA( I ) = PNAM( 1:IOVLEN3 )
                END IF

            ENDIF      !  if fip zero, or nn000, or not.

        ENDDO          !  end read loop

22      CONTINUE       !  exit from loop reading FDEV

        NDAT = I

        IF( NDAT .GT. MXDAT ) THEN
            WRITE( MESG,94010 ) 
     &             'ERROR: Number of pollutant records :', NDAT, 
     &             CRLF() // BLANK10 // 'Memory allocated :', MXDAT
            CALL M3MSG2( MESG )

            MESG = 'Insufficient memory allocated for ' //
     &             'codes/names file'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

        ELSE IF( NDAT .EQ. 0 ) THEN
            MESG ='No entries in codes/names file.'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

        ENDIF

C.........  Sort by number for checking for duplicates
        CALL SORTI1( NDAT, INDX1A, CODESA )

C.........  Check for duplicate numbers - delete duplicates
        LCOD = -9
        DO I = 1, NDAT

            J = INDX1A( I )
            COD = CODESA( J )

            IF( COD .EQ. LCOD ) THEN
                WRITE( MESG,94010 )
     &                 'WARNING: Duplicate code "', CODES( I ),
     &                 '" at line', J, 'in codes/names file.' //
     &                 CRLF() // BLANK5 // 'Skipping record.'
                CALL M3MESG( MESG )

                CODESA( J ) = 0   ! flag for skipping

            ENDIF

            LCOD = COD

        ENDDO

C.........  Sort by name for checking for duplicates
        CALL SORTIC( NDAT, INDX2A, NAMESA )

C.........  Check for duplicate names - rename
        LNAM = EMCMISS3
        DO I = 1, NDAT

            J = INDX2A( I )
            PNAM = NAMESA( J )

            IF( PNAM .EQ. LNAM .AND. CODESA( J ) .NE. 0 ) THEN

                L2 = TRIMLEN( PNAM )
                WRITE( MESG,94010 )
     &                 'WARNING: Duplicate name "' // PNAM( 1:L2 ) //
     &                 '" at line', J, 'in codes/names file.' //
     &                 CRLF() // BLANK5 // 'Skipping record.'
                CALL M3MESG( MESG )

                CODESA( J ) = 0   ! flag for skipping

            ENDIF

            LNAM = PNAM

        ENDDO

C.........  Store valid entries in original (unsorted) order and convert to 
C           uppercase
        J = NDATSAV
        DO I = 1, NDAT

            IF( CODESA( I ) .NE. 0 ) THEN
                J = J + 1
                CODES( J ) = CODESA( I )
                NAMES( J ) = NAMESA( I )
            ENDIF

        ENDDO

        NDAT    = J
        NDATSAV = J

C.........  Rewind file

        REWIND( FDEV )

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Formatted file I/O formats............ 93xxx

93000   FORMAT( A )

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

        END SUBROUTINE RDCODNAM
