
        SUBROUTINE RDSIPOLS( FDEV, NPOLS, CODES, NAMES )

C***********************************************************************
C  subroutine body starts at line 92
C
C  DESCRIPTION:
C     Reads the master inventory pollutants list, concatenates the names if
C     they are longer than 16 characters, checks for duplicates, renames if
C     there are duplicates, and converts all names to uppercase. The subroutine
C     returns the names in their original, unsorted order.
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
        INTEGER       FDEV              !  pollutants file unit number
        INTEGER       NPOLS             !  number of pollutants in file
        INTEGER       CODES( NPOLS )    !  pollutant codes
        CHARACTER(LEN=IOVLEN3) NAMES( NPOLS )    !  pollutant names

C...........   Unsorted pollutant records
        INTEGER       INDX1A( NPOLS )
        INTEGER       INDX2A( NPOLS )
        INTEGER       CODESA( NPOLS )
        CHARACTER(LEN=IOVLEN3) NAMESA( NPOLS )    !  pollutant names

C...........   Other local variables
        INTEGER         COD     !  tmp for pollutant code
        INTEGER         I, J    !  counters and indices
        INTEGER         IOS     !  i/o status
        INTEGER         IREC    !  record counter
        INTEGER         L, L2   !  length indices
        INTEGER         LCOD    !  previous pollutant code

        CHARACTER*300   LNAM    !  previous pollutant name
        CHARACTER*300   PNAM    !  tmp for pollutant name (NOT correct length)
        CHARACTER*300   LINE    !  line buffer
        CHARACTER*300   MESG    !  message buffer

        CHARACTER*16 :: PROGNAME = 'RDSIPOLS' ! program name

C***********************************************************************
C   begin body of subroutine RDSIPOLS

        I    = 0
        IREC = 0
        DO

            READ( FDEV, 93000, END=22, IOSTAT=IOS ) LINE
            IREC = IREC + 1

            IF( IOS .NE. 0 ) THEN
                WRITE( MESG, 94010 ) 
     &                 'Error', IOS,  'reading pollutant names ' // 
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
     &                 'WARNING: Pollutant name too long for ' //
     &                 'I/O API variable names at line', IREC,
     &                 CRLF() // BLANK5 //
     &                 'in pollutant file. Truncating to "' //
     &                 PNAM( 1:IOVLEN3 ) // '"'
                CALL M3MESG( MESG )

            ELSEIF( COD .LE. 0 ) THEN
                WRITE( MESG,94010 )
     &                 'Blank, alphabetic, or 0 pollutant code at line',
     &                 IREC, 'in pollutant file. Skipping record.'
                CALL M3MESG( MESG )
                CYCLE

            ELSE                                        !  county-specific zone

                I = I + 1
                IF( I .LE. NPOLS ) THEN
                    INDX1A( I ) = I
                    INDX2A( I ) = I
                    CODESA( I ) = COD
                    NAMESA( I ) = PNAM( 1:IOVLEN3 )
                END IF

            ENDIF      !  if fip zero, or nn000, or not.

        ENDDO          !  end read loop

22      CONTINUE       !  exit from loop reading FDEV

        IF( I .GT. NPOLS ) THEN
            WRITE( MESG,94010 ) 
     &        'Number of pollutant records :', I, CRLF()// BLANK5//
     &        '           Memory allocated :', NPOLS
            CALL M3MSG2( MESG )

            MESG = 'ERROR: Insufficient memory allocated for ' //
     &             'pollutant codes and names file'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        ENDIF

        NPOLS = I

C.........  Sort by number for checking for duplicates
        CALL SORTI1( NPOLS, INDX1A, CODESA )

C.........  Check for duplicate numbers - delete duplicates
        LCOD = -9
        DO I = 1, NPOLS

            J = INDX1A( I )
            COD = CODESA( J )

            IF( COD .EQ. LCOD ) THEN
                WRITE( MESG,94010 )
     &                 'Duplicate pollutant code "', CODES( I ),
     &                 '" at line', J, 'in pollutant file.' //
     &                 CRLF() // BLANK5 // 'Skipping record.'
                CALL M3MESG( MESG )

                CODESA( J ) = 0   ! flag for skipping

            ENDIF

            LCOD = COD

        ENDDO

C.........  Sort by name for checking for duplicates
        CALL SORTIC( NPOLS, INDX2A, NAMESA )

C.........  Check for duplicate names - rename
        LNAM = EMCMISS3
        DO I = 1, NPOLS

            J = INDX2A( I )
            PNAM = NAMESA( J )

            IF( PNAM .EQ. LNAM .AND. CODESA( J ) .NE. 0 ) THEN

                L2 = TRIMLEN( PNAM )
                WRITE( MESG,94010 )
     &                 'Duplicate pollutant name "' // PNAM( 1:L2 ) //
     &                 '" at line', J, 'in pollutant file.' //
     &                 CRLF() // BLANK5 // 'Skipping record.'
                CALL M3MESG( MESG )

                CODESA( J ) = 0   ! flag for skipping

            ENDIF

            LNAM = PNAM

        ENDDO

C.........  Store valid entries in original (unsorted) order and convert to 
C           uppercase
        J = 0
        DO I = 1, NPOLS

            IF( CODESA( I ) .NE. 0 ) THEN
                J = J + 1
                CODES( J ) = CODESA( I )
                NAMES( J ) = NAMESA( I )
            ENDIF

        ENDDO

        NPOLS = J

C.........  Rewind file

        REWIND( FDEV )

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Formatted file I/O formats............ 93xxx

93000   FORMAT( A )

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

        END
