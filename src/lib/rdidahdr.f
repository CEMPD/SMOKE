
        SUBROUTINE RDIDAHDR( FDEV, JREC, FILTYP, INVYR, MPOL )

C***********************************************************************
C  subroutine body starts at line 
C
C  DESCRIPTION:
C      This subroutine reads a IDA header record.
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C      Created by J. Vukovich 1/99
C
C****************************************************************************/

        IMPLICIT NONE

C...........   INCLUDES

        INCLUDE 'EMCNST3.EXT'   !  emissions constat parameters

C...........   EXTERNAL FUNCTIONS and their descriptions:
        
        INTEGER         GETNLIST
        INTEGER         TRIMLEN

        EXTERNAL        GETNLIST, TRIMLEN

C...........   SUBROUTINE ARGUMENTS
        INTEGER       FDEV        !  unit number of input file (in) 
        CHARACTER*(*) FILTYP      !  description of source category (in) 
        INTEGER       INVYR       !  inventory year (out)  
        INTEGER       MPOL        !  number of pollutants (out) 

C...........   Other local variables
        INTEGER         I, J, L   !  counters and indices

        INTEGER         IOS         !  i/o status
        INTEGER      :: IREC    = 0 !  input line counter
        INTEGER         FILFMT      !  file format code
        INTEGER      :: NLINES  = 0 !  number of lines
        INTEGER      :: NLINEBP = 0 !  number of lines times pollutants
        INTEGER      :: NPOL    = 0 !  number of pollutants at line in file

        CHARACTER*300   BUFFER      !  temporary buffer
        CHARACTER*300   LINE        !  input file line buffer
        CHARACTER*300   MESG        !  message buffer

        CHARACTER*16 :: PROGNAME = 'RDIDAHDR' ! program name

C***********************************************************************
C   begin body of function RDIDAHDR

        DO  J = 1, 5 ! 5 lines in header 

            JLNUM = JREC + J

C.............  Read in header line 
            READ( FDEV,93000, END=111, IOSTAT=IOS ) LINE
 
C.............  Check I/O error status
            IF( IOS .GT. 0 ) THEN
                WRITE( MESG, 94010 )
     &                 'Error', IOS,  'reading IDA input file ' //
     &                 'as character strings at line', IREC
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

            ENDIF

C.............  Scan to make sure its a header line
            IF( LINE( 1:1 ) .EQ. '#' ) THEN

C.................  Scan for pollutant header field
                L = TRIMLEN( LINE )
                ITYP = INDEX1( LINE( 2:L ), MXHDRTYP, HDRTYPES )

C.................  Make sure header type is valid
                IF ( ITYP .LE. 0 ) THEN 
                   WRITE( MESG,94010 )  'Error header type : ' //
     &                 "' // LINE( 2:L) //
     &                 '" in emission file at line', JLNUM,
     &                 CRLF() // BLANK5 //
     &                 'is not in master header list'
                   CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                ENDIF

C.................  
                IF( ITYP .EQ. 1 ) THEN
                    I = I + 5
                    BUFFER = LINE( I:L )
                    L = L - I - 1
  
                    CALL UPCASE( BUFFER )
                    NPOL = GETNLIST( L, BUFFER )

                ENDIF

                CYCLE  ! to end of loop

C.............  Check to ensure header was there!
            ELSEIF( NPOL .EQ. 0 ) THEN

                WRITE( MESG,94010 ) 'No #POLID header in IDA file', FDEV
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

C.............  Otherwise, count the lines and lines times pollutants
            ELSE

                NLINES = NLINES + 1
                NLINEBP = NLINEBP + NPOL

            ENDIF

        ENDDO

111     CONTINUE  ! Exit from read loop

        REWIND( FDEV )

        IF( OUTTYPE .EQ. 1 ) THEN
            GETIDASZ = NLINES

        ELSEIF( OUTTYPE .EQ. 2 ) THEN
            GETIDASZ = NLINEBP

        ELSE
            MESG = 'INTERNAL ERROR: Bad output type in call to ' // 
     &             'subroutine ' // PROGNAME     
            CALL M3MSG2( MESG )
            CALL M3EXIT( PROGNAME, 0, 0, ' ', 2 )

        ENDIF

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Formatted file I/O formats............ 93xxx

93000   FORMAT( A )

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

        END
