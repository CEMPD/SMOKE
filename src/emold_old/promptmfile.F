
C Version "@(#)$Id$ $Source$ $Date$ 

        CHARACTER*16 FUNCTION PROMPTMFILE( PROMPT, FMODE, 
     &                                     DEFAULT, CALLER )

C***********************************************************************
C  function body starts at line 80
C
C       Prompts user for logical file name, then opens the I/O API file
C       associated with it, using the indicated file mode (FSREAD3, 
C       FSRDWR3, FSNEW3, FSUNKN3)
C
C       Provided that '"NONE"'occurs within the prompt, if name entered 
C       is 'NONE', does not attempt to open the file (but still returns 
C       'NONE' as Dthe function value).
C
C  RETURNS:
C       logical name of file opened
C
C  PRECONDITIONS REQUIRED:
C       "setenv <lname> <pathname>" for the file before program launch
C       file description set in FDESC3.EXT structures if appropriate
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C       GETYN, TRIMLEN, OPEN3
C
C  REVISION  HISTORY:
C       prototype 6/95 by CJC
C	Revised  10/95 by CJC:  more robust treatment of 'NONE'
C       Modified  8/96 by CJC:  ! is a comment-designator for input
C***********************************************************************

      IMPLICIT NONE

        INCLUDE 'PARMS3.EXT'
        INCLUDE 'IODECL3.EXT'


C...........   ARGUMENTS and their descriptions:
        
        CHARACTER*(*) PROMPT         !  prompt for user
        INTEGER       FMODE          !  file opening-mode
        CHARACTER*(*) DEFAULT        !  default logical file name
        CHARACTER*(*) CALLER         !  caller-name for logging messages


C...........   PARAMETER

        CHARACTER*16    BLANK16
        CHARACTER*16    NONE16

        PARAMETER     ( BLANK16 = ' ' , 
     &                  NONE16  = 'NONE' )


C...........   EXTERNAL FUNCTIONS and their descriptions:

        INTEGER         TRIMLEN
        LOGICAL         ENVYN
        LOGICAL         GETYN
        EXTERNAL        ENVYN, GETYN, TRIMLEN


C...........   SCRATCH LOCAL VARIABLES and their descriptions:

        INTEGER         I            !  position at which "!" found
        CHARACTER*16    BUFFER       !  buffer
        CHARACTER*16    LNAME        !  logical file name
        INTEGER         LUNIT        !  log-file unit number
        INTEGER         IOS          !  I/O status
        CHARACTER*120   MESG         !  messages
        CHARACTER*80    PBUF         !  prompt buffer
        LOGICAL         NFLAG        !  "NONE" is in the prompt
        LOGICAL         PROMPTON     !  Actually prompt or open default
 
        LOGICAL  FIRSTIME     !  true if first call this execution
        DATA     FIRSTIME / .TRUE. /

        SAVE     FIRSTIME, PROMPTON
 
C***********************************************************************
C   begin body of function  PROMPTMFILE

        IF( FIRSTIME ) THEN
 
            PROMPTON = ENVYN( 'PROMPTFLAG', 'Prompt for input flag',
     &                      .TRUE., IOS )
            FIRSTIME = .FALSE.
 
        ENDIF

C.......   Decide whether 'NONE' is a valid response

        NFLAG = ( INDEX( PROMPT, '"NONE"' ) .GT. 0 )

        IF( PROMPTON ) THEN
C.......   Construct actual prompt; Loop:  get file name until file opens
        
            WRITE( PBUF,94000 ) 
     &              PROMPT ( 1: TRIMLEN( PROMPT  ) ) , ' [',
     &              DEFAULT( 1: TRIMLEN( DEFAULT ) ) , '] >> '

11          CONTINUE
        
                LNAME = ' '
                WRITE( *,95000 ) PBUF( 1:TRIMLEN( PBUF) + 1 )
                READ ( *,93010, IOSTAT=IOS ) LNAME
            
                IF ( IOS .GT. 0 ) THEN

                    MESG = 'Could not read your response'
                    WRITE( *,92000 ) MESG
                    IF ( GETYN( 'Try again?', .TRUE. ) ) THEN
                        GO TO  11
                    ELSE
                        MESG = 'Could not read logical name for file'
                        CALL M3EXIT( 'CTLAMAT', 0, 0, MESG, 2 )
                    END IF

                ELSE IF ( IOS .LT. 0 ) THEN

                    MESG = 'Ending program ' //
     &                     CALLER( 1 : TRIMLEN( CALLER ) ) // '".'
                    CALL M3EXIT( CALLER, 0, 0, MESG, 2 )

                END IF      !  if could not read response

                I = INDEX( LNAME, '!' )
                IF ( I .GT. 0 ) LNAME( I : LEN( LNAME ) ) = ' '

                IF ( LNAME .EQ. BLANK16 ) THEN

                    LNAME = DEFAULT

                END IF

                IF ( NFLAG .AND. ( LNAME .EQ. NONE16 ) ) THEN

                    PROMPTMFILE = NONE16
                    RETURN

                END IF

                LUNIT = INIT3()
            
                IF ( .NOT. OPEN3( LNAME, FMODE, CALLER ) ) THEN !  failure to open

                    MESG = 'Could not open file "' //
     &                     LNAME( 1 : TRIMLEN( LNAME ) ) // '".'
                    WRITE( *,92000 )  ' ', MESG
                    IF ( GETYN( 'Try again?', .TRUE. ) ) THEN
                        GO TO  11
                    ELSE
                        MESG = 'Ending program "' //
     &                          CALLER( 1 : TRIMLEN( CALLER ) ) // '".'
                        CALL M3EXIT( CALLER, 0, 0, MESG, 2 )
                    END IF

                END IF      !  if open3() failed

        ELSE  ! Do not prompt for output

            LNAME = DEFAULT

            IF ( NFLAG )  THEN
 
C           ..  Check if logical name is set in order to permit
C           ..  Study Planner to skip file without having to input "NONE"
 
                CALL ENVSTR( LNAME, 'Input file name', ' ',
     &                       BUFFER, IOS )
 
                IF( IOS   .EQ. -1     .OR.
     &              IOS   .EQ. -2     .OR.
     &              LNAME .EQ. NONE16      ) THEN
 
                    PROMPTMFILE = NONE16
                    RETURN

                ENDIF
            END IF

            LUNIT = INIT3()
 
            IF ( .NOT. OPEN3( LNAME, FMODE, CALLER ) ) THEN !  failure to open
 
                MESG = 'Could not open file "' //
     &                 LNAME( 1 : TRIMLEN( LNAME ) ) // '".'
                WRITE( *,92000 )  ' ', MESG
                MESG = 'Ending program "' //
     &                  CALLER( 1 : TRIMLEN( CALLER ) ) // '".'
                CALL M3EXIT( CALLER, 0, 0, MESG, 2 )

            END IF      !  if open3() failed

        ENDIF

        PROMPTMFILE = LNAME
        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Informational (LOG) message formats... 92xxx
 
92000   FORMAT( /5X, A )


C...........   Formatted file I/O formats............ 93xxx
 
93010   FORMAT( A16 )


C...........   Internal buffering formats ............94xxx
 
94000   FORMAT( 12( A, : ) )
 

C...........   Miscellaneous formats................. 95xxx

95000   FORMAT ( /5X , A , $ )          !  generic prompt format.


        END

