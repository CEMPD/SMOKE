
C Version "@(#)$Id$ $Source$ $Date$ 

        INTEGER FUNCTION PROMPTFFILE( PROMPT, 
     &                                RDONLY, FMTTED, DEFAULT, CALLER )

C***********************************************************************
C  function body starts at line 76
C
C       Prompts user for logical file name, then opens the Fortran file
C       associated with it, for read-only or not, and formatted or not,
C       as indicated by RDONLY and FMTTED
C
C  RETURNS:
C       unit number for the file opened, or 
C       -1 for failure, 
C       -2 for 'NONE', provided that '"NONE"' occurs within the prompt; or
C       -3 for 'ALL',  provided that '"ALL"'  occurs within the prompt
C
C  PRECONDITIONS REQUIRED:
C       "setenv <lname> <pathname>" for the file before program launch
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C       GETEFILE
C
C  REVISION  HISTORY:
C       prototype 6/95 by CJC
C	Revised  10/95 by CJC:  more robust treatment of 'NONE'
C       Modified  8/96 by CJC:  ! is a comment-designator for input
C***********************************************************************

      IMPLICIT NONE

C...........   ARGUMENTS and their descriptions:
        
        CHARACTER*(*) PROMPT         !  prompt for user
        LOGICAL       RDONLY         !  TRUE iff file is input-only
        LOGICAL       FMTTED         !  TRUE iff file should be formatted
        CHARACTER*(*) DEFAULT        !  default logical file name
        CHARACTER*(*) CALLER         !  caller-name for logging messages


C...........   EXTERNAL FUNCTIONS and their descriptions:

        INTEGER         GETEFILE, TRIMLEN
        LOGICAL         GETYN
        EXTERNAL        GETEFILE, GETYN, TRIMLEN


C...........   PARAMETER

        CHARACTER*16    BLANK16, NONE16, ALL16

        PARAMETER     ( BLANK16 = ' ' , 
     &                  NONE16 = 'NONE',
     &                  ALL16  = 'ALL' )


C...........   SCRATCH LOCAL VARIABLES and their descriptions:

        CHARACTER*16    LNAME        !  logical file name
        INTEGER         I            !  position at which "!" found
        INTEGER         IOS          !  I/O error status
        INTEGER         IDEV         !  unit number
        CHARACTER*80    PBUF         !  prompt buffer
        CHARACTER*120   MESG         !  messages
        LOGICAL         AFLAG        !  "ALL"  is in the prompt
        LOGICAL         NFLAG        !  "NONE" is in the prompt


C***********************************************************************
C   begin body of function  PROMPTFFILE

C.......   Get file name; open input control definition file

        AFLAG = ( INDEX( PROMPT, '"ALL"'  ) .GT. 0 )
        NFLAG = ( INDEX( PROMPT, '"NONE"' ) .GT. 0 )

        WRITE( PBUF,94000 ) 
     &          PROMPT ( 1: TRIMLEN( PROMPT  ) ) , ' [',
     &          DEFAULT( 1: TRIMLEN( DEFAULT ) ) , '] >> '

11      CONTINUE
        
            WRITE( *,95000 ) PBUF( 1:TRIMLEN( PBUF) + 1 )
            READ ( *,93010, IOSTAT=IOS ) LNAME
            
            IF ( IOS .GT. 0 ) THEN

                MESG = 'Could not read your response'
                WRITE( *,92000 ) MESG
                IF ( GETYN( 'Try again?', .TRUE. ) ) THEN
                    GO TO  11
                ELSE
                    MESG = 'Could not read logical name for file'
                    CALL M3EXIT( CALLER, 0, 0, MESG, 2 )
                END IF

            ELSE IF ( IOS .LT. 0 ) THEN

                MESG = 'Ending program ' //
     &                 CALLER( 1 : TRIMLEN( CALLER ) ) // '".'
                CALL M3EXIT( CALLER, 0, 0, MESG, 2 )

            END IF      !  if could not read response

            I = INDEX( LNAME, '!' )
            IF ( I .GT. 0 ) LNAME( I : LEN( LNAME ) ) = ' '

            IF ( LNAME .EQ. BLANK16 )  THEN
                LNAME = DEFAULT
            END IF

            IF ( AFLAG .AND. ( LNAME .EQ. ALL16 ) )  THEN
                PROMPTFFILE = -3
                RETURN
            ELSE IF ( NFLAG .AND. ( LNAME .EQ. NONE16 ) )  THEN
                PROMPTFFILE = -2
                RETURN
            END IF

            IDEV = GETEFILE( LNAME, RDONLY, FMTTED, CALLER )
            IF ( IDEV .LT. 0 ) THEN     !  failure to open

                MESG = 'Could not open input file "' //
     &                 LNAME( 1 : TRIMLEN( LNAME ) ) // '".'
                WRITE( *,92000 )  ' ', MESG
                IF ( GETYN( 'Try again?', .TRUE. ) ) THEN
                    GO TO  11
                ELSE
                    MESG = 'Ending program "' //
     &                      CALLER( 1 : TRIMLEN( CALLER ) ) // '".'
                    CALL M3EXIT( CALLER, 0, 0, MESG, 2 )
                END IF

            END IF      !  if getefile() failed

      PROMPTFFILE = IDEV
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

