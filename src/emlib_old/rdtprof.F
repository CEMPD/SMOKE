
C Version "@(#)$Id$ $Source$ $Date$ 

        INTEGER FUNCTION RDTPROF( DEV, HEAD, M, N, CODES, PFACS )

C***********************************************************************
C  function body starts at line  63
C
C  FUNCTION:
C       Reads the indicated segment from DEV, sorts it, and returns 
C       number of records read.        
C
C  PRECONDITIONS REQUIRED:
C       this segment of the input file sorted
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C       TRIMLEN, M3ERR
C
C  REVISION  HISTORY:
C       prototype 3/95 by CJC
C       Version  11/95 sorts the input.
C***********************************************************************

      IMPLICIT NONE

C...........   ARGUMENTS and their descriptions:
        
        INTEGER         DEV     !  unit number for profiles file
        CHARACTER*(*)   HEAD    !  name of this segment type ('MONTHLY', etc.)
        INTEGER         M       !  # of factors per code
        INTEGER         N       !  max number of codes
        INTEGER         CODES( N )      !  list of legal codes
        REAL            PFACS( M, N )   !  returned table of factors


C...........   EXTERNAL FUNCTIONS and their descriptions:
        
        INTEGER         STR2INT         !  read (unjustified) int from str
        INTEGER         TRIMLEN         !  string length w/o trailing blanks
        
        EXTERNAL        STR2INT, TRIMLEN


C...........   SAVED local variables

        INTEGER         REC     !  current line number in profiles file
        DATA            REC / 0 /


C...........   SCRATCH LOCAL VARIABLES and their descriptions:

        INTEGER         CNT
        INTEGER         ID, LID
        INTEGER         IOS
        INTEGER         SUM             !  input row-total
        INTEGER         I, J, K         !  subscripts
        INTEGER         W, WT( 1000 )   !  weight-factors
        REAL            DIV             !  scratch divisor

        CHARACTER*120   LINE    !  character-string buffer for profiles input
        CHARACTER*16    BUF     !  message buffer
        INTEGER         INDXA(     1000 )  !  list of subscripts
        INTEGER         CODEA(     1000 )  !  list of legal codes
        REAL            PFACA( 25, 1000 )  !  returned table of factors
        LOGICAL         EFLAG   !  input error flag
        LOGICAL         SORTED  !  is input in order already?

C***********************************************************************
C   begin body of subroutine  RDTPROF

        REC = REC + 1
        READ( DEV, 93000, IOSTAT=IOS ) LINE

        IF      ( IOS .LT. 0 ) THEN
            BUF = HEAD
            WRITE( *,91010 ) 
     &          'No "' // BUF( 1:TRIMLEN( BUF ) ) // 
     &          '" segment in TEMPORAL PROFILES file'
            RDTPROF = 0
            RETURN
        ELSE IF ( IOS .GT. 0 ) THEN
            WRITE( *,91010 ) 
     &          'Line number ', REC, 
     &          'I/O status  ', IOS,  ' '
            CALL M3ERR( 
     &           'RDTPROF', 0, 0, 
     &          'Error reading month header in TEMPORAL PROFILES file',
     &          .TRUE. )
        ELSE IF ( INDEX( LINE, HEAD ) .EQ. 0 ) THEN
            BUF = HEAD
            WRITE( *,91010 ) 
     &          'No "' // BUF( 1:TRIMLEN( BUF ) ) // 
     &          '" segment in TEMPORAL PROFILES file'
            RDTPROF = 0
            RETURN
        END IF
        
        LID    = -1
        CNT    =  0
        EFLAG  = .FALSE.
        SORTED = .TRUE.
        
11      CONTINUE        !  loop:  read this segment from file

            REC = REC + 1
            READ( DEV, 93000, IOSTAT=IOS ) LINE
            IF ( IOS .LT. 0 ) THEN
                WRITE( *,91010 ) 
     &              'Line number ', REC
                CALL M3WARN( 'RDTPROF', 0, 0, 
     &              'Unexpected END OF FILE in TEMPORAL PROFILES' )
                EFLAG = .TRUE.
                RDTPROF = CNT
                GO TO 100
            ELSE IF ( IOS .GT. 0 ) THEN
                WRITE( *,91010 ) 
     &              'Line number ', REC, 
     &              'I/O status  ', IOS
                CALL M3WARN( 'RDTPROF', 0, 0, 
     &              'Error reading line in TEMPORAL PROFILES file' )
                EFLAG = .TRUE.
                GO TO 11
            ELSE IF ( INDEX( LINE , '/END/' ) .NE. 0 ) THEN
                RDTPROF = CNT
                GO TO  100
            ELSE IF ( CNT .GE. N ) THEN
                WRITE( *,91010 ) 
     &              'Line number ', REC
                CALL M3WARN( 'RDTPROF', 0, 0, 
     &              'Maximum profile count exceeded.' )
                EFLAG = .TRUE.
                GO TO 11
            END IF

            IF ( CNT .GE. N ) THEN
                WRITE( *,91010 ) 
     &              'Line number  ', REC,
     &              'Profile count', CNT
                CALL M3WARN( 'RDTPROF', 0, 0, 
     &              'Maximum profile count exceeded.' )
                EFLAG = .TRUE.
                GO TO 11
            END IF

            ID     = STR2INT( LINE( 1:5 ) )
            SORTED = ( SORTED .AND. ( ID .GT. LID ) )
            LID    = ID
            CNT    = CNT + 1
            CODEA( CNT ) = ID
            INDXA( CNT ) = CNT

            J   = 6
            K   = 9	!  j:k spans 4 characters
            SUM = 0
            DO  22  I = 1, M
                W       = STR2INT( LINE( J:K ) )
                SUM     = SUM + W
                WT( I ) = W
                J = J + 4
                K = K + 4
22          CONTINUE

            DIV = 1.0 / FLOAT( SUM )
            DO  33  I = 1, M
                PFACA( I,CNT ) = DIV * FLOAT( WT( I ) )
33          CONTINUE

            GO TO  11   !  to head of read-segment loop
   
C...........   Codes read.  Check for errors, then sort:

100     CONTINUE
        
        IF ( EFLAG ) THEN
            CALL M3EXIT( 'RDTPROF', 0, 0, 
     &                   'Fatal errors encountered', 2 )
        END IF
        
        IF ( .NOT. SORTED ) CALL SORTI1( CNT, INDXA, CODEA )

        DO  122  I = 1, CNT

            K = INDXA( I )
            CODES( I ) = CODEA( K )

            DO  111 J = 1, M

               PFACS( J,I ) = PFACA( J,K )

111         CONTINUE

122     CONTINUE


C******************  FORMAT  STATEMENTS   ******************************

C...........   Error and warning message formats..... 91xxx

91010   FORMAT ( 5X , A, I6 )


C...........   Formatted file I/O formats............ 93xxx

93000   FORMAT( A )


C...........   Internal buffering formats............ 94xxx

94000   FORMAT( I5 )

94010   FORMAT( 25 I4 )

94020   FORMAT( 3 ( A, I4, :, 1X ) )


        END

