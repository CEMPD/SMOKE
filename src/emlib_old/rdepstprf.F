
        INTEGER FUNCTION RDEPSTPRF( DEV, HEAD, M, N, CODES, PFACS )

C***********************************************************************
C  function body starts at line  102
C
C  DESCRIPTION:
C       Reads the indicated segment from DEV, sorts it, and returns 
C       number of records read.  For EPS2 temporal profiles that do
C       not equal the sum on the far right of each record.   
C       For example, if the WEEKLY packet had a profile of
C          IWKCD   WKF1 WKF2 WKF3 WKF4 WKF5 WKF6 WKF7 WKTOT
C           99      1    1    1    1    1    0    0     7
C           
C       The WKTOT is not the sum of the WKF1-7.  Thus RDTPROF
C       will not work and the user must use this routine to
C       temporally allocated in the desired manner.  This is how
C       EPS2 actually allocates emissions.  If your WKF1-7 does
C       equal the WKTOT ( and the same for the MONTHLY and DIURNAL
C       packet cases ), then use the RDTPROF routine if so desired.
C             
C
C  PRECONDITIONS REQUIRED:
C       this segment of the input file sorted
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C       TRIMLEN, M3ERR
C
C  REVISION  HISTORY:
C       prototype 6/99 by JMV
C***********************************************************************
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
C****************************************************************************

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
        INTEGER         FSUM, SUM             !  input row-total
        INTEGER         I, J, K         !  subscripts
        INTEGER         W, WT( 1000 )   !  weight-factors
        REAL            DIV             !  scratch divisor

        CHARACTER*120   LINE    !  character-string buffer for profiles input
        CHARACTER*256   MESG 
        CHARACTER*16    BUF     !  message buffer
        INTEGER         INDXA(     1000 )  !  list of subscripts
        INTEGER         CODEA(     1000 )  !  list of legal codes
        REAL            PFACA( 25, 1000 )  !  returned table of factors
        LOGICAL         EFLAG   !  input error flag
        LOGICAL         SORTED  !  is input in order already?

C***********************************************************************
C   begin body of subroutine  RDEPSTPRF

        REC = REC + 1
        READ( DEV, 93000, IOSTAT=IOS ) LINE

        IF      ( IOS .LT. 0 ) THEN
            BUF = HEAD
            WRITE( *,91010 ) 
     &          'No "' // BUF( 1:TRIMLEN( BUF ) ) // 
     &          '" segment in TEMPORAL PROFILES file'
            RDEPSTPRF = 0
            RETURN
        ELSE IF ( IOS .GT. 0 ) THEN
            WRITE( *,91010 ) 
     &          'Line number ', REC, 
     &          'I/O status  ', IOS,  ' '
            CALL M3ERR( 
     &           'RDEPSTPRF', 0, 0, 
     &          'Error reading month header in TEMPORAL PROFILES file',
     &          .TRUE. )
        ELSE IF ( INDEX( LINE, HEAD ) .EQ. 0 ) THEN
            BUF = HEAD
            WRITE( *,91010 ) 
     &          'No "' // BUF( 1:TRIMLEN( BUF ) ) // 
     &          '" segment in TEMPORAL PROFILES file'
            RDEPSTPRF = 0
            RETURN
        END IF

        LID    = -1
        CNT    =  0
        EFLAG  = .FALSE.
        SORTED = .TRUE.

11      CONTINUE        !  loop:  read this segment from file

            REC = REC + 1
            READ( DEV, 93000, IOSTAT=IOS ) LINE
            IF ( IOS .NE. 0 ) THEN
                GO TO 100
            ELSE IF ( INDEX( LINE , '/END/' ) .NE. 0 ) THEN
                GO TO  101
            END IF

            CNT = CNT + 1
            IF ( CNT .LE. N ) THEN

                ID     = STR2INT( LINE( 1:5 ) )
                SORTED = ( SORTED .AND. ( ID .GT. LID ) )
                LID    = ID
                CODEA( CNT ) = ID
                INDXA( CNT ) = CNT

                J   = 6
                K   = 9	!  j:k spans 4 characters
                SUM = 0
                DO  22  I = 1, M
                    W       = STR2INT( LINE( J:K ) )
                    SUM     = SUM + W
                    WT( I ) = W
                    J       = J + 4
                    K       = K + 4
22              CONTINUE

                K = K + 2

                FSUM = STR2INT( LINE( J:K ) )
                IF ( SUM .NE. FSUM ) THEN
                   WRITE( MESG, 94040 )
     &              'WARNING: Profile ',ID,'  does not sum correctly',
     &              ' for packet ', HEAD
                   CALL M3MESG( MESG )

                ENDIF
 

                DIV = 1.0 / FLOAT( FSUM )                  
                DO  33  I = 1, M
                    PFACA( I,CNT ) = DIV * FLOAT( WT( I ) )
33              CONTINUE

            END IF

            GO TO  11   !  to head of read-segment loop

100     CONTINUE


C...........   Codes read.  Check for errors, then sort:

        IF ( IOS .LT. 0 ) THEN
            WRITE( *,91010 )  'Line number ', REC
            CALL M3EXIT( 'RDEPSTPRF', 0, 0, 
     &          'Unexpected END OF FILE in TEMPORAL PROFILES' , 2 )
        ELSE IF ( IOS .GT. 0 ) THEN
            WRITE( *,91010 ) 
     &          'Line number ', REC, 
     &          'I/O status  ', IOS
            CALL M3EXIT( 'RDEPSTPRF', 0, 0, 
     &          'Error reading line in TEMPORAL PROFILES file' , 2 )
        END IF

101     CONTINUE

        IF ( CNT .GT. N ) THEN
            WRITE( LINE,94030 ) 
     &          'TPROF overflow--dim:', N,
     &          'actual:', CNT, HEAD
            CALL M3EXIT( 'RDEPSTPRF', 0, 0, LINE, 2 )
        END IF

        IF ( .NOT. SORTED ) CALL SORTI1( CNT, INDXA, CODEA )

        DO  122  I = 1, CNT

            K = INDXA( I )
            CODES( I ) = CODEA( K )

            DO  111 J = 1, M
               PFACS( J,I ) = PFACA( J,K )
111         CONTINUE

122     CONTINUE
        
        RDEPSTPRF = CNT


C******************  FORMAT  STATEMENTS   ******************************

C...........   Error and warning message formats..... 91xxx

91010   FORMAT ( 5X , A, I6 )


C...........   Formatted file I/O formats............ 93xxx

93000   FORMAT( A )


C...........   Internal buffering formats............ 94xxx

94000   FORMAT( I5 )

94010   FORMAT( 25 I4 )

94020   FORMAT( 3 ( A, I4, :, 1X ) )

94030   FORMAT( 10 ( A, :, I7, :, 1X ) )

94040   FORMAT( A, I5, 3 ( A ) )

        END

