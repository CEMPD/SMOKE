
C Version "@(#)$Id$ $Source$ $Date$ 

        SUBROUTINE  LUSTR( STRING )

C***********************************************************************
C  subroutine body starts at line
C
C  FUNCTION:  left-justify and upcase contents of STRING
C
C  PRECONDITIONS REQUIRED: none
C
C  SUBROUTINES AND FUNCTIONS CALLED: none
C
C  REVISION  HISTORY:
C       Prototype 6/95 by CJC
C
C***********************************************************************

      IMPLICIT NONE

C...........   ARGUMENTS and their descriptions:

        CHARACTER*(*)   STRING


C...........   PARAMETERS and their descriptions:

        CHARACTER*1     BLANK
        INTEGER         IDIFF
        
        PARAMETER     ( BLANK = ' ',
     &                  IDIFF = 65 - 97 ) ! = ichar( 'A' ) -  ichar( 'a' )


C...........   SCRATCH LOCAL VARIABLES and their descriptions:

        INTEGER         ILEN
        INTEGER         I, J, IBEG
        CHARACTER*1     CH


C***********************************************************************
C   begin body of subroutine  LUSTR


        ILEN = LEN( STRING )
        DO 10 I = 1,ILEN
            IF( STRING(I:I) .NE. BLANK ) THEN
                IBEG = I
                GO TO  11
            ENDIF
10      CONTINUE
           
C.......   If you get to here:  string all blanks.  No action necessary.
 
        RETURN
 
11      CONTINUE

C.......   Go thru rest of string, replacing lower-case by corresponding upper
 
        J = 1
        DO 20 I = IBEG, ILEN
            CH = STRING( I:I ) 
            IF ( CH .GE. 'a'  .AND.  CH .LE. 'z' ) THEN
                STRING(J:J) = CHAR( ICHAR( CH ) + IDIFF )
            ELSE
                STRING(J:J) = CH
            END IF
            J = J + 1
20      CONTINUE

C.......    pad trailing section of string with blank
        
        DO 30 I = J, ILEN
            STRING(I:I) = BLANK
30      CONTINUE

        RETURN
        END

