
        INTEGER FUNCTION GETNLIST( ILENGTH, STRING )

C***********************************************************************
C  function body starts at line 
C
C  DESCRIPTION:
C      This function counts the number of free-formatted strings in a list
C      of string that may or may not have quotes.  This is used to help
C      when a string is available for reading a series of values, but
C      no indication is available for the number of entries to read.  It
C      accounts for comments appended with a "!" and skips blank lines.
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C      Created by M. Houyoux 1/99
C
C**************************************************************************
C
C Project Title: EDSS Tools Library
C File: @(#)$Id$
C
C COPYRIGHT (C) 2004, Environmental Modeling for Policy Development
C All Rights Reserved
C
C Carolina Environmental Program
C University of North Carolina at Chapel Hill
C 137 E. Franklin St., CB# 6116
C Chapel Hill, NC 27599-6116
C
C smoke@unc.edu
C
C Pathname: $Source$
C Last updated: $Date$ 
C
C***************************************************************************

        IMPLICIT NONE

C...........   EXTERNAL FUNCTIONS
        INTEGER    FINDC
        EXTERNAL   FINDC

C...........   SUBROUTINE ARGUMENTS
        INTEGER       ILENGTH     !  length of string
        CHARACTER(*)  STRING      !  description of source category

C...........   Local parameters
        INTEGER,   PARAMETER :: NDELIM = 4
        CHARACTER, PARAMETER :: DELIMLST( NDELIM ) = 
     &                         (/ ',', ' ', ';', '	' /)
             
C...........   Array of 1-char strings for processing
        CHARACTER     ARRSTR( 5120 )  ! 256 * 20

C...........  Arrays for sorting non-delimiters on a per-machine basis
        INTEGER            NDINDX  ( NDELIM )
        CHARACTER, SAVE :: DELIMSRT( NDELIM )

C...........   Other local variables
        INTEGER         I, J, L, L1, L2  !  counters and indices
        INTEGER         IXP              !  index to non-delimeters
        INTEGER      :: NCNT             !  count of fields

        LOGICAL      :: ALPHA            !  true when within alpha-numeric 
        LOGICAL      :: DELIM            !  true when within or past delimiter 
        LOGICAL      :: FIRSTIME = .TRUE.!  true first time routine is called
        LOGICAL      :: PREVDELIM = .TRUE. !  true when last char was a delim
        LOGICAL      :: NUMBER           !  true when within number in string 
        LOGICAL      :: QUOTED           !  true when within quotes in string
        LOGICAL      :: THISNMBR         !  true when current iteration is numbr

        CHARACTER       CBUF             !  temporary buffer
        CHARACTER    :: DOUBLEQ = '"'
        CHARACTER    :: SINGLEQ = "'"  
        CHARACTER    :: PERIOD  = '.' 
        CHARACTER       QUOTVAL          !  value of starting quote 

        CHARACTER(300)  MESG             ! message buffer

        CHARACTER(16) :: PROGNAME = 'GETNLIST' ! program name

C***********************************************************************
C   begin body of function GETNLIST

C.........  The first time the routine is called, sort the list of delimiters
        IF( FIRSTIME ) THEN
            DO I = 1, NDELIM 
                NDINDX( I ) = I
            END DO

            CALL SORTIC( NDELIM, NDINDX, DELIMLST )

            DO I = 1, NDELIM 
                J = NDINDX( I )
                DELIMSRT( I ) = DELIMLST( J )
            END DO

            FIRSTIME = .FALSE.

        END IF

        L2 = LEN_TRIM( STRING )

C.........  Check for comments, and use to set the end of the line
        L = INDEX( STRING( 1:L2 ), '!' )

        IF( L .LE. 0 ) THEN
            L = L2
        ELSE
            L = L - 1
        END IF

C.........  Skip blank lines
        IF( L .EQ. 0 ) THEN
            GETNLIST = 0
            RETURN
        END IF

C.........  Initialize count and flags
        NCNT    = 0
        ALPHA   = .FALSE.
        DELIM   = .TRUE.
        NUMBER  = .FALSE.
        QUOTED  = .FALSE.

C.........  Process LINE 1-character at a time
        DO I = 1, L

            CBUF = STRING( I:I )

C.............  Look for character in delimiters
            IXP = FINDC( CBUF, NDELIM, DELIMSRT )

C.............  Evaluate the current character for number or not
            THISNMBR = ( CBUF .GE. '0' .AND. CBUF .LE. '9' )

C.............  Waiting for next field...
            IF( DELIM ) THEN

                NUMBER = THISNMBR
                ALPHA  = ( .NOT. NUMBER .AND. IXP .LE. 0 )

                IF( CBUF .EQ. SINGLEQ ) THEN
                    QUOTED  = .TRUE.
                    DELIM   = .FALSE.
                    QUOTVAL = SINGLEQ
                    PREVDELIM = .FALSE.
                    L1     = I + 1
                    NCNT    = NCNT + 1

                ELSE IF( CBUF .EQ. DOUBLEQ ) THEN
                    QUOTED  = .TRUE.
                    DELIM   = .FALSE.
                    QUOTVAL = DOUBLEQ
                    PREVDELIM = .FALSE.
                    L1      = I + 1
                    NCNT    = NCNT + 1

                ELSE IF( ALPHA ) THEN
                    DELIM = .FALSE.
                    PREVDELIM = .FALSE.
                    L1    = I
                    NCNT  = NCNT + 1

                ELSE IF( NUMBER ) THEN
                    DELIM  = .FALSE.
                    PREVDELIM = .FALSE.
                    L1     = I
                    NCNT   = NCNT + 1

C...............  If another delimeter, then another field, but last
C                 field was blank UNLESS delim is a space
                ELSE IF( CBUF .NE. DELIMLST( 2 ) ) THEN
                    
                    IF( PREVDELIM ) THEN
                        NCNT = NCNT + 1
                    ELSE
                        PREVDELIM = .TRUE.
                    END IF

                END IF  ! Else its a space delimiter

C.............  In a quoted field, skip everything unless it is an end quote
            ELSE IF( QUOTED ) THEN

                IF( CBUF .EQ. QUOTVAL ) THEN
                    QUOTED  = .FALSE.
                    DELIM   = .TRUE.
                    PREVDELIM = .FALSE.
                    L2      = I - 1
                  
                END IF

C.............  If start of field was a number, but adjacent character is not
C               a delimiter, then turn field into an alpha
            ELSE IF( NUMBER .AND. .NOT. THISNMBR .AND. IXP .LE. 0 ) THEN
                ALPHA  = .TRUE.
                NUMBER = .FALSE.

C.............  If start of field was a number or alpha, and this is a 
C               delimiter, then end of number has been reached
            ELSE IF( IXP .GT. 0 ) THEN
                ALPHA = .FALSE.
                NUMBER = .FALSE.
                DELIM  = .TRUE.
                PREVDELIM = .TRUE.
                L2     = I - 1

            END IF

        END DO

        GETNLIST = NCNT

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

        END FUNCTION GETNLIST
