
        SUBROUTINE PARSLINE( LINE, N, SEGMENT )

C***********************************************************************
C  subroutine body starts at line 
C
C  DESCRIPTION:
C      This subroutine separates a "list-formatted" line of strings in which
C      the segments may or may not have quotes.  Although fortran requires
C      the quotes for true list-formatting, this subroutine can be used when
C      the quotes are only present to enclose a character (such as space, comma,
C      or semi-colon) that would otherwise be a delimiter.  If an "!" is 
C      encountered, everything after it is treated as a comment.
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C      Created by M. Houyoux 3/99
C
C****************************************************************************/
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
        CHARACTER(2) CRLF
        INTEGER      FINDC

        EXTERNAL     CRLF, FINDC

C...........   SUBROUTINE ARGUMENTS
        CHARACTER(*), INTENT (IN) :: LINE         ! character string to parse
        INTEGER     , INTENT (IN) :: N            ! maximum array length
        CHARACTER(*), INTENT(OUT) :: SEGMENT( N ) ! parsed string

C...........   Local parameters
        INTEGER,   PARAMETER :: NDELIM = 4
        CHARACTER, PARAMETER :: DELIMLST( NDELIM ) = 
     &                         (/ ',', ' ', ';', '	' /)
             
C...........   Array of 1-char strings for processing
        CHARACTER   ARRSTR( 5120 )  ! 256 * 20

C...........  Arrays for sorting non-delimiters on a per-machine basis
        INTEGER            NDINDX  ( NDELIM )
        CHARACTER, SAVE :: DELIMSRT( NDELIM )

C...........   Other local variables
        INTEGER         I, J, K, L, L1, L2  !  counters and indices
        INTEGER         IXP              !  index to non-delimeters
        INTEGER      :: NCNT             !  count of fields

        LOGICAL      :: ALPHA            !  true when within alpha-numeric 
        LOGICAL      :: DELIM            !  true when within or past delimiter 
        LOGICAL, SAVE:: FIRSTIME = .TRUE.!  true first time routine is called
        LOGICAL      :: PREVDELIM = .TRUE. !  true when last char was a delim
        LOGICAL      :: NOSPACDLIM = .FALSE. !  true when encountered a non-space delimiter since the previous field
        LOGICAL      :: NUMBER           !  true when within number in string 
        LOGICAL      :: QUOTED           !  true when within quotes in string
        LOGICAL      :: THISNMBR         !  true when current iteration is numbr

        CHARACTER       CBUF             !  temporary buffer
        CHARACTER    :: DOUBLEQ = '"'
        CHARACTER    :: SINGLEQ = "'"  
        CHARACTER    :: PERIOD  = '.' 
        CHARACTER    :: QUOTVAL = ' '    !  value of starting quote 

        CHARACTER(300)  MESG             ! message buffer

        CHARACTER(16) :: PROGNAME = 'PARSLINE' ! program name

C***********************************************************************
C   begin body of subroutine PARSLINE

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

        L2 = LEN_TRIM( LINE )

C.........  Check for comments, and use to set the end of the line
        L = INDEX( LINE( 1:L2 ), '!' )

        IF( L .LE. 0 ) THEN
            L = L2
        ELSE
            L = L - 1
        END IF

C.........  Skip blank lines
        IF( L .EQ. 0 ) RETURN

C.........  Initialize count, flags, and segments (npte, initializing in
C           the variable definitions is insufficient)
        NCNT    = 0
        SEGMENT = ' ' ! array
        ALPHA   = .FALSE.
        DELIM   = .TRUE.
        NUMBER  = .FALSE.
        QUOTED  = .FALSE.

C.........  Process LINE 1-character at a time
        DO I = 1, L

            CBUF = LINE( I:I )

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
                    NOSPACDLIM = .FALSE.
                    L1     = I + 1
                    NCNT    = NCNT + 1
                    K = 1

                ELSE IF( CBUF .EQ. DOUBLEQ ) THEN
                    QUOTED  = .TRUE.
                    DELIM   = .FALSE.
                    QUOTVAL = DOUBLEQ
                    PREVDELIM = .FALSE.
                    NOSPACDLIM = .FALSE.
                    L1      = I + 1
                    NCNT    = NCNT + 1
                    K = 1

                ELSE IF( ALPHA ) THEN
                    DELIM = .FALSE.
                    PREVDELIM = .FALSE.
                    NOSPACDLIM = .FALSE.
                    L1    = I
                    NCNT  = NCNT + 1
                    K = 1

                ELSE IF( NUMBER ) THEN
                    DELIM  = .FALSE.
                    PREVDELIM = .FALSE.
                    NOSPACDLIM = .FALSE.
                    L1     = I
                    NCNT   = NCNT + 1
                    K = 1

C...............  If hit another non-space delimiter without having
C                 hit another field with contents, then iterate the
C                 field count to create a blank space.
                ELSE IF( CBUF .NE. DELIMLST( 2 ) .AND.
     &                   NOSPACDLIM                   ) THEN
                    NCNT = NCNT + 1
                    PREVDELIM = .TRUE.
                    K = 1

C...............  If first col of line is a delimiter, fill a blank for the first col
                ELSE IF( I == 1 .AND. IXP .GT. 0 .AND.
     &                   CBUF .NE. DELIMLST( 2 )     ) THEN 
                    ALPHA = .FALSE.
                    NUMBER = .FALSE.
                    DELIM  = .TRUE.
                    PREVDELIM = .TRUE.
                    NOSPACDLIM = .TRUE.

                    NCNT = NCNT + 1
                    SEGMENT( NCNT ) = ' '

                END IF  ! Else its a space delimiter

C.............  In a quoted field, skip everything unless it is an end quote
            ELSE IF( QUOTED ) THEN

                IF( CBUF .EQ. QUOTVAL ) THEN
                    QUOTED  = .FALSE.
                    PREVDELIM = .FALSE.
                    K = 2
                END IF

C.............  If start of field was a number, but adjacent character is not
C               a delimiter, then turn field into an alpha
            ELSE IF( NUMBER .AND. .NOT. THISNMBR .AND. IXP .LE. 0 ) THEN
                ALPHA  = .TRUE.
                NUMBER = .FALSE.
                K = 1

C.............  If start of field was a number or alpha, and this char is a 
C               delimiter, then end of the field has been reached
            ELSE IF( IXP .GT. 0 ) THEN
                ALPHA = .FALSE.
                NUMBER = .FALSE.
                DELIM  = .TRUE.
                PREVDELIM = .TRUE.
                IF( CBUF .NE. DELIMLST( 2 ) ) NOSPACDLIM = .TRUE.
                L2     = I - K

                CALL STORE_SEGMENT

            END IF

        END DO

C.........  Store final segment
        IF( CBUF .EQ. QUOTVAL ) L = L - 1
        L2 = L

        IF( IXP .LE. 0 ) CALL STORE_SEGMENT

        RETURN

C******************  FORMAT  STATEMENTS   ************************************

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

C******************  INTERNAL SUBPROGRAMS  *****************************

        CONTAINS

C.............  This subprogram stores the segment from the input string
            SUBROUTINE STORE_SEGMENT

            IF( NCNT .LE. N ) THEN

                SEGMENT( NCNT ) = ADJUSTL( LINE( L1:L2 ) )

            ELSE

                MESG = 'ERROR: Overflow prevented while '//
     &                 'parsing line ' // PROGNAME
                CALL M3MSG2( MESG )
                MESG = 'First 200 characters of line contents are:'
                CALL M3MSG2( MESG )
                MESG = LINE( 1:200 )
                CALL M3MSG2( MESG )

                MESG = 'Formatting problem.'
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                
            END IF

            END SUBROUTINE STORE_SEGMENT

        END SUBROUTINE PARSLINE
