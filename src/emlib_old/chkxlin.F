
        INTEGER FUNCTION CHKXLIN( LINE, IOS, LINENO, PTYPE )

C***********************************************************************
C  subroutine body starts at line 75
C
C  FUNCTION: Checks a line of the temporal Xref file to see if needed 
C            for type of processing (area, point, or mobile).
C
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:  STR2INT, M3MESG, M3EXIT
C
C  REVISION  HISTORY:
C       prototype 4/96 by M Houyoux
C
C***********************************************************************
C
C Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
C                System
C File: @(#)$Id$
C
C COPYRIGHT (C) 1996, MCNC--North Carolina Supercomputing Center
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

C...........   INCLUDES:
        INCLUDE 'PARMS3.EXT'

C...........   ARGUMENTS and their descriptions:

        INTEGER       IOS              ! I/O status from ASCII read
        INTEGER       LINENO           ! Number of current line
        CHARACTER*(*) LINE             ! ASCII LINE from X-ref file
        CHARACTER*(*) PTYPE            ! Indicates the type desired records
                                       !   ( AREA, POINT, or MOBILE )

C...........    EXTERNAL functions and subroutines:

        INTEGER       STR2INT
        EXTERNAL      STR2INT

C...........   LOCAL VARIABLES their descriptions:

        INTEGER       BADVAL       ! Return value for LINE not okay
        INTEGER       OKVAL        ! Return value for LINE okay
        PARAMETER   ( OKVAL  = 0 )
        PARAMETER   ( BADVAL = 1 )

        INTEGER       IMON         ! Monthly profile code
        INTEGER       IWEK         ! Weekly  profile code
        INTEGER       IDIU         ! Diurnal profile code

        CHARACTER*256 MESG         ! Message buffer

C***********************************************************************
C   begin body of subroutine CHKXLIN

        CHKXLIN = OKVAL

C.........  Check ASCII input for errors
        IF ( IOS .NE. 0 ) THEN
            WRITE( MESG,94010 )
     &            'I/O error', IOS,
     &            'reading TEMPORAL XREF file at line', LINENO
            CALL M3MESG( MESG )
            CHKXLIN = BADVAL
            RETURN
        END IF

C.........  Area xref records have non-char SCC and blank plant ID
        IF( PTYPE .EQ. 'AREA' ) THEN
            IF( LINE( 21:21) .LE. '9' .AND. LINE( 7:11 ) .LE. ' ') THEN
                GO TO 111
            ELSE
                CHKXLIN = BADVAL
                RETURN
            ENDIF

C.........  Point xref records have non-char SCC and non-blank plant ID
        ELSEIF( PTYPE .EQ. 'POINT' ) THEN
            IF( LINE( 21:21) .LE. '9' .AND. LINE( 7:11 ) .GT. ' ') THEN
                GO TO 111
            ELSE
                CHKXLIN = BADVAL
                RETURN
            ENDIF

C.........  Mobile xref records have char SCC
        ELSEIF( PTYPE .EQ. 'MOBILE' ) THEN
            IF( LINE( 21:21) .GT. '9' ) THEN
                GO TO 111
            ELSE
                CHKXLIN = BADVAL
                RETURN
            ENDIF

        ELSE
            CALL M3EXIT( CHKXLIN, 0, 0, 
     &           'XREF processing not available for type ' // PTYPE, 2 )

        ENDIF

 111    CONTINUE

C.........  Check x-ref codes input for bad values
        IMON = STR2INT( LINE( 38:40 ) )
        IWEK = STR2INT( LINE( 42:44 ) )
        IDIU = STR2INT( LINE( 45:48 ) )
        
        IF ( IMON .LE. 0 ) THEN       ! Bad monthly profile code

            CHKXLIN = BADVAL
            WRITE( MESG,94010 )
     &             'Bad IMON=', IMON,
     &             'at line', LINENO,
     &             'for FIP ' // LINE( 1:5 ) //
     &             'and SCC ' // LINE( 21:30 ) 
            CALL M3MESG( MESG )

        ENDIF

        IF ( IWEK .LE. 0 ) THEN  ! Bad weekly profile code

            CHKXLIN = BADVAL
            WRITE( MESG,94010 )
     &             'Bad IWEK=', IWEK,
     &             'at line', LINENO,
     &             'for FIP ' // LINE( 1:5 ) //
     &             'and SCC ' // LINE( 21:30 ) 
            CALL M3MESG( MESG )
 
        ENDIF

        IF ( IDIU .LE. 0 ) THEN  ! Bad diurnal profile code
 
            CHKXLIN = BADVAL
            WRITE( MESG,94010 )
     &            'Bad IDIU=', IDIU,
     &            'at line', LINENO,
     &            'for FIP ' // LINE( 1:5 ) //
     &            'and SCC ' // LINE( 21:30 )
            CALL M3MESG( MESG )

        ENDIF

        RETURN

C...........   Internal buffering formats.............94xxx
 
94010   FORMAT( 10( A, :, I10, :, 2X ) )

        END
