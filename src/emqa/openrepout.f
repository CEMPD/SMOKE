
        SUBROUTINE OPENREPOUT( FILNAM, FDEV )

C***********************************************************************
C  subroutine body starts at line
C
C  DESCRIPTION:
C    The OPENREPOUT routine open output file FILNAM, which can be a physical
C    or logical file name.  Returns FDEV to calling program.
C
C  PRECONDITIONS REQUIRED:
C    FILNAM is defined as a physical or logical file name.
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C
C***********************************************************************
C  
C Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
C                System
C File: @(#)$Id$
C  
C COPYRIGHT (C) 2000, MCNC--North Carolina Supercomputing Center
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
C***********************************************************************

        IMPLICIT NONE

C...........   INCLUDES:
        
        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters

C...........   EXTERNAL FUNCTIONS and their descriptions:

        CHARACTER*2   CRLF
        INTEGER       GETEFILE
        INTEGER       JUNIT

        EXTERNAL    CRLF, GETEFILE, JUNIT

C...........   SUBROUTINE ARGUMENTS
        CHARACTER(*), INTENT (IN) :: FILNAM   ! physical or logical file name
        INTEGER     , INTENT(OUT) :: FDEV     ! unit number of file

C...........   Local variables
        INTEGER         L               ! counters and indices
        INTEGER         IOS             ! i/o status

        CHARACTER*16    VARBUF          ! logical file name buffer
        CHARACTER*300   MESG            ! message buffer
        CHARACTER*300   PNAME           ! physical file name on ENVSTR test

        CHARACTER*16 :: PROGNAME = 'OPENREPOUT' ! program name

C***********************************************************************
C   begin body of subroutine OPENREPOUT

        IOS = -1
C.........  If file name is less than 16 characters, check if file name is a
C           defined environment variable
        L = LEN_TRIM( FILNAM )
        IF( L .LE. 16 ) THEN
            MESG = 'Check of file name for logical status'
            VARBUF = FILNAM
            CALL ENVSTR( VARBUF, MESG, ' ', PNAME, IOS )
        END IF

C.........  If it is a defined environment variable, open as logical file name
        IF( IOS .EQ. 0 ) THEN

            FDEV = GETEFILE( FILNAM, .FALSE., .TRUE., PROGNAME )

C.........  If it is not a defined environment variable, open as physical file
C           name
        ELSE

            FDEV = JUNIT()
            OPEN( FDEV,ERR=1006,FILE=FILNAM,STATUS='UNKNOWN',RECL=2500 )

        END IF

        RETURN

C.........  Error opening raw input file
1006    WRITE( MESG,94010 ) ' Could not open file:' //
     &         CRLF() // BLANK5 // FILNAM( 1:LEN_TRIM( FILNAM ) )
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I10, :, 1X ) )

        END SUBROUTINE OPENREPOUT
        
 
