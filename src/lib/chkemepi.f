
        LOGICAL FUNCTION CHKEMEPI( JDATE, JTIME, TIMESTEP, JRUNLEN )

C***********************************************************************
C  function body starts at line
C
C  DESCRIPTION:
C      This subroutine checks the episode parameters to make sure they 
C      are consistent with what SMOKE expects.
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C
C**************************************************************************
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
C***************************************************************************

        IMPLICIT NONE

C...........   INCLUDES:

        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters

C...........   EXTERNAL FUNCTIONS:
        CHARACTER*2     CRLF
        EXTERNAL        CRLF

C...........   SUBROUTINE ARGUMENTS
        INTEGER         JDATE       ! Julian date (DDDYYYY) (in)
        INTEGER         JTIME       ! time (HHMMSS) (in)
        INTEGER         TIMESTEP    ! duration of time step (HHMMSS) (in)
        INTEGER         JRUNLEN     ! episode duration (HHMMSS) (in)

C...........   Other local variables
        REAL            REMNDR  ! REMNDRder for checking time units

        CHARACTER*300   MESG   ! message buffer 

        CHARACTER*16 :: PROGNAME = 'CHKEMEPI' ! program name

C***********************************************************************
C   begin body of function CHKEMEPI

        CHKEMEPI = .TRUE.

        IF( JDATE .LT. 1970001 .OR. JDATE .GT. 2100365 ) THEN

            CHKEMEPI = .FALSE.
            MESG = 'Date is out of acceptable modeling range:' //
     &             CRLF() // BLANK16 // '1970001 to 2100365.'
            CALL M3WARN( PROGNAME, 0, 0, MESG )
            
        ENDIF
        
        IF( JTIME .LT. 0 .OR. JTIME .GT. 235959 ) THEN

            CHKEMEPI = .FALSE.
            MESG = 'Time is out of acceptable modeling range:' //
     &             CRLF() // BLANK16 // '0 to 235959.'
            CALL M3WARN( PROGNAME, 0, 0, MESG )
            
        ENDIF

        REMNDR = MOD( REAL( JTIME ), 10000.0 )

        IF( REMNDR .GT. 0 ) THEN

            CHKEMEPI = .FALSE.
            MESG = 'Time may only be in units of hours.'
            CALL M3WARN( PROGNAME, 0, 0, MESG )
            
        ENDIF

        IF( TIMESTEP .NE. 10000 ) THEN

            CHKEMEPI = .FALSE.
            WRITE( MESG,94010 )
     &             'Timestep may only be one hour, but it is', TIMESTEP,
     &             '(HHMMSS).'
            CALL M3WARN( PROGNAME, 0, 0, MESG )
            
        ENDIF

        REMNDR = MOD( REAL( JRUNLEN ), 10000.0 )

        IF( REMNDR .GT. 0 ) THEN

            CHKEMEPI = .FALSE.
            MESG = 'Run duration may only be in units of hours.'
            CALL M3WARN( PROGNAME, 0, 0, MESG )
            
        ENDIF

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Formatted file I/O formats............ 93xxx

93000   FORMAT( A )

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

        END FUNCTION CHKEMEPI

