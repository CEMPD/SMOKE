
        SUBROUTINE GETM3EPI( TZONE, SDATE, STIME, NSTEPS )

C***********************************************************************
C  subroutine body starts at line
C
C  DESCRIPTION: 
C     This subroutine retrieves the models-3 episode environment
C     variables, using the defaults provided as input.  The time step for
C     emissions must be one hour, so this ends up being one no matter what.
C
C  PRECONDITIONS REQUIRED:
C     
C
C  SUBROUTINES AND FUNCTIONS CALLED:  M3EXIT
C
C  REVISION  HISTORY:
C       Created 3/99 by M Houyoux
C
C***********************************************************************
C
C Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
C                System
C File: @(#)$Id$
C
C COPYRIGHT (C) 1998, MCNC--North Carolina Supercomputing Center
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
        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters

C...........   ARGUMENTS and their descriptions:
        INTEGER, INTENT(IN    ) :: TZONE     ! Time zone or < 0 if unknown
        INTEGER, INTENT(IN OUT) :: SDATE     ! Start date (YYYYDDD)
        INTEGER, INTENT(IN OUT) :: STIME     ! Start time (HHMMSS)
        INTEGER, INTENT(IN OUT) :: NSTEPS    ! Number of time steps
     
C...........   EXTERNAL FUNCTIONS:
        CHARACTER*2     CRLF
        INTEGER         ENVINT  
        INTEGER         GETDATE  
        INTEGER         GETNUM  
        CHARACTER*14    MMDDYY
        INTEGER         SECSDIFF
        INTEGER         STR2INT

        EXTERNAL   CRLF, ENVYN, GETDATE, GETNUM, MMDDYY, SECSDIFF,
     &             STR2INT

C...........   LOCAL VARIABLES their descriptions:

        INTEGER       EDATE        ! ending date from arguments 
        INTEGER       ETIME        ! ending time from arguments 
        INTEGER       G_EDATE      ! episode end date (YYYYDDD)
        INTEGER       G_ETIME      ! episode end time (HHMMSS)
        INTEGER       G_NSTEPS     ! episode number of time steps
        INTEGER       G_SDATE      ! episode start date from E.V. (YYYYDDD)
        INTEGER       G_STIME      ! episode start time from E.V. (HHMMSS)
        INTEGER       G_TSTEP      ! episode time step from E.V. (HHMMSS)
        INTEGER       ISECS        ! tmp duration in seconds
        INTEGER       IOS          ! tmp i/o status
        INTEGER       JRUNLEN      ! episode duration from E.V. (HHMMSS)
        INTEGER       N            ! tmp duration HHMMSS 
        INTEGER       TSTEP        ! time step (HHMMSS) 

        CHARACTER*300 MESG         ! Message buffer

        CHARACTER*14    DTBUF        ! date buffer
        CHARACTER*16 :: PROGNAME = 'GETM3EPI'    ! Program name

C***********************************************************************
C   begin body of subroutine GETM3EPI

C.........  Hardcoded actual time step
        TSTEP = 10000

C.........  Compute ending time from subroutine arguments, using 1-hour time
C           step assumption
        EDATE = SDATE
        ETIME = STIME
        CALL NEXTIME( EDATE, ETIME, NSTEPS * 10000 )

C.........  Get episode settings from the environment (using defaults from
C           actual files in case the environment variables are not set)
        G_SDATE = ENVINT( 'G_STDATE', 'Start date (YYYYDDD)', 
     &                    SDATE, IOS )
        G_STIME = ENVINT( 'G_STTIME', 'Start time (HHMMSS)', 
     &                    STIME, IOS )
        G_TSTEP = ENVINT( 'G_TSTEP', 'Time step (HHMMSS)', 
     &                    TSTEP, IOS )

        N = NSTEPS * TSTEP
        JRUNLEN= ENVINT( 'G_RUNLEN', 'Duration (HHMMSS)' , N, IOS )
        G_NSTEPS = JRUNLEN / G_TSTEP

C.........  Compare environment-based episode settings to those of the files
C.........  If the environment settings are more restrictive, then reset.
        G_EDATE = G_SDATE
        G_ETIME = G_STIME
        CALL NEXTIME( G_EDATE, G_ETIME, G_NSTEPS * G_TSTEP )

        ISECS = SECSDIFF( SDATE, STIME, G_SDATE, G_STIME )
        IF( ISECS .GT. 0 ) THEN  ! environment settings are later
            SDATE = G_SDATE
            STIME = G_STIME
        ELSEIF( ISECS .LT. 0 ) THEN
            MESG = 'Episode start of episode must be later ' //
     &             'than is set by the environment ' //
     &             CRLF() // BLANK5 // 'because of the input files.'
            CALL M3WARN( PROGNAME, 0, 0, MESG )

        ENDIF

        ISECS = SECSDIFF( EDATE, ETIME, G_EDATE, G_ETIME )
        IF( ISECS .LT. 0 ) THEN  ! environment settings are earlier
            EDATE = G_EDATE
            ETIME = G_ETIME

        ELSEIF( ISECS .GT. 0 ) THEN
            MESG = 'End of episode must be earlier than is set ' //
     &             'by the environment ' //
     &             CRLF() // BLANK5 // 'because of the input files.'
            CALL M3WARN( PROGNAME, 0, 0, MESG )

        ENDIF

        NSTEPS = SECSDIFF( SDATE, STIME, EDATE, ETIME ) / 3600 ! assume 1 hour

        SDATE = GETDATE( SDATE,
     &         'Enter simulation starting date (YYYYDDD)|(YYYYMMDD)' )

        STIME = GETNUM( 0, 235959, STIME, 
     &                  'Enter simulation starting time (HHMMSS)' )
  
        NSTEPS= GETNUM( 1, 999999, NSTEPS,
     &                  'Enter output duration (hours)' )

C.........  Write message about episode.  Time zone might be unknown
        IF( TZONE .GE. 0 ) THEN

            DTBUF = MMDDYY( SDATE )
            WRITE( MESG,94050 )
     &        'Output Time Zone:', TZONE,           CRLF() // BLANK5 //
     &        '      Start Date:', DTBUF( 1:LEN_TRIM( DTBUF ) ) //
     &                                              CRLF() // BLANK5 //
     &        '      Start Time:', STIME,'HHMMSS'// CRLF() // BLANK5 //
     &        '       Time Step:', 1    ,'hour'  // CRLF() // BLANK5 //
     &        '        Duration:', NSTEPS, 'hours'
 
        ELSE

            DTBUF = MMDDYY( SDATE )
            WRITE( MESG,94052 )
     &        '      Start Date:', DTBUF( 1:LEN_TRIM( DTBUF ) ) //
     &                                              CRLF() // BLANK5 //
     &        '      Start Time:', STIME,'HHMMSS'// CRLF() // BLANK5 //
     &        '       Time Step:', 1    ,'hour'  // CRLF() // BLANK5 //
     &        '        Duration:', NSTEPS, 'hours'
 
        END IF

        CALL M3MSG2( MESG )

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats.............94xxx
 
94010   FORMAT( 10( A, :, I6, :, 2X ) )

94050   FORMAT( A, 1X, I2.2, A, 1X, A, 1X, I6.6, 1X,
     &          A, 1X, I3.3, 1X, A, 1X, I3.3, 1X, A   )

94052   FORMAT( A, 1X, A, 1X, I6.6, 1X,
     &          A, 1X, I3.3, 1X, A, 1X, I3.3, 1X, A   )


        END SUBROUTINE GETM3EPI
