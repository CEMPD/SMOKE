
        SUBROUTINE GETM3EPI( TZONE, SDATE, STIME, TSTEP, NSTEPS )

C***********************************************************************
C  subroutine body starts at line
C
C  DESCRIPTION: 
C     This subroutine retrieves the models-3 episode environment
C     variables, using the defaults provided as input.  The time step for
C     emissions must be one hour, so this ends up being one no matter what.
C     If the TZONE or NSTEPS arguments are negative, the subroutine will
C     ignore the steps that involve these arguments.
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
C****************************************************************************

        IMPLICIT NONE

C...........   INCLUDES:
        INCLUDE 'IOCNST3.EXT'   !  emissions constant parameters

C...........   ARGUMENTS and their descriptions:
        INTEGER, INTENT(IN    ) :: TZONE     ! Time zone or < 0 if unknown
        INTEGER, INTENT(IN OUT) :: SDATE     ! Start date (YYYYDDD)
        INTEGER, INTENT(IN OUT) :: STIME     ! Start time (HHMMSS)
        INTEGER, INTENT(IN OUT) :: TSTEP     ! Time step (HHMMSS)
        INTEGER, INTENT(IN OUT) :: NSTEPS    ! Number of time steps
     
C...........   EXTERNAL FUNCTIONS:
        CHARACTER(2)    CRLF
        INTEGER         ENVINT  
        INTEGER         GETDATE  
        INTEGER         GETNUM  
        CHARACTER(14)   MMDDYY
        INTEGER         SECSDIFF
        INTEGER         STR2INT

        EXTERNAL   CRLF, ENVINT, GETDATE, GETNUM, MMDDYY, SECSDIFF,
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
        INTEGER       L            ! tmp string length
        INTEGER       N            ! tmp duration HHMMSS 

        CHARACTER(300) MESG        ! Message buffer
        CHARACTER(14)  DTBUF       ! date buffer
        
        CHARACTER(16) :: PROGNAME = 'GETM3EPI'    ! Program name

C***********************************************************************
C   begin body of subroutine GETM3EPI

C.........  Hardcoded actual time step
        TSTEP = 10000

C.........  Get episode settings from the environment (using defaults from
C           actual files in case the environment variables are not set)
        G_SDATE = ENVINT( 'G_STDATE', 'Start date (YYYYDDD)', 
     &                    SDATE, IOS )
        G_STIME = ENVINT( 'G_STTIME', 'Start time (HHMMSS)', 
     &                    STIME, IOS )
        G_TSTEP = ENVINT( 'G_TSTEP', 'Time step (HHMMSS)', 
     &                    TSTEP, IOS )
        TSTEP = G_TSTEP

        IF( NSTEPS .GT. 0 ) THEN  ! only reset if argument is positive.
            N = NSTEPS * TSTEP
            JRUNLEN = ENVINT( 'G_RUNLEN', 'Duration (HHMMSS)' , N, IOS )
            IF( G_TSTEP .NE. 0 ) THEN
                G_NSTEPS = JRUNLEN / G_TSTEP
            ELSE
                G_NSTEPS = 0
            END IF
        END IF

C.........  Check if date and time are already set (by a file) and if so, make
C           comparisons to the environment values.
        IF( SDATE .GT. 0 ) THEN

C.............  Compute ending time from subroutine arguments, using 1-hour time
C               step assumption
            EDATE = SDATE
            ETIME = STIME
            CALL NEXTIME( EDATE, ETIME, NSTEPS * 10000 )

C.............  Compare environment-based episode settings to those of the files
C.............  If the environment settings are more restrictive, then reset.
            G_EDATE = G_SDATE
            G_ETIME = G_STIME
            CALL NEXTIME( G_EDATE, G_ETIME, G_NSTEPS * G_TSTEP )
  
            ISECS = SECSDIFF( SDATE, STIME, G_SDATE, G_STIME )
            IF( ISECS .GT. 0 ) THEN  ! environment settings are later
                SDATE = G_SDATE
                STIME = G_STIME
            ELSE IF( ISECS .LT. 0 ) THEN
                MESG = 'Start of episode must be later ' //
     &                 'than is set by the environment ' //
     &                 CRLF() // BLANK5 // 'because of the input files.'
                CALL M3WARN( PROGNAME, 0, 0, MESG )

            END IF

            ISECS = SECSDIFF( EDATE, ETIME, G_EDATE, G_ETIME )
            IF( ISECS .LT. 0 ) THEN  ! environment settings are earlier
                EDATE = G_EDATE
                ETIME = G_ETIME

            ELSE IF( ISECS .GT. 0 ) THEN
                MESG = 'End of episode must be earlier than is set ' //
     &                 'by the environment ' //
     &                 CRLF() // BLANK5 // 'because of the input files.'
                CALL M3WARN( PROGNAME, 0, 0, MESG )

            END IF

C.........  If date and time are not already set, use the environment variable
C           values as defaults
        ELSE

            SDATE  = G_SDATE
            STIME  = G_STIME            
            EDATE  = SDATE
            ETIME  = STIME

            IF( NSTEPS .GT. 0 ) THEN
                NSTEPS = G_NSTEPS
                CALL NEXTIME( EDATE, ETIME, G_NSTEPS*G_TSTEP )
            END IF

        END IF

        SDATE = GETDATE( SDATE,
     &         'Enter simulation starting date (YYYYDDD)|(YYYYMMDD)' )

        STIME = GETNUM( 0, 235959, STIME, 
     &                  'Enter simulation starting time (HHMMSS)' )
  
C.........  Assuming NSTEPS is requested by the calling program...
C.........  Recompute number of steps assuming 1 hour time steps
        IF( NSTEPS .GT. 0 ) THEN
            N = SECSDIFF( SDATE, STIME, EDATE, ETIME ) / 3600 

C.............  Double check the number of time steps to ensure that the input 
C               files are not out of range from the environment settings
            IF( N .LE. 0 ) THEN
                MESG= 'Date/time from input file(s) are inconsistent '//
     &                'with the environment ' // CRLF() // BLANK10 //
     &                'episode settings.' 
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF

            NSTEPS= GETNUM( 1, G_NSTEPS, N,
     &                      'Enter output duration (hours)' )

            IF ( NSTEPS > N ) THEN
                 MESG= 'Number of time steps selected is most '//
     &                 'likely greater than those available ' //
     &                 'from input data.'
                 CALL M3WARN( PROGNAME, 0, 0, MESG ) 
            ENDIF

        END IF

C.........  Get date buffer field
        DTBUF = MMDDYY( SDATE )
        L = LEN_TRIM( DTBUF )

C.........  Write message about episode.  Time zone might be unknown
        IF( TZONE .GE. 0 .AND. NSTEPS .GT. 0 ) THEN

            WRITE( MESG,94050 )
     &        'Output Time Zone:', TZONE,           CRLF() // BLANK5 //
     &        '      Start Date:', DTBUF( 1:L ) //  CRLF() // BLANK5 //
     &        '      Start Time:', STIME,'HHMMSS'// CRLF() // BLANK5 //
     &        '       Time Step:', 1    ,'hour'  // CRLF() // BLANK5 //
     &        '        Duration:', NSTEPS, 'hours'
 
        ELSE IF( NSTEPS .GT. 0 ) THEN

            WRITE( MESG,94052 )
     &        '      Start Date:', DTBUF( 1:L ) //  CRLF() // BLANK5 //
     &        '      Start Time:', STIME,'HHMMSS'// CRLF() // BLANK5 //
     &        '       Time Step:', 1    ,'hour'  // CRLF() // BLANK5 //
     &        '        Duration:', NSTEPS, 'hours'

        ELSE IF( TZONE .GE. 0 ) THEN

            WRITE( MESG,94050 )
     &        'Output Time Zone:', TZONE,           CRLF() // BLANK5 //
     &        '      Start Date:', DTBUF( 1:L ) //  CRLF() // BLANK5 //
     &        '      Start Time:', STIME,'HHMMSS'

        ELSE

            WRITE( MESG,94052 )
     &        '      Start Date:', DTBUF( 1:L ) //  CRLF() // BLANK5 //
     &        '      Start Time:', STIME,'HHMMSS'
 
        END IF

        CALL M3MSG2( MESG )

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats.............94xxx
 
94010   FORMAT( 10( A, :, I6, :, 2X ) )

94050   FORMAT( A, 1X, I2.2, A, 1X, A, 1X, I6.6, 1X,
     &          A, 1X, I5.3, 1X, A, 1X, I5.3, 1X, A   )

94052   FORMAT( A, 1X, A, 1X, I6.6, 1X,
     &          A, 1X, I5.3, 1X, A, 1X, I5.3, 1X, A   )


        END SUBROUTINE GETM3EPI
