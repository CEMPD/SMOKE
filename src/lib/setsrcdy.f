
        SUBROUTINE SETSRCDY( NSRC, SDATE, TZONES, LDAYSAV, MFLAG,
     &                       DAYBEGT, DAYENDT )

C***********************************************************************
C  subroutine SETSRCDY body starts at line < >
C
C  DESCRIPTION:
C
C  PRECONDITIONS REQUIRED:
C     Sets the start and end of the day for all sources. 
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION HISTORY:
C
C***************************************************************************
C 
C Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
C                System
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

C.........  MODULES for public variables
C.........  This module contains the information about the source category
        USE MODINFO, ONLY: CATEGORY
        
        IMPLICIT NONE

C...........   EXTERNAL FUNCTIONS 
        INTEGER    ENVINT
        LOGICAL    ISDSTIME

        EXTERNAL   ENVINT, ISDSTIME

C...........   SUBROUTINE ARGUMENTS
        INTEGER,      INTENT (IN) :: NSRC             ! no. sources
        INTEGER,      INTENT (IN) :: SDATE            ! start julian date in GMT
        INTEGER,      INTENT (IN) :: TZONES  ( NSRC ) ! time zones per source
        LOGICAL,      INTENT (IN) :: LDAYSAV ( NSRC ) ! true: use daylight time
        LOGICAL,      INTENT (IN) :: MFLAG            ! true: setting for MOVES. false: for MOBILE6
        INTEGER,      INTENT(OUT) :: DAYBEGT( NSRC )  ! start time of SDATE
        INTEGER,      INTENT(OUT) :: DAYENDT( NSRC )  ! end time of SDATE

C...........   Other local variables

        INTEGER       EDATE        ! end julian date of day for source
        INTEGER       ETIME        ! end time HHMMSS
        INTEGER       IOS          ! status of ENVINT
        INTEGER       JDATE        ! tmp Julian date
        INTEGER       S            ! source no. 
        INTEGER       STIME        ! start time HHMMSS 
        INTEGER       STIME_SET    ! start time HHMMSS when all srcs get same

        LOGICAL       :: DAYLIT = .FALSE. ! true: date is daylight savings
        LOGICAL, SAVE :: FIRSTIME = .TRUE.! true: first time routine called
        LOGICAL, SAVE :: UFLAG  = .FALSE. ! true: all srcs use same day start

        CHARACTER(300) MESG          ! message buffer

        CHARACTER(16) :: PROGNAME = 'SETSRCDY' ! program name

C***********************************************************************
C   begin body of subroutine SETSRCDY

C.........  Get environment variable for old-style processing in which
C           all of the sources use the same start and end of the day
        IF( FIRSTIME ) THEN

            MESG = 'Start time for using uniform start time'
            STIME_SET = ENVINT ( 'UNIFORM_STIME', MESG, -1, IOS )

            UFLAG = ( STIME_SET > 0 )

            IF( UFLAG ) THEN

                WRITE( MESG,94010 ) 'NOTE: A daily start time of ',
     &                 STIME_SET, 'is being used for all sources.'
                CALL M3MSG2( MESG )

            END IF

            FIRSTIME = .FALSE.

        END IF

C.........  Set start and end time for all sources and return when uniform
C           start time is being used
        IF( UFLAG ) THEN

            EDATE = SDATE
            ETIME = STIME_SET
            CALL NEXTIME( EDATE, ETIME, 230000 )

            DAYBEGT = STIME_SET   ! array
            DAYENDT = ETIME       ! array

            RETURN

        END IF

C.........  Determine if this date is in the range for daylight savings
         
        DAYLIT = ISDSTIME( SDATE )

C.........  Loop through sources and set start and end date in GMT for
C           each source.  Note that the time zone, adjusted for daylight
C           savings, is the same as the start hour of the day in GMT.0
        DO S = 1, NSRC

            STIME = TZONES( S ) * 10000
            JDATE = SDATE

C.............  If this date is during daylight savings, and if this
C               source is affected by daylight savings
            IF( DAYLIT .AND. LDAYSAV( S ) ) THEN

                CALL NEXTIME( JDATE, STIME, -10000 )
                
            END IF

C.............  Set start time to  6 A.M. local time for MOBILE6 processing
            IF( CATEGORY == 'MOBILE' .AND. .NOT. MFLAG ) THEN
                CALL NEXTIME( JDATE, STIME, 60000 )
            END IF

C.............  Store start time
            DAYBEGT( S ) = STIME

C.............  Compute end date and time
            EDATE = JDATE
            ETIME = STIME
            CALL NEXTIME( EDATE, ETIME, 230000 )

C.............  Store end time
            DAYENDT( S ) = ETIME

        END DO

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I9, :, 1X ) )
 
        END SUBROUTINE SETSRCDY
