
        SUBROUTINE CHKFULLDY( NSRC, SDATE, STIME, EDATE, ETIME, 
     &                        TZONES, LDAYSAV, MODELNAM )
   
C***********************************************************************
C  subroutine CHKFULLDY body starts at line < >
C
C  DESCRIPTION:
C      This routine preprocesses the dates, times, and time zones and looks
C      for which sources do not have complete days at the beginning and
C      end of the period of interest.  For these sources, it tallies information
C      about them, and reports it to the log file.  This routine is needed
C      when min/max information is computed per day, but entire days might not
C      be present.
C
C  PRECONDITIONS REQUIRED:
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

C...........   INCLUDES

        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters

C...........   EXTERNAL FUNCTIONS 
        CHARACTER*2  CRLF

        EXTERNAL     CRLF

C...........   SUBROUTINE ARGUMENTS
        INTEGER, INTENT (IN) :: NSRC            ! no. sources
        INTEGER, INTENT (IN) :: SDATE           ! start julian date
        INTEGER, INTENT (IN) :: STIME           ! start time HHMMSS
        INTEGER, INTENT (IN) :: EDATE           ! end julian date
        INTEGER, INTENT (IN) :: ETIME           ! end time HHMMSS
        INTEGER, INTENT (IN) :: TZONES ( NSRC ) ! time zones per source
        LOGICAL, INTENT (IN) :: LDAYSAV( NSRC ) ! true: use daylight time
        CHARACTER(*), INTENT (IN) :: MODELNAM   ! emission factor model name

C...........  Variables dimensioned by subroutine arguments
        INTEGER       SRTDAYHR( NSRC )
        INTEGER       ENDDAYHR( NSRC )

C...........  Arrays to report by time zone
        INTEGER       NMISSBEG( 0:23 )  ! no. srcs missing hours on first day
        INTEGER       NMISSEND( 0:23 )  ! no. srcs missing hours on last day
        REAL          AVHRBEG ( 0:23 )  ! ave. no. missing hours
        REAL          AVHREND ( 0:23 )  ! ave. no. missing hours
        LOGICAL       ZONEFLAG( 0:23 )  ! true: time zone is in inventory

C...........   Other local variables
        INTEGER       S, Z         ! counters and indices
        INTEGER       HRDIFF       ! tmp hr difference
        INTEGER       ODATE        ! output date (needed only for subrtn call)
        INTEGER       OTIME        ! output time (needed only for subrtn call)

        CHARACTER*300 MESG          ! message buffer

        CHARACTER*16 :: PROGNAME = 'CHKFULLDY' ! program name

C***********************************************************************
C   begin body of subroutine CHKFULLDY

C.........  Initialize information to report missing hours for the first and
C           last day of processing
        NMISSBEG = 0    ! array
        NMISSEND = 0    ! array
        AVHRBEG  = 0.   ! array
        AVHREND  = 0.   ! array

C.........  Determine which time zones are in the inventory

        ZONEFLAG = .FALSE.   ! array

        DO S = 1, NSRC
            Z = TZONES( S )
            ZONEFLAG( Z ) = .TRUE.
        END DO

C.........  Create arrays that contain the start hour and end hour of a day 
C           (in GMT) for each source (it changes by day because of daylight
C           savings) for the first day.
        CALL SETSRCDY( NSRC, SDATE, TZONES, LDAYSAV, 
     &                 SRTDAYHR, ENDDAYHR, MODELNAM )
     
C.........  Count how many sources are missing part of the first day of
C           meteorology data for each time zone that has sources.  This is done
C           assuming that the start time given is in GMT.
        DO S = 1, NSRC

            Z = TZONES( S )
            
            HRDIFF = ( STIME - SRTDAYHR( S ) ) / 10000
            IF( HRDIFF > 0 ) THEN
                NMISSBEG( Z ) = NMISSBEG( Z ) + 1
                AVHRBEG ( Z ) = AVHRBEG ( Z ) + REAL( HRDIFF ) ! tmp sum
            END IF    

        END DO

C.........  Get the day start and end hour for the last day of processing
        CALL SETSRCDY( NSRC, EDATE, TZONES, LDAYSAV, 
     &                 SRTDAYHR, ENDDAYHR, MODELNAM )

C.........  Count how many sources are missing part of the last day of
C           meteorology data for each time zone that has sources.  The end time
C           is assumed to be in GMT.
        DO S = 1, NSRC

            Z = TZONES( S )
            
            HRDIFF = ( ETIME - ENDDAYHR( S ) ) / 10000
            IF( HRDIFF < 0 ) THEN
                NMISSEND( Z ) = NMISSEND( Z ) + 1
                AVHREND ( Z ) = AVHREND ( Z ) - REAL( HRDIFF ) ! tmp sum
            END IF    

        END DO

C.........  Compute average number of days per time zone.  This is averaged
C           in case some sources in a given time zone use daylight savings
C           and some do not.
        DO Z = 0, 23

             IF( NMISSBEG( Z ) > 0 ) THEN
                 AVHRBEG( Z ) = AVHRBEG( Z ) / REAL( NMISSBEG( Z ) )
             END IF

             IF( NMISSEND( Z ) > 0 ) THEN
                 AVHREND( Z ) = AVHREND( Z ) / REAL( NMISSEND( Z ) )
             END IF

        END DO

C.........  Write to log file...

        MESG = 'NOTE: The following statistics indicate how well '//
     &         'minimum and maximum ' // CRLF() // BLANK10 //
     &         'values can be computed using the start and end ' //
     &         'times provided with ' // CRLF() // BLANK10 //
     &         'input meteorology data.'
        CALL M3MESG( MESG )

        MESG = BLANK5 // 'These values should be reviewed to ensure '//
     &         'satisfactory minimum' // CRLF() // BLANK10 //
     &         'and maximum values computation.'

        CALL M3MESG( MESG )

        MESG = '         No. sources affected   Ave. missing hours' //
     &         CRLF() // BLANK5 //
     &         ' Zone        Start      End        Start      End'
        CALL M3MESG( MESG )

        DO Z = 0, 23

            IF( ZONEFLAG( Z ) ) THEN

                WRITE( MESG, 94160 ) Z, NMISSBEG( Z ), NMISSEND( Z ), 
     &                                  AVHRBEG ( Z ), AVHREND ( Z )
                CALL M3MESG( MESG )

            END IF

        END DO

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I9, :, 1X ) )

94160   FORMAT( 2X, I2.2, 6X, I8, 1X, I8, 5X, F8.2, 1X, F8.2 )

        END SUBROUTINE CHKFULLDY
