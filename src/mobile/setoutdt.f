
        SUBROUTINE SETOUTDT( NSRC, SDATE, TZONES, DAYBEGT, DAYENDT,
     &                       OUTDATE, OUTTIME )

C***********************************************************************
C  subroutine SETOUTDT body starts at line < >
C
C  DESCRIPTION:
C
C  PRECONDITIONS REQUIRED:
C     Sets the output date and time based on the time zones and ending hours 
C     of the sources 
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

C...........   EXTERNAL FUNCTIONS 
        INTEGER    SECSDIFF

        EXTERNAL   SECSDIFF

C...........   SUBROUTINE ARGUMENTS
        INTEGER, INTENT (IN) :: NSRC             ! no. sources
        INTEGER, INTENT (IN) :: SDATE            ! comparison julian date in GMT
        INTEGER, INTENT (IN) :: TZONES ( NSRC )  ! time zones per source
        INTEGER, INTENT (IN) :: DAYBEGT( NSRC )  ! start time of SDATE for src
        INTEGER, INTENT (IN) :: DAYENDT( NSRC )  ! end time of SDATE for src
        INTEGER, INTENT(OUT) :: OUTDATE          ! julian date to output min/max
        INTEGER, INTENT(OUT) :: OUTTIME          ! time to output min/max

C...........   Other local variables

        INTEGER       EDATE        ! end julian date of day for source
        INTEGER       ETIME        ! end time HHMMSS
        INTEGER       DSTEP        ! a differenced time step in HHMMSS
        INTEGER       ISDIFF       ! diff in end time and julian date in secs
        INTEGER       JDATE        ! tmp julian time
        INTEGER       JTIME        ! tmp time HHMMSS
        INTEGER       LSDIFF       ! previous ISDIFF of loop
        INTEGER       S            ! source no. 

        CHARACTER*300 MESG          ! message buffer

        CHARACTER*16 :: PROGNAME = 'SETOUTDT' ! program name

C***********************************************************************
C   begin body of subroutine SETOUTDT

C.........  Loop through sources, determine the ending date, and using the
C           ending date and time, compute the difference between the ending
C           dates and times and the reference ones.
        LSDIFF = 0
        DO S = 1, NSRC

C.............  Handle special case of daylight time in time zone zero
            IF( TZONES( S )  .EQ. 0      .AND. 
     &          DAYBEGT( S ) .EQ. 230000       ) THEN
                JDATE = SDATE
                CALL NEXTIME( JDATE, 0, -10000 )
            ELSE
                JDATE = SDATE
            END IF

C.............  Re-calculate the end date in GMT (note that daybegt and
C               dayendt are in GMT)
            EDATE = JDATE
            JTIME = DAYBEGT( S )
            CALL NEXTIME( EDATE, JTIME, 230000 )
            
C.............  Compute the difference in second between the julian date and
C               time in GMT and the sources ending date and time

            ETIME = DAYENDT( S )
            ISDIFF = SECSDIFF( SDATE, 0, EDATE, ETIME )
            IF( ISDIFF .GT. LSDIFF ) THEN
                OUTDATE = EDATE
                OUTTIME = ETIME
                LSDIFF  = ISDIFF
            END IF

        END DO

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I9, :, 1X ) )
 
        END SUBROUTINE SETOUTDT
