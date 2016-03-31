
C.........................................................................
C Version "@(#)$Header$"
C EDSS/Models-3 I/O API.  Copyright (C) 1992-2002 MCNC
C Distributed under the GNU LESSER GENERAL PUBLIC LICENSE version 2.1
C See file "LGPL.txt" for conditions of use.
C.........................................................................

      SUBROUTINE NEXTIME  ( JDATE , JTIME, DTIME )

C.....................................................................
C       subroutine body starts approx line 62
C
C  FUNCTION:
C
C       Subroutine to add  DTIME  to the current time, and then
C       update the clock  JTIME  and calendar  JDATE  accordingly.  
C       Output is fully normalized.
C
C   JDATE is stored in the form   YYYYDDD = YEAR*1000  +  DAY
C   JTIME is stored in the form   HHMMSS  = HOUR*10000 +  MINS*100  +  SECS
C   DTIME is stored in the form   HHMMSS  = HOUR*10000 +  MINS*100  +  SECS
C
C  REVISION HISTORY: 
C       prototype 10/1990 by Carlie J. Coats, Jr., MCNC Environmental Programs
C       Version    3/1993 by CJC for CRAY, etc.
C       Bugfix     1/2002 by CJC:  stepping _backwards_ across a year
C       transition
C       Bugfix     5/2002 by CJC:  stepping _backwards_ by 24 hours.
C
C  CALLS:
C      none
C
C........................................................................

        IMPLICIT  NONE

C.......   ARGUMENTS:

        INTEGER         JDATE           !  date (encoded YYYYDDD)
        INTEGER         JTIME           !  time (encoded  HHMMSS)
        INTEGER         DTIME           !  time increment (encoded HHMMSS)


C.......   LOCAL VARIABLES:  day and year components of date

        INTEGER         YEAR            !  year-component    of JDATE
        INTEGER         DAYS            !  day-component     of JDATE
        INTEGER         HOUR            !  hours-component   of JTIME
        INTEGER         MINS            !  minutes-component of JTIME
        INTEGER         SECS            !  seconds-component of JTIME
        INTEGER         SCR             !  scratch accumulator
        INTEGER         SIGN            !  sign of DTIME
        INTEGER         ATIME           !  absolute value of DTIME


C........................................................................
C       begin  NEXTIME
C.......   Increment seconds part of JTIME by DTIME (secs), and
C.......   re-normalize minute, hour, day part of JDATE:JTIME

        IF ( DTIME .GE. 0 ) THEN
            ATIME = DTIME
            SIGN  = 1
        ELSE
            ATIME = - DTIME
            SIGN  = - 1
        END IF

        SECS  =  MOD ( JTIME , 100 )  +  SIGN * MOD ( ATIME , 100 )
        IF ( SECS .GE. 0 ) THEN
            SCR = SECS / 60
        ELSE
            SCR = - ( 60 - SECS ) / 60
        END IF

        SECS = SECS  -  SCR * 60
        MINS = SCR  +  MOD ( JTIME / 100 , 100 )
     &              +  SIGN * MOD ( ATIME / 100 , 100 )

        IF ( MINS .GE. 0 ) THEN
            SCR = MINS / 60
        ELSE
            SCR = - ( 60 - MINS ) / 60
        END IF

        MINS  =  MINS  -  SCR * 60
        HOUR  =  JTIME / 10000  +  SIGN * ( ATIME / 10000 )  +  SCR

        DAYS  =  MOD ( JDATE , 1000 )

        IF  ( HOUR .LT. -23 )  THEN
            SCR  = -HOUR / 24
            HOUR =  HOUR  +  SCR * 24
            DAYS =  DAYS  -  SCR
        END IF

        IF  ( HOUR .LT. 0 )  THEN
            SCR   =  ( 24 - HOUR ) / 24
            HOUR  =  HOUR  +  SCR * 24
            DAYS  =  DAYS  -  SCR
        ELSE IF  ( HOUR .GT. 23 )  THEN
            SCR   =  HOUR / 24
            HOUR  =  HOUR  -  SCR * 24
            DAYS  =  DAYS  +  SCR
        END IF

        JTIME  =  10000 * HOUR  +  100 * MINS  +  SECS


C...........   Update JDATE:
C...........   Note that each year must be treated individually

        YEAR  =  JDATE / 1000

100     CONTINUE        !  loop normalizing negative day numbers

            IF ( DAYS .LE. 0 ) THEN

                IF (           ( MOD (YEAR,4)   .NE. 1 )        !  nonleap year
     &             .OR. (      ( MOD (YEAR,100) .EQ. 1 )
     &                   .AND. ( MOD (YEAR,400) .NE. 1 ) ) ) THEN

                    DAYS  =  365   +  DAYS
                    YEAR  =  YEAR  -  1

                ELSE            !  leap-year case

                    DAYS  =  366   +  DAYS
                    YEAR  =  YEAR  -  1

                END IF

                GO TO  100

            END IF


200     CONTINUE        !  loop normalizing day numbers > 365,366

            IF ( DAYS .GE. 366 ) THEN

                IF (           ( MOD (YEAR,4)   .NE. 0 )        !  nonleap year
     &              .OR. (     ( MOD (YEAR,100) .EQ. 0 )
     &                   .AND. ( MOD (YEAR,400) .NE. 0 ) ) ) THEN

                    DAYS  =  DAYS   - 365
                    YEAR  =  YEAR  +   1

                    GO TO  200

                ELSE IF ( DAYS .GE. 367 ) THEN           !  leap year case

                    DAYS  =  DAYS   - 366
                    YEAR  =  YEAR  +   1

                    GO TO  200

                END IF

            END IF      !  end DAYS > 365,366 date-normalization loop


        JDATE  =   1000 * YEAR  +  DAYS


        RETURN
        END
