
C Version "@(#)$Id$ $Source$ $Date$ 

        SUBROUTINE  MKTMAT( NSRC,  JDATE, MON, DAY, TZONE, ZONES,
     &                      MFAC,  MDEX, 
     &                      WFAC,  WDEX, 
     &                      DFAC,  DDEX, 
     &                      NDAY,  INDD, TMAT )

C***********************************************************************
C  subroutine body starts at line  93
C
C  FUNCTION:
C       Construct temporal-coefficient-transform matrices for 
C       program TMPPOINT
C
C  PRECONDITIONS REQUIRED:
C       Temporal profile arrays for monthly, weekly, diurnal profiles.
C       MDEX, WDEX entries set to zero for month- or week-independent
C       source records.
C
C  SUBROUTINES AND FUNCTIONS CALLED:  none
C
C  REVISION  HISTORY:
C       Prototype 3/95 by CJC
C       Revised   8/95 by CJC:  daylight savings time correction
C       Revised  10/95 by CJC:  treatment of profile-independence,
C               substitution of day-specific data NDAY,INDD
C***********************************************************************

      IMPLICIT NONE

C...........   INCLUDES:

        INCLUDE 'TMDIMS3.EXT'   !  temporal-allocation dimensioning parameters
        INCLUDE 'CHDIMS3.EXT'   !  chemical-species    dimensioning parameters


C...........   ARGUMENTS and their descriptions:

        INTEGER         NSRC                    ! number of sources
        INTEGER         JDATE                   ! date YYYYDDD
        INTEGER         MON                     ! 1 ... 12
        INTEGER         DAY                     ! 1 ... 7
        INTEGER         TZONE                   ! time zone (5 for Eastern)
        INTEGER         ZONES( NSRC )           ! source time zones
        REAL            MFAC( 12, NMCOD )       ! Monthly profile factors
        REAL            WFAC(  7, NWCOD )       ! Weekly    "       "
        REAL            DFAC( 24, NDCOD  )      ! Diurnal   "       "
        INTEGER         MDEX( NSRC, NIPOL )     ! Monthly profile codes
        INTEGER         WDEX( NSRC, NIPOL )     ! weekly  profile codes
        INTEGER         DDEX( NSRC, NIPOL )     ! diurnal profile codes
        INTEGER         NDAY                    ! # of day-specific obs
        INTEGER         INDD( * )     ! subscripts: day-specific obs
        REAL            TMAT( NSRC, NIPOL, 24 ) ! temporal-profile coeffs

C...........   EXTERNAL FUNCTIONS:

        INTEGER         FIND1
        LOGICAL         ISDSTIME        !  true iff daylight savings time( date)

        EXTERNAL        FIND1, ISDSTIME

C...........   Scratch Local variables:

        INTEGER         H               !  hour-loop counter
        INTEGER         K               !  diurnal-profile index
        INTEGER         L               !  diurnal-xref index
        INTEGER         S               !  source counter
        INTEGER         V               !  chemistry variable counter
        INTEGER         HCORR           !  daylight savings time correction
        REAL            FAC             !  partial matrix factor
        REAL            YR2DAY          !  year-to-day conversion factor


C.......   SAVED Local Variables

        REAL            MWFAC( 12 )
        DATA            MWFAC / 31.0, 28.0, 31.0, 30.0,
     &                          31.0, 30.0, 31.0, 31.0,
     &                          30.0, 31.0, 30.0, 31.0 /

        LOGICAL         FIRSTIME
        DATA            FIRSTIME / .TRUE. /

        SAVE            MWFAC, FIRSTIME


C***********************************************************************
C   begin body of subroutine  PTMAT

C.......   Firstime:  renormalize MFAC, WFAC

        IF ( FIRSTIME ) THEN

            FIRSTIME = .FALSE.

            DO  13 S = 1, NMCOD

                FAC = 0.0
                DO  11 K = 1, 12
                    FAC = FAC + MWFAC( K ) * MFAC( K,S )
11              CONTINUE

                FAC = 1.0 / FAC

                DO  12 K = 1, 12
                    MFAC( K,S ) = FAC * MFAC( K,S )
12              CONTINUE

13          CONTINUE

            DO  23 S = 1, NWCOD
            DO  22 K = 1, 7
                WFAC( K,S ) = 7.0 * WFAC( K,S )
22          CONTINUE
23          CONTINUE

        END IF  !  if firstime


C.......   Compute correct year-to-day conversion factor:

        K = JDATE / 1000        !  year-number
        IF ( MOD( K,4 ) .NE. 0 ) THEN           !  nonleap years
            YR2DAY = 1.0 / 365.0
        ELSE IF ( MOD( K,100 ) .NE. 0 ) THEN    !  noncentury leap years
            YR2DAY = 1.0 / 366.0
        ELSE IF ( MOD( K,400 ) .NE. 0 ) THEN    !  century nonleap years
            YR2DAY = 1.0 / 365.0
        ELSE                                    !  leap centuries
            YR2DAY = 1.0 / 366.0
        END IF


C.......   Compute index correction (offset by 1 because of
C.......   1 + MOD(...) needed below

        HCORR = TZONE + 25

C.......   Compute TMAT:

        DO  133  V = 1, NIPOL
        DO  122  S = 1, NSRC

            IF ( FIND1( S, NDAY, INDD ) .GT. 0 ) THEN

                FAC = 1.0

            ELSE IF ( MDEX( S,V ) .GT. 0 ) THEN

                FAC = MFAC( MON,MDEX( S,V ) ) * WFAC( DAY,WDEX( S,V ) )

            ELSE IF ( WDEX( S,V ) .GT. 0 ) THEN

                FAC = YR2DAY * WFAC( DAY,WDEX( S,V ) )

            ELSE

                FAC = YR2DAY

            END IF

            L = DDEX( S,V )

            DO  111  H = 1, 24

                K = 1 + MOD( H + HCORR - ZONES( S ), 24 )
                TMAT( S,V,H ) = FAC * DFAC( K,L )

111         CONTINUE

122     CONTINUE
133     CONTINUE

        RETURN
        END

