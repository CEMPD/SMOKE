
        SUBROUTINE  MKTMATSI( NSRC,  JDATE, MON, DAY, TZONE, ZONES, TPF,
     &                      MFAC,  MDEX, 
     &                      WFAC,  WDEX, 
     &                      DFAC,  DDEX, 
     &                      NDAY,  INDD, TMAT )

C***********************************************************************
C  subroutine body starts at line  118
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
C       Version   9/96 by CJC for OTAG:  profiles are species-independent
C       Version  11/96 by CJC for OTAG:  new weekday-only normalization
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
C***************************************************************************

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
        INTEGER         TPF  ( NSRC )           ! source temporal treatment flag
        REAL            MFAC( 12, NMCOD )       ! Monthly profile factors
        REAL            WFAC(  7, NWCOD )       ! Weekly    "       "
        REAL            DFAC( 24, NDCOD  )      ! Diurnal   "       "
        INTEGER         MDEX( NSRC )     ! Monthly profile codes
        INTEGER         WDEX( NSRC )     ! weekly  profile codes
        INTEGER         DDEX( NSRC )     ! diurnal profile codes
        INTEGER         NDAY                    ! # of day-specific obs
        INTEGER         INDD( * )     ! subscripts: day-specific obs
        REAL            TMAT( NSRC, 24 ) ! temporal-profile coeffs

C...........   EXTERNAL FUNCTIONS:

        INTEGER         FIND1
        LOGICAL         ISDSTIME        !  true iff daylight savings time( date)

        EXTERNAL        FIND1, ISDSTIME

C...........   Scratch Local variables:

        INTEGER         H               !  hour-loop counter
        INTEGER         K               !  diurnal-profile index
        INTEGER         L               !  diurnal-xref index
        INTEGER         S               !  source counter
        INTEGER         HCORR           !  daylight savings time correction
        REAL            FAC             !  partial matrix factor
        REAL            YR2DAY          !  year-to-day conversion factor
        LOGICAL         EFLAG
        CHARACTER*255   MESG


C.......   SAVED Local Variables

        REAL            MWFAC( 12 )
        DATA            MWFAC / 31.0, 28.0, 31.0, 30.0,
     &                          31.0, 30.0, 31.0, 31.0,
     &                          30.0, 31.0, 30.0, 31.0 /
        REAL            XFAC(  7, NWCOD )       ! Week-day profiles

        LOGICAL         FIRSTIME
        DATA            FIRSTIME / .TRUE. /

        SAVE            MWFAC, XFAC, FIRSTIME


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

13          CONTINUE			!  end loop normalizing month-codes

            DO  23 S = 1, NWCOD

                FAC  = WFAC( 1,S ) + WFAC( 2,S ) + WFAC( 3,S ) +
     &                 WFAC( 4,S ) + WFAC( 5,S )

                IF ( FAC .GT. 0.0 ) THEN

                    FAC = 5.0 / FAC

                ELSE					!  weekend-only profile

                    FAC = WFAC( 6,S ) + WFAC( 7,S )
                    IF ( FAC .GT. 0.0 ) THEN
                        FAC = 2.0 / FAC
                    ELSE
                        FAC = 0.0			!  zero profile
                    END IF

                END IF

                DO  22 K = 1, 7
                    XFAC( K,S ) = FAC * WFAC( K,S )	! for weekday-normalized
                    WFAC( K,S ) = 7.0 * WFAC( K,S )	! for week-normalized
22              CONTINUE

23          CONTINUE		!  end loop normalizing day-of-week codes

C...............     Check cross-reference arrays:

            EFLAG = .FALSE.
            DO  33 S = 1, NSRC
                IF ( MOD( TPF( S ), MTPRFAC ) .EQ. 0 ) THEN
                    IF ( MDEX( S ) .LE. 0 ) THEN
                        WRITE( MESG,94010 ) 
     &                  'Bad MDEX=', MDEX( S ), 'for source', S
                        CALL M3MESG( MESG )
                        EFLAG = .TRUE.
                    END IF
                END IF
                IF ( MOD( TPF( S ), WTPRFAC ) .EQ. 0 ) THEN
                    IF ( WDEX( S ) .LE. 0 ) THEN
                        WRITE( MESG,94010 ) 
     &                  'Bad WDEX=', WDEX( S ), 'for source', S
                        CALL M3MESG( MESG )
                        EFLAG = .TRUE.
                    END IF
                ELSE IF ( MOD( TPF( S ), WDTPFAC ) .EQ. 0 ) THEN
                    IF ( WDEX( S ) .LE. 0 ) THEN
                        WRITE( MESG,94010 ) 
     &                  'Bad WDEX=', WDEX( S ), 'for source', S
                        CALL M3MESG( MESG )
                        EFLAG = .TRUE.
                    END IF
                END IF
33          CONTINUE
            IF ( EFLAG ) THEN
                CALL M3EXIT( 'MKTMATSI', 0, 0, 
     &                       'Bad crosss references', 2 )
            END IF

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

        HCORR = TZONE + 23

C.......   Compute TMAT:

        DO  122  S = 1, NSRC

C.............  Use day-specific data (no adjustments)

            IF ( FIND1( S, NDAY, INDD ) .GT. 0 ) THEN

                FAC = 1.0

C.............  Adjust for year-normal data

            ELSE IF ( MOD( TPF( S ), MTPRFAC ) .EQ. 0 ) THEN

                FAC = MFAC( MON,MDEX( S ) ) * WFAC( DAY,WDEX( S ) )

C.............  Adjust for week-normal data assuming whole week normalizer
                
            ELSE IF ( MOD( TPF( S ), WTPRFAC ) .EQ. 0 ) THEN

               FAC = YR2DAY * WFAC( DAY,WDEX( S ) )

C.............  Adjust for week-normal data assuming week-days normalizer

            ELSE IF ( MOD( TPF( S ), WDTPFAC ) .EQ. 0 ) THEN

                FAC = YR2DAY * XFAC( DAY,WDEX( S ) )

            ELSE

                FAC = YR2DAY

            END IF

            L = DDEX( S )

            DO  111  H = 1, 24

                K = 1 + MOD( H + HCORR - ZONES( S ), 24 )
                TMAT( S,H ) = FAC * DFAC( K,L )

111         CONTINUE

122     CONTINUE

        RETURN
        
94010   FORMAT( 20( A, :, I9, :, 1X ) )

        END

