
        SUBROUTINE  MKTMAT( NSRC, JDATE, MONTH, DAYOW, TZONE, ZONES,
     &                      TPF, MFAC, NMON, MDEX, WFAC, NWEK, WDEX,
     &                      DFAC, NDAY, DDEX, NDSP, INDD, TMAT       )

C***********************************************************************
C  subroutine body starts at line  119
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
C               substitution of day-specific data NDSP,INDD
C       Revised   4/97 by MRH for SMRAQ: weekly/monthly account for time zone
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
C***********************************************************************

      IMPLICIT NONE

C...........   INCLUDES:

        INCLUDE 'TMDIMS3.EXT'   !  temporal-allocation dimensioning parameters
        INCLUDE 'CHDIMS3.EXT'   !  chemical-species    dimensioning parameters


C...........   ARGUMENTS and their descriptions:

        INTEGER         NSRC                    ! number of sources
        INTEGER         JDATE                   ! date YYYYDDD
        INTEGER         MONTH( 24, 0:23 )       ! source time zone's 1 ... 12
        INTEGER         DAYOW( 24, 0:23 )       ! source time zone's 1 ... 7
        INTEGER         TZONE                   ! time zone (5 for Eastern)
        INTEGER         ZONES( NSRC )           ! source time zones
        INTEGER         TPF  ( NSRC )           ! source temporal treatment flag
        INTEGER         NMON                    ! actual number of monthly profs
        INTEGER         NWEK                    ! actual number of weekly profs
        INTEGER         NDAY                    ! actual number of diurnal profs
        REAL            DFAC( 24, NDAY )        ! Diurnal   "       "
        REAL            MFAC( 12, NMON )        ! Monthly profile factors
        REAL            WFAC(  7, NWEK )        ! Weekly    "       "
        INTEGER         MDEX( NSRC, NIPOL )     ! Monthly profile codes
        INTEGER         WDEX( NSRC, NIPOL )     ! weekly  profile codes
        INTEGER         DDEX( NSRC, NIPOL )     ! diurnal profile codes
        INTEGER         NDSP                    ! # of day-specific obs
        INTEGER         INDD( * )     ! subscripts: day-specific obs
        REAL            TMAT( NSRC, NIPOL, 24 ) ! temporal-profile coeffs

C...........   EXTERNAL FUNCTIONS:

        INTEGER         FIND1
        LOGICAL         ISDSTIME       !  true iff daylight savings time( date)
        REAL            YR2DAY

        EXTERNAL        FIND1, ISDSTIME, YR2DAY

C...........   Scratch Local variables:

        INTEGER         DAY             !  day for source and hour pointer
        INTEGER         H               !  hour-loop counter
        INTEGER         HCORR           !  daylight savings time correction
        INTEGER         K               !  diurnal-profile index
        INTEGER         L               !  diurnal-xref index
        INTEGER         MON             !  month for source and hour pointer
        INTEGER         S               !  source counter
        INTEGER         V               !  chemistry variable counter

        REAL            FAC             !  partial matrix factor
        REAL            XFAC( 7, NWCOD) !  Week-day profiles
        REAL            YRFAC           !  year to day factor


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

            DO  13 S = 1, NMON

                FAC = 0.0
                DO  11 K = 1, 12
                    FAC = FAC + MWFAC( K ) * MFAC( K,S )
11              CONTINUE

                FAC = 1.0 / FAC

                DO  12 K = 1, 12
                    MFAC( K,S ) = FAC * MFAC( K,S )
12              CONTINUE

13          CONTINUE    !  end loop normalizing month-codes

            DO  23 S = 1, NWEK

                FAC  = WFAC( 1,S ) + WFAC( 2,S ) + WFAC( 3,S ) +
     &                 WFAC( 4,S ) + WFAC( 5,S )

                IF ( FAC .GT. 0.0 ) THEN
 
                    FAC = 5.0 / FAC
 
                ELSE                                 !  weekend-only profile

                    FAC = WFAC( 6,S ) + WFAC( 7,S )
                    IF ( FAC .GT. 0.0 ) THEN
                        FAC = 2.0 / FAC
                    ELSE
                        FAC = 0.0                    !  zero profile
                    END IF
 
                END IF

                DO  22 K = 1, 7
                    XFAC( K,S ) = FAC * WFAC( K,S )  ! for weekday-normalized
                    WFAC( K,S ) = 7.0 * WFAC( K,S )  ! for week-normalized
22              CONTINUE

23          CONTINUE            !  end loop normalizing day-of-week code

        END IF  !  if firstime


C.......   Compute correct year-to-day conversion factor:

        YRFAC = YR2DAY( JDATE / 1000 )

C.......   Compute index correction (offset by 1 because of
C.......   1 + MOD(...) needed below

        HCORR = TZONE + 23

C.......   Compute TMAT:

        DO  133  V = 1, NIPOL
        DO  122  S = 1, NSRC

            L = DDEX( S,V )

C.............  Use day-specific data (no adjustments)
            IF ( FIND1( S, NDSP, INDD ) .GT. 0 ) THEN

                DO 44   H = 1, 24

                    K = 1 + MOD( H + HCORR - ZONES( S ), 24 )
                    TMAT( S,V,H ) = DFAC( K,L )

44              CONTINUE

C.............  Adjust for year-normal data
            ELSE IF ( MOD( TPF( S ), MTPRFAC ) .EQ. 0 ) THEN

                DO 55   H = 1, 24

                    MON = MONTH( H, ZONES( S ) )
                    DAY = DAYOW( H, ZONES( S ) )
                    FAC = MFAC( MON,MDEX( S,V ) ) * 
     &                    WFAC( DAY,WDEX( S,V ) )
                    K = 1 + MOD( H + HCORR - ZONES( S ), 24 )
                    TMAT( S,V,H ) = FAC * DFAC( K,L )

55              CONTINUE

C.............  Adjust for week-normal data assuming whole week normalizer
            ELSE IF ( MOD( TPF( S ), WTPRFAC ) .EQ. 0 ) THEN

                DO 66   H = 1, 24

                    DAY = DAYOW( H, ZONES( S ) )
                    FAC = YRFAC * WFAC( DAY,WDEX( S,V ) )
                    K = 1 + MOD( H + HCORR - ZONES( S ), 24 )
                    TMAT( S,V,H ) = FAC * DFAC( K,L )

66              CONTINUE

C.............  Adjust for week-normal data assuming week-days normalizer
            ELSE IF ( MOD( TPF( S ), WDTPFAC ) .EQ. 0 ) THEN

                DO 77 H = 1, 24
 
                    DAY = DAYOW( H, ZONES( S ) )
                    FAC = YRFAC * XFAC( DAY,WDEX( S,V ) )
                    K = 1 + MOD( H + HCORR - ZONES( S ), 24 )
                    TMAT( S,V,H ) = FAC * DFAC( K,L )
 
77              CONTINUE

            ELSE

                DO 88 H = 1, 24

                    K = 1 + MOD( H + HCORR - ZONES( S ), 24 )
                    TMAT( S,V,H ) = YRFAC * DFAC( K,L )

88              CONTINUE

            END IF

122     CONTINUE
133     CONTINUE

        RETURN
        END

