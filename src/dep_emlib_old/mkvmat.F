
        SUBROUTINE  MKVMAT( NSRC, NETYPE, JDATE, MON, DAY, TZONE, ZONES,
     &                      MFAC,  NMON, MDEX, 
     &                      WFAC,  NWEK, WDEX, 
     &                      DFAC,  NDAY, DDEX, VMAT, TMAT )

C***********************************************************************
C  subroutine body starts at line  116
C
C  FUNCTION:
C       Construct temporal-coefficient-transform matrices for 
C       program TMPMOBIL
C
C  PRECONDITIONS REQUIRED:
C       Temporal profile arrays for monthly, weekly, diurnal profiles.
C       MDEX, WDEX entries set to zero for month- or week-independent
C       source records.
C
C  SUBROUTINES AND FUNCTIONS CALLED:  none
C
C  REVISION  HISTORY:
C       Prototype 1/96 by CJC adapted from point-source MKTMAT
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
C****************************************************************************

      IMPLICIT NONE

C...........   INCLUDES:

        INCLUDE 'TMDIMS3.EXT'   !  temporal-allocation dimensioning parameters
        INCLUDE 'MBDIMS3.EXT'   !  mobile-source       dimensioning parameters


C...........   ARGUMENTS and their descriptions:

        INTEGER    NSRC                    ! number of sources
        INTEGER    NETYPE                  ! number of emission types
        INTEGER    NMON                    ! number of monthly profiles
        INTEGER    NWEK                    ! number of weekly profiles
        INTEGER    NDAY                    ! number of diurnal profiles
        INTEGER    JDATE                   ! date YYYYDDD
        INTEGER    MON                     ! 1 ... 12
        INTEGER    DAY                     ! 1 ... 7
        INTEGER    TZONE                   ! time zone (5 for Eastern)
        INTEGER    ZONES( NSRC )           ! source time zones
        REAL       MFAC( 12, NMON )        ! Monthly profile factors
        REAL       WFAC(  7, NWEK )        ! Weekly    "       "
        REAL       DFAC( 24, NDAY )        ! Diurnal   "       "
        INTEGER*2  MDEX( NSRC, NVTYPE, MXEMIS )    ! Monthly profile codes
        INTEGER*2  WDEX( NSRC, NVTYPE, MXEMIS )    ! weekly  profile codes
        INTEGER*2  DDEX( NSRC, NVTYPE, MXEMIS )    ! diurnal profile codes
        REAL       VMAT( NSRC, NVTYPE )            ! VMT adjusted by VMTmix
        REAL       TMAT( NSRC, NVTYPE, MXEMIS, 24 ) ! temporal-profile coeffs

C...........   EXTERNAL FUNCTIONS:

        INTEGER         FIND1
        LOGICAL         ISDSTIME        !  true iff daylight savings time( date)

        EXTERNAL        FIND1, ISDSTIME

C...........   Scratch Local variables:

        INTEGER         H               !  hour-loop counter
        INTEGER         I               !  counter
        INTEGER         K               !  diurnal-profile index
        INTEGER         L               !  diurnal-xref index
        INTEGER         S               !  source counter
        INTEGER         V               !  chemistry variable counter
        INTEGER         E               !  emissions mode counter
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
        LOGICAL         EFLAG

        SAVE            MWFAC, FIRSTIME, YR2DAY

        CHARACTER*256   MESG

C***********************************************************************
C   begin body of subroutine  MKVMAT

C.......   Firstime:  renormalize MFAC, WFAC

        IF ( FIRSTIME ) THEN

            FIRSTIME = .FALSE.
            EFLAG    = .FALSE.

C...........   Compute correct year-to-day conversion factor:

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

            DO  13 I = 1, NMON

                FAC = 0.0
                DO  11 K = 1, 12
                    FAC = FAC + MWFAC( K ) * MFAC( K,I )
11              CONTINUE

                IF( FAC .LE. 0. ) THEN
                    EFLAG = .TRUE.
                    WRITE( MESG,94010 )
     &                  'Monthly profiles at count', I, 
     &                  'found to be zero!'
                    CALL M3MESG( MESG )

                ELSE
                    FAC = 1.0 / ( FAC * YR2DAY )

                ENDIF

                DO  12 K = 1, 12
                    MFAC( K,I ) = FAC * MFAC( K,I )
12              CONTINUE

13          CONTINUE

            IF( EFLAG ) THEN
                MESG = 'Error renormalizing temporal profiles'
                CALL M3EXIT( 'MKVMAT', 0, 0, MESG, 2 )
            ENDIF

            DO  23 I = 1, NWEK
            DO  22 K = 1, 7
                WFAC( K,I ) = 7.0 * WFAC( K,I )
22          CONTINUE
23          CONTINUE


        END IF  !  if firstime


C.......   Compute index correction (offset by 1 because of
C.......   1 + MOD(...) needed below

        HCORR = TZONE + 25

C.......   Compute TMAT:

        DO  144  E = 1, NETYPE
        DO  133  V = 1, NVTYPE
        DO  122  S = 1, NSRC

            IF ( MDEX( S,V,E ) .GT. 0 ) THEN

                FAC = VMAT( S,V ) * MFAC( MON,MDEX( S,V,E ) ) 
     &                            * WFAC( DAY,WDEX( S,V,E ) )

            ELSE IF ( WDEX( S,V,E ) .GT. 0 ) THEN

                FAC = VMAT( S,V ) * WFAC( DAY,WDEX( S,V,E ) )

            ELSE

                FAC = VMAT( S,V )

            END IF

            L = DDEX( S,V,E )

            DO  111  H = 1, 24

                K = 1 + MOD( H + HCORR - ZONES( S ), 24 )
                TMAT( S,V,E,H ) = FAC * DFAC( K,L )

111         CONTINUE

122     CONTINUE
133     CONTINUE
144     CONTINUE

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10 ( A, :, I5, :, 2X ) )

        END

