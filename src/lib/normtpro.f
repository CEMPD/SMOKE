
        SUBROUTINE NORMTPRO

C***********************************************************************
C  subroutine body starts at line 
C
C  DESCRIPTION:
C      Normalizes the temporal profiles for use in generating temporal 
C      emissions.  It allocates memory for and computes weekday normalized 
C      weekly factors as well. It gets an environment variable to determine
C      whether or not to renormalize profiles to 1 or not.
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C       Created M Houyoux 1/99
C
C***********************************************************************
C
C Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
C                System
C File: @(#)$Id$
C
C COPYRIGHT (C) 2002, MCNC Environmental Modeling Center
C All Rights Reserved
C
C See file COPYRIGHT for conditions of use.
C
C Environmental Modeling Center
C MCNC
C P.O. Box 12889
C Research Triangle Park, NC  27709-2889
C
C smoke@emc.mcnc.org
C
C Pathname: $Source$
C Last updated: $Date$ 
C
C****************************************************************************

C.........  MODULES for public variables
C.........  For temporal profiles
        USE MODTMPRL, ONLY: HRLFAC, XWKFAC, WEKFAC, MONFAC,
     &                      NHRL, NWEK, NMON

        IMPLICIT NONE

C...........   EXTERNAL FUNCTIONS:
        LOGICAL   ENVYN
        EXTERNAL  ENVYN

C...........   Local parameters:
        REAL, PARAMETER :: MWFAC( 12 ) =
     &                   ( / 31.0, 28.0, 31.0, 30.0, 31.0, 30.0, 
     &                       31.0, 31.0, 30.0, 31.0, 30.0, 31.0 / )

C...........   SCRATCH LOCAL VARIABLES and their descriptions:

        INTEGER         D, I, K          !  counters and indices

        INTEGER         IOS           !  i/o status

        REAL            DIV           !  scratch divisor
        REAL            FAC           !  scratch factor
        REAL            TOT           !  scratch total

        LOGICAL         RENORM   ! true temporal profiles are to be renormalized

        CHARACTER*300   MESG    !  message buffer

        CHARACTER*16  :: PROGNAME = 'NORMTPRO'  ! program name

C***********************************************************************
C   begin body of subroutine  NORMTPRO

C.........  Get information from the environment about renormalization
        MESG = 'Renormalize temporal profiles'
        RENORM = ENVYN ( 'RENORM_TPROF', MESG, .TRUE., IOS )

C.........  Allocate memory for weekday-normalized weekly emissions
        ALLOCATE( XWKFAC( 7,NWEK ), STAT=IOS )
        CALL CHECKMEM( IOS, 'XWKFAC', PROGNAME )

        XWKFAC = 1.0  ! Array
 
C.........  Renormalize temporal profiles, if needed
        IF( RENORM ) THEN

            DO I = 1, NMON   ! Monthly

                TOT = 0.0
                DO K = 1, 12
                    TOT = TOT + MONFAC( K,I )
                END DO

                IF( TOT .GT. 0. ) THEN
                    DIV = 1.0 / TOT
                ELSE
                    DIV = 0.
                END IF

                DO K = 1, 12
                    MONFAC( K,I ) = DIV * MONFAC( K,I )
                END DO

            END DO

            DO I = 1, NWEK   ! Weekly

                TOT = 0.0
                DO K = 1, 7
                    TOT = TOT + WEKFAC( K,I )
                END DO

                IF( TOT .GT. 0. ) THEN
                    DIV = 1.0 / TOT
                ELSE
                    DIV = 0.
                END IF

                DO K = 1, 7
                    WEKFAC( K,I ) = DIV * WEKFAC( K,I )
                END DO

            END DO

            DO D = 1, 7

                DO I = 1, NHRL   ! Diurnal

                    TOT = 0.0
                    DO K = 1, 24
                        TOT = TOT + HRLFAC( K,I,D )
                    END DO

                    IF( TOT .GT. 0. ) THEN
                        DIV = 1.0 / TOT
                    ELSE
                        DIV = 0.
                    END IF

                    DO K = 1, 24
                        HRLFAC( K,I,D ) = DIV * HRLFAC( K,I,D )
                    END DO

                END DO

            END DO

        END IF

C.........  Weight monthly profiles by the number of days in a month
        DO I = 1, NMON

            FAC = 0.0
            TOT = 0.0
            DO K = 1, 12
                TOT = TOT + MONFAC( K,I )
                FAC = FAC + MWFAC( K ) * MONFAC( K,I )
            END DO

            IF( FAC .GT. 0 ) FAC = TOT / FAC

            DO K = 1, 12
                MONFAC( K,I ) = FAC * MONFAC( K,I )
            END DO

        END DO    !  end loop normalizing month-codes

C.........  Weight weekly profiles for weekly and day-of-week adjustments
        DO I = 1, NWEK

            FAC  = WEKFAC( 1,I ) + WEKFAC( 2,I ) + WEKFAC( 3,I ) +
     &             WEKFAC( 4,I ) + WEKFAC( 5,I )

            IF ( FAC .GT. 0.0 ) THEN
 
                FAC = 5.0 / FAC
 
            ELSE                                 !  weekend-only profile

                FAC = WEKFAC( 6,I ) + WEKFAC( 7,I )
                IF ( FAC .GT. 0.0 ) THEN
                    FAC = 2.0 / FAC
                ELSE
                    FAC = 0.0                    !  zero profile
                END IF
 
            END IF

            DO K = 1, 7
                XWKFAC( K,I ) = FAC * WEKFAC( K,I )  ! for weekday-normalized
                WEKFAC( K,I ) = 7.0 * WEKFAC( K,I )  ! for week-normalized
            END DO

        END DO            !  end loop normalizing day-of-week code

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

        END SUBROUTINE NORMTPRO

