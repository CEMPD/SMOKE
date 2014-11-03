
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
C     Created M Houyoux 1/99
C
C     Version 07/2014 by C.Coats for  new GENTPRO CSV profiles and
C     cross-references.  Be sure to normalize "profile-zero" case
C     used to handle the "profile-not-supplied" cases.
C
C***********************************************************************
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
C           For temporal profiles

        USE MODTMPRL, ONLY: HRLFAC, XWKFAC, WEKFAC, MONFAC, MONFAC_ORG,
     &                      DOMFAC, NHRL, NWEK, NMON, NDOM

        IMPLICIT NONE

C.........  INCLUDES

        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters

C.........   EXTERNAL FUNCTIONS:

        LOGICAL, EXTERNAL :: ENVYN

C.........   Local parameters:

        REAL, PARAMETER :: MWFAC( 12 ) = FLOAT( MON_DAYS )

        CHARACTER(16), PARAMETER :: PROGNAME = 'NORMTPRO'  ! program name

C.........   SCRATCH LOCAL VARIABLES and their descriptions:

        INTEGER         I, J, K  !  counters and indices
        INTEGER         IOS         !  i/o status

        REAL            DIV         !  scratch divisor
        REAL(8)         FAC         !  scratch factor
        REAL(8)         TOT         !  scratch total

        LOGICAL       :: RENORM = .FALSE.  ! true temporal profiles are to be renormalized

        CHARACTER(300)  MESG    !  message buffer

C***********************************************************************
C   begin body of subroutine  NORMTPRO

C.........  Get information from the environment about renormalization
        MESG = 'Renormalize temporal profiles'
        RENORM = ENVYN ( 'RENORM_TPROF', MESG, .TRUE., IOS )

C.........  Allocate memory for weekday-normalized weekly emissions

        ALLOCATE( XWKFAC(  7,0:NWEK ), 
     &        MONFAC_ORG( 12,0:NMON ), STAT=IOS )
        CALL CHECKMEM( IOS, 'XWKFAC,MONFAC_ORIG', PROGNAME )

        XWKFAC = 1.0  ! Array
 
C.........  Renormalize temporal profiles, if needed

        IF( RENORM ) THEN

            DO I = 0, NMON   ! Monthly

                TOT = 0.0D0
                DO K = 1, 12
                    TOT = TOT + MONFAC( K,I )
                END DO

                IF( TOT .GT. 0.0D0 ) THEN
                    DIV = 1.0D0 / TOT
                ELSE
                    DIV = 0.0
                END IF

                DO K = 1, 12
                    MONFAC( K,I ) = DIV * MONFAC( K,I )
                END DO

            END DO          !  end loop on month-profiles I

            DO I = 0, NWEK   ! Weekly

                TOT = 0.0D0
                DO K = 1, 7
                    TOT = TOT + WEKFAC( K,I )
                END DO

                IF( TOT .GT. 0.0D0 ) THEN
                    DIV = 1.0D0 / TOT
                ELSE
                    DIV = 0.0D0
                END IF

                DO K = 1, 7
                    WEKFAC( K,I ) = DIV * WEKFAC( K,I )
                END DO

            END DO          !  end loop on week-profiless I

            DO I = 0, NHRL   ! Diurnal

                TOT = 0.0D0
                DO K = 1, 24
                    TOT = TOT + HRLFAC( K,I )
                END DO

                IF( TOT .GT. 0. ) THEN
                    DIV = 1.0 / TOT
                ELSE
                    DIV = 0.0
                END IF

                DO K = 1, 24
                    HRLFAC( K,I ) = DIV * HRLFAC( K,I )
                END DO

            END DO          !  end loop on diurnal profiles I

        END IF

C.........  Backup original MONFAC to MONFAC_ORG before applying monthly adjustment

        IF ( NMON .GT. 0 )   MONFAC_ORG = MONFAC

C.........  Weight monthly profiles by the number of days in a month

        DO I = 0, NMON

            FAC = 0.0D0
            TOT = 0.0D0
            DO K = 1, 12
                TOT = TOT + MONFAC( K,I )
                FAC = FAC + MWFAC( K ) * MONFAC( K,I )
            END DO

            IF( FAC .GT. 0.0D0 ) FAC = TOT / FAC

            DO K = 1, 12
                MONFAC( K,I ) = FAC * MONFAC( K,I )
            END DO

        END DO    !  end loop normalizing month-codes

C.........  Weight weekly profiles for weekly and day-of-week adjustments

        DO I = 0, NWEK

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

C.........  Normalize day-of-month profiles:
        
        DO I = 0, NDOM

            DO J = 1, 12

                FAC = 0.0D0
                DO K = 1, MON_DAYS(J)
                    FAC = FAC + DOMFAC( K,J,I )
                END DO

                IF ( FAC .GT. 0.0D0 ) THEN      !  profile available
                    DIV = 1.0D0 / FAC
                ELSE                            !  profile was missing:  zero it out
                    DIV = 0.0
                END IF

                DO K = 1, MON_DAYS(J)
                    DOMFAC( K,J,I ) = DIV * DOMFAC( K,J,I )
                END DO

            END DO      !  end loop on months J
            
            !  Hack for February:
            
            IF( DOMFAC( 29,2,I ) .LE. 0.0 ) DOMFAC( 29,2,I ) = DOMFAC( 28,2,I )

        END DO      !  end loop on day-of-month profiles I

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

        END SUBROUTINE NORMTPRO

