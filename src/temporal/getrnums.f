
        SUBROUTINE GETRNUMS

C***********************************************************************
C
C  DESCRIPTION:
C     Subroutine GETPROBS generates the random numbers associated with
C     the specified uncertainty sources
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C     Created 3/02 by G. Cano
C
C****************************************************************************/
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
C***************************************************************************

C.........  MODULES for public variables
C...........   This module contains the uncertainty arrays
        USE MODUNCERT

        IMPLICIT NONE

C..........  EXTERNAL FUNCTIONS and their descriptions:
        INTEGER         ENVINT
        EXTERNAL    ENVINT


C.........  Local variables
        INTEGER       H, I, J, K, S, T   ! counters and indices
        INTEGER       IOS                ! i/o status
        INTEGER       ISEED              ! seed for random number generator

        LOGICAL    :: EFLAG = .FALSE.    ! true: error found
        LOGICAL    :: FIRSTIME = .TRUE.  ! true: first time the routine called

        CHARACTER*300 MESG                    ! message buffer
        CHARACTER*16 :: PROGNAME = 'GETRNUMS' ! program name

C***********************************************************************
C   begin body of subroutine GETRNUMS

        IF ( FIRSTIME ) THEN

            ISEED = ENVINT( 'SMK_RANDOMSEED', 
     &                  'Seed for random number genration', 27195, IOS )
            CALL RANDOM_SEED( ISEED )

        END IF

        DO I = 1, UNSRC
            DO J = 1, UNIPPA

               SELECT CASE( METHOD( I, J ) )

               CASE( 0 ) ! Empirical

                   SELECT CASE( EPTYP( I, J ) ) ! Empirical statistics
        
                   CASE( 1 ) ! Stepwise 
        
                   CASE( 2 ) ! Linear Interpolation 

                   CASE DEFAULT
C....................  Error case.  invalid method entry
                       MESG = 'Invalid empirical Method entry'
                       CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                   END SELECT ! end type select

C .............  Error case.  invalid entry                    

               CASE( 1 ) ! Parametric

                   SELECT CASE( EPTYP( I, J ) ) ! Parametric distribution
        
                   CASE( 1 ) ! Normal
        
                       SAMPGEN( I,J ) = RNORM( PARMS( UNCIDX( I,J ), 1),
     &                                         PARMS( UNCIDX( I,J ), 2))

                   CASE( 2 ) ! Lognormal

                       SAMPGEN( I,J ) = EXP(RNORM(PARMS(UNCIDX(I,J),1),
     &                                            PARMS(UNCIDX(I,J),2))) 

                   CASE( 3 ) ! Gamma
C BUG................  code exists but not linked in
                       SAMPGEN( I,J ) = RGAMMA(PARMS(UNCIDX(I,J), 1),
     &                                         PARMS(UNCIDX(I,J), 2) )

                   CASE( 4 ) ! Weibull
C BUG................  code exists but not linked in
                       SAMPGEN( I,J ) = RGAMMA(PARMS(UNCIDX(I,J), 1),
     &                                         PARMS(UNCIDX(I,J), 2) )

                   CASE( 5 ) ! Beta
C BUG................  code exists but not linked in
                       SAMPGEN( I,J ) = RGAMMA(PARMS(UNCIDX(I,J), 1),
     &                                         PARMS(UNCIDX(I,J), 2) )

                   CASE DEFAULT
C....................  Error case.  invalid method entry
                       MESG = 'Invalid parametric Method entry'
                       CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                   END SELECT ! end type select

               CASE DEFAULT
C................  Error case.  invalid type entry
                   MESG = 'Invalid Type entry'
                   CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
               END SELECT ! end method select

            END DO ! J loop
        END DO ! I loop

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C.........  Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

C.........  test buffering formats............ 94xxx

98010   FORMAT( 10( A, :, F8.4, :, 1X ) )

C******************  INTERNAL SUBPROGRAMS  *****************************

        CONTAINS

            REAL FUNCTION RNORM( MU, SIGMA )

C.............  function parameters
            REAL,             INTENT(IN) :: MU     ! distribution Mean 
            REAL,             INTENT(IN) :: SIGMA  ! distribution deviation

C.............  Local variables
            LOGICAL,     SAVE :: X2_SET = .FALSE.! Normal deviate saved

            REAL                            FAC    ! -2*(log(RSQ)/RSQ)
            REAL                         :: NEARZERO = 0.001! value close to 0 
            REAL                            RSQ    ! 
            REAL                            U1, U2 ! uniform (0,1) RNs
            REAL                            V1, V2 ! RNs in the unit circle

C**********************************************************************

            IF ( SIGMA .LE. NEARZERO ) THEN

                CALL M3EXIT( PROGNAME, 0, 0, 'Negative Normal Sigma ' //
     &                       'parameter.  Should be > 0', 2 )
            END IF

            RSQ = 0.0
            DO WHILE ( RSQ .EQ. 0.0 .OR. RSQ .GT. 1.0 ) 
                CALL RANDOM_NUMBER( U1 )
                V1 = 2.0 * U1 - 1.0
                CALL RANDOM_NUMBER( U2 )
                V2 = 2.0 * U2 - 1.0
                RSQ = V1*V1 + V2*V2
            END DO 
            
            FAC = SQRT( -2.0 * log( RSQ ) / RSQ )
            RNORM = (V1 * FAC) * SIGMA + MU  ! transform N(0,1) to N(MU,SIGMA)
            
            RETURN

            END FUNCTION RNORM

C**********************************************************************
            REAL FUNCTION RGAMMA( ALPHA, BETA )

C.............  function parameters
            REAL,             INTENT(IN) :: ALPHA  ! distribution Mean 
            REAL,             INTENT(IN) :: BETA   ! distribution Mean 

C.............  Local variables
            REAL                         :: NEARZERO = 0.001 ! near 0 constant
            REAL                            U      ! uniform (0,1) RNs

C**********************************************************************

C BUG: needs to be tested

            IF ( ALPHA .LT. NEARZERO .AND. ALPHA .GT. - NEARZERO .OR.
     &           BETA .LT. NEARZERO .AND. BETA .GT. - NEARZERO ) THEN

                CALL M3EXIT( PROGNAME, 0, 0, 'Negative Alpha or Beta '//
     &                       'parameters.  Both should be > 0', 2 )

            END IF


c            IF ( ALPHA .LT. 1.0 + NEARZERO .AND. 
c     &           ALPHA .GT. 1.0 - NEARZERO ) THEN ! 
C BUG: insert Gamma code here
                 RGAMMA = REXPONENTIAL( BETA )

c            ELSE
c compute a Gamma deviate
c            END IF

            RETURN

            END FUNCTION RGAMMA

C**********************************************************************
            REAL FUNCTION REXPONENTIAL( BETA )


C.............  function parameters
            REAL,             INTENT(IN) :: BETA   ! distribution Mean 

C.............  Local variables
            REAL                            U      ! uniform (0,1) RNs

C**********************************************************************

C BUG: needs to be tested

            IF ( BETA .LE. 0.0 ) THEN
                CALL M3EXIT( PROGNAME, 0, 0, 'Negative Exponetial ' //
     &                       'parameter.  Should be > 0', 2 )
            END IF

            U = 0.0
            DO WHILE ( U .LE. 0.0 ) 
                CALL RANDOM_NUMBER( U )
            END DO 

            REXPONENTIAL = -1.0 * BETA * LOG( 1.0 - U )

            RETURN

            END FUNCTION REXPONENTIAL

C**********************************************************************

            SUBROUTINE RSAMPLE

C.............  Local variables
            CHARACTER*9, PARAMETER :: PDF( 5 ) = ( / 'NORMAL   ',
     &                                               'LOGNORMAL',
     &                                               'GAMMA    ',
     &                                               'WEIBULL  ',
     &                                               'BETA     ' / )

            INTEGER              I, J      !  counters
            INTEGER       :: ISEED = 27195 !  seed for random number (RN)

            REAL, ALLOCATABLE :: RNUMS( : )!  array for uncertainty RNs

            REAL                 S, X      !  temporary variables
            REAL                 RBETA     !  temporary Beta RN
            REAL                 REXP      !  temporary Exponential RN
            REAL                 RGAMMA    !  temporary Gamma RN
            REAL                 RWEIBULL  !  temporary Weibull RN
            REAL              :: MU = 0.5  !  dummy mean
            REAL              :: SIGMA=1.2 !  dummy variance
            REAL,        SAVE :: URAND01   !  quasi Uniform RN

C**********************************************************************


C.............  Get Weibull random number
            CALL RANDOM_NUMBER( URAND01 )    !  get quasi uniform (0,1) RN
            REXP = -LOG( URAND01 ) * MU      !  transform to an exponential RN
            RWEIBULL = REXP ** ( 1 / SIGMA ) !  tranform exp to a weibull RN

C.............  Get Gamma random number 
            X = 1.0
            DO I = 1, 5
                CALL RANDOM_NUMBER( URAND01 )    !  get quasi uniform (0,1) RN
                X = X * URAND01
            END DO
            RGAMMA = REXP ** ( 1 / SIGMA)    !  transform to a Gamma RN

C.............  Get Beta random number 
            RBETA = RGAMMA
            DO I = 1, 6                      !  Get a Gamma RN
                CALL RANDOM_NUMBER( URAND01 )    !  get quasi uniform (0,1) RN
                RGAMMA = REXP ** ( 1 / SIGMA)    !  transform to a Gamma RN
            END DO
            RBETA = RBETA / ( RBETA + RGAMMA )

C.............  Get Beta random number 
            RBETA = RGAMMA
            DO I = 1, 6                      !  Get a Gamma RN
                CALL RANDOM_NUMBER( URAND01 )    !  get quasi uniform (0,1) RN
                RGAMMA = REXP ** ( 1 / SIGMA)    !  transform to a Gamma RN
            END DO
            RBETA = RBETA / ( RBETA + RGAMMA )

            END SUBROUTINE RSAMPLE

        END SUBROUTINE GETRNUMS
