
        SUBROUTINE GETRNUMS( ITIME, ISRC, IPA )

C***********************************************************************
C
C  DESCRIPTION:
C     Subroutine GETRNUMS generates the random numbers associated with
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

C...........   INCLUDES
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters

C..........  EXTERNAL FUNCTIONS and their descriptions:
        INTEGER         ENVINT

        EXTERNAL    ENVINT

C...........   SUBROUTINE ARGUMENTS
        INTEGER     , INTENT (IN)    :: IPA       ! index of pollutant of activity
        INTEGER     , INTENT (IN)    :: ISRC      ! source index in uncertainty array
        INTEGER     , INTENT (IN)    :: ITIME     ! current processing time

C.........  Local variables
        INTEGER       IOS                         ! i/o status
        INTEGER       ISEED                       ! seed for RN generator

        LOGICAL    :: FIRSTIME = .TRUE.  ! true: first time the routine called

        CHARACTER*300 MESG                    ! message buffer
        CHARACTER*16 :: PROGNAME = 'GETRNUMS' ! program name

C***********************************************************************
C   begin body of subroutine GETRNUMS

        IF ( FIRSTIME ) THEN

            ISEED = ENVINT( 'SMK_RANDOMSEED', 
     &                  'Seed for random number genration', 27195, IOS )
            CALL RANDOM_SEED( ISEED )

            FIRSTIME = .FALSE.

        END IF

C.........  Determine Parametric or epirical distributoin
        SELECT CASE( METHOD( ISRC, IPA ) )
        
C.........  Empirical sampling
        CASE( 0 )
        
C.............  Stepwise sampling method
            IF ( EPTYP( ISRC, IPA ) .EQ. 6 ) THEN 

                SAMPGEN( ISRC,IPA ) = REMPSTEP( APRCH( ISRC,IPA ),
     &                                          ISRC, IPA         )

C.............  Linear interpolation sampling method
            ELSE IF ( EPTYP( ISRC, IPA ) .EQ. 7 ) THEN 

c                SAMPGEN( ISRC,IPA ) = RLINTRPL( APRCH( ISRC,IPA ),
                SAMPGEN( ISRC,IPA ) = REMPSTEP( APRCH( ISRC,IPA ),
     &                                          ISRC, IPA         )

C.............  Error case.  invalid empirical method
            ELSE

                WRITE( MESG, 94010)
     &             'Invalid Empirical type for source ',
     &             ISRC, 'and pollutant/activity ', IPA
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

            END IF
        
C.........  Parametric sampling method        
        CASE( 1 )

C.............  Distributions to sample from
            SELECT CASE( EPTYP( ISRC,IPA ) )

C.............  Normal random deviate needed
            CASE( 1 ) 
        
                SAMPGEN( ISRC,IPA ) = 
     &                   RNORM( PARMS( UNCIDX( ISRC,IPA ), 1 ),
     &                          PARMS( UNCIDX( ISRC,IPA ), 2 ) )
   

            CASE( 2 ) ! Lognormal
        
                SAMPGEN( ISRC,IPA ) = 
     &                   EXP( RNORM( PARMS( UNCIDX( ISRC,IPA ), 1 ),
     &                               PARMS( UNCIDX( ISRC,IPA ), 2) ) )
  
            CASE( 3 ) ! Gamma
C.........  known BUG: code not tested
                SAMPGEN( ISRC,IPA ) = 
     &                   RGAMMA( PARMS( UNCIDX( ISRC,IPA ), 1 ),
     &                           PARMS( UNCIDX( ISRC,IPA ), 2 ) )
        
            CASE( 4 ) ! Weibull
C.........  known BUG: code not tested
                SAMPGEN( ISRC,IPA ) = 
     &                   RWEIBULL( PARMS( UNCIDX( ISRC,IPA ), 1 ),
     &                             PARMS( UNCIDX( ISRC,IPA ), 2 ) )
 
            CASE( 5 ) ! Beta
C.........  known BUG: code not tested
                SAMPGEN( ISRC,IPA ) = 
     &                   RBETA( PARMS( UNCIDX( ISRC,IPA ), 1 ),
     &                          PARMS( UNCIDX( ISRC,IPA ), 2 ) )
   
            CASE DEFAULT
C.............  Error case.  invalid method entry
                WRITE( MESG, 94010)
     &             'Invalid Parametric type for source ',
     &             ISRC, 'and pollutant/activity ', IPA
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END SELECT ! end type select
        
C.........  Error case.  invalid type entry
        CASE DEFAULT
            WRITE( MESG, 94010)
     &             'Invalid Methods type for source ', 
     &             ISRC, 'and pollutant/activity ', IPA
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END SELECT ! end method select

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C.........  Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

C.........  test buffering formats............ 94xxx

98010   FORMAT( 10( A, :, F8.4, :, 1X ) )

C******************  INTERNAL SUBPROGRAMS  *****************************

        CONTAINS

C**********************************************************************

            REAL FUNCTION REMPSTEP( APR , S, PA )

C.............  function parameters
            INTEGER,          INTENT(IN) :: APR
            INTEGER,          INTENT(IN) :: S
            INTEGER,          INTENT(IN) :: PA

C.............  Local variables
            INTEGER                         I     ! counters
            INTEGER                         IOS   ! i/o status
            INTEGER                         IDX   ! index for empirical dist.

            INTEGER, SAVE                :: LTIME ! used to track changing time

            LOGICAL, SAVE                :: FIRSTIME = .TRUE. ! first call

            REAL                            U      ! uniform (0,1) RNs
	    REAL, SAVE, ALLOCATABLE      :: SSAMP ( : )  ! holds S method RNs
	    REAL, SAVE, ALLOCATABLE      :: STSAMP ( : ) ! holds ST method RNs 

C**********************************************************************

            IF( FIRSTIME ) THEN

                LTIME = ITIME 

                ALLOCATE( SSAMP( NUCPCKT ), STAT=IOS )
                CALL CHECKMEM( IOS, 'SSAMP', PROGNAME )
                SSAMP = AMISS3  ! array

                ALLOCATE( STSAMP( NUCPCKT ), STAT=IOS )
                CALL CHECKMEM( IOS, 'STSAMP', PROGNAME )
                STSAMP = AMISS3  ! array

                FIRSTIME = .FALSE.

            END IF

C.............  Unconditionally sample if no value exists
C               otherwise determine approach for sampling.
C.............  Approach for sampling means what the conditions 
C               are for drawing a random number given a packet
C               for a source and hour.  The condition will be either
C               the same or diffferent for each

C.............  (S) Sources-same, hourly-same
            IF( APR .EQ. 1 ) THEN

                IF( SSAMP( UNCIDX( S,PA ) ) .EQ. AMISS3 ) THEN

                    CALL RANDOM_NUMBER( U )
                    IDX = INT( U*NUMEP( S,PA ) + 1 )
                    SSAMP( UNCIDX( S,PA ) ) 
     &                  = EMFVAL( IDX, UNCIDX( S,PA ) )

                END IF

                REMPSTEP = SSAMP( UNCIDX( S,PA ) )

C.............  (I) Sources-different, hourly-same
            ELSE IF( APR .EQ. 2 ) THEN

                IF( PA .EQ. 1 ) THEN

                    CALL RANDOM_NUMBER( U )
                    IDX = INT( U*NUMEP( S,PA ) + 1 )
                    REMPSTEP = EMFVAL( IDX, UNCIDX( S,PA ) )

                ELSE

                    REMPSTEP = SAMPGEN( S,1 )

                END IF

C.............  (ST) Sources-same, hourly-different
            ELSE IF( APR .EQ. 3 ) THEN

                IF( LTIME .NE. ITIME ) STSAMP = AMISS3  ! array

                IF( STSAMP( UNCIDX( S,PA ) ) .EQ. AMISS3 ) THEN

                    CALL RANDOM_NUMBER( U )
                    IDX = INT( U*NUMEP( S,PA ) + 1 )
                    STSAMP( UNCIDX( S,PA ) ) 
     &                  = EMFVAL( IDX, UNCIDX( S,PA ) )

                    WRITE( MESG, 94010 ) 
     &                     'IDX= ', IDX,
     &                     'UNCIDX( S,PA ) = ', UNCIDX( S,PA ),
     &                     'LTIME =', LTIME,
     &                     'ITIME =',ITIME
                    CALL M3MSG2( MESG )

                END IF

                REMPSTEP = STSAMP( UNCIDX( S,PA ) )

C.............  (IT) Sources-same, hourly-different
            ELSE IF( APR .EQ. 4 ) THEN

                CALL RANDOM_NUMBER( U )
                IDX = INT( U*NUMEP( S,PA ) + 1 )
                REMPSTEP = EMFVAL( IDX, UNCIDX( S,PA ) )

            ELSE

                WRITE( MESG, 94010 ) 
     &                 'Invalid empirical stepwise approach for ' //
     &                 'source ', S, 'and pollutant/activity ', PA
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

            END IF

            LTIME = ITIME


C.........  Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

            END FUNCTION REMPSTEP

C**********************************************************************

            REAL FUNCTION RWEIBULL( ALPHA, BETA )

C.............  function parameters
            REAL,             INTENT(IN) :: ALPHA  
            REAL,             INTENT(IN) :: BETA   

C.............  Local variables
            REAL                            U      ! uniform (0,1) RNs

C**********************************************************************

C BUG: needs to be tested

            CALL RANDOM_NUMBER( U )

            RWEIBULL = ALPHA*( -LOG( U ) ) ** ( 1.0 / BETA )

            RETURN

            END FUNCTION RWEIBULL
C**********************************************************************

            REAL FUNCTION RBETA( ALPHA, BETA )

C.............  function parameters
            REAL,             INTENT(IN) :: ALPHA  
            REAL,             INTENT(IN) :: BETA   

C.............  Local variables
            REAL                            Y1, Y2 ! uniform (0,1) RNs

C**********************************************************************

C BUG: needs to be tested

            Y1 = RGAMMA( ALPHA, 1.0 )
            Y2 = RGAMMA( BETA, 1.0 )
            RBETA = Y1 / ( Y1 + Y2 )

            RETURN

            END FUNCTION RBETA

C**********************************************************************

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
            LOGICAL           NOT_FOUND            ! flag to signal RN found

            REAL              AA, BB, D, P, Q, V, W, Y, Z ! temporay reals

            REAL           :: NEARZERO = 0.001     ! near 0 constant
            REAL           :: LAMBDA = 0.0         ! model shift value
            REAL              U1, U2               ! uniform (0,1) RNs

C**********************************************************************

            IF( ALPHA .GT. 1.0 - NEARZERO .AND.
     &          ALPHA .LT. 1.0 + NEARZERO      ) THEN 

                RGAMMA = REXPONENTIAL( BETA )

            ELSE IF ( ALPHA .GT. 1.0 ) THEN

                NOT_FOUND = .TRUE.
                AA = 1/SQRT( 2*ALPHA - 1.0 )
                BB = ALPHA - LOG( 4.0 )
                Q  = ALPHA + 1.0/AA
                D  = 1.0 + LOG( 4.5 )
                DO WHILE( NOT_FOUND )

                    CALL RANDOM_NUMBER( U1 )
                    CALL RANDOM_NUMBER( U2 )
                    V = AA*LOG( U1/( 1.0 - U1 ) ) 
                    Y = ALPHA*EXP( V )
                    Z = U1*U1*U2
                    W = BB + Q*V - Y
                    P = LOG( Z )
                    IF( ( W + D - 4.5*Z) .GE. 0.0 ) THEN

                        RGAMMA = Y*BETA + LAMBDA
                        NOT_FOUND = .FALSE.

                    ELSE IF( W .GE. P ) THEN

                        RGAMMA = Y*BETA + LAMBDA
                        NOT_FOUND = .FALSE.

                    END IF
                END DO

            ELSE IF ( ALPHA .LT. 1.0 ) THEN

                NOT_FOUND = .TRUE.
                BB = ( EXP( 1.0 ) + ALPHA ) / EXP( 1.0 )
                DO WHILE( NOT_FOUND )
                    
                    CALL RANDOM_NUMBER( U1 )
                    CALL RANDOM_NUMBER( U2 )
                    IF( (BB*U1) .GT. 1.0 ) THEN

                        Y = -LOG(( BB - BB*U1) / ALPHA )
                        IF( U2 .LE. Y**(ALPHA - 1.0) ) THEN
                            RGAMMA = Y*BETA + LAMBDA
                            NOT_FOUND = .FALSE.
                        END IF

                    ELSE

                        Y = ( BB*U1 ) ** ( 1.0 / ALPHA )
                        IF( U2 .LE. EXP( Y ) ) THEN
                            RGAMMA = Y*BETA + LAMBDA
                            NOT_FOUND = .FALSE.
                        END IF

                    END IF
                END DO

            END IF

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

        END SUBROUTINE GETRNUMS
