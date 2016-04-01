
        SUBROUTINE HRBEIS( JDATE, JTIME, NX, NY, MSPCS, PX_VERSION,
     &                     INITIAL_HOUR, COSZEN, SEMIS,
     &                     GROWAGNO, NGROWAGNO, NONAGNO, SLAI,
     &                     TA, SOILM, SOILT, ISLTYP, RAIN, TSOLAR,
     &                     PRES, PTYPE, PULSEDATE, PULSETIME, EMPOL )

C***********************************************************************
C  subroutine body starts at line  143
C
C  DESCRIPTION:
C
C     Uses PAR and sfc temperature data to calculate
C     biogenic ISOP and MBO emissions.  Other emissions are
C     calculated using the temperature data only.
C
C  PRECONDITIONS REQUIRED:
C     PAR and Surface Temperature
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C    4/01 : Prototype by JMV
C    6/05 : updates for BEIS3.3 by D. Schwede (BEIS3.13)
C    8/05 : additional diagnostic messages for PAR out of bounds (G. Pouliot)
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
C***********************************************************************

        IMPLICIT NONE

C.........  INCLUDES
        INCLUDE 'B3V14DIMS3.EXT'  ! biogenic-related constants

C.........  ARGUMENTS and their descriptions
        INTEGER, INTENT (IN)  :: JDATE   !  current simulation date (YYYYDDD)
        INTEGER, INTENT (IN)  :: JTIME   !  current simulation time (HHMMSS)
        INTEGER, INTENT (IN)  :: NX      !  no. columns
        INTEGER, INTENT (IN)  :: NY      !  no. rows
        INTEGER, INTENT (IN)  :: MSPCS   !  no. of output species

        LOGICAL, INTENT (IN)  :: PX_VERSION    ! true: using PX version of MCIP
        LOGICAL, INTENT (IN)  :: INITIAL_HOUR  ! true:

        REAL, INTENT (IN)  ::  COSZEN   ( NX, NY )    !  cosine of zenith angle
        REAL, INTENT (IN)  ::  SEMIS    ( NX, NY, NSEF ) ! norm emissions
        REAL, INTENT (IN)  ::  GROWAGNO ( NX, NY )    ! growing season NO emissions
        REAL, INTENT (IN)  ::  NGROWAGNO( NX, NY )    ! non-growing season NO emissions
        REAL, INTENT (IN)  ::  NONAGNO  ( NX, NY )    ! non-agriculuture NO emissions
        REAL, INTENT (IN)  ::  SLAI  ( NX, NY, NLAI ) ! leaf area indices
        REAL, INTENT (IN)  ::  TA    ( NX, NY )       ! air temperature (K)
        REAL, INTENT (IN)  ::  SOILM ( NX, NY )       ! soil moisture
        REAL, INTENT (IN)  ::  SOILT ( NX, NY )       ! soil temperature
        REAL, INTENT (IN)  ::  ISLTYP( NX, NY )       ! soil type
        REAL, INTENT (IN)  ::  RAIN  ( NX, NY)        ! rainfall rate (cm/ 24hr)
        REAL, INTENT (IN)  ::  TSOLAR( NX, NY )       ! PAR
        REAL, INTENT (IN)  ::  PRES  ( NX, NY )       ! surface pressure (mb)

        INTEGER, INTENT (IN OUT) :: PTYPE(NX, NY)      ! 'pulse' type
        INTEGER, INTENT (IN OUT) :: PULSEDATE (NX, NY) ! date of pulse start
        INTEGER, INTENT (IN OUT) :: PULSETIME (NX, NY) ! date of pulse end

        REAL, INTENT (OUT)  ::  EMPOL( NX, NY, NSEF ) !  output pol emissions

C.........  SCRATCH LOCAL VARIABLES and their descriptions
        INTEGER         R, C, L, I      !  counters
        INTEGER         IAFTER

        REAL            CFOTHR       !  isop corr fac -- non-forest
        REAL            CFCLAI       !  ISOP CORR FAC -- LAI
        REAL            CFNO         !  NO correction factor
        REAL            CFOVOC       !  non-isop corr fac
	REAL            CFSESQT      !  sesquiterpene corr fac
        REAL            PAR          !  photo. actinic flux (UE/M**2-S)
        REAL            CT, DT       !  temperature correction
        REAL            TAIR         !  surface temperature
        REAL            RK           !  k from Geron and Guenther
        REAL            CSUBL        !  C sub l
        REAL            TLAI         !  temporary storage of LAI
        REAL            SOLTMP       !  temporary storage of radiation
        REAL            PSFC         !  temporary storage of sfc pressure (mb)
        REAL            ZEN          !  zenith angle
        REAL            PARDB        !  par direct beam
        REAL            PARDIF       !  par diffuse

        CHARACTER(5)    BTMP         !  temporary variable name
        CHARACTER(256)  MESG         !  message buffer

        CHARACTER(16) :: PROGNAME = 'HRBEIS'   !  program name

C***********************************************************************
C   begin body of subroutine HRBEIS

C.........  Loop through cells
        DO R = 1, NY
            DO C = 1, NX

                TAIR = TA( C, R )         ! unit in degree K

C..................  Check max and min bounds for temperature
C                    Note we no longer cap temperature for isoprene
                IF (TAIR .LT. 200.0) THEN
                    WRITE( MESG, 94010 ) 'TAIR=', TAIR,
     &                  'out of range at (C,R)=', C, R
                    CALL M3EXIT( PROGNAME, JDATE, JTIME, MESG, 2 )
                END IF

                 IF (TAIR .GT. 315.0 ) THEN
                     WRITE( MESG, 94020 ) 'TAIR=', TAIR,
     &                  'out of range at (C,R)=', C, R,
     &                  ' resetting to 315K'
                     CALL M3WARN( PROGNAME, JDATE, JTIME, MESG )
                     TAIR = 315.0
                 END IF

C.................  Calculate temperature correction term
                DT = 28668.514 / TAIR
                CT = EXP( 37.711 - 0.398570815 * DT ) /
     &                  (1.0 + EXP( 91.301 - DT ) )

                SOLTMP = TSOLAR( C, R )

C.................  Cosine of zenith angle to zenith angle (radians)
                ZEN =  ACOS( COSZEN( C, R ) )
                PSFC = PRES( C, R )

                CALL GETPAR( SOLTMP, PSFC, ZEN, PARDB, PARDIF )

                PAR = PARDB + PARDIF

C.................  Check max/min bounds of PAR and calculate
C                   biogenic ISOP
                IF ( PAR .LT. 0.00 .OR. PAR .GT. 2600.0 ) THEN
		
                    WRITE( MESG, 94030 ) 'PAR=', PAR,
     &                  'out of range at (C,R)=', C, R,
     &                  'PARDB  = ', PARDB,
     &                  'PARDIF = ', PARDIF,
     &                  'SOLTMP = ', SOLTMP,
     &                  'PSFC   = ', PSFC,
     &                  'ZEN    = ', ZEN
    
                     CALL M3MSG2(MESG)
!                    CALL M3EXIT( PROGNAME, JDATE, JTIME, MESG, 2 )
                END IF

C.................  Compute ISOP and MBO and METH emissions first
C                   Note assumption that these are the first 3
C                   species in LAITYPE and BIOTYPE arrays
                DO I = 1, NLAI

                    BTMP = LAITYPES( I )
                    TLAI = SLAI( C, R, I )

C.....................  Adjust methanol based on T. Pierce recommendation (1-16-03)
                    IF( TRIM( BTMP ) == 'METH' ) THEN
                        TLAI = MAX( 3.0, TLAI )
                    END IF

                    IF ( TLAI .GT. 10.0 ) THEN
                        WRITE( MESG, 94010 ) 'LAI=', TLAI,
     &                  'out of range at (C,R)=', C, R
                        CALL M3EXIT( PROGNAME, JDATE, JTIME, MESG, 2 )
                    END IF

C.....................  Initialize csubl
                    CSUBL = 0.0

                    IF ( PAR .LE. 0.01 .OR.
     &                   COSZEN( C,R ) .LE. 0.02079483 ) THEN

                        EMPOL( C, R, I ) = 0.0

                    ELSE
                        IF ( TLAI .GT. 0.1 ) THEN
                           CSUBL = CLNEW( ZEN, PARDB, PARDIF, TLAI )
                        ELSE  ! keep this or not?
                            CSUBL  = CGUEN( PAR )
                        END IF

                        EMPOL( C, R, I ) = SEMIS( C,R, I ) * CT * CSUBL
                    END IF

                END DO ! end ISOP and MBO calculations loop

C.................  Calculate other biogenic emissions except NO
C                   Note not speciated here
C                   Limit temerature to 315 K for monoterpenes and other VOCs
               IF (TAIR .GT. 315.0 ) THEN
                    WRITE( MESG, 94020 ) 'TAIR=', TAIR,
     &               'out of range at (C,R)=', C, R,
     &               ' resetting to 315K for monoterpene and other VOCs'
                    CALL M3WARN( PROGNAME, JDATE, JTIME, MESG )
                    TAIR = 315.0
                END IF

                CFOVOC  = EXP( 0.09 * ( TAIR - 303.0 ) )
                CFSESQT = EXP( 0.17 * ( TAIR - 303.0 ) )
		
                IAFTER = NLAI + 1
                DO I = IAFTER, NSEF - 2
                    EMPOL( C,R,I ) = SEMIS( C,R,I ) * CFOVOC
                END DO
		DO I = NSEF, NSEF  ! Sesquiterpene emissions
		    EMPOL( C,R,I ) = SEMIS( C,R,I ) * CFSESQT
		END DO

            END DO ! end loop over columns
        END DO ! end loop over rows

C.........  Calculate NO emissions
        CALL HRNO( JDATE, JTIME, NX, NY,  TA, SOILM, SOILT,
     &             ISLTYP, RAIN, GROWAGNO, NGROWAGNO, NONAGNO,
     &             PX_VERSION, INITIAL_HOUR, PTYPE, PULSEDATE,
     &             PULSETIME, EMPOL )

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( A, F10.2, 1X, A, I3, ',', I3 )
94020   FORMAT( A, F10.2, 1X, A, I3, ',', I3, A )
94030   FORMAT( A, F10.2, 1X, A, I3, ',', I3, 5(A, F10.2))

C***************** CONTAINS ********************************************

        CONTAINS

C.............  Function to calculate csubl based on zenith angle, par, and lai
            REAL FUNCTION CLNEW( ZEN, PARDB, PARDIF, TLAI )
C******** Reference:CN98
C      Campbell, G.S. and J.M. Norman. 1998. An Introduction to Environmental Biophysics, 
C      Springer-Verlag, New York.
C
C
            IMPLICIT NONE

C.............  Function arguments
            REAL, INTENT (IN) :: PARDB    ! direct beam PAR( umol/m2-s)
            REAL, INTENT (IN) :: PARDIF   ! diffuse PAR ( umol/m2-s)
            REAL, INTENT (IN) :: ZEN      ! solar zenith angle (radians)
            REAL, INTENT (IN) :: TLAI     ! leaf area index for grid cell

C.............  Local variables
            REAL ALPHA              ! leave absorptivity
            REAL KBE                ! extinction coefficient for direct beam
            REAL KD                 ! extinction coefficient for diffuse radiation
            REAL CANPARSCAT         ! exponentially wtd scattered PAR (umol/m2-s)
            REAL CANPARDIF          ! exponentially wtd diffuse PAR (umol/m2-s)
            REAL PARSHADE           ! PAR on shaded leaves (umol/m2-s)
            REAL PARSUN             ! PAR on sunlit leaves (umol/m2-s)
            REAL LAISUN             ! LAI that is sunlit
            REAL FRACSUN            ! fraction of leaves that are sunlit
            REAL FRACSHADE          ! fraction of leaves that are shaded
            REAL SQALPHA            ! square root of alpha

C-----------------------------------------------------------------------------
            ALPHA = 0.8
            SQALPHA = SQRT(0.8)
            KD = 0.68

C.............  CN98 - eqn 15.4, assume x=1
            KBE = 0.5 * SQRT(1. + TAN( ZEN ) * TAN( ZEN ))

C.............  CN98 - p. 261 (this is usually small)
            CANPARSCAT = 0.5 * PARDB * (EXP(-1.* SQALPHA * KBE * TLAI) -
     &                   EXP(-1.* KBE * TLAI))

C.............  CN98 - p. 261 (assume exponentially wtd avg)
            CANPARDIF  = PARDIF * (1. - EXP(-1. * SQALPHA * KD * TLAI))
     &                 /(SQALPHA * KD * TLAI)

C.............  CN98 - p. 261 (for next 3 eqns)
C               note that we use the incoming (not absorbed) PAR
            PARSHADE   = CANPARDIF + CANPARSCAT
            PARSUN     = KBE * PARDB + PARSHADE
	      LAISUN     = (1. - EXP(-1. * KBE * TLAI))/KBE
	      FRACSUN    = LAISUN/TLAI
	      FRACSHADE  = 1. - FRACSUN

C...........  cguen is guenther's eqn for computing light correction as a
C             function of PAR...fracSun should probably be higher since
C             sunlit leaves tend to be thicker than shaded leaves.  But
C             since we need to make crude asmptns regarding leave
C             orientation (x=1), will not attempt to fix at the moment.

            CLNEW = FRACSUN * CGUEN( PARSUN ) +
     &              FRACSHADE * CGUEN( PARSHADE )

            RETURN

            END FUNCTION CLNEW

C-----------------------------------------------------------------------------

C.............  Function to calculate Guenther's equation for computing
C               light correction
C    Reference:   Guenther, A., B. Baugh, G. Brasseur, J. Greenberg, P. Harley, L. Klinger,
C   D. Serca, and L. Vierling, 1999: Isoprene emission estimates and uncertainties 
C   for the Central African EXPRESSO Study domain. J. Geophys. Res., 104, 30625-30639.
C
C
C
            REAL FUNCTION CGUEN( PARTMP )

            IMPLICIT NONE

C.............  Function arguments
            REAL, INTENT (IN) :: PARTMP
            REAL, PARAMETER :: ALPHA = 0.001
            REAL, PARAMETER :: CL = 1.42

C-----------------------------------------------------------------------------

            IF ( PARTMP .LE. 0.01) THEN
                CGUEN = 0.0
            ELSE
                CGUEN = (ALPHA * CL * PARTMP) /
     &                  SQRT(1. + ALPHA * ALPHA * PARTMP * PARTMP)
            END IF

            RETURN

            END FUNCTION CGUEN

C-----------------------------------------------------------------------------

        END SUBROUTINE HRBEIS

