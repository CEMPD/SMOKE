
        SUBROUTINE HRBEIS4( JDATE, JTIME, NX, NY, MSPCS, 
     &                PX_VERSION,INITIAL_HOUR, COSZEN, SEMIS,
     &                GROWAGNO, NGROWAGNO, NONAGNO, TA,
     &                SOILM, SOILT, WSAT, WRF_WSAT, ISLTYP, RAIN,PRES,
     &                RSOLAR,LAI, SLAI, WRF_LAI, USTAR,RSTOM,RATM,
     &                Q2,TEMPG,PTYPE, PULSEDATE, PULSETIME, EMPOL )
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
C     HRNO_BEIS4, GETPARB
C
C  REVISION  HISTORY:
C    4/22 : Prototype by JMV based on CMAQ inline logic with updates
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
        LOGICAL, INTENT( IN ) :: WRF_LAI
        LOGICAL, INTENT( IN ) :: WRF_WSAT

        REAL, INTENT (IN)  ::  COSZEN   ( NX, NY )    !  cosine of zenith angle
        REAL, INTENT (IN)  ::  SEMIS    ( NX, NY, NSEF ) ! norm emissions
        REAL, INTENT (IN)  ::  GROWAGNO ( NX, NY )    ! growing season NO emissions
        REAL, INTENT (IN)  ::  NGROWAGNO( NX, NY )    ! non-growing season NO emissions
        REAL, INTENT (IN)  ::  NONAGNO  ( NX, NY )    ! non-agriculuture NO emissions
        REAL, INTENT (IN)  ::  TA    ( NX, NY )       ! air temperature (K)
        REAL, INTENT (IN)  ::  SOILM ( NX, NY )       ! soil moisture
        REAL, INTENT (IN)  ::  SOILT ( NX, NY )       ! soil temperature
        REAL, INTENT (IN)  ::  WSAT ( NX, NY )       ! WSAT 

        REAL, INTENT (IN)  ::  ISLTYP( NX, NY )       ! soil type
        REAL, INTENT (IN)  ::  RAIN  ( NX, NY)   ! rainfall rate (cm/ 24hr)

        REAL,    INTENT( IN ) :: RSOLAR( NX,NY )   ! surface radiation [w/m**2]
        REAL,    INTENT( IN ) :: PRES  ( NX,NY )  ! surface pressure [Pa]  (new)
        REAL,    INTENT (IN ) :: SLAI  ( NX, NY, NLAI ) ! leaf area indices
        REAL,    INTENT( IN ) :: LAI   ( NX,NY )   ! WRF LAI [m**2/m**2]
        REAL,    INTENT( IN ) :: USTAR ( NX,NY )   ! sfc friciton velocity [m/s]
        REAL,    INTENT( IN ) :: RSTOM ( NX,NY )  ! Stomatal resistance [s/m]
        REAL,    INTENT( IN ) :: RATM  ( NX,NY )  ! Aerodynamic resistance [s/m]
        REAL,    INTENT( IN ) :: Q2    ( NX,NY )  ! 2-m water vapor mix ratio [kg/kg]
        REAL,    INTENT( IN ) :: TEMPG ( NX,NY ) ! Ground temperature [K]
      
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


        REAL            CSUBL        !  C sub l

        REAL            SOLTMP       !  temporary storage of radiation




!!!! NEW FOR VERSION 3.50 BEGIN
      REAL           CT_SUN       ! temperature correction
      REAL           DT_SUN       ! temperature correction
      REAL           CT_SHADE     ! temperature correction       
      REAL           DT_SHADE     ! temperature correction      
      REAL           TAIR         ! local 2 meter temperature
      REAL           DTLEAF_SUN   ! Difference between mean canopy leaf and ambient temperature [K]
      REAL           DTLEAF_SHADE ! Difference between mean canopy leaf and ambient temperature [K]      
      REAL           TLEAF_SUN    ! Mean canopy leaf temperature [K]
      REAL           TLEAF_SHADE  ! Mean canopy leaf temperature [K]
      REAL           RBW          ! Quasi-laminar boundary layer resistance for water vapor [s/m]
      REAL           RBH          ! Quasi-laminar boundary layer resistance for heat [s/m]     
      REAL           RH           ! Relative humidity [ratio 0-1]
      REAL           ES           ! Saturation vapor pressure for 2 meter T  [Pa]          
      REAL           SHF          ! Soil heat flux [W/m**2]    
      REAL           DVAP         ! vapor pressure deficit [Pa/Pa]  
      REAL           SSVP         ! Slope of the saturation vapor pressure curve over P [1/K] 
      REAL           GVAP         ! canopy water vapor conuctance m/s
      REAL           GHT          ! canopy heat conductance m/s      
      REAL           CPAIR        ! specific heat of air
      REAL           LHV          ! Latent heat of vaporization       
      REAL           CPOT         ! potential temperature conversion 
      REAL           DENS         ! Dry air density kg/m**3
      REAL           LHSH_DIV     ! W/m**2 to K units conversion 
      REAL           LHSH_COMP    ! latent/sensible heat flux component of leaf energy bal
      REAL           RK           ! k from Geron and Guenther
      REAL           CSUBL_SUN    ! C sub l
      REAL           CSUBL_SHADE  ! C sub l
      REAL           FRACSUN      ! Fraction sun
      REAL           FRACSHADE    ! Fraction shade            
      REAL           TLAI         ! local LAI
      REAL           SOLRAD       ! local solar radiation [W/m**2]
      REAL           PSFC         ! local sfc pressure (mb)
      REAL           ZEN          ! zenith angle
      REAL           PARDB        ! PAR direct beam
      REAL           PARDIF       ! PAR diffuse
      REAL           COSZ         ! local cosine of zenith angle

      REAL, PARAMETER :: SVP2       = 17.67   ! from MM5 and WRF PX
      REAL, PARAMETER :: SVP3       = 29.65   ! from MM5 and WRF PX
      REAL, PARAMETER :: CV         = 8.0e-6  ! Resistance to soil heat conductance under vegetation [s/m]
      REAL, PARAMETER :: RRAD       = 230.0   ! Atmospheric radiative resistance Monteith 1973 [s/m]
      REAL, PARAMETER :: PR         = 0.709   ! prandtl number [dim'less]
      REAL, PARAMETER :: KVIS       = 0.132   ! kinimatic viscosity of air
      REAL, PARAMETER :: DWAT       = 0.2178  ! Diffusivity of water vapor in air
      REAL, PARAMETER :: SCW        = KVIS/DWAT ! schmidt number for water vapor
      REAL, PARAMETER :: TWOTHIRDS  = 2.0 / 3.0
      REAL, PARAMETER :: REFLDV     = 0.057 ! visible light reflection coefficient from MEGAN 2.10 


!   The following constants come from CMAQ code !!! (should this be an include or not?
! vapor press of water at 0 C [ Pa ] Source: CRC76 pp. 6-15
      REAL, PARAMETER :: VP0 = 611.29

      REAL,      PARAMETER :: PI = 3.14159265
! length of a sidereal day [ sec ]
! FSB: Source: CRC76 pp. 14-6 
      REAL, PARAMETER :: SIDAY = 86164.09

! Standard Temperature [ K ]
      REAL, PARAMETER :: STDTEMP = 273.15
! latent heat of vaporization of water at 0 C [ J/kg ]
      REAL, PARAMETER :: LV0 = 2.501E6
! universal gas constant [ J/mol-K ]
      REAL, PARAMETER :: RGASUNIV = 8.314510

! mean molecular weight for dry air [ g/mol ]
! FSB: 78.06% N2, 21% O2, and 0.943% A on a mole 
! fraction basis ( Source : Hobbs, 1995) pp. 69-70
      REAL, PARAMETER :: MWAIR = 28.9628

! dry-air gas constant [ J / kg-K ]
      REAL, PARAMETER :: RDGAS = 1.0E3 * RGASUNIV / MWAIR   ! 287.07548994
! specific heat of dry air at constant pressure [ J/kg-K ]
      REAL, PARAMETER :: CPD = 7.0 * RDGAS / 2.0            ! 1004.7642148 

! standard atmosphere  [ Pa ]
      REAL, PARAMETER :: STDATMPA = 101325.0

! end standard CMAQ constants

      
      
        CHARACTER(5)    BTMP         !  temporary variable name
        CHARACTER(256)  MESG         !  message buffer

        CHARACTER(16) :: PROGNAME = 'HRBEIS4'   !  program name

C***********************************************************************
C   begin body of subroutine HRBEIS

C.........  Loop through cells
        DO R = 1, NY
          DO C = 1, NX

            TAIR = TA( C, R )         ! unit in degree K

!CCCC BEBIN UPDATED CODE FOR REVISED TEMP LEAF CALCULATION

            COSZ = COSZEN( C,R ) 

C..................  Check max and min bounds for temperature
C                    Note we no longer cap temperature for isoprene
            IF (TAIR .LT. 200.0) THEN
                 WRITE( MESG, 94010 ) 'TAIR=', TAIR,
     &                'out of range at (C,R)=', C, R
                    CALL M3EXIT( PROGNAME, JDATE, JTIME, MESG, 2 )
            END IF

            SOLRAD = RSOLAR( C,R )
C Cosine of zenith angle to zenith angle (radians)
            ZEN =  ACOS( COSZ )
            PSFC = PRES( C,R )

C atmospheric water vapor variables used for leaf latent heat flux
            IF ( TAIR .LE. STDTEMP ) THEN
               ES = VP0 * EXP( 22.514 - (6.15e3 / TAIR) )
            ELSE
               ES = VP0 * EXP( SVP2 * (TAIR - STDTEMP) /
     &                                (TAIR - SVP3) )
            END IF
            RH   = Q2( C,R )*(PSFC-ES)/(0.622*ES)
            RH   = MIN(1.0,RH)
            DVAP = ES*(1.0-RH)/PSFC
            SSVP = (SVP2*(STDTEMP-SVP3)*ES/(TAIR-SVP3)**2)/
     &              PSFC

C calculate the soil heaf flux under a canopy folowing WRF3.4.1 PX
            SHF = -2.0*PI/SIDAY*(TEMPG( C,R ) - TAIR)/CV

C calculate the heat and water vapor quasilaminar boundary layer resistance
            RBH = 5.0/USTAR( C,R )
            RBW = RBH*(SCW/PR)**TWOTHIRDS
  
C calculate the specific heat of air and latent heat of vaporization
            CPAIR = CPD * (1.0 + 0.84 * Q2( C,R ))
            LHV   = LV0 - 2370.0 * (TAIR - STDTEMP)

C calculate the leaf  water vapor and heat conductance 
            GHT  = 1 / (RATM( C,R ) + RBH  ) + 1/RRAD
            GVAP = 1 / (RATM( C,R ) + RBW + RSTOM( C,R ) * LAI( C,R ))

C Calculate the potential temperature conversion from the sensible heat flux in WRF 3.4.1
            CPOT = (STDATMPA/PSFC)**(RDGAS/CPAIR)

C Calculate the solar radiation and soil heat flux of the leaf energy budget
            DENS     =  PSFC /( RDGAS * TAIR )
            LHSH_DIV =  DENS * CPOT * CPAIR 
     &                  * (GHT + 1 / (RATM( C,R ) + RBH  )) + 
     &                  DENS * LHV * SSVP * GVAP

C calculate the latent heat flux portion of the leaf energy budget 
            LHSH_COMP = SHF - LHV * DENS * GVAP * (ES-RH*ES)/ PSFC

C Direct and diffuse photosynthetically active radiation
            CALL GETPARB( SOLRAD, PSFC, COSZ, PARDB, PARDIF )

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
               IF ( WRF_LAI ) THEN
                TLAI = LAI( C, R )
               ELSE
                TLAI = SLAI( C, R, I )
               ENDIF
C.....................  Adjust methanol based on T. Pierce recommendation (1-16-03)
               IF( TRIM( BTMP ) == 'METH' ) THEN
                 TLAI = MAX( 3.0, TLAI )
               END IF

               IF ( TLAI .GT. 10.0 ) THEN
                 WRITE( MESG, 94010 ) 'LAI=', TLAI,
     &                  'out of range at (C,R)=', C, R
                 CALL M3EXIT( PROGNAME, JDATE, JTIME, MESG, 2 )
               END IF

C Initialize csubl
               CSUBL_SUN   = 0.0
               CSUBL_SHADE = 0.0       

               IF ( PARDB + PARDIF .EQ. 0.0 ) THEN
                 EMPOL( C,R,I ) = 0.0
               ELSE
                 CALL CLNEW_SUB( ZEN, PARDB, PARDIF, TLAI, LHSH_DIV,
     &                           LHSH_COMP, DTLEAF_SUN, DTLEAF_SHADE, 
     &                           CSUBL_SUN, CSUBL_SHADE, FRACSUN, 
     &                           FRACSHADE, SOLRAD, REFLDV )

                 TLEAF_SUN   = DTLEAF_SUN   + TAIR
                 TLEAF_SHADE = DTLEAF_SHADE + TAIR
C Calculate temperature correction term
                 DT_SUN   = 28668.514 / TLEAF_SUN
                 DT_SHADE = 28668.514 / TLEAF_SHADE 
                 CT_SUN   = EXP( 37.711 - 0.398570815 * DT_SUN ) /
     &                         ( 1.0 + EXP( 91.301 - DT_SUN ) )
                 CT_SHADE = EXP( 37.711 - 0.398570815 * DT_SHADE ) /
     &                         ( 1.0 + EXP( 91.301 - DT_SHADE ) )     
                 EMPOL( C,R,I ) = SEMIS( C,R,I )*( FRACSUN * CT_SUN * 
     &             CSUBL_SUN + FRACSHADE * CT_SHADE * CSUBL_SHADE )

               END IF

            END DO ! end ISOP and MBO calculations loop

C Only estiamte BVOC emissions for vegatation
C Estimate emissions for OVOC and SESQT 
 
            IF( TLAI .GT. 0.0 ) THEN

              CALL CLNEW_SUB( ZEN, PARDB, PARDIF, LAI(C,R), LHSH_DIV,
     &                      LHSH_COMP, DTLEAF_SUN, DTLEAF_SHADE, 
     &                      CSUBL_SUN, CSUBL_SHADE, FRACSUN, FRACSHADE,
     &                      SOLRAD, REFLDV )
     
              TLEAF_SUN   = TAIR + DTLEAF_SUN
              TLEAF_SHADE = TAIR + DTLEAF_SHADE
C Calculate other biogenic emissions except NO
C Note not speciated here
C Limit temerature to 315 K for monoterpenes and other VOCs
              TLEAF_SUN   = MIN( TLEAF_SUN, 315.0 )
              TLEAF_SHADE = MIN( TLEAF_SHADE, 315.0 )    

              CFOVOC  = FRACSUN * EXP( 0.09 * ( TLEAF_SUN   - 303.0 ) )+
     &                FRACSHADE * EXP( 0.09 * ( TLEAF_SHADE - 303.0 ) ) 
              CFSESQT = FRACSUN * EXP( 0.17 * ( TLEAF_SUN   - 303.0 ) )+
     &                FRACSHADE * EXP( 0.17 * ( TLEAF_SHADE - 303.0 ) ) 
            ELSE
C If LAI = 0 zero out emission factors
              CFOVOC  = 0.0
              CFSESQT = 0.0
            END IF

            DO I = NLAI + 1, NSEF - 2
               EMPOL( C,R,I ) = SEMIS( C,R,I ) * CFOVOC
            END DO
C
            DO I = NSEF, NSEF  ! Sesquiterpene emissions
               EMPOL( C,R,I ) = SEMIS( C,R,I ) * CFSESQT
            END DO

          END DO ! end loop over columns
        END DO ! end loop over rows

C.........  Calculate NO emissions
        CALL HRNO_BEIS4( JDATE, JTIME, NX, NY,  TA, SOILM, SOILT,
     &            WSAT, WRF_WSAT, ISLTYP, RAIN, GROWAGNO, NGROWAGNO, 
     &            NONAGNO, PX_VERSION, INITIAL_HOUR, PTYPE, PULSEDATE,
     &            PULSETIME, EMPOL )


        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( A, F10.2, 1X, A, I3, ',', I3 )
94020   FORMAT( A, F10.2, 1X, A, I3, ',', I3, A )
94030   FORMAT( A, F10.2, 1X, A, I3, ',', I3, 5(A, F10.2))

C***************** CONTAINS ********************************************

        CONTAINS

C Function to calculate csubl based on zenith angle, par, and lai
         SUBROUTINE CLNEW_SUB( ZEN, PARDB, PARDIF, TLAI, LHSH_DIV,
     &                           LHSH_COMP, DTLSUN, DTLSHADE, 
     &                           CSUBL_SUN, CSUBL_SHADE, FRACSUN, 
     &                           FRACSHADE, SOLRAD, REFLDV )

         IMPLICIT NONE

C Function arguments:
         REAL, INTENT( IN )  :: PARDB    ! direct beam PAR( umol/m2-s)
         REAL, INTENT( IN )  :: PARDIF   ! diffuse PAR ( umol/m2-s)
         REAL, INTENT( IN )  :: ZEN      ! solar zenith angle (radians)
         REAL, INTENT( IN )  :: TLAI     ! leaf area index for grid cell
         REAL, INTENT( IN )  :: LHSH_DIV  
         REAL, INTENT( IN )  :: LHSH_COMP
         REAL, INTENT( IN )  :: SOLRAD
         REAL, INTENT( IN )  :: REFLDV
         REAL, INTENT( OUT ) :: CSUBL_SUN
         REAL, INTENT( OUT ) :: CSUBL_SHADE
         REAL, INTENT( OUT ) :: DTLSUN           ! Sun leaf temperature [K]
         REAL, INTENT( OUT ) :: DTLSHADE         ! Sun leaf temperature [K]
         REAL, INTENT( OUT ) :: FRACSUN          ! fraction of leaves that are sunlit
         REAL, INTENT( OUT ) :: FRACSHADE        ! fraction of leaves that are shaded
C Parameters:
         REAL, PARAMETER :: ALPHA = 0.8 ! leaf absorptivity
         REAL, PARAMETER :: KD = 0.68   ! extinction coefficient for diffuse radiation
C Local variables:
         REAL, SAVE :: SQALPHA ! square root of alpha
         REAL KBE              ! extinction coefficient for direct beam
         REAL CANPARSCAT       ! exponentially wtd scattered PAR (umol/m2-s)
         REAL CANPARDIF_SUN    ! exponentially wtd diffuse PAR at the top of the canopy (umol/m2-s)
         REAL CANPARDIF_SHADE  ! exponentially wtd diffuse PAR in the shaded part of the canopy (umol/m2-s)
         REAL PARSHADE         ! PAR on shaded leaves (umol/m2-s)
         REAL PARSUN           ! PAR on sunlit leaves (umol/m2-s)
         REAL SOLSUN           ! RS transmitted to sunlit leaves W/m**2
         REAL SOLSHADE         ! RS transmitted to shaded leaves W/m**2
         REAL LAISUN           ! LAI that is sunlit
         REAL LAISHADE         ! LAI that is shaded


         LOGICAL, SAVE :: FIRSTIME = .TRUE.

C-----------------------------------------------------------------------
         IF ( FIRSTIME ) THEN
            FIRSTIME = .FALSE.
            SQALPHA = SQRT( ALPHA )
         END IF
C CN98 - eqn 15.4, assume x=1 (can use a table or atributes to change this)
C Set a ceiling for KBE to prevent a blow up at high zenith angles. This has
C little impact on the results because direct PAR is low under these conditions
         IF( ZEN .GE. 1.57 ) THEN
            KBE = 627.9
         ELSE
            KBE = 0.5 * SQRT( 1.0 + TAN( ZEN )**2 )
         END IF
         IF ( TLAI .GT. 0.1 ) THEN
            IF ( PARDB + PARDIF .GT. 0.0 ) THEN

C CN98 p-259 Sun and shaded areas of the canopy
               LAISUN     = ( 1.0 - EXP( -1.0 * KBE * TLAI ) ) / KBE
               LAISHADE   = MAX( TLAI - LAISUN, 0.0 )
               FRACSUN    = LAISUN / TLAI             
               FRACSHADE  = 1.0 - FRACSUN

C CN98 - p. 261 (this is usually small)
               CANPARSCAT = 0.5 * PARDB * ( EXP( -1.0 * SQALPHA * KBE *
     &           TLAI ) - EXP( -1.0 * KBE * TLAI ) )

C CN98 - p. 261 (assume exponentially wtd avg)
               CANPARDIF_SUN    = PARDIF * ( 1.0 - EXP( -1.0 * SQALPHA 
     &           * KD * LAISUN ) ) / ( SQALPHA * KD * LAISUN )

               CANPARDIF_SHADE  = CANPARDIF_SUN * ( EXP( -1.0 * SQALPHA
     &           * KD * LAISUN ) - EXP( -1.0 * SQALPHA * KD * TLAI ) ) 
     &             / ( SQALPHA * KD * (TLAI - LAISUN) )

C CN98 - p. 261 (for next 3 eqns)
C note that we use the incoming (not absorbed) PAR
               PARSHADE   = CANPARDIF_SHADE + CANPARSCAT
               PARSUN     = KBE * PARDB + CANPARDIF_SUN + CANPARSCAT
     
C calculate the leaf temperature following Campbel and Norman 1998 eq 14.6 
C with the addition of incomming atmospheric long wave irradiation resulting 
C in the cacelation of the long wave radiation budget
               SOLSUN    = SOLRAD * PARSUN / ( PARSUN + PARSHADE )
               SOLSHADE  = SOLRAD * PARSHADE / ( PARSUN + PARSHADE )
C Because PARSUN + PARSHADE != PARDB+PARDIF as the first considers a
C hemispheric sphere (crown canopy) and the later a horizontal plane (the
C ground). This changes the peak leaf temperature a bit, usually about 1
C hour earlier.
C               SOLSUN    = SOLRAD * PARSUN / ( PARDB + PARDIF )
C               SOLSHADE  = SOLRAD * PARSHADE / ( PARDB + PARDIF )
               DTLSUN    = ((1.0 - REFLDV) * SOLSUN + LHSH_COMP ) / 
     &                      LHSH_DIV
               DTLSHADE  = ((1.0 - REFLDV) * SOLSHADE + LHSH_COMP ) / 
     &                      LHSH_DIV
               DTLSUN    = MIN(DTLSUN,  10.0) 
               DTLSUN    = MAX(DTLSUN, -10.0) 
               DTLSHADE  = MIN(DTLSHADE,  10.0) 
               DTLSHADE  = MAX(DTLSHADE, -10.0) 

C cguen is Guenther's eqn for computing light correction as a function of
C PAR...fracSun should probably be higher since sunlit leaves tend to be
C thicker than shaded leaves. But since we need to make crude assumptions
C regarding leaf orientation (x=1), we will not attempt to fix at the moment.

               CSUBL_SUN   = CGUEN( PARDB + PARDIF, 0.0, LAISUN, KBE )
C By definition diffusive radiation, use the diffusive attenuation
C coefficient and Diffusive par at the bottom of the sunlit layer
               CSUBL_SHADE = CGUEN( CANPARDIF_SUN, LAISUN, TLAI, KD )
C               CSUBL_SHADE = CGUEN( PARDB + PARDIF, LAISUN, TLAI, KBE )
            
            ELSE ! to prevent divide by 0 when there is no solar rad
               CSUBL_SUN   = 0.0
               CSUBL_SHADE = 0.0
               FRACSUN     = 0.2
               FRACSHADE   = 0.8
               DTLSUN      = LHSH_COMP / LHSH_DIV
               DTLSHADE    = LHSH_COMP / LHSH_DIV
               DTLSUN      = MIN(DTLSUN,  10.0) 
               DTLSUN      = MAX(DTLSUN, -10.0) 
               DTLSHADE    = MIN(DTLSHADE,  10.0) 
               DTLSHADE    = MAX(DTLSHADE, -10.0) 
            END IF       
    
         ELSE 
            CSUBL_SUN   = CGUEN( PARDB + PARDIF, 0.0, TLAI, KBE )
            CSUBL_SHADE = 0.0
            FRACSUN     = 1.0
            FRACSHADE   = 0.0
            DTLSUN  = ((1.0 - REFLDV) * SOLRAD + LHSH_COMP ) / LHSH_DIV
            DTLSHADE    = 0.0
            DTLSUN      = MIN(DTLSUN,  10.0) 
            DTLSUN      = MAX(DTLSUN, -10.0)      
         END IF

         END SUBROUTINE CLNEW_SUB



         REAL FUNCTION CGUEN( PAR, LAI1, LAI2, KBE )

C 11/14 J. Bash - Updated to Niinemets et al. 2010 doi:10.1029/2010JG001436 
C                 Big leaf model which updates Guenther et al. 1993 doi:10.1029/93JD00527 for 
C                 in-canopy gradients

         IMPLICIT NONE

C Function arguments:
         REAL, INTENT( IN ) :: PAR
         REAL, INTENT( IN ) :: LAI1 ! top of the layer LAI
         REAL, INTENT( IN ) :: LAI2 ! bottom of the layer LAI
         REAL, INTENT( IN ) :: KBE  ! light extenction coefficient

C Parameters:
C         REAL, PARAMETER :: ALPHA = 0.0027 ! Guenther et al. 1993
C         REAL, PARAMETER :: CL    = 1.066  ! Guenther et al. 1993
C Consistent with Niinemets et al. 2010a
C Mean of reported Quercus rubra and Liquidambar styraciflua

         REAL, PARAMETER :: ALPHA = 0.0015 
         REAL, PARAMETER :: CL    = 1.2716 

C  Niinemets et al. 2010b to return 1 at standard conditions
C  standard conditions (PAR=1000, KBE = 0.68 )

C---------------------------------------------------------------------
C-----------------------------------------------------------------------
         IF ( PAR .LE. 0.01 ) THEN
            CGUEN = 0.0
         ELSE
C Niinemets et al. 2010 equation A9 integrated from LAI1 to LAI2
           CGUEN = CL * ( SQRT(1+ALPHA**2 * PAR**2 * EXP(-2*LAI1*KBE)) -
     &                  SQRT(1+ALPHA**2 * PAR**2 * EXP(-2*LAI2*KBE)) ) /
     &                  ( ALPHA * KBE * PAR )
         END IF

         RETURN

         END FUNCTION CGUEN

C-----------------------------------------------------------------------------
        END SUBROUTINE HRBEIS4
