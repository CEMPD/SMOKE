 
        SUBROUTINE  HRBEIS3( JDATE, JTIME, NX, NY, COSZEN,
     &                      SISOP, SMONO, SOVOC, SNO, SLAI,
     &                      TA, TSOLAR, PRES, EMPOL )

C***********************************************************************
C  subroutine body starts at line  143
C
C  DESCRIPTION:
C  
C     Uses PAR and sfc temperature data to calculate
C     biogenic ISOP emissions.  OVOC, MONO and NO emissions are
C     calculated using the temperature data only.  OVOC and MONO
C     are not speciated in this rountine. 
C
C  PRECONDITIONS REQUIRED:
C     PAR and Surface Temperature
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C    3/01 : Prototype by JMV
C 
C***********************************************************************
C
C Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
C                System
C File: @(#)$Id$
C
C COPYRIGHT (C) 2001, MCNC--North Carolina Supercomputing Center
C All Rights Reserved
C
C See file COPYRIGHT for conditions of use.
C
C MCNC-Environmental Programs Group
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

        INCLUDE 'PARMS3.EXT'      ! I/O API constants
        INCLUDE 'FDESC3.EXT'      ! I/O API file description data structure
        INCLUDE 'IODECL3.EXT'     ! I/O API function declarations
        INCLUDE 'EMCNST3.EXT'     ! Emissions constants
        INCLUDE 'CONST3.EXT'      ! More constants
        INCLUDE 'B3DIMS3.EXT'     ! biogenic-related constants

C...........   EXTERNAL FUNCTIONS and their descriptions:

        INTEGER         INDEX1
        EXTERNAL        INDEX1

C...........   ARGUMENTS and their descriptions:

        INTEGER, INTENT (IN)  :: JDATE   !  current simulation date (YYYYDDD)
        INTEGER, INTENT (IN)  :: JTIME   !  current simulation time (HHMMSS)
        INTEGER, INTENT (IN)  :: NX      !  no. columnse
        INTEGER, INTENT (IN)  :: NY      !  no. rows

        REAL, INTENT (IN)  ::  TA    ( NX, NY )    !  air temperature (K)
        REAL, INTENT (IN)  ::  TSOLAR( NX, NY )    !  PAR
        REAL, INTENT (IN)  ::  SISOP ( NX, NY )    ! norm ISOP emissions
        REAL, INTENT (IN)  ::  SMONO ( NX, NY )    !  norm MONO emissions
        REAL, INTENT (IN)  ::  SOVOC ( NX, NY )    !  norm OVOC emissions
        REAL, INTENT (IN)  ::  SNO   ( NX, NY )    !  nor NO emissions
        REAL, INTENT (IN)  ::  SLAI  ( NX, NY )    !  leaf area index
        REAL, INTENT (IN)  ::  COSZEN( NX, NY )    !  cosine of zenith angle
        REAL, INTENT (IN)  ::  PRES  ( NX, NY )    !  surface pressure (mb)

        REAL, INTENT (OUT)  ::  EMPOL ( NX, NY, BSPCS )      !  output pol emissions

C...........   SCRATCH LOCAL VARIABLES and their descriptions:

        INTEGER         R, C, L      !  counters
        REAL            CFOTHR       !  isop corr fac -- non-forest
        REAL            CFCLAI       !  ISOP CORR FAC -- LAI
        REAL            CFNO         !  NO correction factor
        REAL            CFOVOC       !  non-isop corr fac
        REAL            PAR          !  photo. actinic flux (UE/M**2-S)
        REAL            CT, DT       !  temperature correction
        REAL            TAIR         ! surface temperature
        REAL            RK           !  k from Geron and Guenther
        REAL            CSUBL        !  C sub l
        REAL            TLAI         !  temporary storage of LAI
        REAL            SOLTMP       !  temporary storage of radiation
        REAL            PSFC         !  temporary storage of sfc pressure (mb)
        REAL            ZEN          !  zenith angle
        REAL            PARDB        !  par direct beam
        REAL            PARDIF       !  par diffuse

        CHARACTER*16 :: PROGNAME = 'HRBEIS3'   !  program name

        CHARACTER*256   MESG

C***********************************************************************
C   begin body of subroutine  HRBIOS

C........... Find indices from BIOSPC array

            NISO  = INDEX1 ('ISOP', BSPCS, BIOSPC)
            NNO   = INDEX1 ( 'NO' , BSPCS, BIOSPC) 
            NMONO = INDEX1 ('MONO', BSPCS, BIOSPC) 
            NOVOC = INDEX1 ('OVOC', BSPCS, BIOSPC) 

             
C...........   loop thru cells

            DO R = 1, NY
               DO C = 1, NX
                
                  TAIR = TA( C, R )         ! unit in degree K


C..........    Perform checks on max and min bounds for temperature

                  IF (TAIR .LT. 200.0) THEN

                      WRITE( MESG, 94010 )
     &                 'TAIR=', TAIR,
     &                 'out of range at (C,R)=',  C, R
                      CALL M3EXIT( PROGNAME, JDATE, JTIME, MESG, 2 )

                  END IF

                  IF (TAIR .GT. 315.0 ) THEN
                      WRITE( MESG, 94020 )
     &                 'TAIR=', TAIR,
     &                 'out of range at (C,R)=',  C, R,
     &                 ' resetting to 315K'
                      CALL M3WARN( PROGNAME, JDATE, JTIME, MESG )
                      TAIR = 315.0
                  ENDIF

C..............  Calculate temperature correction term

                  DT   = 28668.514 / TAIR
                  CT   = EXP( 37.711 - 0.398570815 * DT ) /
     &                   (1.0 + EXP( 91.301 - DT ) )

                  SOLTMP   = TSOLAR( C, R )

C...................   cosine of zenith angle to zenith angle (radians)

                  ZEN =  ACOS( COSZEN( C, R ) ) 
                  TLAI = SLAI( C, R )

                  IF ( TLAI .GT. 10.0 ) THEN
                      WRITE( MESG, 94010 )
     &                 'LAI=', TLAI, 
     &                 'out of range at (C,R)=',  C, R
                      CALL M3EXIT( PROGNAME, JDATE, JTIME, MESG, 2 )
                  ENDIF

                  PSFC = PRES( C, R )

                  CALL GETPAR( SOLTMP, PSFC, ZEN, PARDB, PARDIF )
          
                  PAR = PARDB + PARDIF

C................... Check max/min bounds of PAR and calculate
C................... biogenic ISOP          

                  IF ( PAR .LT. 0.00 .OR. PAR .GT. 2500.0 ) THEN

                      WRITE( MESG, 94010 )
     &                 'PAR=', PAR, 
     &                 'out of range at (C,R)=',  C, R
                      CALL M3EXIT( PROGNAME, JDATE, JTIME, MESG, 2 )

                  ENDIF

C................... Initialize csubl 
                  CSUBL = 0.0

                  IF ( PAR .LE. 0.01 .OR. COSZEN( C,R) .LE. 0.02079483 ) 
     &                THEN

                      EMPOL( C,R, NISO ) = 0.0

                  ELSE
                     
                      IF ( TLAI .GT. 0.1 ) THEN 

                         CSUBL = CLNEW( ZEN, PARDB, PARDIF, TLAI )

                      ELSE  ! keep this or not?

                         CSUBL  = CGUEN( PAR ) 
                         
                      ENDIF
 
                      EMPOL( C,R, NISO ) = SISOP( C,R ) * CT * CSUBL

                  ENDIF

C..............  calculate OVOC and MONO emissions ; not speciated here
           
                  CFOVOC = EXP( 0.09 * ( TAIR - 303.0 ) )

                  EMPOL( C,R, NMONO ) = SMONO( C,R ) * CFOVOC
                  EMPOL( C,R, NOVOC ) = SOVOC( C,R ) * CFOVOC
  
C............. calculate NO emissions by going thru temperature cases

                  IF ( TAIR .GT. 303.00 ) TAIR = 303.00

                  IF ( TAIR .GT. 268.8690 ) THEN  !  most frequent case first

                       CFNO = EXP( 0.05112 * TAIR  -  15.68248 ) !  agriculture
                       EMPOL( C,R, NNO ) =  SNO( C,R ) * CFNO 
                  ELSE
                       EMPOL( C,R, NNO ) = 0.00

                  END IF      

               ENDDO
            ENDDO

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Informational (LOG) message formats... 92xxx


C...........   Internal buffering formats............ 94xxx


94010   FORMAT( A, F10.2, 1X, A, I3, ',', I3 )
94020   FORMAT( A, F10.2, 1X, A, I3, ',', I3, A )

C***************** CONTAINS ********************************************
            CONTAINS

            REAL FUNCTION CLNEW( ZEN, PARDB, PARDIF, TLAI )

C........  Function to calculate csubl based on zenith angle, par, and lai

            IMPLICIT NONE

            REAL, INTENT (IN) :: PARDB    ! direct beam PAR( umol/m2-s)
            REAL, INTENT (IN) :: PARDIF   ! diffuse PAR ( umol/m2-s)
            REAL, INTENT (IN) :: ZEN      ! solar zenith angle (radians)
            REAL, INTENT (IN) :: TLAI     ! leaf area index for grid cell
            REAL KBE                ! extinction coefficient for direct beam
            REAL CANPARSCAT         ! exponentially wtd scattered PAR (umol/m2-s)
            REAL CANPARDIF          ! exponentially wtd diffuse PAR (umol/m2-s)
            REAL PARSHADE           ! PAR on shaded leaves (umol/m2-s)
            REAL PARSUN             ! PAR on sunlit leaves (umol/m2-s)
            REAL LAISUN             ! LAI that is sunlit
            REAL FRACSUN            ! fraction of leaves that are sunlit
            REAL FRACSHADE          ! fraction of leaves that are shaded

C...........  CN98 - eqn 15.4, assume x=1

            KBE = 0.5 * SQRT(1. + TAN( ZEN ) * TAN( ZEN ))
 
C..........   CN98 - p. 261 (this is usually small)

            CANPARSCAT = 0.5 * PARDB * (EXP(-0.894 * KBE * TLAI) - 
     &               EXP(-1.* KBE * TLAI))

C..........   CN98 - p. 261 (assume exponentially wtd avg)

            CANPARDIF  = PARDIF * (1. - EXP(-0.61 * TLAI))/(0.61 * TLAI)

C.........    CN98 - p. 261 (for next 3 eqns)

            PARSHADE   = CANPARDIF + CANPARSCAT
	    PARSUN     = KBE * (PARDB + PARDIF) + PARSHADE
	    LAISUN     = (1. - EXP( -KBE * TLAI))/KBE
	    FRACSUN    = LAISUN/TLAI
	    FRACSHADE  = 1. - FRACSUN

C..........  cguen is guenther's eqn for computing light correction as a 
C..........  function of PAR...fracSun should probably be higher since 
C..........  sunlit leaves tend to be thicker than shaded leaves.  But 
C..........  since we need to make crude asmptns regarding leave 
C..........  orientation (x=1), will not attempt to fix at the moment.

            CLNEW =FRACSUN * CGUEN(PARSUN) + FRACSHADE * CGUEN(PARSHADE)

            RETURN 
            END FUNCTION CLNEW

            REAL FUNCTION CGUEN( PARTMP ) 

C..........  Guenther's equation for computing light correction

            IMPLICIT NONE
            REAL, INTENT (IN) :: PARTMP
            REAL, PARAMETER :: ALPHA2 = 0.00000729 

            IF ( PARTMP .LE. 0.01) THEN
               CGUEN = 0.0
            ELSE
               CGUEN = (0.0028782 * PARTMP) /
     &             SQRT(1. + ALPHA2 * PARTMP * PARTMP)
            ENDIF
 
            RETURN
            END FUNCTION CGUEN


        END SUBROUTINE HRBEIS3

