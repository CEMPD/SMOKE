 
        SUBROUTINE  HRBIOS( JDATE, JTIME, NX, NY, 
     &                      PINE, DECD, CONF, AGRC, LAI, OTHR, AVLAI,
     &                      NORNO, TA, TSOLAR, EMPOL )

C***********************************************************************
C  subroutine body starts at line  143
C
C  DESCRIPTION:
C  
C     Uses PAR and sfc temperature data to calculate
C     biogenic ISOP emissions.  OVOC, TERP and NO emissions are
C     calculated using the temperature data only.  OVOC and TERP
C     are not speciated in this rountine. 
C
C  PRECONDITIONS REQUIRED:
C     PAR and Surface Temperature
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C     m3io 
C
C  REVISION  HISTORY:
C     11/99: by Jeff Vukovich taken from v4.3 hrbios.F SMOKE prototype
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
C***********************************************************************

C...........   Modules for public variables
C...........   This module contains the global variables for the 3-d grid
        USE MODGRID

        IMPLICIT NONE

C...........   INCLUDES:

        INCLUDE 'PARMS3.EXT'      ! I/O API constants
        INCLUDE 'FDESC3.EXT'      ! I/O API file description data structure
        INCLUDE 'IODECL3.EXT'     ! I/O API function declarations
        INCLUDE 'EMCNST3.EXT'     !
        INCLUDE 'BIODIMS3.EXT'    ! biogenic-related constants

C...........   EXTERNAL FUNCTIONS and their descriptions:

        INTEGER         INDEX1

        EXTERNAL        INDEX1

C...........   ARGUMENTS and their descriptions:

        INTEGER, INTENT (IN)  :: JDATE   !  current simulation date (YYYYDDD)
        INTEGER, INTENT (IN)  :: JTIME   !  current simulation time (HHMMSS)
        INTEGER, INTENT (IN)  :: NX      !  no. met columns
        INTEGER, INTENT (IN)  :: NY      !  no. met rows

        REAL, INTENT (IN) :: TA    ( NX, NY )  !  air temperature (K)
        REAL, INTENT (IN) :: TSOLAR( NX, NY )  !  PAR
        REAL, INTENT (IN) :: PINE  ( NCOLS, NROWS, BSPCS-1 ) ! nor VOC emissions
        REAL, INTENT (IN) :: DECD  ( NCOLS, NROWS, BSPCS-1 ) ! nor VOC emissions
        REAL, INTENT (IN) :: CONF  ( NCOLS, NROWS, BSPCS-1 ) ! nor VOC emissions
        REAL, INTENT (IN) :: AGRC  ( NCOLS, NROWS, BSPCS-1 ) ! nor VOC emissions
        REAL, INTENT (IN) :: LAI   ( NCOLS, NROWS, BSPCS-1 ) ! nor VOC emissions
        REAL, INTENT (IN) :: OTHR  ( NCOLS, NROWS, BSPCS-1 ) ! nor VOC emissions
        REAL, INTENT (IN) :: AVLAI ( NCOLS, NROWS )          ! average LAI
        REAL, INTENT (IN) :: NORNO ( NCOLS, NROWS, LUSES )   ! nor NO  emissions

        REAL, INTENT(OUT) :: EMPOL ( NCOLS, NROWS, BSPCS )   ! output pol emissions


C...........   PARAMETERS and their descriptions:

      INTEGER        CLAYS              !  dim:  canopy layers
      INTEGER        LPINE, LDECD, LCONF
      REAL           CN, PRES0, DPRES0
      REAL           ALPHA

      PARAMETER    ( CLAYS  =    5, 
     &               LPINE  =    1 ,
     &               LDECD  =    2 ,
     &               LCONF  =    3 , 
     &               CN     =    1.0,
     &               PRES0  = 1013.0 ,
     &               DPRES0 = 1.0 / PRES0 ,
     &               ALPHA = 0.00000729  )
 

C...........   SCRATCH LOCAL VARIABLES and their descriptions:
        INTEGER         I, J      !  indices
        INTEGER         R, C, L

        REAL            CFPINE    !  isop corr fac -- pine
        REAL            CFDECD    !  isop corr fac -- deciduous
        REAL            CFCONF    !  isop corr fac -- coniferous
        REAL            CFOTHR    !  isop corr fac -- non-forest
        REAL            CFCLAI    !  ISOP CORR FAC -- LAI
        REAL            CFNOG     !  NO   corr fac for land use GRAS 
        REAL            CFNOF     !  NO   corr fac for land use FORE
        REAL            CFNOW     !  NO   corr fac for land use WETL
        REAL            CFNOA     !  NO   corr fac for land use AGRI
        REAL            CFOVOC    !  non-isop corr fac
        REAL            PAR             !  photo. actinic flux (UE/M**2-S)
        REAL            PARZ, SQPARZ, CT, DT
        REAL            TAIR
        REAL            FPINE
        REAL            FDECD
        REAL            FCONF
        REAL            FCLAI
        REAL            BISOP, BMONO, BOVOC        !  biogenic VOC spcs emis.
        CHARACTER*240   MESG

        CHARACTER*16 :: PROGNAME = 'HRBIOS'   !  program name

C...........   SAVED LOCAL VARIABLES and their descriptions:

        REAL         CANATN( CLAYS, TREETY )    !  Canopy attenuation factors
        REAL         BIOFRA( CLAYS )            !  biomass fractions

        DATA   CANATN  / 0.88, 0.69, 0.53, 0.41, 0.32,    !  pine
     &                   0.81, 0.53, 0.35, 0.23, 0.15,    !  deciduous
     &                   0.75, 0.41, 0.23, 0.13, 0.07 /   !  coniferous
      
        DATA   BIOFRA  / 0.27, 0.21, 0.18, 0.17, 0.17 /

        SAVE   CANATN, BIOFRA

C***********************************************************************
C   begin body of subroutine  HRBIOS

C........... Find indices from BIOSPC array

            NISO = INDEX1 ( 'ISOP', BSPCS, BIOSPC)
            NNO  = INDEX1 ( 'NO' , BSPCS, BIOSPC) 
            NTERP = INDEX1 ('TERP', BSPCS, BIOSPC) 
            NOVOC = INDEX1 ('OVOC', BSPCS, BIOSPC) 
 
C...........   loop thru cells

            DO  R = 1, NY
            DO  C = 1, NX
                
C.................  Adjust for subgrid
                I = C - XOFF
                J = R - YOFF

                IF( I .LE. 0 .OR. I .GT. NCOLS .OR.
     &              J .LE. 0 .OR. J .GT. NROWS       ) CYCLE

                TAIR = TA(C,R)         ! unit in degree K

C..........    Perform checks on max and min bounds for temperature

                IF (TAIR .LT. 200.0) THEN

                    WRITE( MESG, 94010 )
     &               'TAIR=', TAIR,
     &               'out of range at (C,R)=',  C, R
                    CALL M3EXIT( PROGNAME, JDATE, JTIME, MESG, 2 )

                END IF

                IF (TAIR .GT. 315.0 ) THEN
                    WRITE( MESG, 94020 )
     &               'TAIR=', TAIR,
     &               'out of range at (C,R)=',  C, R,
     &               ' resetting to 315K'
                    CALL M3WARN( PROGNAME, JDATE, JTIME, MESG )
                    TAIR = 315.0
                ENDIF

C..............  Calculate temperature correction term

                DT   = 28668.514 / TAIR
                CT   = EXP( 37.711 - 0.398570815 * DT ) /
     &                 (1.0 + EXP( 91.301 - DT ) )

C...................   Compute biogenic VOC emissions, due to
C...................   photosynthetic activity:

                PAR   = TSOLAR( C, R )

C................... Check max/min bounds of PAR and calculate
C................... biogenic ISOP          
 
                IF (PAR .LT. 0.0 .OR. PAR .GT. 2500.0) THEN

                    WRITE( MESG, 94010 )
     &               'PAR=', PAR, 
     &               'out of range at (C,R)=',  C, R
                    CALL M3EXIT( PROGNAME, JDATE, JTIME, MESG, 2 )

                ELSE IF ( PAR .GT. 0.0 ) THEN
                
                    FPINE = 0.0
                    FDECD = 0.0
                    FCONF = 0.0
                    FCLAI = 0.0
C................ cycle thru canopy layers
               
                    DO  L = 1, CLAYS
                        PARZ   = PAR * CANATN( L,LPINE )
                        SQPARZ = PARZ * PARZ
                        FPINE  = FPINE + 0.2 * 0.002878 * PARZ 
     &                                  / SQRT( 1.0 + ALPHA * SQPARZ )
               
                        PARZ   = PAR * CANATN( L,LDECD )
                        SQPARZ = PARZ * PARZ
                        FDECD  = FDECD + BIOFRA(L) * 0.002878 * PARZ 
     &                                  / SQRT( 1.0 + ALPHA * SQPARZ )
               
                        PARZ   = PAR * CANATN( L,LCONF )
                        SQPARZ = PARZ * PARZ
                        FCONF  = FCONF + 0.2 * 0.002878 * PARZ 
     &                                  / SQRT( 1.0 + ALPHA * SQPARZ )
               
                        PARZ   = PAR * EXP(-0.042 * AVLAI( I,J )
     &                                            * FLOAT( 2 * L - 1 ) )
                        SQPARZ = PARZ * PARZ
                        FCLAI  = FCLAI + 0.2 * 0.002878 * PARZ
     &                                  / SQRT( 1.0 + ALPHA*SQPARZ )
               
                    ENDDO
                    
                    CFPINE = CT * FPINE
                    CFDECD = CT * FDECD
                    CFCONF = CT * FCONF
                    CFCLAI = CT * FCLAI
                    CFOTHR = CT * 0.002878 * PAR 
     &                          / SQRT( 1.0 + ALPHA * PAR * PAR )

                    BISOP = PINE( I,J,ISOP ) * CFPINE +
     &                      DECD( I,J,ISOP ) * CFDECD +
     &                      CONF( I,J,ISOP ) * CFCONF +
     &                      LAI ( I,J,ISOP ) * CFCLAI +
     &                      AGRC( I,J,ISOP ) * CFOTHR + 
     &                      OTHR( I,J,ISOP ) * CFOTHR
                
                    EMPOL( I,J, NISO ) = BISOP
                    
                ELSE
                
                    EMPOL( I,J, NISO ) = 0.0

                END IF

C..............  calculate OVOC and MONO emissions ; not speciated here
           
                CFOVOC = EXP( 0.09 * ( TAIR - 303.0 ) )

                BMONO = ( PINE( I,J,MONO ) +
     &                    DECD( I,J,MONO ) +
     &                    CONF( I,J,MONO ) +
     &                    AGRC( I,J,MONO ) +  ! sl
     &                    LAI ( I,J,MONO ) +
     &                    OTHR( I,J,MONO )  ) * CFOVOC

                BOVOC = ( PINE( I,J,OVOC ) +
     &                    DECD( I,J,OVOC ) +
     &                    CONF( I,J,OVOC ) +
     &                    AGRC( I,J,OVOC ) +  ! sl
     &                    LAI ( I,J,OVOC ) +
     &                    OTHR( I,J,OVOC )  ) * CFOVOC

                EMPOL( I,J, NTERP ) = BMONO
                EMPOL( I,J, NOVOC ) = BOVOC

C............. calculate NO emissions by going thru temperature cases

                IF ( TAIR .GT. 268.8690 ) THEN  !  most frequent case first

                    CFNOG = EXP( 0.04686 * TAIR  -  14.30579 ) !  grass
                    CFNOF = EXP( 0.05964 * TAIR  -  18.16535 ) !  forest
                    CFNOW = EXP( 0.06532 * TAIR  -  19.66061 ) !  wetland
                    CFNOA = EXP( 0.05112 * TAIR  -  15.68248 ) !  agriculture

                    EMPOL( I,J, NNO ) =  NORNO( I,J,GRAS ) * CFNOG +
     &                                   NORNO( I,J,FORE ) * CFNOF +
     &                                   NORNO( I,J,WETL ) * CFNOW +
     &                                   NORNO( I,J,AGRI ) * CFNOA

                ELSE IF ( TAIR .GT. 268.3804 ) THEN     !  no forest NO

                    CFNOG = EXP( 0.04686 * TAIR  -  14.30579 ) !  grass
                    CFNOW = EXP( 0.06532 * TAIR  -  19.66061 ) !  wetland
                    CFNOA = EXP( 0.05112 * TAIR  -  15.68248 ) !  agriculture

                    EMPOL( I,J, NNO ) =  NORNO( I,J,GRAS ) * CFNOG +
     &                                   NORNO( I,J,WETL ) * CFNOW +
     &                                   NORNO( I,J,AGRI ) * CFNOA

                ELSE IF ( TAIR .GT. 265.11111 ) THEN    !  no forest, wet NO

                    CFNOG = EXP( 0.04686 * TAIR  -  14.30579 ) !  grass
                    CFNOA = EXP( 0.05112 * TAIR  -  15.68248 ) !  agriculture

                    EMPOL( I,J, NNO ) =  NORNO( I,J,GRAS ) * CFNOG +
     &                                   NORNO( I,J,AGRI ) * CFNOA

                ELSE IF ( TAIR .GT. 259.8333 ) THEN     ! NO from grass only

                    CFNOG = EXP( 0.04686 * TAIR  -  14.30579 ) !  grass

                    EMPOL( I,J, NNO ) = NORNO( I,J,GRAS ) * CFNOG

                ELSE       !  else tair <= 259.8333

                    EMPOL( I,J, NNO ) = 0.0

                END IF          !  5-way conditional on TAIR

            ENDDO
            ENDDO

        RETURN
94010   FORMAT( A, F10.2, 1X, A, I3, ',', I3 )
94020   FORMAT( A, F10.2, 1X, A, I3, ',', I3, A )

        END

