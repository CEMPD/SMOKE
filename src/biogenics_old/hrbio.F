
C Version "@(#)$Id$ $Source$ $Date$ 

        SUBROUTINE  HRBIO( JDATE, JTIME, LAT, LON, PRES, TA, 
     &                     PINE, DECD, CONF, LAI, OTHR, AVLAI, NORNO,
     &                     EMIS )

C***********************************************************************
C  subroutine body starts at line  160
C
C  FUNCTION:
C
C   SOLBIO part:
C     Computes the radiation terms needed for biogenic calculations.
C     Two terms are computed: TOTAL SOLAR RADIATION (LANGLEYS/MIN)
C     and PAR (UE/M**2-S).  The methodology for this calculation
C     was taken from IQBAL, M., 1983. AN INTRODUCTION TO SOLAR
C     RADIATION. ACADEMIC PRESS, NEW YORK, PP. 202-210.
C
C     Development of SOLBIO was prompted by the need for a
C     horizontal rather than actinic flux calculation which had
C     been performed by soleng.  Furthermore, soleng computed total
C     radiation only out to the near-IR spectrum.  This program
C     is designed only for approximate radiation estimates to be used
C     for biogenic emission calculations.
C
C   HRBIO part:
C     Computes hourly correction factors for various emitted species,
C     for program TMPBIO
C
C   HRLYEM part:
C     Computes hourly gridded natural emissions from normalized
C     emission inputs by applying correction factors for temperature
C     and solar energy.
C
C  PRECONDITIONS REQUIRED:
C     surface pressure PRES
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C	CZANGLE, DAYMON, CLDATEN
C
C  REVISION  HISTORY:
C       Prototype 12/95 by CJC adapted from UAM BEIS2 subroutines SOLBIO, 
C       CORRECT, and HRLYEM.
C
C***********************************************************************

      IMPLICIT NONE

C...........   INCLUDES:

      INCLUDE 'PARMS3.EXT'      ! I/O API constants
      INCLUDE 'FDESC3.EXT'      ! I/O API file description data structure
      INCLUDE 'IODECL3.EXT'     ! I/O API function declarations
      INCLUDE 'CONST3.EXT'      ! physical and mathematical constants
      INCLUDE 'BIODIMS3.EXT'    ! biogenic-related constants
      INCLUDE 'GRDIMS3.EXT'     ! grid-constants


C...........   ARGUMENTS and their descriptions:

        INTEGER         JDATE   !  current simulation date (YYYYDDD)
        INTEGER         JTIME   !  current simulation time (HHMMSS)
        REAL            LAT  ( NCOLS, NROWS )  !  lat (deg) -90 <= LAT <= 90
        REAL            LON  ( NCOLS, NROWS )  !  lon (deg) -180 <= LON <= 180 
        REAL            PRES ( NCOLS, NROWS )  !  pressure (Pa)
        REAL            TA   ( NCOLS, NROWS )  !  air temperature (K)
        REAL            PINE ( NCOLS, NROWS, BSPCS-1 )     !  nor VOC emissions
        REAL            DECD ( NCOLS, NROWS, BSPCS-1 )     !  nor VOC emissions
        REAL            CONF ( NCOLS, NROWS, BSPCS-1 )     !  nor VOC emissions
        REAL            LAI  ( NCOLS, NROWS, BSPCS-1 )     !  nor VOC emissions
        REAL            OTHR ( NCOLS, NROWS, BSPCS-1 )     !  nor VOC emissions
        REAL            AVLAI( NCOLS, NROWS )              !  average LAI
        REAL            NORNO( NCOLS, NROWS, LUSES )       !  nor NO  emissions
        REAL            EMIS ( NCOLS, NROWS, MSPCS )       !  output  emissions


C...........   PARAMETERS and their descriptions:

      INTEGER        CLAYS              !  dim:  canopy layers
      INTEGER        LPINE, LDECD, LCONF
      REAL           WM2LY, LY2UE, CN, PRES0, DPRES0, SOLFAC
      REAL           D28, D29, D30, D31
      REAL           ALPHA

      PARAMETER    ( CLAYS  =    5, 
     &               LPINE  =    1 ,
     &               LDECD  =    2 ,
     &               LCONF  =    3 , 
     &               WM2LY  =    0.001433,
     &               LY2UE  = 2916.0,
     &               CN     =    1.0,
     &               PRES0  = 1013.0 ,
     &               SOLFAC = 0.5 * LY2UE * WM2LY * CN ,
     &               DPRES0 = 1.0 / PRES0 ,
     &               D28    = 1.0 / 28.0 ,
     &               D29    = 1.0 / 29.0 ,
     &               D30    = 1.0 / 30.0 ,
     &               D31    = 1.0 / 31.0 ,
     &               ALPHA = 0.00000729  )
 

C...........   SCRATCH LOCAL VARIABLES and their descriptions:

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
        INTEGER         JDAY, JMON, R, C, L
        LOGICAL         GETATN
        REAL            AA, BB, CC
        REAL            P, Q
        REAL            EXPA
        REAL            COSZEN( NCOLS, NROWS )
        REAL            CLDATN( NCOLS, NROWS )
        REAL            PAR             !  photo. actinic flux (UE/M**2-S)
        REAL            PARZ, SQPARZ, CT, DT
        REAL            TAIR
        REAL            FPINE
        REAL            FDECD
        REAL            FCONF
        REAL            FCLAI
        REAL            BISOP, BMONO, BOVOC        !  biogenic VOC spcs emis.

C...........   SAVED LOCAL VARIABLES and their descriptions:

        LOGICAL      FIRSTIME
        REAL         DMON(   12 )       !  inverse length of month
        REAL         ADAY( 0:13 )
        REAL         BDAY( 0:13 )
        REAL         CDAY( 0:13 )
        REAL         CANATN( CLAYS, TREETY )    !  Canopy attenuation factors
        REAL         BIOFRA( CLAYS )            !  biomass fractions

        DATA FIRSTIME  / .TRUE. /

        DATA DMON
     & / D31, D28, D31, D30, D31, D30, 
     &   D31, D31, D30, D31, D30, D31 /
 
        DATA ADAY
     & /1203.0, 1202.0, 1187.0, 1164.0, 1130.0, 1106.0, 1092.0,
     &  1093.0, 1107.0, 1136.0, 1136.0, 1190.0, 1204.0, 1203.0 /
 
        DATA BDAY
     & / 0.141, 0.141,  0.142,  0.149,  0.164,  0.177,  0.185,
     &   0.186, 0.182,  0.165,  0.152,  0.144,  0.141,  0.141 /
 
        DATA CDAY
     & / 0.103, 0.103,  0.104,  0.109,  0.120,  0.130,  0.137,
     &   0.138, 0.134,  0.121,  0.111,  0.106,  0.103,  0.103 /

        DATA CANATN  / 0.88, 0.69, 0.53, 0.41, 0.32,    !  pine
     &                 0.81, 0.53, 0.35, 0.23, 0.15,    !  deciduous
     &                 0.75, 0.41, 0.23, 0.13, 0.07 /   !  coniferous
      
        DATA BIOFRA  / 0.27, 0.21, 0.18, 0.17, 0.17 /

        SAVE   FIRSTIME, ADAY, BDAY, CDAY, CANATN, BIOFRA
 

C***********************************************************************
C   begin body of subroutine  HRBIO

C.......   Fold sfac into ADAY, dpres0 into BDAY

        IF( FIRSTIME ) THEN

            FIRSTIME = .FALSE.

            DO  11  JMON = 1, 14
                ADAY( JMON ) = SOLFAC * ADAY( JMON )
                BDAY( JMON ) = DPRES0 * BDAY( JMON )
11          CONTINUE

        END IF          !  if firstime:  fold sfac into ADAY

C.......   Begin processing for this hour:  look at cos of 
C.......   zenith angle, and whether sun is above the horizon:

        CALL CZANGLE( JDATE, JTIME, LAT, LON, COSZEN, GETATN )

        IF ( GETATN ) THEN !  sun above horizon:; compute direct radiation: 

C.......   Interpolate [ A | B | C ]DAY to this JDATE:

            CALL DAYMON( JDATE, JMON, JDAY )
            P  = DMON( JMON ) * FLOAT( JDAY - 1 )
            Q  = 1.0 - P
            AA = Q * ADAY( JMON-1 )  +  P * ADAY( JMON )
            BB = Q * BDAY( JMON-1 )  +  P * BDAY( JMON )
            CC = Q * CDAY( JMON-1 )  +  P * CDAY( JMON )

            CALL CLDATEN( JDATE, JTIME, LAT, LON, CLDATN )

C...........   compute direct horizontal radiation, attenuation:
C...........   first, perform the table look up

            DO  122  R = 1, NROWS
            DO  111  C = 1, NCOLS
                
                TAIR = TA( C,R )

                DT   = 28668.514 / TAIR
                CT   = EXP( 37.711 - 0.398570815 * DT ) /
     &                 (1.0 + EXP( 91.301 - DT ) )

                IF ( COSZEN( C,R ) .LT. 0.02 ) THEN     !  effectively no PAR

C...................   Biogenic ISOP emissions (due to photo. activity) zero

                    EMIS( C,R, NISO ) = 0.0

                ELSE

C...................   Compute biogenic VOC emissions, due to
C...................   photosynthetic activity:

                    EXPA  = AA * EXP( -BB*PRES( C,R ) / COSZEN( C,R ) )
                    PAR   = EXPA*CLDATN( C,R ) * ( COSZEN( C,R ) + CC )
                    FPINE = 0.0
                    FDECD = 0.0
                    FCONF = 0.0
                    FCLAI = 0.0
                    
                    DO  101  L = 1, CLAYS

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

                        PARZ   = PAR * EXP(-0.042 * AVLAI( C,R )
     &                                            * FLOAT( 2 * L - 1 ) )
                        SQPARZ = PARZ * PARZ
                        FCLAI  = FCLAI + PARZ/SQRT( 1.0 + ALPHA*SQPARZ )

101                 CONTINUE
                    
                    CFPINE = CT * FPINE
                    CFDECD = CT * FDECD
                    CFCONF = CT * FCONF
                    CFCLAI = CT * 0.002878 * 0.2 * FCLAI
                    CFOTHR = CT * 0.002878 * PAR 
     &                          / SQRT( 1.0 + ALPHA * PAR * PAR )

                    BISOP = PINE( C,R,ISOP ) * CFPINE +
     &                      DECD( C,R,ISOP ) * CFDECD +
     &                      CONF( C,R,ISOP ) * CFCONF +
     &                      LAI ( C,R,ISOP ) * CFCLAI +
     &                      OTHR( C,R,ISOP ) * CFOTHR

                    EMIS( C,R, NISO ) = DENISO   * BISOP

                END IF                  !  if coszen > 0.02 or not
                    
                CFOVOC = EXP( 0.09 * ( TAIR - 303.0 ) )

                BMONO = ( PINE( C,R,MONO ) +
     &                    DECD( C,R,MONO ) +
     &                    CONF( C,R,MONO ) +
     &                    LAI ( C,R,MONO ) +
     &                    OTHR( C,R,MONO )  ) * CFOVOC

                BOVOC = ( PINE( C,R,OVOC ) +
     &                    DECD( C,R,OVOC ) +
     &                    CONF( C,R,OVOC ) +
     &                    LAI ( C,R,OVOC ) +
     &                    OTHR( C,R,OVOC )  ) * CFOVOC

                EMIS( C,R, NPAR ) = MON2PAR  * BMONO +  OTH2PAR * BOVOC
                EMIS( C,R, NOLE ) = MON2OLE  * BMONO +  OTH2OLE * BOVOC
                EMIS( C,R, NALD ) = MON2ALD2 * BMONO


C...........   Compute correction factors for NO, NO emissions:
C...........   In the following, simultaneously solve inequalities for 
C...........   TS in all 4 land use categories; then do substitution of 
C...........   the various A * TAIR + B into the  0.071 * (TS(I) - 303.16)
C...........   formula for CFNO(i) to arrive at the final version:
C
C                TS(GRAS) = 0.66*TAIR + 101.67
C                TS(FORE) = 0.84*TAIR + 47.31
C                TS(WETL) = 0.92*TAIR + 26.25
C                TS(AGRI) = 0.72*TAIR + 82.28
C                DO 80 I= 1, LUSES
C                  IF (TS(I) .LE. 273.16) THEN
C                     CFNO(I) = 0.0
C                  ELSE
C                     CFNO(I) = EXP(0.071 * (TS(I) - 303.16))
C                  END IF
C80              CONTINUE
               
                IF ( TAIR .GT. 268.8690 ) THEN  !  most frequent case first

                    CFNOG = EXP( 0.04686 * TAIR  -  14.30579 ) !  grass
                    CFNOF = EXP( 0.05964 * TAIR  -  18.16535 ) !  forest
                    CFNOW = EXP( 0.06532 * TAIR  -  19.66061 ) !  wetland
                    CFNOA = EXP( 0.05112 * TAIR  -  15.82480 ) !  agriculture

                    EMIS( C,R, NNO ) = ( NORNO( C,R,GRAS ) * CFNOG +
     &                                   NORNO( C,R,FORE ) * CFNOF +
     &                                   NORNO( C,R,WETL ) * CFNOW +
     &                                   NORNO( C,R,AGRI ) * CFNOA
     &                                 ) * DENISO

                ELSE IF ( TAIR .GT. 268.3804 ) THEN     !  no forest NO

                    CFNOG = EXP( 0.04686 * TAIR  -  14.30579 ) !  grass
                    CFNOW = EXP( 0.06532 * TAIR  -  19.66061 ) !  wetland
                    CFNOA = EXP( 0.05112 * TAIR  -  15.82480 ) !  agriculture

                    EMIS( C,R, NNO ) = ( NORNO( C,R,GRAS ) * CFNOG +
     &                                   NORNO( C,R,WETL ) * CFNOW +
     &                                   NORNO( C,R,AGRI ) * CFNOA
     &                                 ) * DENISO

                ELSE IF ( TAIR .GT. 265.11111 ) THEN    !  no forest, wet NO

                    CFNOG = EXP( 0.04686 * TAIR  -  14.30579 ) !  grass
                    CFNOA = EXP( 0.05112 * TAIR  -  15.82480 ) !  agriculture

                    EMIS( C,R, NNO ) = ( NORNO( C,R,GRAS ) * CFNOG +
     &                                   NORNO( C,R,AGRI ) * CFNOA
     &                                 ) * DENISO

                ELSE IF ( TAIR .GT. 259.8333 ) THEN     ! NO from grass only

                    CFNOG = EXP( 0.04686 * TAIR  -  14.30579 ) !  grass

                    EMIS( C,R, NNO ) = NORNO( C,R,GRAS ) * CFNOG*DENISO

                ELSE       !  else tair <= 259.8333

                    EMIS( C,R, NNO ) = 0.0

                END IF          !  5-way conditional on TAIR

111         CONTINUE
122         CONTINUE

        ELSE    !  not GETATN:
        
            DO 222 R = 1,NROWS
            DO 211 C = 1,NCOLS
                   
C...........   Sun below horizon; set photosynthetic emissions factors 
C...........   (and therefore biogenic ISOP emissions) to zero

                TAIR   = TA( C,R )

                CFOVOC = EXP( 0.09 * ( TAIR - 303.0 ) )

                BMONO  = ( PINE( C,R,MONO ) +
     &                     DECD( C,R,MONO ) +
     &                     CONF( C,R,MONO ) +
     &                     LAI ( C,R,MONO ) +
     &                     OTHR( C,R,MONO )  ) * CFOVOC

                BOVOC  = ( PINE( C,R,OVOC ) +
     &                     DECD( C,R,OVOC ) +
     &                     CONF( C,R,OVOC ) +
     &                     LAI ( C,R,OVOC ) +
     &                     OTHR( C,R,OVOC )  ) * CFOVOC

                EMIS( C,R, NPAR ) = MON2PAR  * BMONO +  OTH2PAR * BOVOC
                EMIS( C,R, NOLE ) = MON2OLE  * BMONO +  OTH2OLE * BOVOC
                EMIS( C,R, NALD ) = MON2ALD2 * BMONO
                EMIS( C,R, NISO ) = 0.0

C...........   Compute correction factors for NO (see note, above)

                IF ( TAIR .GT. 268.8690 ) THEN  !  most frequent case first

                    CFNOG = EXP( 0.04686 * TAIR  -  14.30579 ) !  grass
                    CFNOF = EXP( 0.05964 * TAIR  -  18.16535 ) !  forest
                    CFNOW = EXP( 0.06532 * TAIR  -  19.66061 ) !  wetland
                    CFNOA = EXP( 0.05112 * TAIR  -  15.82480 ) !  agriculture

                    EMIS( C,R, NNO ) = ( NORNO( C,R,GRAS ) * CFNOG +
     &                                   NORNO( C,R,FORE ) * CFNOF +
     &                                   NORNO( C,R,WETL ) * CFNOW +
     &                                   NORNO( C,R,AGRI ) * CFNOA
     &                                 ) * DENISO

                ELSE IF ( TAIR .GT. 268.3804 ) THEN     !  no forest NO

                    CFNOG = EXP( 0.04686 * TAIR  -  14.30579 ) !  grass
                    CFNOW = EXP( 0.06532 * TAIR  -  19.66061 ) !  wetland
                    CFNOA = EXP( 0.05112 * TAIR  -  15.82480 ) !  agriculture

                    EMIS( C,R, NNO ) = ( NORNO( C,R,GRAS ) * CFNOG +
     &                                   NORNO( C,R,WETL ) * CFNOW +
     &                                   NORNO( C,R,AGRI ) * CFNOA
     &                                 ) * DENISO

                ELSE IF ( TAIR .GT. 265.11111 ) THEN    !  no forest, wetl NO

                    CFNOG = EXP( 0.04686 * TAIR  -  14.30579 ) !  grass
                    CFNOA = EXP( 0.05112 * TAIR  -  15.82480 ) !  agriculture

                    EMIS( C,R, NNO ) = ( NORNO( C,R,GRAS ) * CFNOG +
     &                                   NORNO( C,R,AGRI ) * CFNOA
     &                                 ) * DENISO

                ELSE IF ( TAIR .GT. 259.8333 ) THEN     !  only grass NO

                    CFNOG = EXP( 0.04686 * TAIR  -  14.30579 ) !  grass

                    EMIS( C,R, NNO ) = NORNO( C,R,GRAS ) * CFNOG*DENISO

                ELSE            !  tair <= 259.8333; no NO

                    EMIS( C,R, NO ) = 0.0

                END IF          !  5-way conditional on TAIR

211         CONTINUE
222         CONTINUE

        END IF          !  if getatn (sun above horizon) or not

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Error and warning message formats..... 91xxx

91000   FORMAT ( //5X , '*** ERROR ABORT in subroutine HRBIO ***',
     &            /5X , A ,
     &           // )        !  generic error message format


C...........   Informational (LOG) message formats... 92xxx

92000   FORMAT ( 5X , A )


C...........   Formatted file I/O formats............ 93xxx


C...........   Internal buffering formats............ 94xxx


C...........   Miscellaneous formats................. 95xxx

95000   FORMAT ( /5X , A , $ )          !  generic prompt format.


        END

