 
        SUBROUTINE  HRBIOS( JDATE, JTIME,
     &                      PINE, DECD, CONF, AGRC, LAI, OTHR, AVLAI,
     &                      NORNO, TSFC, TSOLAR, EMIS )

C***********************************************************************
C  subroutine body starts at line  167
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
C     PAR and Surface Temperature
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C     m3io 
C
C  REVISION  HISTORY:
C       Prototype 12/95 by CJC adapted from UAM BEIS2 subroutines SOLBIO, 
C       CORRECT, and HRLYEM.
C
C       Version 8/96 by SL and CJC to use Dyntel-style ASCII solar radiation
C       and surface temperature file instead of MCIP meteorology.
C       Name change HRBIO ~~> HRBIOS to accompany arg-list change.
C
C***********************************************************************
C
C Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
C                System
C File Version: @(#)$Id$
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
C***************************************************************************/

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
        REAL            PINE ( NCOLS, NROWS, BSPCS-1 )     !  nor VOC emissions
        REAL            DECD ( NCOLS, NROWS, BSPCS-1 )     !  nor VOC emissions
        REAL            CONF ( NCOLS, NROWS, BSPCS-1 )     !  nor VOC emissions
        REAL            LAI  ( NCOLS, NROWS, BSPCS-1 )     !  nor VOC emissions
        REAL            AGRC ( NCOLS, NROWS, BSPCS-1 )     !  nor VOC emissions
        REAL            OTHR ( NCOLS, NROWS, BSPCS-1 )     !  nor VOC emissions
        REAL            AVLAI( NCOLS, NROWS )              !  average LAI
        REAL            NORNO( NCOLS, NROWS, LUSES )       !  nor NO  emissions
        REAL            EMIS ( NCOLS, NROWS, MSPCS )       !  output  emissions
        REAL            TSFC  ( NCOLS, NROWS )		   !  surface temp (C)
        REAL            TSOLAR( NCOLS, NROWS )		   !  PAR


C...........   PARAMETERS and their descriptions:

      INTEGER        CLAYS              !  dim:  canopy layers
      INTEGER        LPINE, LDECD, LCONF
      REAL           WM2LY, LY2UE, CN, PRES0, DPRES0
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
        INTEGER         R, C, L
        REAL            PAR             !  photo. actinic flux (UE/M**2-S)
        REAL            PARZ, SQPARZ, CT, DT
        REAL            TAIR
        REAL            FPINE
        REAL            FDECD
        REAL            FCONF
        REAL            FCLAI
        REAL            BISOP, BMONO, BOVOC        !  biogenic VOC spcs emis.
        CHARACTER*240   MESG

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

C...........   first, perform the table look up

            DO  122  R = 1, NROWS
            DO  111  C = 1, NCOLS
                
                TAIR = TSFC(C,R)         ! unit in degree K
                DT   = 28668.514 / TAIR
                CT   = EXP( 37.711 - 0.398570815 * DT ) /
     &                 (1.0 + EXP( 91.301 - DT ) )

C...................   Compute biogenic VOC emissions, due to
C...................   photosynthetic activity:

                PAR   = TSOLAR( C, R )
 
                IF (TAIR .LT. 200.0  .OR. TAIR .GT. 315.0) THEN

                    WRITE( MESG, 94010 ) 
     &               'TAIR=', TAIR, 
     &               'out of range at (C,R)=',  C, R
                    CALL M3EXIT( 'TMPBIO', JDATE, JTIME, MESG, 2 )

                END IF
 
                IF (PAR .LT. 0.0 .OR. PAR .GT. 2500.0) THEN

                    WRITE( MESG, 94010 )
     &               'PAR=', PAR, 
     &               'out of range at (C,R)=',  C, R
                    CALL M3EXIT( 'TMPBIO', JDATE, JTIME, MESG, 2 )

                ELSE IF ( PAR .GT. 0.0 ) THEN
                
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
                        FCLAI  = FCLAI + 0.2 * 0.002878 * PARZ
     &                                  / SQRT( 1.0 + ALPHA*SQPARZ )
               
101                 CONTINUE
                    
                    CFPINE = CT * FPINE
                    CFDECD = CT * FDECD
                    CFCONF = CT * FCONF
                    CFCLAI = CT * FCLAI
                    CFOTHR = CT * 0.002878 * PAR 
     &                          / SQRT( 1.0 + ALPHA * PAR * PAR )
               
                    BISOP = PINE( C,R,ISOP ) * CFPINE +
     &                      DECD( C,R,ISOP ) * CFDECD +
     &                      CONF( C,R,ISOP ) * CFCONF +
     &                      LAI ( C,R,ISOP ) * CFCLAI +
     &                      AGRC( C,R,ISOP ) * CFOTHR +  ! sl from hrlyem2.f
     &                      OTHR( C,R,ISOP ) * CFOTHR
                
                    EMIS( C,R, NISO ) = DENISO   * BISOP

                ELSE
                
                    EMIS( C,R, NISO ) = 0.0

                END IF
           
                CFOVOC = EXP( 0.09 * ( TAIR - 303.0 ) )

                BMONO = ( PINE( C,R,MONO ) +
     &                    DECD( C,R,MONO ) +
     &                    CONF( C,R,MONO ) +
     &                    AGRC( C,R,MONO ) +  ! sl
     &                    LAI ( C,R,MONO ) +
     &                    OTHR( C,R,MONO )  ) * CFOVOC

                BOVOC = ( PINE( C,R,OVOC ) +
     &                    DECD( C,R,OVOC ) +
     &                    CONF( C,R,OVOC ) +
     &                    AGRC( C,R,OVOC ) +  ! sl
     &                    LAI ( C,R,OVOC ) +
     &                    OTHR( C,R,OVOC )  ) * CFOVOC

                EMIS( C,R, NPAR ) = MON2PAR  * BMONO +  OTH2PAR * BOVOC
                EMIS( C,R, NOLE ) = MON2OLE  * BMONO +  OTH2OLE * BOVOC
                EMIS( C,R, NALD ) = MON2ALD2 * BMONO

                IF ( TAIR .GT. 268.8690 ) THEN  !  most frequent case first

                    CFNOG = EXP( 0.04686 * TAIR  -  14.30579 ) !  grass
                    CFNOF = EXP( 0.05964 * TAIR  -  18.16535 ) !  forest
                    CFNOW = EXP( 0.06532 * TAIR  -  19.66061 ) !  wetland
                    CFNOA = EXP( 0.05112 * TAIR  -  15.68248 ) !  agriculture

                    EMIS( C,R, NNO ) = ( NORNO( C,R,GRAS ) * CFNOG +
     &                                   NORNO( C,R,FORE ) * CFNOF +
     &                                   NORNO( C,R,WETL ) * CFNOW +
     &                                   NORNO( C,R,AGRI ) * CFNOA
     &                                 ) * DENNO     ! sl

                ELSE IF ( TAIR .GT. 268.3804 ) THEN     !  no forest NO

                    CFNOG = EXP( 0.04686 * TAIR  -  14.30579 ) !  grass
                    CFNOW = EXP( 0.06532 * TAIR  -  19.66061 ) !  wetland
                    CFNOA = EXP( 0.05112 * TAIR  -  15.68248 ) !  agriculture

                    EMIS( C,R, NNO ) = ( NORNO( C,R,GRAS ) * CFNOG +
     &                                   NORNO( C,R,WETL ) * CFNOW +
     &                                   NORNO( C,R,AGRI ) * CFNOA
     &                                 ) * DENNO     ! sl

                ELSE IF ( TAIR .GT. 265.11111 ) THEN    !  no forest, wet NO

                    CFNOG = EXP( 0.04686 * TAIR  -  14.30579 ) !  grass
                    CFNOA = EXP( 0.05112 * TAIR  -  15.68248 ) !  agriculture

                    EMIS( C,R, NNO ) = ( NORNO( C,R,GRAS ) * CFNOG +
     &                                   NORNO( C,R,AGRI ) * CFNOA
     &                                 ) * DENNO     ! sl

                ELSE IF ( TAIR .GT. 259.8333 ) THEN     ! NO from grass only

                    CFNOG = EXP( 0.04686 * TAIR  -  14.30579 ) !  grass

                    EMIS( C,R, NNO ) = NORNO( C,R,GRAS ) * CFNOG*DENNO  ! sl

                ELSE       !  else tair <= 259.8333

                    EMIS( C,R, NNO ) = 0.0

                END IF          !  5-way conditional on TAIR

111         CONTINUE
122         CONTINUE

        RETURN
94010   FORMAT( A, F10.2, 1X, A, I3, ',', I3 )
        END

