
      SUBROUTINE CLDATEN_CALC(JDATE,JTIME,LAT,LON,CLOUDTYPE,CLDATN)

C***********************************************************************
C  program body starts at line 161
C
C  DESCRIPTION:
C       Computes gridded cloud attenuation factor CLDATN for SMOKE BEIS2
C       Based on the cloud type string in CLDTYPE, reads cloud cover data
C       from Met-Chem Interface Processor for either the KUO (ANTHES-KUO)
C       cloud recalculation, or the McHenry-Kain (KAIN-) cloud recalculation.
C       If CLDTYPE= 'NO CUM PARAM', it is assumed that MM5 ran a simulation
C       without a deep convection parameterization. This is typically done
C       only for urban-scale problems where the grid-resolution is 4km or
C       less. 
C
C       The file containing the information from the Met-Chem Interface Processor
C       is MET_CRO_2D.
C
C  PRECONDITIONS REQUIRED:
C       Must have MET_CRO_2D file from MCIP 
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C
C        Based on structure and simplified approach over the related 
C        routine CLDATEN_OBS. (CLDATEN_OBS computes gridded cloud atten.
C        factors based on observed surface met data...).      
C
C        MET_CRO_2D contains data on cloud fractional coverage, cld base, 
C        cld top, and cloud liquid water content (which are used to calculate
C        attenuated actinic flux in the chemical transport model). Thus, for 
C        this routine we have data for only 1 "combined" cloud, rather than 
C        the multiple layers (up to 3) contained in a typical station 
C        observation file. The "combined" cloud layer may actually represent 
C        multiple layers of cloud, and each layer may represent multiple model 
C        levels. 
C
C        Prototype 1/97 by JM
C
C        Revised and Validated 2/97 by SL
C
C        Revised by M Houyoux - see SCCS update documentation
C
C        Revised for urban-scale applications by John McHenry, 8/22/97.
C
C        FUTURE UPGRADES: Would enable access to the multiple layers
C                         and level-resolved fractional coverages
C                         of the clouds calculated within the McHenry-
C                         Kain (KAIN-) recalculation in MCIP for a
C                         more detailed treatment.... 
C
C***********************************************************************
C  
C Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
C                System
C File Version: @(#)$Id$
C  
C COPYRIGHT (C) 1998, MCNC--North Carolina Supercomputing Center
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
C**********************************************************************

      IMPLICIT NONE

C...........   INCLUDES:

      INCLUDE 'CONST3.EXT'      ! physical and mathematical constants
      INCLUDE 'PARMS3.EXT'      ! I/O API constants
      INCLUDE 'FDESC3.EXT'      ! I/O API file description data structure
      INCLUDE 'IODECL3.EXT'     ! I/O API function declarations
      INCLUDE 'GRDIMS3.EXT'     ! grid-related constants
      INCLUDE 'BIODIMS3.EXT'    ! biogenic-related constants


C...........   ARGUMENTS and their descriptions:

      INTEGER         JDATE           !  current simulation date (YYYYDDD)
      INTEGER         JTIME           !  current simulation time (HHMMSS)
      INTEGER         MDATE           !  current simulation date 
      INTEGER         MTIME           !  current simulation time + 0.5 hr
      REAL            CLDATN( NCOLS, NROWS )  !  cloud attenuation factor
      REAL            LAT  ( NCOLS, NROWS )  ! latitude of cross pt. grid
      REAL            LON  ( NCOLS, NROWS )  ! longitude of cross pt. grid
      CHARACTER*16    CLOUDTYPE
      CHARACTER*5     CLTYPE

C...........   PARAMETERS and their descriptions:
C
      REAL THIN
      PARAMETER(THIN = 200.)  ! thin cld cover <= 200m thick
      REAL            SIGA, SDEC, D100, D60, ROTDAY

      PARAMETER (
     &            SIGA = 279.9348  ,
     &            SDEC = 0.39784984 , ! SIN (23^26'37.8") the declination angle
     &            D100 = 1.0 /100.0,
     &             D60 = 1.0 / 60.0,
     &          ROTDAY = 360.0 / 365.242 )   ! rotation per day


C.......   Variables for zenith angle calculation:

      REAL            DADM, SNDAY, CSDAY, SNDA2, CSDA2
      REAL            AMM, DESIN, DECOS
      REAL            SIGMA, CZ
      REAL            H0, H1
C
C -- empirical factors for attenuation related to height of the
C    cloud base:  position 1: for obscured sky (not used here);
C                 position 2: for cld bases < 1000m
C                 position 3: for cld bases < 3000m
C                 position 4: for cld bases < 7000m
C                 position 5: for cld bases above 7000m
C
      REAL            AA( 5 ), BB( 5 )

      DATA            AA / 0.165, 0.18,  0.365,  0.47,  0.89 /
      DATA            BB / 0.005, 0.02, -0.015, -0.01, -0.04 /

C
C...........   LOCAL VARIABLES and their descriptions:
C
      REAL            CLDT ( NCOLS, NROWS )  ! ave. cloud top in meters
      REAL            CLDB ( NCOLS, NROWS )  ! ave. cloud bottom in meters
      REAL            CFRAC( NCOLS, NROWS )  ! fractional cloud coverage
      REAL            CBOT, CTOP
      REAL            FRAC, CTHK
      REAL            TTYPE 
      REAL            T

      INTEGER         R,C,ROW,COL, I, NCLAY

      LOGICAL         LTHIN
      LOGICAL         FIRSTIME
      DATA            FIRSTIME / .TRUE. /
      SAVE            FIRSTIME

      CHARACTER*16    CLDTVAR   ! Name of input variable
      CHARACTER*16    CLDBVAR   ! Name of input variable
      CHARACTER*16    CFRACVAR  ! Name of input variable
      CHARACTER*16    PNAME
      DATA            PNAME / 'CLDATEN_CALC' /
C
C***********************************************************************
C   begin body of program CLDATEN_CALC

      IF ( FIRSTIME ) THEN

          FIRSTIME = .FALSE.
          IF(CLOUDTYPE(1:5).EQ.'ANTHE') THEN
            CALL M3MESG( 'Using Kuo-package cld data for CLDATEN' )
            CLDTVAR = 'CLDT_KUO'
            CLDBVAR = 'CLDB_KUO'
            CFRACVAR= 'CFRAC_KUO'

          ELSE IF(CLOUDTYPE(1:5).EQ.'KAIN-') THEN
            CALL M3MESG( 'Using KF-package cld data for CLDATEN' )
            CLDTVAR = 'CLDT_KF'
            CLDBVAR = 'CLDB_KF'
            CFRACVAR= 'CFRAC_KF'

          ELSE IF(CLOUDTYPE(1:5).EQ.'NO CU') THEN
            CALL M3MESG( 'Using explicit w/out deep cu cld data 
     &                    for CLDATEN' )
            CLDTVAR = 'CLDT_URB'
            CLDBVAR = 'CLDB_URB'
            CFRACVAR = 'CFRAC_URB'

          ELSE
            CALL M3EXIT( PNAME, JDATE, JTIME,
     &                   'Unsupported or unknown cloud type', 2 )

          END IF
          CLTYPE(1:5) = CLOUDTYPE(1:5)

      END IF                  !  if firstime
C
      IF(CLTYPE(1:5).NE.CLOUDTYPE(1:5)) THEN
        PRINT*,'CLOUD TYPE MISMATCH FOUND...'
            CALL M3EXIT( PNAME, JDATE, JTIME,
     &                   'Unsupported or unknown cloud type', 2 )
      END IF
C
C --- first compute (from spherical geometry) factors needed to 
C     get the cosine of the zenith angle dependent on time of day.
C
C
      DADM  = ROTDAY * FLOAT( MOD( JDATE, 1000 ) - 1 )
      SNDAY = SIN( PI180 * DADM )
      CSDAY = COS( PI180 * DADM )
      SNDA2 = 2.0 * SNDAY * CSDAY
      CSDA2 = CSDAY * CSDAY - SNDAY * SNDAY
      SIGMA  =   SIGA + DADM
     &             + 1.914827 * SNDAY  - 0.079525 * CSDAY
     &             + 0.019938 * SNDA2 - 0.001620 * CSDA2
      DESIN  =  SDEC * SIN( PI180 * SIGMA )
      DECOS  =  SQRT( 1.0 - DESIN * DESIN )

      AMM   =  12.0 + 0.123570 * SNDAY - 0.004289 * CSDAY
     &                + 0.153809 * SNDA2 + 0.060783 * CSDA2
      H0     =  15.0 * ( FLOAT( JTIME / 10000 ) +
     &                     D60 * ( FLOAT( MOD( JTIME / 100, 100 ) )
     &                           + D60 * FLOAT( MOD( JTIME , 100 ) ) )
     &                   - AMM )
C
C...........   Read CFRAC,CLDT,CLDB.
C

C Read CLDT
      MDATE = JDATE
      MTIME = JTIME
C     CALL NEXTIME ( MDATE, MTIME, 3000 )

      IF ( .NOT. READ3( 'SFCCRO', CLDTVAR , ALLAYS3,
     &                    MDATE, MTIME, CLDT         ) ) THEN

            CALL M3EXIT( PNAME, MDATE, MTIME,
     &                  'Could not read CLDT from SFCCRO', 2 )

      END IF  !  if READ3 failed

C Read CLDB

      IF ( .NOT. READ3( 'SFCCRO', CLDBVAR , ALLAYS3,
     &                  MDATE, MTIME, CLDB           ) ) THEN

            CALL M3EXIT( PNAME, MDATE, MTIME,
     &                  'Could not read CLDB from SFCCRO', 2 )

      END IF  !  if READ3 failed

C Read CFRAC

      IF ( .NOT. READ3( 'SFCCRO', CFRACVAR , ALLAYS3,
     &                  MDATE, MTIME, CFRAC           ) ) THEN

            CALL M3EXIT( PNAME, MDATE, MTIME,
     &                  'Could not read CFRAC from SFCCRO', 2 )

      END IF  !  if READ3 failed

C --- compute cld attenuation factors for all cols,rows in domain

      DO 22 ROW = 1,NROWS
        DO 33 COL = 1,NCOLS
          H1  = PI180 * ( H0 + LON( COL, ROW ) ) 
          CZ  = DESIN * SIN (PI180 * LAT(COL,ROW))
     &          + DECOS * COS (PI180 * LAT(COL,ROW)) * COS(H1)
          IF( CFRAC( COL, ROW ) .LE. 0.001) THEN
             CLDATN( COL, ROW ) = 1.
          ELSE
             CZ = 1.0 / MAX( CZ, 0.2 )
             FRAC = CFRAC(COL,ROW)
             CBOT = CLDB(COL,ROW)
             CTOP = CLDT(COL,ROW)
             T = 1.0
             CTHK = CTOP - CBOT
             IF(CTHK.GT.THIN) THEN
               LTHIN = .FALSE.
             ELSE
               LTHIN = .TRUE.
             END IF
C
C -- per consistency with BEIS2 cld attenuation using station data,
C    calculate number of cloud layers with divisions at 1000, 3000,
C    and 7000m by knowing composite column cloud base and top, which
C    have already accounted for multiple layers and model levels...
C
             IF(CBOT.LT.1000.) THEN
                IF(CTOP.LT.1000.) NCLAY = 1
                IF(CTOP.GE.1000..AND.CTOP.LT.3000.) NCLAY = 2
                IF(CTOP.GE.3000.) NCLAY = 3
             ELSE IF(CBOT.GE.1000..AND.CBOT.LT.3000.) THEN
                IF(CTOP.LT.3000.) NCLAY = 1
                IF(CTOP.GE.3000..AND.CTOP.LT.7000.) NCLAY = 2
                IF(CTOP.GE.7000.) NCLAY = 3
             ELSE IF(CBOT.GE.3000..AND.CBOT.LT.7000.) THEN
                IF(CTOP.GE.3000..AND.CTOP.LT.7000.) NCLAY = 1
                IF(CTOP.GE.7000.) NCLAY = 2
             ELSE IF(CBOT.GE.7000.) THEN
                NCLAY = 1
             END IF
C
             DO 11  I = 1,NCLAY
               IF(CBOT.LT.1000.0) THEN
                 TTYPE = AA(I+1) + BB(I+1) * CZ
               ELSE IF(CBOT.GE.1000..AND.CBOT.LT.3000. ) THEN
                 TTYPE = AA(I+2) + BB(I+2) * CZ
               ELSE IF(CBOT.GE.3000..AND.CBOT.LT.7000. ) THEN
                 TTYPE = AA(I+3) + BB(I+3) * CZ
               ELSE
                 TTYPE = AA(5) + BB(5) * CZ
               END IF
               T = T * ( 1.0 - ( 1.0 - TTYPE ) * FRAC)
11           CONTINUE

             IF(LTHIN) THEN
               CLDATN(COL,ROW) = 0.5 + 0.5 * T
             ELSE
               CLDATN(COL,ROW) = T
             END IF
          END IF

33      CONTINUE
22    CONTINUE
C

      DO  55  R = 1, NROWS
      DO  44  C = 1, NCOLS
         T = CLDATN( C,R )
         IF ( T .LT. 0.05 ) THEN
             CLDATN( C,R ) = 0.05
         ELSE IF ( T .GT.  1.0 ) THEN
             CLDATN( C,R ) = 1.0
         END IF
44    CONTINUE
55    CONTINUE

155   CONTINUE

      RETURN
      END
