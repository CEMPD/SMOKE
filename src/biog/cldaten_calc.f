
      SUBROUTINE CLDATEN_CALC( JDATE, JTIME, NX, NY, LAT ,LON ,
     &                         CLOUDTYPE ,CLDATN )

C***********************************************************************
C  program body starts at line 138
C
C  DESCRIPTION:
C       Computes gridded cloud attenuation factor CLDATN for SMOKE BEIS2
C       Based on the cloud type string in CLDTYPE, reads cloud cover data
C       from Met-Chem Interface Processor for either the KUO (ANTHES-KUO)
C       cloud recalculation, or the McHenry-Kain (KAIN-) cloud recalculation.
C       If CLDTYPE= 'No deep convection param', it is assumed that MM5 ran 
C       a simulation without a deep convection parameterization. This is 
C       typically done  only for urban-scale problems where the 
C       grid-resolution is 4km or less. 
C
C  PRECONDITIONS REQUIRED:
C
C       NOT TESTED FULLY YET!  Use BG_CLOUD_TYPE set to 1 or 5 at this
C       time. 
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C
C      11/99: by Jeff Vukovich taken from v4.4 of SMOKE prototype
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

C...........   INCLUDES:

        INCLUDE 'PARMS3.EXT'      ! I/O API constants
        INCLUDE 'FDESC3.EXT'      ! I/O API file description data structure
        INCLUDE 'IODECL3.EXT'     ! I/O API function declarations
        INCLUDE 'CONST3.EXT'      ! I/O API constants
        INCLUDE 'EMCNST3.EXT'     ! emissions constants
        INCLUDE 'BIODIMS3.EXT'    ! biogenic-related constants

C...........   ARGUMENTS and their descriptions:

        INTEGER, INTENT (IN)  :: JDATE   ! current simulation date (YYYYDDD)
        INTEGER, INTENT (IN)  :: JTIME   ! current simulation time (HHMMSS)
        INTEGER, INTENT (IN)  :: NX      ! no. cols
        INTEGER, INTENT (IN)  :: NY      ! no. rows
        REAL, INTENT (IN)  ::  LAT  ( NX, NY )  !  lat (deg) -90 <= LAT <= 90
        REAL, INTENT (IN)  ::  LON  ( NX, NY )  !  lon (deg) -180 <= LON <= 180
        CHARACTER(16), INTENT(IN) :: CLOUDTYPE

        REAL, INTENT (OUT)  ::  CLDATN( NX, NY ) !  cloud attenuation factor 
      
C...........   PARAMETERS and their descriptions:
C
      INTEGER   IOS                    ! iostat variable

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
      REAL, ALLOCATABLE, SAVE  :: CLDT ( :, : )  ! ave. cloud top in meters
      REAL, ALLOCATABLE, SAVE  :: CLDB ( :, : )  ! ave. cloud bottom in meters
      REAL, ALLOCATABLE, SAVE  :: CFRAC( :, : )  ! fractional cloud coverage
      REAL            CBOT, CTOP
      REAL            FRAC, CTHK
      REAL            TTYPE 
      REAL            T

      INTEGER         R,C,ROW,COL, I, NCLAY

      LOGICAL         LTHIN
      LOGICAL         FIRSTIME
      DATA            FIRSTIME / .TRUE. /
      SAVE            FIRSTIME

      CHARACTER(5)    CLTYPE
      CHARACTER(16)   CLDTVAR   ! Name of input variable
      CHARACTER(16)   CLDBVAR   ! Name of input variable
      CHARACTER(16)   CFRACVAR  ! Name of input variable
      SAVE            CLTYPE, CLDTVAR, CLDBVAR, CFRACVAR

      CHARACTER(16)   PNAME
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

C
C........... Allocate memory for cloud variables
         ALLOCATE( CLDT( NX, NY ), STAT=IOS )
         CALL CHECKMEM( IOS, 'CLDT', PNAME )

         ALLOCATE( CLDB( NX, NY ), STAT=IOS )
         CALL CHECKMEM( IOS, 'CLDB', PNAME )

         ALLOCATE( CFRAC( NX, NY ), STAT=IOS )
         CALL CHECKMEM( IOS, 'CFRAC', PNAME )

      END IF                  !  if firstime
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

      IF ( .NOT. READ3( 'MET_CRO_3D', CLDTVAR , ALLAYS3,
     &                    JDATE, JTIME, CLDT         ) ) THEN

            CALL M3EXIT( PNAME, JDATE, JTIME,
     &                  'Could not read CLDT from SFCCRO', 2 )

      END IF  !  if READ3 failed

C Read CLDB

      IF ( .NOT. READ3( 'MET_CRO_3D', CLDBVAR , ALLAYS3,
     &                  JDATE, JTIME, CLDB           ) ) THEN

            CALL M3EXIT( PNAME, JDATE, JTIME,
     &                  'Could not read CLDB from SFCCRO', 2 )

      END IF  !  if READ3 failed

C Read CFRAC

      IF ( .NOT. READ3( 'MET_CRO_3D', CFRACVAR , ALLAYS3,
     &                  JDATE, JTIME, CFRAC           ) ) THEN

            CALL M3EXIT( PNAME, JDATE, JTIME,
     &                  'Could not read CFRAC from SFCCRO', 2 )

      END IF  !  if READ3 failed

C --- compute cld attenuation factors for all cols,rows in domain

      DO  ROW = 1,NY
        DO  COL = 1,NX
          H1  = PI180 * ( H0 + LON( COL, ROW ) ) 
          CZ  = DESIN * SIN ( PI180 * LAT( COL, ROW ))
     &          + DECOS * COS ( PI180 * LAT( COL, ROW )) * COS( H1 )
          IF( CFRAC( COL, ROW ) .LE. 0.001) THEN
             CLDATN( COL, ROW ) = 1.
          ELSE
             CZ = 1.0 / MAX( CZ, 0.2 )
             FRAC = CFRAC( COL, ROW )
             CBOT =  CLDB( COL, ROW )
             CTOP =  CLDT( COL, ROW )
             T = 1.0
             CTHK = CTOP - CBOT
             IF( CTHK .GT. THIN ) THEN
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
             IF( CBOT .LT. 1000. ) THEN
                IF( CTOP .LT. 1000. ) NCLAY = 1
                IF( CTOP .GE. 1000. .AND. CTOP .LT. 3000. ) NCLAY = 2
                IF( CTOP .GE. 3000. ) NCLAY = 3
             ELSE IF( CBOT .GE. 1000. .AND. CBOT .LT. 3000. ) THEN
                IF( CTOP .LT. 3000. ) NCLAY = 1
                IF( CTOP .GE. 3000. .AND. CTOP .LT. 7000. ) NCLAY = 2
                IF( CTOP .GE. 7000. ) NCLAY = 3
             ELSE IF( CBOT .GE. 3000. .AND. CBOT .LT. 7000. ) THEN
                IF( CTOP .GE. 3000. .AND. CTOP .LT. 7000. ) NCLAY = 1
                IF( CTOP .GE. 7000. ) NCLAY = 2
             ELSE IF( CBOT .GE. 7000. ) THEN
                NCLAY = 1
             END IF
C
             DO 11  I = 1,NCLAY
               IF( CBOT .LT. 1000.0) THEN
                 TTYPE = AA(I+1) + BB(I+1) * CZ
               ELSE IF( CBOT .GE. 1000. .AND. CBOT .LT. 3000. ) THEN
                 TTYPE = AA(I+2) + BB(I+2) * CZ
               ELSE IF( CBOT .GE. 3000. .AND. CBOT .LT. 7000. ) THEN
                 TTYPE = AA(I+3) + BB(I+3) * CZ
               ELSE
                 TTYPE = AA(5) + BB(5) * CZ
               END IF
               T = T * ( 1.0 - ( 1.0 - TTYPE ) * FRAC)
11           CONTINUE

             IF( LTHIN ) THEN
               CLDATN(COL,ROW) = 0.5 + 0.5 * T
             ELSE
               CLDATN(COL,ROW) = T
             END IF
          END IF

        ENDDO
      ENDDO
C

      DO  R = 1, NY
      DO  C = 1, NX
         T = CLDATN( C,R )
         IF ( T .LT. 0.05 ) THEN
             CLDATN( C,R ) = 0.05
         ELSE IF ( T .GT.  1.0 ) THEN
             CLDATN( C,R ) = 1.0
         END IF
      ENDDO
      ENDDO

155   CONTINUE

      RETURN
      END
