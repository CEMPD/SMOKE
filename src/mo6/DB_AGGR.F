      SUBROUTINE DB_AGGR(ICY,INERR) 
C 
C     This subroutine outputs a aggregated daily evaporative, 
C     exhaust and refueling emissions to the database output file. 
C 
C     Called by:  PRCSCN 
C 
C     Calls:  None 
C 
C     CHANGES  :  (LAST CHANGE FIRST) 
C 
C  06 Nov 01 AIR Task 27: Updated fuel economy calculations via PMFLECON. 
C  19 Oct 01 AIR Task: Put in code (again) so that all vehicles have MPG output 
C  07 Sep 01 @EPA-djb Bug 333, fixed the summation for ZSTARTS and ZENDS. 
C  27 Aug 01 @EPA-djb Added MEVMYR weighting to ZSTARTS and ZENDS. 
C  07 Aug 01 @EPA-djb Bug330, Changed RFE weighting to TF instead of MEVMYR. 
C  10 Jul 01 @EPA-djb Wish130, Expanded STARTSPERDAY array to include age. 
C  16 Oct 00 @EPA-djb Task W2, Created to calculate aggregated 
C            database output option. 
C 
C  Input on call: 
C 
C    Constants:                None 
C    argument list:            ICY 
C    common blocks:            (See include files below) 
C 
C  Output on return:           None 
C 
C  Local array subscripts: 
C  common blocks:              (See include files below) 
C  Local variable dictionary:  (See declarations below) 
C 
      USE DATABASE, ONLY : DBUNIT,TABCHAR, 
     &                     DBSELPOL,DBSELVEH,DBSELEFT, 
     &                     DBNFILE,DBNRUN,DBNSCEN 
C 
      IMPLICIT NONE 
C 
      INCLUDE  'CEVBMY.I'          ! BMYMPD 
      INCLUDE  'EVAHS1.I'          ! HSRED 
      INCLUDE  'FLAGS3.I'          ! RLFLAG 
      INCLUDE  'IVTYPE.I'          ! VTLDGV, VVGAS, VTGASBUS, VTMC 
      INCLUDE  'MAXIMA.I'          ! MAXPOL, MAXYRS, MAXVEH 
      INCLUDE  'MYRCAL.I'          ! MEVMYR, TF 
      INCLUDE  'RESUL2.I'          ! G_MI_EVP 
      INCLUDE  'RESUL4.I'          ! RUNEF,STARTEF 
      INCLUDE  'STDIST.I'          ! STARTSPERDAY 
      INCLUDE  'VMXCOM.I'          ! VMTMIX 
      INCLUDE  'WEEKLY.I'          ! IWEEK 
C 
      REAL, EXTERNAL :: PMFLECON 
C 
C Input arguments 
C 
      INTEGER, INTENT(IN) :: ICY 
      INTEGER, INTENT(INOUT) :: INERR 
C 
C Local variables 
C 
      INTEGER             :: ETYPE     ! Emission type index 
      INTEGER             :: EVAPT     ! Evaporative emission type index 
      INTEGER             :: IDX       ! Reverse vehicle age index 
      INTEGER             :: IP        ! Pollutant index 
      INTEGER             :: IV        ! Vehicle class index 
      INTEGER             :: JDX       ! Vehicle age index 
      INTEGER             :: MY        ! Model year 
C 
      REAL                :: HSRFAC    ! Fraction of all trips with hot soaks (trip ends) 
      REAL                :: RFE       ! Onroad fuel economy 
      REAL                :: ZENDS     ! SUM of DAILY ENDS 
      REAL                :: ZGMDAY    ! SUM of DAILY GRAMS/DAY 
      REAL                :: ZGMMILE   ! Stores DAILY GRAMS/MILE 
      REAL                :: ZMILES    ! SUM of DAILY MILES 
      REAL                :: ZRFE      ! All model year fuel economy value 
      REAL                :: ZSTARTS   ! SUM of DAILY STARTS 
C 
      CHARACTER(22)       ::  DBHEADER !Initial part of database output record. 
C 
C     Begin emission type and pollutant looping 
C 
      DO IP=1,MAXPOL 
        IF(DBSELPOL(IP).EQ.1) CYCLE 
        DO IV=1,MAXVEH 
          IF(DBSELVEH(IV).EQ.1) CYCLE 
C 
          IF (IV == VTLDGV) THEN 
            HSRFAC = HSRED(1) 
          ELSE 
            HSRFAC = HSRED(2) 
          ENDIF 
C 
          ZSTARTS = 0.0 
          ZENDS   = 0.0 
          DO JDX=1,25 
            ZSTARTS = ZSTARTS + 
     &                STARTSPERDAY(JDX,IV,IWEEK)*MEVMYR(JDX,IV) 
            ZENDS   = ZENDS + 
     &                STARTSPERDAY(JDX,IV,IWEEK)*HSRFAC*MEVMYR(JDX,IV) 
          END DO 
C 
C     Construct the header section of the record 
C 
      WRITE(DBHEADER,110) DBNFILE,TABCHAR,DBNRUN,TABCHAR, 
     &      DBNSCEN,TABCHAR,ICY,TABCHAR,IP,TABCHAR,IV,TABCHAR 
  110 FORMAT(I3,A1,I3,A1,I3,A1,I4,A1,I1,A1,I2,A1) 
C 
C     Sum up emissions, miles and fuel economy 
C 
          ZGMMILE = 0.0 
          ZGMDAY  = 0.0 
          DO ETYPE=1,8 
            IF(DBSELEFT(ETYPE).EQ.1) CYCLE 
            IF(ETYPE.GE.3 .AND. IP.GT.1) CYCLE 
            IF(ETYPE.EQ.8.AND.RLFLAG.EQ.5) CYCLE 
            IF(ETYPE.LE.3) EVAPT = 1 
            IF(ETYPE.EQ.4) EVAPT = 2 
            IF(ETYPE.EQ.5) EVAPT = 5 
            IF(ETYPE.EQ.6) EVAPT = 3 
            IF(ETYPE.EQ.7) EVAPT = 4 
            IF(ETYPE.GE.8) EVAPT = 6 
C           Add appropriate emissions to output value 
            IF(ETYPE.EQ.1) THEN 
              ZGMMILE = ZGMMILE + RUNEF(IP,IV) 
            ELSEIF(ETYPE.EQ.2) THEN 
              ZGMMILE = ZGMMILE + STARTEF(IP,IV) 
            ELSE 
              ZGMMILE = ZGMMILE + G_MI_EVP(EVAPT,IV) 
            ENDIF 
            ZMILES  = 0.0 
            ZRFE    = 0.0 
C 
            DO JDX=1,MAXYRS 
C 
C             Daily miles traveled 
              ZMILES = ZMILES + BMYMPD(JDX,IV)*MEVMYR(JDX,IV) 
C 
C             Fleet average fuel economy (miles per gallon) 
              MY=ICY-JDX+1 
              RFE = PMFLECON(MY,IV,INERR) ! Compute fuel economy for all vehicle classes 
C             Convert mi/gal to gal/mi, weight and sum. (harmonic mean) 
              IF(RFE.GT.0.0) THEN 
                 IDX=(MAXYRS+1)-JDX 
                 ZRFE = ZRFE + (1.0/RFE)*TF(IDX,IV) 
              ENDIF 
C 
            END DO  !JDX 
            IF(ZRFE.LE.0.0) THEN 
              ZRFE = 0.0 
            ELSE 
              ZRFE = 1.0 / ZRFE   !Convert gal/mi units to miles per gallon 
            ENDIF 
          END DO    !ETYPE 
C         Calculate grams per day per vehicle emissions 
          ZGMDAY  = ZGMMILE * ZMILES 
          IF(VMTMIX(IV).LE.0.0) THEN 
            ZSTARTS=0.0 
            ZENDS=0.0 
          ENDIF 
C 
C         Write the aggregated emissions to the database output file. 
          WRITE (DBUNIT(1),220) DBHEADER, 
     &       ZGMMILE, TABCHAR, ZGMDAY, TABCHAR, ZSTARTS, TABCHAR,  
     &       ZENDS, TABCHAR, ZMILES, TABCHAR, ZRFE, TABCHAR, 
     &       VMTMIX(IV) 
  220     FORMAT(A22,F9.4,A1,F10.3,A1,F8.4,A1,F8.4,A1,F9.4,A1, 
     &       F7.2,A1,F9.6) 
C 
        END DO  !IV 
      END DO    !IP 
C 
      RETURN 
      END SUBROUTINE DB_AGGR 
