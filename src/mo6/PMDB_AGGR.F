      SUBROUTINE PMDB_AGGR(ICY) 
C 
C     This subroutine outputs a aggregated daily PM emissions 
C     to the database output file. 
C 
C     CHANGES  :  (LAST CHANGE FIRST) 
C 
C  04 Sep 01 AIR Task 15-22: Applied Bug Fix 331. 
C  06 Mar 01 AIR Task 11: New MOBILE6 routine. 
C 
C     Called by:  PMOUTDB 
C 
C     Calls:  None 
C 
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
     &                     DBSELPART,DBSELVEH,DBSELEFT, 
     &                     DBNFILE,DBNRUN,DBNSCEN 
C 
      IMPLICIT NONE 
C 
      INCLUDE  'CEVBMY.I'          ! BMYMPD 
      INCLUDE  'EVAHS1.I'          ! HSRED 
      INCLUDE  'FLAGS3.I'          ! RLFLAG 
      INCLUDE  'IVTYPE.I'          ! VTLDGV, VVGAS, VTGASBUS, VTMC 
      INCLUDE  'MAXIMA.I'          ! MAXYRS, MAXVEH 
      INCLUDE  'MYRCAL.I'          ! MEVMYR 
      INCLUDE  'PART1.I'           ! MAXIPPM 
      INCLUDE  'PART3.I'           ! EFPM 
      INCLUDE  'RESUL2.I'          ! G_MI_EVP 
      INCLUDE  'RESUL4.I'          ! RUNEF,STARTEF 
      INCLUDE  'STDIST.I'          ! STARTSPERDAY 
      INCLUDE  'VMXCOM.I'          ! VMTMIX 
      INCLUDE  'WEEKLY.I'          ! IWEEK 
C 
C Input arguments 
C 
      INTEGER, INTENT(IN) :: ICY 
C 
C  Declare external functions. 
C 
      REAL, EXTERNAL :: PMFLECON 
C 
C Local variables 
C 
      INTEGER             :: ETYPE     ! Emission type index 
      INTEGER             :: IDX       ! Year index 
      INTEGER             :: IP        ! MOBILE6 extended pollutant identifier 
      INTEGER             :: IPPM      ! PM pollutant identifier 
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
      CHARACTER(23)       :: DBHEADER  !Initial part of database output record. 
C 
C     Begin PM emission type and pollutant looping, skip if not selected. 
C 
      DO IPPM=1,MAXIPPM 
        IF(DBSELPART(IPPM).EQ.1) CYCLE 
C 
        IP=IPPM+6 
C 
        DO IV=1,MAXVEH 
          IF(DBSELVEH(IV).EQ.1) CYCLE 
C 
C     Skip if GASPM for DIESELS, OCARBON/ECARBON for GAS 
C 
          IF(IPPM.EQ.5.AND.VVDSL(IV).EQ.1 .OR. 
     *      IPPM.GE.2.AND.IPPM.LE.3.AND.VVGASMC(IV).EQ.1) CYCLE 
C 
C     Set up 
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
  110 FORMAT(I3,A1,I3,A1,I3,A1,I4,A1,I2,A1,I2,A1) 
C 
C     Sum up emissions, miles and fuel economy 
C 
          ZGMMILE = 0.0 
          ZGMDAY  = 0.0 
C 
          DO ETYPE=1,10 
            IF(DBSELEFT(ETYPE).EQ.1) CYCLE 
C           Add appropriate emissions to output value 
            ZGMMILE = EFPM(IPPM,IV) 
C 
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
              RFE = PMFLECON(MY,IV) 
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
  220     FORMAT(A23,F9.4,A1,F10.3,A1,F8.4,A1,F8.4,A1,F9.4,A1, 
     &       F7.2,A1,F9.6) 
C 
        END DO  !IV 
      END DO    !IPPM 
C 
      RETURN 
C 
      END SUBROUTINE PMDB_AGGR 
