      SUBROUTINE CO2DB_AGGR(ICY,INERR) 
C 
C     This subroutine outputs aggregated daily CO2 emissions 
C     to the database output file. The output units are g/mi. 
C 
C     CHANGES  :  (LAST CHANGE FIRST) 
C 
C  04 Dec 01 AIR Task 27: New routine. 
C 
C     Called by:  CO2OUTDB 
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
     &                     DBSELPOL,DBSELVEH,DBSELEFT, 
     &                     DBNFILE,DBNRUN,DBNSCEN 
C 
      IMPLICIT NONE 
C 
      INCLUDE  'RESUL1.I'          ! EFFTP 
      INCLUDE  'CEVBMY.I'          ! BMYMPD 
      INCLUDE  'EVAHS1.I'          ! HSRED 
      INCLUDE  'FLAGS3.I'          ! RLFLAG 
      INCLUDE  'IVTYPE.I'          ! VTLDGV, VVGAS, VTGASBUS, VTMC, VVDSL 
      INCLUDE  'MAXIMA.I'          ! MAXYRS, MAXVEH 
      INCLUDE  'MYRCAL.I'          ! MEVMYR 
      INCLUDE  'PART1.I'           ! MAXIPPM 
      INCLUDE  'STDIST.I'          ! STARTSPERDAY 
      INCLUDE  'VMXCOM.I'          ! VMTMIX 
      INCLUDE  'WEEKLY.I'          ! IWEEK 
C 
C Input arguments 
C 
      INTEGER, INTENT(IN) :: ICY 
      INTEGER, INTENT(INOUT) :: INERR 
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
C     Return if CO2 is not selected. 
C 
      IP=4 
C    
        IF(DBSELPOL(IP).EQ.1) RETURN 
C 
C  CO2 is exhaust only (running + start)  
C 
          ETYPE=1 
C 
          IF(DBSELEFT(ETYPE).EQ.1) RETURN 
C 
C     Loop for each vehicle class, print if selected. 
C 
        DO IV=1,MAXVEH 
          IF(DBSELVEH(IV).EQ.1) CYCLE 
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
C           Add appropriate emissions to output value.  
C 
            ZGMMILE = EFFTP(IP,IV) 
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
              RFE = PMFLECON(MY,IV,INERR) 
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
  220     FORMAT(A23,F9.1,A1,F10.1,A1,F8.4,A1,F8.4,A1,F9.4,A1, 
     &       F7.2,A1,F9.6) 
C 
          END DO    !IV 
C 
      RETURN 
C 
      END SUBROUTINE CO2DB_AGGR 
