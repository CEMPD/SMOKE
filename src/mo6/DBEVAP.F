      SUBROUTINE DBEVAP(IEVAP,IV,JDX,GRAMS_MI,GRAMS_HR,IROAD,MY,INERR) 
C 
C     This subroutine outputs the evaporative emissions, daily or hourly, 
C       to the Database. 
C 
C 
C     Called by:  EVP_HRLY 
C 
C     Calls:  DB_DAILY 
C 
C     CHANGES  :  (LAST CHANGE FIRST) 
C 
C  06 Nov 01 AIR Task 27: Updated fuel economy calculations via PMFLECON. 
C            Added INERR processing. 
C  10 Jul 01 @EPA-djb Wish130, Expanded STARTSPERDAY array to include age. 
C  13 Feb 01 @EPA-elg increased the size of FVMT arrays to include vehicle class 
C  24 Oct 00 @EPA-BG  New labels for database age, model year and hour selection 
C 19 Oct 00 @EPA-djb Task W2, Added DBAGGR condition on output. 
C 20 Oct 00 AIR Task 03: Added MYR to daily and hourly output 
C 25 Jul 00 @EPA-bag bug52  No refueling losses for motorcycles 
C 10 Jul 00 @EPA-ddj Bug187 Removed 2nd parameter from list of 8 in 
C           call to DB_DAILY. 
C 06 Jul 00 @EPA-bag bug52  No refueling losses for motorcycles 
C 23 Jun 00 @EPA-bag bug67  Added diurnal temperature, DIRNL_TEMP, for output 
C 23 Jun 00 @EPA-bag bug52  No running losses for motorcycles 
C 30 May 00 @EPA-elg Replaced BDSPEED1.I with SPEED9.I 
C  05 MAY 00 @EPA -HMC; CHANGED OUTPUT FORMATS, AND ADDED DAILY OUTPUT 
C 
C    common blocks: 
C 
C    /SPEED9/   FVMT,HVMT 
C    /CEVBMY/   BMYMPD 
C    /IVTYPE/   VVSTARTS 
C    /MYRCAL/   MEVMYR 
C    /STDIST/   STARTDIST, STARTSPERDAY 
C    /TEMPS/    TEMHRLY 
C    /VMXCOM/   VCOUNT 
C    /WEEKLY/   IWEEK 
C 
      USE DATABASE, ONLY : DBSELPOL,DBFLAG,DBSELVEH,DBSELEFT,DBSELFAC, 
     &                     DBNFILE,DBNRUN,DBNSCEN,DBDAILY,DBAGGR, 
     &                     DBUNIT,TABCHAR,DBSELHR,DBSELAGE 
      USE DITEMDAT, ONLY : DIRNL_TEMP  ! Source for Diurnal Temperatures 
C 
      IMPLICIT NONE 
C 
      INCLUDE  'SPEED9.I'          ! FVMT,HVMT 
      INCLUDE  'CEVBMY.I'          ! BMYMPD 
      INCLUDE  'EVAHS1.I'          ! HSRED, HS_HRLY_TRP_FRCN 
      INCLUDE  'IVTYPE.I'          ! VVSTARTS 
      INCLUDE  'MYRCAL.I'          ! MEVMYR 
      INCLUDE  'STDIST.I'          ! STARTDIST, STARTSPERDAY 
      INCLUDE  'TEMPS.I'           ! TEMHRLY 
      INCLUDE  'VMXCOM.I'          ! DBVCOUNT 
      INCLUDE  'WEEKLY.I'          ! IWEEK 
C 
C Externals 
C 
      REAL, EXTERNAL :: PMFLECON 
C 
C Input arguments 
C 
      REAL, INTENT(IN)    :: GRAMS_MI(24) 
      REAL, INTENT(IN)    :: GRAMS_HR(24) 
      INTEGER, INTENT(IN) :: IV 
      INTEGER, INTENT(IN) :: JDX 
      INTEGER, INTENT(IN) :: IROAD 
      INTEGER, INTENT(IN) :: IEVAP 
      INTEGER, INTENT(IN) :: MY 
      INTEGER, INTENT(INOUT) :: INERR 
C 
C Local variables 
C 
      INTEGER             :: IH       ! HOUR INDEX 
      INTEGER, PARAMETER  :: IMOTORCYCLE = 24 
      INTEGER, PARAMETER  :: IGASBUS = 25 
      INTEGER, SAVE       :: UPPER_LMT=2 
      INTEGER, SAVE       :: LOWER_LMT=1 
      LOGICAL             :: IVOK     ! FLAG FOR VEHICLES TO PROCESS 
      LOGICAL             :: AGE_OK   ! FLAG FOR AGES/MODEL YEARS TO PROCESS 
      LOGICAL             :: NO_DBOUT ! FLAG SET IF NO DATABASE OUTPUT 
                                      ! REQUESTED. 
      REAL                :: RFE      ! MPG 
      REAL                :: FACVMT 
      REAL                :: HSRFAC 
      REAL, DIMENSION(24) :: DTEMP 
C 
      CHARACTER(LEN=25) :: DBHEADER 
C 
C IF THIS OUTPUT IS NOT SELECTED IN DATABASE OPTIONS DO NOTHING 
C 
      IVOK=( 
     &      IV<=13 .OR. IV==IGASBUS .OR. 
     &      (IV==IMOTORCYCLE.AND.(IEVAP/=6.AND.IEVAP/=8))  ! no (database) running AND 
     &     )                                               ! refueling losses for motorcycles 
C 
      AGE_OK = ( 
     &          DBSELAGE(LOWER_LMT) <= (JDX-1) .AND. 
     &            (JDX-1)<= DBSELAGE(UPPER_LMT) 
     &         ) 
C 
      NO_DBOUT=( 
     &          DBSELPOL(1) /=2 .OR. 
     &          DBFLAG /= 2 .OR. 
     &          DBSELVEH(IV) /= 2 .OR. 
     &          DBSELEFT(IEVAP) /=2 .OR. 
     &          DBSELFAC(IROAD) /=2 .OR. 
     &          .NOT.IVOK  .OR. 
     &          .NOT.AGE_OK 
     &         ) 
C 
      IF (IV == VTLDGV) THEN 
        HSRFAC = HSRED(1) 
      ELSE 
        HSRFAC = HSRED(2) 
      ENDIF 
C 
      IF (NO_DBOUT) THEN 
        RETURN 
      ELSE 
C 
        IF(IV.EQ.VTMC) THEN 
          RFE=0.0 
        ELSE 
          RFE=PMFLECON(MY,IV,INERR) 
        END IF 
C 
C       PREPARE HEADER INFORMATION FOR EITHER FORM OF OUTPUT 
C 
        WRITE(DBHEADER,210) DBNFILE,TABCHAR,DBNRUN,TABCHAR, 
     &     DBNSCEN,TABCHAR,1,TABCHAR,IV, 
     &     TABCHAR,IEVAP,TABCHAR,IROAD,TABCHAR, 
     &     JDX-1,TABCHAR 
C 
        IF (DBDAILY) THEN 
C 
           CALL DB_DAILY(GRAMS_MI,HSRFAC,DBHEADER,RFE,IROAD,IV,JDX,MY) 
C 
        ELSEIF (.NOT.DBAGGR) THEN 
C 
C          WRITE HOURLY DATABASE OUTPUT 
C 
           DO IH=DBSELHR(LOWER_LMT),DBSELHR(UPPER_LMT) 
               SELECT CASE (IROAD) 
                   CASE (1:4) 
                       FACVMT=FVMT(IROAD,IH,IV) 
                   CASE (5) 
                       FACVMT=1.0 
               END SELECT 
               DTEMP(1:19)=DIRNL_TEMP(6:24) 
               DTEMP(20:24)=DIRNL_TEMP(1:5) 
               WRITE(DBUNIT(1),215) DBHEADER,IH, 
     &             TABCHAR, GRAMS_MI(IH), 
     &             TABCHAR, GRAMS_HR(IH), 
     &             TABCHAR, STARTSPERDAY(JDX,IV,IWEEK)* 
     &                      STARTDIST(IH,IWEEK), 
     &             TABCHAR, STARTSPERDAY(JDX,IV,IWEEK)* 
     &                      HSRFAC*HS_HRLY_TRP_FRCN(IH,IWEEK), 
     &             TABCHAR, BMYMPD(JDX,IV) * HVMT(IH), 
     &             TABCHAR, RFE,TABCHAR,HVMT(IH),TABCHAR,FACVMT, 
     &             TABCHAR, MEVMYR(JDX,IV),TABCHAR,DBVCOUNT(IV), 
     &             TABCHAR, TEMHRLY(IH),TABCHAR,DTEMP(IH), 
     &             TABCHAR, MY 
           END DO 
        ENDIF 
C 
      ENDIF 
  210 FORMAT(I3,A1,I3,A1,I3,A1,I1,A1,I2,A1,I2,A1,I1,A1,I2,A1) 
  215 FORMAT(A25,I2,A1,F9.4,A1,F9.4,A1,F8.4,A1,F8.4,A1,F9.4,A1, 
     &       F7.2,A1,F9.4,A1,F9.4,A1,F8.4,A1,F10.4,A1,F6.1,A1,F6.1, 
     &       A1,I4) 
C 
      RETURN 
      END SUBROUTINE 
