      SUBROUTINE ATCALX(ICY,INERR) 
C 
C  ATCALX computes the exhaust and evap Air Toxic emission factors. 
C 
      USE DATABASE, ONLY : DBSELTOX 
C 
      IMPLICIT NONE 
C 
      INCLUDE 'ATOX1.I' 
      INCLUDE 'ATOX2.I' 
      INCLUDE 'ATOX3.I' 
      INCLUDE 'ATOX4.I' 
      INCLUDE 'ATOX5.I'   ! 
      INCLUDE 'CITRV1.I'  ! RVPICY 
      INCLUDE 'EVAPACT.I' ! SAWFLAG 
      INCLUDE 'FLAGS2.I'  ! EVIMFLAG 
      INCLUDE 'FLAGS3.I'  ! CAAFLG, CAAMY 
      INCLUDE 'FLAGS5.I'  ! SPD_FLG,SOAK_FLAG 
      INCLUDE 'IVTYPE.I'  ! MAXVEH,VTMC,VVLDGAS 
      INCLUDE 'MAXIMA.I'  ! MAXIH, MAXYRS 
      INCLUDE 'MYRCAL.I'  ! TF 
      INCLUDE 'SPEED9.I'  ! FVMT, HVMT 
      INCLUDE 'VMXCOM.I'  ! VMTMIX 
C 
C  Declare parameter list. 
C 
      INTEGER ICY 
      INTEGER INERR 
C 
C  Declare local variables/arrays. 
C 
      INTEGER ACTUAL_MY 
      INTEGER I 
      INTEGER IAT 
      INTEGER IB 
      INTEGER ICHK 
      INTEGER IDX 
      INTEGER IE 
      INTEGER IEVAP 
      INTEGER IH 
      INTEGER INH 
      INTEGER IOXY 
      INTEGER IROAD 
      INTEGER IV 
      INTEGER IVGAS 
      INTEGER IVT 
      INTEGER IXV 
      INTEGER JDX 
      INTEGER MY 
      REAL    ATEVAP(5,2) 
      REAL    CAT(3) 
      REAL    MTBEVOL 
      REAL    P3W 
      REAL    POX 
      REAL    RVP 
      REAL    SUMD 
      REAL    SUMIH 
      REAL    SUMINH 
      REAL    SUMJDX 
      REAL    SUMN 
      REAL    TCAT 
 
CDFK    OPEN(99,FILE='ATJUNK') 
C 
C  Check that all required air toxic data has been entered. 
C 
      DO ICHK=1,6 
        IF(TOXCHK(ICHK).EQ.1) THEN 
          CALL QUITER(0.0,ICHK,116,INERR) 
        END IF 
      END DO 
C 
      IF(INERR.GT.0) RETURN 
C 
C ********************** Set Up Exhuast ******************************* 
C 
C  Initialize the Air Toxic Exhaust EF arrays 
C 
      DBAT = 0.0 
      EFAT = 0.0 
       VAT = 0.0 
C 
C  Determine the market share weighted oxygen content. The oxygen content 
C  must be in wt%. 
C 
      SUMN=0. 
      SUMD=0. 
      DO I=1,4                                     !I=1:MTBE, 2:ETBE, 3:ETHANOL, 4:TAME 
        SUMN=SUMN+ATOXY(1,I)*ATOXY(2,I)            !Sum of OXY*MARKETSHARE, in wt% 
        SUMD=SUMD+ATOXY(2,I)                       !Sum of MARKETSHARE 
      END DO 
C 
      ATWTOXY=0. 
      IF(SUMD.GT.0.) ATWTOXY=SUMN/SUMD   !Market share weighted oxygen content, in wt% 
C 
C  Set up the Complex Model target fuel parameters. Note that the Sulfur 
C  value is set to 330 ppm, NOT the in-use fuel value from MOBILE6. 
C  This lets MOBILE6 adjust for sulfur, not the outdated routines in 
C  the Complex Model. The RVP is the MOBILE6 weathered RVP. 
C 
      ATTRGT( 2)=330.          !Sulfur level, held at 330 ppm for Complex Model  
      ATTRGT( 3)=RVPICY        !Weathered RVP 
      ATTRGT( 4)=E200 
      ATTRGT( 5)=E300 
      ATTRGT( 6)=AROM 
      ATTRGT( 7)=OLEF 
      ATTRGT( 8)=BENZ 
C 
C  Call the soak function, SOAKFCN, to set up the array, SOAKFUN, for HC (IP=1). 
C  These values are used to determine hourly soak time correction factors for use 
C  in ATHRLOOP. They are dependent only on soak time, pollutant, and type of catalyst. 
C 
      IF(SOAK_FLAG .EQ. 1) CALL SOAKFCN(1)   ! (IP=1) 
C 
C ********************** Set Up Evap ******************************* 
C 
C  Initialize the Air Toxic Evaporative Conversion Fractions 
C  (Currently, the conversion fractions are independent of vehicle 
C   type and model year) 
C 
      ATEVAP = 0. 
C 
C  Compute the evap RVP. Limit is set to 8.7 psi to avoid 
C  extrapolation. 
C 
      RVP=MIN(8.7,RVPICY) 
C 
C  Compute the Benzene evap fractions. Oxygen content (ATWTOXY) is in wt% 
C  However, the BENZ values are in vol%. 
C 
      IF(BENZ.GT.0.) THEN 
      ATEVAP(1,1) = (- 0.03420  * ATWTOXY  
     *               - 0.080274 * RVP  
     *               + 1.4448)  * BENZ / 100. !Hot Soak 
      ATEVAP(2,1) = (- 0.02895  * ATWTOXY  
     *               - 0.080274 * RVP 
     *               + 1.3758)  * BENZ / 100. !Diurnal 
      ATEVAP(3,1) = (- 0.03420  * ATWTOXY  
     *               - 0.080274 * RVP  
     *               + 1.4448)  * BENZ / 100. !Running 
      ATEVAP(4,1) = (- 0.02895  * ATWTOXY  
     *               - 0.080274 * RVP 
     *               + 1.3758)  * BENZ / 100. !Resting 
      ATEVAP(5,1) = (- 0.02955  * ATWTOXY 
     *               - 0.081507 * RVP 
     *               + 1.3972)  * BENZ / 100. !Refueling 
      END IF 
C 
C  Compute the MTBE evap fractions. Oxygen content (ATWTOXY) is in wt%. 
C  However, the MTBE values must be in vol% 
C 
      IF(ATOXY(1,1).GT.0.) THEN 
      MTBEVOL=ATOXY(1,1)/VOLTOWT(1) 
      ATEVAP(1,2) = (24.2050 - 1.7460 * RVP) * MTBEVOL / 1000. !Hot Soak     
      ATEVAP(2,2) = (22.1980 - 1.7460 * RVP) * MTBEVOL / 1000. !Diurnal     
      ATEVAP(3,2) = (17.8538 - 1.6622 * RVP) * MTBEVOL / 1000. !Running     
      ATEVAP(4,2) = (22.1980 - 1.7460 * RVP) * MTBEVOL / 1000. !Resting     
      ATEVAP(5,2) = 1.743 
     *            * (-0.02995 * ATWTOXY - 0.081507 * RVP + 1.3972) 
     *            * MTBEVOL / 100.                            !Refueling 
      END IF 
C 
C  Compute the hourly refueling losses. 
C 
      CALL REFUEL(ICY,INERR) 
C 
C  Initialize the diurnal fractions 
C 
      CALL DIRNLINI() 
C 
C  Initialize gas vehicle index, IVGAS, used for evaporative 
C  emissions calculations 
C 
      IVGAS = 0 
C 
C   To avoid problems with EVAPCALC call below, trap case where user has entered 
C   Evap I/M programs but has also specified NO EVP SAWTOOTH. Reset SAWFLAG to 1 
C   (sawtooth enabled) whenever this happens. 
C 
      IF(EVIMFLAG.GT.1 .AND. SAWFLAG.EQ.0) THEN 
        CALL QUITER(0.,0,563,INERR) 
        SAWFLAG = 1 
      END IF 
C 
C ********************** Check and Set Up HAPS ********************** 
C 
      IF(ATHAPS) CALL ATHCHK(ICY,INERR) 
C 
C ********************** Begin Vehicle Loop ************************* 
C 
      DO IV=1,MAXVEH 
C 
C  Determine the evap vehicle index 
C 
      IVGAS = IVGAS + VVGASMC(IV) 
C 
C  Skip if VMTMIX is 0. 
C 
      IF (VMTMIX(IV) .EQ. 0.) CYCLE 
C 
C  Calculate the speed correction factor for 2.5-65 m/h for HC (IP=1) 
C 
      IF(VVLDGAS(IV).NE.1.AND.IV.NE.VTMC.AND.SPD_FLAG.EQ.1) 
     &                       CALL SCFTAB(IV,1) 
C 
C   Set up the Evap weighting factors, without and with I/M, 
C   and without and with OBD, for all ages. Only call this 
C   routine for all gasoline powered vehicles. 
C 
      IF(VVGASMC(IV).EQ.1) CALL EVAPCALC(ICY,IV) 
C 
C ********************** Begin Yearly Loop ************************* 
C 
      DO IDX=1,MAXYRS 
C 
C  Compute JDX (age) 
C 
      JDX=(MAXYRS+1)-IDX 
C 
C  If travel fraction is 0, skip 
C 
      IF(TF(IDX,IV).EQ.0.) CYCLE 
C 
C  Determine model year, modified for no Clean Air Act 
C 
      ACTUAL_MY=ICY-JDX+1 
      MY=ACTUAL_MY 
C 
      IF(CAAFLG.EQ.1 .AND. ICY.GE.CAAMY) MY = MIN(CAAMY,ACTUAL_MY) 
CDFK 
CDFK      WRITE(99,91) IV,MY,TF(IDX,IV) 
CDFK   91 FORMAT(/2I5,F10.5) 
C 
C  Calculate catalyst fractions; only gas vehicles have catalysts. 
C 
      IF(VVGAS(IV) .EQ. 1) THEN 
         CALL P3WPOX(MY,IV,P3W,POX) 
         CAT(1) = 1.0 - POX - P3W   !NoCat 
         CAT(2) = POX               !Ox 
         CAT(3) = P3W               !3-way 
       ELSE 
         CAT(1) = 1.0  !NoCat/Diesel 
         CAT(2) = 0.0  !Ox 
         CAT(3) = 0.0  !3-way 
      ENDIF 
C 
      TCAT = CAT(2) + CAT(3)   !All Catalysts 
C 
C *************** Compute Exhaust Fractions ************************* 
C 
C  Initalize the Air Toxic Conversion Array 
C 
      ATCONV = 0. 
C 
C  Calculate the MOBILE6 VOC and TOG values by hour and road type. 
C  Also calculate the MOBILE6 normal/high emitter fractions. 
C 
      CALL ATHRLOOP(IV,ICY,MY,ACTUAL_MY,IDX,TCAT,INERR) 
CDFK 
CDFK      ATHEF=1. 
CDFK 
CDFK      WRITE(99,92) ATHEF 
CDFK   92 FORMAT(24F10.5) 
C 
CDFK      WRITE(99,94) ATNHWT 
CDFK   94 FORMAT('*',F9.5,F10.5) 
C 
C  Calculate the Air Toxic Technology Group Weighting Factors 
C  for 3-way catalyst LDGV/LDGT1-4 only. 
C 
      IF(VVLDGAS(IV) .EQ. 1   .AND.                     !LDGV/LDGT1-4 
     *        CAT(3) .GT. 0.) THEN                      !3-way 
       CALL ATGROUP(MY,IV,INERR) 
CDFK 
CDFK      WRITE(99,7776) ICY,IV,MY,GROUP 
CDFK 7776 FORMAT('GROUP',I5,I3,I5,10F10.6) 
C 
      END IF 
C 
C  Loop over the various oxygenates. The Complex Model can only 
C  handle one oxygenate type at a time. 
C 
      DO IOXY=1,5  !1:MTBE, 2:ETBE, 3:ETHANOL, 4:TAME, 5:No Oxygenate 
C 
C  Use the Complex Model to obtain EF fractions for 
C  LDGV/T1-4 with 3-way catalysts only. 
C 
      IF(VVLDGAS(IV) .EQ. 1   .AND.                     !LDGV/LDGT1-4 
     *        CAT(3) .GT. 0.  .AND.                     !3-way 
     * ATOXY(2,IOXY) .GT. 0.)                           !Non-zero market share 
     * CALL ATCOMPLEX(IV,MY,IOXY,CAT) 
CDFK 
CDFK      WRITE(99,7777) ICY,IV,MY,IOXY, 
CDFK     &    BENZRATIO,BUTARATIO,ACETRATIO,FORMRATIO,MTBERATIO 
CDFK 7777 FORMAT('RATIO',I5,I3,I5,I2, 
CDFK     &       10F10.6/20X,10F10.6/20X,10F10.6/20X,10F10.6/20X,10F10.6) 
C 
C  Use the alternate equations for all other 
C  vehicle/catalyst combinations. 
C 
      CALL ATCALX2(CAT,IV,IOXY) 
C 
      END DO   !IOXY 
C 
C  Calculate the Acrolein EF fractions 
C 
      CALL ATACRO(CAT,IV) 
CDFK 
CDFK      WRITE(99,93) ICY,IV,MY,ATCONV 
CDFK   93 FORMAT(I5,I3,I5,24F10.5) 
C 
C  Compute the Exhaust Air Toxic emission factors by hour and road type. 
C  Store the values in the data base arrays in g/mi. 
C  Store the values in the descriptive output arrays in mg/mi. 
C  Set IXV=1 since these are exhaust only emissions. 
C  Note that starts and running are added together for the composite EFAT. 
C  (Programming note: The ATCONV values have already been weighted by 
C   catalyst and oxygenate fractions, so they can be multiplied by the 
C   VOC and TOG efs and added directly.) 
C 
      IXV = 1 
C 
      DO IAT=1,MAXIAT !1=Benz, 2=MTBE, 3=Buta, 4=Form, 5=Acet, 6=Acrolein 
C 
       IF(DBSELTOX(IAT) .EQ. 1) CYCLE    !Skip if toxic is not selected 
C 
        SUMJDX=0. 
C 
        DO IH=1,MAXIH 
          SUMIH=0. 
C 
          DO IROAD=1,5          ! 5 = Start Emissions, 1-4 are Road Types 
            SUMINH=0. 
C 
            DO INH=1,2          ! 1 = Normal, 2 = High 
              DO IVT=1,2        ! 1 = VOC, 2 = TOG 
                SUMINH = SUMINH  
     *                 + ATHEF(IH,IROAD,IVT) 
     *                 * ATCONV(INH,IVT,IAT)  
     *                 * ATNHWT(INH) 
              END DO !IVT 
            END DO !INH 
C 
            DBAT(IDX,IH,IROAD,IXV,IAT,IV) = SUMINH  ! g/mi 
C 
             IF(IROAD.LE.4) THEN  
               SUMIH = SUMIH + SUMINH * FVMT(IROAD,IH,IV) ! Road 
             ELSE 
               SUMIH = SUMIH + SUMINH                     ! Start 
             END IF   
            END DO !IROAD 
C 
            SUMJDX = SUMJDX + SUMIH * HVMT(IH)   
           END DO !IH 
C 
           EFAT(IXV,IAT,IV) = EFAT(IXV,IAT,IV)  
     *                      + SUMJDX * TF(IDX,IV) * 1000.  ! mg/mi   
          END DO !IAT 
C 
C ********************* Compute Evap Fractions ************************* 
C 
C  Only perform Evap Air Toxic calculations for gasoline vehicles 
C 
      IF (VVGASMC(IV).EQ.1) THEN  
C 
C  Compute the MOBILE6 hourly resting, running, diurnal, and  
C  hot soak emission factors. 
C 
      CALL ATEVP_HRLY(ICY,IVGAS,IDX,IV,JDX,MY,ACTUAL_MY) 
CDFK 
CDFK 
CDFK      IF(IV.EQ.VTMC) THEN 
CDFK      WRITE(99,771) IV,MY 
CDFK  771 FORMAT(I2,I5) 
CDFK      WRITE(99,772) ATHEV 
CDFK  772 FORMAT(24F10.6) 
CDFK      END IF 
C 
C  Compute the Evap Air Toxic emission factors by hour and road type. 
C  (Currently, the only evap values are Benzene and MTBE.) 
C  Store the values in the data base arrays in g/mi. 
C  Store the values in the descriptive output arrays in mg/mi. 
C 
      DO IAT=1,2 !1=Benzene, 2=MTBE 
C 
      DO IEVAP=1,5 !1=Hot Soak, 2=Diurnal, 3=Running, 4=Resting, 5=Refueling 
C 
      IXV = IEVAP + 1 ! IXV=1 for Exhaust only, 2-6 for the five Evaps 
C 
      IF(IEVAP.EQ.3) THEN  ! Only Running Evap is road dependent 
        IB = 1 
        IE = 4   ! Do not include IROAD=5 to avoid double counting 
      ELSE 
        IB = 5 
        IE = 5 
      END IF 
         
        SUMJDX=0. 
C 
        DO IH=1,MAXIH 
C 
          SUMIH=0. 
C 
          DO IROAD=IB,IE  ! 5 = Road Independent, 1-4 are Road Types 
C 
             DBAT(IDX,IH,IROAD,IXV,IAT,IV) = ATHEV(IH,IROAD,IEVAP)  
     *                                     * ATEVAP(IEVAP,IAT) ! g/mi 
             IF(IROAD.LE.4) THEN 
               SUMIH = SUMIH  
     *               + DBAT(IDX,IH,IROAD,IXV,IAT,IV) * FVMT(IROAD,IH,IV) 
             ELSE 
               SUMIH = SUMIH + DBAT(IDX,IH,IROAD,IXV,IAT,IV) 
             END IF 
C 
          END DO !IROAD 
C 
          SUMJDX = SUMJDX + SUMIH * HVMT(IH)   
         END DO !IH 
C 
        EFAT(IXV,IAT,IV) = EFAT(IXV,IAT,IV) 
     *                   + SUMJDX * TF(IDX,IV) * 1000.  ! mg/mi   
        END DO !IEVAP 
       END DO !IAT 
C 
      END IF !Gas Vehicles 
C 
      END DO !IDX  
C 
C ********************** End of Yearly Loop ************************* 
C 
C  Compute All Vehicle composites for Exhaust and Evap 
C 
      DO IAT=1,MAXIAT 
        DO IXV=1,6  ! IXV=1 for Exhaust, 2-6 for the five Evaps 
          VAT(IXV,IAT) = VAT(IXV,IAT) + EFAT(IXV,IAT,IV) * VMTMIX(IV)  
        END DO  ! IXV 
      END DO  ! IAT 
C 
      END DO !IV 
C 
C ********************** End of Vehicle Loop ************************* 
C 
C  Output the descriptive results. 
C 
      CALL ATOUTDSC(ICY) 
C 
C  Output the database results. 
C 
      CALL ATOUTDB(ICY) 
C 
C  Output the spreadsheet results. 
C 
      CALL ATOUTSS(ICY) 
C 
      RETURN 
      END 
 
       
