      SUBROUTINE PMLEAD(ICY,INERR) 
C 
C  Calculates Lead Particulate Emission Factors, based on PART5. 
C 
C Changes: 
C 
C 15 Nov 01 AIR Task 27: Added INERR support to PMLEADED. 
C 09 Feb 01 AIR Task 11: Modified and merged into MOBIL6.1 Explicit 
C                        use of PART5 vehicle class IVPART. 
C 17 Nov 00 AIR Task 06: Removed Dead Code. 
C 10 Nov 00 AIR Task 06: Updated code to MOBILE6 standards. 
C 04 Mar 94 @ CSC-bsg Request 340 this was added to overhaul 
C           lead emission factor calculation 
C 17 Jun 94 @ CSC-bsg Request 430 MC use lead instead of gas carbon 
C 24 Jun-94 @ CSC-bsg Request 430 MC use lead and carbon 
C 
C  Called by PMCALX. 
C 
C  Calls LEADED, PMLEDPTR, and PMTAMPER. 
C 
C  Input on call: 
C 
C    parameter list: ICY 
C 
C    common blocks:  See Includes below 
C 
C  Output on return: EFPM, VPM, DBPM 
C 
C  Local variable dictionary: 
C 
C   Name   Type              Description 
C  ------  ----  ------------------------------------------------------------- 
C  CADJ     R    Combined ZEV and NDV adjustment 
C  CATFCT   R    Catalyst fraction of the gasoline vehicle population 
C  CATS     R    Fraction of vehicles that have not had the catalyst removed 
C  EFLDRT   R    Lead emission rate   
C  FRACA    R    Catalyst fractions: 1: Ox No Air, 2: Ox Air, 3: 3w No Air, 4: 3w Air 
C  IDX      I    IDX=(MAXYRS+1)-JDX 
C  IH       I    Hour 
C  IROAD    I    Road Type 
C  IV       I    Vehicle Class 
C  IYPT     I    Discretionary fuel switching index variable 
C  JDX      I    Age 
C  MY       I    Model Year 
C  NGVADJ   R    Adjustment for NGVs 
C  NGVPT    I    NGV pointer 
C  PLNOCT   R    Particulate emission factor for non-catalyst vehicles using leaded fuel 
C  PLYSCT   R    Particulate emission factor for catalyst equipped vehicles using leaded fuel 
C  PUNOCT   R    Particulate emission factor for non-catalyst vehicles using unleaded fuel 
C  PUYSCT   R    Particulate emission factor for catalyst equipped vehicles using unleaded fuel 
C  SUM      R    Summation variable 
C  SUMIH    R    Summation over IH variable 
C  VLNOCT   R    Number of non-catalyst vehicles using leaded fuel 
C  VLYSCT   R    Number of catalyst equipped vehicles using leaded fuel 
C  VUNOCT   R    Number of non-catalyst vehicles using unleaded fuel 
C  VUYSCT   R    Number of catalyst equipped vehicles using unleaded fuel 
C  ZEVADJ   R    Adjustment for ZEVs 
C 
C 
      IMPLICIT NONE 
C 
      INCLUDE 'DFSBLK.I'  ! DFS 
      INCLUDE 'IVTYPE.I'  ! MAXVEH, VVGASMC, VM6M5B, VTMC 
      INCLUDE 'MAXIMA.I'  ! MAXIH, MAXYRS 
      INCLUDE 'MYRCAL.I'  ! TF 
      INCLUDE 'NGV.I'     ! NGVPHIN, NGV_START, NGVFLG1 
      INCLUDE 'PART2.I'   ! TAMFRC, RMIS 
      INCLUDE 'PART3.I'   ! EFPM, DBPM, VPM 
      INCLUDE 'SPEED9.I'  ! HVMT, FVMT 
      INCLUDE 'VMXCOM.I'  ! VMTMIX 
C 
C  Declare parameter list. 
C 
      INTEGER ICY 
      INTEGER INERR 
C 
C  Declare external functions. 
C 
      REAL, EXTERNAL :: PCTLEV 
C 
C  Declare local variables/arrays. 
C 
      INTEGER IDX 
      INTEGER IH 
      INTEGER IROAD 
      INTEGER IV 
      INTEGER IVPART 
      INTEGER IYPT 
      INTEGER JDX 
      INTEGER MY 
      INTEGER NGVPT 
      REAL    CADJ 
      REAL    CATFCT 
      REAL    CATS 
      REAL    EFLDRT 
      REAL    NGVADJ 
      REAL    PLNOCT 
      REAL    PLYSCT 
      REAL    PUNOCT 
      REAL    PUYSCT 
      REAL    SUM 
      REAL    SUMIH 
      REAL    VLNOCT 
      REAL    VLYSCT 
      REAL    VUNOCT 
      REAL    VUYSCT 
      REAL    ZEVADJ 
C 
C 
C  Loop for each MOBILE6 vehicle class 
C 
      DO IV=1,MAXVEH 
C 
C  Skip Diesels. 
C 
      IF(VVGASMC(IV).EQ.0) CYCLE 
C 
C  Skip if VMTMIX is zero 
C 
      IF(VMTMIX(IV).EQ.0.) CYCLE 
C 
C  Map MOBILE6 vehicle class into PART5 space 
C 
      IVPART=VM6M5B(IV) 
C 
C  Loop for all ages 
C 
      DO JDX=1,MAXYRS 
C 
C  Set IDX 
C 
      IDX=(MAXYRS+1)-JDX 
C 
C  If the travel fraction for this age/vehicle combination is zero, skip. 
C 
      IF(TF(IDX,IV) .EQ. 0.) CYCLE 
C 
C  Set model year 
C 
      MY=ICY-JDX+1 
C 
C  Determine the ZEV adjustment for LDGV, LDGT1-LDGT4. 
C 
      ZEVADJ = 1.0 - PCTLEV(MY,IV,11) 
C 
C  Compute the NGV adjustment 
C 
      IF(NGVFLG1.EQ.1) THEN 
         IF(MY.LT.NGV_START) THEN 
            NGVADJ = 1.0 
         ELSE 
            NGVPT = MY - NGV_START + 1 
            NGVADJ = 1.0 - NGVPHIN(NGVPT,IV) 
         END IF 
      ELSE 
         NGVADJ = 1.0 
      END IF 
C 
C  Combine the NGV and ZEV adjustments (neither ever emit lead). 
C  This form assumes that ZEVs are removed first, and the NGVs 
C  are from the remainder. 
C 
      CADJ = ZEVADJ * NGVADJ 
C 
C  Determine the lead particulate emission rates 
C 
      CALL PMLEADED(MY,ICY,IV,PLNOCT,PUNOCT,PLYSCT,PUYSCT,INERR) 
C 
C  MC use lead and carbon 
C 
      IF(IV.EQ.VTMC) THEN 
        EFLDRT = PLNOCT + PUNOCT + PLYSCT + PUYSCT 
C 
C  Compute lead values for remaining classes 
C 
      ELSE 
C 
C  Get the catalyst fraction of the gasoline vehicle population. 
C  (Use MOBILE6 computed values instead of PART5's FRACA.) 
C 
      CALL PMFRACA(MY,IV,INERR) 
       
      CATFCT = FRACA(1) + FRACA(2) + FRACA(3) + FRACA(4) 
C 
C  Determine discretionary fuel switching index variable. 
C  Use PART5 vehicle class. 
C 
      CALL PMLEDPTR(IYPT,IV,MY) 
C 
C  determine the tampering and misfueling rates. 
C  (Must pass M6 and PART5 vehicle classes) 
C 
      CALL PMTAMPER(ICY,MY,IV) 
C 
C  determine the fraction of all vehicles which 
C  have NOT had the catalyst removed. 
C 
      CATS = CATFCT * (1.-TAMFRC) 
C 
C  determine the fraction of all vehicles which are 
C  non-catalyst vehicles fueled with leaded gasoline. 
C 
      VLNOCT = (1.-CATS) * (1.-DFS(IYPT,IVPART)) 
C 
C  determine the fraction of all vehicles which are 
C  non-catalyst vehicles fueled with unleaded gasoline. 
C 
      VUNOCT = (1.-CATS) * DFS(IYPT,IVPART) 
C 
C  determine the fraction of all vehicles which are 
C  catalyst equipped vehicles fueled with leaded gasoline. 
C 
      VLYSCT = CATFCT * (RMIS - TAMFRC) 
C 
C  determine the fraction of all vehicles which are 
C  catalyst equipped vehicles fueled with unleaded gasoline. 
C 
      VUYSCT = (1.-VLNOCT-VUNOCT-VLYSCT) 
C 
C  determine the combined lead particulate emission factor 
C 
      EFLDRT = 
     &        PLNOCT * VLNOCT 
     &      + PUNOCT * VUNOCT 
     &      + PLYSCT * VLYSCT 
     &      + PUYSCT * VUYSCT 
C 
      ENDIF 
C 
C  Adjust for ZEVs/NGVs 
C 
      EFLDRT = EFLDRT * CADJ 
C 
C  Loop for each Hour and Road Type. Compute composites. 
C 
          SUMIH=0. 
C 
          DO IH=1,MAXIH 
            SUM=0. 
C 
            DO IROAD=1,4 
              DBPM(5,IROAD,IH,JDX,IV)=EFLDRT 
              SUM = SUM + DBPM(5,IROAD,IH,JDX,IV) * FVMT(IROAD,IH,IV) 
            END DO !IROAD 
C 
            SUMIH = SUMIH + SUM * HVMT(IH) 
C 
          END DO   !IH 
C 
          EFPM(5,IV)=EFPM(5,IV)+SUMIH*TF(IDX,IV) 
C 
        END DO     !JDX 
C 
        VPM(5)=VPM(5)+EFPM(5,IV)*VMTMIX(IV) 
C 
      END DO       !IV 
C 
      RETURN 
      END 
