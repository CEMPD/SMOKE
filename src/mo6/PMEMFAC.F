      SUBROUTINE PMEMFAC(ICY,INERR) 
C 
C  Computes the Total PM (PM = RCP + SOF + SO4) based on the EPA methodology. 
C  The output values are at particle size PM30, with 50 ppm sulfur for gasoline, 
C  and 500 ppm suflur for diesel. 
C 
C  Changes : 
C 
C  07 Feb 01 AIR Task 15-23c: Fixed bug where BERSULF wasn't being initialized 
C                          to the HDR 2007 rule properly. 
C  28 Aug 01 AIR Task 23c: Implement DR1, DR2 and PMDRAGE 
C  24 Aug 01 AIR Task 23b: Implement EPA-based ZMLs/DRs and methodology. 
C  06 Jun 01 AIR Task 17: Added CA diesel aromatic control program code. 
C  23 May 01 AIR Task 17: Made No-Speed-Corrections as default. Removed 4-Mode 
C            ZMLs/DRs and made MY start in 1950- in stead of 1975-. 
C  04 May 01 AIR Task 12: Tier2/2007 Rule 
C  17 Apr 01 AIR Task 11: Changed code to handle by model year 
C            sulfur correction for gasoline fueled vehicles and 
C            expanded FVMT dimensions. 
C  01 Feb 01 AIR Task 11: New Module 
C 
C     Called by:  PMCALX 
C 
C     Calls: PMSCF 
C 
C  Input on call: 
C 
C    Parameter list: ICY 
C 
C    Common blocks:  See Includes below 
C 
C  Output on return:  
C    /PART3.I/ EFPM,VPM,DBPM 
C 
C  Local arrays: 
C 
C  FCAT(2)       - FCAT (ICAT)     
C  PMBER(2,2)    - PMBER(ICAT) 
C 
C  Local variable / array dictionary: 
C 
C   Name        Type              Description 
C  ------       ----  ----------------------------------------------------- 
C  BERMAX        R    Maximum BER allowed by 2007 HD Rule 
C  BODOM         R    Base odometer mileage 
C  FCAT          R    Fraction without (1) and with (2) catalyst 
C  IBIN          I    BIN number 
C  ICAT          I    Loop for catalyst-equipped, 1=without, 2=with 
C  ICY           I    Calendar Year 
C  IDX           I    IDX=(MAXYRS+1)-JDX 
C  IG5           I    Technology pointer 
C  IH            I    Hour index 
C  IMY           I    Temporary model year 
C  IROAD         I    Road index 
C  IV            I    Vehicle Class 
C  IVTEMP        I    temporary IV 
C  JDX           I    Age 
C  KODOM         R    Mileage after kink 
C  MY            I    Model year 
C  MYPM          I    model year pointer into PM_BER arrays, 1950-=1, 2020+=71 
C  P3W           R    Fraction of 3-way catalysts 
C  PMBER         R    Basic emmission rates without (1) and with (2) catalyst 
C  POX           R    Fraction of oxidation catalysts 
C  SCOR          R    Sulfur correction value for diesel fuel at 500 ppm, MY dependent 
C  SUM           R    Temperory summation 
C  SUMIH         R    Temporary summation over IH 
C  XPM           R    Temporary PM 
C 
      IMPLICIT NONE 
C 
      INCLUDE 'CUMCOM.I'   ! CUMMIL 
      INCLUDE 'FLAGS5.I'   ! SOAK_FLG,TIER2_FLAG,HDR2007 
      INCLUDE 'IVPCOM.I'   ! IVPTRT 
      INCLUDE 'IVTYPE.I'   ! MAXVEH,VTHDG2B,VTHDG8B,VTHDD2B,VTHDD8B,VTGASBUS,VTSCHOOL,VVLDDSL 
      INCLUDE 'MAXIMA.I'   ! MAXYRS, MAXIH 
      INCLUDE 'MYRCAL.I'   ! TF 
      INCLUDE 'PART1.I'    ! PRTCHK, D_PPM, BERSULF 
      INCLUDE 'PART3.I'    ! EFPM, DBPM, VPM 
      INCLUDE 'PART4.I'    ! PMZML, PMDR1, PMDR2, BSULFLVL 
      INCLUDE 'PART5.I'    ! PMSCFC,PMSPDFLG 
      INCLUDE 'SPEED9.I'   ! FVMT, HVMT 
      INCLUDE 'SULFUR.I'   ! S_PPM 
      INCLUDE 'TIER2PM.I'  ! T2PM_50K, T2PM_FRAC 
      INCLUDE 'VMXCOM.I'   ! VMTMIX 
C 
C  Declare external functions. 
C 
      INTEGER, EXTERNAL :: IEVPTR 
      REAL   , EXTERNAL :: HDCPIC 
      REAL   , EXTERNAL :: PMSCF 
C 
C  Declare parameter list 
C 
      INTEGER ICY 
      INTEGER INERR 
C 
C  Declare local variables/arrays 
C 
      INTEGER IAGE 
      INTEGER IBIN 
      INTEGER ICAT 
      INTEGER IDX 
      INTEGER IH 
      INTEGER IMY 
      INTEGER IROAD 
      INTEGER IV 
      INTEGER IVTEMP 
      INTEGER JDX 
      INTEGER MY 
      INTEGER MYPM 
      REAL    BERMAX 
      REAL    BODOM 
      REAL    FCAT(2) 
      REAL    KODOM 
      REAL    P3W 
      REAL    PMBER(2) 
      REAL    POX 
      REAL    SUM 
      REAL    SUMIH 
      REAL    XPM 
C 
C 
C  Initialize the BER associated suflur levels. 
C 
      BERSULF = 0. 
C 
C  If HDD2007 rule is in effect and the calendar year diesel 
C  sulfur is above 15 ppm, warn user. Do Not reset the sulfur level. 
C  This warning is issued because the diesel ECARBON and OCARBON 
C  values may become zero. 
C 
      IF(HDR2007.EQ.1.AND.ICY.GE.2007.AND.D_PPM.GT.15.) 
     *   CALL QUITER(D_PPM,0,111,INERR) 
C 
C  Loop for each vehicle class 
C 
      DO IV=1,MAXVEH 
C 
C  Skip if VMTMIX is 0. 
C 
      IF(VMTMIX(IV).EQ.0.) CYCLE 
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
C  Determine the age and apply diesel DR age limits as follows: 
C   12 years for HDDT2B-3 
C    5 years for HDDT4-8B, HDDBT, HDDBS 
C  The diesel limits take into account engine rebuilds. 
C 
      IAGE = JDX       
C 
      IF((IV.EQ.VTHDDV2B .OR. IV.EQ.VTHDDV3) .AND. 
     &    JDX.GT.12) IAGE = 12 
C 
      IF((IV.GE.VTHDDV4 .AND. IV.LE.VTHDDV8B) .OR. 
     &   (IV.EQ.VTURBAN . OR. IV.EQ.VTSCHOOL) .AND. 
     &    JDX.GT.5) IAGE = 5 
C 
C  Determine the base and kink mileage, in 10K. 
C 
      IF(IAGE.LE.PMDRAGE(IV)) THEN 
        BODOM=CUMMIL(IAGE,IV)/10000. 
        KODOM=0. 
      ELSE 
        IF(PMDRAGE(IV).EQ.0) THEN 
          BODOM=0. 
        ELSE 
          BODOM=CUMMIL(PMDRAGE(IV),IV)/10000. 
        END IF 
        KODOM=CUMMIL(IAGE,IV)/10000. - BODOM 
      END IF       
C 
C  Determine the model year pointer into the PM_BER arrays. 
C 
      MYPM=MY-1950+1 
      IF(MYPM.LT. 1) MYPM= 1 
      IF(MYPM.GT.71) MYPM=71 
C 
C  COMPUTE THE BERS 
C 
C  Loop for catalyst 
C 
      DO ICAT=1,2 
C 
C  Compute the basic emission rates (BER) in g/mi. 
C 
        PMBER(ICAT) = PMZML(MYPM,ICAT,IV) 
     *              + PMDR1(MYPM,ICAT,IV) * BODOM 
     *              + PMDR2(MYPM,ICAT,IV) * KODOM 
C 
      END DO ! ICAT 
C 
C  Determine the fraction of vehicles without/with catalysts 
C 
      CALL P3WPOX(MY,IV,P3W,POX) 
C 
      FCAT(2)=P3W+POX 
      FCAT(1)=1.0-FCAT(2) 
C 
C  Determine the associated fuel sulfur content of the BER. (For gas, the  
C  BERs are based on the input sulfur for CY1999- and CY2000+. For diesel, 
C  the BERs are based on a single input sulfur level.) Store the sulfur values 
C  for use in determining the amount of associated SO4 to remove later. 
C  BSULFLVL now contains the base sulfur levels input by the user. Its array 
C  order is 1=Gas CY1999-, 2=Gas CY2000+, 3=Diesel 
C 
      IF(VVDSL(IV).EQ.1) THEN 
        BERSULF(JDX,IV)=BSULFLVL(3) 
      ELSE 
        IF(ICY.LE.1999) THEN 
          BERSULF(JDX,IV)=BSULFLVL(1) 
        ELSE 
          BERSULF(JDX,IV)=BSULFLVL(2) 
        END IF 
      ENDIF 
C 
C  Gas BINs are applicable to Tier2 LDGV/T MY2004+ when the TIER2 flag is on. 
C  Diesel BINs are applicable to LDDV/T MY2007+ when the HDV 2007 rule flag 
C  is on and the diesel fuel sulfur level is 15 ppm or less. The BIN weighted 
C  PM emission factor limit is applied to the BERs. 
C  Note: The associated fuel sulfur content for TIER2 BERs is 30 ppm for 
C  gas and 8 ppm diesel. 
C 
      IF((TIER2_FLAG .EQ. 1         .AND.        ! LDGT/LDGT applicability 
     *            MY .GE. 2004      .AND. 
     *            IV .LE. VTLDGT4) 
     *                              .OR. 
     *      (HDR2007 .EQ. 1         .AND.        ! LDDV/LDDT applicability 
     *            MY .GE. 2007      .AND. 
     *         D_PPM .LE. 15.       .AND. 
     *   VVLDDSL(IV) .EQ. 1)) THEN 
C 
        IVTEMP=IV                               ! Map IV into BINs 
        IF(IV.EQ.VTLDDV .OR. 
     *     IV.EQ.VTLDDT12) IVTEMP=IV-8 
        IF(IV.EQ.VTLDDT34) IVTEMP=8 
C 
        IMY=MY-2004+1                           ! Map MY into BINs 
        IF(IMY.GT.12) IMY=12 
C 
        BERMAX=0. 
C 
        DO IBIN = 1, 12                         ! Compute BIN Limit 
          BERMAX = BERMAX  
     *           + T2PM_50K(IVTEMP,IBIN) 
     *           * T2PM_FRAC(IMY,IBIN,IVTEMP) 
        END DO 
C 
        BERMAX = BERMAX * HDCPIC(MY,IV)         ! Convert BIN to g/mi 
C 
        SUM=0. 
C 
          DO ICAT=1,2                           ! No Cat + Cat 
            SUM=SUM+PMBER(ICAT)*FCAT(ICAT) 
          END DO 
C 
        IF(SUM.GT.BERMAX) THEN 
            DO ICAT=1,2                         ! Apply Limit 
              PMBER(ICAT) = PMBER(ICAT) * BERMAX/SUM 
            END DO 
        END IF 
C 
          IF(VVDSL(IV).EQ.1) THEN               ! Store associated sulfur level 
            BERSULF(JDX,IV)=8. 
          ELSE 
            BERSULF(JDX,IV)=30. 
          ENDIF 
C 
      END IF 
C 
C  Determine if 2007 HDV Rule is in effect. If so, and diesel fuel 
C  level is 15 ppm or less, apply limit of 0.0092 g/bhp-hr to PMBER 
C  (start + running).  Apply same limit to HDGV, regardless of  
C  gasoline sulfur level.  
C  Note: The associated fuel sulfur content for 2007 HDV BERs is 30 ppm for 
C  gas and 8 ppm diesel. 
C  (Programming Note: currently the deterioration rate for the 0.0092 
C   g/bhp-hr standard is set to 0.0.) 
C 
      IF(HDR2007.EQ.1.AND.MY.GE.2007) THEN             ! HDV 2007 applicability 
C 
      IF(VVHDDV(IV).EQ.1.AND.D_PPM.LE.15. .OR.         ! Fuel sulfur limit 
     *   VVHDGV(IV).EQ.1) THEN 
C 
        BERMAX = 0.0092 * HDCPIC(MY,IV)                ! Compute 2007 HDV BER limit in g/mi. 
C 
        SUM=0. 
C 
          DO ICAT=1,2                                  ! No Cat + Cat  
            SUM=SUM+PMBER(ICAT)*FCAT(ICAT) 
          END DO 
C 
        IF(SUM.GT.BERMAX) THEN 
            DO ICAT=1,2                                ! Apply limit 
              PMBER(ICAT) = PMBER(ICAT) * BERMAX/SUM 
            END DO 
        END IF 
C 
          IF(VVDSL(IV).EQ.1) THEN                      ! Store associated sulfur level 
            BERSULF(JDX,IV)=8. 
          ELSE 
            BERSULF(JDX,IV)=30. 
          ENDIF 
C 
      END IF 
      END IF 
C 
C  BUILD BY-MODEL-YEAR AND COMPOSITE EMISSION FACTORS 
C 
      SUMIH = 0. 
C 
C  Loop for each hour 
C 
      DO IH = 1,MAXIH 
C 
       SUM=0. 
C 
C  Loop for each road type 
C 
       DO IROAD = 1,4 
C 
C  Loop for catalyst 
C 
         DO ICAT=1,2 
C 
C  Apply the MOBILE6-based Cat/No Cat fractions. 
C 
            XPM = PMBER(ICAT) * FCAT(ICAT) 
C 
C  Build the Total PM here. 
C 
            DBPM(10,IROAD,IH,JDX,IV) = DBPM(10,IROAD,IH,JDX,IV) + XPM 
C 
           END DO  !ICAT 
C 
          SUM = SUM + DBPM(10,IROAD,IH,JDX,IV) * FVMT(IROAD,IH,IV) 
C 
         END DO    !IROAD 
C 
         SUMIH = SUMIH + SUM * HVMT(IH) 
C 
        END DO     !IH 
C 
        EFPM(10,IV)= EFPM(10,IV) + SUMIH * TF(IDX,IV) 
C 
       END DO ! JDX 
C 
       VPM(10)= VPM(10) + EFPM(10,IV) * VMTMIX(IV) 
C 
      END DO ! IV 
C 
      RETURN 
      END 
