      SUBROUTINE DIRNLDAY(MY,ACTUAL_MY,ICY,IVGAS) 
C 
C  This subroutine will calculate the daily diurnal emission 
C 
C  Called by EVP_HRLY 
C  Calls EVWTPRE, EVWTENH 
C 
C  Changes:  (Last change first) 
C 
C  27 Mar 01 @EPA-BAg bug321 added LESS_THN10_FCTR to motorcycle diurnals 
C  27 Feb 01 @EPA-BAg bug314 revised tier2 effects on heavy duty vehicles 
C  31 Jan 01 @EPA-bag added high altitude adjustment; no longer applying the factor 
C            to liquid leaks 
C  18 Jan 01 @EPA-bag 2007 rule for heavy duty vehicles 
C  13 Dec 00 @EPA-bag t-x6 MY and ACTUAL_MY added to parameter list for CAA switch 
C  31 Oct 00 @EPA-djb Bug290, Capped pass/pass diurnal emissions at failing levels, 
C            such that passing <= fail purge <= fail pressure <= gross liquid leaks. 
C  01 Nov 00 @EPA-bag BUG296 precontrol vehicles; new logic for pre-control 
C            vehicles 
C  07 Aug 00 @EPA bug252&253 1996-98 evaporative weighting for liquid leaks; vehicle 
C            type diurnal categorizations 
C  22 Aug 00 @EPA-bag  bug25 1972-1978 ldgt3&4 should use hdgv2b&3 parameters 
C  21 Jul 00 @EPA-elg  Added code to account for Tier2 Evap on IV=6 vehicles (2b's). 
C  13 Aug 00 @EPA-bag for clarity substituted DIU1_PRE72 for DIURNL24A in first 
C            portion of pre-1972 diurnal calculation 
C  07 Aug 00 @EPA-bag bug215 Added DIURNL24B array; they are weighted 1999 portion 
C            of the full day diurnals for the model years 1996-1998; DIURNL24A are 
C            are used for all model years and are also the 1995 portion of the full 
C            day diurnals for the model years 1996-1998 
C  07 Jul 00 @EPA-ddj    e-15  DIU1 and DIU2 arrays are no longer weighted 
C  07 Jul 00 @EPA-ddj    e-15  moved DIU1 and DIU2 arrays into DIRNLDAT module 
C                        for multiday diurnals 
C  07 Jan 00 @Dyntel-bg 1-005 updated diurnal code(new) 
C 
C  Input on call: 
C 
C    Constants:  
C 
C    argument list: MY,ICY,IVGAS 
C    functions :    EVWTENH, EVWTPRE 
C    modules   : 
C    DITEMDAT - LESS_THN10_FCTR,VPPR 
C    DIRNLDAT - COMC,COLDV,DIRNL_DAY_LQD,WA,WB,COLDVA,COLDVB, 
C               COLDVC,COLDVD,DIURNL24A,DIURNL24B,LDGV_FCTR, 
C               HDGV_FCTR1,HDGV_FCTR2,DIU1,DIU2,DV5 
C 
C    common blocks: 
C    /EVAL/      MEVAL 
C    /EVAPHD/    HDWGT 
C    /EVAPAR/    FINJ 
C    /EVAHIA/    HIADJ 
C    /FLAGS5/    TIER2_FLAG 
C    /IVTYPE/    LDGV,LDGT1,LDGT2,LDGT3,LDGT4 
C    /REGION/    IREJN 
C    /TIER2/     T2_EVAP,T2_RED 
C 
C 
C  Output on return: 
C 
C    argument list: 
C    modules      : 
C    DIRNLDAT - DIURNL24 
C    common blocks: 
C 
C  Local variable / array dictionary: 
C 
C   Name          Type              Description 
C  ------         ----  ------------------------------------------------------- 
C  HOT_DIURN_RUN   I    diurnal parameter for EVWT, HOT_DIURN_RUN=2 
C  PRE1972         I    model year range parameter, PRE1972=0 
C  NO_SPECL_CAS    I    model year range parameter, NO_SPECL_CAS=1 
C  BTWEEN_96_98    I    model year range parameter, BTWEEN_96_98=2 
C  DIU1_PRE72      R    temporary storage of unweighted daily diurnal 
C                       which is used for comparison with 1972 carb. fail 
C                       pressure result 
C  IYES            I    used to compare flag values; it has the value of 1 
C  INO             I    used to compare flag values; it has the value of 0 
C  IVGAS_MC        I    used to distinguish for IVGAS = 14, i.e., motorcycles 
C  AGE             I    vehicle age, AGE=ICY-MY 
C  IFUEL           I    DO loop index for carbureted and fuel injection 
C                       IFUEL mappings : carbureted vehicles,    IFUEL=1 
C                                        fuel injected vehicles, IFUEL=2 
C  IPPS            I    DO loop index for the three pass/fail pressure and purge tests 
C                       categories 
C                       IPPS mappings : fail press(IPPS=1) 
C                                       fail purge(IPPS=2) 
C                                       pass both (IPPS=3) 
C  IPPS_NUM        I    upper limit of IPPS DO loop for the three pass/fail pressure and purge tests 
C                       categories 
C  IFUEL_SYS       I    maps IFUEL and the pass/fail pressure and purge test 
C                       categories into a single index, 
C                       IFUEL_SYS = IPPS + I_0OR3(IFUEL)    [I_0OR3(1)=0 & I_0OR3(2)=3] 
C                                   1: carbureted     fail pressure 
C                                   2:   "            fail pruge 
C                                   3:   "            pass both 
C                                   4: fuel injection fail pressure 
C                                   5:   "     "      fail purge 
C                                   6:   "     "      pass both 
C  I_0OR3          I    used to map IFUEL_SYS index for daily diurnals to the loop index, 
C                       IPPS, i.e., IFUEL_SYS = IPPS + I_0OR3(IFUEL) 
C  IFUEL_NUM1      I    sets the lower limit to the IFUEL DO loop index, allows 
C                       for loop to be set for carbureted and fuel injection, 
C                       carbureted only, or fuel injection only 
C  IFUEL_NUM2      I    sets the upper limit to the IFUEL DO loop index, allows 
C                       for loop to be set for carbureted and fuel injection, 
C                       carbureted only, or fuel injection only 
C  IMY             I    used with the model year ranges to set the model year index 
C                       on the parameters for the daily diurnals 
C  JKAGE           I    vehicle age passed to evaporative weighting function, JKAGE=ICY-MY+1 
C  MOTRCYCL        I    has value of IYES=1 or INO=0 depending on whether or not the vehicle 
C                       type is motorcycle; used to skip over calculations not used for motorcycles 
C  MYRNG           I    model year range parameter which l 
C  HEAVY_FCTR      R    multiplicative factor for heavy duty vehicles 
C  ALPHA           R    temporary storage for results of daily diurnal parameterization 
C                       used to check for results less than zero 
C  PP_WT           R    pass/fail weighting 
C  FDS_WT          R    fuel delivery system weighting 
C                       needed for mixed years, 1996-1998 
C  ZERO_LIMIT      R      parameter to limit differences which are close to zero 
C 
C  Notes:  diurnal emissions from liquid leaks is determined in DIRNLINI 
C 
      USE DITEMDAT, ONLY : LESS_THN10_FCTR,VPPR 
      USE DIRNLDAT, ONLY : COMC,COLDV,W95,W99,COLDVA,COLDVB,COLDVC, 
     *                     COLDVD,DIURNL24A,DIURNL24B,LDGV_FCTR, 
     *                     DIU1,DIU2,DV5 
      IMPLICIT NONE 
C 
      INCLUDE 'EVAL.I' 
      INCLUDE 'EVAPHD.I' 
      INCLUDE 'EVAHIA.I' 
      INCLUDE 'EVAPAR.I' 
      INCLUDE 'FLAGS5.I' 
      INCLUDE 'IVTYPE.I' 
      INCLUDE 'INJECT.I' 
      INCLUDE 'REGION.I' 
C 
C     Declare parameters. 
C 
      INTEGER, INTENT(IN) :: ICY 
      INTEGER, INTENT(IN) :: MY        ! Model year (corrected if NO CLEAN AIR ACT (NO CAA) is enabled) 
      INTEGER, INTENT(IN) :: ACTUAL_MY ! Model year (uncorrected; ACTUAL_MY=ICY+IDX-MAXYRS) 
      INTEGER, INTENT(IN) :: IVGAS 
C 
C     Declare called functions 
C 
      REAL :: EVWTENH 
      REAL :: EVWTPRE 
      REAL :: TR2FCTR 
C 
C     Declare local variables and constants. 
C 
      INTEGER, SAVE :: IDIURNAL=2 
      INTEGER, SAVE :: PRE1972=0,NO_SPECL_CAS=1,BTWEEN_96_98=2 
      INTEGER, SAVE :: IYES=1,INO=0 
      REAL,    SAVE :: ZERO_LIMIT=0.001 
      REAL    :: DIU1_PRE72 
      REAL    :: AGE 
      INTEGER :: IFUEL 
      INTEGER :: IFUEL_SYS 
      INTEGER :: IFUEL_NUM1 
      INTEGER :: IFUEL_NUM2 
      INTEGER :: IMY 
      INTEGER :: IPPS 
      INTEGER :: IPPS_NUM 
      INTEGER :: JKAGE 
      INTEGER :: MOTRCYCL 
      INTEGER :: MYRNG 
      REAL    :: HEAVY_FCTR 
      REAL    :: ALPHA 
      REAL, DIMENSION(4) :: EVP_WT 
      REAL, DIMENSION(4) :: EVP96_98 
      REAL, DIMENSION(2) :: FDS_WT 
      INTEGER, DIMENSION(2), SAVE :: I_0OR3 = (/0,3/) 
C 
C     initializations 
C 
      ALPHA  = 0.0 
      EVP_WT = 0.0 
      EVP96_98  = 0.0 
      FDS_WT(1) = 1.0  !initialize all vehicles to carbureted vehicles 
      FDS_WT(2) = 0.0  !initialize all vehicles to carbureted vehicles 
      DIU1(1:6) = 0.0  !liquid leaks, DIU1(7), are initialized in DIRNLINI.FOR 
      DIU2(1:6) = 0.0  !liquid leaks, DIU2(7), are initialized in DIRNLINI.FOR 
      DIU1_PRE72 = 0.0 
      DIURNL24A = 0.0 
      DIURNL24B = 0.0 
      HEAVY_FCTR=LDGV_FCTR 
      MOTRCYCL=INO 
      IMY=-1 
C 
      IF( LESS_THN10_FCTR < ZERO_LIMIT ) GO TO 99 ! no diurnals for zero 
                                                  ! max. - min temperature 
                                                  ! differences 
C 
C 
C     For runs based on January first 
      IF(MEVAL==1) THEN 
        AGE=ICY-ACTUAL_MY 
C 
C     For runs based on July first 
      ELSE 
        AGE=ICY-ACTUAL_MY+0.5 
C 
      END IF 
C 
C     Calculate the weighting factor 
C     PURGE/Pressure   Fail Press(1)  Fail Purge(2)  Pass Both(3) 
C                      gross liquid leakers(4) 
      JKAGE=ICY-ACTUAL_MY+1 
      IF(MY < 1996) THEN 
        EVP_WT(1)=EVWTPRE(MY,ACTUAL_MY,JKAGE,3,IDIURNAL) 
     *        +EVWTPRE(MY,ACTUAL_MY,JKAGE,4,IDIURNAL) 
        EVP_WT(2)=EVWTPRE(MY,ACTUAL_MY,JKAGE,2,IDIURNAL) 
        EVP_WT(3)=EVWTPRE(MY,ACTUAL_MY,JKAGE,1,IDIURNAL) 
        EVP_WT(4)=EVWTPRE(MY,ACTUAL_MY,JKAGE,5,IDIURNAL) 
        EVP96_98(1)=0.0 
        EVP96_98(2)=0.0 
        EVP96_98(3)=0.0 
        EVP96_98(4)=0.0 
      ELSEIF(MY>=1996 .AND. MY<=1998) THEN 
        EVP_WT(1)=W95(MY-1995)*(EVWTPRE(MY,ACTUAL_MY,JKAGE,3,IDIURNAL) 
     *        +EVWTPRE(MY,ACTUAL_MY,JKAGE,4,IDIURNAL)) 
        EVP_WT(2)=W95(MY-1995)*EVWTPRE(MY,ACTUAL_MY,JKAGE,2,IDIURNAL) 
        EVP_WT(3)=W95(MY-1995)*EVWTPRE(MY,ACTUAL_MY,JKAGE,1,IDIURNAL) 
        EVP_WT(4)=W95(MY-1995)*EVWTPRE(MY,ACTUAL_MY,JKAGE,5,IDIURNAL) 
        EVP96_98(1)=W99(MY-1995)*(EVWTENH(MY,ACTUAL_MY,JKAGE,3,IDIURNAL) 
     *        +EVWTPRE(MY,ACTUAL_MY,JKAGE,4,IDIURNAL)) 
        EVP96_98(2)=W99(MY-1995)*EVWTENH(MY,ACTUAL_MY,JKAGE,2,IDIURNAL) 
        EVP96_98(3)=W99(MY-1995)*EVWTENH(MY,ACTUAL_MY,JKAGE,1,IDIURNAL) 
        EVP96_98(4)=W99(MY-1995)*EVWTENH(MY,ACTUAL_MY,JKAGE,5,IDIURNAL) 
      ELSEIF(MY > 1998) THEN 
        EVP_WT(1)=EVWTENH(MY,ACTUAL_MY,JKAGE,3,IDIURNAL) 
     *        +EVWTENH(MY,ACTUAL_MY,JKAGE,4,IDIURNAL) 
        EVP_WT(2)=EVWTENH(MY,ACTUAL_MY,JKAGE,2,IDIURNAL) 
        EVP_WT(3)=EVWTENH(MY,ACTUAL_MY,JKAGE,1,IDIURNAL) 
        EVP_WT(4)=EVWTENH(MY,ACTUAL_MY,JKAGE,5,IDIURNAL) 
        EVP96_98(1)=0.0 
        EVP96_98(2)=0.0 
        EVP96_98(3)=0.0 
        EVP96_98(4)=0.0 
      ENDIF 
C 
C     Calculate the weighting factors 
C     for carbureted (1) and fuel injection (2) 
C 
      FDS_WT(1)=1.0-FINJ(1)-FINJ(2) 
      IF(FDS_WT(1)<=0.00001) FDS_WT(1) = 0.0 
      FDS_WT(2)=FINJ(1)+FINJ(2) 
      IF(FDS_WT(2)<=0.00001) FDS_WT(2) = 0.0 
C 
C 
      IFUEL_NUM1=1 
      SELECT CASE (IVGAS) 
          CASE (:5)        ! light duty gas 
            HEAVY_FCTR=LDGV_FCTR 
            IPPS_NUM=3 
            IFUEL_NUM2=2 
          CASE (6:7)       ! heavy duty gas 2b,3 
            HEAVY_FCTR=HDWGT(1) 
            IPPS_NUM=1     !use fail press(IPPS=1) parameters 
            IFUEL_NUM2=1   !use carbureted vehicles only 
          CASE (8:13,15)   ! heavy duty gas 4,5,6,7,8a,8b,gas bus 
            HEAVY_FCTR=HDWGT(2) 
            IPPS_NUM=1     !use fail press(IPPS=1) parameters 
            IFUEL_NUM2=1   !use carbureted vehicles only 
          CASE (14)        ! motorcycle 
            MOTRCYCL=IYES 
            DIU1(1)=(COMC(1)+COMC(2)*VPPR+COMC(3)*AGE)*LESS_THN10_FCTR 
            IF(DIU1(1) < 0.00001) DIU1(1)=0.0 
            DIU1(2:6) = 0.0 
            DIURNL24A(1) = EVP_WT(1)*DIU1(1) 
            DIURNL24A(2) = EVP_WT(2)*DIU1(1) 
            DIURNL24A(3) = EVP_WT(3)*DIU1(1) 
            DIURNL24A(4) = 0.0 
            DIURNL24A(5) = 0.0 
            DIURNL24A(6) = 0.0 
            DIURNL24A(7) = EVP_WT(4)*DIU1(7)*LESS_THN10_FCTR ! DIU1(7) = COGLL(1)+COGLL(2)*DIFF_TEMP 
C 
            DIURNL24B(1) = EVP96_98(1)*DIU1(1) 
            DIURNL24B(2) = EVP96_98(2)*DIU1(1) 
            DIURNL24B(3) = EVP96_98(3)*DIU1(1) 
            DIURNL24B(4) = 0.0 
            DIURNL24B(5) = 0.0 
            DIURNL24B(6) = 0.0 
            DIURNL24B(7) = EVP96_98(4)*DIU1(7)*LESS_THN10_FCTR ! DIU1(7) = COGLL(1)+COGLL(2)*DIFF_TEMP 
      END SELECT 
C 
C 
      IF(MOTRCYCL==IYES) GOTO 99 ! motorcycle daily diurnal calculation is done 
C 
C 
C     determine the model year range index for the 
C     diurnal parameters 
      MYRNG=NO_SPECL_CAS 
C 
      MODEL_YR_RANG : SELECT CASE (MY) 
C 
          CASE (:1971) MODEL_YR_RANG 
            IMY = 1        !need to compare with 1972 emission 
            IFUEL_NUM2 = 1 !use carbureted vehicles only 
            IPPS_NUM=3     !reset for heavy duty vehicles 
            MYRNG = PRE1972 
            DIU1_PRE72 = COLDV(1)+COLDV(2)*VPPR+COLDV(3)*AGE 
            IF(DIU1_PRE72 <= 0.00001) DIU1_PRE72=0.0 
            IF(IVGAS==4 .OR. IVGAS==5) THEN ! heavy duty 2b and 3 factors 
              HEAVY_FCTR=HDWGT(1)           ! for 1972-1978 ldgt3&4 
            ENDIF 
C 
          CASE (1972:1979) MODEL_YR_RANG 
            IMY = 1       !1972-1979 daily diurnal parameters 
            IFUEL_NUM2=1  !use carbureted vehicles only 
            IPPS_NUM=3    !reset for heavy duty vehicles 
            IF(MY/=1979 .AND. (IVGAS==4 .OR. IVGAS==5)) THEN ! heavy duty 2b and 3 factors 
              HEAVY_FCTR=HDWGT(1)                            ! for 1972-1978 ldgt3&4 
            ENDIF 
C 
          CASE (1980:1985) MODEL_YR_RANG 
            IMY = 2       !1980-1985 daily diurnal parameters 
            IF(IVGAS>5 .AND. MY==1985) THEN 
              IPPS_NUM=3    !reset for heavy duty vehicles 
              IFUEL_NUM2=2  !reset for heavy duty vehicles 
            ENDIF 
C 
          CASE (1986:1995) MODEL_YR_RANG 
            IMY = 3       !1986-1995 daily diurnal parameters 
            IPPS_NUM=3    !reset for heavy duty vehicles 
            IFUEL_NUM2=2  !reset for heavy duty vehicles 
C 
          CASE (1996:1998) MODEL_YR_RANG 
            MYRNG = BTWEEN_96_98  !mix of IMY=4 and IMY=3 for light duty vehicles 
            IMY = 3 
            IPPS_NUM=3    !reset for heavy duty vehicles 
            IFUEL_NUM2=2  !reset for heavy duty vehicles 
C 
          CASE (1999:) MODEL_YR_RANG 
            IMY = 4       !1999 daily diurnal parameters 
            IPPS_NUM=3    !reset for heavy duty vehicles 
            IFUEL_NUM2=2  !reset for heavy duty vehicles 
C 
C           Tier 2 Modification 
            IF(TIER2_FLAG == IYES) THEN 
              EVP_WT(3)=EVP_WT(3)*TR2FCTR(MY,ACTUAL_MY,IVGAS,1) !apply tier2 reductions 
              EVP_WT(2)=EVP_WT(2)*TR2FCTR(MY,ACTUAL_MY,IVGAS,2) !apply tier2 reductions 
              EVP_WT(1)=EVP_WT(1)*TR2FCTR(MY,ACTUAL_MY,IVGAS,3) !apply tier2 reductions 
            ENDIF 
C 
      END SELECT MODEL_YR_RANG 
C 
C     Calculate Full-Day's Diurnal for the different types of fuel delivery 
C     categories 
      DO IPPS = 1,IPPS_NUM  !IPPS mappings : fail press(IPPS=1) 
                            !                fail purge(IPPS=2) 
                            !                pass both (IPPS=3) 
        DO IFUEL= IFUEL_NUM1,IFUEL_NUM2  !IFUEL mappings : carb(IFUEL=1) 
                                         !                 FI(IFUEL=2) 
C 
          IFUEL_SYS = IPPS + I_0OR3(IFUEL)   ! 1: carbureted     fail pressure 
                                             ! 2:   "            fail purge 
                                             ! 3:   "            pass both 
                                             ! 4: fuel injection fail pressure 
                                             ! 5:   "     "      fail purge 
                                             ! 6:   "     "      pass both 
          ALPHA = COLDVA(IMY, IPPS, IFUEL) 
     *                   + COLDVB(IMY, IPPS, IFUEL)*VPPR 
     *                   + COLDVC(IMY, IPPS, IFUEL)*(VPPR**2)/1000 
     *                   + COLDVD(IMY, IPPS, IFUEL)*(VPPR**3)/1000000 
          DIU1(IFUEL_SYS) = ALPHA 
          IF(DIU1(IFUEL_SYS).LE.0.00001) DIU1(IFUEL_SYS)=0.0 
C 
          IF(MYRNG == PRE1972) THEN 
C           check pre-1972 parametrization against 1972 carb. fail pressure result 
            DIU1(IFUEL_SYS) = AMAX1(DIU1_PRE72,DIU1(IFUEL_SYS)) 
          ENDIF 
C 
          IF(MYRNG == BTWEEN_96_98) THEN  !for light duty vehicles only; see above 
            ALPHA = COLDVA(4, IPPS, IFUEL) 
     *                   + COLDVB(4, IPPS, IFUEL)*VPPR 
     *                   + COLDVC(4, IPPS, IFUEL)*(VPPR**2)/1000 
     *                   + COLDVD(4, IPPS, IFUEL)*(VPPR**3)/1000000 
            DIU2(IFUEL_SYS) = ALPHA 
            IF(DIU2(IFUEL_SYS).LE.0.00001) DIU2(IFUEL_SYS)=0.0 
C 
          ENDIF 
C 
        END DO !fuel system types 
      END DO   !pass/fail purge/pressure types 
C 
C         After the last calculation, check to see that the diurnal 
C         emissions of vehicles failing pressure is less than or equal to the 
C         diurnal emissions of gross liquid leakers. Both arrays containing 
C         diurnal emissions are checked and capped. 
C 
      IF(DIU1(1).GT.DIU1(7)) DIU1(1)=DIU1(7)   !Check carb failing pressure 
      IF(DIU1(4).GT.DIU1(7)) DIU1(4)=DIU1(7)   !Check fi failing pressure 
      IF(DIU2(1).GT.DIU2(7)) DIU2(1)=DIU2(7)   !Check carb failing pressure 
      IF(DIU2(4).GT.DIU2(7)) DIU2(4)=DIU2(7)   !Check fi failing pressure 
C 
C         After the last calculation, check to see that the diurnal 
C         emissions of vehicles failing purge is less than or equal to the 
C         diurnal emissions of pressure failing vehicles. Both arrays containing 
C         diurnal emissions are checked and capped. 
C 
      IF(DIU1(2).GT.DIU1(1)) DIU1(2)=DIU1(1)   !Check carb failing purge 
      IF(DIU1(5).GT.DIU1(4)) DIU1(5)=DIU1(4)   !Check fi failing purge 
      IF(DIU2(2).GT.DIU2(1)) DIU2(2)=DIU2(1)   !Check carb failing purge 
      IF(DIU2(5).GT.DIU2(4)) DIU2(5)=DIU2(4)   !Check fi failing purge 
C 
C         After the last calculation, check to see that the diurnal 
C         emissions of passing vehicles is less than or equal to the 
C         diurnal emissions of purge fail vehicles.  Both arrays containing 
C         diurnal emissions are checked and capped. 
C 
      IF(DIU1(3).GT.DIU1(2)) DIU1(3)=DIU1(2)   !Check carb passing 
      IF(DIU1(6).GT.DIU1(5)) DIU1(6)=DIU1(5)   !Check fi passing 
      IF(DIU2(3).GT.DIU2(2)) DIU2(3)=DIU2(2)   !Check carb passing 
      IF(DIU2(6).GT.DIU2(5)) DIU2(6)=DIU2(5)   !Check fi passing 
C 
      IF(MY<1980 .AND. IVGAS>=4) THEN    !for heavy duty vehicles find 
        DIU1(1) = AMAX1(DIU1(1),DIU1(2)) !the largest diurnal emission 
        DIU1(1) = AMAX1(DIU1(1),DIU1(3)) !amongst the pressure/purge 
      ENDIF                              !test categories 
C 
      SELECT CASE (MY) 
C 
        CASE (:1971) 
C 
          DIU1(2:3) = DIU1(1) ! carbureted fail pressure 
                              ! parameters for all vehicles 
          DIURNL24A(1) = EVP_WT(1)*FDS_WT(1)*LESS_THN10_FCTR* 
     *                    HEAVY_FCTR*DIU1(1) 
          DIURNL24A(2) = EVP_WT(2)*FDS_WT(1)*LESS_THN10_FCTR* 
     *                    HEAVY_FCTR*DIU1(2) 
          DIURNL24A(3) = EVP_WT(3)*FDS_WT(1)*LESS_THN10_FCTR* 
     *                    HEAVY_FCTR*DIU1(3) 
          DIURNL24A(4) = 0.0 
          DIURNL24A(5) = 0.0 
          DIURNL24A(6) = 0.0 
          DIURNL24A(7) = LESS_THN10_FCTR*EVP_WT(4)*DIU1(7) 
          DIURNL24B    = 0.0 
C 
        CASE (1972:1979) 
          IF((IVGAS>=4.AND.MY/=1979).OR.IVGAS>=6) DIU1(2:3) = DIU1(1) ! carbureted fail pressure parameters 
                                                                      ! for all heavy duty vehicles 
                                                                      ! ldgt3&4 were considered heavy duty trucks 
                                                                      ! for model years earlier than 1979 
          DIURNL24A(1) = EVP_WT(1)*FDS_WT(1)*LESS_THN10_FCTR* 
     *                    HEAVY_FCTR*DIU1(1) 
          DIURNL24A(2) = EVP_WT(2)*FDS_WT(1)*LESS_THN10_FCTR* 
     *                    HEAVY_FCTR*DIU1(2) 
          DIURNL24A(3) = EVP_WT(3)*FDS_WT(1)*LESS_THN10_FCTR* 
     *                    HEAVY_FCTR*DIU1(3) 
          DIURNL24A(4) = 0.0 
          DIURNL24A(5) = 0.0 
          DIURNL24A(6) = 0.0 
          DIURNL24A(7) = LESS_THN10_FCTR*EVP_WT(4)*DIU1(7) 
          DIURNL24B    = 0.0 
C 
        CASE (1996:1998) 
          DO IFUEL_SYS =1,6 
            IPPS=IFUEL_SYS 
            IFUEL=1 
            IF(IFUEL_SYS > 3) THEN 
              IPPS=IFUEL_SYS-3 
              IFUEL=2 
            ENDIF 
            DIURNL24A(IFUEL_SYS) = LESS_THN10_FCTR*HEAVY_FCTR* 
     *                  EVP_WT(IPPS)*FDS_WT(IFUEL)*DIU1(IFUEL_SYS) 
C 
            DIURNL24B(IFUEL_SYS) = LESS_THN10_FCTR*HEAVY_FCTR* 
     *                    EVP96_98(IPPS)*FDS_WT(IFUEL)*DIU2(IFUEL_SYS) 
          END DO 
C         liquid leaks: 
          DIURNL24A(7) = LESS_THN10_FCTR*EVP_WT(4)*DIU1(7) 
          DIURNL24B(7) = DIURNL24A(7) 
C 
        CASE DEFAULT 
          IF(IVGAS>5 .AND. IMY==2 .AND. MY/=1985) DIU1(2:3) = DIU1(1) ! 1980-1985 heavy duty vehicles 
          DO IFUEL_SYS =1,6 
            IPPS=IFUEL_SYS 
            IFUEL=1 
            IF(IFUEL_SYS > 3) THEN 
              IPPS=IFUEL_SYS-3 
              IFUEL=2 
            ENDIF 
            DIURNL24A(IFUEL_SYS) = EVP_WT(IPPS)*FDS_WT(IFUEL)* 
     *                        LESS_THN10_FCTR*HEAVY_FCTR*DIU1(IFUEL_SYS) 
C 
          END DO                 
C         liquid leaks: 
          DIURNL24A(7) = EVP_WT(4)*LESS_THN10_FCTR*DIU1(7) 
          DIURNL24B    = 0.0 
C 
      END SELECT 
C 
C     Adjust if the user requests a high altitude region. 
  99  IF(IREJN==2) THEN 
        DIURNL24A(1:6) = HIADJ*DIURNL24A(1:6) 
        DIURNL24B(1:6) = HIADJ*DIURNL24B(1:6) 
      ENDIF 
C 
      RETURN 
      END SUBROUTINE DIRNLDAY 
