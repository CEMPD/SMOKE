      REAL FUNCTION S_CF_LD(MY,IP,IV,INH,ISR,STD,AX_PPM) 
C 
C     S_CF computes the fuel correction factor for the passed gasoline 
C     sulfur content (AX_PPM). 
C 
C 
C     Called by SUL_COR 
C 
C 
C     Calls CORR_S and P3WPOX 
C 
C 
C  Changes: (Last change first) 
C 
C  04 Dec 00 @EPA-bag Task X6, removed MY, JDX, and IDX from MYCODE.I; removed 
C            MYCODE.I and added model year, MY, to the parameter list of P3WPOX 
C     30 Aug 00 @EPA-djb Bug 242, Set IVI same for LDGV and LDGT1. 
C     28 Jun 00 @EPA-ddj initialized SULFUR to 0.0 rather than 1.0 
C     21 Jun 00 @EPA-ddj Bug157 Split correction into LDGVs (S_CF_LD) 
C               and HDGVs (S_CF_HD). 
C     07 Jun 00 @EPA-elg NOx High emitter sulfur correction for Tier1 
C               and later changed.  60% of Normal correction is applied to the 
C               final High correction rather than separately to the numerator 
C               and denominator of the final High correction factor. 
C     25 May 00 @EPA-ddj Bug 108, Changed normalization of sulfur 
C               correction from 50 ppm to 30 ppm. 
C     20 Apr 00 @EPA-ddj E-7  Changed name to S_CF to conform to 8 
C               character naming convention. 
C     28 Oct 99 @DynTel-ddj #1-003, Set SULFUR equal to zero for 
C               vehicle type VTMC (i.e. motorcycles). 
C     21 Jul 99 @DynTel-ddj #2-694, Changes necessary to embed the IH 
C               loop inside the JSTD loop.  Also fixed HC correction 
C               to always use NMHC correction.  THC correction is 
C               retained in code, in case EPA has a change of heart, 
C               but is bypassed in calculation. 
C      5 May 99 @DynTel-ddj #2-699, New function 
C 
C 
C   Input on CALL: 
C 
C     from INCLUDE file 
C         'BASEQ9.I'   LEVIMP(STD,IYR,IV) 
C         'IVTYPE.I'   VTMC 
C         'SULFUR.I'   X_PPM, SUL_FLAG, S_COEF(IPI,ISR,JSTD,INH) 
C 
C     Argument List:  MY, IP, IV, INH, ISR, STD, AX_PPM 
C 
C 
C   Output on RETURN: 
C 
C     SULFUR_CF, fuel correction factor for gasoline sulfur content. 
C 
C 
C   Local Array Table: 
C 
C     AX(I)    Array containing input to CORR_S for TIER0 
C              and relevant sulfur content. 
C     AX330(I) Array containing input to CORR_S for TIER0 and 
C              330 ppm sulfur content. 
C     EX(I)    Array containing input to CORR_S for relevant standard 
C              and sulfur content. 
C     PPM(I)   Array containing the gasoline sulfur content in ppm for: 
C                 I == 1, scenario gasoline sulfur content. 
C                   == 2, correction factor based sulfur level (30 ppm). 
C 
C 
C 
      IMPLICIT NONE 
C 
      INCLUDE 'BASEQ9.I' 
      INCLUDE 'IVTYPE.I' 
      INCLUDE 'SULFUR.I' 
      INCLUDE 'IOUCOM.I' 
C 
      INTEGER, INTENT(IN)      :: IP 
      INTEGER, INTENT(IN)      :: IV 
      INTEGER, INTENT(IN)      :: INH 
      INTEGER, INTENT(IN)      :: ISR 
      INTEGER, INTENT(IN)      :: STD 
      INTEGER, INTENT(IN)      :: MY 
      INTEGER I 
      INTEGER IPI 
      INTEGER IVI 
      INTEGER JSTD 
C 
      REAL SULFUR 
      REAL P3W 
      REAL POX 
      REAL X_PPM 
      REAL AX_PPM 
C 
      REAL, DIMENSION(2)    :: AX 
      REAL, DIMENSION(2)    :: AX30 
      REAL, DIMENSION(2)    :: AX330 
      REAL, DIMENSION(2)    :: EX 
      REAL, DIMENSION(2)    :: PPM 
C 
C  Compute the correction factor for gasoline sulfur content (X_PPM). 
C 
      CALL P3WPOX(MY,IV,P3W,POX) 
C 
      IF (IP.EQ.3) THEN 
C 
C      Sulfur correction factor not used for NON catalyst vehicles and 
C      NOx emissions not adjusted for oxidation catalyst vehicles 
C 
        POX = 0.0 
      END IF 
C 
      IF (IV.EQ.VTLDGV .OR. IV.EQ.VTLDGT1) THEN 
         IVI = 1 
      ELSE 
         IVI = 2 
      END IF 
C 
C     Set pollutant code: 
C          IPI = 1, HC 
C                2, NMHC 
C                3, CO 
C                4, NOx 
C 
      IPI = IP + 1 
C 
C     Per 7/21/99 meeting with EPA, the non-methane hydrocarbons NMHC 
C     (IPI = 2) are to be used.  Differentiation between THC (IPI=1) 
C     and NMHC will not be necessary. 
C 
C     Array PPM(I) contains the gasoline sulfur content in ppm for: 
C           I == 1, gasoline sulfur content. 
C             == 2, correction factor based sulfur level (30 ppm). 
C 
C     This code sets a maximum of 600 ppm Sulfur.  This was the extent 
C     of the underlying data. 
C 
      X_PPM  = AX_PPM 
C 
      IF(AX_PPM .GE. 600.) THEN 
        X_PPM = 600. 
      ENDIF 
C 
      PPM(1) = X_PPM 
      PPM(2) = 30.0 
C 
      IF(STD .GT.4) THEN 
         JSTD = 3 
      ELSE 
         JSTD = STD 
      END IF 
C 
C      Sulfur correction factor not used for NON catalyst vehicles 
C      (i.e. motorcycles). 
C 
      SULFUR = 0.0 
C 
      IF (VVGAS(IV).EQ.1) THEN   !No sulfur correction for motorcycles 
C 
C SUL_FLAG       Flag indicating source of sulfur correction factors: 
C                    1 Use Default built-in Sulfur Levels. 
C                    2 User supplied Average and Maximum Sulfur Levels. 
C                    3 User supplied Average and Maximum Sulfur Levels 
C                       and sulfur correction factors. 
C 
        IF (SUL_FLAG.EQ.1 .OR. SUL_FLAG.EQ.2) THEN 
C 
C     For post 1981 model years. 
C 
          IF (MY.GT.1980) THEN 
             IF (STD.LT.6) THEN   !Not a ZEV all other cases computed. 
C 
C    Compute the Sulfur correction factor by vehicle standard. 
C 
               SELECT CASE (JSTD) 
C 
C   *******************  TIER0 SULFUR CORRECTION  ******************* 
C 
                 CASE (1) 
                   IF (INH.EQ.1) THEN 
C 
C     For TIER0 normal emitters. 
C 
                     DO I = 1, 2 
                       EX(I) = EXP(S_COEF(IPI,ISR,JSTD,INH,IVI) 
     &                               *ALOG(PPM(I))) 
                     END DO 
                   ELSE 
C 
C     For TIER0 high emitters. 
C 
                     DO I = 1, 2 
                       EX(I) = EXP(S_COEF(IPI,ISR,JSTD,INH,IVI) 
     &                               *PPM(I)) 
                     END DO 
                   ENDIF 
C 
C   *******************  TIER1 SULFUR CORRECTION  ******************* 
C 
                 CASE (2) 
C 
C     For TIER 1, NORMAL emitters are a log-linear regression, 
C     HC and CO HIGH emitters us TIER 0 log-linear regression, 
C     and NOx uses 60% of NORMAL emitters. 
C 
                   IF (X_PPM.GT.330.0) THEN 
C 
C     If sulfur content greater than 330 ppm, TIER1 uses an 
C     adjustment based on TIER0. 
C 
                     PPM(1) = 330.0 
                     IF (INH.EQ.1) THEN 
C 
C     Compute correction to normal TIER1 using TIER0 normal emitters 
C     with 330 ppm. 
C 
                         AX330(2) = EXP(S_COEF(IPI,ISR,1,INH,IVI) 
     &                                    *ALOG( PPM(2) )) 
                         AX(2) = AX330(2) 
                         AX30(2) = AX330(2) 
C 
                         AX330(1) = EXP(S_COEF(IPI,ISR,1,INH,IVI) 
     &                                    *ALOG( PPM(1) )) 
                         AX(1) = EXP(S_COEF(IPI,ISR,1,INH,IVI) 
     &                                 * ALOG(X_PPM) ) 
                         AX30(1) = EXP(S_COEF(IPI,ISR,1,INH,IVI) 
     &                                 *ALOG(30.0) ) 
                     ELSE 
C 
C     Compute TIER1 high emitter correction using TIER0 high emitters 
C     with 330 ppm. 
C 
                         AX330(2) = EXP(S_COEF(IPI,ISR,1,INH,IVI) 
     &                                     *PPM(2) ) 
                         AX(2) = AX330(2) 
                         AX30(2) = AX330(2) 
C 
                         AX330(1) = EXP(S_COEF(IPI,ISR,1,INH,IVI) 
     &                                    *PPM(1) ) 
                         AX(1) = EXP(S_COEF(IPI,ISR,1,INH,IVI) 
     &                                 *X_PPM ) 
                         AX30(1) = EXP(S_COEF(IPI,ISR,1,INH,IVI)*30.0) 
                     END IF   ! NORMAL AND HIGH CONTROL 
                   END IF 
C 
C     Compute TIER1 gasoline sulfur content correction factor. 
C 
                   IF (INH.EQ.1) THEN 
                     DO I = 1, 2 
                       EX(I) = EXP(S_COEF(IPI,ISR,JSTD,INH,IVI) 
     &                               *PPM(I) ) 
                     END DO 
                   ELSE 
                     IF (IPI.NE.4) THEN 
                       DO I = 1, 2 
                         EX(I) = EXP(S_COEF(IPI,ISR,1,INH,IVI) 
     &                                 *PPM(I) ) 
                       END DO 
                     ELSE 
                       DO I = 1, 2 
                         EX(I) = EXP(S_COEF(IPI,ISR,JSTD,1,IVI) 
     &                                   *PPM(I) ) 
                       END DO 
                     END IF 
                   END IF 
C 
C     Reset gasoline sulfur content 
C 
                   PPM(1) = X_PPM 
C 
                 CASE DEFAULT ! APPLIES to LEV, ULEV, and 
                              ! Tier2 Standards. 
C 
C     For LEV and ULEV, NORMAL emitters are a log-log regression, 
C     HC and CO HIGH emitters us TIER 0 log-linear regression, 
C     and NOx uses 60% of NORMAL emitters. 
C 
                   IF (INH.EQ.1) THEN 
                     DO I = 1, 2 
                       EX(I) = EXP(S_COEF(IPI,ISR,JSTD,INH,IVI) 
     &                               *ALOG(PPM(I))) 
                     END DO 
                   ELSE 
                     IF (IPI.NE.4) THEN 
                       DO I = 1, 2 
                         EX(I) = EXP(S_COEF(IPI,ISR,1,INH,IVI) 
     &                                 *PPM(I) ) 
                       END DO 
                     ELSE 
                       DO I = 1, 2 
                         EX(I) = EXP(S_COEF(IPI,ISR,JSTD,1,IVI) 
     &                               *ALOG(PPM(I))) 
                       END DO 
                     END IF 
                   END IF 
C 
               END SELECT 
C 
C          ***************  END OF CASE STATEMENT  ************** 
C 
               IF (JSTD.EQ.2.AND.X_PPM.GT.330.0) THEN 
C 
C     If TIER1 and sulfur content greater than 330 ppm, apply 
C     TIER0 correction (p.15, M6.FUL.001, 1/5/99). 
C 
                 SULFUR = (( CORR_S(AX(1),AX(2)) - 
     &                       CORR_S(AX30(1),AX30(2)) ) / 
     &                     ( CORR_S(AX330(1),AX330(2)) - 
     &                       CORR_S(AX30(1),AX30(2)) )) 
     &                     * CORR_S(EX(1),EX(2)) * (P3W + POX) 
                 IF(INH.EQ.2 .AND. IP.EQ.3 .AND. JSTD.GE.2) THEN 
                   SULFUR = SULFUR*0.60 
                 ENDIF 
C 
               ELSE 
                 SULFUR = CORR_S(EX(1),EX(2)) * (P3W + POX) 
                 IF(INH.EQ.2 .AND. IP.EQ.3 .AND. JSTD.GE.2) THEN 
                   SULFUR = SULFUR*0.60 
                 ENDIF 
               END IF 
C 
             END IF           ! STD less than 6 IF statement. 
C 
          ELSE                ! ELSE the pre-1981 vehicles are computed 
C 
C     For pre-1981 vehicles with Oxidation cataysts, the TIER0 
C     sulfur effects for HIGH emitters is applied to HC and CO. 
C 
            IF (IP.LT.3) THEN 
C 
C      NOx emissions not adjusted for oxidation catalyst vehicles. 
C      (EPA memo from D. Brzezinski 1/27/99) 
C 
              DO I = 1, 2 
                 EX(I) = EXP(S_COEF(IPI,ISR,1,2,IVI)*PPM(I) ) 
              END DO 
              SULFUR = CORR_S(EX(1),EX(2))*POX 
            END IF 
          END IF              ! End IF for Pre-Post 1981 vehicles 
        ELSE                  ! IF SUL_FLAG 
C 
C     User Supplied Sulfur Correction Factors (USSCF) applied in 
C     place of the regression coefficients. 
C 
           SULFUR = (USSCF(JSTD,ISR,INH,IPI,IVI) - 1.0) * (POX + P3W) 
        END IF                ! IF SUL_FLAG 
      END IF                  ! IF IV not VTMC 
C 
      S_CF_LD = SULFUR 
C 
      RETURN 
C 
      CONTAINS 
         REAL FUNCTION CORR_S(A,B) 
C 
C        Define INTERNAL Function CORR_S 
C        Calculates the 
C 
         REAL A 
         REAL B 
C 
         CORR_S = (A - B) / B 
C 
         END FUNCTION CORR_S 
C 
      END 
