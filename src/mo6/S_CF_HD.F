      REAL FUNCTION S_CF_HD(MY,IP,IV,STD,AX_PPM) 
C 
C     S_CF computes the sulfur fuel correction factor for the HDGVs 
C     using the passed gasoline sulfur content (AX_PPM). 
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
C            MYCODE.I; MY passed through parameter list 
C     03 Jul 00 @EPA-ddj Bug157 Initialized SULFUR to 0.0 rather than 1.0 
C                       (toward end of subroutine) 
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
C         'SULFUR.I'   X_PPM, SUL_FLAG 
C 
C     Argument List:  MY, IP, IV, STD, AX_PPM 
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
      INTEGER, INTENT(IN)      :: STD 
      INTEGER, INTENT(IN)      :: MY 
C 
      INTEGER I 
      INTEGER IPI 
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
      REAL, DIMENSION(3,3), SAVE  :: S_COMP = RESHAPE( 
     &                    (/ 0.055,    0.0308,   0.076, 
     &                       7.22E-04, 3.02E-04, 6.30E-04, 
     &                       0.09,     0.146,    0.151 / ), 
     &                       (/3,3/) ) 
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
C     Initialize Sulfur Correction 
C 
      SULFUR = 0.0 
C 
      IF (VVHDGV(IV).EQ.1) THEN 
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
          IF (MY .GE. 2005) THEN 
C 
             DO I = 1, 2 
                EX(I) = EXP(S_COMP(IP,3)*ALOG(PPM(I))) 
             END DO 
C 
          ELSEIF(MY .GE. 1998) THEN 
             IF(IP .LT. 3) THEN     !Use Tier1 composite for 3-way plus 
                DO I = 1, 2         !Tier0 ccomposite for oxidation cat. 
                   EX(I) = P3W*EXP(S_COMP(IP,2)*PPM(I) ) + 
     &                     POX*EXP(S_COMP(IP,1)*ALOG(PPM(I)) ) 
                END DO 
             ELSE                   !If NOx, then Tier1 for 3-way only. 
                DO I = 1, 2 
                   EX(I) = P3W*EXP(S_COMP(IP,2)*PPM(I) ) 
                END DO 
             END IF 
C 
          ELSE    ! pre-1998 HDGVs 
             IF(IP .LT. 3) THEN     !Use Tier0 composite for both 3-way 
                DO I = 1, 2         !and oxidation cat. (HC and CO). 
                   EX(I) = (P3W + POX)*EXP(S_COMP(IP,1)*ALOG(PPM(I)) ) 
                END DO 
             ELSE                   !Use Tier0 composite for NOx 3-way 
                DO I = 1, 2         !only. 
                   EX(I) = P3W*EXP(S_COMP(IP,2)*ALOG(PPM(I))) 
                END DO 
             END IF 
          END IF 
C 
          IF(EX(2) .GT. 0.0) THEN 
            IF (JSTD.EQ.2.AND.X_PPM.GT.330.0) THEN 
C 
C     If TIER1 and sulfur content greater than 330 ppm, apply 
C     TIER0 correction (p.15, M6.FUL.001, 1/5/99). 
C 
C     Compute correction to normal TIER1 using TIER0 normal emitters 
C     with 330 ppm. 
C 
              AX330(2) = EXP(S_COMP(IP,1) * ALOG( PPM(2) )) 
              AX(2) = AX330(2) 
              AX30(2) = AX330(2) 
C 
              AX330(1) = EXP(S_COMP(IP,1) * ALOG( PPM(1) )) 
              AX(1) = EXP(S_COMP(IP,1) * ALOG(X_PPM) ) 
              AX30(1) = EXP(S_COMP(IP,1) * ALOG(30.0) ) 
C 
              SULFUR = (( CORR_S(AX(1),AX(2)) - 
     &                    CORR_S(AX30(1),AX30(2)) ) / 
     &                  ( CORR_S(AX330(1),AX330(2)) - 
     &                    CORR_S(AX30(1),AX30(2)) )) 
     &                * ( CORR_S(EX(1),EX(2)) ) 
C 
            ELSE 
              SULFUR = CORR_S(EX(1),EX(2)) 
            END IF 
          END IF 
        ELSE 
C 
C     User Supplied Sulfur Correction Factors (USSCF) applied in 
C     place of the regression coefficients. 
C 
C     Set pollutant code: 
C          IPI = 1, HC 
C                2, NMHC 
C                3, CO 
C                4, NOx 
C 
          IPI = IP + 1 
C 
             SULFUR = (USSCF(JSTD,1,1,IPI,3) - 1.0) * (POX + P3W) 
        END IF                ! IF SUL_FLAG 
      END IF                  ! IF IV not VTMC 
C 
      S_CF_HD = SULFUR 
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
