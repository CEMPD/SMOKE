      REAL FUNCTION SUL_COR(IP,STD,ICY,MY,IV,INH,ISR) 
C 
C SULFUR computes the emissions correction factor associated 
C with the sulfur content of gasoline. 
C 
C Changes : 
C 
C  21 Nov 00 @EPA-elg Added No Clean Air Act Sulfur initializations 
C  08 Aug 00 @EPA-EG BUG#216; sulfur content initialization 
C  26 Jun 00 @EPA-EG BUG RFG Phase 1 RFG is always 300 ppm Sulfur; so added 
C                 conditional for phase 1 (ICY.LE.1999 .AND. RFGFLG.EQ.2) 
C  21 Jun 00 @EPA-ddj Bug 157, S_CF.FOR split for heavy duty vehicles, 
C                     S_CF_HD.FOR, and light duty vehicles, S_CF_LD.FOR 
C  12 Jun 00 @EPA-ddj Bug 127, Added conditional to restrict the 
C            application of irreversibility to LDGVs only. 
C  07 Jun 00 @EPA-elg Changed Sulfur Irreversibility code so that it 
C                     is not applied to pre-2004 MYs. 
C  25 May 00 @EPA-ddj Bug 108,  Comments 
C  24 May 00 @EPA-ddj BUG107, Revised assignment of S_PPM from function 
C            of model year to function of calendar year. 
C  16 May 00 @EPS-ddj RFG BUG, Revised selection of RFG sulfur 
C  20 Apr 00 @EPA-ddj Revised function to handle short-term, long-term 
C            and irreversibility effect of sulfur content in gasoline. 
C     28 Oct 99 @DynTel-ddj #1-003, Set SULFUR equal to zero for 
C               vehicle type VTMC (i.e. motorcycles). 
C     21 Jul 99 @DynTel-ddj #2-694, Changes necessary to embed the IH 
C               loop inside the ISTD loop.  Also fixed HC correction 
C               to always use NMHC correction.  THC correction is 
C               retained in code, in case EPA has a change of heart, 
C               but is bypassed in calculation. 
C      5 May 99 @DynTel-ddj #2-699, New function 
C 
C 
C Called by: HRLOOP 
C 
C 
C Calls: S_CF_LD, S_CF_HD and STANDARD. 
C 
C 
C Input on call: 
C 
C   argument list : 
C       IP, STD, ICY, MY, IV, INH, ISR 
C 
C   common blocks : 
C       LONGTERM(IP) 
C 
C Output on return: 
C 
C   argument list : 
C 
C   common blocks : 
C       LONGTERM(IP) 
C 
C   function      : 
C       SULFUR 
C 
C 
C Local variable / array dictionary : 
C Name            Type         Description 
C -----          ------    ----------------------------------------- 
C CURRENT           R      Short term sulfur effect of current sulfur 
C                          content of gasoline for calendar year ICY. 
C MAXIMUM           R      Short term sulfur effect of due to maximum 
C                          sulfur content exposed for model year MY. 
C JSTD              I      Index indicating emissions standard. 
C 
C 
C Notes : 
C 
C 
      IMPLICIT NONE 
C 
      INCLUDE 'FLAGS3.I' 
      INCLUDE 'FLAGS5.I' 
      INCLUDE 'IVTYPE.I' 
      INCLUDE 'RFORM1.I' 
      INCLUDE 'SULFUR.I' 
C 
      INTEGER, INTENT(IN)    :: ICY 
      INTEGER, INTENT(IN)    :: INH 
      INTEGER, INTENT(IN)    :: IP 
      INTEGER, INTENT(IN)    :: ISR 
      INTEGER, INTENT(IN)    :: IV 
      INTEGER, INTENT(IN)    :: MY 
      INTEGER, INTENT(IN)    :: STD 
C 
      REAL CURRENT 
      REAL MAXIMUM 
      REAL, DIMENSION(6)     :: IR =  ! 0.425 for Tier2, 0.15 all 0thers 
     & (/0.15, 0.15, 0.15, 0.15, 0.425, 0.15/) 
C        TIER0,T1INT,TIER1,LEV,TIER2,ZEV 
C 
      INTEGER IYR 
      INTEGER IMY 
      INTEGER IT2 
      INTEGER JSTD 
C 
C     FUNCTION CALLS 
C 
      INTEGER STANDARD 
      REAL S_CF_LD 
      REAL S_CF_HD 
C 
C     Initialize function 
C 
      SUL_COR = 1.0 
C 
C     If Clean Air Act command is set (CAAFLG = 1) 
C 
      IF(CAAFLG.EQ.1 .AND. ICY.GE.CAAMY) THEN 
        IT2 = 1 
      ELSE 
        IT2 = TIER2_FLAG + 1 
      ENDIF 
C 
                                          ! Pre-2000, use default/input 
      IF(ICY .GT. 2015) THEN 
         S_PPM = S_FUEL(1,FUEL_P,16,IT2)  ! Post-2015, use matrix 
      ELSE IF(ICY.GE.2000) THEN 
         IYR = ICY - 1999 
         S_PPM = S_FUEL(1,FUEL_P,IYR,IT2) ! Post-1999, use matrix 
      END IF 
C 
      JSTD = STANDARD(STD)                ! Convert emissions standard 
C 
C   Compute Short-term effect for current sulfur level. 
C 
      IF(VVLDGAS(IV) .EQ. 1) THEN 
         CURRENT = S_CF_LD(MY,IP,IV,INH,ISR,JSTD,S_PPM) 
C 
         IF(JSTD .GE. 3 .AND. JSTD .LE. 5) THEN 
C 
C    Apply long-term affect to LEV and cleaner vehicles. 
C    Normalizing to base level (30 ppm). 
C 
           CURRENT = CURRENT*LONGTERM(IP) 
           IF(MY .GT. 2015) THEN 
              IMY = 16 
           ELSE 
              IF(MY .LT. 2000) THEN 
                 IMY = 1 
              ELSE 
                 IMY = MY - 1999 
              END IF 
           END IF 
C                                   ! MAXIMUM is exposure effect. 
           MAXIMUM = 
     &         S_CF_LD(MY,IP,IV,INH,ISR,JSTD,S_FUEL(2,FUEL_P,IMY,IT2)) 
C 
C     Weight current short-term effect and maximum exposure effect 
C     for each model year. 
C 
           IF(MY.GE.2004) THEN 
              SUL_COR = 1.0 + CURRENT*(1.0 - IR(JSTD) ) + 
     &                        MAXIMUM*IR(JSTD) 
           ELSE 
              SUL_COR = 1.0 + CURRENT 
           ENDIF 
C 
         ELSEIF(JSTD.GE.1 .AND. JSTD.LE.2) THEN 
           SUL_COR = 1.0 + CURRENT    !Current is applied to pre-2000 MY LDGVs. 
         ELSE 
           SUL_COR = 1.0 
         END IF                       !JSTD ENDIF Condition. 
C 
      ELSEIF(VVHDGV(IV) .EQ. 1) THEN 
         CURRENT = S_CF_HD(MY,IP,IV,JSTD,S_PPM) 
         IF(MY .GE. 2005) THEN   !HDG sulfur longterm affect, post-2004 
           CURRENT = CURRENT * LONGTERM(IP) 
         END IF 
         SUL_COR = 1.0 + CURRENT 
      ENDIF 
C 
      RETURN 
      END 
