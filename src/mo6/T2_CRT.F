      SUBROUTINE T2_CRT(ICY) 
C 
C 
C     This subroutine builds the T2CRT50K array used in EF_DSPTCH to 
C     compute the Tier2 LDDV and LDDTs basic emission rates and in 
C     EF_TIER2 to compute the Tier2 LDGVs/LDGTs basic emission rates. 
C 
C 
C     Called by:  EFCALX. 
C 
C 
C     CALLS:  none 
C 
C 
C  Changes: (Last change first) 
C 
C  02 Feb 00 @DynTel-ddj 1-032, New subroutine. 
C 
C 
C  Input on call or from calls: 
C 
C    argument list: ICY 
C 
C 
C  common blocks: 
C    /TIER2/  TIER2_50K (ISTD, IP, IV) 
C             T2_BIN_FRAC(IMY, IBIN, IP, IVTEMP) 
C 
C  Output on return: 
C    /TIER2/  T2CRT50K(IV, IP, JDX) 
C 
C 
C    common blocks: 
C 
C  Local variable / array dictionary: 
C 
C   Name   Type              Description 
C  ------  ----  ----------------------------------------------------- 
C 
C 
C 
C  Notes: 
C 
C 
      IMPLICIT NONE 
C 
      INCLUDE 'MAXIMA.I' 
      INCLUDE 'TIER2.I' 
C 
      INTEGER, INTENT(IN)      :: ICY 
C 
      INTEGER JDX 
      INTEGER IBIN 
      INTEGER IP 
      INTEGER IMY 
      INTEGER IVTEMP 
      INTEGER MY 
C 
      REAL CRT50K 
C 
      DO IVTEMP = 1, 5 ! 1,LDV/LDDV, 2 LDT1, 3 LDT2/LDDT12, 4 LDT3, 5 LDT4/LDDT34 
         DO IP = 1, 3 
            DO JDX = 1, MAXYRS 
               MY = ICY - JDX + 1 
               IF(MY .GT. 2004) THEN 
                  IF(MY .LT. 2015) THEN 
                     IMY = MY - 2003 
                  ELSE 
                     IMY = 12 
                  END IF 
               ELSE 
                  IMY = 1 
               END IF 
C 
               CRT50K = 0.0 
C 
C      M6.EXH.011 states that for LDDT12 emission standards, use 
C      a 50/50 weighting of LDT1 and LDT2.  However, EPA has 
C      directed for task 1-032 to phase in standards using LDT2 
C      for LDDT12. 
C 
               DO IBIN = 1, 12 
                  CRT50K = CRT50K + TIER2_50K(IVTEMP,IBIN,IP)* 
     &                              T2_BIN_FRAC(IMY,IBIN,IP,IVTEMP) 
               END DO 
               T2CRT50K(IVTEMP,IP,JDX) = CRT50K 
            END DO 
         END DO 
      END DO 
C 
      RETURN 
      END SUBROUTINE 
