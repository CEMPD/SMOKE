      REAL FUNCTION EF_TIER2(IV,IP,MY,JDX,INH,EF_T1) 
C 
C  TIER2 computes the weighted Tier2 basic emission rate. 
C  Initially only the light duty gasoline classes are affected. 
C  This may change at a later date to include the heavy duty 
C  classes. 
C 
C 
C  Called by EF_DSPTCH. 
C 
C 
C  Calls None. 
C 
C 
C  Changes :  
C 
C   8 Jun 00 @EPA-ddj Bug73, Moved CRT50 ratio for TIER2 emissions 
C            to EF_TIER2. 
C  15 Mar 00 @DynTel-ddj 1-032  Implemented array T2CRT50K. 
C  30 Dec 99 @DYNTEL-ddj 1-021  New Function. 
C 
C 
C  Input on call: 
C 
C    Constants:  
C 
C 
C    argument list: IV, IP, MY, EF_T1 
C 
C    common blocks: 
C    /TIER2/    TIER2_50K, T2_BIN_FRAC 
C    /FLAGS5/   TIER2_FLAG 
C    /LDGBER/   CRT050K 
C 
C 
C  Output on return: 
C 
C    EF_TIER2     Tier2 basic emission rate for vehicle type IV, 
C                 pollutant IP and model year MY. 
C 
C    common blocks: 
C 
C 
C  Local array subscripts: 
C   TIER2_50K     (IV, IBIN, IP) 
C   T2_BIN_FRAC   (JDX, IBIN, IP, IV) 
C   CRT50K        (ISTD, IP, IV) 
C 
C 
C  Local variable / array dictionary: 
C 
C   Name         Type              Description 
C  ------        ----  ------------------------------------------------------- 
C  IBIN           I    Tier2 bin index. 
C  TIER2_50K      R    Tier2 50K certification standards. 
C  T2_BIN_FRAC    R    Tier2 bin phase-in fraction. 
C 
C  Notes: 
C 
C 
      IMPLICIT NONE 
C 
      INCLUDE 'FLAGS5.I' 
      INCLUDE 'LDGBER.I' 
      INCLUDE 'IVTYPE.I' 
      INCLUDE 'TIER2.I' 
C 
      INTEGER, INTENT(IN)   :: JDX 
      INTEGER, INTENT(IN)   :: INH 
      INTEGER, INTENT(IN)   :: IV 
      INTEGER, INTENT(IN)   :: IP 
      INTEGER, INTENT(IN)   :: MY 
      REAL,    INTENT(IN)   :: EF_T1 
C 
      INTEGER, PARAMETER :: NORMAL=1 
      INTEGER, PARAMETER :: HC=1, CO=2, NOX=3 
      INTEGER, PARAMETER :: TIER0=1, TIER1=3 
C 
      REAL    :: C50K 
C 
      IF (TIER2_FLAG .EQ. 1 .AND. MY .GT. 2003) THEN 
C 
C     Tier2 50K certification standard (T2CRT50K) is the weighted sum 
C     of the 12 individual Tier2 bins.  T2CRT50K is computed in 
C     subroutine T2_CRT. 
C 
           IF(IP.EQ.HC) THEN 
             C50K = T2CRT50K(IV,HC,JDX) / CRT050K(TIER0,HC,VTLDGV) 
           ELSE IF(IP.EQ.CO) THEN 
             C50K = T2CRT50K(IV,CO,JDX) / CRT050K(TIER1,CO,IV) 
           ELSE           ! IP = NOx 
             C50K = T2CRT50K(IV,NOX,JDX) / CRT050K(TIER1,NOX,VTLDGV) 
           END IF 
C 
         IF (INH .EQ. NORMAL) THEN 
C 
C    NORMAL Tier2 emission factor obtained by multiplying the Tier1 
C    emission factors by the ratio of the weighted Tier2 to Tier1 50K 
C    certification standards. 
C 
            EF_TIER2 = EF_T1 * C50K 
C 
         ELSE  ! HIGH Emitter 
C 
C    HIGH Tier2 emission factor equal to the AVERAGE of the Tier1 
C    emission factor multiplied by the ratio of the 50K Certification 
C    Standards AND the original Tier1 emission factor. 
C 
           EF_TIER2 = (EF_T1*C50K + EF_T1) / 2.0 
         END IF 
C 
      ELSE ! TIER2 standard not applicable 
         EF_TIER2 = EF_T1 
      END IF 
C 
      RETURN 
      END 
