      SUBROUTINE CALAER 
C 
C  CALAER calculates the effect of an anti-tamping program on 
C  tampering rates and stores the results in the EVPAER AND AER arrays. 
C 
C  Called by GETATP and TAMPER 
C 
C  Changes : 
C 
C  08 Jun 00 AIR Task 03: Removed type declarations for intrinsic functions. 
C  20 Aug 98 @DynTel-ZK 2-000 explicit Typing 
C 
C 
C  Input on call: 
C 
C    common blocks: 
C    /ATPAR3/ AIR,CAT,ECK,TNK,EGR 
C    /ATPAR4/ CAN,PCV,CAP 
C 
C  Output on return: 
C 
C    common blocks: 
C    /EGSCAL/ AER 
C    /EVPDAT/ EVPAER 
C 
C  Local variable dictionary: 
C 
C   Name   Type              Description 
C  ------  ----  ------------------------------------------------------- 
C  AIREGR   R    temporary array to hold either AIR or EGR depending 
C                on IPG 
C  DIFF     R    temporary variable to hold overlap difference 
C  DIFFA    R    temporary variable to hold another overlap difference 
C  DIFFB    R    temporary variable to hold another overlap difference 
C  DIFFC    R    temporary variable to hold another overlap difference 
C 
C  Notes: 
C 
C  None. 
C 
      IMPLICIT NONE 
      INCLUDE 'ATPAR3.I' 
      INCLUDE 'ATPAR4.I' 
      INCLUDE 'EGSCAL.I' 
      INCLUDE 'EVPDAT.I' 
      INTEGER IAY, IPG, NCKIAY 
      REAL AIREGR, DIFF, DIFFA, DIFFB, DIFFC 
      DIMENSION AIREGR(2) 
C 
C 
      DO 56 IAY=1,2 
C 
      DO 54 IPG=1,2 
C 
C  Choose Air Pump or EGR effectiveness on pollutant group (IPG). 
C 
      IF(IPG.EQ.1) AIREGR(IAY)=AIR(IAY) 
      IF(IPG.EQ.2) AIREGR(IAY)=EGR(IAY) 
C 
C  Region 1 calculations 
C  Overlap between Catalyst Removal and Air Pump Tampering or 
C  Overlap between Catalyst Removal and EGR      Tampering 
C 
      DIFF=AIREGR(IAY)-CAT(IAY) 
      IF(DIFF.GE.0.0) GOTO 10 
C 
      AER(1,8,IPG,IAY)=ABS(DIFF) 
      AER(1,1,IPG,IAY)=1.00-CAT(IAY) 
      GOTO 12 
C 
   10 AER(1,9,IPG,IAY)=DIFF 
      AER(1,1,IPG,IAY)=1.00-AIREGR(IAY) 
C 
C  Region 2 calculations 
C  Overlap between Air Pump Tampering & Fuel Inlet Tampering or 
C  Overlap between EGR      Tampering & Fuel Inlet Tampering 
C 
   12 NCKIAY=IAY+1 
C 
      DIFF=AIREGR(IAY)-ECK(NCKIAY) 
      IF(DIFF.GE.0.0) GOTO 14 
C 
      AER(2,8,IPG,IAY)=ABS(DIFF) 
      AER(2,2,IPG,IAY)=1.00-ECK(NCKIAY) 
      GOTO 16 
C 
   14 AER(2,10,IPG,IAY)=DIFF 
      AER(2,2,IPG,IAY)=1.00-AIREGR(IAY) 
C 
C  Region 3 calculations 
C  Overlap between Air Pump Tampering & Other Misfueling or 
C  Overlap between EGR      Tampering & Other Misfueling 
C 
   16 DIFF=AIREGR(IAY)-TNK(IAY) 
      IF(DIFF.GE.0.0) GOTO 18 
C 
      AER(3,8,IPG,IAY)=ABS(DIFF) 
      AER(3,3,IPG,IAY)=1.00-TNK(IAY) 
      GOTO 20 
C 
   18 AER(3,11,IPG,IAY)=DIFF 
      AER(3,3,IPG,IAY)=1.00-AIREGR(IAY) 
C 
C  Region 4 calculations 
C  Overlap between Air Pump, Fuel Inlet & Catalyst Removal or 
C  Overlap between EGR,      Fuel Inlet & Catalyst Removal 
C 
   20 NCKIAY=1 
      IF(IAY.EQ.2) NCKIAY=3 
C 
      DIFFA=AIREGR(IAY)-CAT(IAY) 
      DIFFB=AIREGR(IAY)-ECK(NCKIAY) 
      DIFFC=CAT(IAY)-ECK(NCKIAY) 
C 
      IF(DIFFA.GE.0.0) GOTO 26 
C 
      IF(DIFFB.LT.0.0) GOTO 22 
C 
C  CAT > AIR > ECK  .OR.  CAT > EGR > ECK 
C 
      AER(4,4,IPG,IAY)=1.00-CAT(IAY) 
      AER(4,2,IPG,IAY)=ABS(DIFFA) 
      AER(4,10,IPG,IAY)=DIFFB 
      GOTO 32 
C 
   22 IF(DIFFC.GE.0.0) GOTO 24 
C 
C  ECK > CAT > AIR  .OR.  ECK > CAT > EGR 
C 
      AER(4,4,IPG,IAY)=1.00-ECK(NCKIAY) 
      AER(4,1,IPG,IAY)=ABS(DIFFC) 
      AER(4,8,IPG,IAY)=ABS(DIFFA) 
      GOTO 32 
C 
C  CAT > ECK > AIR  .OR.  CAT > ECK > EGR 
C 
   24 AER(4,4,IPG,IAY)=1.00-CAT(IAY) 
      AER(4,2,IPG,IAY)=DIFFC 
      AER(4,8,IPG,IAY)=ABS(DIFFB) 
      GOTO 32 
C 
   26 IF(DIFFB.GE.0.0) GOTO 28 
C 
C  ECK > AIR > CAT  .OR.  ECK > EGR > CAT 
C 
      AER(4,4,IPG,IAY)=1.00-ECK(NCKIAY) 
      AER(4,1,IPG,IAY)=ABS(DIFFB) 
      AER(4,9,IPG,IAY)=DIFFA 
      GOTO 32 
C 
   28 IF(DIFFC.GE.0.0) GOTO 30 
C 
C  AIR > ECK > CAT  .OR.  EGR > ECK > CAT 
C 
      AER(4,4,IPG,IAY)=1.00-AIREGR(IAY) 
      AER(4,6,IPG,IAY)=DIFFB 
      AER(4,9,IPG,IAY)=ABS(DIFFC) 
      GOTO 32 
C 
C  AIR > CAT > ECK  .OR.  EGR > CAT > ECK 
C 
   30 AER(4,4,IPG,IAY)=1.00-AIREGR(IAY) 
      AER(4,6,IPG,IAY)=DIFFA 
      AER(4,10,IPG,IAY)=DIFFC 
C 
C  Region 5 calculations 
C  Overlap between Air Pump, Other Misfueling & Catalyst Removal or 
C  Overlap between EGR,      Other Misfueling & Catalyst Removal 
C 
   32 DIFFA=AIREGR(IAY)-CAT(IAY) 
      DIFFB=AIREGR(IAY)-TNK(IAY) 
      DIFFC=CAT(IAY)-TNK(IAY) 
C 
      IF(DIFFA.GE.0.0) GOTO 38 
C 
      IF(DIFFB.LT.0.0) GOTO 34 
C 
C  CAT > AIR > TNK  .OR.  CAT > EGR > TNK 
C 
      AER(5,5,IPG,IAY)=1.00-CAT(IAY) 
      AER(5,3,IPG,IAY)=ABS(DIFFA) 
      AER(5,11,IPG,IAY)=DIFFB 
      GOTO 44 
C 
   34 IF(DIFFC.GE.0.0) GOTO 36 
C 
C  TNK > CAT > AIR  .OR.  TNK > CAT > EGR 
C 
      AER(5,5,IPG,IAY)=1.00-ECK(IAY) 
      AER(5,1,IPG,IAY)=ABS(DIFFC) 
      AER(5,8,IPG,IAY)=ABS(DIFFA) 
      GOTO 44 
C 
C  CAT > ECK > AIR  .OR.  CAT > ECK > EGR 
C 
   36 AER(5,5,IPG,IAY)=1.00-CAT(IAY) 
      AER(5,3,IPG,IAY)=DIFFC 
      AER(5,8,IPG,IAY)=ABS(DIFFB) 
      GOTO 44 
C 
   38 IF(DIFFB.GE.0.0) GOTO 40 
C 
C  ECK > AIR > CAT  .OR.  ECK > EGR > CAT 
C 
      AER(5,5,IPG,IAY)=1.00-TNK(IAY) 
      AER(5,1,IPG,IAY)=ABS(DIFFB) 
      AER(5,9,IPG,IAY)=DIFFA 
      GOTO 44 
C 
   40 IF(DIFFC.GE.0.0) GOTO 42 
C 
C  AIR > TNK > CAT  .OR.  EGR > TNK > CAT 
C 
      AER(5,5,IPG,IAY)=1.00-AIREGR(IAY) 
      AER(5,7,IPG,IAY)=DIFFB 
      AER(5,9,IPG,IAY)=ABS(DIFFC) 
      GOTO 44 
C 
C  AIR > CAT > TNK  .OR.  EGR > CAT > TNK 
C 
   42 AER(5,5,IPG,IAY)=1.00-AIREGR(IAY) 
      AER(5,7,IPG,IAY)=DIFFA 
      AER(5,11,IPG,IAY)=DIFFC 
C 
C  Region 6 calculations 
C  Overlap between Catalyst Removal and Fuel Inlet Tampering 
C 
   44 NCKIAY=1 
      IF(IAY.EQ.2) NCKIAY=3 
C 
      DIFF=CAT(IAY)-ECK(NCKIAY) 
      IF(DIFF.GE.0.0) GOTO 46 
C 
      AER(6,9,IPG,IAY)=ABS(DIFF) 
      AER(6,6,IPG,IAY)=1.00-ECK(NCKIAY) 
      GOTO 48 
C 
   46 AER(6,10,IPG,IAY)=DIFF 
      AER(6,6,IPG,IAY)=1.00-CAT(IAY) 
C 
C  Region 7 Calculations 
C  Overlap between Catalyst Removal and Other Misfueling 
C 
   48 DIFF=CAT(IAY)-TNK(IAY) 
      IF(DIFF.GE.0.0) GOTO 50 
C 
      AER(7,9,IPG,IAY)=ABS(DIFF) 
      AER(7,7,IPG,IAY)=1.00-TNK(IAY) 
      GOTO 52 
C 
   50 AER(7,11,IPG,IAY)=DIFF 
      AER(7,7,IPG,IAY)=1.00-CAT(IAY) 
C 
C  Regions 8,9,10,&11 Calculations 
C 
C  8 : Air Pump or EGR Tampering Only 
C  9 : Catalyst Removal Only 
C  10 : Fuel Inlet Tampering Only 
C  11 : Other Misfueling Only 
C 
   52 NCKIAY=IAY+1 
      AER(8,8,IPG,IAY)=1.00-AIREGR(IAY) 
      AER(9,9,IPG,IAY)=1.00-CAT(IAY) 
      AER(10,10,IPG,IAY)=1.00-ECK(NCKIAY) 
      AER(11,11,IPG,IAY)=1.00-TNK(IAY) 
C 
   54 CONTINUE 
C 
C  Evap. Canister, PCV, & Gas Cap calculations 
C 
      EVPAER(1,IAY)=1.00-CAN(IAY) 
      EVPAER(2,IAY)=1.00-PCV(IAY) 
      EVPAER(3,IAY)=1.00-CAP(IAY) 
C 
   56 CONTINUE 
C 
      RETURN 
C 
      END 
