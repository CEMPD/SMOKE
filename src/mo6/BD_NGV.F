      BLOCK DATA BD_NGV 
C 
C     BLOCK BD_NGV initializes the NGV variables 
C 
C  Changes: (Last change first) 
C 
C   24 apr 00 @EPA ELG   Increased size of NGVPHIN 
C   24 FEB 00 @EPA ELG   New data block 
C 
C 
C  Common block dictionary: 
C 
C   Name    Type              Description 
C  ------   ----  ---------------------------------------------- 
C  NGVDET    R    Deterioration rate (10,000K miles) of NGV vehicles. 
C  NGVFLG1   I    Toggle for Alternate NGV Phase-in Fractions 0-DEFAULT 1-Alternate. 
C  NGVFLG2   I    Toggle for Alternate NGV Emission Factors 0-DEFAULT 1-Alternate. 
C  NGVHIGH   R    High Emitter Emission Level for NGVs. 
C  NGVPHIN   R    NGV Phase-in Fractions BY model year and vehicles TYPE. 
C  NGV_START I    Earliest NGV Vehicle that can be modeled.  Default = 1994 MY 
C  NGVZML    R    Zero mile Emission Level of NGV vehicles. 
C 
C 
      IMPLICIT NONE 
C 
      INCLUDE 'NGV.I' 
C 
      DATA NGVFLG1, NGVFLG2 / 2*0 / 
      DATA NGV_START / 1994 / 
C 
      DATA NGVPHIN / 1596*0.00/ 
C 
      DATA NGVZML   / 168*-1.0 / 
      DATA NGVDET   / 168*-1.0 / 
      DATA NGVHIGH  / 168*-1.0 / 
C 
C     HC is NMHC. 
C                       HC     CO    NOX 
      DATA NGV_NZML / 0.016, 0.665, 0.133,       !LDV 
     *                0.016, 0.665, 0.133,       !LDT1 
     *                0.020, 0.916, 0.265,       !LDT2 
     *                0.020, 0.916, 0.265,       !LDT3 
     *                0.046, 1.040, 0.398 /      !LDT4 
C 
C                        HC     CO      NOX 
      DATA NGV_NDET / 0.0018, 0.2340, 0.0070,    !LDV 
     *                0.0018, 0.2340, 0.0070,    !LDT1 
     *                0.0023, 0.2800, 0.0150,    !LDT2 
     *                0.0023, 0.2800, 0.0150,    !LDT3 
     *                0.0053, 0.2800, 0.0220 /   !LDT4 
C 
C                       HC      CO    NOX 
      DATA NGV_HIGH / 0.068,  5.800, 0.740,      !LDV 
     *                0.068,  5.800, 0.740,      !LDT1 
     *                0.070,  5.900, 0.990,      !LDT2 
     *                0.070,  5.900, 0.990,      !LDT3 
     *                0.080,  5.900, 1.230 /     !LDT4 
C 
C 
C                       HC     CO    NOX 
      DATA NGV_HDZML/ 0.049, 1.116, 0.427,       !LIGHT 
     *                0.620, 0.340, 1.840,       !MEDIUM 
     *                0.440, 0.428, 1.840 /      !HEAVY 
C 
C                        HC     CO      NOX 
      DATA NGV_HDDET/ 0.0050, 0.3180, 0.0450,    !LIGHT 
     *                0.0020, 0.0040, 0.0010,    !MEDIUM 
     *                0.0020, 0.0020, 0.0030 /   !HEAVY 
C 
C 
C 
C     Ratios Taken from !From M6.EXH.007 (HC and NOx) and M6.EXH.009 (CO) 
C 
C                        HC     CO      NOX 
      DATA NGV_P    / 0.2300, 0.3380, 0.9000,    !ISR = RUNNING 
     *               11.0000, 4.1490, 1.2700 /   !ISR = START 
C 
C 
C 
C 
      END 
