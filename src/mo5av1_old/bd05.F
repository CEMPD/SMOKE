C Version @(#)$Id$ $Source$ $Date$ 

      BLOCK DATA BD05
C
C  BLOCK DATA Subprogram 05: /BASEQ8/,/BASEQ9/,/TIER1S/
C
C  BLOCK DATA 05 initializes the zero mile levels and deterioration rates
C  (intercepts & slopes) of the base LEV - Low Emitting Vehicles
C  exhaust emission rate equations.
C
C  Common block array subscripts:
C
C  LEVIMP(8,12,5)    -  LEVIMP ( LEV-TYPE, LMY, IV-CLASS )
C
C  LEVSTD(3,3,5,3)   -  LEVSTD ( LEN,TYPE,RATE,IP )
C
C  HDVZML(25,2)      -  HDVZML ( MY,   IREJN )
C  HDVBDR(25,2)      -  HDVBDR ( MY,   IREJN )
C  IHDVMY(25,2)      -  IHDVMY ( MY,   IREJN )
C  IPMOD(3)          -  IPMOD  ( IP )
C  TR1SL1(6,3,3)     -  TR1SL1 ( MY,   IP,   IV )
C
C  Common block dictionary:
C
C   Name   Type              Description
C  ------  ----  --------------------------------------------
C
C
C  /BASEQ8/:
C
C  The three following HDVXXX arrays contain the Old information to
C  correctly compute the Idle CO for HDV.
C
C  IHDVMY   R    HDV CO Idle Model Year
C  HDVZML   R    HDV CO Idle Zero Mile ER equation
C  HDVBDR   R    HDV CO Idle Base Deterioration rates
C
C  Local array subscripts:
C
C
C  /BASEQ9/:
C
C  LEVIMP   R    California LEV's % of Light Duty Fleet from
C                1992 to 2005
C  LEVSTD   R    Zero Mile, Deterioration Rates 1 and 2
C
C
C
C  Notes:
C
C  BD05 was modified to include kinked EFs for LDGT in MOBILE5.
C  BD05 was modified for post release MOBILE5 for LEV phase in
C  correction.
C
      REAL LEVSTD, LEVIMP, LDTAB
C
      COMMON /BASEQ8/ IHDVMY(23,2),HDVZML(23,2),HDVBDR(23,2)
      COMMON /BASEQ9/ LEVSTD(3,6,5,3),LEVIMP(8,12,5),LDTAB(2,2)
      COMMON /BASE10/ ZEVTF(25,3),ZEVTOT(3)
      COMMON /BASE12/ CAPIV1, CAPIV2
      COMMON /TIER1S/ TR1SL1(6,3,3)
C
      REAL LEVHC,  LEVCO,  LEVNO
      REAL LDGVIP, GT1AIP, GT1BIP, GT2AIP, GT2BIP
C
      DIMENSION LEVHC(3,6,5), LEVCO(3,6,5), LEVNO(3,6,5)
C
      DIMENSION LDGVIP(8,12), GT1AIP(8,12), GT1BIP(8,12)
      DIMENSION GT2AIP(8,12), GT2BIP(8,12)
C
      EQUIVALENCE  (LEVIMP(1,1,1),LDGVIP(1,1))
      EQUIVALENCE  (LEVIMP(1,1,2),GT1AIP(1,1))
      EQUIVALENCE  (LEVIMP(1,1,3),GT1BIP(1,1))
      EQUIVALENCE  (LEVIMP(1,1,4),GT2AIP(1,1))
      EQUIVALENCE  (LEVIMP(1,1,5),GT2BIP(1,1))
C
      EQUIVALENCE  (LEVSTD(1,1,1,1),LEVHC(1,1,1))
      EQUIVALENCE  (LEVSTD(1,1,1,2),LEVCO(1,1,1))
      EQUIVALENCE  (LEVSTD(1,1,1,3),LEVNO(1,1,1))
C
C  MOBILE5 Processes the California Emission Standards for
C  the LEV, Low Emission Vehicles.
C
C  Hardcoded FTP ER equation zero mile levels and deterioration (per 10000
C  miles) rates are blocked by vehicle class.  Within each block, the
C  continuation code is the pollutant type: H = HC, C = CO and N = NOx.
C  There are up to 2 lines and 20 cases per pollutant per vehicle class per
C  region per equation parameter.
C
C  /TIER1S/
C  Tier 1 LDGV's Slopes for I/M Special Credit
C
        DATA TR1SL1 /
C
C  LDGV
C   HC 1994  1995  1996  1997  1998  1999
     * .036, .034, .031, .029, .028, .028,
C   CO
     * .251, .251, .251, .251, .251, .251,
C   NOx
     * .082, .057, .044, .044, .044, .044,
C
C  LDGT1
     * .060, .044, .034, .032, .030, .030,
C   CO
     * .931, .622, .394, .319, .282, .282,
C   NOx
     * .073, .064, .059, .059, .059, .059,
C
C  LDGT2
C   HC
     * .076, .076, .074, .072, .072, .072,
C   CO
     *1.448,1.448,1.448,1.448,1.448,1.448,
C   NOx
     * .083, .083, .083, .083, .083, .083 /
C
C  /BASEQ8/
C
C  HDV Idle CO fix - in BEFIDL HDV(IV=4) gets the same BEF as LDGV2(IV=3)
C  IHDVMY HDV Model Year Groups ( same for Low & High Altitude )
C
        DATA IHDVMY/
C
C  HDV  Model Year Groups  Low altitude
C
     C 1969, 1973, 1978, 1980, 1983, 1984, 1985, 1986, 1987, 1988,
     C 1989, 1990, 1991, 2020, 9*0,
C
C  HDV  Model Year Groups  - High altitude
C
     C 1969, 1973, 1978, 1980, 1981, 1983, 1984, 1985, 1986, 1987,
     C 1988, 1989, 1990, 1991, 2020, 8*0/
C
        DATA HDVZML/
C
C  HDV Zero Mile Low altitude
     C 93.98, 60.08, 60.08, 12.28, 12.58, 9.43, 7.03, 5.76, 4.42,
     C 5.29, 5.26, 5.22, 5.21, 5.07, 9*0.,
C
C  HDV Zero Mile High altitude
C
     C141.35, 107.72, 107.72, 44.25, 44.25, 30.16, 23.35, 9.84, 8.06,
     C 6.19, 7.41, 7.36, 7.31, 7.29, 7.10, 8*0./
C
C  HDV Low & High altitude Deterioration rates
C
        DATA HDVBDR/
C  HDV Basic Deterioration Rates Low altitude
     C 2.25, 2.55, 2.44, 2.43, 1.46, 1.46, 0.73, 0.73, 0.73, 0.73,
     C 0.73, 0.73, 0.73, 0.73, 9*0.,
C
C  HDV Basic Deterioration Rates High Altitude
     C 2.25, 2.55, 2.44, 2.43, 2.43, 1.46, 1.46, 0.73, 0.73, 0.73,
     C 0.73, 0.73, 0.73, 0.73, 0.73, 8*0./
C
C     Arrays for LDGV, LDGT1, and LDGT2 Vehicles as
C     low emission vehicles ( LEV's ) equivalenced into LEVIMP
C
C     LDGV   -> LDGVIP(8,12)
C     LDGT1A -> GT1AIP(8,12)
C     LDGT1B -> GT1BIP(8,12)
C     LDGT2A -> GT2AIP(8,12)
C     LDGT2B -> GT2BIP(8,12)
C
C    LDGVIP(8,12)  LDGV - LEV - IMPLEMENTATION Phase In
C
      DATA LDGVIP /
C MYR  TIER1/ TLEV / TLEV  / LEV  / LEV  / ULEV / ULEV / ZEV
C1994       / inter/       / inter/      / inter/      /
     * 0.900, 0.100, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
C1995
     * 0.850, 0.150, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
C1996
     * 0.800, 0.000, 0.200, 0.000, 0.000, 0.000, 0.000, 0.000,
C1997
     * 0.730, 0.000, 0.000, 0.250, 0.000, 0.020, 0.000, 0.000,
C1998
     * 0.480, 0.000, 0.000, 0.480, 0.000, 0.020, 0.000, 0.020,
C1999
     * 0.230, 0.000, 0.000, 0.000, 0.730, 0.000, 0.020, 0.020,
C2000
     * 0.000, 0.000, 0.000, 0.000, 0.960, 0.000, 0.020, 0.020,
C2001
     * 0.000, 0.000, 0.000, 0.000, 0.900, 0.000, 0.050, 0.050,
C2002
     * 0.000, 0.000, 0.000, 0.000, 0.850, 0.000, 0.100, 0.050,
C2003
     * 0.000, 0.000, 0.000, 0.000, 0.750, 0.000, 0.150, 0.100,
C2004
     * 0.000, 0.000, 0.000, 0.000, 0.750, 0.000, 0.150, 0.100,
C2005
     * 0.000, 0.000, 0.000, 0.000, 0.750, 0.000, 0.150, 0.100/
CC
CC
CC
CC    GT1AIP(8,12)  LDGT1A - LEV - IMPLEMENTATION Phase In
CC
        DATA GT1AIP /
C MYR  TIER1/ TLEV / TLEV  / LEV  / LEV  / ULEV / ULEV / ZEV
C1994       / inter/       / inter/      / inter/      /
     * 0.900, 0.100, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
C1995
     * 0.850, 0.150, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
C1996
     * 0.800, 0.000, 0.200, 0.000, 0.000, 0.000, 0.000, 0.000,
C1997
     * 0.730, 0.000, 0.000, 0.250, 0.000, 0.020, 0.000, 0.000,
C1998
     * 0.480, 0.000, 0.000, 0.480, 0.000, 0.020, 0.000, 0.020,
C1999
     * 0.230, 0.000, 0.000, 0.000, 0.730, 0.000, 0.020, 0.020,
C2000
     * 0.000, 0.000, 0.000, 0.000, 0.960, 0.000, 0.020, 0.020,
C2001
     * 0.000, 0.000, 0.000, 0.000, 0.900, 0.000, 0.050, 0.050,
C2002
     * 0.000, 0.000, 0.000, 0.000, 0.850, 0.000, 0.100, 0.050,
C2003
     * 0.000, 0.000, 0.000, 0.000, 0.750, 0.000, 0.150, 0.100,
C2004
     * 0.000, 0.000, 0.000, 0.000, 0.750, 0.000, 0.150, 0.100,
C2005
     * 0.000, 0.000, 0.000, 0.000, 0.750, 0.000, 0.150, 0.100/
CC
CC
CC
CC    GT1BIP(8,12)  LDGT1B - LEV - IMPLEMENTATION Phase In
CC
        DATA GT1BIP /
C MYR  TIER1/ TLEV / TLEV  / LEV  / LEV  / ULEV / ULEV / ZEV
C1994       / inter/       / inter/      / inter/      /
     * 0.900, 0.100, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
C1995
     * 0.850, 0.150, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
C1996
     * 0.800, 0.000, 0.200, 0.000, 0.000, 0.000, 0.000, 0.000,
C1997
     * 0.730, 0.000, 0.000, 0.250, 0.000, 0.020, 0.000, 0.000,
C1998
     * 0.480, 0.000, 0.000, 0.500, 0.000, 0.000, 0.020, 0.000,
C1999
     * 0.230, 0.000, 0.000, 0.000, 0.750, 0.000, 0.020, 0.000,
C2000
     * 0.000, 0.000, 0.000, 0.000, 0.980, 0.000, 0.020, 0.000,
C2001
     * 0.000, 0.000, 0.000, 0.000, 0.950, 0.000, 0.050, 0.000,
C2002
     * 0.000, 0.000, 0.000, 0.000, 0.900, 0.000, 0.100, 0.000,
C2003
     * 0.000, 0.000, 0.000, 0.000, 0.850, 0.000, 0.150, 0.000,
C2004
     * 0.000, 0.000, 0.000, 0.000, 0.850, 0.000, 0.150, 0.000,
C2005
     * 0.000, 0.000, 0.000, 0.000, 0.850, 0.000, 0.150, 0.000/
CC
CC
CC    GT2AIP(8,12)  LDGT2A - LEV - IMPLEMENTATION Phase In
CC
        DATA GT2AIP /
C MYR  TIER1/ TLEV / TLEV  / LEV  / LEV  / ULEV / ULEV / ZEV
C1994       / inter/       / inter/      / inter/      /
     * 1.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
C1995
     * 1.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
C1996
     * 1.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
C1997
     * 1.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
C1998
     * 0.730, 0.250, 0.000, 0.000, 0.000, 0.020, 0.000, 0.000,
C1999
     * 0.480, 0.500, 0.000, 0.000, 0.000, 0.020, 0.000, 0.000,
C2000
     * 0.230, 0.000, 0.000, 0.000, 0.750, 0.000, 0.020, 0.000,
C2001
     * 0.000, 0.000, 0.000, 0.000, 0.950, 0.000, 0.050, 0.000,
C2002
     * 0.000, 0.000, 0.000, 0.000, 0.900, 0.000, 0.100, 0.000,
C2003
     * 0.000, 0.000, 0.000, 0.000, 0.850, 0.000, 0.150, 0.000,
C2004
     * 0.000, 0.000, 0.000, 0.000, 0.850, 0.000, 0.150, 0.000,
C2005
     * 0.000, 0.000, 0.000, 0.000, 0.850, 0.000, 0.150, 0.000/
CC
CC
CC
CC    GT2BIP(8,12)  LDGT2B - LEV - IMPLEMENTATION Phase In
CC
        DATA GT2BIP /
C MYR  TIER1/ TLEV / TLEV  / LEV  / LEV  / ULEV / ULEV / ZEV
C1994       / inter/       / inter/      / inter/      /
     * 1.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
C1995
     * 1.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
C1996
     * 1.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
C1997
     * 1.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
C1998
     * 0.730, 0.000, 0.000, 0.250, 0.000, 0.000, 0.020, 0.000,
C1999
     * 0.480, 0.000, 0.000, 0.500, 0.000, 0.000, 0.020, 0.000,
C2000
     * 0.230, 0.000, 0.000, 0.000, 0.750, 0.000, 0.020, 0.000,
C2001
     * 0.000, 0.000, 0.000, 0.000, 0.950, 0.000, 0.050, 0.000,
C2002
     * 0.000, 0.000, 0.000, 0.000, 0.900, 0.000, 0.100, 0.000,
C2003
     * 0.000, 0.000, 0.000, 0.000, 0.850, 0.000, 0.150, 0.000,
C2004
     * 0.000, 0.000, 0.000, 0.000, 0.850, 0.000, 0.150, 0.000,
C2005
     * 0.000, 0.000, 0.000, 0.000, 0.850, 0.000, 0.150, 0.000/
CC
CC
C     initilize the TF for ZEV, Zero Emission Vehicles
C
      DATA ZEVTF / 75*0.0 /
      DATA ZEVTOT/  3*0.0 /
CC
CC    LEVSTD(3,6,5,3) - LEV Standards for Special I/M Test
CC     Intermediate TLEV      TLEV                 Intermediate LEV
CC     ZML   DET1    DET2     ZML   DET1    DET2   ZML   DET1    DET2
CC     LEV                    Intermediate ULEV    ULEV
CC     ZML   DET1    DET2     ZML   DET1    DET2   ZML   DET1    DET2
CC
CC    - HC -
      DATA LEVHC /
C     LDGV
     *0.134,0.0175,0.0175, 0.089,0.0116,0.0116, 0.075, 0.0098, 0.0098,
     *0.056,0.0073,0.0073, 0.043,0.0057,0.0055, 0.030, 0.0039, 0.0039,
C     LDGT1a
     *0.134,0.0175,0.0175, 0.089,0.0116,0.0116, 0.075, 0.0098, 0.0098,
     *0.056,0.0073,0.0073, 0.043,0.0057,0.0055, 0.030, 0.0039, 0.0039,
C     LDGT1b
     *0.169,0.0221,0.0221, 0.114, 0.0149, 0.0149, 0.096, 0.0125,0.0125,
     *0.075,0.0098,0.0098, 0.056, 0.0073, 0.0073, 0.037, 0.0049,0.0049,
C     LDGT2a
     *0.169,0.0221,0.0221, 0.114, 0.0149, 0.0149, 0.096, 0.0125,0.0125,
     *0.075,0.0098,0.0098, 0.056, 0.0073, 0.0073, 0.037, 0.0049,0.0049,
C     LDGT2b
     *0.169,0.0221,0.0221, 0.114, 0.0149, 0.0149, 0.096, 0.0125,0.0125,
     *0.075,0.0098,0.0098, 0.056, 0.0073, 0.0073, 0.037, 0.0049,0.0049 /
CC
CC    - CO -
      DATA LEVCO /
C      LDGV
     *2.482,0.2896,0.2896, 2.482,0.2896,0.2896, 0.620,0.4758,0.4758,
     *0.620,0.4758,0.4758, 0.620,0.4758,0.4758, 0.620,0.4758,0.4758,
C     LDGT1a
     *2.482,0.2896,0.2896, 2.482,0.2896,0.2896, 0.620,0.4758,0.4758,
     *0.620,0.4758,0.4758, 0.620,0.4758,0.4758, 0.620,0.4758,0.4758,
C     LDGT1b
     *3.211,0.3748,0.3748, 3.211,0.3748,0.3748, 0.803,0.6157,0.6157,
     *0.803,0.6157,0.6157, 0.803,0.6157,0.6157, 0.803,0.6157,0.6157,
C     LDGT2a
     *3.211,0.3748,0.3748, 3.211,0.3748,0.3748, 0.803,0.6157,0.6157,
     *0.803,0.6157,0.6157, 0.803,0.6157,0.6157, 0.803,0.6157,0.6157,
C     LDGT2b
     *3.211,0.3748,0.3748, 3.211,0.3748,0.3748, 0.803,0.6157,0.6157,
     *0.803,0.6157,0.6157, 0.803,0.6157,0.6157, 0.803,0.6157,0.6157 /
CC
CC    - NO -
      DATA LEVNO /
C      LDGV
     * .174,.0433,.0433, .174,.0433,.0433, .131,.0325,.0325,
     * .087,.0217,.0217, .131,.0325,.0325, .087,.0217,.0217,
C      LDGT1a
     * .174,.0433,.0433, .174,.0433,.0433, .131,.0325,.0325,
     * .087,.0217,.0217, .131,.0325,.0325, .087,.0217,.0217,
C      LDGT1b
     * .305,.0758,.0758, .305,.0758,.0758, .218,.0542,.0542,
     * .174,.0433,.0433, .218,.0542,.0542, .174,.0433,.0433,
C      LDGT2a
     * .305,.0758,.0758, .305,.0758,.0758, .218,.0542,.0542,
     * .174,.0433,.0433, .218,.0542,.0542, .174,.0433,.0433,
C      LDGT2b
     * .305,.0758,.0758, .305,.0758,.0758, .218,.0542,.0542,
     * .174,.0433,.0433, .218,.0542,.0542, .174,.0433,.0433 /
CC
CC    LDTAB  R  Weighting fraction for light and heavy duty
CC              LDGT1 and LDGT2
CC
CC                 (1,1) (2,1) (1,2) (2,2)
      DATA LDTAB / 0.57, 0.65, 0.43, 0.35 /
CC
CC    CAPIV1 Upper Limit of the VMTAGE LDGV IV=1
CC    CAPIV2 Upper Limit of the VMTAGE LTGT1/2 IV=2/3 in LEVEF & BEFUN
CC    in SUBROUTINES LEVEF & BEFUN
CC
      DATA CAPIV1 / 10.0 /
      DATA CAPIV2 / 12.0 /
CC
CC
      END
