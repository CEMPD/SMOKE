      BLOCK DATA BD11 
C 
C  BLOCK DATA Subprogram 11: /EVAHS1/EVAPAR/,/EVAPHD/,/EVASTD/, 
C                            /EVAPGR/,/ETPPAR/,/EVAPACT/ 
C 
C  BLOCK DATA 11 initializes the parameters and rates affecting both hot soak 
C  and diurnal or hot soak alone (Hot soak HC + Diurnal HC = Evaporative HC). 
C  BD11 also intializes the Evaportive Test Procedure (ETP) parameters. 
C 
C  Changes: (Last change first) 
C 
C  27 Feb 01 @EPA-BAg bug314 removed HSIV5 and RLV5 variables. 
C  31 Oct 00 @EPA-elg Fixed Bug 282 - Hot soak coefficents were changed. 
C  06 Oct 00 @EPA-BG  BUG155 removed NEWETP flag previously used for user input 
C            of ETP phase in 
C  15 Sep 00 AIR Task 08: Updated the ETP phase in rates. 
C  20 Jul 00 @EPA-elg E10 Tier2 MDPV; added HSIV5(hot soak) and RLV5(resting loss) 
C            initialization. 
C  11 Jul 00 AIR Task 06: Added initialization of Hot Soak Activity Fractions 
C            and Fixed P/P Rates. 
C  27 Jul 00 @EPA-djb Bub 208, removed MAXSTD. 
C   2 Jun 00 @EPA-ddj Bug 86, Added back HDGVs values in array KSTD81. 
C  22 Jan 00 @DynTel-bag 1-005  Removed diurnal parameters 
C  13 Dec 99 @DynTel-bag 1-006  Removed partial diurnal arrays 
C  15 Jan 99 @DynTel-bag 2-664  Added hot soak fraction array, HSFRAC(IMIN), 
C            rather than calculating them in HOTSOK 
C  06 Nov 98 @DynTel-bag 2-000  Added temporary storage variables and EQUIVLENCE 
C            statements for DIUFRC0 - DIUFRC7 and DIUACT0 - DIUACT7 for the 
C            diurnal parameters, DIURNL_FRCTN and NVHCLES_SOAKING, respectively, 
C            to limit continuation lines used in DATA statements.(LF90 compiling) 
C  30 Oct 98 @DynTel-bag 2-677  Added parameters for hourly diurnal calculations. 
C   8 Oct 98 @DynTel-MLA 2-670  Revisd HSCOEF coefficients  
C  24 Sep 98 @DynTel-MLA 2-670  Added HSCOEF, coefficients for the hourly 
C            hot soak distribution equations. 
C 
C  11 Aug 98 @DynTel-MLA 2-000  Changed HDG and HDD names to conform to 
C            new naming conventions. 
C 
C  Common block array subscripts: 
C 
C  EFHS(3,2,5,15)   -  EFHS   ( KEYEQ  , IRVP , ISTD, IVGAS ) 
C  ETP(7,3)         -  ETP    ( IETP, KEYEQ ) 
C  ETPPIV(4)        -  ETPPIV ( IETPY ) 
C  EVP(4)           -  EVP    ( IEVP ) 
C  FINJ(2)          -  FINJ   ( 2 ) 
C  GREVP(6,29)      -  GREVP  ( IEVP , IV ) 
C  HDWGT(2)         -  HDWGT  ( by weight class; IHDSAL ) 
C  HIHS(3,2)        -  HIHS   ( KEYEQ, IRVP ) 
C  HSCOEF(4,14,2)   -  HSCOEF ( ICOEF, IH, IWEEK ) 
C  HSEQ(3,2,9)      -  HSEQ   ( IFDS  , IVV  , IEVEQ ) 
C  HSFRACT(60)      -  HSFRACT( IMIN ) 
C  IETPYR(4)        -  IETPYR ( IETPY ) 
C  IVVMAP(15)       -  IVVMAP ( IVGAS ) 
C  KSTD81(5,3)      -  KSTD81 ( IVTAM , IFDS ) 
C  MYGSTD(5,15)     -  MYGSTD ( IG    , IVGAS ) 
C  HSACTUSER(60,14) -  HSACTUSER( IMIN  , IHR ) 
C  DSACTUSER(60,14) -  DSACTUSER( IH  , IHR_SOAK ) 
C 
C  Common block dictionary: 
C 
C   Name             Type              Description 
C  ------            ----  --------------------------------------------- 
C  /EVAHIA/: 
C 
C  HIHS               R    High alt, 1972-76, LDGV & LDGT12, HS rates on 2 RVP fuels 
C  HIADJ              R    High alt adjustment to low alt based evap rate 
C 
C  /EVAHS1/: 
C 
C  HSEQ               R    1981+ LDGT FI AND TBI HS equations parameters 
C  HSEQ_81            R    1981-85 MY LDGV & LDGT (carb only) HS equation parameters 
C  HSEQ_81            R    1986+ MY LDGV & LDGT (carb only) HS equation parameters 
C  HSFRACT            R    hot soak fractions for hourly hot soak calculations 
C  EFHS               R    EF HS rates for IVGAS vehicles on Pass/Fail and 2 RVP fuels 
C  HS_HRLY_TRP_FRCN   R    hot soak hourly trip fractions 
C  HSIV5              R    hot soak for IV=5 used to calculate HS for IV=6 Tier2 vehicles. 
C 
C  /EVAPAR/: 
C 
C  IFDS               I    See Subscript Dictionary 
C  ISTDC              I    Carbureted engine emission standard case 
C  ISTDT              I    TBI engine emission standard case 
C  ISTDP              I    PFI engine emission standard case 
C  FINJ               R    Fuel-injected fraction of vehicle class IVGAS built in MY, 
C                          index is 1 = TBI, 2 = PBI 
C 
C  /EVAPHD/: 
C 
C  HDWGT              R    HDGV weighting factor on value returned from EV81 
C  IVVMAP             I    IVGAS to IVV mapping vector 
C 
C  /EVASTD/: 
C 
C  MYGSTD             I    Model year groups for emission standards by vehicle type 
C  KSTD81             I    Emission standard key for IFDS sensitive 1981+ LDGV/T 
C 
C  /EVAPGR/: 
C 
C  EVP                R    Gram evaporative component for a model year 
C  GREVP              R    Gram evaporative component for a vehicle type 
C  RSLOSS             R    incremental resting loss HC 
C 
C  /ETPPAR/: 
C 
C  ETP                R    Evaporative Test Procedure rates for Hot Soak, Full 
C                          Diurnal, Multiple Diurnal 2-3, Multiple Diurnal 4+, 
C                          Partial Diurnal, Running Loss and Resting Loss for 
C                          Problem Free, Purge Failure and Pressure Failure 
C  IETPYR             I    The phase-in years of the ETP program. 
C  ETPPIV             I    The phase-in values (percentages) of the ETP program. 
C  NETPYR             I    The number of ETP phase-in years 
C 
C  /EVAPACT/: 
C 
C  SAWFLAG            I    Flag to disable/enable evap sawtooth weighting 
C  HSACTUSER          R    User-supplied hot soak activity rates 
C  HSACTFLAG          I    Flag to trigger user-supplied hot soak activity rates 
C  FPPRUSER           R    User-supplied fixed P/P rates (pass/pass, pass pressure/fail purge, 
C                          fail pressure/pass purge, fail/fail) 
C 
C  Local array subscripts: 
C 
C  LDHS(3,2,5)     -  LDHS     ( KEYEQ, IRVP , ISTD  ) 
C  T1AHS(3,2,5)    -  T1AHS    ( KEYEQ, IRVP , ISTD  ) 
C  T1BHS(3,2,5)    -  T1BHS    ( KEYEQ, IRVP , ISTD  ) 
C  T2AHS(3,2,5)    -  T2AHS    ( KEYEQ, IRVP , ISTD  ) 
C  T2BHS(3,2,5)    -  T2BHS    ( KEYEQ, IRVP , ISTD  ) 
C  HD2BHS(3,2,5)   -  HD2BHS   ( KEYEQ, IRVP , ISTD  ) 
C  HD3HS(3,2,5)    -  HD3HS    ( KEYEQ, IRVP , ISTD  ) 
C  HD4HS(3,2,5)    -  HD4HS    ( KEYEQ, IRVP , ISTD  ) 
C  HD5HS(3,2,5)    -  HD5HS    ( KEYEQ, IRVP , ISTD  ) 
C  HD6HS(3,2,5)    -  HD6HS    ( KEYEQ, IRVP , ISTD  ) 
C  HD7HS(3,2,5)    -  HD7HS    ( KEYEQ, IRVP , ISTD  ) 
C  HD8AHS(3,2,5)   -  HD8AHS   ( KEYEQ, IRVP , ISTD  ) 
C  HD8BHS(3,2,5)   -  HD8BHS   ( KEYEQ, IRVP , ISTD  ) 
C  CLHS(3,2,5)     -  CLHS     ( KEYEQ, IRVP , ISTD  ) 
C  GBUSHS(3,2,5)   -  GBUSHS   ( KEYEQ, IRVP , ISTD  ) 
C  HSEQ1(3,2,5)    -  HSEQ1    ( IFDS ,  IVV , IEVEQ ) 
C  HSEQ2(3,2,4)    -  HSEQ2    ( IFDS ,  IVV , IEVEQ ) 
C 
C  Local array dictionary: 
C 
C   Name   Type                Description 
C  ------  ----  ----------------------------------------------- 
C  LDHS      R    EF HS rates for pre-81 LDGV 
C  T1AHS     R    EF HS rates for pre-81 LDGT1 
C  T1BHS     R    EF HS rates for pre-81 LDGT2 
C  T2AHS     R    EF HS rates for pre-81 LDGT3 
C  T2BHS     R    EF HS rates for pre-81 LDGT4 
C  HD2BHS    R    EF HS rates for pre-85 HDGV2B 
C  HD3HS     R    EF HS rates for pre-85 HDGV3 
C  HD4HS     R    EF HS rates for pre-85 HDGV4 
C  HD5HS     R    EF HS rates for pre-85 HDGV5 
C  HD6HS     R    EF HS rates for pre-85 HDGV6 
C  HD7HS     R    EF HS rates for pre-85 HDGV7 
C  HD8AHS    R    EF HS rates for pre-85 HDGV8A 
C  HD8BHS    R    EF HS rates for pre-85 HDGV8B 
C  CLHS      R    EF HS rates for MC 
C  GBUSHS    R    EF HS rates for GASBUS 
C  HSEQ1     R    1981+ LDGV and LDGT1-4 HS equations parameters 
C                        for IEVEQ = 1-5 
C  HSEQ2     R    1981+ LDGV and LDGT1-4 HS temperature equations parameters 
C                        for IEVEQ = 6-10 
C 
C  Notes: 
C 
C  These data structures are in support of the Release MOBILE4.1 evaporative 
C  emissions calculation methodology. 
C 
C  MOBILE4.1 changes the units of the evaporative component variables 
C  under HCFLAG=3 from G/Day to G/Mi (HS) and from Grams to G/Mi (Du, 
C  Mult & all Partials). Crankcase was already in G/Mi. 
C  Block Data 11 was modified for MOBILE5 to make Hot Soak TCF for 
C  LDGT=LDGV 
C 
C 18 Nov 96 DynTel@MLA 2-622  Expanded the arrays to accomodate the new  
C           Heavy Duty Gasoline, Heavy Duty Diesel, and LDGT vehicle types. 
C 
C  4 Apr 97 DynTel@MLA 2-622  Expanded arrays EFHS, GREVP, GRPD, KSTD81, 
C           and MYGSTD to accomodate the new BUS vehicle types. 
C 
C 11 Jun 97 DynTel@MLA 2-622  It appears from reviewing the code that 
C           KSTD81 is used only by LDGV, LDGT1, and LGDT2. The first 
C           dimension of this array therefore has been reduced to 5. 
C 
C 10 Mar 98 DynTel@RJD 2-622 Added integer labels for emission factors to 
C           facilitate the efficient calculation of emission factors for 
C           composite vehicles in EFCVT. 
C 
C 23 Jul 98 @DynTel-MLA 2-000  Changed LDGT names to conform to 
C           certification naming conventions. 
C 
      IMPLICIT NONE 
C 
      INCLUDE 'ETPPAR.I' 
      INCLUDE 'EVAHIA.I' 
      INCLUDE 'EVAHS1.I' 
      INCLUDE 'EVAPACT.I' 
      INCLUDE 'EVAPAR.I' 
      INCLUDE 'EVAPGR.I' 
      INCLUDE 'EVAPHD.I' 
      INCLUDE 'EVASTD.I' 
      INCLUDE 'PARAME.I' 
C 
      REAL LDHS 
      REAL T1AHS,T1BHS,T2AHS, T2BHS,HD2BHS,HD3HS,HD4HS, 
     *          HD5HS,HD6HS,HD7HS,HD8AHS,HD8BHS,CLHS,GBUSHS 
C 
      DIMENSION LDHS(3,2,5),T1AHS(3,2,5),T1BHS(3,2,5),T2AHS(3,2,5), 
     *          T2BHS(3,2,5),HD2BHS(3,2,5),HD3HS(3,2,5),HD4HS(3,2,5), 
     *          HD5HS(3,2,5),HD6HS(3,2,5),HD7HS(3,2,5),HD8AHS(3,2,5), 
     *          HD8BHS(3,2,5),CLHS(3,2,5),GBUSHS(3,2,5) 
C 
      REAL  HSEQ1,HSEQ2,HSEQ3,HSEQ4,HSEQ5,HSEQ6 
C 
      DIMENSION HSEQ1(3,2,5),HSEQ2(3,2,4) 
      DIMENSION HSEQ3(3,2,5),HSEQ4(3,2,4) 
      DIMENSION HSEQ5(3,2,5),HSEQ6(3,2,4) 
 
C 
      EQUIVALENCE 
     *     (EFHS(1,1,1, 1),LDHS(1,1,1)), 
     *     (EFHS(1,1,1, 2),T1AHS(1,1,1)), 
     *     (EFHS(1,1,1, 3),T1BHS(1,1,1)), 
     *     (EFHS(1,1,1, 4),T2AHS(1,1,1)), 
     *     (EFHS(1,1,1, 5),T2BHS(1,1,1)), 
     *     (EFHS(1,1,1, 6),HD2BHS(1,1,1)), 
     *     (EFHS(1,1,1, 7),HD3HS(1,1,1)), 
     *     (EFHS(1,1,1, 8),HD4HS(1,1,1)), 
     *     (EFHS(1,1,1, 9),HD5HS(1,1,1)), 
     *     (EFHS(1,1,1,10),HD6HS(1,1,1)), 
     *     (EFHS(1,1,1,11),HD7HS(1,1,1)), 
     *     (EFHS(1,1,1,12),HD8AHS(1,1,1)), 
     *     (EFHS(1,1,1,13),HD8BHS(1,1,1)), 
     *     (EFHS(1,1,1,14),CLHS(1,1,1)), 
     *     (EFHS(1,1,1,15),GBUSHS(1,1,1)), 
     *     (HSEQ(1,1,1),HSEQ1(1,1,1)), 
     *     (HSEQ(1,1,6),HSEQ2(1,1,1)) 
C 
      EQUIVALENCE 
     *     (HSEQ_81(1,1,1),HSEQ3(1,1,1)), 
     *     (HSEQ_81(1,1,6),HSEQ4(1,1,1)), 
     *     (HSEQ_86(1,1,1),HSEQ5(1,1,1)), 
     *     (HSEQ_86(1,1,6),HSEQ6(1,1,1)) 
C 
C 
C  /EVAHIA/: 
C 
C                      At 9.0 RVP                 At 11.5 RVP 
C                       Failed   Failed            Failed   Failed 
C                Pass    Purge Pressure     Pass    Purge Pressure 
      DATA HIHS/ 
     1          10.03,   14.52,   19.07,   18.87,   27.30,   36.74/ 
C 
      DATA HIADJ/1.30/ 
C 
C  /EVAHS1/: 
C 
C     Hot soak emissions hourly trip fractions 
C      
      DATA HS_HRLY_TRP_FRCN/ 
C   weekdays : 
C  IH =  1       2        3       4       5       6       7      8 
     D 0.0233, 0.0605, 0.0630, 0.0462, 0.0508, 0.0632, 0.0780, 0.0732, 
C 
C******************************************************************* 
C  IH =  9       10      11      12       13     14 
C     D 0.0787, 0.0863, 0.0871, 0.0799, 0.0588, 0.151, 
C******************************************************************* 
C 
C  IH =  9       10      11      12       13     14      15      16 
     D 0.0787, 0.0863, 0.0871, 0.0799, 0.0588, 0.0137, 0.0137, 0.0137, 
C  IH =  17      18      19      20       21     22      23      24 
     D 0.0137, 0.0137, 0.0137, 0.0137, 0.0137, 0.0137,0.0137, 0.0137, 
C******************************************************************* 
C   weekends : 
C  IH =  1       2        3       4       5       6       7      8 
     E 0.0099, 0.0226, 0.0338, 0.0641, 0.0698, 0.088 , 0.0923, 0.0740, 
C 
C******************************************************************* 
C  IH =  9       10      11      12       13     14 
C     E 0.081 , 0.0662, 0.0803, 0.0691, 0.0627, 0.186 / 
C******************************************************************* 
C 
C  IH =  9       10      11      12       13     14      15      16 
     E 0.081 , 0.0662, 0.0803, 0.0691, 0.0627, 0.0169, 0.0169, 0.0169, 
C  IH =  17      18      19      20       21     22      23      24 
     E 0.0169, 0.0169, 0.0169, 0.0169, 0.0169, 0.0169, 0.0169, 0.0169 / 
C 
C     Hot soak fractions for times of 1 to 60 minutes. 
C     Set the hot soak fraction at 60 minutes to 1.0 ('full hot soak'). 
C         HSFRACT(IMIN) = 0.0258*IMIN - 0.000156*IMIN*IMIN 
C 
      DATA HSFRACT/ 2.56E-02,5.10E-02,7.60E-02,0.101, 0.125, 0.149, 
     *              0.173,   0.196,   0.220,   0.242, 0.265, 0.287, 
     *              0.309,   0.331,   0.352,   0.373, 0.394, 0.414, 
     *              0.434,   0.454,   0.473,   0.492, 0.511, 0.529, 
     *              0.547,   0.565,   0.583,   0.600, 0.617, 0.634, 
     *              0.650,   0.666,   0.682,   0.697, 0.712, 0.727, 
     *              0.741,   0.755,   0.769,   0.782, 0.796, 0.808, 
     *              0.821,   0.833,   0.845,   0.857, 0.868, 0.879, 
     *              0.890,   0.900,   0.910,   0.920, 0.929, 0.938, 
     *              0.947,   0.956,   0.964,   0.972, 0.979, 1.000/ 
C 
C  HSEQ provides the following coefficients and constants: 
C 
C  For MOBILE4.1, lines a to e are changed to reflect 
C  LDGV and LDGT hot soak rates for pass vehicles. 
C 
C  HS81 = a + b*RVP                 for RVP <  9.0 
C  HS81 = c + d*RVP + e*RVP**2      for RVP >= 9.0 
C 
C  Continuation lines f,g,h,i are the 1981+ LDGV and LDGT1-4 EF HS 
C  temperature correction factor equation coefficients: 
C 
C      TCFHS81 = (f + g*HSTEMP + h*HSTEMP**2)/i 
C 
C                 CARB        TBI         PFI 
      DATA HSEQ1/ 
     a            0.25593   , 0.258327  ,-0.40673  , 
     a           -0.164070  , 0.078327  , 0.078327 , 
     b            0.13823   , 0.041297  , 0.10297  , 
     b            0.13823   , 0.041297  , 0.041297 , 
     c           -1.326657  ,-4.60971   ,18.49588  , 
     c           -5.196600  ,-4.789710  ,-4.789710 , 
     d            0.0       , 0.58219   ,-4.2538   , 
     d            0.69740   , 0.58219   , 0.58219  , 
     e            0.034897  , 0.0       , 0.25072  , 
     e            0.0       , 0.0       , 0.0      / 
C 
      DATA HSEQ2/ 
     f           -2.4636    ,-2.4636    , 0.0      , 
     f           -2.4636    ,-2.4636    , 0.0      , 
     g            0.0       , 0.0       , 0.0055541, 
     g            0.0       , 0.0       , 0.0055541, 
     h            0.00056161, 0.00056161, 0.0      , 
     h            0.00056161, 0.00056161, 0.0      , 
     i            1.31      , 1.31      , 0.46     , 
     i            1.31      , 1.31      , 0.46     / 
C 
C 
C..HSEQ3 and HSEQ4 are for 1981-85 MY cars and CARB Trucks. 
C 
C                 CARB        TBI         PFI 
      DATA HSEQ3/ 
     a            -1.13591  , -0.52111  ,-0.058967 , 
     a             1.29368  , 0.078327  , 0.078327 , 
     b            0.39098   , 0.159322  , 0.100658 , 
     b            0.08904   , 0.041297  , 0.041297 , 
     c           -1.326657  ,-4.60971   ,18.49588  , 
     c           -5.196600  ,-4.789710  ,-4.789710 , 
     d            0.0       , 0.58219   ,-4.2538   , 
     d            0.69740   , 0.58219   , 0.58219  , 
     e            0.034897  , 0.0       , 0.25072  , 
     e            0.0       , 0.0       , 0.0      / 
C 
      DATA HSEQ4/ 
     f           -2.4636    ,-2.4636    , 0.0      , 
     f           -2.4636    ,-2.4636    , 0.0      , 
     g            0.0       , 0.0       , 0.0055541, 
     g            0.0       , 0.0       , 0.0055541, 
     h            0.00056161, 0.00056161, 0.0      , 
     h            0.00056161, 0.00056161, 0.0      , 
     i            2.081      , 1.898     , 0.749    , 
     i            2.541      , 1.31      , 0.46     / 
C 
C 
C 
C..HSEQ5 and HSEQ6 are for 1986+ MY cars and CARB Trucks. 
C 
C..Values for HSEQ5 (a) were modified so that the equations 
C  for RVP < 9.0 match the equations for RVP >= 9.0 at 9.0 RVP. 
C 
C 
C                 CARB        TBI         PFI 
      DATA HSEQ5/ 
     a            -2.5687   , -1.96640  ,-0.226800 , 
     a            -2.8719   , -1.15190  , 0.008200 , 
     b            0.45214   , 0.288530  , 0.08281  , 
     b            0.43908   , 0.178030  , 0.04906  , 
     c           -1.326657  ,-4.60971   ,18.49588  , 
     c           -5.196600  ,-4.789710  ,-4.789710 , 
     d            0.0       , 0.58219   ,-4.2538   , 
     d            0.69740   , 0.58219   , 0.58219  , 
     e            0.034897  , 0.0       , 0.25072  , 
     e            0.0       , 0.0       , 0.0      / 
C 
      DATA HSEQ6/ 
     f           -2.4636    ,-2.4636    , 0.0      , 
     f           -2.4636    ,-2.4636    , 0.0      , 
     g            0.0       , 0.0       , 0.0055541, 
     g            0.0       , 0.0       , 0.0055541, 
     h            0.00056161, 0.00056161, 0.0      , 
     h            0.00056161, 0.00056161, 0.0      , 
     i            2.041      , 2.748      , 0.651   , 
     i            2.527      , 2.596      , 0.805   / 
C 
C 
C 
C  Continuation line code for LDHS, T1AHS, T1BHS, ..., and CLHS 
C  are the MY standard groups. 
C 
C                      At 9.0 RVP                 At 11.5 RVP 
C                       Failed   Failed            Failed   Failed 
C                Pass    Purge Pressure     Pass    Purge Pressure 
      DATA LDHS/ 
     1          14.67,   14.67,   14.67,   22.45,   22.45,   22.45, 
     2          10.63,   14.67,   14.67,   15.68,   22.45,   22.45, 
     3           7.79,   14.50,   14.67,   11.96,   22.19,   22.45, 
     4           2.20,    5.83,    5.89,    4.17,   17.69,   17.90, 
     *           0.00,    0.00,    0.00,    0.00,    0.00,    0.00/ 
C 
      DATA T1AHS/ 
     1          14.67,   14.67,   14.67,   22.45,   22.45,   22.45, 
     2          10.63,   14.67,   14.67,   15.68,   22.45,   22.45, 
     3           7.79,   14.50,   14.67,   11.96,   22.19,   22.45, 
     4           2.20,    5.83,    5.89,    4.17,   17.69,   17.90, 
     *           0.00,    0.00,    0.00,    0.00,    0.00,    0.00/ 
C       
      DATA T1BHS/ 
     1          14.67,   14.67,   14.67,   22.45,   22.45,   22.45, 
     2          10.63,   14.67,   14.67,   15.68,   22.45,   22.45, 
     3           7.79,   14.50,   14.67,   11.96,   22.19,   22.45, 
     4           2.20,    5.83,    5.89,    4.17,   17.69,   17.90, 
     *           0.00,    0.00,    0.00,    0.00,    0.00,    0.00/ 
C 
      DATA T2AHS/ 
     1          18.08,   18.08,   18.08,   27.97,   27.97,   27.97, 
     2           2.20,    5.83,    5.89,    4.17,   17.69,   17.90, 
     *           0.00,    0.00,    0.00,    0.00,    0.00,    0.00, 
     *           0.00,    0.00,    0.00,    0.00,    0.00,    0.00, 
     *           0.00,    0.00,    0.00,    0.00,    0.00,    0.00/ 
C 
      DATA T2BHS/ 
     1          18.08,   18.08,   18.08,   27.97,   27.97,   27.97, 
     2           2.20,    5.83,    5.89,    4.17,   17.69,   17.90, 
     *           0.00,    0.00,    0.00,    0.00,    0.00,    0.00, 
     *           0.00,    0.00,    0.00,    0.00,    0.00,    0.00, 
     *           0.00,    0.00,    0.00,    0.00,    0.00,    0.00/ 
C 
      DATA HD2BHS/ 
     1          18.08,   18.08,   18.08,   27.97,   27.97,   27.97, 
     2           1.69,    9.10,    9.21,    4.41,   27.64,   27.97, 
     *           0.00,    0.00,    0.00,    0.00,    0.00,    0.00, 
     *           0.00,    0.00,    0.00,    0.00,    0.00,    0.00, 
     *           0.00,    0.00,    0.00,    0.00,    0.00,    0.00/ 
C 
      DATA HD3HS/ 
     1          18.08,   18.08,   18.08,   27.97,   27.97,   27.97, 
     2           1.69,    9.10,    9.21,    4.41,   27.64,   27.97, 
     *           0.00,    0.00,    0.00,    0.00,    0.00,    0.00, 
     *           0.00,    0.00,    0.00,    0.00,    0.00,    0.00, 
     *           0.00,    0.00,    0.00,    0.00,    0.00,    0.00/ 
C 
      DATA HD4HS/ 
     1          18.08,   18.08,   18.08,   27.97,   27.97,   27.97, 
     2           1.69,    9.10,    9.21,    4.41,   27.64,   27.97, 
     *           0.00,    0.00,    0.00,    0.00,    0.00,    0.00, 
     *           0.00,    0.00,    0.00,    0.00,    0.00,    0.00, 
     *           0.00,    0.00,    0.00,    0.00,    0.00,    0.00/ 
C 
      DATA HD5HS/ 
     1          18.08,   18.08,   18.08,   27.97,   27.97,   27.97, 
     2           1.69,    9.10,    9.21,    4.41,   27.64,   27.97, 
     *           0.00,    0.00,    0.00,    0.00,    0.00,    0.00, 
     *           0.00,    0.00,    0.00,    0.00,    0.00,    0.00, 
     *           0.00,    0.00,    0.00,    0.00,    0.00,    0.00/ 
C 
      DATA HD6HS/ 
     1          18.08,   18.08,   18.08,   27.97,   27.97,   27.97, 
     2           1.69,    9.10,    9.21,    4.41,   27.64,   27.97, 
     *           0.00,    0.00,    0.00,    0.00,    0.00,    0.00, 
     *           0.00,    0.00,    0.00,    0.00,    0.00,    0.00, 
     *           0.00,    0.00,    0.00,    0.00,    0.00,    0.00/ 
C 
      DATA HD7HS/ 
     1          18.08,   18.08,   18.08,   27.97,   27.97,   27.97, 
     2           1.69,    9.10,    9.21,    4.41,   27.64,   27.97, 
     *           0.00,    0.00,    0.00,    0.00,    0.00,    0.00, 
     *           0.00,    0.00,    0.00,    0.00,    0.00,    0.00, 
     *           0.00,    0.00,    0.00,    0.00,    0.00,    0.00/ 
C 
      DATA HD8AHS/ 
     1          18.08,   18.08,   18.08,   27.97,   27.97,   27.97, 
     2           1.69,    9.10,    9.21,    4.41,   27.64,   27.97, 
     *           0.00,    0.00,    0.00,    0.00,    0.00,    0.00, 
     *           0.00,    0.00,    0.00,    0.00,    0.00,    0.00, 
     *           0.00,    0.00,    0.00,    0.00,    0.00,    0.00/ 
C 
      DATA HD8BHS/ 
     1          18.08,   18.08,   18.08,   27.97,   27.97,   27.97, 
     2           1.69,    9.10,    9.21,    4.41,   27.64,   27.97, 
     *           0.00,    0.00,    0.00,    0.00,    0.00,    0.00, 
     *           0.00,    0.00,    0.00,    0.00,    0.00,    0.00, 
     *           0.00,    0.00,    0.00,    0.00,    0.00,    0.00/ 
C 
      DATA CLHS/ 
     1           4.01,    0.00,    0.00,    6.14,    0.00,    0.00, 
     2           9.01,    0.00,    0.00,   13.79,    0.00,    0.00, 
     3           9.64,    0.00,    0.00,   14.75,    0.00,    0.00, 
     4           9.95,    0.00,    0.00,   15.22,    0.00,    0.00, 
     5           9.90,    0.00,    0.00,   15.15,    0.00,    0.00/ 
C 
      DATA GBUSHS/ 
     1          18.08,   18.08,   18.08,   27.97,   27.97,   27.97, 
     2           1.69,    9.10,    9.21,    4.41,   27.64,   27.97, 
     *           0.00,    0.00,    0.00,    0.00,    0.00,    0.00, 
     *           0.00,    0.00,    0.00,    0.00,    0.00,    0.00, 
     *           0.00,    0.00,    0.00,    0.00,    0.00,    0.00/ 
C 
C  There are 4 coefficeints for each of the 14 hourly intervals defined 
C  during the day. First, the coefficients for week days. 
C 
      DATA HSCOEF / 1143.5,  20.261,  -4.028,  -0.001095, 
     &              1749.6,  24.665,  -4.259,  -0.001225, 
     &              2483.7,  29.051,  -4.449,  -0.000981, 
     &              3212.9,  32.712,  -4.589,  -0.001003, 
     &              4010.7,  36.230,  -4.709,  -0.000929, 
     &              2985.7,  31.546,  -4.552,  -0.001310, 
     &              3208.4,  32.605,  -4.590,  -0.001202, 
     &              4042.0,  36.357,  -4.714,  -0.0009702, 
     &              3066.0,  31.957,  -4.565,  -0.001189, 
     &              3207.6,  32.627,  -4.590,  -0.001167, 
     &              2957.4,  31.546,  -4.549,  -0.001149, 
     &              2435.8,  28.726,  -4.440,  -0.001239, 
     &              2096.7,  26.827,  -4.361,  -0.001445, 
     &              1906.5,  25.712,  -4.306,  -0.000900 , 
C 
C  The second group of coefficients is for the weekends. 
C 
     &              46.7959,47.3457,   0.4006,  0.7685, 
     &              55.31  , 0.0000633,-13.695,  -0.0325, 
     &              2732.31, 0.1819,  -9.620,  -0.000651, 
     &              2208.97, 0.1645,  -9.507,  -0.000750, 
     &              2706.57, 0.1784,  -9.631,  -0.000628, 
     &              2432.10, 0.1674,  -9.591,  -0.000862, 
     &              1824.07, 0.1564,  -9.364,  -0.000857, 
     &              1930.61, 0.1506,  -9.464,  -0.000904, 
     &              2424.95, 0.1761,  -9.531,  -0.000692, 
     &              1921.98, 0.1496,  -9.464,  -0.000861, 
     &              2129.29, 0.1602,  -9.498,  -0.000819, 
     &              1292.06, 0.1333,  -9.193,  -0.001161, 
     &               178.02, 0.0327,  -8.586,  -0.007661, 
     &               520.28, 0.0902,  -8.653,  -0.001729 / 
C 
C  /EVAPACT/: These are the user-supplied Hot Soak Activity Fractions 
C             and Fixed P/P rates and their associated flags. 
C 
      DATA HSACTUSER,DSACTUSER/840*0.0,1296*0.0/ 
      DATA FPPRUSER/4*0./ 
      DATA SAWFLAG,HSACTFLAG,DSACTFLAG,FPPRFLAG/1,0,0,0/ 
C 
C  /EVAPAR/: These are all parameters set during execution by EVP_HRLY & IEVPTR. 
C            Initialization is not necessary, so zeroes are used here. New 
C            parameters (as of Oct 14, 1998) are used in hourly diurnal calculations. 
C 
      DATA IFDS,ISTDC,ISTDP,ISTDT/4*0/,FINJ/2*0.0/ 
C 
C  /EVAPHD/: These values are used in both HS and DU calculations. 
C 
      DATA HDWGT/1.5,2.0/ 
      DATA IVVMAP/1,2,2,2,2,2,2,2,2,2,2,2,2,1,2/ 
C 
C  /EVASTD/: The continuation line codes are the old Mobile5 IVGAS numbers 
C            for MYGSTD (1=LDGV, 2=LDGT12, 3=LDGT34, 4=HDGV, and 5=MC). The 
C            code B is used for the new Gas Bus category. The continuation  
C            code for KSTD81 is IFDS. 
C 
      DATA MYGSTD/ 
     1            1970, 1971, 1977, 1980, 0000, 
     2            1970, 1971, 1977, 1980, 0000, 
     2            1970, 1971, 1977, 1980, 0000, 
     3            1978, 1980, 0000, 0000, 0000, 
     3            1978, 1980, 0000, 0000, 0000, 
     4            1984, NEYR, 0000, 0000, 0000, 
     4            1984, NEYR, 0000, 0000, 0000, 
     4            1984, NEYR, 0000, 0000, 0000, 
     4            1984, NEYR, 0000, 0000, 0000, 
     4            1984, NEYR, 0000, 0000, 0000, 
     4            1984, NEYR, 0000, 0000, 0000, 
     4            1984, NEYR, 0000, 0000, 0000, 
     4            1984, NEYR, 0000, 0000, 0000, 
     5            1977, 1979, 1981, 1984, NEYR, 
     B            1984, NEYR, 0000, 0000, 0000/ 
C 
      DATA KSTD81/ 
     1              5,   5,5,  3,3,  2,2,2,2,2,2,2,2,     2, 
     2              6,   6,6,  4,4,  2,2,2,2,2,2,2,2,     2, 
     3              7,   7,7,  5,5,  2,2,2,2,2,2,2,2,     2/ 
C 
C  /EVAPGR/:  These values are calculated in HOTSOK and DIURNL. 
C 
      DATA GREVP/174*0.0/ 
      DATA IGREVP    /20/ 
C 
C 
C  /ETPPAR/:  Evaporative Test Procedure parameters. 
C 
C  Evap Test Procedure Reduction Percentages 
C  1st line: Problem Free, 2nd: Purge Failure, 3rd: Pressure Failure 
C  (Note: For fail purge and pressure, the hot soak and running loss 
C  reductions are computed in ETPCAL since they are dependent upon 
C  RVP and temperature.) 
C 
C                HS   FDU MDU23  MDU4   PDU  RUNL RESTL 
      DATA ETP/ 50.,  50.,  75.,  40.,  50.,  80.,  75., 
     *          30.,  00.,  00.,  00.,  00.,  30.,  75., 
     *          30.,  00.,  00.,  00.,  00.,  30.,  75./ 
      DATA NETPYR/4/ 
      DATA IETPYR/1996,1997,1998,1999/ 
      DATA ETPPIV/ 30., 55., 90.,100./ 
C 
      END 
