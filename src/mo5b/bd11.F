C Version "@(#)$Id$ $Source$ $Date$ 

      BLOCK DATA BD11
C
C  BLOCK DATA Subprogram 11: /EVAHIA/,/EVAHS1/
C                            /EVAPAR/,/EVAPHD/,/EVASTD/,
C                            /EVAPGR/,/ETPPAR/
C
C  BLOCK DATA 11 initializes the parameters and rates affecting both hot soak
C  and diurnal or hot soak alone (Hot soak HC + Diurnal HC = Evaporative HC).
C  BD11 also intializes the Evaportive Test Procedure (ETP) parameters.
C
C  Common block array subscripts:
C
C  EFHS(3,2,5,5)-  EFHS   ( KEYEQ  , IRVP , ISTD, IVGAS )
C  ETP(7,3)     -  ETP    ( IETP, KEYEQ )
C  ETPPIV(5)    -  ETPPIV ( IETPY )
C  EVP(4)       -  EVP    ( IEVP )
C  FINJ(2)      -  FINJ   ( 2 )
C  GREVP(4,9)   -  GREVP  ( IEVP , IV )
C  GRPD(3,9)    -  GRPD   ( IPARTL, IVLDG )
C  HDWGT(2)     -  HDWGT  ( IHDSAL )
C  HIHS(3,2)    -  HIHS   ( KEYEQ, IRVP )
C  HSEQ(3,2,9)  -  HSEQ   ( IFDS  , IVV  , IEVEQ )
C  IETPYR(5)    -  IETPYR ( IETPY )
C  IVVMAP(5)    -  IVVMAP ( IVGAS )
C  KSTD81(4,3)  -  KSTD81 ( IVTAM , IFDS )
C  MYGSTD(5,5)  -  MYGSTD ( IG    , IVGAS )
C  PARD(3)      -  PARD   ( IPARTL )
C  VGREVP(4)    -  VGREVP ( IEVP )
C
C  Common block dictionary:
C
C   Name   Type              Description
C  ------  ----  ---------------------------------------------
C  /EVAHIA/:
C
C  HIHS     R    High alt, 1972-76, LDGV & LDGT1, HS rates on 2 RVP fuels
C  HIADJ    R    High alt adjustment to low alt based evap rate
C
C  /EVAHS1/:
C
C  HSEQ     R    1981+ LDGV and LDGT1/2 HS equations parameters
C  EFHS     R    EF HS rates for IVGAS vehicles on Pass/Fail and 2 RVP fuels
C
C  /EVAPAR/:
C
C  IFDS     I    See Subscript Dictionary
C  ISTDC    I    Carbureted engine emission standard case
C  ISTDT    I    TBI engine emission standard case
C  ISTDP    I    PFI engine emission standard case
C  FINJ     R    Fuel-injected fraction of vehicle class IVGAS built in MY,
C                index is 1 = TBI, 2 = PBI
C
C  /EVAPHD/:
C
C  HDSAL    R    Sales fraction for HDGV weighing < 14000 lbs
C  HDWGT    R    HDGV weighting factor on value returned from EV81
C  IVVMAP   I    IVGAS to IVV mapping vector
C
C  /EVASTD/:
C
C  MYGSTD   I    Model year groups for emission standards by vehicle type
C  MAXSTD   I    Maximum number of myg's for emission standards, given IV
C  KSTD81   I    Emission standard key for IFDS sensitive 1981+ LDGV/T
C
C  /EVAPGR/:
C
C  EVP      R    Gram evaporative component for a model year
C  GREVP    R    Gram evaporative component for a vehicle type
C  VGREVP   R    Gram evaporative component for all vehicles combined
C  PARD     R    G/Mi partial diurnal case for a model year
C  GRPD     R    G/Mi partial diurnal case for a vehicle type (LDGV/T1/T2)
C
C  /ETPPAR/:
C
C  ETP      R    Evaporative Test Procedure rates for Hot Soak, Full
C                Diurnal, Multiple Diurnal 2-3, Multiple Diurnal 4+,
C                Partial Diurnal, Running Loss and Resting Loss for
C                Problem Free, Purge Failure and Pressure Failure
C  NEWETP   I    A flag which tells if the user input new ETP rates (NEWETP=2)
C  IETPYR   I    The phase-in years of the ETP program.
C  ETPPIV   I    The phase-in values (percentages) of the ETP program.
C  NETPYR   I    The number of ETP phase-in years
C
C  Local array subscripts:
C
C  LDHS(3,2,5)  -  LDHS    ( KEYEQ, IRVP , ISTD  )
C  T1HS(3,2,5)  -  T1HS    ( KEYEQ, IRVP , ISTD  )
C  T2HS(3,2,5)  -  T2HS    ( KEYEQ, IRVP , ISTD  )
C  HDHS(3,2,5)  -  HDHS    ( KEYEQ, IRVP , ISTD  )
C  CLHS(3,2,5)  -  CLHS    ( KEYEQ, IRVP , ISTD  )
C  HSEQ1(3,2,5) -  HSEQ1   ( IFDS  , IVV  , IEVEQ )
C  HSEQ2(3,2,4) -  HSEQ2   ( IFDS  , IVV  , IEVEQ )
C
C  Local array dictionary:
C
C   Name   Type                Description
C  ------  ----  -----------------------------------------------
C  LDHS     R    EF HS rates for pre-81 LDGV
C  T1HS     R    EF HS rates for pre-81 LDGT1
C  T2HS     R    EF HS rates for pre-81 LDGT2
C  HDHS     R    EF HS rates for pre-85 HDGV
C  CLHS     R    EF HS rates for MC
C  HSEQ1    R    1981+ LDGV and LDGT1/2 HS equations parameters
C                for IEVEQ = 1-5
C  HSEQ2    R    1981+ LDGV and LDGT1/2 HS temperature equations parameters
C                for IEVEQ = 6-10
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
C
      INCLUDE 'ETPPAR.I'
      INCLUDE 'EVAHIA.I'
      INCLUDE 'EVAHS1.I'
      INCLUDE 'EVAPAR.I'
      INCLUDE 'EVAPGR.I'
      INCLUDE 'EVAPHD.I'
      INCLUDE 'EVASTD.I'
      INCLUDE 'PARAME.I'
C
      REAL LDHS
C
      DIMENSION LDHS(3,2,5),T1HS(3,2,5),
     *          T2HS(3,2,5),HDHS(3,2,5),CLHS(3,2,5)
      DIMENSION HSEQ1(3,2,5),HSEQ2(3,2,4)
C
      EQUIVALENCE
     *     (EFHS(1,1,1,1),LDHS(1,1,1)),
     *     (EFHS(1,1,1,2),T1HS(1,1,1)),
     *     (EFHS(1,1,1,3),T2HS(1,1,1)),
     *     (EFHS(1,1,1,4),HDHS(1,1,1)),
     *     (EFHS(1,1,1,5),CLHS(1,1,1)),
     *     (HSEQ(1,1,1),HSEQ1(1,1,1)),
     *     (HSEQ(1,1,6),HSEQ2(1,1,1))
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
C  HSEQ provides the following coefficients and constants:
C
C  For MOBILE4.1, lines a to e are changed to reflect
C  LDGV and LDGT hot soak rates for pass vehicles.
C
C  HS81 = a + b*RVP                 for RVP <  9.0
C  HS81 = c + d*RVP + e*RVP**2      for RVP >= 9.0
C
C  Continuation lines f,g,h,i are the 1981+ LDGV and LDGT1/2 EF HS
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
C  Continuation line code for LDHS, T1HS, T2HS, HDHS, and CLHS
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
      DATA T1HS/
     1          14.67,   14.67,   14.67,   22.45,   22.45,   22.45,
     2          10.63,   14.67,   14.67,   15.68,   22.45,   22.45,
     3           7.79,   14.50,   14.67,   11.96,   22.19,   22.45,
     4           2.20,    5.83,    5.89,    4.17,   17.69,   17.90,
     *           0.00,    0.00,    0.00,    0.00,    0.00,    0.00/
C
      DATA T2HS/
     1          18.08,   18.08,   18.08,   27.97,   27.97,   27.97,
     2           2.20,    5.83,    5.89,    4.17,   17.69,   17.90,
     *           0.00,    0.00,    0.00,    0.00,    0.00,    0.00,
     *           0.00,    0.00,    0.00,    0.00,    0.00,    0.00,
     *           0.00,    0.00,    0.00,    0.00,    0.00,    0.00/
C
      DATA HDHS/
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
C  /EVAPAR/: These are all parameters set during execution by CCEVRT & IEVPTR.
C            Initialization is not necessary, so zeroes are used here.
C
      DATA IFDS,ISTDC,ISTDP,ISTDT/4*0/,FINJ/2*0.0/
C
C  /EVAPHD/: These values are used in both HS and DU calculations.
C
      DATA HDSAL/0.875/, HDWGT/1.5,2.0/, IVVMAP/1,2,2,2,1/
C
C  /EVASTD/: The continuation line codes are IVGAS for MYGSTD and IFDS for
C            KSTD81.
C
      DATA MYGSTD/
     1            1970, 1971, 1977, 1980, 0000,
     2            1970, 1971, 1977, 1980, 0000,
     3            1978, 1980, 0000, 0000, 0000,
     4            1984, NEYR, 0000, 0000, 0000,
     5            1977, 1979, 1981, 1984, NEYR/
C
      DATA MAXSTD/5/
C
C                 LDGV LDGT1 LDGT2 HDGV
      DATA KSTD81/
     1              5,   5,    3,   2,
     2              6,   6,    4,   0,
     3              7,   7,    5,   0/
C
C  There is no fuel injected case for HDGV presently, thus the 0 in KSTD81.
C
C
C  /EVAPGR/:  These values are calculated in HOTSOK and DIURNL.
C
      DATA EVP/4*0.0/, GREVP/36*0.0/, VGREVP/4*0.0/
      DATA PARD/3*0.0/,GRPD/27*0.0/
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
      DATA NEWETP,NETPYR/1,4/
      DATA IETPYR/1996,1997,1998,1999/
      DATA ETPPIV/ 20., 40., 90.,100./
C
      END
