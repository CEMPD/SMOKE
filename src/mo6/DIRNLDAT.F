      MODULE DIRNLDAT 
C 
C    This module initializes the variables for hourly diurnal 
C    calculations 
C       
C    This data includes the vehicle activity fractions and parameters 
C    for diurnal fractions, daily diurnals and coefficients for 
C    equations to calculate ratios R12/R23 
C 
C 
C  Changes:  (Last change first) 
C 
C  27 Jan 01 @EPA-BAG bug314 removed HDGV_FCTR1 & HDGV_FCTR2; diurnals using factors 
C            from EVAPHD.I 
C  21 Sep 00 @EPA-BAG 263 Updated evaporative test procedure (ETP) phase in fractions 
C            for 1996 and 1997 
C  11 Sep 00 @EPA-BAG 252&253 Removed partial day diurnal fraction 
C            normalization variable NRMLIZ_ID(IFUEL_SYS,ISOAK) 
C  20 Jul 00 @EPA-elg E10 Tier2 MDPV; added DV5 initialization. 
C  07 Aug 00 @EPA-bag bug215 Added DIURNL24B array; they are weighted 1999 portion 
C            of the full day diurnals for the model years 1996-1998; DIURNL24A are 
C            are used for all model years and are also the 1995 portion of the full 
C            day diurnals for the model years 1996-1998 
C  24 JuL 00 @EPA-BAG D15 Added multi-day ratios from diuratio.for 
C  28 Jun 00 @EPA-ddj D15 Changed dimensions for R12COEFF and R23COEFF; put 
C            DIU1 and DIU2 arrays in this module for multi-day diurnals 
C            from 5,4 to 4,5. 
C  07 Jan 00 @Dyntel-bg 1-005 updated diurnal code(new) 
C 
C 
C  MODULE dictionary: 
C 
C     Name           Type              Description 
C    ------          ----  --------------------------------------------- 
C  NVHCLES_SOAKING    R    The activity weight/fraction is a normalized 
C                          distribution of vehicles soaking for a given hour 
C                          of the day and the length of the soak time. 
C  IFUEL              I    fuel system type :    carbureted,    IFUEL=1 
C                                                fuel injected, IFUEL=2 
C  LDGVREDUCT         R    tier2 phase in 
C  LDGTREDUCT         R    tier2 phase in 
C  LDGV_FCTR          R    tier2 phase in 
C  HDGV_FCTR1         R    heavy duty factor for HDGV2B and HDGV3 (=1.5) 
C  HDGV_FCTR2         R    heavy duty factor for HDGV4,5,6,7,8a,8b, and HDGB (=2.0) 
C  DIURNL_FRCTN       R    diurnal fractions 
C  DIURNL24A          R    weighted 24hr/daily diurnal 
C  DIURNL24B          R    weighted 24hr/daily diurnal needed for mixed years, 1996-1998 
C  R12COEFF           R    equation coefficients for R12 
C  R23COEFF           R    equation coefficients for R23 
C  DIU1               R    temporary storage for unweighted daily diurnal results 
C  DIU2               R    temporary storage for unweighted daily diurnal results 
C                          needed for mixed years, 1996-1998 
C  DIURFRA            R    parameter for diurnal fractions 
C  DIURFRB            R    parameter for diurnal fractions 
C  DIURFRC            R    parameter for diurnal fractions 
C  DIURFRD            R    parameter for diurnal fractions 
C  DIURFRE            R    parameter for diurnal fractions 
C  DIURFRF            R    parameter for diurnal fractions 
C  DIURFRG            R    parameter for diurnal fractions 
C  DIURFRH            R    parameter for diurnal fractions 
C  DIURFRI            R    parameter for diurnal fractions 
C  DIURFRJ            R    parameter for diurnal fractions 
C  DIURFRK            R    parameter for diurnal fractions 
C  DIURFRL            R    parameter for diurnal fractions 
C  COMC               R    parameters for the motorcycle daily diurnal equation 
C  COGLL              R    parameters for liquid leak daily diurnal equation 
C  COLDV              R    parameters for pre-1972 daily diurnal equation 
C  COLDVA             R    parameters for the daily diurnal equation 
C  COLDVB             R    parameters for the daily diurnal equation 
C  COLDVC             R    parameters for the daily diurnal equation 
C  COLDVD             R    parameters for the daily diurnal equation 
C  W95                R    for 1996<=MY<=1998 weight of 1995 daily diurnal 
C  W99                R    for 1996<=MY<=1998 weight of 1999 daily diurnal 
C  LDGV_FCTR          R    daily diurnal multiplicative factor for light duty gas vehicles, 
C                          LDGV_FCTR=1.0 
C  HDGV_FCTR1         R    daily diurnal multiplicative factor for heavy duty gas vehicles 
C                          2b and 3, HDGV1_FCTR=1.5 
C  HDGV_FCTR2         R    daily diurnal multiplicative factor fro heavy duty gas vehicles 
C                          4,5,6,7,8a,8b, and gas bus, HDGV2_FCTR=2.0 
C  R12                R    Ratio of the 2nd day's diurnal 
C                          emissions to the 1st day's 
C  R23                R    Ratio of the 3rd day's diurnal 
C                          emissions to the 2nd day's 
C  R12Y96A            R    First ratio of the 2nd day's diurnal 
C                          emissions to the 1st day's for phase-in 
C                          model years 1996-98 
C  R12Y96B            R    Second ratio of the 2nd day's diurnal 
C                          emissions to the 1st day's for phase-in 
C                          model years 1996-98 
C  R23Y96A            R    First ratio of the 3rd day's diurnal 
C                          emissions to the 2nd day's for phase-in 
C                          model years 1996-98 
C  R23Y96B            R    Second ratio of the 3rd day's diurnal 
C                          emissions to the 2nd day's for phase-in 
C                          model years 1996-98 
C 
C  Local array subscripts: 
C 
C  COLDVNEWPR(2) 
C  COLDVNEWPU(2) 
C  W95(3)                  -    WA(MY-1995); 1996<=MY<=1998 
C  W99(3)                  -    WB(MY-1995); 1996<=MY<=1998 
C  COMC(3)                -    3 parameters for the motorcycle daily diurnal equation 
C  COGLL(2)               -    2 parameters for liquid leak daily diurnal equation 
C  COLDV(3)               -    3 parameters for pre-1972 daily diurnal equation 
C  COLDVA(4,3,2)          -    (MYRNG,     IPPS,  IFUEL) 
C  COLDVB(4,3,2)          -    (MYRNG,     IPPS,  IFUEL) 
C  COLDVC(4,3,2)          -    (MYRNG,     IPPS,  IFUEL) 
C  COLDVD(4,3,2)          -    (MYRNG,     IPPS,  IFUEL) 
C  COLDVNEWPP(3) 
C  DIURFRA(7)             -    ( IFUEL_SYS             ) 
C  DIURFRB(7)             -    ( IFUEL_SYS             ) 
C  DIURFRC(7)             -    ( IFUEL_SYS             ) 
C  DIURFRD(7)             -    ( IFUEL_SYS             ) 
C  DIURFRE(7)             -    ( IFUEL_SYS             ) 
C  DIURFRF(7)             -    ( IFUEL_SYS             ) 
C  DIURFRG(7)             -    ( IFUEL_SYS             ) 
C  DIURFRH(7)             -    ( IFUEL_SYS             ) 
C  DIURFRI(7)             -    ( IFUEL_SYS             ) 
C  DIURFRJ(7)             -    ( IFUEL_SYS             ) 
C  DIURFRK(7)             -    ( IFUEL_SYS             ) 
C  DIURFRL(7)             -    ( IFUEL_SYS             ) 
C  DIURNL24(7)            -    ( IFUEL_SYS             ) 
C  DIURNL_FRCTN(19,20,7)  -    ( IH,       ISOAK, IFUEL_SYS) 
C  NVHCLES_SOAKING(18,72) -    ( IH   ,    IHR_SOAK    ) 
C  R12(7)                 -    ( IFUEL_SYS             ) 
C  R23(7)                 -    ( IFUEL_SYS             ) 
C  R12Y96A(7)             -    ( IFUEL_SYS             ) 
C  R12Y96B(7)             -    ( IFUEL_SYS             ) 
C  R23Y96A(7)             -    ( IFUEL_SYS             ) 
C  R23Y96B(7)             -    ( IFUEL_SYS             ) 
C 
      IMPLICIT NONE 
C 
C******* NVHCLES_SOAKING(IHR=1-18,IHR_SOAK=1-72) *********** 
C 
C Normalized vehicle activity distributions for diurnal emissions  
C used for both vapor leakers and liquid leakers (see report M6.EVP.002 
C and associated spreadsheets, EVP_DIU2.XLS, from which this data was 
C taken) 
C 
C      NVHCLES_SOAKING(IHR,IHR_SOAK) : 
      REAL, DIMENSION(18,72) :: NVHCLES_SOAKING 
C 
C The coefficient for hourly diurnal fraction 
C 
      REAL, DIMENSION(7), SAVE :: DIURFRA = 
C 
C         FPRES     FPURG       PBOTH 
     * (/0.010549,  0.006724,  0.007032,  ! CARB 
     *   0.006515,  0.007882,  0.008001,  ! FI 
     *   0.021349/)                       ! liquid leaks 
C 
      REAL, DIMENSION(7), SAVE :: DIURFRB = 
C 
C          FPRES      FPURG      PBOTH 
     *  (/0.001758,  0.003966,  0.003586, ! CARB 
     *   0.001963,  -0.001172,  0.001961, ! FI 
     *   0.002065/)                      ! liquid leaks 
C 
      REAL, DIMENSION(7), SAVE :: DIURFRC = 
C 
C          FPRES       FPURG      PBOTH 
     * (/ 0.000000, -0.000018,  0.000000,  ! CARB 
     *    0.000000,  0.000084,  0.000056,  ! FI 
     *    0.000000/)                       ! liquid leaks 
C 
      REAL, DIMENSION(7), SAVE :: DIURFRD = 
C 
C          FPRES     FPURG      PBOTH 
     * (/0.000000,   0.000000,  0.000000,  ! CARB 
     *   0.000000,   0.000000,  0.000000,  ! FI 
     *   0.010137/)                        ! liquid leaks 
C 
      REAL, DIMENSION(7), SAVE :: DIURFRE = 
C 
C          FPRES      FPURG      PBOTH 
     * (/0.000000,   0.000000,  0.000000,  ! CARB 
     *   0.000000,   0.000000,  0.005964,  ! FI 
     *   0.000000/)                        ! liquid leaks 
C 
      REAL, DIMENSION(7), SAVE :: DIURFRF = 
C 
C          FPRES      FPURG      PBOTH 
     * (/0.000000,   0.000000,  0.000000,  ! CARB 
     *   0.000000,   0.000118,  0.000000,  ! FI 
     *   0.000000/)                        ! liquid leaks 
C 
      REAL, DIMENSION(7), SAVE :: DIURFRG = 
C 
C          FPRES      FPURG      PBOTH 
     * (/0.001138,   0.000000,  0.000000, ! CARB 
     *   0.001194,   0.000855,  0.000535, ! FI 
     *   0.000000/)                       ! liquid leaks 
C  
      REAL, DIMENSION(7), SAVE :: DIURFRH = 
C 
C          FPRES      FPURG      PBOTH 
     * (/0.000000,   0.000023,  0.000023, ! CARB 
     *   0.000000,   0.000000,  0.000000, ! FI 
     *   0.000000/)                       ! liquid leaks 
C 
      REAL, DIMENSION(7), SAVE :: DIURFRI = 
C 
C          FPRES     FPURG      PBOTH 
     * (/ 0.000000, -0.001122, -0.001111, ! CARB 
     *    0.000574,  0.000825,  0.000000, ! FI 
     *    0.000000/)                     ! liquid leaks 
C 
      REAL, DIMENSION(7), SAVE :: DIURFRJ = 
C 
C         FPRES     FPURG      PBOTH 
     * (/0.000000, 0.000000,  0.000000, ! CARB 
     *   0.000000, 0.000000, -0.000060, ! FI 
     *   0.000000/)                     ! liquid leaks 
CC 
      REAL, DIMENSION(7), SAVE :: DIURFRK = 
C 
C         FPRES      FPURG      PBOTH 
     * (/0.001765,  0.000000,  0.000000, ! CARB 
     *   0.001329,  0.006960,  0.000000, ! FI 
     *   0.000000/)                       ! liquid leaks 
C 
      REAL, DIMENSION(7), SAVE :: DIURFRL = 
C 
C         FPRES     FPURG      PBOTH 
     * (/0.000000,  0.000019,  0.000000,  ! CARB 
     *   0.000000, -0.000160,  0.000000,  ! FI 
     *   0.000000/)                       ! liquid leaks 
C 
C 
C     Daily Diurnal Emission Coefficients 
C 
C     liquid leaks, all model years, all vehicle types 
      REAL, DIMENSION(2), SAVE :: COGLL = (/ 40.10367, 2.616201/) 
C 
C     pre-1972; all fuel system categories 
      REAL, DIMENSION(3), SAVE :: COLDV = (/-40.67512,0.05115,1.41114/) 
C 
C     MY(MYRNG)             72-79(MYRNG=1)      80-85(MYRNG=2)      86-95(MYRNG=3)    99+(MYRNG=4) 
C     PURGE/Pressure(IPPS)  Fail Press(IPPS=1)  Fail Purge(IPPS=2)  Pass Both(IPPS=3) 
C     FUEL SYSTEM(IFUEL)    CARB(IFUEL=1)       FI(IFUEL=2) 
C     COLDVA(MYRNG, IPPS,IFUEL) 
      REAL, DIMENSION(4,3,2), SAVE :: COLDVA = RESHAPE( 
C         1972-1979  1980-1985  1986-1995  1996- 
C 
C  IPPS 
     1  (/6.90895,   9.71190,   1.51716,   0.47846,   !  \ 
     2   -4.58719,   6.88852,   4.11975,   3.25800,   !   carbureted vehicles(IFUEL=1) 
     3   -8.09426,   2.97845,   1.42298,   0.19415,   !  / 
C 
     1    6.90895,  -3.14389,   0.47846,   0.47846,   !  \ 
     2   -4.58719,   2.39612,   3.25800,   3.25800,   !   fuel injected vehicles(IFUEL=2) 
     3   -8.09426,  -1.29432,   0.38830,   0.19415/), !  / 
     & (/4,3,2/)) 
C 
      REAL, DIMENSION(4,3,2), SAVE :: COLDVB = RESHAPE( 
C 
     1  (/0.00000,   0.00000,   0.00000,   0.00000,   !  \ 
     2    0.00000,   0.00000,   0.00000,   0.00000,   !   carbureted vehicles(IFUEL=1) 
     3    0.00000,   0.00000,   0.00000,   0.00000,   !  / 
C 
     1    0.00000,   0.03337,   0.00000,   0.00000,   !  \ 
     2    0.00000,   0.00000,   0.00000,   0.00000,   !   fuel injected vehicles(IFUEL=2) 
     3    0.00000,   0.00000,   0.00000,   0.00000/), !  / 
     & (/4,3,2/)) 
C 
      REAL, DIMENSION(4,3,2), SAVE ::  COLDVC = RESHAPE( 
C 
     1 (/0.02955,    0.02156,   0.02156,   0.01497,    ! \ 
     2   0.04825,    0.00000,   0.00000,   0.00941,    !  carbureted vehicles(IFUEL=1) 
     3   0.04825,    0.00000,   0.00000,   0.00252,    ! / 
C 
     1   0.02955,    0.00000,   0.01497,   0.01497,    ! \ 
     2   0.04825,    0.00954,   0.00941,   0.00941,    !  fuel injected vehicles(IFUEL=2) 
     3   0.04825,    0.00954,   0.00503,   0.00252/),  ! / 
     & (/4,3,2/)) 
C 
      REAL, DIMENSION(4,3,2), SAVE ::  COLDVD = RESHAPE( 
C 
     1 (/0.00000,    0.00000,   0.00000,   0.00000,    ! \ 
     2   0.00000,    0.01896,   0.01896,   0.00000,    !  carbureted vehicles(IFUEL=1) 
     3   0.00000,    0.01464,   0.01464,   0.00000,    ! / 
C 
     1   0.00000,    0.00000,   0.00000,   0.00000,    ! \ 
     2   0.00000,    0.00000,   0.00000,   0.00000,    !  fuel injected vehicles(IFUEL=2) 
     3   0.00000,    0.00000,   0.00000,   0.00000/),  ! / 
     & (/4,3,2/)) 
C 
C     coefficients for 1996 <= MY <= 1999; weights of 1995 and 1999 daily diurnals 
      REAL, DIMENSION(3), SAVE :: W95 = (/0.70,0.45,0.10/) 
C 
      REAL, DIMENSION(3), SAVE :: W99 = (/0.30,0.55,0.90/) 
C 
C     coefficients for motorcycle daily diurnals 
      REAL, DIMENSION(3), SAVE :: COMC=(/-36.7971, 0.058251, 0.855491/) 
C 
C Equation coefficients for ratios R12 and R23 
C 
C                   Coefficients by Strata 
C 
C               Failing   Passing   Failing   Passing 
C               either     both     either     both 
C                purge     purge     purge     purge 
C               pressure  pressure  pressure  pressure 
C 
C                 CARB      CARB       FI        FI 
C 
      REAL, DIMENSION(4,5), SAVE :: R12COEFF = RESHAPE( 
     A          (/50.19    ,49.06     ,48.61    ,47.48, 
     B            -0.70    ,-0.70     ,-0.70    ,-0.70, 
     C            -8.11    ,-8.11     ,-8.11    ,-8.11, 
     D             0.12    , 0.12     , 0.12    , 0.12, 
     E             0.64    , 0.74     , 0.64    , 0.74/), 
     & (/4,5/)) 
C 
      REAL, DIMENSION(4,5), SAVE :: R23COEFF = RESHAPE( 
     A          (/12.25    ,12.25     ,12.25    ,12.25, 
     B            -0.21    ,-0.21     ,-0.21    ,-0.21, 
     C            -2.61    ,-2.61     ,-2.61    ,-2.61, 
     D             0.04    , 0.04     , 0.04    , 0.04, 
     E             0.81    , 0.81     , 0.81    , 0.81/), 
     & (/4,5/)) 
C 
      REAL, SAVE :: LDGVREDUCT = 0.525 
      REAL, SAVE :: LDGTREDUCT = 0.400 
      REAL, SAVE :: LDGV_FCTR =  1.0 
C 
      REAL, DIMENSION(7) :: DIURNL24A = 0.0 
      REAL, DIMENSION(7) :: DIURNL24B = 0.0 
      REAL, DIMENSION(7) :: DIU1 = 0.0 
      REAL, DIMENSION(7) :: DIU2 = 0.0 
      REAL, DIMENSION(6:24,20,7) :: DIURNL_FRCTN = 0.0 
      REAL, DIMENSION(7,20) :: NRMLIZ_FD = 0.0 
C 
      REAL, DIMENSION(7) :: R12 = 0.0 
      REAL, DIMENSION(7) :: R23 = 0.0 
      REAL, DIMENSION(7) :: R12Y96A = 0.0 
      REAL, DIMENSION(7) :: R12Y96B = 0.0 
      REAL, DIMENSION(7) :: R23Y96A = 0.0 
      REAL, DIMENSION(7) :: R23Y96B = 0.0 
      REAL, DIMENSION(25,2,3)  :: DV5 = 0.0 
C 
      END MODULE DIRNLDAT 
 
