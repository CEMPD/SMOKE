      BLOCK DATA BD_IM1
C
C  BLOCK DATA IM initializes data used in the calculation of IM credits
C  for 1981-1993 model year light duty gasoline vehicles and trucks.
C
C  BLOCK DATA Subprogram IM :  /IMPAR1/,/IMPAR2/,/IM240P/
C
C   Changes.
C  15 Oct 02 @EPA-BAG Bug399 : OBDIM changed back to 0.99; mistakenly changed
C            to 0.9 via task 32
C  23 Sep 02 AIR Task 35B: Added OBDIMFLG to prevent using OBD EXH RESPONSE
C            whenever an I/M program is specified. Also corrected comments.
C  26 Jul 02 AIR Task 32: OBDIM chanegd to 0.9 (was 0.99)
C  04 Apr 02 AIR Task 32: added OBD EXHAUST RESPONSE and OBD EVAP RESPONSE
C            default values.
C  04 Jan 02 @EPA-BAG        smooth BIENRATIO versus vehicle age
C  17 Dec 01 @EPA-BAG        Remove sawtooth
C  17 Oct 01 @EPA-BAG bug331 Biennial/ annual I&M program discrepancy
C   25 May 01 @EPA-ELG I&M parameters for loaded idle
C   06 Oct 00 @EPA-BAG BUG155 Remove OBDPC
C   15 Sep 00 AIR Task 08: Increased all I/M related arrays from 5 to 7
C             to allow for Evap I/M programs. Included/changed initialization
C             of EXEMPTAGE to MAXYRS, IMGRACE to 1, CRIM to 85.0, WAVPRE81
C             and WAV81PLUS to 0.05, WAVRDC to 0.80, WMIL TO 0.85, OBDIM to 0.99,
C             OBDF to 0.9,0.1,0.0, OBDREP to 1.5, OBDPC to 0.
C   07 Jul 00 @EPA-elg           Added ASMRAT array.
C   14 Jun 00 @EPA-elg           Updated array R_S_CONV with new values.
C   03 Aug 99 @DynTel-BG  2-663  Renamed BD43.FOR to BD_IM1.FOR and included
C             I&M parameters from BD27.FOR. Fixed typo in RSDCSE(2,3,1); it
C             was 2% now 7%.  Commments
C   29 Apr 99 @DynTel-ddj 2-694, Added data for array R_S_CONV.
C   22 Sep 98 @DynTel-RJD 2-663  Created BD43 to initialize 81-93 IM
C             benefits data.
C  Common block array subscripts:
C
C    /IM240P/:
C    DSIZE(25,14)      -  DSIZE  ( IDX, IVTAM )
C
C    /IMPAR1/:
C    ICYIM(2,7)       -  ICYIM  ( 1st/last, I/M Program )
C    ISTRIN(7)        -  ISTRIN ( I/M Program )
C    MODYR(2,7)       -  MODYR  ( 1st/ Last, I/M Program )
C    WAVPRE81(7)      -  WAVPRE81  ( IMPGM )
C    WAV81PLUS(7)     -  WAV81PLUS ( IMPGM )
C    EXEMPTAGE(7)     -  EXEMPTAGE ( IMPGM )
C    IMGRACE(7)       -  IMGRACE   ( IMPGM )
C    CRIM(7)          -  CRIM      ( IMPGM )
C
C    /IMPAR2/:
C    ILDT(14,7)        -  ILDT   ( IVIM, IM PROGRAM )
C    ITEST(7)          -  ITEST  ( I/M Program )
C
C     /LDGIMC/:
C     ASEBZML(13,3,5,3) - ASEBZML(IGY,IP,IV,IFDS)
C     ASMCP1 (4,3,6)    - ASMCP1 (ICOEF,IP,ITEST-5)
C     ASMCP2 (4,3,6)    - ASMCP2 (ICOEF,IP,ITEST-5)
C     IDRID1 (3,2,3)    - IDRID1 (IFDS,ISR,IP)
C     IDRID2 (3,2,3)    - IDRID2 (IFDS,ISR,IP)
C     IMPEFF (7,3)      - IMPEFF (IMPGM,IP)
C
C     /LDGOBD/:
C     R_S_CONV(11,3,2)  - R_S_CONV(ISTD, IP, ISR)
C     OBDHINIM(25,3,5)     - OBDHINIM   (JDX,IP,IV)
C     OBDHI_IM(25,3,5)     - OBDHI_IM   (JDX,IP,IV)
C     OBDXM(2)
C     OBDXF(3)
C     OBDVM(3)
C     OBDVF(3)
C
C     /LDGSIM/:
C     RSDEFF (15,3)      - RSDEFF (ICT,IP)
C     RSDCSE (4,3,2)     - RSDCSE (ICT,IP,IMSL)
C     RSDCUT (7)         - RSDCUT (IMPGM)
C     RSDNOX (7)         - RSDNOX (IMPGM)
C
C
C  Common block dictionary:
C
C  /IMPAR1/:
C
C  EXEMPTAGE  I    I&M programs can have vehicle ages at which they no longer are required
C                  to undergo testing; EXEMPTAGE(IPGM) is this vehicle age when it is no longer
C                  required to be tested. It can have values greater than or equal to 15 years.
C  ICYIM      I    I/M program start and end years
C  IMGRACE    I    I&M programs for new vehicles have a variable age at which the vehicle
C                  is required to start being tested for high emissions. IMGRACE(IMPGM) is
C                  that age.
C  ISTRIN     I    Technology I/II stringency level of I/M program
C  MODYR      I    first (earliest) and last (latest) year included in the I/M pro
C  WAIVER     R    waiver rate - percent of reviewed fleet that is exempted from
C                  I/M testing
C  CRIM       R    compliance rate/enforcement level - percent of fleet reviewed for I/M
C                  testing WAIVER is based on CRIM: the percent of the total fleet
C                  tested = (1.0 - WAIVER(IGW)) * (CRIM * fleetsize),
C                  IGW = 1 for pre-1981 and IGW = 2 for 1981+ MY vehicles.
C  EVFBEN     R    Evap fractional benefit
C
C  /IMPAR2/:
C
C  ILDT     I    inclusion vector - vehicle classes to be covered by I/M:
C                    1 = omit    2 = included
C  ITEST    I    I&M test type: (1-19 are exhaust, 20+ are evap)
C                1 = idle test
C                2 = 2500 rpm idle test
C                3 = loaded idle test
C                4 = IM240 test
C                5 = ASM test procedure 
C                6 = ASM 2525 PHASE-IN
C                7 = ASM 2525 FINAL
C                8 = ASM 5015 PHASE-IN
C                9 = ASM 5015 FINAL
C               10 = ASM 2525/5015 PHASE-IN
C               11 = ASM 2525/5015 FINAL
C               12 = OBD I/M
C            13-19 = for future exhaust I/M use
C               20 = EVAP OBD
C               21 = EVAP OBD & GC
C               22 = EVAP FP & GC
C               23 = EVAP GC
C
C  /IM240P/:
C  DSIZE    R    non-ATP tampering rates captured in EVPEME
C
C /LDGIMC/:
C  ASEBZML   R    after repair start zero mile level for high emitters
C                 The after repair emissions line has zero slope. So,
C                 the 'deterioration' is a constant, ASEBZML(IGY,IP,IV,IFDS),
C                 over mileage. (See table 6 from M6.IM.001, p.53.)
C  ASMCP1    R    coefficients used in the calculation of the high
C                 emitter identification rate for vehicle ages less
C                 than or equal to 10 years
C  ASMCP2    R    coefficients used in the calculation of the high
C                 emitter identification rate for vehicle ages greater
C                 than 10 years
C  IDRID1    R    high emitter identification rate for idle test
C  IDRID2    R    high emitter identification rate for idle/2500 test
C
C /LDGSIM/
C  RSDCSE    R    RSD clean screening effectiveness fractions
C  RSDCUT    R    RSD CO cutpoint
C  RSDEFF    R    probability that a remote sensing device will detect a
C                 high emitting vehicle that passes over it.
C  RSDNOX    R    NOX cutpoint for the RSD clean screening program
C
C /LDGOBD/
C  R_S_CONV  R    Conversion factors from FTP to running/start emissions
C                 rates.
C  OBDHINIM  R    fraction of high emitting vehicles with OBD and without I&M program
C  OBDHI_IM  R    fraction of high emitting vehicles with OBD and I&M program
C  OBDXM     R    OBD owner exhaust mileage inflection points (in 10,000 miles)
C  OBDXF     R    OBD owner exhaust response fractions
C                 OBDXF(1) holds owner response fraction for 
C                          vehicle mileage =< OBDXM(1)
C                 OBDXF(2) holds owner response fraction when 
C                          when OBDXM(1) < vehicle mileage =< OBDXM(2)
C                 OBDXF(3) holds owner response fraction when 
C                          vehicle mileage > OBDXM(2)
C  OBDXFLG   L    flag to determine if OBD EXH RESPONSE has been properly entered
C  OBDXLIN   L    flag to determine if Step (.false.) or Linear (.true.) function used
C                 for mileage vs fraction
C  OBDVA     I    OBD owner evap age inflection points (in years)
C  OBDVF     R    OBD owner evap response fractions
C                 OBDVF(1) holds owner response fraction for 
C                          vehicle age =< OBDVM(1)
C                 OBDVF(2) holds owner response fraction when 
C                          when OBDVM(1) < vehicle age =< OBDVM(2)
C                 OBDVF(3) holds owner response fraction when 
C                          vehicle age > OBDVM(2)
C  OBDVFLG   L    flag to determine if OBD EVAP RESPONSE has been properly entered
C  OBDVLIN   L    flag to determine if Step (.false.) or Linear (.true.) function used
C                 for age vs fraction
C  WMIL      R    fraction of vehicles with on board emissions malfunction
C                 indicator light ON; default in REINIT.FOR is set at 85% of
C                 high emitting vehicles will have properly functioning
C                 emission malfunction indicators
C  OBDIM     R    OBD I&M program owner response fractions; default value set
C                 at 90% in REINIT.FOR
C  OBDREP    R    Emission level after an OBD-induced repair; the default
C                 value is 1.5, i.e., 1.5 time the applicable 50 kmile
C                 certification standard (See IMPROG.FOR)
C  OBDIMFLG  L    Flag which aids in preventing the OBD EXH RESPONSE command from
C                 being used whenver an I/M program has been entered.
C
      IMPLICIT NONE
      INCLUDE 'IMPAR1.I'
      INCLUDE 'IMPAR2.I'
      INCLUDE 'IM240P.I'
      INCLUDE 'LDGIMC.I'
      INCLUDE 'LDGOBD.I'
      INCLUDE 'LDGSIM.I'
      INCLUDE 'PARAME.I'
C
      REAL AZSBC(13,3,5), AZSBT(13,3,5), AZSBP(13,3,5)
C
      EQUIVALENCE
     *  (ASEBZML(1,1,1,1),AZSBC(1,1,1)),
     *  (ASEBZML(1,1,1,2),AZSBT(1,1,1)),
     *  (ASEBZML(1,1,1,3),AZSBP(1,1,1))
C
C     /IM240P/:
      DATA DSIZE/350*0.0/
C     /IMPAR2/:
      DATA ILDT,ITEST /105*0/
C     /IMPAR1/:
      DATA ICYIM,MODYR/28*NEYR/
      DATA ISTRIN     /7*0/
      DATA WAVPRE81   /7*0.05/
      DATA WAV81PLUS  /7*0.05/
      DATA CRIM       /7*85.0/
      DATA EXEMPTAGE  /7*25/
      DATA IMGRACE    /7*1/
      DATA EVFBEN     /28*1.0/
      DATA OLDFRAC    /0.75/
      DATA NEWFRAC    /0.25/
C     /LDGIMC/:
      DATA IMPEFF     /21*1.0/
      DATA WAVRDC     /7*0.80/
C     /LOGOBD/
      DATA WMIL       /0.85/
      DATA OBDIM      /0.99/
      DATA OBDXM      /3.6,8.0/
      DATA OBDXF      /0.9,0.1,0.0/
      DATA OBDXFLG    /.FALSE./
      DATA OBDXLIN    /.FALSE./
      DATA OBDVA      /3,7/
      DATA OBDVF      /0.9,0.1,0.0/
      DATA OBDVFLG    /.FALSE./
      DATA OBDVLIN    /.FALSE./
      DATA OBDIMFLG   /.FALSE./
      DATA OBDREP     /1.5/
C     /LDGSIM/:
      DATA RSDCUT     /7*1.0/
C
C High emitter identification rates IDR for idle tests
C
      DATA IDRID1/
C      RUNN.  RUNN.  RUNN.  START  START  START
C      CARB   TBI    PFI    CARB   TBI    PFI
     H 0.540, 0.580, 0.580, 0.250, 0.350, 0.350,
     C 0.540, 0.580, 0.580, 0.230, 0.300, 0.300,
     N 0.000, 0.000, 0.000, 0.000, 0.000, 0.000/
C
C High emitter identification rates IDR for idle/2500 tests
C
      DATA IDRID2/
C      RUNN.  RUNN.  RUNN.  START  START  START
C      CARB   TBI    PFI    CARB   TBI    PFI
     H 0.700, 0.600, 0.600, 0.300, 0.360, 0.360,
     C 0.650, 0.600, 0.600, 0.270, 0.320, 0.320,
     N 0.000, 0.000, 0.000, 0.000, 0.000, 0.000/
C
C  Hardcoded zero mile exhaust levels for starts emissions after
C  repair.
C
C  Carbureted vehicles 
C
        DATA AZSBC/
C  LDGV
     H 5*2.70,   8*3.11,
     C 5*28.33,  8*30.05,
     N 5*1.84,   8*1.49,
C  LDGT1
     H 5*2.70,   8*3.11,
     C 5*28.33,  8*30.05,
     N 5*1.84,   8*1.49,
C  LDGT2
     H 5*2.70,   8*3.11,
     C 5*28.33,  8*30.05,
     N 5*1.84,   8*1.49,
C  LDGT3
     H 5*2.70,   8*3.11,
     C 5*28.33,  8*30.05,
     N 5*1.84,   8*1.49,
C  LDGT4
     H 5*2.70,   8*3.11,
     C 5*28.33,  8*30.05,
     N 5*1.84,   8*1.49/
C
C  Ported Fuel-Injected vehicles
C
        DATA AZSBP/
C  LDGV
     H 5*2.70,   4*3.11,   4*2.60,
     C 5*28.33,  4*30.05,  4*19.0,
     N 5*1.84,   4*1.49,   4*1.48,
C  LDGT1
     H 5*2.70,   4*3.11,   4*2.60,
     C 5*28.33,  4*30.05,  4*19.0,
     N 5*1.84,   4*1.49,   4*1.48,
C  LDGT2
     H 5*2.70,   4*3.11,   4*2.60,
     C 5*28.33,  4*30.05,  4*19.0,
     N 5*1.84,   4*1.49,   4*1.48,
C  LDGT3
     H 5*2.70,   4*3.11,   4*2.60,
     C 5*28.33,  4*30.05,  4*19.0,
     N 5*1.84,   4*1.49,   4*1.48,
C  LDGT4
     H 5*2.70,   4*3.11,   4*2.60,
     C 5*28.33,  4*30.05,  4*19.0,
     N 5*1.84,   4*1.49,   4*1.48/
C
C  Throttle-Body Fuel-Injected vehicles
C
        DATA AZSBT/
C  LDGV
     H 5*2.70,   4*3.11,   4*2.60,
     C 5*28.33,  4*30.05,  4*19.0,
     N 5*1.84,   4*1.49,   4*1.48,
C  LDGT1
     H 5*2.70,   4*3.11,   4*2.60,
     C 5*28.33,  4*30.05,  4*19.0,
     N 5*1.84,   4*1.49,   4*1.48,
C  LDGT2
     H 5*2.70,   4*3.11,   4*2.60,
     C 5*28.33,  4*30.05,  4*19.0,
     N 5*1.84,   4*1.49,   4*1.48,
C  LDGT3
     H 5*2.70,   4*3.11,   4*2.60,
     C 5*28.33,  4*30.05,  4*19.0,
     N 5*1.84,   4*1.49,   4*1.48,
C  LDGT4
     H 5*2.70,   4*3.11,   4*2.60,
     C 5*28.33,  4*30.05,  4*19.0,
     N 5*1.84,   4*1.49,   4*1.48/
C
C  Coefficients used in the calculation of higher emitter identification
C  rates for ASM tests.(See tables 7a,7b,7c,& 7d pages 59 & 60 of M6.IM.001)
C
      DATA ASMCP1/
C  1. ASM 2525 with phase in cutpoints (vehicle ages less than or equal to 10)
     H -0.0301,   0.002811, -9.764e-4, 0.7324,
     C  0.00171,  6.151e-4,  0.00390,  0.2676,
     N  0.0289,  -0.001775,  0.015844,-1.0368,
C  2. ASM 2525 with final cutpoints (vehicle ages less than or equal to 10)
     H  0.001759, 1.588e-4 ,-0.001383, 0.9655,
     C  0.007035,-1.826e-4 , 0.001893, 0.6641,
     N -0.05733,  0.002875 ,-0.03234,  4.1179,
C  3. ASM 5015 with phase in cutpoints (vehicle ages less than or equal to 10)
     H -0.03507,  0.003147,  0.002789, 0.4213,
     C  1.775e-5, 7.131e-4,  0.006458, 0.05215,
     N -0.07537,  0.003535,  5.805e-4, 0.9142,
C  4. ASM 5015 with final cutpoints (vehicle ages less than or equal to 10)
     H -0.00494,  5.766e-4 ,-0.002577, 1.0894,
     C  0.003682, 5.1103e-5, 0.001625, 0.6787,
     N -0.12004,  0.006165, -0.02997,  4.1908,
C  5. ASM 2 mode with phase in cutpoints (vehicle ages less than or equal to 10)
     H -0.03397,  0.003077, -2.039e-4, 0.6986,
     C -3.874e-4, 7.258e-4,  0.004660, 0.2228,
     N -0.3024,   0.01462,  -0.10688,  11.890,
C  6. ASM 2 mode with final cutpoints (vehicle ages less than or equal to 10)
     H -0.006005, 6.0291e-4,-0.001904, 1.0791,
     C  9.809e-4, 1.5345e-4, 0.002478, 0.6573,
     N -0.1461,   0.007589, -0.036311, 4.9515/
C
      DATA ASMCP2/
C  1. ASM 2525 with phase in cutpoints (vehicle ages greater than 10)
     H -0.01390,  3.1742e-4,-0.001254, 0.8387,
     C -0.00747,  1.698e-4,  0.00331,  0.4557,
     N -0.00118,  7.0546e-5, 0.00833, -0.2571,
C  2. ASM 2525 with final cutpoints (vehicle ages greater than 10)
     H -0.005209, 1.161e-4, -0.001458, 1.0459,
     C -0.002926, 6.517e-5,  0.001498, 0.7751,
     N -0.001853, 3.9017e-5,-0.006412, 1.5271,
C  3. ASM 5015 with phase in cutpoints (vehicle ages greater than 10)
     H -0.01281,  2.945e-4,  0.003960, 0.3707,
     C -0.007068, 1.6312e-4, 0.005987, 0.2188,
     N -0.005603, 1.7994e-4, 0.01033, -0.3290,
C  4. ASM 5015 with final cutpoints (vehicle ages greater than 10)
     H -0.005936, 1.3247e-4,-0.002217, 1.1102,
     C -0.003783, 8.5342e-5, 0.001490, 0.7627,
     N -0.004784, 1.1196e-4, 5.749e-4, 0.9120,
C  5. ASM 2 mode with phase in cutpoints (vehicle ages greater than 10)
     H -0.01242, 2.8523e-4,  6.932e-4, 0.6740,
     C -0.00714, 1.6475e-4,  0.004418, 0.3669,
     N -0.01342, 3.5207e-4, -0.02177,  2.8099,
C  6. ASM 2 mode with final cutpoints (vehicle ages greater than 10)
     H -0.004063, 9.1144e-5,-0.001374, 1.0614,
     C -0.002706, 6.0543e-5,  0.002150, 0.7326,
     N -0.005176, 1.1762e-4,-0.002785, 1.2899/
C
C  RSD Effectiveness: the ability of a remote sensing device to
C  to recognize a high emitter if one passes over it..
C
        DATA RSDEFF/
C
     H 0.570, 0.433, 0.387, 0.348, 0.319, 0.262, 0.217, 0.182,
     H 0.150, 0.109, 0.071, 0.060, 0.046, 0.039, 0.028,
     C 0.596, 0.499, 0.442, 0.396, 0.352, 0.278, 0.213, 0.178,
     C 0.133, 0.107, 0.072, 0.053, 0.044, 0.034, 0.017,
     N 0.283, 0.178, 0.122, 0.091, 0.059, 0.054, 0.042, 0.018,
     N 0.015, 0.009, 0.006, 0.003, 0.003, 0.003, 0.003/
C
        DATA RSDCSE/
C lenient I/M program
     H 0.02, 0.01, 0.01, 0.02,
     C 0.07, 0.00, 0.00, 0.07,
     N 0.23, 0.07, 0.11, 0.12,
C strict  I/M program
     H 0.09, 0.04, 0.05, 0.06,
     C 0.07, 0.01, 0.01, 0.05,
     N 0.28, 0.07, 0.12, 0.15/
C
C   Repair running/start conversion factors
C
      DATA R_S_CONV/
C running      Tier1   LEV
     H  0.0, 3*0.230, 7*0.230,
     C  0.0, 3*0.338, 7*0.338,
     N  0.0, 3*0.900, 7*0.900,
C start
     H  0.0,3*11.000, 7*11.000,
     C  0.0, 3*4.149, 7* 4.149,
     N  0.0, 3*1.370, 7* 1.370/
C
      DATA ASMRAT / 1350*0.0 /
C
C AGE:  1       2      3      4      5      6      7      8      9
C      10      11     12     13     14     15     16     17     18
C      19      20     21     22     23     24     25
      DATA BIENRATIO /
     H 0.0000,0.4966,0.5877,0.6900,0.7400,0.7773,0.8000,0.8356,0.8740,
     H 0.8914,0.9200,0.9393,0.9468,0.9532,0.9595,0.9648,0.9689,0.9729,
     H 0.9755,0.9776,0.9794,0.9810,0.9828,0.9844,0.9852,
     C 0.0000,0.4976,0.5991,0.7100,0.7600,0.8000,0.8300,0.8640,0.8943,
     C 0.9083,0.9300,0.9469,0.9530,0.9589,0.9632,0.9673,0.9709,0.9744,
     C 0.9769,0.9788,0.9813,0.9829,0.9836,0.9849,0.9864,
     N 0.0000,0.5167,0.6136,0.7000,0.7500,0.7804,0.8100,0.8372,0.8730,
     N 0.8966,0.9134,0.9246,0.9353,0.9439,0.9515,0.9568,0.9615,0.9670,
     N 0.9720,0.9741,0.9757,0.9781,0.9793,0.9815,0.9826/
C
      END
