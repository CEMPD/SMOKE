C Version @(#)$Id$ $Source$ $Date$ 

      BLOCK DATA BD20
C
C  BLOCK DATA Subprogram 20: /MYRSAV/,/VMXCOM/
C
C  BLOCK DATA 20 initializes fleet characterization data:
C  (1) annual mileage accumulation rates by relative (to evaluation year)
C      model year (JDX) and vehicle type (IV)
C  (2) July 1 registration distributions, also by JDX and IV
C  (3) vehicle miles travelled (vmt) mix, by IV
C  (4) registration (reg) mix, by IV
C  (5) vehicle count (VCOUNT) , by ICY and IV
C
C  Common block array subscripts:
C
C  AMAR(25,8)    -  AMAR   ( JDX, IV )
C  JULMYR(25,8)  -  JULMYR ( JDX, IV )
C  REGMIX(8)     -  REGMIX ( IV )
C  TFNORM(8)     -  TFNORM ( IV )
C  VCOUNT(39,8)  -  VCOUNT ( ICY, IV )
C  VMTMIX(8)     -  VMTMIX ( IV )
C
C  Common block dictionary:
C
C   Name   Type              Description
C  ------  ----  ------------------------------------------------
C  /MYRSAV/:
C
C  NEWCUM   I    cumulative model year mileage recalculation flag
C  AMAR     R    annual mileage accumulation rates
C  JULMYR   R    July 1 registration distributions
C
C  /VMXCOM/:
C
C  REGMIX   R    normalized registration mix (registration percentage by IV)
C  TFNORM   R    for each IV, the sum across its 20 model year window of
C                registrations times mileage accrual rate
C  VCOUNT   R    vehicle count by calender year (ICY) and vehicle type (IV)
C  VMTMIX   R    normalized vehicle miles traveled mix (vmt percentage by IV)
C
C  Local array subscripts:
C
C  V12MAR(25,2)  -  V12MAR (JDX, IV )
C  V34MAR(25,2)  -  V34MAR (JDX, IV )
C  V56MAR(25,2)  -  V56MAR (JDX, IV )
C  V78MAR(25,2)  -  V78MAR (JDX, IV )
C
C  Local array dictionary:
C
C   Name   Type                Descriptions
C  ------  ----  ------------------------------------------------------
C  V12MAR   R    annual mileage accumulation rates for LDGVs and LDGT1s
C  V34MAR   R    annual mileage accumulation rates for LDGT2s and HDGVs
C  V56MAR   R    annual mileage accumulation rates for LDDVs and LDDTs
C  V78MAR   R    annual mileage accumulation rates for HDDVs and MCs
C
C
C  Notes:
C
C  The model year mileage recalculation flag "NEWCUM" is initialized here.
C
C  HDDV MAR values are initialized as all zeros, since they are
C  calculated by subroutine HDDMAR or supplied by the user.
C  Block data 20 was modified for MOBILE5v13 to include new
C  mileage accumulation rates.
C
C
      REAL JULMYR
C
      COMMON /MYRSAV/ AMAR(25,8),JULMYR(25,8),NEWCUM
      COMMON /VMXCOM/ REGMIX(8),TFNORM(8),VMTMIX(8),VCOUNT(39,8)
C
      DIMENSION V12MAR(25,2),V34MAR(25,2),V56MAR(25,2),V78MAR(25,2)
      DIMENSION V12COU(39,2),V34COU(39,2),V56COU(39,2),V78COU(39,2)
C
      EQUIVALENCE
     *(AMAR(1,1),V12MAR(1,1)),
     *(AMAR(1,3),V34MAR(1,1)),
     *(AMAR(1,5),V56MAR(1,1)),
     *(AMAR(1,7),V78MAR(1,1)),
C
     *(VCOUNT(1,1),V12COU(1,1)),
     *(VCOUNT(1,3),V34COU(1,1)),
     *(VCOUNT(1,5),V56COU(1,1)),
     *(VCOUNT(1,7),V78COU(1,1))
C
C  Initialize NEWCUM.
C
      DATA NEWCUM /1/
C
C  Initialize REGMIX and TFNORM.
C
      DATA REGMIX,TFNORM /16*0.0/
C
C  The vmt mix below represents a typical distribution, but is NOT used during
C  execution, where the mix is either computed in TFCALX or supplied by user:
C
      DATA VMTMIX/
     *     .7824, .0825, .0468, .0419, .0021, .0004, .0354, .0085/
C
C  Annual average mileage accumulation rates (divided by 100000)
C
      DATA V12MAR/
C  LDGV
     *.14390, .13612, .12875, .12180, .11522,
     *.10899, .10310, .09751, .09225, .08726,
     *.08254, .07807, .07386, .06987, .06608,
     *.06251, .05913, .05594, .05291, .05005,
     *.04735, .04478, .04237, .04007, .03790,
C  LDGT1
     *.15442, .14508, .13631, .12807, .12032,
     *.11305, .10621, .09979, .09376, .08809,
     *.08276, .07776, .07306, .06864, .06449,
     *.06059, .05693, .05348, .05025, .04721,
     *.04436, .04168, .03916, .03679, .03456/
      DATA V34MAR/
C  LDGT2
     *.14779, .14259, .13758, .13275, .12809,
     *.12359, .11924, .11505, .11101, .10711,
     *.10335, .09972, .09621, .09283, .08957,
     *.08642, .08339, .08046, .07763, .07490,
     *.07227, .06973, .06728, .06492, .06264,
C  HDGV
     *.17251, .16185, .15185, .14246, .13365,
     *.12539, .11764, .11037, .10355, .09715,
     *.09114, .08551, .08022, .07526, .07061,
     *.06625, .06215, .05831, .05471, .05132,
     *.04815, .04517, .04238, .03976, .03730/
      DATA V56MAR/
C  LDDV
     *.17825, .16478, .15233, .14081, .13017,
     *.12033, .11124, .10283, .09506, .08788,
     *.08123, .07509, .06942, .06417, .05932,
     *.05484, .05069, .04686, .04332, .04005,
     *.03702, .03422, .03163, .02924, .02703,
C  LDDT
     *.21004, .19125, .17415, .15858, .14440,
     *.13149, .11973, .10902, .09927, .09040,
     *.08231, .07495, .06825, .06215, .05659,
     *.05153, .04692, .04272, .03890, .03543,
     *.03226, .02937, .02675, .02435, .02218/
      DATA V78MAR/
C  HDDV
C  Warning: There are no default HDDV MAR values.  Either the user supplies
C           them or they are computed by subroutine HDDMAR.  See BD26.
     * .000, .000, .000, .000, .000, .000, .000, .000, .000, .000,
     * .000, .000, .000, .000, .000, .000, .000, .000, .000, .000,
     * .000, .000, .000, .000, .000,
C  MC
     *.04786, .04475, .04164, .03853, .03543,
     *.03232, .02921, .02611, .02300, .01989,
     *.01678, .01368, .00000, .00000, .00000,
     *.00000, .00000, .00000, .00000, .00000,
     *.00000, .00000, .00000, .00000, .00000/
C
C  Vehicle registration distributions (enter unmodified data)
C
      DATA JULMYR/
C  LDGV
     *.049,.079,.083,.082,.084,.081,.077,.056,.050,.051,.050,.054,.047,
     *.037,.024,.019,.014,.015,.011,.008,.006,.005,.004,.003,.010,
C  LDGT1
     *.063,.084,.084,.084,.084,.069,.059,.044,.036,.031,.030,.053,.047,
     *.046,.036,.028,.017,.022,.017,.014,.009,.008,.008,.005,.025,
C  LDGT2
     *.054,.072,.072,.072,.072,.052,.050,.034,.054,.031,.028,.080,.084,
     *.049,.039,.030,.018,.023,.018,.015,.009,.008,.009,.006,.026,
C  HDGV
     *.023,.047,.047,.047,.047,.038,.033,.021,.026,.029,.034,.064,.054,
     *.058,.051,.038,.043,.041,.035,.029,.021,.022,.022,.014,.117,
C  LDDV
     *.049,.079,.083,.082,.084,.081,.077,.056,.050,.051,.050,.054,.047,
     *.037,.024,.019,.014,.015,.011,.008,.006,.005,.004,.003,.010,
C  LDDT
     *.063,.084,.084,.084,.084,.069,.059,.044,.036,.031,.030,.053,.047,
     *.046,.036,.028,.017,.022,.017,.014,.009,.008,.008,.005,.025,
C  HDDV
     *.034,.067,.067,.067,.067,.073,.061,.040,.041,.051,.053,.066,.055,
     *.057,.045,.019,.023,.028,.024,.016,.011,.009,.007,.005,.016,
C  MC
     *.144,.168,.135,.109,.088,.070,.056,.045,.036,.029,.023,.097,.000,
     *.000,.000,.000,.000,.000,.000,.000,.000,.000,.000,.000,.000/
C
C VCOUNT - Vehicle Count by IV class, scaled by / 10**5; Note LDGV & LDDV
C counts are the identical sum of their registrations and ditto with LDGT1
C and LDDT. Hence MOBILE4 needs GSFVCT to split the sums properly, for the
C other IV's GSFVCT is 1.0 and has no effect.
C
C Calender Years run from 1982 to 2020, 39 Calender Years total.
C 1982 Vcounts are used for all pre-1982 Vehicles Counts
C IV = 8 (MC) Motor Cycles is a constant 4.219.
C  Modifed Release MOBILE4.1 [ ST99:D0V.S-TEST   20-AUG-1991 ]
C
C  LDGV
      DATA V12COU/
     * 106.867, 108.960, 112.018, 114.662, 117.268, 119.849, 121.519,
     * 122.758, 124.658, 126.556, 128.535, 130.591, 132.720, 134.917,
     * 137.030, 139.203, 141.431, 143.712, 146.044, 148.353, 150.707,
     * 153.104, 155.542, 158.020, 160.496, 163.006, 165.554, 168.132,
     * 170.745, 173.389, 176.062, 178.765, 181.497, 184.257, 186.638,
     * 189.037, 191.454, 193.887, 196.338,
C
C  LDGT1
C
     * 20.785, 21.026, 21.930, 23.102, 24.574, 26.002, 27.460,
     * 29.264, 30.558, 32.164, 33.745, 35.315, 36.822, 38.280,
     * 39.672, 41.045, 42.412, 43.768, 45.109, 46.427, 47.746,
     * 49.067, 50.383, 51.701, 53.017, 54.301, 55.605, 56.923,
     * 58.254, 59.594, 60.947, 62.309, 63.684, 65.072, 66.381,
     * 67.701, 69.039, 70.389, 71.759/
C
C  LDGT2
      DATA V34COU/
C
     * 10.175, 10.444, 10.942, 11.207, 11.886, 12.381, 13.372,
     * 13.940, 14.073, 14.503, 14.925, 15.330, 15.775, 16.246,
     * 16.704, 17.165, 17.624, 18.082, 18.551, 19.008, 19.464,
     * 19.919, 20.380, 20.844, 21.298, 21.791, 22.272, 22.749,
     * 23.225, 23.704, 24.183, 24.666, 25.152, 25.641, 26.047,
     * 26.452, 26.857, 27.261, 27.664,
C
C  HDGV
C
     *  4.522,  4.978,  5.300,  5.806,  5.827,  6.091,  6.190,
     *  6.443,  6.499,  6.680,  6.853,  7.026,  7.207,  7.396,
     *  7.569,  7.743,  7.922,  8.112,  8.303,  8.496,  8.700,
     *  8.907,  9.002,  9.225,  9.456,  9.694,  9.941, 10.182,
     * 10.394, 10.653, 10.920, 11.190, 11.463, 11.739, 11.976,
     * 12.217, 12.457, 12.695, 12.934/
C
C  LDDV
      DATA V56COU/
     * 106.867, 108.960, 112.018, 114.662, 117.268, 119.849, 121.519,
     * 122.758, 124.658, 126.556, 128.535, 130.591, 132.720, 134.917,
     * 137.030, 139.203, 141.431, 143.712, 146.044, 148.353, 150.707,
     * 153.104, 155.542, 158.020, 160.496, 163.006, 165.554, 168.132,
     * 170.745, 173.389, 176.062, 178.765, 181.497, 184.257, 186.638,
     * 189.037, 191.454, 193.887, 196.338,
C
C  LDDT
C
     * 20.785, 21.026, 21.930, 23.102, 24.574, 26.002, 27.460,
     * 29.264, 30.558, 32.164, 33.745, 35.315, 36.822, 38.280,
     * 39.672, 41.045, 42.412, 43.768, 45.109, 46.427, 47.746,
     * 49.067, 50.383, 51.701, 53.017, 54.301, 55.605, 56.923,
     * 58.254, 59.594, 60.947, 62.309, 63.684, 65.072, 66.381,
     * 67.701, 69.039, 70.389, 71.759/
C
C  HDDV
C
      DATA V78COU/
     *  1.505,  1.696,  1.972,  2.275,  2.541,  2.870,  3.202,
     *  3.554,  3.807,  4.144,  4.496,  4.853,  5.211,  5.566,
     *  5.920,  6.274,  6.630,  6.979,  7.329,  7.665,  7.996,
     *  8.328,  8.668,  8.988,  9.296,  9.601,  9.897, 10.202,
     * 10.494, 10.783, 11.067, 11.349, 11.629, 11.909, 12.148,
     * 12.382, 12.616, 12.847, 13.078,
C
C  MC
C
     * 39*4.219 /
C
      END
