      BLOCK DATA BD40 
C 
C   BLOCK DATA 40 contains coefficients used in determining the start 
C   and running basic emission factors. 
C  
C  Changes: (Last change first) 
C 
C  21 Sep 00 @EPA-jwk Bug235, Changed TLEV stds for LDGV/LDGT in CRT050K array. 
C  20 Sep 00 @EPA-djb Bug262, Changed diesel stds in CRT050K array. 
C  31 May 00 @EPA-ddj Bug114, Changed LDT1/2 HC Starts from 6*0.0964 to 
C            6*0.00964. 
C  25 May 00 @EPA-JK  Fixed Tier2 values in CRT050K for light duty gas 
C                     vehicles (cars and trucks). 
C  12 May 00 @EPA-djb Set T2 values in CRT050K indentical to LEV values. 
C  10 May 00 @epa-RAG    BUG#52 revised zero mile values and slope of 
C                        revised emission line 
C   2 Feb 00 @DynTel-ddj 1-032  Changed LDDT to LDDT12 and created 
C            LDDT34 by copying LDDT12 data to place holder for 
C            vehicle type 28. 
C    8 Mar 00 @EPA-eg            hrloop revisions 
C    2 Nov 99 @EPA-ELG, replaced values in arrays LDGVRC and LDGTRC. 
C   16 Jun 99 @DynTel-ddj 2-700, Fixed data entry error for LDGT1 CO NDSBP.   
C             [Start emissions deterioration rate/ slopes of start emissions 
C              curves as a function of mileage.] 
C   21 May 99 @DynTel-ddj 2-700, Modified CRIT2 in arrays LDGVRC and 
C             LDGTRC and made corrections to typos in NDSBC, NDSBT, 
C             NZRBP, NDRBC, NDRBP and NDRBT. 
C    4 May 99 @DynTel-ddj 2-694, Updated CRT050K array to reflect values 
C             provided in M6.EXH.007 and M6.EXH.009.  Removed CRT100K. 
C   10 Mar 99 @DynTel-RJD 2-694  Corrected data entry error in NSEBDR. 
C   01 Feb 99 @DynTel-RJD 2-694  Introduced certification standards for 
C             1994 and later model year LDGV/T vehicles. 
C   21 Sep 98 @DynTel-RJD 2-663  Moved the after repairs data to BD41. 
C   14 Sep 98 @DynTel-RJD 0-632  Revised data for cars and introduced 
C             data for trucks. 
C   27 Aug 98 @DynTel-RJD 0-632  Filled in blanks in basic emission 
C             factor arrays. Removed the MINNHY and MAXNHY arrays since 
C             the elements now all have the same value. 
C   11 Aug 98 @DynTel-RJD 2-663  Changed IT variable to IFDS. 
C 
C   Common block array subscripts: 
C 
C    /LDGBER/: 
C 
C     NSEBZML(13,3,5,3) - NSEBZML(IGY,IP,IV,IFDS) 
C     NSEBDR(13,3,5,3)  - NSEBDR(IGY,IP,IV,IFDS) 
C     HSEBZML(13,3,5,3) - HSEBZML(IGY,IP,IV,IFDS) 
C     NREBZML(13,3,5,3) - NREBZML(IGY,IP,IV,IFDS) 
C     NREBDR(13,3,5,3)  - NREBDR(IGY,IP,IV,IFDS) 
C     HREBZML(13,3,5,3) - HREBZML(IGY,IP,IV,IFDS) 
C     LDGRC(6,7,3,2)    - LDGRC(ICOEF,ITG,IP,IV) 
C     ITLDG(13,3,2)     - ITLDG(IGY,IFDS,IV) 
C     CRT050K(11,3,5)   - CRT050K(ISTD,IP,IV) 
C     CRT100K(11,3,5)   - CRT100K(ISTD,IP,IV) 
C      
C   Common block array dictionary: 
C 
C   Name   Type              Description 
C  ------  ----  --------------------------------------------- 
C 
C /LDGBER/: 
C 
C  NSEBZML  R     start zero mile level for normal emitters       
C  NSEBDR   R     start deterioration rate for normal emitters 
C  HSEBZML  R     start zero mile level for high emitters       
C  NREBZML  R     running zero mile level for normal emitters 
C  NREBDR   R     running deterioration rate for normal emitters 
C  HREBZML  R     running zero mile level for high emitters       
C  LDGRC    R     regression coefficients for average running emissions 
C  ITLDG    I     technology group number for 81-93 LDGV/Ts 
C  CRT050K  R     certification standards for 1994+ LDGV/Ts under 
C                 50000 miles. 
C  CRT100K  R     certification standards for 1994+ LDGV/Ts between 
C                 50000 and 100000 (120000 for trucks 3/4) miles. 
C 
C  Local array subscripts: 
C 
C  JKLn(13,3,5)  -  JKLn ( IGY, IP, IV ) 
C 
C  where J = N (normal), H (high) 
C        K = Z (zero mile level),  D (deterioration rate) 
C        L = S (start),  R (running) 
C        n = technology group number for fuel delivery system 
C 
C 
C  Local array dictionary: 
C 
C   Name      Type                Descriptions 
C  ---------  ----  ---------------------------------------------------- 
C  Jnm         R     as defined above 
C 
C  Notes: 
C 
C 01 Jul 98 @DynTel-RJD 2-660 Created BD40 to initialize coefficients 
C            for the computation of separate basic start and running 
C            emissions for light duty gasolene vehicles and trucks. 
C 02 Jul 98 @DynTel-RJD 2-663 Included after repair data for IM credit 
C            calculations. 
C 23 Jul 98 @DynTel-MLA 2-000  Changed LDGT names to conform to 
C           certification naming conventions. 
C 
C 
      IMPLICIT NONE 
      INCLUDE 'LDGBER.I' 
C 
      INTEGER ITLDGV(13,3), ITLDGT(13,3) 
      REAL NZSBC(13,3,5), NZSBT(13,3,5), NZSBP(13,3,5) 
      REAL NDSBC(13,3,5), NDSBT(13,3,5), NDSBP(13,3,5) 
      REAL HZSBC(13,3,5), HZSBT(13,3,5), HZSBP(13,3,5) 
      REAL NZRBC(13,3,5), NZRBT(13,3,5), NZRBP(13,3,5) 
      REAL NDRBC(13,3,5), NDRBT(13,3,5), NDRBP(13,3,5) 
      REAL HZRBC(13,3,5), HZRBT(13,3,5), HZRBP(13,3,5) 
      REAL LDGVRC(6,7,3), LDGTRC(6,7,3) 
C 
      EQUIVALENCE 
     *  (NSEBZML(1,1,1,1),NZSBC(1,1,1)), 
     *  (NSEBZML(1,1,1,2),NZSBT(1,1,1)), 
     *  (NSEBZML(1,1,1,3),NZSBP(1,1,1)) 
C 
      EQUIVALENCE 
     *  (NSEBDR(1,1,1,1), NDSBC(1,1,1)), 
     *  (NSEBDR(1,1,1,2), NDSBT(1,1,1)), 
     *  (NSEBDR(1,1,1,3), NDSBP(1,1,1)) 
C 
      EQUIVALENCE 
     *  (HSEBZML(1,1,1,1),HZSBC(1,1,1)), 
     *  (HSEBZML(1,1,1,2),HZSBT(1,1,1)), 
     *  (HSEBZML(1,1,1,3),HZSBP(1,1,1)) 
C 
      EQUIVALENCE 
     *  (NREBZML(1,1,1,1),NZRBC(1,1,1)), 
     *  (NREBZML(1,1,1,2),NZRBT(1,1,1)), 
     *  (NREBZML(1,1,1,3),NZRBP(1,1,1)) 
C 
      EQUIVALENCE 
     *  (NREBDR(1,1,1,1), NDRBC(1,1,1)), 
     *  (NREBDR(1,1,1,2), NDRBT(1,1,1)), 
     *  (NREBDR(1,1,1,3), NDRBP(1,1,1)) 
C 
      EQUIVALENCE 
     *  (HREBZML(1,1,1,1),HZRBC(1,1,1)), 
     *  (HREBZML(1,1,1,2),HZRBT(1,1,1)), 
     *  (HREBZML(1,1,1,3),HZRBP(1,1,1)) 
C 
      EQUIVALENCE 
     *  (LDGRC(1,1,1,1),LDGVRC(1,1,1)), 
     *  (LDGRC(1,1,1,2),LDGTRC(1,1,1)) 
 
 
      EQUIVALENCE 
     *  (ITLDG(1,1,1),ITLDGV(1,1)), 
     *  (ITLDG(1,1,2),ITLDGT(1,1)) 
C 
C  /LDGBER/: 
C 
C  Hardcoded zero mile exhaust levels for starts emissions in grams/start 
C 
C  Carbureted vehicles  
C 
        DATA NZSBC/ 
C  LDGV 
     H 2*2.121,  3*1.589,  8*1.493, 
     C 2*28.637, 3*24.442, 8*24.698, 
     N 2*1.6018, 3*0.748,  8*1.405, 
C  LDGT1 
     H 3*6.817,  10*3.916, 
     C 3*98.432, 10*78.286, 
     N 3*1.082,  10*0.143, 
C  LDGT2 
     H 3*6.817,  10*3.916, 
     C 3*98.432, 10*78.286, 
     N 3*1.082,  10*0.143, 
C  LDGT3 
     H 3*6.817,  10*3.916, 
     C 3*98.432, 10*78.286, 
     N 3*1.082,  10*0.143, 
C  LDGT4 
     H 3*6.817,  10*3.916, 
     C 3*98.432, 10*78.286, 
     N 3*1.082,  10*0.143/ 
C 
C  Ported Fuel-Injected vehicles 
C 
        DATA NZSBP/ 
C  LDGV 
     H 2*2.354,  5*2.3589, 6*1.9987, 
     C 2*20.038, 5*19.949, 6*18.972, 
     N 2*1.53,   5*1.461,  6*1.444, 
C  LDGT1 
     H 7*2.599,  6*2.873, 
     C 7*23.497, 6*32.178, 
     N 7*1.384,  6*1.597, 
C  LDGT2 
     H 7*2.599,  6*2.873, 
     C 7*23.497, 6*32.178, 
     N 7*1.384,  6*1.597, 
C  LDGT3 
     H 7*2.599,  6*2.873, 
     C 7*23.497, 6*32.178, 
     N 7*1.384,  6*1.597, 
C  LDGT4 
     H 7*2.599,  6*2.873, 
     C 7*23.497, 6*32.178, 
     N 7*1.384,  6*1.597/ 
C 
C  Throttle-Body Fuel-Injected vehicles 
C 
        DATA NZSBT/ 
C  LDGV 
     H 2*2.354,  5*2.3589, 6*1.9019, 
     C 2*20.038, 5*19.949, 6*19.233, 
     N 2*1.53,   5*1.461,  6*2.3, 
C  LDGT1 
     H 7*2.599,  6*2.599, 
     C 7*23.497, 6*42.456, 
     N 7*1.384,  6*4.294, 
C  LDGT2 
     H 7*2.599,  6*2.599, 
     C 7*23.497, 6*42.456, 
     N 7*1.384,  6*4.294, 
C  LDGT3 
     H 7*2.599,  6*2.599, 
     C 7*23.497, 6*42.456, 
     N 7*1.384,  6*4.294, 
C  LDGT4 
     H 7*2.599,  6*2.599, 
     C 7*23.497, 6*42.456, 
     N 7*1.384,  6*4.294/ 
C 
C  Hardcoded deterioration rates (emissions increase with the 
C  total miles travelled by a particular vehicle class) for normal 
C  starts emissions (units of grams/mile) 
C     basic_emissions(grams/start) = deter_rate(grams/(mile*start))*miles + zero_mile_exhaust(grams/start) 
C 
C  Carbureted vehicles 
C 
        DATA NDSBC/ 
C  LDGV 
     H 2*0.01361, 3*0.009408, 8*0.018238, 
     C 2*0.22673, 3*0.10577,  8*0.10947, 
     N 2*0.0,     3*0.00524,  8*0.0, 
C  LDGT1 
     H 3*0.00154, 10*0.00854, 
     C 3*0.324,   10*0.2564, 
     N 3*0.0,     10*0.00436, 
C  LDGT2 
     H 3*0.00154, 10*0.00854, 
     C 3*0.324,   10*0.2564, 
     N 3*0.0,     10*0.00436, 
C  LDGT3 
     H 3*0.00154, 10*0.00854, 
     C 3*0.324,   10*0.2564, 
     N 3*0.0,     10*0.00436, 
C  LDGT4 
     H 3*0.00154, 10*0.00854, 
     C 3*0.324,   10*0.2564, 
     N 3*0.0,     10*0.00436/ 
C 
C  Ported Fuel-Injected vehicles 
C 
        DATA NDSBP/ 
C  LDGV 
     H 2*0.00853, 5*0.00139, 6*0.00683, 
     C 2*0.22673, 5*0.0,     6*0.00703, 
     N 2*0.00059, 5*0.00141, 6*0.0022, 
C  LDGT1 
     H 7*0.00964, 6*0.0, 
     C 7*0.0613,  6*0.0168, 
     N 13*0.0, 
C  LDGT2 
     H 7*0.00964, 6*0.0, 
     C 7*0.0613,  6*0.0168, 
     N 13*0.0, 
C  LDGT3 
     H 7*0.00964, 6*0.0, 
     C 7*0.0613,  6*0.0168, 
     N 13*0.0, 
C  LDGT4 
     H 7*0.00964, 6*0.0, 
     C 7*0.0613,  6*0.0168, 
     N 13*0.0/ 
C 
C  Throttle-Body Fuel-Injected vehicles 
C 
        DATA NDSBT/ 
C  LDGV 
     H 2*0.00853, 5*0.00139, 6*0.002679, 
     C 2*0.22673, 5*0.0,     6*0.0, 
     N 2*0.00059, 5*0.00141, 6*0.0, 
C  LDGT1 
     H 7*0.00964, 6*0.00964, 
     C 7*0.0613,  6*0.1411, 
     N 7*0.0,     6*0.00324, 
C  LDGT2 
     H 7*0.00964, 6*0.00964, 
     C 7*0.0613,  6*0.1411, 
     N 7*0.0,     6*0.00324, 
C  LDGT3 
     H 7*0.00964, 6*0.01309, 
     C 7*0.0613,  6*0.1411, 
     N 7*0.0,     6*0.00324, 
C  LDGT4 
     H 7*0.00964, 6*0.01309, 
     C 7*0.0613,  6*0.1411, 
     N 7*0.0,     6*0.00324/ 
C 
C  Hardcoded zero mile exhaust levels for high starts emissions 
C 
C  Carbureted vehicles  
C 
        DATA HZSBC/ 
C  LDGV 
     H 13*10.52, 
     C 13*92.82, 
     N 7*0.748,6*1.405, 
C  LDGT1 
     H 3*17.865,  10*9.406, 
     C 3*179.549, 10*162.115, 
     N 3*1.082,10*0.143, 
C  LDGT2 
     H 3*17.865,  10*9.406, 
     C 3*179.549, 10*162.115, 
     N 3*1.082,10*0.143, 
C  LDGT3 
     H 3*17.865,  10*9.406, 
     C 3*179.549, 10*162.115, 
     N 3*1.082,10*0.143, 
C  LDGT4 
     H 3*17.865,  10*9.406, 
     C 3*179.549, 10*162.115, 
     N 3*1.082,10*0.143/ 
C 
C  Ported Fuel-Injected vehicles 
C 
      DATA HZSBP/ 
C  LDGV 
     H 7*5.313,  6*4.829, 
     C 2*92.82,  5*65.31, 6*38.06, 
     N 7*0.0,6*1.444, 
C  LDGT1 
     H 7*5.826,  6*5.212, 
     C 7*60.319, 6*83.862, 
     N 7*0.0,6*1.597, 
C  LDGT2 
     H 7*5.826,  6*5.212, 
     C 7*60.319, 6*83.862, 
     N 7*0.0,6*1.597, 
C  LDGT3 
     H 7*5.826,  6*5.212, 
     C 7*60.319, 6*83.862, 
     N 7*0.0,6*1.597, 
C  LDGT4 
     H 7*5.826,  6*5.212, 
     C 7*60.319, 6*83.862, 
     N 7*0.0,6*1.597/ 
C 
C  Throttle-Body Fuel-Injected vehicles 
C 
        DATA HZSBT/ 
C  LDGV 
     H 7*5.313,  6*3.293, 
     C 2*92.82,  5*65.31, 6*27.16, 
     N 2*1.53,5*1.461,6*2.3, 
C  LDGT1 
     H 7*5.826,  6*5.212, 
     C 7*60.319, 6*83.862, 
     N 7*1.384,6*4.294, 
C  LDGT2 
     H 7*5.826,  6*5.212, 
     C 7*60.319, 6*83.862, 
     N 7*1.384,6*4.294, 
C  LDGT3 
     H 7*5.826,  6*5.212, 
     C 7*60.319, 6*83.862, 
     N 7*1.384,6*4.294, 
C  LDGT4 
     H 7*5.826,  6*5.212, 
     C 7*60.319, 6*83.862, 
     N 7*1.384,6*4.294/ 
C 
C  Hardcoded zero mile exhaust levels for normal running emissions 
C 
C  Carbureted vehicles  
C 
        DATA NZRBC/ 
C  LDGV 
     H 2*0.1539, 3*0.1266, 8*0.0774, 
     C 2*1.3932, 3*0.7276, 8*0.5666, 
     N 2*0.5834, 3*0.5555, 8*0.4959, 
C  LDGT1 
     H 3*0.4918, 10*0.2684, 
     C 3*7.4202, 10*1.3553, 
     N 3*1.6159, 10*1.2872, 
C  LDGT2 
     H 3*0.4918, 10*0.2684, 
     C 3*7.4202, 10*1.3553, 
     N 3*1.6159, 10*1.2872, 
C  LDGT3 
     H 3*0.4918, 10*0.2684, 
     C 3*7.4202, 10*1.3553, 
     N 3*1.6159, 10*1.2872, 
C  LDGT4 
     H 3*0.4918, 10*0.2684, 
     C 3*7.4202, 10*1.3553, 
     N 3*1.6159, 10*1.2872/ 
C 
C  Ported Fuel-Injected vehicles 
C 
        DATA NZRBP/ 
C  LDGV 
     H 2*0.097,  5*0.0942, 6*0.0214, 
     C 2*1.5762, 5*1.4448, 6*0.4588, 
     N 2*0.4597, 5*0.4798, 6*0.2006, 
C  LDGT1 
     H 7*0.1338, 6*0.02989, 
     C 7*1.6222, 6*0.4927, 
     N 7*0.3150, 6*0.3024, 
C  LDGT2 
     H 7*0.1338, 6*0.02989, 
     C 7*1.6222, 6*0.4927, 
     N 7*0.3150, 6*0.3024, 
C  LDGT3 
     H 7*0.1338, 6*0.02989, 
     C 7*1.6222, 6*0.4927, 
     N 7*0.3150, 6*0.3024, 
C  LDGT4 
     H 7*0.1338, 6*0.02989, 
     C 7*1.6222, 6*0.4927, 
     N 7*0.3150, 6*0.3024/ 
C 
C  Throttle-Body Fuel-Injected vehicles 
C 
        DATA NZRBT/ 
C  LDGV 
     H 2*0.0970, 5*0.0942, 6*0.0042, 
     C 2*1.5762, 5*1.4448, 6*0.0, 
     N 2*0.4597, 5*0.4798, 6*0.2253, 
C  LDGT1 
     H 7*0.1338, 6*0.04664, 
     C 7*1.6222, 6*0.7663, 
     N 7*0.3150, 6*0.3150, 
C  LDGT2 
     H 7*0.1338, 6*0.04664, 
     C 7*1.6222, 6*0.7663, 
     N 7*0.3150, 6*0.3150, 
C  LDGT3 
     H 7*0.1338, 6*0.04664, 
     C 7*1.6222, 6*0.7663, 
     N 7*0.3150, 6*0.3150, 
C  LDGT4 
     H 7*0.1338, 6*0.04664, 
     C 7*1.6222, 6*0.7663, 
     N 7*0.3150, 6*0.3150/ 
C 
C  Hardcoded deterioration rates for normal running emissions 
C 
C  Carbureted vehicles 
C 
        DATA NDRBC/ 
C  LDGV 
     H 2*0.00127, 3*0.00121,  8*8.12E-4, 
     C 2*0.01389, 3*0.01691,  8*0.01371, 
     N 2*0.00233, 3*0.00273,  8*0.00170, 
C  LDGT1 
     H 3*0.006485, 10*0.002701, 
     C 3*0.03293,  10*0.0666, 
     N 3*0.000025, 10*0.0001, 
C  LDGT2 
     H 3*0.006485, 10*0.002701, 
     C 3*0.03293,  10*0.0666, 
     N 3*0.000025, 10*0.0001, 
C  LDGT3 
     H 3*0.006485, 10*0.002701, 
     C 3*0.03293,  10*0.0666, 
     N 3*0.000025, 10*0.0001, 
C  LDGT4 
     H 3*0.006485, 10*0.002701, 
     C 3*0.03293,  10*0.0666, 
     N 3*0.000025, 10*0.0001/ 
C 
C  Ported Fuel-Injected vehicles 
C 
        DATA NDRBP/ 
C  LDGV 
     H 2*0.00225, 5*0.00144, 6*0.001385, 
     C 2*0.02150, 5*0.01959, 6*0.02293, 
     N 2*0.00633, 5*0.00188, 6*0.00376, 
C  LDGT1 
     H 7*0.00328,  6*0.002376, 
     C 7*0.04311,  6*0.02678, 
     N 7*0.003171, 6*0.003904, 
C  LDGT2 
     H 7*0.00328,  6*0.002376, 
     C 7*0.04311,  6*0.02678, 
     N 7*0.003171, 6*0.003904, 
C  LDGT3 
     H 7*0.00328,  6*0.002376, 
     C 7*0.04311,  6*0.02678, 
     N 7*0.003171, 6*0.003904, 
C  LDGT4 
     H 7*0.00328,  6*0.002376, 
     C 7*0.04311,  6*0.02678, 
     N 7*0.003171, 6*0.003904/ 
C 
C  Throttle-Body Fuel-Injected vehicles 
C 
        DATA NDRBT/ 
C  LDGV 
     H 2*0.00225, 5*0.00144, 6*1.70E-3, 
     C 2*0.02150, 5*0.01959, 6*0.01990, 
     N 2*0.00633, 5*0.00188, 6*0.00381, 
C  LDGT1 
     H 7*0.00328,  6*0.002998, 
     C 7*0.04311,  6*0.03442, 
     N 7*0.003171, 6*0.003171, 
C  LDGT2 
     H 7*0.00328,  6*0.002998, 
     C 7*0.04311,  6*0.03442, 
     N 7*0.003171, 6*0.003171, 
C  LDGT3 
     H 7*0.00328,  6*0.002998, 
     C 7*0.04311,  6*0.03442, 
     N 7*0.003171, 6*0.003171, 
C  LDGT4 
     H 7*0.00328,  6*0.002998, 
     C 7*0.04311,  6*0.03442, 
     N 7*0.003171, 6*0.003171/ 
C 
C  Hardcoded zero mile exhaust levels for high running emissions 
C 
C  Carbureted vehicles  
C 
        DATA HZRBC/ 
C  LDGV 
     H 2*2.372,  11*1.845, 
     C 2*37.933, 11*27.653, 
     N 2*2.951,  11*2.872, 
C  LDGT1 
     H 3*3.71,   10*2.012, 
     C 3*80.726, 10*39.415, 
     N 3*5.014,  10*4.988, 
C  LDGT2 
     H 3*3.71,   10*2.012, 
     C 3*80.726, 10*39.415, 
     N 3*5.014,  10*4.988, 
C  LDGT3 
     H 3*3.71,   10*2.012, 
     C 3*80.726, 10*39.415, 
     N 3*5.014,  10*4.988, 
C  LDGT4 
     H 3*3.71,   10*2.012, 
     C 3*80.726, 10*39.415, 
     N 3*5.014,  10*4.988/ 
C 
C  Ported Fuel-Injected vehicles 
C 
        DATA HZRBP/ 
C  LDGV 
     H 7*2.372,  6*1.74, 
     C 7*37.933, 6*36.106, 
     N 7*2.951,  6*2.846, 
C  LDGT1 
     H 7*2.446,  6*2.12, 
     C 7*43.87,  6*33.283, 
     N 7*2.846,  6*2.846, 
C  LDGT2 
     H 7*2.446,  6*2.12, 
     C 7*43.87,  6*33.283, 
     N 7*2.846,  6*2.846, 
C  LDGT3 
     H 7*2.446,  6*2.12, 
     C 7*43.87,  6*33.283, 
     N 7*2.846,  6*2.846, 
C  LDGT4 
     H 7*2.446,  6*2.12, 
     C 7*43.87,  6*33.283, 
     N 7*2.846,  6*2.846/ 
C 
C  Throttle-Body Fuel-Injected vehicles 
C 
        DATA HZRBT/ 
C  LDGV 
     H 7*2.372,  6*3.394, 
     C 7*37.933, 6*46.527, 
     N 7*2.951,  6*2.872, 
C  LDGT1 
     H 7*2.446,  6*3.241, 
     C 7*43.87,  6*33.283, 
     N 7*2.846,  6*2.846, 
C  LDGT2 
     H 7*2.446,  6*3.241, 
     C 7*43.87,  6*33.283, 
     N 7*2.846,  6*2.846, 
C  LDGT3 
     H 7*2.446,  6*3.241, 
     C 7*43.87,  6*33.283, 
     N 7*2.846,  6*2.846, 
C  LDGT4 
     H 7*2.446,  6*3.241, 
     C 7*43.87,  6*33.283, 
     N 7*2.846,  6*2.846/ 
C 
C   The data for LDGV's is divided into 7 technology groups: 
C 
C    Group1     1988-93     PFI 
C    Group2     1988-93     TBI 
C    Group3     1983-87     FI 
C    Group4     1986-93     Carb 
C    Group5     1983-85     Carb 
C    Group6     1981-82     FI 
C    Group7     1981-82     Carb 
C 
C Regression Coefficients for LDGV average running Emissions 
C 
      DATA LDGVRC/    
C HC 
C Group   RZML     DET1      DET2      DET3    CRIT1    CRIT2 
     1  0.05158, 0.001300, 0.003600, 0.003600, 20.03,   100.0, 
     2  0.08427, 0.001300, 0.003300, 0.003300, 34.39,   100.0, 
     3  0.14791, 0.0,      0.007800, 0.0059,   18.8900, 81.3761, 
     4  0.08151, 0.003900, 0.005800, 0.005800, 19.83,   100.0, 
     5  0.16909, 0.0003,   0.009800, 0.009800, 25.2400, 100.0, 
     6  0.124,   0.009400, 0.013200, 0.0131,   11.2865, 70.5546, 
     7  0.2108,  0.004800, 0.015800, 0.015800, 10.1788, 100.0, 
C CO 
     1  0.7983,  0.031,    0.0707,   0.0707,   13.78,       100.0, 
     2  2.5684,  0.031,    0.031,    0.031,    13.78,   100.0, 
     3  2.14158, 0.0,      0.1091,   0.0846,   19.0400, 69.7826, 
     4  0.69103, 0.07270,  0.1034,   0.1034,   21.1292, 100.0, 
     5  1.09825, 0.00000,  0.1537,   0.1537,   25.6800, 100.0, 
     6  1.72695, 0.1817,   0.2401,   0.2401,   16.6175,     100.0, 
     7  2.93608, 0.1414,   0.2908,   0.2873,   8.7946,  15.0178, 
C NOX 
     1  0.2583,  0.00103,  0.005793, 0.005793, 18.58,   100.0, 
     2  0.2931,  0.00103,  0.005771, 0.005771, 21.55,   100.0, 
     3  0.59755, 0.002281, 0.006432, 0.006432, 34.2524, 100.0, 
     4  0.55222, 0.002111, 0.004458, 0.004458, 26.124,  100.0, 
     5  0.56136, 0.000313, 0.006248, 0.006248, 12.5204, 100.0, 
     6  0.63699, 0.0,      0.006897, 0.006897, 30.6571, 100.0, 
     7  0.61213, 0.000300, 0.006593, 0.0063,   8.7946,  16.9951/ 
C 
C                        CARB             TBI             PFI 
        DATA ITLDGV/ 2*7, 3*5, 8*4,  2*6, 5*3, 6*2,  2*6, 5*3, 6*1/ 
C 
C   The data for LDGT's is divided into 5 technology groups: 
C     
C    Group1     1988-93   PFI 
C    Group2     1988-93   TBI 
C    Group3     1981-87   FI(both PFI and TBI) 
C    Group4     1984-93   CARB 
C    Group5     1981-83   CARB 
C 
C Regression Coefficients for LDGT average running Emissions 
C 
        DATA LDGTRC/    
C HC 
C Group RZML      DET1      DET2      DET3     CRIT1    CRIT2 
     1 0.09321, 0.001300, 0.003800, 0.003800, 23.4042, 100.0, 
     2 0.07834, 0.001300, 0.005600, 0.0055,   16.2429, 55.1613, 
     3 0.29272, 0.0,      0.009900, 0.009900, 40.5800, 100.0, 
     4 0.24951, 0.0,      0.008300, 0.008300, 36.0100, 100.0, 
     5 0.65869, 0.001800, 0.012700, 0.012700, 15.9948, 100.0, 
     6 0.00000, 0.000000, 0.000000, 0.000000, 0.00000, 0.00000, 
     7 0.00000, 0.000000, 0.000000, 0.000000, 0.00000, 0.00000, 
C CO 
     1 0.90173, 0.032600, 0.06830,  0.0623,   16.8003, 58.6841, 
     2 1.14387, 0.032600, 0.08170,  0.08170,  17.5399, 100.0, 
     3 5.23367, 0.054500, 0.11900,  0.11900,  55.0316, 100.0, 
     4 1.53835, 0.0,      0.13270,  0.13270,  28.9000, 100.0, 
     5 9.07038, 0.1040,   0.16750,  0.16750,  18.8569, 100.0, 
     6 0.00000, 0.00000,  0.00000,  0.00000,  0.00000, 0.00000, 
     7 0.00000, 0.00000,  0.00000,  0.00000,  0.00000, 0.00000, 
C NOx 
     1 0.37824, 0.000208, 0.0046,   0.0046,   21.204,  100.0, 
     2 0.33462, 0.000208, 0.0042,   0.0034,   16.243,  55.1613, 
     3 0.53883, 0.0,      0.0056,   0.0056,   32.211,  100.0, 
     4 1.32337, 0.0,      0.0001,   0.0001,   74.2057, 100.0, 
     5 1.66598, 0.00080,  0.00080,  0.00080,  100.0,   100.0, 
     6 0.00000, 0.000000, 0.000000, 0.000000, 0.00000, 0.00000, 

     7 0.00000, 0.000000, 0.000000, 0.000000, 0.00000, 0.00000/ 
C 
C                     CARB              TBI             PFI 
      DATA ITLDGT/  3*5, 10*4,       7*3, 6*2,        7*3, 6*1/ 
C 
C Certification standards for 1994+ LDGV/Ts under 50000 miles 
C 
      DATA CRT050K/ 
C  LDGV 
C      TR0  ITR1 TR1  TR2   ITLEV  TLEV   ILEV,  LEV    IULEV ULEV  ZEV 
     H 0.41,0.25,0.25,0.25, 0.125, 0.125, 0.075, 0.075, 0.04, 0.04, 0.0, 
     C 3.4, 3.4, 3.4, 3.4,  3.4,   3.4,   3.4,   3.4,   1.7,  1.7,  0.0, 
     N 1.0, 0.4, 0.4, 0.4,  0.4,   0.4,   0.2,   0.2,   0.2,  0.2,  0.0, 
C  LDGT1 
     H 0.41,0.25,0.25,0.25, 0.125, 0.125, 0.075, 0.075, 0.04, 0.04, 0.0, 
     C 3.4, 3.4, 3.4, 3.4,  3.4,   3.4,   3.4,   3.4,   1.7,  1.7,  0.0, 
     N 1.0, 0.4, 0.4, 0.4,  0.4,   0.4,   0.2,   0.2,   0.2,  0.2,  0.0, 
C  LDGT2 
     H 0.41,0.32,0.32,0.32, 0.16,  0.16,  0.1,   0.1,   0.05, 0.05, 0.0, 
     C 4.4, 4.4, 4.4, 4.4,  4.4,   4.4,   4.4,   4.4,   2.2,  2.2,  0.0, 
     N 0.7, 0.7, 0.7, 0.7,  0.7,   0.7,   0.4,   0.4,   0.4,  0.4,  0.0, 
C  LDGT3 
     H 0.41,0.32,0.32,0.32, 0.16,  0.16,  0.16,  0.16,  0.1,  0.1,  0.0, 
     C 4.4, 4.4, 4.4, 4.4,  4.4,   4.4,   4.4,   4.4,   2.2,  2.2,  0.0, 
     N 0.7, 0.7, 0.7, 0.7,  0.4,   0.4,   0.4,   0.4,   0.4,  0.4,  0.0, 
C  LDGT4 
     H 0.41,0.39,0.39,0.39, 0.195, 0.195, 0.195, 0.195, 0.117,0.117,0.0, 
     C 5.0, 5.0, 5.0, 5.0,  5.0,   5.0,   5.0,   5.0,   2.5,  2.5,  0.0, 
     N 1.1, 1.1, 1.1, 1.1,  0.6,   0.6,   0.6,   0.6,   0.6,  0.6,  0.0, 
C  LDDV (50k levels) 
     H 0.34,0.25,0.25,0.0, 0.125, 0.125, 0.075, 0.075, 0.04,  0.04, 0.0, 
     C 3.4, 3.4, 3.4, 0.0,  3.4,   3.4,   3.4,   3.4,   1.7,  1.7,  0.0, 
     N 1.0, 1.0, 1.0, 0.0,  0.4,   0.4,   0.2,   0.2,   0.2,  0.2,  0.0, 
C  LDDT12 (Full Life for LDT2 Only) 
     H 0.67,0.4, 0.4, 0.0,  0.2,   0.2,  0.13,  0.13,  0.07, 0.07,  0.0, 
     C10.0, 5.5, 5.5, 0.0,  5.5,   5.5,   5.5,   5.5,   2.8,  2.8,  0.0, 
     N 1.7,0.97,0.97, 0.0,  0.9,   0.9,   0.5,   0.5,   0.5,  0.5,  0.0, 
C  LDDT34 (Full Life for LDT4 Only) 
     H 0.67,0.56,0.56,0.0, 0.28,  0.28,  0.28,  0.28, 0.167,0.167,  0.0, 
     C10.0, 7.3, 7.3, 0.0,  7.3,   7.3,   7.3,   7.3,   3.7,  3.7,  0.0, 
     N 1.7,1.53,1.53, 0.0,  1.5,   1.5,   1.5,   1.5,   0.8,  0.8,  0.0/ 
C 
      END 
