      BLOCK DATA BD_AC 
C 
C     BLOCK BD_AC initializes the AC variables 
C 
C  Changes: (Last change first) 
C 
C   14 Dec 99 @DynTel ddf 1-013 New data block 
C 
C 
C  Common block dictionary: 
C 
C   Name                 Type              Description 
C  ------                ----  ---------------------------------------------- 
C  /AC_DATA/: 
C 
C  VHCLES_AC_EQPD(IV,MY)   R    fractions of LDGV & LDGT1-4 vehicles 
C                               equipped with a/c 
C  NTFCNNG_AC_UNITS(IV,MY) R    fractions of LDGV & LDGT1-4 vehicles 
C                               equipped with a/c which are not 
C                               functioning 
C  FULL_A(index)           R    first parameter for A/C correction 
C                               factors at Full Usage. 
C  FULL_B(index)                second parameter for A/C compressor 
C                               factors ar Full Usage. 
C  DMND(C,DLY_PRIOD)       R    C = 1, zero order parameter for A/C 
C                               compressor loading 
C                               dependence on the heat index 
C                               C = 2, first order parameter for A/C 
C                               compressor loading 
C                               dependence on the heat index 
C                               C = 3, second order parameter for A/C 
C                               compressor loading 
C                               dependence on the heat index 
C 
C    DEFAULT VALUES are 
C          SUNRISE      ==  1,    6:00 AM 
C          SUNSET       == 16,    9:00 PM 
C          PEAKSUN(1)   ==  5,   10:00 AM 
C          PEAKSUN(2)   == 11,    4:00 PM 
C 
C 
      IMPLICIT NONE 
C 
      INCLUDE 'AC.I' 
C 
      DATA SUNRISE, SUNSET, PEAKSUN/ 1, 16, 5, 11/ 
C 
      DATA CLD_CVR/ 0.0/ 
C 
      DATA FULL_A, FULL_B / 
     H            0.001162, 0.506, 
     C            0.815, 0.104, 0.154, 0.678, 0.119, 
     N            4.867, 1.930, 0.655, 
C     FULL_B 
     H            0.00, 0.00, 
     C            0.05272, 0.00, 0.00, 0.00, 0.00, 
     N            2.296, 0.769, 0.00 / 
C 
C Yearly fractions of vehicles manufactured with air conditioners 
C 
      DATA VHCLES_AC_EQPD/ 
C 197x  LDVs   LDTs 
     2  0.59,  0.29, 
     3  0.72,  0.29, 
     4  0.61,  0.29, 
     5  0.63,  0.29, 
     6  0.67,  0.31, 
     7  0.72,  0.35, 
     8  0.72,  0.38, 
     9  0.69,  0.36, 
     O  0.63,  0.35, 
     1  0.67,  0.39, 
     2  0.70,  0.45, 
     3  0.75,  0.47, 
     4  0.78,  0.52, 
     5  0.80,  0.53, 
     6  0.80,  0.55, 
     7  0.76,  0.60, 
     8  0.79,  0.65, 
     9  0.78,  0.73, 
     O  0.88,  0.76, 
     1  0.88,  0.78, 
     2  0.89,  0.82, 
     3  0.90,  0.84, 
     4  0.92,  0.85, 
     5  0.94,  0.89, 
     6  0.95,  0.91, 
     7  0.96,  0.93, 
     8  0.98,  0.95/ 
C 
C Yearly estimates of fractions of vehicles with air conditioners 
C that are NOT functioning. 
C 
      DATA NTFCNNG_AC_UNITS/ 
C  AGE   Freon  R-134a 
     1   0.000, 0.000, 
     2   0.000, 0.000, 
     3   0.000, 0.000, 
     4   0.018, 0.004, 
     5   0.022, 0.006, 
     6   0.026, 0.008, 
     7   0.030, 0.010, 
     8   0.034, 0.012, 
     9   0.050, 0.025, 
     O   0.050, 0.025, 
     1   0.050, 0.025, 
     2   0.050, 0.025, 
     3   0.050, 0.025, 
     4   0.060, 0.038, 
     5   0.060, 0.038, 
     6   0.060, 0.038, 
     7   0.060, 0.038, 
     8   0.060, 0.038, 
     9   0.069, 0.050, 
     O   0.069, 0.050/ 
C 
C A/C demand factor equation parameters ; 
C      indexing : MORNING/ AFTERNOON = 1; PEAK SUN = 2; 
C                 DAYTIME = 3; NIGHT = 4; ALL = 5 
C 
      DATA DMND/ 
C INDEX    constant       a           b 
     1    -2.930273,  0.059110,  -0.000213, 
     2    -5.307355,  0.113973,  -0.000521, 
     3    -4.101082,  0.086382,  -0.000367, 
     4    -1.257412,  0.006753,   0.000143, 
     5    -3.631541,  0.072465,  -0.000276/ 
C 
      END 
