      BLOCK DATA BD29
C
C  BLOCK DATA Subprogram 29: /HDCCOM/
C
C  BLOCK DATA 29 initializes grams-per-brakehorsepower-hour to grams-per-mile
C  units conversion variables.
C
C  Changes : (Last change first)
C  22 Oct 02 @EPA-EG      Changed 24 HDCFAC values from 0.000 to 1.792 for class HDCGAS8
C  26 Sep 00 @EPA-EG      Added initialization of HD conversion factor arrays.
C  07 march 2000 @DynTel-ZK 1-031 Entered new conversion factor from spread sheets. 
C  07 Oct 98 @DynTel-HQ 2-678 Comment statements
C  29 Sep 98 @DynTel-HQ 2-678 Conversion Factor
C
C  Common block array subscripts:
C
C  HDCFAC(56,19)  -  HDCFAC ( IG, IVHD )
C
C  Common block dictionary:
C
C   Name   Type              Description
C  ------  ----  -------------------------------------------
C  /HDCCOM/:
C
C  HDCFAC   R    HDV gpbhp to gpm conversion factors
C  LOWYR    I    lowest my for which HDCFAC is defined.  All my <= LOWYR are
C                assigned LOWYR's conversion factors.
C  MAXCIG   I    maximum HDCFAC index allowed = 1st dimension of HDCFAC = 56
C                = index of highest my for which HDCFAC is defined, currently
C                2021.  All my >= 2021 are assigned MAXCIG's factors.
C
C  Local array subscripts:
C
C  HDCGAS1(56)  -  HDCGAS ( IG )
C  HDCGAS2(56)  -  HDCGAS ( IG )
C  HDCGAS3(56)  -  HDCGAS ( IG )
C  HDCGAS4(56)  -  HDCGAS ( IG )
C  HDCGAS5(56)  -  HDCGAS ( IG )
C  HDCGAS6(56)  -  HDCGAS ( IG )
C  HDCGAS7(56)  -  HDCGAS ( IG )
C  HDCGAS8(56)  -  HDCGAS ( IG )
C  HDCDIE1(56)  -  HDCDIE ( IG )
C  HDCDIE2(56)  -  HDCDIE ( IG )
C  HDCDIE3(56)  -  HDCDIE ( IG )
C  HDCDIE4(56)  -  HDCDIE ( IG )
C  HDCDIE5(56)  -  HDCDIE ( IG )
C  HDCDIE6(56)  -  HDCDIE ( IG )
C  HDCDIE7(56)  -  HDCDIE ( IG )
C  HDCDIE8(56)  -  HDCDIE ( IG )
C  HDCGBS1(56)  -  HDCGBS ( IG )
C
C  Local array dictionary:
C
C   Name   Type                Description
C  ------  ----  -----------------------------------------
C  HDCGAS   R    HDGV gpbhp to gpm units conversion factors
C  HDCDIE   R    HDDV gpbhp to gpm units conversion factors
C  HDCGBS  R     GAS BUS gpbhp to gpm units conversion factors
C  HDCDBS  R     Diesel BUS gpbhp to gpm units conversion factors
C
C  Notes:
C
C 12 Jan 98 @DynTel-MLA 1-653  Modified the descriptions of LOWYR and
C           MAXCIG in the common block dictionary above to improve the
C           accuracy and clarity of the descriptions.
C
      IMPLICIT NONE
      INCLUDE 'HDCCOM.I'
C
      REAL HDCGAS1(56),HDCGAS2(56),HDCGAS3(56),HDCGAS4(56), 
     *     HDCGAS5(56),HDCGAS6(56),HDCGAS7(56),HDCGAS8(56),
     *     HDCDIE1(56),HDCDIE2(56),HDCDIE3(56),HDCDIE4(56),
     *     HDCDIE5(56),HDCDIE6(56),HDCDIE7(56),HDCDIE8(56),
     *     HDCGBS1(56),HDCDBS1(56),HDCDBS2(56)
C
      EQUIVALENCE
     *   (HDCFAC(1,1),HDCGAS1(1)),
     *   (HDCFAC(1,2),HDCGAS2(1)),
     *   (HDCFAC(1,3),HDCGAS3(1)),
     *   (HDCFAC(1,4),HDCGAS4(1)),
     *   (HDCFAC(1,5),HDCGAS5(1)),
     *   (HDCFAC(1,6),HDCGAS6(1)),
     *   (HDCFAC(1,7),HDCGAS7(1)),
     *   (HDCFAC(1,8),HDCGAS8(1)),
     *   (HDCFAC(1,9),HDCDIE1(1)),
     *   (HDCFAC(1,10),HDCDIE2(1)),
     *   (HDCFAC(1,11),HDCDIE3(1)),
     *   (HDCFAC(1,12),HDCDIE4(1)),
     *   (HDCFAC(1,13),HDCDIE5(1)),
     *   (HDCFAC(1,14),HDCDIE6(1)),
     *   (HDCFAC(1,15),HDCDIE7(1)),
     *   (HDCFAC(1,16),HDCDIE8(1)),
     *   (HDCFAC(1,17),HDCGBS1(1)),
     *   (HDCFAC(1,18),HDCDBS1(1)),
     *   (HDCFAC(1,19),HDCDBS2(1))
C
C  Gpbhp to gpm units conversion factors
C
C  For HDGV
C
      DATA HDCGAS1/
     *    0.960, 0.960, 0.960, 0.960, 0.960, 0.960, 0.960, 0.960, 0.938,
     *    0.917, 0.896, 0.877, 0.858, 0.841, 0.835, 0.819, 0.814, 0.813,
     *    0.813, 0.811, 0.809, 1.068, 1.071, 1.073, 1.076, 1.079, 1.082,
     *    1.085, 1.089, 1.092, 1.096, 1.096, 24*1.096/
      DATA HDCGAS2/
     *    1.353, 1.353, 1.353, 1.353, 1.353, 1.353, 1.353, 1.353, 1.349,
     *    1.345, 1.341, 1.336, 1.332, 1.328, 1.343, 1.339, 1.355, 1.354,
     *    1.353, 1.350, 1.346, 1.152, 1.151, 1.150, 1.150, 1.149, 1.149,
     *    1.149, 1.149, 1.149, 1.150, 1.150, 24*1.150/
      DATA HDCGAS3/
     *    1.353, 1.353, 1.353, 1.353, 1.353, 1.353, 1.353, 1.353, 1.349,
     *    1.345, 1.341, 1.336, 1.332, 1.328, 1.343, 1.339, 1.355, 1.354,
     *    1.353, 1.351, 1.348, 1.164, 1.160, 1.156, 1.152, 1.149, 1.146,
     *    1.143, 1.140, 1.137, 1.134, 1.134, 24*1.134/
      DATA HDCGAS4/
     *    1.353, 1.353, 1.353, 1.353, 1.353, 1.353, 1.353, 1.353, 1.349,
     *    1.345, 1.341, 1.336, 1.332, 1.328, 1.343, 1.339, 1.355, 1.353,
     *    1.351, 1.347, 1.342, 1.291, 1.294, 1.297, 1.301, 1.305, 1.308,
     *    1.312, 1.316, 1.32,  1.324, 1.324, 24*1.324/
      DATA HDCGAS5/
     *    1.513, 1.513, 1.513, 1.513, 1.513, 1.513, 1.513, 1.513, 1.540, 
     *    1.554, 1.551, 1.537, 1.518, 1.492, 1.452, 1.419, 1.358, 1.347,
     *    1.337, 1.327, 1.317, 1.311, 1.310, 1.310, 1.309, 1.309, 1.309,
     *    1.309, 1.309, 1.310, 1.311, 1.311, 24*1.311/
      DATA HDCGAS6/
     *    1.900, 1.900, 1.900, 1.900, 1.900, 1.900, 1.900, 1.900, 1.891,
     *    1.887, 1.879, 1.863, 1.813, 1.751, 1.706, 1.710, 1.714, 1.702,
     *    1.692, 1.680, 1.668, 1.383, 1.389, 1.395, 1.402, 1.409, 1.416,
     *    1.423, 1.430, 1.438, 1.446, 1.446, 24*1.446/
      DATA HDCGAS7/
     *    2.052, 2.052, 2.052, 2.052, 2.052, 2.052, 2.052, 2.052, 1.991,
     *    1.938, 1.973, 1.839, 1.790, 1.747, 1.731, 1.691, 1.674, 1.662,
     *    1.651, 1.639, 1.627, 1.503, 1.507, 1.510, 1.513, 1.517, 1.521,
     *    1.526, 1.530, 1.535, 1.540, 1.540, 24*1.540/
      DATA HDCGAS8/
     *    2.597, 2.597, 2.597, 2.597, 2.597, 2.597, 2.597, 2.597, 2.486,
     *    2.390, 2.296, 2.140, 2.083, 2.033, 2.014, 1.968, 1.948, 1.934,
     *    1.921, 1.907, 1.893, 1.749, 1.754, 1.757, 1.761, 1.765, 1.770,
     *    1.776, 1.780, 1.786, 1.792, 1.792, 24*1.792/
C
C  For HDDV
C
      DATA HDCDIE1/
     *    0.942, 0.942, 0.942, 0.942, 0.942, 0.942, 0.942, 0.942, 0.942,
     *    0.942, 0.942, 0.942, 0.942, 0.942, 0.942, 0.942, 0.923, 0.923,
     *    0.922, 0.921, 0.919, 0.919, 1.103, 1.101, 1.099, 1.097, 1.095,
     *    1.094, 1.093, 1.091, 1.090, 1.090, 24*1.090/
      DATA HDCDIE2/
     *    2.161, 2.161, 2.161, 2.161, 2.161, 2.161, 2.161, 2.161, 2.009,
     *    1.879, 1.762, 1.762, 1.762, 1.762, 1.762, 1.762, 1.762, 1.762,
     *    1.762, 1.762, 1.762, 1.762, 1.253, 1.252, 1.251, 1.250, 1.250,
     *    1.250, 1.250, 1.250, 1.250, 1.250, 24*1.250/
      DATA HDCDIE3/
     *    2.161, 2.161, 2.161, 2.161, 2.161, 2.161, 2.161, 2.161, 2.009,
     *    1.879, 1.762, 1.762, 1.762, 1.762, 1.762, 1.762, 1.762, 1.762,
     *    1.762, 1.762, 1.762, 1.762, 1.496, 1.491, 1.486, 1.481, 1.476,
     *    1.472, 1.467, 1.463, 1.458, 1.458, 24*1.458/
      DATA HDCDIE4/
     *    2.161, 2.161, 2.161, 2.161, 2.161, 2.161, 2.161, 2.161, 2.009,
     *    1.879, 1.762, 1.762, 1.762, 1.762, 1.762, 1.762, 1.762, 1.762,
     *    1.762, 1.762, 1.762, 1.762, 1.676, 1.662, 1.649, 1.636, 1.623,
     *    1.610, 1.597, 1.585, 1.573, 1.573, 24*1.573/
      DATA HDCDIE5/
     *    2.452, 2.452, 2.452, 2.502, 2.502, 2.554, 2.608, 2.608, 2.551,
     *    2.442, 2.342, 2.303, 2.216, 2.135, 2.060, 1.992, 1.927, 1.911,
     *    1.892, 1.877, 1.865, 1.865, 1.979, 1.974, 1.969, 1.964, 1.960,
     *    1.955, 1.951, 1.947, 1.942, 1.942, 24*1.942/
      DATA HDCDIE6/
     *    2.591, 2.645, 2.645, 2.645, 2.701, 2.701, 2.701, 2.760, 2.775,
     *    2.780, 2.760, 2.753, 2.616, 2.412, 2.254, 2.229, 2.193, 2.176,
     *    2.159, 2.143, 2.127, 2.127, 2.392, 2.394, 2.396, 2.398, 2.400,
     *    2.403, 2.405, 2.407, 2.409, 2.409, 24*2.409/
      DATA HDCDIE7/
     *    2.845, 2.904, 2.904, 2.966, 2.966, 3.031, 3.031, 3.031, 3.104,
     *    3.098, 3.156, 3.180, 3.126, 3.085, 3.062, 3.106, 3.089, 3.059,
     *    3.035, 3.010, 2.987, 2.987, 2.946, 2.922, 2.898, 2.874, 2.851,
     *    2.828, 2.806, 2.784, 2.763, 2.763, 24*2.763/
      DATA HDCDIE8/
     *    3.074, 3.138, 3.138, 3.205, 3.205, 3.275, 3.275, 3.275, 3.326,
     *    3.299, 3.353, 3.402, 3.361, 3.307, 3.332, 3.255, 3.152, 3.150,
     *    3.141, 3.138, 3.129, 3.129, 3.263, 3.231, 3.201, 3.171, 3.141,
     *    3.113, 3.085, 3.058, 3.031, 3.031, 24*3.031/
C
C     GAS BUS
C
C
      DATA HDCGBS1/
     *    1.242, 1.242, 1.245, 1.243, 1.243, 1.243, 1.180, 1.182, 1.182,
     *    1.190, 1.185, 1.185, 1.184, 1.185, 1.147, 1.155, 1.161, 1.161,
     *    1.161, 1.161, 1.161, 1.622, 1.628, 1.633, 1.639, 1.645, 1.651,
     *    1.658, 1.664, 1.671, 1.796, 1.796, 24*1.796/
C
C    DIESEL BUS
C
      DATA HDCDBS1/
     *    2.847, 2.847, 2.963, 3.021, 3.000, 3.000, 3.000, 3.000, 3.000,
     *    3.000, 3.000, 2.839, 2.872, 2.875, 2.912, 3.005, 3.066, 3.066,
     *    3.066, 3.066, 3.066, 3.066, 4.009, 4.009, 4.010, 4.010, 4.013,
     *    4.015, 4.017, 4.021, 4.025, 4.025, 24*4.025/
      DATA HDCDBS2/
     *    1.653, 1.653, 1.653, 1.653, 1.653, 1.653, 1.653, 1.653, 1.653,
     *    1.653, 1.653, 1.653, 1.653, 1.653, 1.603, 1.610, 1.615, 1.615,
     *    1.615, 1.615, 1.615, 1.615, 2.673, 2.685, 2.697, 2.708, 2.771,
     *    2.823, 2.877, 2.932, 2.989, 2.989, 24*2.989/
C
C
      DATA LOWYR,MAXCIG/1966,56/
      DATA HDCFAC_FLG/0/
      DATA HDCFAC2/1064*1.0/
C
      END
