C Version @(#)$Id$ $Source$ $Date$ 

      BLOCK DATA BD29
C
C  BLOCK DATA Subprogram 29: /HDCCOM/
C
C  BLOCK DATA 29 initializes grams-per-brakehorsepower-hour to grams-per-mile
C  units conversion variables.
C
C  Common block array subscripts:
C
C  HDCFAC(37,2)  -  HDCFAC ( IG, IVHD )
C
C  Common block dictionary:
C
C   Name   Type              Description
C  ------  ----  -------------------------------------------
C  /HDCCOM/:
C
C  HDCFAC   R    HDV gpbhp to gpm conversion factors
C  LOWYR    I    lowest my for which HDCFAC is defined.  All my < LOWYR are
C                assigned LOWYR's conversion factors.
C  MAXCIG   I    maximum HDCFAC index allowed = 1st dimension of HDCFAC = 37
C                = index of highest my for which HDCFAC is defined, currently
C                1997.  All my > 1997 are assigned MAXCIG's factors.
C
C  Local array subscripts:
C
C  HDCGAS(37)  -  HDCGAS ( IG )
C  HDCDIE(37)  -  HDCDIE ( IG )
C
C  Local array dictionary:
C
C   Name   Type                Description
C  ------  ----  -----------------------------------------
C  HDCGAS   R    HDGV gpbhp to gpm units conversion factors
C  HDCDIE   R    HDDV gpbhp to gpm units conversion factors
C
C  Notes:
C
C  None.
C
C
      COMMON /HDCCOM/ HDCFAC(37,2),LOWYR,MAXCIG
C
      DIMENSION HDCGAS(37),HDCDIE(37)
C
      EQUIVALENCE
     *   (HDCFAC(1,1),HDCGAS(1)),
     *   (HDCFAC(1,2),HDCDIE(1))
C
C  Gpbhp to gpm units conversion factors
C
C  For HDGV
C
      DATA HDCGAS/
     *    1.548, 1.536, 1.536, 1.527, 1.516, 1.518, 1.497, 1.479, 1.449,
     *    1.449, 1.437, 1.419, 1.422, 1.313, 1.244, 1.073, 1.057, 1.022,
     *    0.961, 0.935, 0.912, 0.884, 0.907, 0.896, 0.894, 0.897, 0.895,
     *    0.894, 0.893, 0.893, 0.892, 0.891, 0.890, 0.889, 0.887, 0.886,
     *    0.885/
C
C  For HDDV
C
      DATA HDCDIE/
     *    2.850, 2.858, 2.874, 2.890, 2.900, 2.964, 2.995, 3.074, 3.100,
     *    3.161, 3.197, 3.188, 3.213, 3.146, 3.179, 3.246, 3.187, 2.999,
     *    2.716, 2.698, 2.376, 2.277, 2.406, 2.211, 2.214, 2.167, 2.132,
     *    2.099, 2.066, 2.050, 2.033, 2.033, 2.033, 2.039, 2.039, 2.037,
     *    2.036/
C
      DATA LOWYR,MAXCIG/1962,37/
C
      END
