      MODULE LEAKYDAT 
C 
C    This module initializes the liquid leak flag, the liquid leak parameters 
C    for the liquid leak equations , and phase in parameters for liquid 
C    evaporative emissions calculations. 
C 
C  Changes (Last change first): 
C 
C  28 Aug 00 AIR Task 08: Updated the Phase-In values. Comments. 
C 
C  MODULE dictionary: 
C 
C     Name           Type              Description 
C    ------          ----  --------------------------------------------- 
C  ILIQUID_LEAKS      I    Flag used to enable users to turn off or on the 
C                          addition of liquid leaks in the evaporative 
C                          emissions functions 
C                          values : ILIQUID_LEAKS=0 sets liquid leak fractions 
C                                                   equal to zero, which effectively 
C                                                   eliminates liquid leak emissions 
C                                   ILIQUID_LEAKS=1 sets liquid leak fractions does not 
C                                                   set leak emissions equal to zero 
C                                                   (this is the default value) 
C  GLL_COEF           R   parameters used in the determination of the liquid leak fractions 
C                         (PROGRAMMING NOTE: The C coefficients for MY99+ are half those 
C                          of MY95-. This accounts for extra durability.) 
C  PHASE_IN           R   1995 phase in fractions for model years 1996 thru 1998 
C 
C  Local array subscripts: 
C 
C  GLL_COEF(4,4)          -    (COEFFICIENT #, MYRNG) 
C  PHASE_IN(3)            -    (IMY=1->1996;IMY=2->1997;IMY=3->1998) 
C 
      IMPLICIT NONE 
C 
      INTEGER :: ILIQUID_LEAKS 
C 
      REAL, DIMENSION(4,4) :: GLL_COEF = RESHAPE( 
C DIURNAL 
C   1995 and earlier   A       B       C    D 
     &            (/0.08902,414.613,-0.3684,0.0, 
C   1999 and later     A       B       C    D 
     &              0.08902,414.613,-0.1842,0.0, 
C RUNNING 
C  1995 and earlier   A       B      C     D 
     &              0.06,   120.0, -0.4  ,0.0, 
C  1999 and later     A       B      C     D 
     &              0.06,   120.0, -0.2,  0.0/), 
     & (/4,4/)) 
C 
C Note-these are the unenhanced coeffs, i.e., fractions for 1995 emissions 
      REAL, DIMENSION(3) :: PHASE_IN = (/0.70,0.45,0.10/) 
C 
C 
      END MODULE LEAKYDAT 
 
