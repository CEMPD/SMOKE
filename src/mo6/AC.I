C 
C  Changes: (Last change first) 
C 
C 14 Dec 99 @DynTel-ddj 1-013  Variables for AC correction calculation. 
C 
C The initialization of the data can be found in BD_AC. 
C 
      INTEGER SUNRISE 
      INTEGER SUNSET 
      INTEGER PEAKSUN 
C 
      REAL    CLD_CVR 
      REAL    DMND 
      REAL    FULL_A 
      REAL    FULL_B 
      REAL    NTFCNNG_AC_UNITS 
      REAL    VHCLES_AC_EQPD 
C 
      COMMON /AC_TIMES/ SUNRISE, SUNSET, PEAKSUN(2) 
      COMMON /AC_DATA/  CLD_CVR, FULL_A(10), FULL_B(10), 
     &                  VHCLES_AC_EQPD(2,27), NTFCNNG_AC_UNITS(2,20), 
     &                  DMND(3,5) 
 
