C 
C  Changes : (Last change first) 
C  18 Jan 01 @EPA-bag  2007 rule for heavy duty vehicles 
C  08 Aug 00 @EPA-elg Added NO_EVAPCALC to shut down evap calculations. 
C  20 Jul 00 @EPA-bag running loss: added NO_RSTLOS for zeroing resting losses 
C  19 jan 00 @Dyntel-jwr 1-017 replaced resting loss total 
C             by RSTLSS_EF 
C 
C The data is initialized in BD34. 
C 
      REAL RSTLSS_EF, RSTLOS_FR_RNLOSS 
      INTEGER NO_RSTLOS, NO_EVAPCALC 
      COMMON /RESTLS1/RSTLSS_EF,NO_RSTLOS,NO_EVAPCALC,RSTLOS_FR_RNLOSS 
