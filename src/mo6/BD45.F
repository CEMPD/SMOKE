      BLOCK DATA BD45 
C 
C  EMFAC2000v2.02 PM Speed Correction Factor Coefficients. 
C 
C  Changes : 
C  28 Mar 01 AIR Task 11: Due to Lahey bug, had to replace 
C            BLOCK DATA with GETPMSCFC and PMSCFC.D. 
C  20 Mar 01 AIR Task 07/2K: New Module 
C 
C  PMSCFC(46,4,3,2,2,28) - PMSCFC (MY,ICOEF,IFDS,ICAT,MODE,IVEH) 
C 
C  C and ICOEF are for SCF coefficient, as follows 
C    A, B, C, D 
C 
C  IFDS is for fuel delivery system, as follows 
C    1 = Carb/Diesel 
C    2 = TBI 
C    3 = MPFI 
C 
C  ICAT is for catalyst, as follows: 
C    1 = no catalyst 
C    2 = catalyst 
C 
C  M and Mode are for Mode, as follows: 
C    1 = FTP 
C    2 = Unified Cycle 
C 
C  Common Block Dictionary: 
C 
C   Name    Type                  Description 
C  ------   ----  ----------------------------------------------------- 
C /PART5/ 
C   PMSCFC   R    PM Speed Correction Factor Coefficients 
C 
      IMPLICIT NONE 
C 
      INCLUDE 'PART5.I' 
C 
      DATA PMSCFC   /61824*0./ 
      DATA PMSPDFLG /0/ 
C 
      END 
