C 
C  PART5 contains the  EMFAC2000v2.02 PM Speed Correction Factor Coefficients. 
C 
C  Changes: (Last change first) 
C 
C  23 May 01: AIR Task 17: Added PM spped correction factor flag. 
C  20 Mar 01: AIR Task 11: New INCLUDE. 
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
C   PMSPDFLG I    PM Speed Correctin Factor Flag (0=off (default, 1=on) 
C 
C The data is initialized in BD45. 
C 
      INTEGER PMSPDFLG 
      REAL PMSCFC 
C 
      COMMON /PART5/ PMSCFC(46,4,3,2,2,28),PMSPDFLG 
