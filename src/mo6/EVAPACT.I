C 
C  Changes: (Last change first) 
C  20 Jul 00 AIR Task 06: New Include 
C 
C The data is initialized in BD11. 
C 
      REAL HSACTUSER, DSACTUSER, FPPRUSER 
      INTEGER SAWFLAG, HSACTFLAG, DSACTFLAG, FPPRFLAG 
      COMMON /EVAPACT/ HSACTUSER(60,14), DSACTUSER(18,72), FPPRUSER(4), 
     *                 SAWFLAG, HSACTFLAG, DSACTFLAG, FPPRFLAG                                                             
