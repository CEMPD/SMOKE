C 
C  PPEFFD contains the variable that holds the values of the credits for 
C   pressure and purge test effectiveness 
C 
C  Changes: (Last change first) 
C 
C 
C  The data is initialized in BD_IM2. 
C 
C  Common Block Dictionary: 
C 
C   Name    Type                  Description 
C  ------   ----  ----------------------------------------------------- 
C  PPEFF1    R    holds the default values of the credits for pressure and  
C                 purge test effectiveness 
C  PPEFF     R    holds the values that have to be for the current run. 
C  PPDES     C    descriptive record 
C  RBHPPF    I    RBH flag for Pressure / Purge program  
C  Notes: 
       REAL PPEFF, PPEFF1 
       INTEGER RBHPPF 
       CHARACTER*80 PPDES  
       COMMON /PPEFFD/ PPEFF(25,2),PPEFF1(25,2),PPDES,RBHPPF 
