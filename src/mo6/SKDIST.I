C 
C  Changes: (Last change first) 
C 
C   3 May 99 @DynTel-bag  2-650 Increased the size of the hourly 
C            index on SOAKDIST to 24. 
C 
C The data is initialized in BD41. 
C 
      REAL SOAKDIST 
      COMMON /SKDIST/ SOAKDIST(24,70,2) 
