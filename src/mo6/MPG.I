C 
C  Changes: (Last change first) 
C 
C 06 Nov 01: AIR Task 27: New INCLUDE. 
C 
C The data is initialized in BD48. 
C 
      LOGICAL MPGFLAG 
      REAL    ROADFE 
C 
      COMMON /MPG/ ROADFE(99,28,2), MPGFLAG 
