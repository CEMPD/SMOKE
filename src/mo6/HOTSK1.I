C 
C Changes : (Last change first) 
C 
C 16 Nov 1999 @Dyntel 1-014 added GLL_SWITCH DECLARATION AND TO HOTSK1 
C 
C 
C  The HSK data is initialized in BD27. 
C 
      REAL HSK 
      LOGICAL GLL_SWITCH 
      COMMON /HOTSK1/ HSK(24),GLL_SWITCH 
