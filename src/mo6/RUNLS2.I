C 
C  Changes: (Last change first) 
C 
C 13 Jan 99 - DynTel@RAG 2-664  Moved TLEMI array for pass pressure and purge tests emissions 
C             from RUNLS1.I. 
C 
C The data is initialized in BD34. 
C 
      REAL TLEMI, TLFAIL 
      COMMON /RUNLS2/ TLEMI(6,3,4,4,4,14),TLFAIL(6,4,4) 
