C 
C  Changes: (Last change first) 
C 
C 13 Jan 99 - DynTel@RAG 2-664  Moved TLEMI array for pass pressure and 
C             purge tests emissions to RUNLS2.I; expanded TLVMT array to 
C             include hourly and weekday/end indices. 
C 16 Apr 97 - DynTel@MLA  Expanded array TLEMI to accomodate the new GAS 
C           BUS vehicle type. 
C  6 Dec 96 DynTel@MLA Expanded array TLEMI to accomodate the new Heavy  
C           Duty Gasoline and LDGT vehicle categories. 
C 
C 
C  
C The data is initialized in BD34. 
C 
      INTEGER MAXTL 
      REAL    TLVMT, GLVMTU, TLVMTU 
C 
      COMMON /RUNLS1/ TLVMT(6,14,2),TLVMTU(6,14,2),GLVMTU(6,14,2), 
     *                MAXTL 
