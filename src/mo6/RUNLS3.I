C 
C  Changes: (Last change first) 
C 
C 16 Apr 97 - DynTel@MLA  Expanded array MYGRNL to accomodate the new GAS  
C           BUS vehicle types. 
C  6 Dec 96 DynTel@MLA Expanded array MYGRNL to accomodate the new Heavy 
C           Duty Gasoline and LDGT vehicle categories. 
C 
C 
C The data is initialized in BD34. 
C 
      INTEGER MYGRNL, MAXRNL 
      REAL    RNLRVP, RNLTEM, RNLSPD               
C     
      COMMON /RUNLS3/ MYGRNL(4,14),RNLRVP(4),RNLTEM(4),RNLSPD(3),MAXRNL                               
