C 
C  Changes: (Last change first) 
C 
C  8 Apr 99 @DynTel-DDJ 2-696 Added hourly index to the RVFPCF array and 
C 16 Apr 97 - DynTel@MLA  Expanded arrays MYGRVP and RVPCF to accomodate  
C           the new GAS BUS vehicle types. 
C  3 Dec 96 DynTel@MLA Expanded arrays MYGRVP and RVPCF to accomodate the  
C           new Heavy Duty Gasoline and LDGT vehicle categories. 
C 
C 
C 
C The data is initialized in BD31. 
C 
      INTEGER MYGRVP 
      REAL RVPCF 
      COMMON /RVPEX2/ MYGRVP(4,14),RVPCF(25,3,14,24) 
