C 
C  Changes: (Last change first) 
C 
C 14 Apr 97 - DynTel@MLA  Expanded array DSIZE to accomodate the new 
C             GAS BUS vehicle type. 
C 
C  3 Dec 96 DynTel@MLA Expanded array DSIZE to accomodate the new Heavy  
C           Duty Gasoline and LDGT vehicle categories. 
C 
C The data is initialized in BD_IM1. 
C   
      REAL DSIZE 
      COMMON /IM240P/ DSIZE(25,14)                                                                    
