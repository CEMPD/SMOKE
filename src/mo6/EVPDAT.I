C 
C  Changes: (Last change first) 
C 
C 10 Apr 97 - DynTel@MLA  Expanded arrays EVPTGS and CCEMI to accomodate  
C           the new GAS BUS vehicle type. 
C 
C 26 Nov 97 DynTel@MLA Expanded arrays EVPTGS and CCEMI to accomodate the  
C           new Heavy Duty Gasoline and LDGT vehicle categories. 
C 
C The data is initialized in BD25. 
C 
      REAL EVPAER,EVPTGS,CCEMI 
      COMMON /EVPDAT/ EVPAER(3,2),EVPTGS(3,2,14),CCEMI(7,14)                                           
