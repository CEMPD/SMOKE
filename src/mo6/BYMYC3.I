C 
C  Changes: (Last change first) 
C 
C  9 Dec 96 DynTel@MLA Expanded arrays BYEVAP, BYREFL, BYRSTL, and BYRUNL to  
C           accomodate the new Heavy Duty Gasoline, Heavy Duty Diesel, and  
C           LDGT vehicle categories. 
C 
C 16 Apr 97 - DynTel@MLA  Expanded arrays BYEVAP, BYREFL, BYRSTL, and   
C           BYRUNL to accomodate the new BUS vehicle types. 
C 
C The data is initialized inBD35. 
C   
      REAL BYEVAP,BYRUNL,BYREFL,BYRSTL 
      COMMON /BYMYC3/BYEVAP(25,28),BYRUNL(25,14),BYREFL(25,14), 
     &           BYRSTL(25,28) 
