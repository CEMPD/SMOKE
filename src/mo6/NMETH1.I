C 
C  Changes: (Last change first) 
C 
C  4 Apr 97 - DynTel@MLA Expanded arrays ZDMTH2 and ZDMTH3 to accomodate  
C           the new BUS vehicle types.  
C 
C  9 Dec 96 DynTel@MLA Expanded arrays ZDMTH1, ZDMTH2, and ZDMTH3 to  
C           accomodate the new Heavy Duty Gasoline, Heavy Duty Diesel,  
C           and LDGT vehicle categories. 
C 
C 
C  The data is initialized in BD10. 
C 
      REAL ZDMTH1, ZDMTH2, ZDMTH3 
      COMMON /NMETH1/ ZDMTH1(2,3,2,7,5,2),ZDMTH2(2,7,10,2), 
     &                ZDMTH3(7,13,2) 
