C 
C  Changes: (Last change first) 
C 
C 25 Jan 2000 @EPA-djb reduced size of JULMYR array. 
C  7 Apr 97 - DynTel@MLA  Expanded arrays AMAR and JULMYR to accomodate  
C           the new BUS vehicle types. 
C 
C 20 Nov 96 DynTel@MLA Expanded the arrays to accomodate the new Heavy Duty 
C           Gasoline, Heavy Duty Diesel, and LDGT vehicle categories. 
C 
C 
      INTEGER NEWCUM 
      REAL AMAR, JULMYR  
      COMMON /MYRSAV/ AMAR(25,28),JULMYR(25,16),NEWCUM 
