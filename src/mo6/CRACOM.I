C 
C  Changes: (Last change first) 
C 
C 15 Nov 96 DynTel@MLA Expanded arrays CRAVAL and MYGCRA to accomodate the  
C           new Heavy Duty Gasoline, Heavy Duty Diesel, and LDGT vehicle  
C           categories. 
C 
C 10 Apr 97 - DynTel@MLA  Expanded arrays CRAVAL and MYGCRA to accomodate  
C           the new BUS vehicle types. 
C 
C The data is initialized in BD15. 
C 
      INTEGER MAXCRA,MYGCRA 
      REAL CRAVAL 
      COMMON /CRACOM/ MAXCRA,CRAVAL(3,28,2),MYGCRA(3,28,2)                                             
