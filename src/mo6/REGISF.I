C 
C  Changes: (Last change first) 
C 
C 25 Jan 2000 @EPA-djb expanded arrays. 
C 14 Apr 97 - DynTel@MLA  Expanded array GSFVCT to accomodate the new BUS  
C           vehicle types. 
C  3 Dec 96 DynTel@MLA Expanded array GSFVCT to accomodate the new Heavy 
C           Duty Gasoline, Heavy Duty Diesel, and LDGT vehicle categories. 
C 
C 
C 
C The data is initialized in BD27. 
C 
      REAL GSFVCT, USRGSF, GLUSGF 
      COMMON /REGISF/ GSFVCT(28),USRGSF(25,14),GLUSGF(25,14) 
