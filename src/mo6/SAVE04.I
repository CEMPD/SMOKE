C 
C  Changes: (Last change first) 
C 
C 25 Jan 2000 @EPA-djb reduced size of VA008 array. 
C 16 Apr 97 - DynTel@MLA  Expanded arrays VA007 and VA008 to accomodate  
C           the new BUS vehicle types. 
C 20 Nov 96 DynTel@MLA Expanded arrays VA007 and VA008 to accomodate the  
C           new Heavy Duty Gasoline, Heavy Duty Diesel, and LDGT vehicle  
C           categories. 
C 
C 
C 
C 
C The data is initialized in BD12. 
C 
      REAL VA007, VA008 
      COMMON /SAVE04/ VA007(25,28),VA008(25,16) 
