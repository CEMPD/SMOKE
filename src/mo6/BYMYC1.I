C 
C  Changes: (Last change first) 
C 
C 16 Jun 00 @EPA-BAG  Removed IVPICK. This variable is obsolete in Mobile6. 
C 25 Feb 97 MLA  Removed IMPICK. This variable is obsolete in Mobile6. 
C 16 Apr 97 - DynTel@MLA  Expanded arrays BYCUM and IVPICK to accomodate  
C           the new BUS vehicle types. 
C  9 Dec 96 DynTel@MLA Expanded arrays IVPICK and BYCUM to accomodate the  
C           new Heavy Duty Gasoline, Heavy Duty Diesel, and LDGT vehicle  
C           categories. 
C 
C The data is initialized in BD35. 
C 
      REAL BYCUM 
      COMMON /BYMYC1/ BYCUM(25,28) 
