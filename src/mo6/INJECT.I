C 
C  Changes: (Last change first) 
C 
C  09 May 01 @EPA bag Added fractions and additional model years for post 2004 fuel injection vehicle percentages 
C 16 Apr 97 - DynTel@MLA  Expanded arrays PFI and TBI to accomodate the  
C           new GAS BUS vehicle types. 
C 
C  3 Dec 96 DynTel@MLA Expanded arrays PFI and TBI to accomodate the new  
C           Heavy Duty Gasoline and LDGT vehicle categories. 
C 
C The data is initialized in BD17 and in BD30. 
C 
      REAL TBI,PFI 
      INTEGER MAXFIY,MINFIY 
      COMMON /INJECT/ TBI(17,14),PFI(17,14),MAXFIY,MINFIY 
