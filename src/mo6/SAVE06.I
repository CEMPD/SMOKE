C 
C  Changes: (Last change first) 
C 
C 25 Jan 2000 @EPA-djb reduced array sizes. 
C  9 Dec 96 DynTel@MLA Expanded array VA011 to accomodate the new Heavy 
C           Duty Gasoline, Heavy Duty Diesel, and LDGT vehicle categories. 
C 
C 
C  The data is initialized in BD12. 
C 
      REAL VA011, VA012, VA014 
      INTEGER IVA013 
      COMMON /SAVE06/ VA011(16),VA012(7,3),IVA013(4),VA014(4) 
