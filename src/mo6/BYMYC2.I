C 
C  Changes: (Last change first) 
C 
C   5 Aug 99 @DynTel-ddj 2-699 Removed BYTAM array. 
C  22 Feb 99 @DynTel-ddj 2-696 Added hourly index to BYTAM. 
C  9 Dec 96 DynTel@MLA Expanded arrays BYBEF4, BYFER, and BYTAM to accomodate  
C           the new Heavy Duty Gasoline, Heavy Duty Diesel, and LDGT vehicle  
C           categories. 
C 
C 16 Apr 97 - DynTel@MLA  Expanded arrays BYBEF4, BYFER, and BYTAM to  
C           accomodate the new BUS vehicle types. 
C 
C 23 Apr 98 @DynTel-RJD  Added an extra dimension to the BYBEF4 array to 
C           distinguish between starts and running emissions. 
C 
C The data is initialized in BD35. 
C 
      REAL BYBEF4,BYFER 
      COMMON /BYMYC2/ BYBEF4(3,25,28,2),BYFER(3,25,28) 
