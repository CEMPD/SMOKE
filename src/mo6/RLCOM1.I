C 
C  Changes: (Last change first) 
C 
C 26 Feb 01 @EPA-djb Added HDGV2b to the ORVR arrays (OBIPY & OBPIP). 
C  2 Jun 99 @DynTel-MLA 2-694  Removed obsolete variable OBDF. 
C 14 Apr 97 - DynTel@MLA  Expanded arrays OBDF and ROADFE to accomodate  
C           the new GAS BUS vehicle type. 
C  3 Dec 96 DynTel@MLA Expanded arrays OBDF, OBPIP, OBPIY, and ROADFE to  
C           accomodate the new Heavy Duty Gasoline and LDGT vehicle  
C           categories. 
C 
C 
C The data is initialized in BD28. 
C 
      INTEGER OBPIY 
      REAL ROADFE, DISPL, SPILL, OBED, OBES, OBPIP 
      COMMON /RLCOM1/ ROADFE(32,14),DISPL,SPILL,OBED,OBES, 
     *                OBPIY(3,6), OBPIP(3,6) 
