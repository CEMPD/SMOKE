C 
C  Changes: (Last change first) 
C 
C  9 Dec 96 DynTel@MLA Expanded array BYVTN to accomodate the new Heavy  
C           Duty Gasoline, Heavy Duty Diesel, and LDGT vehicle categories. 
C 
C 16 Apr 97 - DynTel@MLA  Expanded array BYVTN to accomodate the new BUS  
C             vehicle types. 
C 
C The data is initialized in BD35. 
C  
      CHARACTER*26 BYVTN 
      CHARACTER*17 YFM 
      CHARACTER*7 PN,U7 
      CHARACTER*3 PSEC 
      CHARACTER*1 VB 
      COMMON /BYMYC4/ BYVTN(28),YFM,PN(10),U7,PSEC(3),VB 
