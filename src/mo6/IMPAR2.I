C 
C  Changes: (Last change first) 
C 
C 25 Aug 00 AIR Task 08  Expanded arrays ILDT and ITEST to 
C           accomodate 7 IM programs. 
C  8 Dec 98 @DynTel-MLA 2-663  Expanded arrays ILDT and ITEST to 
C           accomodate 5 IM programs. 
C 
C 14 Apr 97 - DynTel@MLA  Expanded array ILDT to accomodate the new 
C             GAS BUS vehicle type. 
C 
C  3 Dec 96 DynTel@MLA Expanded array ILDT to accommodate the new Heavy 
C           Duty Gasoline and LDGT vehicle categories. 
C 
C 14 Apr 97 - DynTel@MLA  Expanded array ILDT to accommodate the new 
C             GAS BUS vehicle type. 
C 
C 
C The data is initialized in BD_IM1. 
C 
      INTEGER ILDT,ITEST 
      COMMON /IMPAR2/ ILDT(14,7),ITEST(7) 
