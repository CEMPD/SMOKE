C 
C  Changes: (Last change first) 
C 
C 27 Jul 00 @EPA-djb Bug 208, removed MAXSTD. 
C  2 Jun 00 @EPA-ddj, Bug 86, Expanded array KSTD81 to accomodate the 
C           HDGVs. 
C  3 Apr 97 - DynTel@MLA  Expanded arrays MYGSTD and KSTD81 to accomodate 
C             the new GAS BUS vehicle type. 
C 
C 15 Nov 96 DynTel@MLA Expanded the arrays MYGSTD and KSTD81 to accomodate 
C           the new Heavy Duty Gasoline and LDGT vehicle categories. 
C 
C 
C The data is initialized in BD11. 
C 
      INTEGER MYGSTD,KSTD81 
      COMMON /EVASTD/ MYGSTD(5,15),KSTD81(14,3) 
