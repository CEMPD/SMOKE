C 
C  Changes: (Last change first) 
C 
C 28 Aug 00 @EPA-djb Bug 233, Removed HDSAL. 
C 15 Nov 96 DynTel@MLA Expanded array IVVMAP to accomodate the new Heavy 
C           Duty Gasoline and LDGT vehicle categories. 
C 
C  7 Apr 97 - DynTel@MLA  Expanded array IVVMAP to accomodate the new   
C             GASBUS vehicle type. 
C 
C 
C  The data is initialized in BD11. 
C 
      REAL HDWGT 
      INTEGER IVVMAP 
      COMMON /EVAPHD/ HDWGT(2),IVVMAP(15) 
