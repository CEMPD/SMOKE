C 
C  Changes: (Last change first) 
C 
C 03 Nov 99 DynTel@BAG 1-006 Removed partial diurnal array, PARD, and 
C           evaporative sum array, EVP 
C 4 Mar 99  DynTel@BAG 2-664 Added resting loss and running loss variables 
C 10 Mar 98 DynTel@RJD 2-622 Added integer labels for emission factors to 
C           facilitate the efficient calculation of emission factors for 
C           composite vehicles in EFCVT. 
C 
C  3 Apr 97 - DynTel@MLA  Expanded arrays GREVP and GRPD to accomodate the  
C             new BUS vehicle types. 
C 
C 15 Nov 96 DynTel@MLA Expanded arrays GREVP and GRPD to accomodate the new  
C           Heavy Duty Gasoline, Heavy Duty Diesel, and LDGT vehicle types. 
C 
C  The data is initialized in BD11.FOR 
C 
      INTEGER IGREVP 
      REAL GREVP 
 
      COMMON /EVAPGR/ GREVP(6,29),IGREVP 
