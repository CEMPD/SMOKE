C 
C  Changes: (Last change first) 
C 
C 20 Jun 00 AIR Task 03: Removed OMTTAM. 
C 02 Jun 00 AIR Task 03: Removed OMTCF and added TCFSTART, TCFRUN, and 
C           TCF12 
C 23 May 00 @EPA-elg   Removed TEMP_ST array. 
C  8 Apr 99 @DynTel-DDJ 2-696 Added hourly index to the OMTCF array and 
C           added DIST_SUM AND TEMP_ST arrays 
C 21 Apr 98 @DynTel-RJD Added an extra dimension to the OMTCF array to 
C            distinguish between start and running emissions. 
C 26 Nov 97 DynTel@MLA Expanded arrays OMTCF and OMTTAM to accomodate the  
C           new Heavy Duty Gasoline, Heavy Duty Diesel, and LDGT vehicle  
C           categories. 
C 10 Apr 97 - DynTel@MLA  Expanded arrays OMTCF and OMTTAM to accomodate  
C           the new GAS BUS vehicle type. 
C 
C 
C The data is initialized in BD25. 
C 
      REAL TCFSTART, TCFRUN, TCF12, DIST_SUM 
      COMMON /OMTCOM/ TCFSTART(69), TCFRUN, TCF12, DIST_SUM(24) 
