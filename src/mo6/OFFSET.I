C 
C  Changes: (Last change first) 
C 
C 02 Jun 00 AIR Task 03: Changed OFFCO to OFFSET. 
C 08 May 00 AIR Task 03: Put OFFCO back in for MOBILE5b TCF useage. 
C 25 Apr 00 @EPA-elg  Added MOFFMTH Array. 
C 14 Oct 99 @DynTel-ddj 2-714  Removed array OFFCO.  OFFCO is now a 
C           local variable in MAINTCF. 
C 28 Sep 99 @DynTel-ddj for EPA  Added ISTD index to OFFMTH. 
C  8 Apr 99 @DynTel-DDJ 2-696 Added hourly index to the OFFCO array 
C 21 Apr 98 @DynTel-RJD Added an extra dimension to the OFFMTH array to 
C           distinguish between starts and running emissions. 
C 14 Apr 97 - DynTel@MLA  Expanded array OFFMTH to accomodate the new BUS  
C             vehicle types. 
C  3 Dec 96 DynTel@MLA Expanded arrays OFFCO and OFFMTH to accomodate the  
C           new Heavy Duty Gasoline, Heavy Duty Diesel, and LDGT vehicle  
C           categories. 
C 
C 
C 
C The data is initialized in BD27. 
C 
      REAL OFFMTH, MOFFMTH, OFFSET 
      COMMON / OFFSET/ OFFMTH(25,28,2,10,3), MOFFMTH(25,28,2,10,3),  
     *                 OFFSET 
