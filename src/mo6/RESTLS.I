C 
C  Changes : (Last change first) 
C 
C  19 JAN 00 @Dyntel-jwr 1-017 Expanded RSTLS1, added RSTLS2 
C 30 Oct 98 @DynTel-bag 2-677 Added hourly resting loss model parameters, 
C           RSTLS1(model year range,pressure /purge pass/fail,fuel delivery system). 
C 16 Apr 97 @DynTel-MLA  Expanded array PTGRST to accomodate the new GAS  
C           BUS vehicle type. 
C  6 Dec 96 @DynTel-MLA Expanded array PTGRST to accomodate the new Heavy  
C           Duty Gasoline and LDGT vehicle categories. 
C 
C 
C The data is initialized in BD34. 
C 
      REAL PTGRST, RESTA, RESTB,RSTLS1,RSTLS2                                                  
      COMMON /RESTLS/ PTGRST(13,15,2),RESTA(2),RESTB,RSTLS1(4,4,2),                                                  
     & RSTLS2(4,4)