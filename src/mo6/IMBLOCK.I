C 
C  Changes: (Last change first) 
C 
 
C 29 Aug 00  @EPA-elg  Created IMBLOCK 
C 
C The initialization of the data can be found in BD_IM2.FOR. 
C 
C     FUELCF   R    Fuel correction Factor 
C 
      REAL, DIMENSION(2,2,10,3,24) :: FUELCF 
C 
      COMMON /IMBLOCK1/ FUELCF 
