C 
C The data is initialized in BD33. 
C 
C  Changes: (Last change first) 
C 
C 22 Jun 99 @DynTel-MLA 2-694  Removed obsolete variables FEET and ALT. 
C  4 Jun 98 @DynTel-MLA 2-659  Added IROLD for Mobile6. This variable 
C           is set in MOBILE and PRCSCN. 
C 
      INTEGER IREJN, INITPR, IROLD 
C 
      COMMON /REGION/ IREJN,IROLD,INITPR                                                       
