C 
C The data is initialized in BD38. 
C 
C  Changes: (Last change first) 
C 
C  4 Jun 98 @DynTel-MLA 2-659  Added MEVOLD for Mobile6. This variable 
C           is set in MOBILE and PRCSCN. 
C 
      INTEGER MEVAL,MEVOLD 
      CHARACTER *4 MONTH 
      REAL FRC,FRN 
      COMMON /EVAL/ MEVAL,MEVOLD,MONTH(12),FRC(28),FRN(28) 
