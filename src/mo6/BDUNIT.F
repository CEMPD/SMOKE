      BLOCK DATA BDUNIT 
C 
C 
C 
C  BLOCK DATA Subprogram UNIT:  /UNIT_PARAM/ 
C 
C  BDUNIT defines and initializes the symbolic names for unit conversions 
C 
C  Changes: (Last change first) 
C 
C 
C  Common block array subscripts: 
C 
C 
C  Notes: 
C 
C 
C  Variable Dictionary: 
C 
C   Name    Type                  Description 
C  ------   ----  ----------------------------------------------------- 
C 
C  /UNIT_PARAM/ 
C 
C GPT       I     grams per trip; has a value of 1 
C GPH       I     grams per hour; has a value of 2 
C GPG       I     grams per gallon; has a value of 3 
C GPM_D     I     grams per mile; has a value of 4 
C GPM       I     daily average grams per mile; has a value of 5 
C GPS       I     grams per start; has a value of 6 
C 
C .................................................................... 
C 
      IMPLICIT NONE 
      INCLUDE 'UNIT_PARAM.I' 
C 
      DATA GPT,GPH,GPG,GPM_D,GPM,GPS/1, 2, 3, 4, 5, 6/ 
C 
      END 
