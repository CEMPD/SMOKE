      INTEGER FUNCTION STANDARD(STD) 
C 
C    STANDARD converts from 1 to 10 STD standards to 1 of 4 ISTD 
C    standards. 
C 
C    Called by SUL_COR 
C 
C  Changes: (Last change first) 
C 
C    5 May 99 @DynTel-ddj #2-699, New function 
C 
C 
C  Input on CALL: 
C 
C     Argument List:  STD 
C 
C 
C   Output on RETURN: 
C 
C     STANDARD (value of 1, 2, 3 or 4) used in sulfur correction 
C     function SUL_COR. 
C 
C    STD to ISTD relationship: 
C 
C    STD 1      >>   ISTD ==   1  for Tier0. 
C    STD 2 - 3  >>   ISTD ==   2  for Tier1 Inter, Tier1. 
C    STD 4      >>   ISTD ==   5  Standard Tier 2. 
C    STD 5 - 8  >>   ISTD ==   3  TLEV Inter, TLEV, LEV and 
C                                 LEV Inter. 
C    STD 9 - 10 >>   ISTD ==   4  ULEV and ULEV Inter. 
C    STD > 10   >>   ISTD ==   6  Default (ZEV, no adjustment). 
C 
C 
      IMPLICIT NONE 
C 
      INTEGER, INTENT(IN)       :: STD 
      INTEGER  ISTD 
C 
      IF(STD.GT.10) THEN 
         ISTD = 6 
      ELSE 
         IF (STD.GE.9) THEN 
            ISTD = 4 
         ELSE 
            IF (STD.GE.5) THEN 
               ISTD = 3 
            ELSE 
               IF (STD.EQ.1) THEN 
                  ISTD = 1 
               ELSE 
                  IF (STD.LT.4) THEN 
                     ISTD = 2 
                  ELSE 
                     ISTD = 5 
                  END IF 
               END IF 
            END IF 
         END IF 
      END IF 
      STANDARD = ISTD 
      RETURN 
      END 
       
 
