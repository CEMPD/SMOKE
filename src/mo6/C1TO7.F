      SUBROUTINE C1TO7(CSFUSE) 
C 
C  C1TO7 calculates overlap categories 1 to 7 for tampering effects. 
C  These cover all the allowable overlap combinations of the 4 base tampering 
C  cases for the bag effects for the given pollutant group. 
C 
C  Called by ADJ12 and EFFGRP. 
C 
C Changes : 
C 
C  19 Aug 98 @DynTel-ZK 2-000 Explicit Typing 
C 
C  Input on call: 
C 
C    argument list: CSFUSE 
C    common blocks: 
C    /LOOKUP/ JPGD 
C    /MYCODE/ IAY 
C    /SIZCAL/ ATR 
C 
C  Output on return: 
C 
C    common block: /SIZCAL/ CSIZE 
C 
C  Notes: 
C 
C  None. 
C 
       
                   
      IMPLICIT NONE 
      INCLUDE 'LOOKUP.I' 
      INCLUDE 'MYCODE.I' 
      INCLUDE 'SIZCAL.I' 
      
C  The first 7 category sizes are computed directly from the size coefficient 
C  and tampering rates. 
      DIMENSION CSFUSE(7) 
      REAL CSFUSE        
      INTEGER IC7 
      DO 10 IC7=1,3 
      CSIZE(IC7,IAY)=CSFUSE(IC7)*ATR(JPGD,IAY) 
   10 CONTINUE 
C 
      DO 20 IC7=4,7 
      CSIZE(IC7,IAY)=CSFUSE(IC7)*ATR(2,IAY) 
   20 CONTINUE 
C 
      RETURN 
      END 
