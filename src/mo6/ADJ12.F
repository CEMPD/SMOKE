      SUBROUTINE ADJ12(CSF1ST) 
C 
C  ADJ12 adjusts all category sizes to eliminate a negative category 12 
C  size. 
C 
C  Called by EFFGRP. 
C 
C  Calls C1TO7 and C8TO12. 
C 
C  Changes : (Last change first) 
C 
C   15 Jul 98 @DynTel-ZK 2-000 Explicit Typing 
C 
C  Input on call: 
C 
C    argument list: CSF1ST 
C    common blocks: 
C    /LOOKUP/ IPG 
C    /MYCODE/ IAY 
C    /SIZCAL/ CSIZE 
C 
C  Output on return: 
C 
C    common block: /SIZCAL/ CSIZE 
C 
C  Local array subscripts: 
C 
C  CSF2ND(7)  -  CSF2ND ( IC7 ) 
C 
C  Local variable / array dictionary: 
C 
C   Name   Type              Description 
C  ------  ----  ------------------------------------------------------- 
C  CSF2ND   R    category size factor for 2nd calculation of CSIZE; 
C                passed to C1TO7 
C  ADJSUM   R    adjustment factor to scale down 1st calculation of CSIZE 
C 
C  Notes: 
C 
C  None. 
C 
C 
      IMPLICIT NONE 
      REAL CSF2ND, ADJSUM, CSF1ST  
      INTEGER IC11, IC7, IC 
C 
      DIMENSION CSF1ST(7,2,3),CSF2ND(7) 
C 
      INCLUDE 'LOOKUP.I' 
      INCLUDE 'MYCODE.I' 
      INCLUDE 'SIZCAL.I' 
C 
C  Adjustment sum is sum of the positive emission effect category sizes. 
C 
      ADJSUM=0.0 
      DO 10 IC11=1,11 
      ADJSUM=ADJSUM+CSIZE(IC11,IAY) 
   10 CONTINUE 
C 
C  Apply adjustment sum to all the 1st calculation category size factors 
C  to get the 2nd calculation size factors. 
C 
      DO 20 IC7=1,7 
      CSF2ND(IC7)=ADJSUM*CSF1ST(IC7,IPG,IGCSF) 
   20 CONTINUE 
C 
C  Now recalculate all 12 category sizes.  category size 12 will go to 0.0 
C  (the sum of the first 11 sizes will become 1.0). 
C 
      CALL C1TO7(CSF2ND) 
      CALL C8TO12 
C 
C  Check if any of the sizes for 8 to 11 went negative on the recalculation. 
C 
      DO 30 IC=8,11 
      IF(CSIZE(IC,IAY).LT.0) CALL A8TO11(IC) 
   30 CONTINUE 
C 
C  2nd adjustments to 8-11 will not drive 12 negative, due to the changed 
C  CSF's.  Hence, we do not need to check 12 again.  Checking again would make 
C  the whole procedure recursive and probably open to infinite loops. 
C 
      RETURN 
      END 
