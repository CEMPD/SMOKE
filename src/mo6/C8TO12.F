      SUBROUTINE C8TO12 
C 
C  C8TO12 calculates category sizes 8 to 12 for tampering effects. 
C  Categories 8 to 11 are the no overlap cases (AIRS only, CATS only, etc.). 
C  Category 12 is no tampering at all case. 
C 
C  Called by ADJ12, A8TO11 and EFFGRP. 
C 
C  Changes : 
C 
C  19 Aug 98 @DynTel-ZK 2-000 Explicit Typing 
C 
C  Input on call: 
C 
C    common blocks: 
C    /MYCODE/ IAY 
C    /SIZPAR/ NTERMS,INDXCS 
C    /SIZCAL/ CSIZE 
C 
C  Output on return: 
C 
C    common block: /SIZCAL/ CSIZE 
C 
C  Local array subscripts: 
C 
C  SUB(4)  -  SUB ( IS ) 
C 
C  Local variable / array dictionary: 
C 
C   Name   Type              Description 
C  ------  ----  ------------------------------------------------------- 
C  C1ST11   R    sum of the first 11 category sizes 
C  NT       I    number of terms in SUB(IS) 
C  SUB      R    sum of category sizes to be subtracted from ATR to get CSIZE, 
C                given IS, for sizes 8 to 11 
C 
C  Notes: 
C 
C  None. 
C 
      IMPLICIT NONE 
      DIMENSION SUB(4) 
        
      INCLUDE 'LOOKUP.I' 
      INCLUDE 'MYCODE.I' 
      INCLUDE 'SIZPAR.I' 
      INCLUDE 'SIZCAL.I' 
C 
C  Categories 8 to 12 are derived from the first 7.  first sum up the terms to 
C  be subtracted. 
      REAL SUB, C1ST11  
      INTEGER NT, ICS, IC, IS, IC11 
      DO 20 IS=1,4 
      NT=NTERMS(IS) 
      SUB(IS)=0.0 
      DO 10 ICS=1,NT 
      SUB(IS)=SUB(IS)+CSIZE(INDXCS(ICS,IS),IAY) 
   10 CONTINUE 
   20 CONTINUE 
C 
      CSIZE(8,IAY)=ATR(JPGD,IAY)-SUB(1) 
      DO 30 IC=9,11 
      IS=IC-7 
      CSIZE(IC,IAY)=ATR(IS,IAY)-SUB(IS) 
   30 CONTINUE 
C 
C  Assign the 12th category size (= no tampering effect) the residual 
C  percentage. 
C 
      C1ST11=0.0 
      DO 40 IC11=1,11 
      C1ST11=C1ST11+CSIZE(IC11,IAY) 
   40 CONTINUE 
      CSIZE(12,IAY)=1.0-C1ST11 
C 
      RETURN 
      END 
