      SUBROUTINE A8TO11(IC) 
C 
C  A8TO11 adjusts effects category sizes when 1 or more of the no overlap 
C  categories 8 to 11 are negative. 
C 
C  Called by EFFGRP. 
C 
C  Calls C8TO12. 
C 
C  Changes : 
C 
C 
C  Input on call: 
C 
C    argument list: IC (flags the negative category size) 
C    common blocks: 
C    /LOOKUP/ JPGD 
C    /MYCODE/ IAY 
C    /SIZCAL/ CSIZE 
C    /SIZPAR/ NTERMS,INDXCS 
C 
C  Output on return: 
C 
C    common block: /SIZCAL/ CSIZE 
C 
C  Local array subscripts: 
C 
C  ADJFAC(2)  -  ADJFAC ( IADJ ) 
C 
C  Local variable / array dictionary: 
C 
C   Name   Type              Description 
C  ------  ----  ------------------------------------------------------- 
C  ADJFAC   R    adjustment factor computed and used to modify CSIZE 
C  DIV      R    divisor in ADJFAC(1) formula = sum of subset of CSIZE 
C  NT       R    number of CSIZE terms in DIV 
C 
C  Notes: 
C 
C  None. 
C 
C 
      IMPLICIT NONE 
C 
      INCLUDE 'LOOKUP.I' 
      INCLUDE 'MYCODE.I' 
      INCLUDE 'SIZCAL.I' 
      INCLUDE 'SIZPAR.I' 
C 
      REAL ADJFAC, DIV 
      INTEGER IS, IC, ICS, JCS, IAD, NT 
C 
      DIMENSION ADJFAC(2) 
C 
C  Sum up the same terms subtracted off to get CSIZE(IC,IAY). 
C 
      IS=IC-7 
      NT=NTERMS(IS) 
      DIV=0.0 
      DO 10 ICS=1,NT 
      JCS=INDXCS(ICS,IS) 
      DIV=DIV+CSIZE(JCS,IAY) 
   10 CONTINUE 
C 
C  Divide the sum into the tamper rate for the no overlap case (no other 
C  tamper types occurring) case that is being computed.  The result is the 
C  first adjustment factor.  One minus the result is the second adjustment 
C  factor. 
C 
      IAD=IS 
      IF(IAD.EQ.1) IAD=JPGD 
      ADJFAC(1)=ATR(IAD,IAY)/DIV 
      ADJFAC(2)=1.0-ADJFAC(1) 
C 
C  Multiply the second adjustment factor times  the appropriate IC overlap 
C  case and add to the related no-IC-occurring overlap case.  The uncorrected 
C  IC values are required for this step.  Hence, ADJFAC(2) is applied before 
C  ADJFAC(1). 
C 
      IF(IC.EQ.8) CSIZE(6,IAY)=CSIZE(6,IAY)+CSIZE(4,IAY)*ADJFAC(2) 
      IF(IC.EQ.8) CSIZE(7,IAY)=CSIZE(7,IAY)+CSIZE(5,IAY)*ADJFAC(2) 
      IF(IC.EQ.9) CSIZE(2,IAY)=CSIZE(2,IAY)+CSIZE(4,IAY)*ADJFAC(2) 
      IF(IC.EQ.9) CSIZE(3,IAY)=CSIZE(3,IAY)+CSIZE(5,IAY)*ADJFAC(2) 
      IF(IC.EQ.10) CSIZE(1,IAY)=CSIZE(1,IAY)+CSIZE(4,IAY)*ADJFAC(2) 
      IF(IC.EQ.11) CSIZE(1,IAY)=CSIZE(1,IAY)+CSIZE(5,IAY)*ADJFAC(2) 
C 
C  Multiply the first adjustment factor times each of IC's overlap cases to 
C  correct their values.  This completes the adjustment of sizes 1 to 7. 
C 
      DO 20 ICS=1,NT 
      JCS=INDXCS(ICS,IS) 
      CSIZE(JCS,IAY)=CSIZE(JCS,IAY)*ADJFAC(1) 
   20 CONTINUE 
C 
C  Recalculate category sizes 8 to 12.  Use the same code as used on the 
C  previous calculation(s). 
C 
      CALL C8TO12 
C 
C  Size adjustment for category IC is finished.  Note that the process 
C  increases the size by decreasing the related overlap category sizes.  The 
C  other non-overlap categories will either also increase in size or not be 
C  affected.  Fixing one negative 8-11 size will never drive another 8-11 size 
C  negative.  The code therefore does not have to go back and check for such 
C  a result. 
C 
      RETURN 
      END 
