C Version "@(#)$Id$ $Source$ $Date$ 

      SUBROUTINE ATPEFF(*)
C
C  ATPEFF applies the anti-tampering effectiveness factors to the 11
C  positive tampering effect categories.
C
C  Called by EFFGRP.
C
C  Input on call:
C
C    common blocks:
C    /EGSCAL/ AER
C    /LOOKUP/ IPG
C    /MYCODE/ IAY
C    /SIZCAL/ CSIZE
C
C  Output on return:
C
C    common block: /SIZCAL/ CSAE
C
C  Local variable dictionary:
C
C   Name   Type              Description
C  ------  ----  -------------------------------------------------------
C  C1ST11   R    sum of first 11 category sizes, assigned as 12th category
C                size
C
C  Notes:
C
C  None.
C
C
      INCLUDE 'EGSCAL.I'
      INCLUDE 'LOOKUP.I'
      INCLUDE 'MYCODE.I'
      INCLUDE 'SIZCAL.I'
C
C  The sum of the 11 positive categories will be subtracted from 1.0 to get
C  the 12th, no effect, category.
C
      C1ST11=0.0
C
C  ATP effectiveness modified size ICOL will be the product of the entire
C  unmodified row size vector times the ICOL'th column of the effectiveness
C  matrix for the given pollutant group and pre/post ATP year case.
C
C  The effectiveness matrices are currently identical.  Each has the major
C  diagonal all 1.0's and all other cells 0.0.  the result of this subroutine
C  is thus to merely echo CSIZE into CSAE with no changes.
C
      DO 20 ICOL=1,11
      CSAE(ICOL)=0.0
C
      DO 10 IROW=1,11
      CSAE(ICOL)=CSAE(ICOL)+CSIZE(IROW,IAY)*AER(IROW,ICOL,IPG,IAY)
   10 CONTINUE
C
      C1ST11=C1ST11+CSAE(ICOL)
   20 CONTINUE
C
C  The no effects case is not used hereafter.  It is only computed for
C  error checking.
C
      CSAE(12)=1.0-C1ST11
      RETURN1
      END
