
C Version "@(#)$Id$ $Source$ $Date$ 

      SUBROUTINE OUTHDR(IOUOUT)
C
C  OUTHDR writes the output header ( = run title ) on unit IOUOUT.
C
C  Called by CONSEC.
C
C  Input on call:
C
C    parameter list: IOUOUT
C    common block: /PROJEC/ PROJID
C
C  Notes:
C
C  None.
C
C
      CHARACTER*4 PROJID
      COMMON /PROJEC/ PROJID(20)
C
cmh   WRITE(IOUOUT,200) PROJID
  200 FORMAT('1',20A4/' ','MOBILE5a (26-Mar-93)')
C
      RETURN
      END
