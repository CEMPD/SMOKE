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
C                  /VDATA / VERSON
C
C  Notes:
C
C  12-Dec-1994 @ CSC-rwk Req. 2-000 Include file VDATA was added to make it
C  easier to change the name of the version of the model.  This was done for quality
C  control.
C
C
      INCLUDE 'PROJEC.I'
      INCLUDE 'VDATA.I'
C
      WRITE(IOUOUT,200) PROJID,VERSON
  200 FORMAT('1',20A4/' ',A80)
C
      RETURN
      END
