C Version @(#)$Id$ $Source$ $Date$ 

      SUBROUTINE OBYHDR(IOUOUT)
C
C  OBYHDR writes the By Model Year Tables page header ( = run + scenario titles)
C  on unit IOUOUT.
C
C  Called by OUTBY1, OUTBY2, OUTBY3, and OUTBY4.
C
C  Input on call:
C
C    parameter list: IOUOUT
C    common blocks:
C    /CITPAR/ SCNAME
C    /PROJEC/ PROJID
C
C  Notes:
C
C  (1) OBYHDR was added for MOBILE4 In-House Version 02 (By Model Year Tables).
C
C
      CHARACTER*4 PROJID,SCNAME
C
      COMMON /CITPAR/ SCNAME(4)
      COMMON /PROJEC/ PROJID(20)
C
      WRITE(IOUOUT,200) PROJID
  200 FORMAT('1',20A4)
C
      WRITE(IOUOUT,210) SCNAME
  210 FORMAT(' ',4A4)
C
      RETURN
      END
