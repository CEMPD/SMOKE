C Version @(#)$Id$ $Source$ $Date$ 

      FUNCTION HDCPIC(MY,IV)
C
C  HDCPIC picks the HDV gpbhp to gpm units conversion factor.
C
C  Called by HDCGPM.
C
C  Input on call:
C
C    parameter list: MY,IV
C    common blocks:
C    /HDCCOM/ HDCFAC,LOWYR,MAXCIG
C
C  Output on return:
C
C    function: HDCPIC
C
C  Notes:
C
C  None.
C
C
      COMMON /HDCCOM/ HDCFAC(37,2),LOWYR,MAXCIG
C
C  There is a different conversion factor for each model year from
C  LOWYR up through MAXCIG's year.  Below LOWYR, use LOWYR's.
C  Above MAXCIG's year, use MAXCIG's year's.
C
      IG=MY-LOWYR+1
      IF(IG.LT.1) IG=1
      IF(IG.GT.MAXCIG) IG=MAXCIG
C
C  IV must be a HDV code (4 or 7).
C
      IF(IV.EQ.4) IVHD=1
      IF(IV.EQ.7) IVHD=2
C
      HDCPIC=HDCFAC(IG,IVHD)
C
      RETURN
      END
