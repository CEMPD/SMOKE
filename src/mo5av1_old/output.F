C Version @(#)$Id$ $Source$ $Date$ 

      SUBROUTINE OUTPUT(ICY)
C
C  OUTPUT writes scenario results on unit IOUREP.
C
C  Called by MOBILE.
C
C  Calls OUTZTT, OUTDT1, OUTDT2, OUTDT3, OUTDT4,
C        OUTHD1, OUTHD2, OUTHD3 and OUTHD4.
C
C  Input on call:
C
C    /LEVBLK/ LEVFLG
C
C    parameter list: ICY
C    common blocks:
C    /EVAL/   MEVAL,MONTH
C    /FLAGS3/ OUTFMT
C    /IOUCOM/ IOUREP
C    /REGION/ INITPR
C
C  Notes:
C
C  By model year table(s) generation routing and print routine call were added
C  for Version 02.  In MOBILE5 OUTFMT=6 for spreadsheets was added.
C
C
      CHARACTER*4 MONTH
      INTEGER ATPFLG,TPDFLG,RLFLAG,TEMFLG,OUTFMT
C
      COMMON /EVAL/ MEVAL,MONTH(12)
      COMMON /FLAGS3/ ATPFLG,TPDFLG,RLFLAG,LOCFLG,TEMFLG,OUTFMT
      COMMON /IOUCOM/ IOUIMD,IOUGEN,IOUREP,IOUERR,IOUASK
      COMMON /LEVBLK/ LEVFLG,IMFLG,LEVYRS
      COMMON /MYRCAL/ XMYM(25,8),JANMYR(25,8),TF(25,8),TFMYM(25,8)
      COMMON /REGION/ FEET(2),IREJN,ALT,INITPR
      COMMON /YEARS4/ IY1941,IY1960,IY1994,IY2020
C
C     Local Variables
C
C     now call outztt, that calculates the ZEV TF percentage
C
      IF(LEVFLG.EQ.2) CALL OUTZTT(ICY)
C
C  Select report format.
C
c  LMOS Modification 04/08/93
c
c  Add output format 7
c
c      GOTO(10,20,30,40,50,60),OUTFMT
      GOTO(10,20,30,40,50,60,70),OUTFMT
c
C
C  OUTFMT=1
C
   10 IF(INITPR.EQ.1) CALL OUTHD1
      WRITE(IOUREP,100) MONTH(MEVAL)
  100 FORMAT('0','Emission factors are as of ',A4,
     *          ' 1st of the indicated calendar year.')
C
      IF(LEVFLG.EQ.2.AND.ICY.GE.IY1994) THEN
        IF (IMFLG.EQ.2) WRITE(IOUREP,205) LEVYRS
  205   FORMAT(' ',
     * 'LEV phase-in begins in ',I4,' using ',
     * '(12/1/92) Guidance Memo Credits')
        IF (IMFLG.EQ.1) WRITE(IOUREP,210) LEVYRS
  210   FORMAT(' ',
     * 'LEV phase-in begins in ',I4,' *NOT* using ',
     * '(12/1/92) Guidance Memo Credits')
      ENDIF
C
      CALL OUTDT1(ICY)
      GOTO 99
C
C  OUTFMT=2
C
   20 IF(INITPR.EQ.1) CALL OUTHD2
      WRITE(IOUREP,100) MONTH(MEVAL)
      IF(LEVFLG.EQ.2.AND.ICY.GE.IY1994) THEN
        IF(IMFLG.EQ.2) WRITE(IOUREP,205)LEVYRS
        IF(IMFLG.EQ.1) WRITE(IOUREP,210)LEVYRS
      ENDIF
C
      CALL OUTDT2(ICY)
      GOTO 99
C
C  OUTFMT=3
C
   30 IF(INITPR.EQ.1) CALL OUTHD3(ICY)
      WRITE(IOUREP,100) MONTH(MEVAL)
      IF(LEVFLG.EQ.2.AND.ICY.GE.IY1994) THEN
        IF(IMFLG.EQ.2) WRITE(IOUREP,205) LEVYRS
        IF(IMFLG.EQ.1) WRITE(IOUREP,210) LEVYRS
      ENDIF
C
      CALL OUTDT3(ICY)
        IF(OUTFMT.EQ.5) THEN
          CALL OUTDT5(ICY)
          IF(LEVFLG.EQ.2.AND.ICY.GE.IY1994) THEN
            IF(IMFLG.EQ.2) WRITE(IOUREP,205) LEVYRS
            IF(IMFLG.EQ.1) WRITE(IOUREP,210) LEVYRS
          ENDIF
        ENDIF
      GOTO 99
C
C  OUTFMT=4
C
   40 IF(INITPR.EQ.1) CALL OUTHD4(ICY)
      WRITE(IOUREP,100) MONTH(MEVAL)
      IF(LEVFLG.EQ.2.AND.ICY.GE.IY1994) THEN
        IF(IMFLG.EQ.2) WRITE(IOUREP,205)LEVYRS
        IF(IMFLG.EQ.1) WRITE(IOUREP,210)LEVYRS
      ENDIF
C
      CALL OUTDT4(ICY)
      GOTO 99
C
C  OUTFMT=5
C
   50 IF(INITPR.EQ.1) CALL OUTHD3(ICY)
      WRITE(IOUREP,100) MONTH(MEVAL)
CCCC  CALL OUTDT3(ICY)
      IF(LEVFLG.EQ.2.AND.ICY.GE.IY1994) THEN
        IF(IMFLG.EQ.2) WRITE(IOUREP,205)LEVYRS
        IF(IMFLG.EQ.1) WRITE(IOUREP,210)LEVYRS
      ENDIF
      CALL OUTDT3(ICY)
      CALL OUTDT5(ICY)
      GOTO 99
C
C  OUTFMT=6
C
   60 CALL OUTHD1
      WRITE(IOUREP,110) MONTH(MEVAL)
  110 FORMAT('"Emission factors are as of","","",',A4,
     *          ',"1st of the indicated calendar year."')
      IF(LEVFLG.EQ.2.AND.ICY.GE.IY1994) THEN
        IF(IMFLG.EQ.2) WRITE(IOUREP,205)LEVYRS
        IF(IMFLG.EQ.1) WRITE(IOUREP,210)LEVYRS
      ENDIF
C
      CALL OUTDT1(ICY)
C
c
c  LMOS Modification 04/08/93
c
c  Add output format 7, OUTFMT=7
c
      GOTO 99
   70 CALL OUTDT7(ICY)
c
C
   99 RETURN
      END
