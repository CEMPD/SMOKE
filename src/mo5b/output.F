C Version "@(#)$Id$ $Source$ $Date$ 

      SUBROUTINE OUTPUT(ICY)
C
C  OUTPUT writes scenario results on unit IOUREP.
C
C  Called by MOBILE.
C
C  Calls OUTZTT, OUTCSV, OUTDT1, OUTDT2, OUTDT3, OUTDT4,
C        OUTHD1, OUTHD2, OUTHD3 and OUTHD4.
C
C  Input on call:
C
C    /LEVBLK/ LEVFLG
C
C    parameter list: ICY
C    common blocks:
C    /EVAL/   MEVAL,MONTH
C    /FLAGS1/ PPTFLG 
C    /FLAGS3/ OUTFMT
C    /IOUCOM/ IOUREP
C    /LEVINF/ LEVNAM
C    /REGION/ INITPR
C
C
C  Notes:
C
C  For MOBILE4.1 by model year was added.
C  For MOBILE5 OUTFMT=6 for spreadsheets was added.
C  Feb 24,1993 changes to OUTPUT were needed to make the one
C  record spreadsheet version.
C  June-14-1993 @ ARC-bsg Subtask 222
C  August 11, 1994 @ CSC-pme request 2-455. The subroutine was modified
C  to print the name of the data file that contains the LEV phase-in
C  schedule and the ZML and DET for LEV standards.
C  August 30, 1994 @ CSC-pme request 2-455. Extra code that checked outdt5
C  value when OUTFMT = 3 was removed.
C  October-4-1994 @CSC-pme request 2-472. Checks for APPLEV flag were added
C  to identify the type of messages that need to be printed whether or not
C  the LEV standards are applied to LDGT2. Also it is checked to identify
C  if the ZEV TF percentages should be calculated.
C  Sep-13-1994 @ bsg-CSC Request 2-000, 2-392 removed call to OUTSIX
C  Jan-17-1995 @CSC-ZP  Request 2-392. When OUTFMT flag equal six,
C  subroutine OUTHD1 will be called without checking INITPR flag.
C
      INCLUDE 'EVAL.I'
      INCLUDE 'FLAGS1.I'
      INCLUDE 'FLAGS3.I'
      INCLUDE 'IOUCOM.I'
      INCLUDE 'LEVBLK.I'
      INCLUDE 'LEVINF.I'
      INCLUDE 'MYRCAL.I'
      INCLUDE 'REGION.I'
      INCLUDE 'YEARS4.I'
C
C     Local Variables
C
C     now call outztt, that calculates the ZEV TF percentage
C
      IF(LEVFLG.EQ.2) CALL OUTZTT(ICY)
C
C  Select report format.
C
      GOTO(10,20,30,40,50,10,60),OUTFMT
C
C  OUTFMT=1 or 6
C
   10        IF(INITPR.EQ.1 .OR. OUTFMT.EQ.6)  CALL OUTHD1
C
      WRITE(IOUREP,210) MONTH(MEVAL)
  210 FORMAT('0','Emission factors are as of ',A4,
     *          ' 1st of the indicated calendar year.')
C
      IF(LEVFLG.EQ.2.AND.ICY.GE.IY1994) THEN
        IF(APPLEV.NE.2)THEN
          IF (IMFLG.EQ.1) WRITE(IOUREP,220) LEVYRS
  220     FORMAT(' ',
     *   'LEV phase-in begins in ',I4,' without using ',
     *   '(4/8/94) Guidance Memo Credits')
          IF (IMFLG.EQ.2) WRITE(IOUREP,230) LEVYRS
  230    FORMAT(' ',
     *  'LEV phase-in begins in ',I4,' using ',
     *  '(4/8/94) Guidance Memo Credits')
        ELSE
          IF (IMFLG.EQ.1) WRITE(IOUREP,240) LEVYRS
  240     FORMAT(' ',
     *    'LEV phase-in begins in ',I4,' including LDGT2 and ',
     *    'without using ',
     *   '(4/8/94) Guidance Memo Credits')
        IF (IMFLG.EQ.2) WRITE(IOUREP,250) LEVYRS
  250     FORMAT(' ',
     *    'LEV phase-in begins in ',I4,' including LDGT2 and ',
     *    'using (4/8/94) Guidance Memo Credits')
        ENDIF
      ENDIF
C
      CALL OUTDT1(ICY)
      GOTO 99
C
C  OUTFMT=2
C
   20 IF(INITPR.EQ.1) CALL OUTHD2
      WRITE(IOUREP,210) MONTH(MEVAL)
      IF(LEVFLG.EQ.2.AND.ICY.GE.IY1994) THEN
        IF(APPLEV.NE.2)THEN
          IF (IMFLG.EQ.1) WRITE(IOUREP,220) LEVYRS
          IF (IMFLG.EQ.2) WRITE(IOUREP,230) LEVYRS
        ELSE
          IF (IMFLG.EQ.1) WRITE(IOUREP,240) LEVYRS
          IF (IMFLG.EQ.2) WRITE(IOUREP,250) LEVYRS
        ENDIF
      ENDIF
C
      CALL OUTDT2(ICY)
      GOTO 99
C
C  OUTFMT=3
C
   30 IF(INITPR.EQ.1) CALL OUTHD3(ICY)
      WRITE(IOUREP,210) MONTH(MEVAL)
      IF(LEVFLG.EQ.2.AND.ICY.GE.IY1994) THEN
        IF(APPLEV.NE.2)THEN
          IF (IMFLG.EQ.1) WRITE(IOUREP,220) LEVYRS
          IF (IMFLG.EQ.2) WRITE(IOUREP,230) LEVYRS
        ELSE
          IF (IMFLG.EQ.1) WRITE(IOUREP,240) LEVYRS
          IF (IMFLG.EQ.2) WRITE(IOUREP,250) LEVYRS
        ENDIF
      ENDIF
C
C     Print the name of file with LEV phase-in data.
C
      IF(PPTFLG.EQ.5.OR.PPTFLG.EQ.6)THEN
        WRITE(IOUREP,200) LEVNAM
  200   FORMAT('0','LEV phase-in data read from file: ',A40)
      ENDIF
C
      CALL OUTDT3(ICY)
      GOTO 99
C
C  OUTFMT=4
C
   40 IF(INITPR.EQ.1) CALL OUTHD4(ICY)
      WRITE(IOUREP,210) MONTH(MEVAL)
      IF(LEVFLG.EQ.2.AND.ICY.GE.IY1994) THEN
        IF(APPLEV.NE.2)THEN
          IF (IMFLG.EQ.1) WRITE(IOUREP,220) LEVYRS
          IF (IMFLG.EQ.2) WRITE(IOUREP,230) LEVYRS
        ELSE
          IF (IMFLG.EQ.1) WRITE(IOUREP,240) LEVYRS
          IF (IMFLG.EQ.2) WRITE(IOUREP,250) LEVYRS
        ENDIF
      ENDIF
C
C     Print the name of file with LEV phase-in data.
C
      IF(PPTFLG.EQ.5.OR.PPTFLG.EQ.6)THEN
       WRITE(IOUREP,200) LEVNAM
      ENDIF
C
      CALL OUTDT4(ICY)
      GOTO 99
C
C  OUTFMT=5
C
   50 IF(INITPR.EQ.1) CALL OUTHD3(ICY)
      WRITE(IOUREP,210) MONTH(MEVAL)
C
      IF(LEVFLG.EQ.2.AND.ICY.GE.IY1994) THEN
        IF(APPLEV.NE.2)THEN
          IF (IMFLG.EQ.1) WRITE(IOUREP,220) LEVYRS
          IF (IMFLG.EQ.2) WRITE(IOUREP,230) LEVYRS
        ELSE
          IF (IMFLG.EQ.1) WRITE(IOUREP,240) LEVYRS
          IF (IMFLG.EQ.2) WRITE(IOUREP,250) LEVYRS
        ENDIF
      ENDIF
C
C     Print the name of file with LEV phase-in data.
C
      IF(PPTFLG.EQ.5.OR.PPTFLG.EQ.6)THEN
        WRITE(IOUREP,200) LEVNAM
      ENDIF
C
      CALL OUTDT3(ICY)
      CALL OUTDT5(ICY)
      GOTO 99
C
C  OUTFMT=7  (Temporary long version CSV)
C
   60 CALL OUTCSV(ICY)
   99 RETURN
C
      END
