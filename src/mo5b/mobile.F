C Version "@(#)$Id$ $Source$ $Date$ 

      SUBROUTINE MOBILE(INERR, *)
C
C  Called by user's program (example: Air Quality Analysis
C  System).
C
C  Calls ADJUST, BDSAVE, CONSEC, EFCALX, ONESEC, OUTNEW, OUTPUT,
C        PARSEC, QUITER, REINIT, RESTOR, REGMOD and SAVER.
C
C  The following variables are passed in argument lists or common
C  blocks and in MOBILE:
C
C    argument lists: ICY,INERR
C    common blocks:
C    /EVAL/   MEVAL
C    /FLAGS3/ OUTFMT
C    /IOUCOM/ IOUERR
C    /REGION/ IREJN
C    /SAVE01/ JCALL
C    /YEARS4/ IYEND
C
C  Output on return:
C
C    parameter list: INERR
C
C  Local variable dictionary:
C
C   Name   Type              Description
C  ------  ----  -------------------------------------------------------
C  ICYOLD   I    calendar year from previous pass of scenario loop
C  IROLD    I    region from previous pass of scenario loop
C                      0 = message has not yet been given
C                      1 = message has been given: not again!
C  IY       I    loop counting variable for interpolation runs
C  LEAST1   I    flag: 0 = OUTPUT has not yet been called.
C                      1 = OUTPUT has been called at least 1 time.
C  MEVOLD   I    calendar month from previous pass of scenario loop
C  NIY      I    number of year interations for by-month interpolation
C  NIYOLD   I    number of year interations from previous pass
C                of scenario loop
C
C
C  Notes:
C
C  For MOBILE4.1 the subroutine MOBILE was created from MAIN.
C  The /BYMYCd/ CBs were added.
C  May-18-1993 @ ARC-bsg Subtask 238 the common block IMPAR8
C  was removed, it contained the variable ICALC
C  June-15-1993 @ ARC-bk Subtask 244 the common block IMPAR5
C  was modified, the HDGV I/M credit was expanded to include
C  an initial value for NOx.
C  24-August-1994 @CSC-pme request 2-446 Include file BASE12.I was removed
C  (24-May-1996) @DynTel-yc Request 2-621 The calendar year was expanded
C                 from year 2020 to year 2051
C
      INCLUDE 'EVAL.I'
      INCLUDE 'FLAGS2.I'
      INCLUDE 'FLAGS3.I'
      INCLUDE 'FLAGS4.I'
      INCLUDE 'IOUCOM.I'
      INCLUDE 'REGION.I'
      INCLUDE 'SAVE01.I'
      INCLUDE 'YEARS4.I'
C
      IF(JCALL.EQ.0) GOTO 10
C
C  Subroutine call, setup run.
C
      CALL REINIT
      IF(JCALL.EQ.1) CALL BDSAVE
      IF(JCALL.EQ.2) CALL RESTOR
      JCALL=2
C
C  Initialize non-COMMON variables (INERR & local variables).
C
   10 INERR=0
      ICYOLD=-1
      IROLD=-1
      MEVOLD=-1
      NIYOLD=-1
      LEAST1=0
C
C  Control and One-Time Data Sections:
C
C  Get the execution control flags and parameters, and then use them to read
C  in data, if any, that is read in once and applied to all scenarios.  Note:
C  alternate RETURN1 is the fatal read error exit: execution stops.
C
      IF(IOUASK.EQ.0.OR.IOUASK.EQ.9) WRITE(IOUASK,110)
  110 FORMAT(' Reading information.')
C
      CALL CONSEC(*99,*999)
      CALL ONESEC(INERR,*99)
C
C  ljn hc=3 3/5/93 Warning MSG in QUITER
C
      IF(HCFLAG.EQ.3.AND.OUTFMT.NE.3.AND.OUTFMT.NE.5)
     *                    CALL QUITER(0.,0,153,INERR)
C
C  only get error diagnostics out of the run.
C
   20 CALL PARSEC(ICY,INERR,*98)
C
C  50+ errors is (arbritrary) run limit => stop execution.  Otherwise,
C  one or more input processing errors => abort current scenario & then do
C  input processing for error diagnostics on the remaining scenarios.  INERR
C  is never reset to 0, but instead accumulates until the run ends or the
C  cutoff value (50) is reached.
C
      IF(INERR.GT.50) GOTO 30
      IF(INERR.GT.0) GOTO 20
C
      IF(IOUASK.EQ.0.OR.IOUASK.EQ.9) WRITE(IOUASK,115)
  115 FORMAT(' Performing calculations.')
C
C  For interpolated (non Jan 1st) emissions, loop through ICY and ICY+1,
C  save intermediate results, then interpolate to the evaluation month.
C
      NIY=2
      IF(MEVAL.EQ.1.OR.ICY.EQ.IYEND) NIY=1
C
      DO 35 IY=1,NIY
      ICY1=ICY+IY-1
      IF(NIY.GT.1 .OR.
     *   NIYOLD.NE.NIY .OR.
     *   ICYOLD.NE.ICY1 .OR.
     *   IROLD.NE.IREJN .OR.
     *   MEVOLD.NE.MEVAL) CALL REGMOD(ICY1,INERR)
      IF(INERR.GT.50) GOTO 30
      IF(INERR.GT.0) GOTO 20
C
      ICYOLD=ICY1
      IROLD=IREJN
      MEVOLD=MEVAL
      NIYOLD=NIY
C
      CALL EFCALX(ICY1,INERR,IY)
      IF(IDLFLG.EQ.2) CALL IDLCAL(ICY)
      IF(INERR.GT.50) GOTO 30
      IF(INERR.GT.0) GOTO 20
C
      IF(NIY.GT.1) CALL SAVER(IY)
C
  35  CONTINUE
C
      IF(NIY.GT.1) CALL ADJUST
C
      IF(IOUASK.EQ.0.OR.IOUASK.EQ.9) WRITE(IOUASK,120)
  120 FORMAT(' Preparing output.')
C
      CALL OUTPUT(ICY)
      LEAST1=1
      GOTO 20
C
   30 CALL QUITER(0.0,0,28,INERR)
C
C  If user is supplying replacement basic FTP ef parameter records, always
C  print the records contents and disposition table, even when input errors
C  have prevented the normal output routines from being called.  The only
C  exception is when a READ error / end-of-file occurs in CONSEC or ONESEC.
C
C  This section can be expanded to print other optional input.
C
   98 IF(LEAST1.EQ.0.AND.
     *  (NEWFLG.EQ.2.OR.NEWFLG.EQ.4.OR.NEWFLG.EQ.6)) CALL OUTNEW
C
   99 RETURN
C
C  End-of-file on 1st read of call => input for call missing entirely.
C  Calling program can branch to attach a new input file or terminate
C  calls to MOBILE.
C
  999 RETURN 1
      END
