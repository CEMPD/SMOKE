C Version @(#)$Id$ $Source$ $Date$ 

      SUBROUTINE GETHOR(FLGERR,*)
C
C  GETHOR reads in and validates the Control Section flags via the horizontal
C  input mode: all flag values entered on a single line (I1,1X,I1,1X,...,I1).
C
C  Called by CONSEC.
C
C  Calls QUITER.
C
C  Input on call:
C
C    common blocks:
C    /FLAGS1/ PROMPT
C    /IOUCOM/ IOUGEN,IOUASK
C
C  Output on return:
C
C    parameter list: FLGERR
C    common blocks:
C    /FLAGS1/ TAMFLG,SPDFLG,VMFLAG
C    /FLAGS2/ MYMRFG,NEWFLG,IMFLAG,ALHFLG,IMTIER
C    /FLAGS3/ ATPFLG,RLFLAG,LOCFLG,TEMFLG,OUTFMT
C    /FLAGS4/ PRTFLG,IDLFLG,NMHFLG,HCFLAG,COLDFG
C    /PRSCH2/ PRSFLG
C    /PRGCH2/ PRGFLG
C
C  Local variable dictionary:
C
C   Name   Type              Description
C  ------  ----  -------------------------------------------------------
C  USERFG   I    User's input of ATPFLG before it is broken down into
C                three flags: ATPFLG, PRSFLG, and PRGFLG.
C
C
C  Notes:
C
C  OUTFMT's range was extended up to 5 for Version 02's BYMY Tables option.
C  HCFLAG's range was extended up to 6 for MOBILE4.1.
C  COLDFG was added, but commented out for future use.
C  COLDFG is set in BD17.
C  Extended range of SPDFLG to 4 to read in TLVMT via WTLRUL for MOBILE4.1 D02
C  NEWFLG's range was extended up to 6 for MOBILE5.
C  OUTFMT's range was extended up to 6 for version 3 of MOBILE5 and allows the
C  spreadsheet readable/machine readable format option.
C  and it will be extended again to increase the range to 7 to allow a combo
C  spreadsheet/machine readable and by model year.
C  USERFG was added in MOBILE5 version 3 break up ATPFLG.
C
      INTEGER PROMPT,TAMFLG,SPDFLG,VMFLAG,OXYFLG,DSFLAG
      INTEGER ALHFLG,PRSFLG,PRGFLG,PRSPGM,PRSFQT
      INTEGER ATPFLG,TPDFLG,RLFLAG,TEMFLG,OUTFMT
      INTEGER PRTFLG,HCFLAG,COLDFG,PRGPGM,PRGFQT
C
C
      COMMON /FLAGS1/ PROMPT,TAMFLG,SPDFLG,VMFLAG,OXYFLG,DSFLAG
      COMMON /FLAGS2/ MYMRFG,NEWFLG,IMFLAG,ALHFLG,IMTIER
      COMMON /FLAGS3/ ATPFLG,TPDFLG,RLFLAG,LOCFLG,TEMFLG,OUTFMT
      COMMON /FLAGS4/ PRTFLG,IDLFLG,NMHFLG,HCFLAG,COLDFG
      COMMON /IOUCOM/ IOUIMD,IOUGEN,IOUREP,IOUERR,IOUASK
      COMMON /PRSCH2/ PRSPGM,PRSFQT,CRPRS,PRSFLG
      COMMON /PRGCH2/ PRGPGM,PRGFQT,CRPRG,PRGFLG
C
      INTEGER FLGERR,USERFG
C
      IF(PROMPT.EQ.2) WRITE(IOUASK,200)
  200 FORMAT('&','Enter TAMFLG, SPDFLG, VMFLAG,',
     *        ' MYMRFG, NEWFLG, IMFLAG, ALHFLG,'/
     *        ' ATPFLG, RLFLAG, LOCFLG, TEMFLG, OUTFMT,',
     *        ' PRTFLG, IDLFLG, NMHFLG, HCFLAG')
C    *        ' PRTFLG, IDLFLG, NMHFLG, HCFLAG & COLDFG')
C
      READ(IOUGEN,100,ERR=99,END=99)
     *                  TAMFLG, SPDFLG, VMFLAG,
     *          MYMRFG, NEWFLG, IMFLAG, ALHFLG,
     *          USERFG, RLFLAG, LOCFLG, TEMFLG, OUTFMT,
     *          PRTFLG, IDLFLG, NMHFLG, HCFLAG
C    *          PRTFLG, IDLFLG, NMHFLG, HCFLAG, COLDFG
  100 FORMAT(I1,17(1X,I1))
C
      IF(TAMFLG.LT.0.OR.TAMFLG.GT.2) CALL QUITER(1.,TAMFLG,1,FLGERR)
      IF(SPDFLG.LT.1.OR.SPDFLG.GT.4) CALL QUITER(2.,SPDFLG,1,FLGERR)
      IF(VMFLAG.LT.1.OR.VMFLAG.GT.3) CALL QUITER(3.,VMFLAG,1,FLGERR)
      IF(MYMRFG.LT.1.OR.MYMRFG.GT.4) CALL QUITER(4.,MYMRFG,1,FLGERR)
      IF(NEWFLG.LT.1.OR.NEWFLG.GT.6) CALL QUITER(5.,NEWFLG,1,FLGERR)
      IF(IMFLAG.LT.1.OR.IMFLAG.GT.5) CALL QUITER(6.,IMFLAG,1,FLGERR)
      IF(ALHFLG.LT.1.OR.ALHFLG.GT.3) CALL QUITER(7.,ALHFLG,1,FLGERR)
      IF(USERFG.LT.1.OR.USERFG.GT.8) CALL QUITER(8.,USERFG,1,FLGERR)
      IF(RLFLAG.LT.1.OR.RLFLAG.GT.5) CALL QUITER(9.,RLFLAG,1,FLGERR)
      IF(LOCFLG.LT.1.OR.LOCFLG.GT.2) CALL QUITER(10.,LOCFLG,1,FLGERR)
      IF(TEMFLG.LT.1.OR.TEMFLG.GT.2) CALL QUITER(11.,TEMFLG,1,FLGERR)
c
c  LMOS Modification 04/08/93
c  Allow a new output format 7
c
c      IF(OUTFMT.LT.1.OR.OUTFMT.GT.6) CALL QUITER(12.,OUTFMT,1,FLGERR)
      IF(OUTFMT.LT.1.OR.OUTFMT.GT.7) CALL QUITER(12.,OUTFMT,1,FLGERR)
      IF(PRTFLG.LT.1.OR.PRTFLG.GT.4) CALL QUITER(13.,PRTFLG,1,FLGERR)
      IF(IDLFLG.LT.1.OR.IDLFLG.GT.2) CALL QUITER(14.,IDLFLG,1,FLGERR)
      IF(NMHFLG.LT.1.OR.NMHFLG.GT.5) CALL QUITER(15.,NMHFLG,1,FLGERR)
      IF(HCFLAG.LT.1.OR.HCFLAG.GT.6) CALL QUITER(16.,HCFLAG,1,FLGERR)
C     IF(COLDFG.LT.1.OR.COLDFG.GT.3) CALL QUITER(17.,COLDFG,1,FLGERR)
CC
CC    IDLE FLAG FIX
      IDLFLG = 1
C
      COLDFG=2
      IF(NEWFLG.GE.5) THEN
         COLDFG=1
         CALL QUITER(0.,0,135,INERR)
      ENDIF
C
C  Break USERFG into three flags
c
      ATPFLG=1
      PRSFLG=1
      PRGFLG=1
      GOTO (10,20,30,40,50,60,70,80) USERFG
   10 RETURN
C
   20        ATPFLG=2
      RETURN
C
   30        PRSFLG=2
      RETURN
C
   40        PRGFLG=2
      RETURN
C
   50        ATPFLG=2
             PRSFLG=2
      RETURN
C
   60        ATPFLG=2
             PRGFLG=2
      RETURN
C
   70        PRSFLG=2
             PRGFLG=2
      RETURN
C
   80        ATPFLG=2
             PRSFLG=2
             PRGFLG=2
      RETURN
C
C  Error or EOF on any attempted read => take alternate
C  return 1 => run aborts.
C
   99 CALL QUITER(0.,0,125,INERR)
C
      RETURN1
      END
