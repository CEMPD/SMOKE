C Version "@(#)$Id$ $Source$ $Date$ 

      SUBROUTINE GETVER(FLGERR,*)
C
C  GETVER reads in and validates the Control Section flags via
C  the vertical input mode: one flag per input line in I1 format.
C
C  Called by CONSEC.
C
C  Calls QUITER.
C
C  Input on call:
C
C    common blocks:
C    /FLAGS1/ PROMPT
C    /FLAGS3/ OUTFMT
C    /NAMES1/ NAMFLG
C    /STRING/ COLON
C
C  Output on return:
C
C    parameter list: FLGERR
C    common blocks:
C    /FLAGS1/ TAMFLG,SPDFLG,VMFLAG
C    /FLAGS2/ MYMRFG,NEWFLG,IMFLAG,ALHFLG
C    /FLAGS3/ ATPFLG,RLFLAG,LOCFLG,TEMFLG,OUTFMT
C    /FLAGS4/ PRTFLG,IDLFLG,NMHFLG,HCFLAG,COLDFG
C    /NFTPDT/ NONFTP
C    /PRGCH2/ PRGFLG
C    /PRSCH2/ PRSFLG
C
C  Local variable dictionary:
C
C   Name   Type              Description
C  ------  ----  -----------------------------------------------
C  USERFG   I    User's input of ATPFLG before it is broken down
C                into three flags: ATPFLG, PRSFLG, and PRGFLG.
C  USERFM   I    User's input of OUTFMT before is determined
C                that OUTFMT should be reassigned to CSV.
C  TMPOUT   I    Temporary storage for OUTFMT
C
C
C  Notes:
C
C  For MOBILE4.1. OUTFMT's range was extended up to 5 for the by
C  model year option.  HCFLAG's range was extended up to 6.
C  COLDFG was added and is set internally.  Extended range of
C  SPDFLG to 4 to read in TLVMT via WTLRUL.
C  For MOBILE5 NEWFLG's range was extended up to 6. OUTFMT's
C  range was extended up to 6 and allows the spreadsheet
C  readable/machine readable format option.  USERFG was
C  added to break up ATPFLG.  GETVER was modified Feb 24, 1993,
C  because the CSV code (to set OUTMFT=6) was moved to the
C  DRIVER.  OUTFMT's range is limited to 6 for user input.
C  June-14-1993 @ ARC-bsg Subtask 222
C 
C  TMM@CSC REQUEST 2-392, comment out code that requires the
C     file extension be .CSV when csv type output is selected.
C  Aug-18-1994 @ CSC-bsg Request 462 extended range of PRTFLG and made common
C  for non-FTP indicators when PRTFLG is 5 or 6
C
C  April-1995 @ CSC-zp Req. 2-495 "Add MOB5a-H to MOBILE5b"
C  Changed the range of the valid values for  IMFLAG
C
C  Jun-20-1996 @DynTel-yc Request 2-624  Non-FTP indicators(PRTFLG=5 or 6) 
C               have been disabled to remove non-FTP effects in MOBILE5b
C
      INCLUDE 'FLAGS1.I'
      INCLUDE 'FLAGS2.I'
      INCLUDE 'FLAGS3.I'
      INCLUDE 'FLAGS4.I'
      INCLUDE 'IOUCOM.I'
      INCLUDE 'NAMES1.I'
      INCLUDE 'NFTPDT.I'
      INCLUDE 'PRSCH2.I'
      INCLUDE 'PRGCH2.I'
      INCLUDE 'STRING.I' 
C
      INTEGER FLGERR,USERFG,USERFM
C
C  Do the in-range check immediately after flag read, since if
C  in interactive mode (PROMPT.EQ.2), user can then receive
C  out-of-range message immediately after entering an erroneous
C  flag value.
C
      IF(PROMPT.EQ.2) WRITE(IOUASK,200)NAMFLG(1),COLON
  200 FORMAT('&','Enter',T7,A8,T14,A1)
      READ(IOUGEN,100,ERR=99,END=99) TAMFLG
  100 FORMAT(I1)
      IF(TAMFLG.LT.0.OR.TAMFLG.GT.2) CALL QUITER(1.,TAMFLG,1,FLGERR)
C
      IF(PROMPT.EQ.2) WRITE(IOUASK,200)NAMFLG(2),COLON
      READ(IOUGEN,100,ERR=99,END=99) SPDFLG
      IF(SPDFLG.LT.1.OR.SPDFLG.GT.4) CALL QUITER(2.,SPDFLG,1,FLGERR)
C
      IF(PROMPT.EQ.2) WRITE(IOUASK,200)NAMFLG(3),COLON
      READ(IOUGEN,100,ERR=99,END=99) VMFLAG
      IF(VMFLAG.LT.1.OR.VMFLAG.GT.3) CALL QUITER(3.,VMFLAG,1,FLGERR)
C
      IF(PROMPT.EQ.2) WRITE(IOUASK,200)NAMFLG(4),COLON
      READ(IOUGEN,100,ERR=99,END=99) MYMRFG
      IF(MYMRFG.LT.1.OR.MYMRFG.GT.4) CALL QUITER(4.,MYMRFG,1,FLGERR)
C
      IF(PROMPT.EQ.2) WRITE(IOUASK,200)NAMFLG(5),COLON
      READ(IOUGEN,100,ERR=99,END=99) NEWFLG
      IF(NEWFLG.LT.1.OR.NEWFLG.GT.6) CALL QUITER(5.,NEWFLG,1,FLGERR)
C
      IF(PROMPT.EQ.2) WRITE(IOUASK,200)NAMFLG(6),COLON
C
      READ(IOUGEN,100,ERR=99,END=99) IMFLAG
C
      IF(IMFLAG.LT.1.OR.IMFLAG.GT.6) CALL QUITER(6.,IMFLAG,1,FLGERR)
C
      IF(PROMPT.EQ.2) WRITE(IOUASK,200)NAMFLG(7),COLON
      READ(IOUGEN,100,ERR=99,END=99) ALHFLG
      IF(ALHFLG.LT.1.OR.ALHFLG.GT.3) CALL QUITER(7.,ALHFLG,1,FLGERR)
C
      IF(PROMPT.EQ.2) WRITE(IOUASK,200)NAMFLG(8),COLON
      READ(IOUGEN,100,ERR=99,END=99) USERFG
      IF(USERFG.LT.1.OR.USERFG.GT.8) CALL QUITER(8.,USERFG,1,FLGERR)
C
      IF(PROMPT.EQ.2) WRITE(IOUASK,200)NAMFLG(9),COLON
      READ(IOUGEN,100,ERR=99,END=99) RLFLAG
      IF(RLFLAG.LT.1.OR.RLFLAG.GT.5) CALL QUITER(9.,RLFLAG,1,FLGERR)
C
      IF(PROMPT.EQ.2) WRITE(IOUASK,200)NAMFLG(10),COLON
      READ(IOUGEN,100,ERR=99,END=99) LOCFLG
      IF(LOCFLG.LT.1.OR.LOCFLG.GT.2) CALL QUITER(10.,LOCFLG,1,FLGERR)
C
      IF(PROMPT.EQ.2) WRITE(IOUASK,200)NAMFLG(11),COLON
      READ(IOUGEN,100,ERR=99,END=99) TEMFLG
      IF(TEMFLG.LT.1.OR.TEMFLG.GT.2) CALL QUITER(11.,TEMFLG,1,FLGERR)
C
      IF(PROMPT.EQ.2) WRITE(IOUASK,200)NAMFLG(12),COLON
      READ(IOUGEN,100,ERR=99,END=99) USERFM
      IF(USERFM.LT.1.OR.USERFM.GT.8) CALL QUITER(12.,USERFM,1,FLGERR)
C
C  Reassign OUTFMT to USERFM if CSV is not selected
C
C     TMM@CSC Request 2-392, comment out code that requires
C         that the spreadsheet output file end with .CSV
C
CC    IF(OUTFMT.EQ.7.AND.USERFM.EQ.6) TMPOUT=USERFM
CC    IF(OUTFMT.EQ.7.AND.USERFM.NE.6) TMPOUT=OUTFMT
C
C  Error because USERFM=6 is not allowed in a non-CSV format
C
CC    IF(OUTFMT.NE.7.AND.USERFM.EQ.6)
CC   *    CALL QUITER(12.,USERFM,74,FLGERR)
CC    IF(OUTFMT.NE.7.AND.USERFM.NE.6) TMPOUT=USERFM
C
C  Error because USERFM=6 is not allowed in a non-CSV format
C
      OUTFMT=USERFM
        IF(OUTFMT.EQ.7) CALL HDRCSV
C
C  Aug-18-1994 @ CSC-bsg Request 462 extended range of PRTFLG and made common
C  for non-FTP indicators when PRTFLG is 5 or 6
C
      IF(PROMPT.EQ.2) WRITE(IOUASK,200)NAMFLG(13),COLON
      READ(IOUGEN,100,ERR=99,END=99) PRTFLG
      IF(PRTFLG.LT.1.OR.PRTFLG.GT.4) CALL QUITER(13.,PRTFLG,1,FLGERR)
C
C  Jun-20-1996 @DynTel-yc Request 2-624  Comment out code related to 
C               non-FTP effects in conjunction with PRTFLG=5 or 6.
C
CC      IF(PRTFLG.LT.5) THEN
        NONFTP=1
CC      ELSEIF(PRTFLG.EQ.5) THEN
CC        NONFTP=2
CC      ELSE
CC        NONFTP=3
CC      ENDIF
C
      IF(PROMPT.EQ.2) WRITE(IOUASK,200)NAMFLG(14),COLON
      READ(IOUGEN,100,ERR=99,END=99) IDLFLG
      IF(IDLFLG.LT.1.OR.IDLFLG.GT.2) CALL QUITER(14.,IDLFLG,1,FLGERR)
C
      IF(PROMPT.EQ.2) WRITE(IOUASK,200)NAMFLG(15),COLON
      READ(IOUGEN,100,ERR=99,END=99) NMHFLG
      IF(NMHFLG.LT.1.OR.NMHFLG.GT.5) CALL QUITER(15.,NMHFLG,1,FLGERR)
C
      IF(PROMPT.EQ.2) WRITE(IOUASK,200)NAMFLG(16),COLON
      READ(IOUGEN,100,ERR=99,END=99) HCFLAG
      IF(HCFLAG.LT.1.OR.HCFLAG.GT.6) CALL QUITER(16.,HCFLAG,1,FLGERR)
C
      COLDFG=2
      IF(NEWFLG.GE.5) THEN
         COLDFG=1
         CALL QUITER(0.,0,135,INERR)
      ENDIF
C
C  Break USERFG into three flags
C
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
   99 CALL QUITER(0.,0,124,INERR)
C
      RETURN1
      END
