C Version @(#)$Id$ $Source$ $Date$ 

      SUBROUTINE ONESEC(INERR,*)
C
C  ONESEC checks if control flags indicate One-time Data. It then reads in,
C  validates, and stores the data.
C
C  Called by MOBILE.
C
C  Calls CALCUM, CEVTPD, GETATP, GETBYM, GETCUM, GETERQ, GETIMP, GETIM2,
C        GETIMC, GETMRN, GETNAT, GETTAM, GETVMT, RDCITY, SAVEPB, SKIPIN,
C        CCKFUN and WTLRUL.
C
C  Input on call:
C
C    parameter list: INERR
C    common blocks:
C    /FLAGS1/ TAMFLG,VMFLAG
C    /FLAGS2/ NEWFLG,IMFLAG
C    /FLAGS3/ ATPFLG,LOCFLG,OUTFMT
C    /FLAGS4/ PRTFLG
C    /PRSCH2/ PRSFLG
C    /PRGCH2/ PRGCH2
C    /TAMPB1/ LPOD
C
C  Output on call:
C
C    parameter list: INERR
C
C  Notes:
C
C  A call to GETBYM (if user asks for BYMY tables) was added for Version 02.
C  Added call to WTLRUL (SPDFLG=4) in MOBILE4.1, D02.
C  ONESEC was modified in MOBILE4.1v6 to call CALCUM and SAVEPB.
C  Modified in MOBILE5.3 to decide if Pressure/Purge checks are made.
C  Modified in MOBILE5.3 to call CHKPGM a replacement subroutine for
C  checking if IM ATP Pressure and/or Purge program parameters are not
C  all equal if applicable.  The old subroutine was CHCKIM.  CHKPGM is
C  called unconditionally.  The old subroutine was only called if
C  there was at least two I/M programs or an ATP and an I/M.
C
      INTEGER PROMPT,TAMFLG,SPDFLG,VMFLAG,OXYFLG,DSFLAG
      INTEGER ALHFLG,PRSFLG,PRGFLG,PRSPGM,PRSFQT
      INTEGER ATPFLG,TPDFLG,RLFLAG,TEMFLG,OUTFMT
      INTEGER PRTFLG,HCFLAG,COLDFG,PRGPGM,PRGFQT
C
      COMMON /FLAGS1/ PROMPT,TAMFLG,SPDFLG,VMFLAG,OXYFLG,DSFLAG
      COMMON /FLAGS2/ MYMRFG,NEWFLG,IMFLAG,ALHFLG,IMTIER
      COMMON /FLAGS3/ ATPFLG,TPDFLG,RLFLAG,LOCFLG,TEMFLG,OUTFMT
      COMMON /FLAGS4/ PRTFLG,IDLFLG,NMHFLG,HCFLAG,COLDFG
      COMMON /IOUCOM/ IOUIMD,IOUGEN,IOUREP,IOUERR,IOUASK
      COMMON /PRSCH2/ PRSPGM,PRSFQT,CRPRS,PRSFLG
      COMMON /PRGCH2/ PRGPGM,PRGFQT,CRPRG,PRGFLG
      COMMON /TAMPB1/ LPOD,PBBTR(3,2,25,4)
C
C  Look for new tampering rates.
C
   10 IF(TAMFLG.EQ.2) CALL GETTAM(INERR,*99)
C
C  Check for user supplied vmt mix to be used for all scenarios of this run.
C
      IF(VMFLAG.EQ.3) CALL GETVMT(INERR,*99)
C
C  Optionally read in and screen model year mileage accrual rates and / or
C  registration distributions.  In all cases, normalize the distributions.
C
      CALL GETMRN(INERR,*99)
C
C  If requested, read, screen and store emission factor modification data.
C
      IF(NEWFLG.EQ.2.OR.NEWFLG.EQ.4.OR.NEWFLG.EQ.6)
     *    CALL GETERQ(INERR,*99)
C
C  If requested, read, screen and store Evap Test Procedure (ETP) data.
C
      IF(NEWFLG.EQ.3.OR.NEWFLG.EQ.4) CALL GETETP(INERR,*99)
C
C  If requested, read, screen and store I/M data.
C
      IF(IMFLAG.NE.1) CALL GETIMP(INERR,*99)
      IF(IMFLAG.EQ.3.OR.IMFLAG.EQ.5) CALL GETIM2(INERR,*99)
      IF(IMFLAG.NE.1) CALL GETIMC(INERR,*99)
C
C  Anti-tampering program => get parameters and effectiveness rates.
C
      IF(ATPFLG.EQ.2) CALL GETATP(INERR,*99)
C
C  Pressure Functional Inspections program
C  get parameters and effectiveness rates.
C
      IF(PRSFLG.EQ.2) CALL GETPRS(INERR,*99)
C
C  Purge Functional Inspections program
C  get parameters and effectiveness rates.
C
      IF(PRGFLG.EQ.2) CALL GETPRG(INERR,*99)
C
C  The CHKPGM calls QUITER subroutine when there are parameters that
C  are not all equal.(ATP,IM,Pressure and/or Purge)
C
      CALL CHKPGM(INERR)
C
C  TAMFLG, MYMRFG, IMFLAG, and ATPFLG control inputs affecting the lead
C  phaseout date tampering feature.  Now set up the lead phaseout save area.
C
   60 CALL CALCUM(LPOD)
      CALL GETCUM
      CALL SAVEPB
C
C  Assign trips per day / miles per day.
C
C  CEVTPD calculates trips per day by model year group.  The
C  corresponding miles will be calculated in CEVMPD everytime REGMOD is
C  called.
C
      CALL CEVTPD
C
C  If user requests either Stage II or Onboard controls or both,
C  get user's (HC refueling loss control) program(s) specifications.
C  These parameters apply nationwide.  They do not vary with local conditions
C  and therefore are entered in One-Time Data instead of by scenario.
C
      IF(RLFLAG.GE.2.AND.RLFLAG.LE.4) CALL GETNAT(INERR,*99)
C
C  If the user asks for 1 local area parameters record to apply to all
C  scenarios on the run, read it in here.
C
      IF(LOCFLG.EQ.2) CALL RDCITY(INERR,*99)
C
C  If the user asks for 1 %VMT trip length record to apply to all
C  scenarios on the run, read it it here.
C
      IF(SPDFLG.EQ.4) CALL WTLRUL(INERR,*99)
C
C  If the user asks for by model year emission factor results tables, find
C  out which tables are to be printed.
C
      IF(OUTFMT.EQ.5) CALL GETBYM(INERR,*99)
C
C  Normal return: Either there are no one-time data errors ( => + no other
C  errors => evaluate scenarios) OR there are errors ( => check for errors,
C  but do not process the scenario(s) or + number of errors > 50 =>
C  terminate the run).
C
      RETURN
C
C  EOF on any attempted read => take alternate return 1 => run aborts.
C
   99 RETURN1
      END
