C Version @(#)$Id$ $Source$ $Date$ 

      SUBROUTINE GETNAT(INERR,*)
C
C  GETNAT reads in and processes refueling loss related HC parameters
C  that are to apply nationwide - that is, to all scenarios of the run.
C
C  Called by ONESEC.
C
C  Calls YRTEST.
C
C  Input on call:
C
C    common blocks:
C    /FLAGS1/ PROMPT
C    /FLAGS3/ RLFLAG
C    /IOUCOM/ IOUGEN,IOUASK
C    /YEARS4/ IY2020
C
C  Output on return:
C
C    parameter list: INERR
C    common blocks:
C    /RLCOM2/ IOBMY,IVOB,IS2SY,NPHASE,S2EFF
C
C  Local variable / array dictionary:
C
C   Name   Type              Description
C  ------  ----  --------------------------------------------------------------
C  LOWVYR   I    lowest allowed S2 or ob vrs installation start year
C  PCTEL    I    percent efficiency of S2 vrs in reducing LDGV/T refueling loss
C  PCTEH    I    percent efficiency of S2 vrs in reducing HDGV   refueling loss
C
C  Notes:
C
C  None.
C
C
      INTEGER PROMPT,TAMFLG,SPDFLG,VMFLAG,OXYFLG,DSFLAG
      INTEGER ATPFLG,TPDFLG,RLFLAG,TEMFLG,OUTFMT
C
      COMMON /FLAGS1/ PROMPT,TAMFLG,SPDFLG,VMFLAG,OXYFLG,DSFLAG
      COMMON /FLAGS3/ ATPFLG,TPDFLG,RLFLAG,LOCFLG,TEMFLG,OUTFMT
      COMMON /IOUCOM/ IOUIMD,IOUGEN,IOUREP,IOUERR,IOUASK
      COMMON /RLCOM2/ IOBMY,IVOB(4),IS2SY,NPHASE,S2EFF(4)
      COMMON /YEARS4/ IY1941,IY1960,IY1994,IY2020
C
      DATA LOWVYR/1989/
C
C  If uncontrolled losses or no refueling losses are selected, then skip
C  the next section.
C
      IF(RLFLAG.EQ.1.OR.RLFLAG.EQ.5) GOTO 30
C
C  Vapor Recovery Systems (vrs) input parameters section:
C
C  Stage II vrs input: start year, no. of phase-in years & percent efficiency
C  for LDG & HDG are entered.
C
      IF(RLFLAG.EQ.3) GOTO 10
C
      IF(PROMPT.EQ.2) WRITE(IOUASK,200)
  200 FORMAT(' ','Enter Stage II VRS parameters -'/
     *      ' ','start year, # phase-in years, LDG & HDG % efficiency:')
      READ(IOUGEN,100,ERR=99,END=99) IS2SY,NPHASE,PCTEL,PCTEH
  100 FORMAT(I2,1X,I1,2(1X,F3.0))
C
      CALL YRTEST(IS2SY,87,LOWVYR,IY2020,INERR)
      IF(NPHASE.EQ.0) CALL QUITER(0.,NPHASE,86,INERR)
      IF(PCTEL.LT.0.0.OR.PCTEL.GT.100.0) CALL QUITER(PCTEL,1,99,INERR)
      IF(PCTEH.LT.0.0.OR.PCTEH.GT.100.0) CALL QUITER(PCTEH,2,99,INERR)
C
      S2EFF(1)=PCTEL/100.
      S2EFF(2)=S2EFF(1)
      S2EFF(3)=S2EFF(2)
      S2EFF(4)=PCTEH/100.
      IF(RLFLAG.EQ.2) GOTO 30
C
C  Users specifies onboard vrs controlled refueling losses.  The user must
C  supply the start year for controls and which vehicle classes are to be
C  covered.
C
   10 IF(PROMPT.EQ.2) WRITE(IOUASK,210)
  210 FORMAT(' ','Enter onboard vrs first model year & vehicle classes',
     *           ' (I2,1X,4I1):')
      READ(IOUGEN,110,ERR=99,END=99) IOBMY,IVOB
  110 FORMAT(I2,1X,4I1)
C
      CALL YRTEST(IOBMY,87,LOWVYR,IY2020,INERR)
C
      DO 20 IVTAM=1,4
      IF(IVOB(IVTAM).LT.1.OR.IVOB(IVTAM).GT.2)
     *      CALL QUITER(0.,IVOB(IVTAM),84,INERR)
   20 CONTINUE
C
C  RLFLAG one-time data processing is completed.
C
   30 RETURN
C
C  Error or EOF on any attempted read => take alternate
C  return 1 => run aborts.
C
   99 CALL QUITER(0.,0,130,INERR)
C
      RETURN1
      END
