C Version "@(#)$Id$ $Source$ $Date$ 

      SUBROUTINE GETALH(INERR,*)
C
C  GETALH reads in and validates the optional correction factors (BIGALH)
C  parameters record.
C
C  Called by PARSEC.
C
C  Calls QUITER and YRTEST.
C
C  Input on call:
C
C    parameter list: INERR
C    common blocks:
C    /FLAGS1/ PROMPT
C    /FLAGS2/ ALHFLG
C    /IOUCOM/ IOUGEN,IOUASK
C    /MAXIMA/ MAXLTW
C
C  Output on return:
C
C    parameter list: INERR
C    common blocks:
C    /ALUHIN/ AC,XLOAD,TRAILR,ABSHUM,DB,WB
C
C  Notes:
C
C  None.
C
C
      INCLUDE 'ALUHIN.I'
      INCLUDE 'FLAGS1.I'
      INCLUDE 'FLAGS2.I'
      INCLUDE 'IOUCOM.I'
      INCLUDE 'MAXIMA.I'
C
C  GETALH is called once for each scenario of the run.
C
C  Assign db/wb default temperatures for a/c usage (90% for defaults).
C
      DB=85.
      WB=75.
C
C  Read in BIGCFX parameters only if user flags their presence with ALHFLG.
C  ALHFLG also determines the choice of inputs and the corresponding format.
C
      IF(ALHFLG.EQ.1) GOTO 50
      IF(ALHFLG.EQ.3) GOTO 20
C
      IF(PROMPT.EQ.2) WRITE(IOUASK,200)
  200 FORMAT(' ','Enter AC,XLOAD(3),TRAILR(1),ABSHUM')
      READ(IOUGEN,100,ERR=99,END=99) AC,XLOAD,TRAILR(1),ABSHUM
  100 FORMAT(5F4.3,F4.1)
C
      DO 10 IVP=2,MAXLTW
      TRAILR(IVP)=TRAILR(1)
   10 CONTINUE
      GOTO 30
C
   20 IF(PROMPT.EQ.2) WRITE(IOUASK,210)
  210 FORMAT(' ','Enter AC,XLOAD(3),TRAILR(3),ABSHUM,DBTEMP,WBTEMP')
      READ(IOUGEN,110,ERR=99,END=99) AC,XLOAD,TRAILR,ABSHUM,DB,WB
  110 FORMAT(7F4.3,3F4.1)
C
C  Validate inputs for BIGALH.
C
   30 IF(AC.LT.0.0.OR.AC.GT.1.) CALL QUITER(AC,0,44,INERR)
      IF(AC.GT.0.0) CALL QUITER(0.,0,56,INERR)
C
      DO 40 IVP=1,MAXLTW
      IF(XLOAD(IVP).LT.0.0.OR.XLOAD(IVP).GT.1.)
     *   CALL QUITER(XLOAD(IVP),0,45,INERR)
      IF(TRAILR(IVP).LT.0.0.OR.TRAILR(IVP).GT.1.)
     *   CALL QUITER(TRAILR(IVP),0,46,INERR)
   40 CONTINUE
C
      IF(ABSHUM.LT.20.0.OR.ABSHUM.GT.140.0)
     *    CALL QUITER(ABSHUM,0,47,INERR)
C
      IF(DB.LT.0.0.OR.DB.GT.110.0) CALL QUITER(DB,0,40,INERR)
      IF(WB.LT.0.0.OR.WB.GT.110.0) CALL QUITER(WB,0,40,INERR)
      IF(WB.GT.DB) CALL QUITER(0.,0,57,INERR)
   50 RETURN
C
C  Error or EOF on any attempted read => take alternate
C  return 1 => run aborts.
C
   99 CALL QUITER(0.,0,133,INERR)
C
      RETURN1
      END
