C Version @(#)$Id$ $Source$ $Date$ 

      SUBROUTINE GETETP(INERR,*)
C
C  GETETP gets and verifies the evaporative test procedure parameters.
C
C  Called by ONESEC.
C
C  Calls YRTEST.
C
C  Input on call:
C
C    parameter list: INERR
C    common Blocks:
C    /FLAGS1/ PROMPT
C    /IOUCOM/ IOUGEN,IOUASK
C    /YEARS4/ IY1960,IY2020
C
C  Output on return:
C
C    Parameter list: INERR
C    common blocks:
C    /ETPPAR/ ETP,NEWETP,IETPYR,ETPPIV
C
C  Local array subscripts :
C
C  ETPNAM(3)  -  ETPNAM ( KEYEQ )
C
C  Local variable / array dictionary:
C
C  Name   Type              Description
C  ------ ----  ----------------------------------------------------
C  ETPNAM  C    names of the ETP groups: Problem Free, Purge Failure,
C               and Pressure Failure
C
C  Notes:
C
C  GETETP was added for MOBILE5.
C
C
      INTEGER PROMPT,TAMFLG,SPDFLG,VMFLAG,OXYFLG,DSFLAG
      CHARACTER*17 ETPNAM(3)
C
      COMMON /FLAGS1/ PROMPT,TAMFLG,SPDFLG,VMFLAG,OXYFLG,DSFLAG
      COMMON /ETPPAR/ ETP(7,3),NEWETP,IETPYR(4),ETPPIV(4),NETPYR
      COMMON /IOUCOM/ IOUIMD,IOUGEN,IOUREP,IOUERR,IOUASK
      COMMON /YEARS4/ IY1941,IY1960,IY1994,IY2020
C
      DATA ETPNAM/'Problem Free:    ',
     *            'Purge Failure:   ',
     *            'Pressure Failure:'/
C
C  Obtain number of evaporative test procedure phase in years and
C  new rates flag. Prompt only for number of phase-in years.
C
      IF(PROMPT.EQ.2) WRITE(IOUASK,200)
  200 FORMAT(' ','Enter the number of evaporative ',
     *           'test procedure phase-in years:')
      READ(IOUGEN,100,END=99,ERR=99) NETPYR,NEWETP
  100 FORMAT(I1,I2)
C
      IF(NETPYR.LT.1.OR.NETPYR.GT.4) CALL QUITER(0.,NETPYR,115,INERR)
C
      IF(NEWETP.LT.1) NEWETP=1
      IF(NEWETP.GT.2) NEWETP=2
C
C  Acquire and verify phase-in years and corresponding percentages.
C
      DO 5 IETPY=1,NETPYR
        IF(PROMPT.EQ.2) WRITE(IOUASK,205) IETPY
  205   FORMAT(' ','Enter the',I2,' phase-in year and percentage:')
        READ(IOUGEN,105,END=99,ERR=99) IETPYR(IETPY),ETPPIV(IETPY)
  105   FORMAT(I2,F5.0)
C
        CALL YRTEST(IETPYR(IETPY),139,1995,IY2020,INERR)
        IF(ETPPIV(IETPY).LT.0..OR.ETPPIV(IETPY).GT.100.)
     *    CALL QUITER(ETPPIV(IETPY),0,116,INERR)
    5 CONTINUE
C
C  Obtain the alternate ETP rates.
C
      IF(NEWETP.EQ.2) THEN
C
      DO 10 KEYEQ=1,3
      IF(PROMPT.EQ.2) WRITE(IOUASK,210) ETPNAM(KEYEQ)
  210 FORMAT(' ','Enter the Hot Soak, Full Diurnal, Multiple Diurnal'/
     *       ' ','2-3, Multiple Diurnal 4+, Partial Diurnal, Running'/
     *       ' ','Loss, and Resting Loss rates (%) for ',A17)
      READ(IOUGEN,110,END=99,ERR=99) (ETP(IETP,KEYEQ),IETP=1,7)
  110 FORMAT(7F5.0)
   10 CONTINUE
C
C  Verify ETP parameters:
C  Values must lie between 0 and 100%, inclusive.
C  The full and partial diurnal rates must be equal.
C  The running loss purge and pressure values must be equal.
C  All the resting loss values must be the same.
C
      DO 20 KEYEQ=1,3
      DO 20 IETP=1,7
      IF(ETP(IETP,KEYEQ).LT.0..OR.ETP(IETP,KEYEQ).GT.100.)
     *    CALL QUITER(ETP(IETP,KEYEQ),0,140,INERR)
   20 CONTINUE
C
      IF(ETP(2,1).NE.ETP(5,1).OR.
     *   ETP(2,2).NE.ETP(5,2).OR.
     *   ETP(2,3).NE.ETP(5,3)) CALL QUITER(0.,0,141,INERR)
      IF(ETP(6,2).NE.ETP(6,3)) CALL QUITER(0.,0,142,INERR)
      IF(ETP(7,1).NE.ETP(7,2).OR.ETP(7,1).NE.ETP(7,3))
     *    CALL QUITER(0.,0,143,INERR)
C
      ENDIF
C
      RETURN
C
C  EOF or ERR on any attempted read => take alternate return 1 => run aborts.
C
   99 CALL QUITER(0.,0,144,INERR)
C
      RETURN1
      END
