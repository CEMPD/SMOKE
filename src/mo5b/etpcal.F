C Version "@(#)$Id$ $Source$ $Date$ 

      FUNCTION ETPCAL(IEVTYP,MY,KEYEQ,ETPRVP,ETPTMP)
C
C  ETPCAL computes the Evaporative Test Procedure rate benefit.
C
C  Called by HOTSOK, DIURNL, RSTLOS and RULOSS.
C
C  Calls ETPRED.
C
C  Input on call:
C
C    parameter list: IEVTYP,MY,KEYEQ,ETPRVP,ETPTMP
C    common blocks:
C    /ETPPAR/ ETP,IETPYR,ETPPIV,NETPYR
C    /FLAGS2/ NEWFLG
C
C  Output on return:
C
C    funtion: ETPCAL
C
C  Local variable / array dictionary:
C
C   Name   Type              Description
C  ------  ----  -------------------------------------------------------
C  ETPRVP   R    Evaporative Test Procedure RVP
C  ETPTMP   R    Evaporative Test Procedure temperature
C  ETPEFF   R    Evaporative Test Procedure effectiveness
C
C  Notes:
C
C  ETPCAL was added for MOBILE5v3.
C
      INCLUDE 'ETPPAR.I'
      INCLUDE 'FLAGS2.I'
C
C  Initialize.
C
      ETPCAL=0.0
C
C  Return if MY is before ETP start year or if
C  user has specified not to include ETP.
C
      IF(MY.LT.IETPYR(1).OR.NEWFLG.GE.5) RETURN
C
C  Find phase-in rate.
C
      DO 10 IETPY=1,NETPYR
      IF(MY.GE.IETPYR(IETPY)) ETPCAL=ETPPIV(IETPY)/100.
   10 CONTINUE
C
C  Find ETP effectiveness based on evaporative (hot soaks,
C  diurnals, etc) type and pass/purge/pressure status.
C
      ETPEFF=ETP(IEVTYP,KEYEQ)
C
C  However, if fail and hot soak or running loss, calculate the
C  effectiveness rate based on RVP, temperature and user input,
C  if supplied.
C
      IF(KEYEQ.GE.2.AND.(IEVTYP.EQ.1.OR.IEVTYP.EQ.6))
     *  ETPEFF=ETPRED(IEVTYP,KEYEQ,ETPRVP,ETPTMP)
C
C  Determine the final ETP rate.
C
      ETPCAL=ETPCAL*ETPEFF/100.
C
      RETURN
      END
