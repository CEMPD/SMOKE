      REAL FUNCTION ETPCAL(IEVTYP,MY,KEYEQ,ETPRVP,ETPTMP) 
C 
C  ETPCAL computes the Evaporative Test Procedure rate benefit. 
C 
C  Called by EVP_HRLY, HOTSOK, DIUR24, and RNLOSS. 
C 
C  Calls ETPRED. 
C 
C Changes : (Last change first) 
C 
C  06 Oct 00 @EPA-BG  BUG155 removed NEWETP flag previously used for user input 
C            of ETP phase in ; removed code asociated with NEWETP and user supplied 
C            phase in, i.e., removed KEYEQ and IEVTYP from argument list of ETPRED 
C  11 Mar 1999 @DynTel-bag 2-664 Explicitly named a REAL FUNCTION 
C  23 sept 1998 @DynTel-ZK 2-000 Explicit Typing 
C 
C  29 NOV 1999  @ELG-EPA Added new local variable PHASEIN 
C               Phase-in of ETPCAL using ETPPIV is not done for Hot Soak. 
C               It is now done in HOTSOK. 
C 
C  Input on call: 
C 
C    argument list: IEVTYP,MY,KEYEQ,ETPRVP,ETPTMP 
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
C  PHASEIN  R    ETP Phase-in Variable. 
C 
C  Notes: 
C 
C  ETPCAL was added for MOBILE5v3. 
C 
      IMPLICIT NONE 
      INCLUDE 'ETPPAR.I' 
      INCLUDE 'FLAGS2.I' 
      INTEGER  MY, IETPY, IEVTYP, KEYEQ 
      REAL ETPEFF, ETPRVP, ETPTMP, ETPRED 
      REAL PHASEIN 
 
C  Initialize. 
C 
      ETPCAL=0.0 
      PHASEIN = 1.0 
C 
C  Return if MY is before ETP start year or if 
C  user has specified not to include ETP. 
C 
      IF(MY.LT.IETPYR(1).OR.NEWFLG.GE.5) RETURN 
C 
C  Find phase-in rate for all evap types except Hot Soak. 
C  Hot Soak ETP Phase-in is calculated in HOTSOK 
C 
      IF(IEVTYP.NE.1) THEN 
        DO 10 IETPY=1,NETPYR 
          IF(MY.GE.IETPYR(IETPY)) PHASEIN=ETPPIV(IETPY)/100. 
   10   CONTINUE 
      ELSE 
        PHASEIN = 1.0 
      ENDIF 
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
     *  ETPEFF=ETPRED(ETPRVP,ETPTMP) 
C 
C  Determine the final ETP rate. 
C 
      ETPCAL=PHASEIN*ETPEFF/100. 
C 
      RETURN 
      END 
