C Version "@(#)$Id$ $Source$ $Date$ 

      SUBROUTINE OUTHD4(ICY)
C
C  OUTHD4 writes the run header on unit IOUREP using output format (OUTFMT) 4.
C  The header consists of the run title and a partial echo of user selected
C  options and data.
C
C  Called by OUTPUT.
C
C  Calls OUTATP, OUTIMP, OUTLAP, OUTMTH, OUTNEW, OUTTAM, USDPTR and OUTDGT.
C
C  Input on call:
C
C    common blocks:
C    parameter list: ICY
C    /ETPPAR/ NEWETP
C    /FLAGS1/ TAMFLG,DTGFLG
C    /FLAGS2/ NEWFLG,IMFLAG
C    /FLAGS3/ ATPFLG
C    /IOUCOM/ IOUREP
C    /PRGCH2/ PRGFLG
C    /PRSCH2/ PRSFLG
C
C  Output on return:
C
C    common block: /REGION/ INITPR
C
C  Notes:
C
C  The call to TAMPER was restored in Version I00.
C
C  September-14-1994 @CSC-pme request 2-464. Include message that indicate 
C  that the effects of detergent gasoline additive rulemaking have been 
C  altered by the user. Subroutine OUTDTG was created to print the  message. 
C  A call to subroutine OUTDTG was added. (One LAP record for all scenarios) 
C
C  08-Feb-1996 @ DynTel Request 2-615 Detergent gasoline effects are disabled
C                when NEWFLG is set to 5 or 6.
C
      INCLUDE 'ETPPAR.I'
      INCLUDE 'FLAGS1.I'
      INCLUDE 'FLAGS2.I'
      INCLUDE 'FLAGS3.I'
      INCLUDE 'IOUCOM.I'
      INCLUDE 'PRGCH2.I'
      INCLUDE 'PRSCH2.I'
      INCLUDE 'REGION.I'
C
      INITPR=0
C
      IF(TAMFLG.EQ.2) CALL OUTTAM
C
      IF(NEWFLG.EQ.2.OR.NEWFLG.EQ.4.OR.NEWFLG.EQ.6) CALL OUTNEW
C
      IF(NEWFLG.EQ.3.OR.NEWFLG.EQ.4) CALL OUTETP
C
      IF(IMFLAG.NE.1) CALL OUTIMP
C
      IF(ATPFLG.EQ.2.OR.PRSFLG.EQ.2.OR.PRGFLG.EQ.2) CALL OUTATP
C
      IF(RLFLAG.EQ.2) CALL OUTSTG
C
      IF(LOCFLG.EQ.2) CALL OUTLAP
C
      IF(DSFLAG.EQ.2) CALL OUTGSF(ICY)
C
      IF(DTGFLG.EQ.2.AND.LOCFLG.EQ.2.AND.NEWFLG.LT.5) CALL OUTDTG
C
      CALL OUTMTH
C
      WRITE(IOUREP,200)
  200 FORMAT('0',80('_'))
C
      CALL USDPTR
C
      RETURN
      END
