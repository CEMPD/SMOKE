C Version @(#)$Id$ $Source$ $Date$ 

      SUBROUTINE BIGSAL(IP,IV,MY,IDX,ISP)
C
C  BIGSAL computes SALHCF, the speed and, if LDGV/T or MC, optional adjustments
C  multiplicative correction factor table and computes SCFIDL.
C  Also computes SCIADJ, the idle speed correction factor adjustment.
C
C  Called by BIGCFX.
C
C  Calls ALHRET, BIGSC1, BIGSC2, and BIGSC3.
C
C  Input on call:
C
C    parameter list: IP,IV,MY,IDX,ISP
C    common blocks:
C    /ALUOUT/ ALHRET
C    /FLAGS2/ ALHFLG
C    /SCENE1/ SPD
C    /SPEED4/ LB1STS,IVA
C    /SPEED6/ SALHCF
C
C  Output on return:
C
C    common blocks:
C    /SPEED6/ SALHCF,SCFIDL,SCIADJ
C
C  Notes:
C
C  BIGSAL was modified for MOBILE4.1 v4 to separate into subroutines the
C  primary algorithm cases.
C
C
      INTEGER ALHFLG
C
      COMMON /ALUOUT/ ALHRET(25,4,3)
      COMMON /FLAGS2/ MYMRFG,NEWFLG,IMFLAG,ALHFLG,IMTIER
      COMMON /SCENE1/ SPD(8),PCCN,PCHC,PCCC
      COMMON /SPEED4/ MAXSP1,LB1STS(3),LBLAST,IVB,IVA,SCUNA0(18,3,8)
      COMMON /SPEED6/ SALHCF(25,3,8),SCFIDL(25,4),SCIADJ(25,3,8)
C
C SCFIDL for BEFIDL'S OFFCO correction
C
      IF(IP.EQ.2) THEN
      IF((MY.GE.1981.AND.IV.EQ.1).OR.(MY.GE.1984.AND.IV.EQ.2) .OR.
     1   (MY.GE.1984.AND.IV.EQ.3)) THEN
         SCFIDL(IDX,IVA)=1.0
         CALL IDLSC3(IP,MY,IDX)
      ENDIF
      ENDIF
C
      SALHCF(IDX,IP,IV)=1.0
      SCIADJ(IDX,IP,IV)=1.0
C
      GOTO(20,20,20,10,10,10,10,20),IV
C
C  SCF for Group B IV: HDGV, LDDV, LDDT, and HDDV.
C
   10 CALL BIGSC1(IP,IV,IDX)
      RETURN
C
C  SCF for Group A IV: LDGV, LDGT1, LDGT2 and MC.
C
C  If optional corrections selected, correct LDGV, LDGT1, LDGT2 and MC
C  for air conditioning, extra load, trailer towing and/or humidity.
C
   20 IF(ALHFLG.GT.1) SALHCF(IDX,IP,IV)=ALHRET(IDX,IVA,IP)
C
C  For Speed cf, either use old MOBILE3 algorithm,
C
C  Note: Split the IF, since at least 1 compiler chokes on LB1STS(8).
C
C  Motor Cycles (IV=8) All Speeds and all IP's
      IF(IV.EQ.8) THEN
      CALL BIGSC2(IP,IV,MY,IDX,ISP)
      RETURN
      ENDIF
C
C  Pre 1980 LDGV, LDGT1, and LDGT2
      IF(MY.LT.LB1STS(IV)) THEN
      CALL BIGSC2(IP,IV,MY,IDX,ISP)
      RETURN
      ENDIF
C
C  Or use new equation form, based on speed and pollutant type.
C
C  Low, Mid and High Speed for all three IP's
      IF(MY.GE.LB1STS(IV)) THEN
       CALL BIGSC3(IP,MY,IDX,ISP)
       RETURN
       ENDIF
C
      RETURN
      END
