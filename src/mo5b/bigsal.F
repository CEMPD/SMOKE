C Version "@(#)$Id$ $Source$ $Date$ 

      SUBROUTINE BIGSAL(IP,IV,MY,IDX,ISP)
C
C  BIGSAL computes SALHCF, the speed and, if LDGV/T or MC, optional adjustments
C  multiplicative correction factor table and computes SCFIDL.
C
C  Called by BIGCFX and IDLCAL.
C
C  Calls BIGSC1, BIGSC2, and BIGSC3.
C
C  Input on call:
C
C    parameter list: IP,IV,MY,IDX,ISP
C    common blocks:
C    /ALUOUT/ ALHRET
C    /FLAGS2/ ALHFLG
C    /SCENE1/ SPD
C    /SPEED4/ LB1STS,IVA
C    /SPEED6/ SALHCF,IDLYES
C
C  Output on return:
C
C    common blocks:
C    /SPEED6/ SALHCF,SCFIDL
C
C  Notes:
C
C  BIGSAL was modified for MOBILE4.1 v4 to separate into subroutines the
C  primary algorithm cases.
C  BIGSAL was modified for MOBILE5a in house version for new idle algorithm.
C  15 July 1993 @ ARC-rwk Task 2-277
C
C
      INCLUDE 'ALUOUT.I'
      INCLUDE 'FLAGS2.I'
      INCLUDE 'SCENE1.I'
      INCLUDE 'SPEED4.I'
      INCLUDE 'SPEED6.I'
C
      IF(IDLYES) THEN
        SCFIDL=1.0
      ELSE
        SALHCF(IDX,IP,IV)=1.0
      ENDIF
C
C  If optional corrections selected, correct LDGV, LDGT1, LDGT2 and MC
C  for air conditioning, extra load, trailer towing and/or humidity.
C
      IF((IV.LE.3.OR.IV.EQ.8).AND.ALHFLG.GT.1) THEN
        IF(IDLYES) THEN
          SCFIDL=ALHRET(IDX,IVA,IP)
        ELSE
          SALHCF(IDX,IP,IV)=ALHRET(IDX,IVA,IP)
        ENDIF
      ENDIF
C
C  SCF for Group A IV: LDGV, LDGT1, LDGT2 and MC.
C
C  LDGV & T
C
      IF(IV.LE.3) THEN
C
C  For Speed cf, either use old MOBILE3 algorithm,
C  Pre 1980 LDGV, LDGT1, and LDGT2
C
        IF(MY.LT.LB1STS(IV)) THEN
          CALL BIGSC2(IP,IV,MY,IDX,ISP)
C
C  Or use new equation form, based on speed and pollutant type.
C  Low, Mid and High Speed for all three IP's
C  Post 1980 LDGV, LDGT1, and LDGT2
C
        ELSE
          CALL BIGSC3(IP,MY,IDX,ISP)
        ENDIF
C
C  Motorcycles use old MOBILE3 algorithm.
C
      ELSEIF(IV.EQ.8) THEN
        CALL BIGSC2(IP,IV,MY,IDX,ISP)
C
C  SCF for Group B IV: HDGV, LDDV, LDDT, and HDDV.
C
      ELSE
        CALL BIGSC1(IP,IV,IDX)
      ENDIF
C
      RETURN
      END
