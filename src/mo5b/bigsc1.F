C Version "@(#)$Id$ $Source$ $Date$ 

      SUBROUTINE BIGSC1(IP,IV,IDX)
C
C  BIGSC1 computes SALHCF for IVB (Group B IV): HDGV, LDDV, LDDT, and HDDV.
C  It also computes the Idle Speed Correction Factor Adjustment for this group.
C
C  Called by BIGSAL.
C
C  Calls EXP.
C
C  Input on call:
C
C    parameter list: IP,IV,IDX
C    common blocks:
C    /SCENE1/ SPD
C    /SPEED3/ SADJ
C    /SPEED4/ IVB
C    /SPEED6/ SALHCF,SCFIDL,IDLYES
C    /SPEED7/ GPBSCO
C
C  Output on return:
C
C    common blocks:
C    /SPEED6/ SALHCF,SCFIDL
C
C  Local variable / array dictionary:
C
C   Name   Type              Description
C  ------  ----  -------------------------------------------------------
C  SB       R    vehicle speed from SPD, given IV, for IVB scf calculation
C  SPDGPB   R    accumulator: stores scf calculations for subgroup IVB
C
C  Notes:
C
C  BIGSC1 added for MOBILE4.1 v4.
C  BIGSC1 was modified for MOBILE4.1.
C  BIGSC1 was modified for MOBILE5a in house version for new idle algorithm.
C  15 July 1993 @ ARC-rwk Task 2-277
C
C
      INCLUDE 'SCENE1.I'
      INCLUDE 'SPEED3.I'
      INCLUDE 'SPEED4.I'
      INCLUDE 'SPEED6.I'
      INCLUDE 'SPEED7.I'
C
      SB=SPD(IV)
      IF(IDLYES) SB=2.5
C
C  Heavy Duty Vehicle classes
C
      IF(IV.EQ.4.OR.IV.EQ.7) THEN
          SPDGPB=GPBSCO(3,IP,IVB)*SB*SB
     *          +GPBSCO(2,IP,IVB)*SB
     *          +GPBSCO(1,IP,IVB)
C
C  Light Duty Diesel classes
C
      ELSEIF(IV.EQ.5.OR.IV.EQ.6) THEN
        SPDGPB=GPBSCO(3,IP,IVB)*(SB*SB-SADJ*SADJ)
     *        +GPBSCO(2,IP,IVB)*(SB-SADJ)
      ENDIF
C
C  Unless HDDV and NOX, exponentiate the 2nd degree polynomial result.
C
      IF(IV.NE.4.OR.IP.NE.3) SPDGPB=EXP(SPDGPB)
C
      IF(IDLYES) THEN
        SCFIDL=SPDGPB
      ELSE
        SALHCF(IDX,IP,IV)=SPDGPB
      ENDIF
C
      RETURN
      END
