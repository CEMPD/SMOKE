C Version @(#)$Id$ $Source$ $Date$ 

      SUBROUTINE BIGSC1(IP,IV,IDX)
C
C  BIGSC1 computes SALHCF for IVB (Group B IV): HDGV, LDDV, LDDT, and HDDV.
C  It also computes the Idle Speed Correction Factor Adjustment for this group.
C
C  Called by BIGSAL.
C
C  calls EXP.
C
C  Input on call:
C
C    parameter list: IP,IV,IDX
C    common blocks:
C    /SCENE1/ SPD
C    /SPEED3/ SADJ
C    /SPEED4/ IVB
C    /SPEED7/ GPBSCO
C
C  Output on return:
C
C    common blocks:
C    /SPEED6/ SALHCF,SCIADJ
C
C  Local variable / array dictionary:
C
C   Name   Type              Description
C  ------  ----  -------------------------------------------------------
C  SB       R    vehicle speed from SPD, given IV, for IVB scf calculation
C  SPD16    R    accumulator: stores scf calculations for idle SCF adj @ 16mph
C  SPDGPB   R    accumulator: stores scf calculations for subgroup IVB
C
C  Notes:
C
C  BIGSC1 added for MOBILE4.1 v4.
C  BIGSC1 was modified for MOBILE4.1.
C
C
      COMMON /SCENE1/ SPD(8),PCCN,PCHC,PCCC
      COMMON /SPEED3/ SADJ,SCUNA1(18,3,8),SCADJ1(3,4),ISPGRP(3,4)
      COMMON /SPEED4/ MAXSP1,LB1STS(3),LBLAST,IVB,IVA,SCUNA0(18,3,8)
      COMMON /SPEED6/ SALHCF(25,3,8),SCFIDL(25,4),SCIADJ(25,3,8)
      COMMON /SPEED7/ GPBSCO(3,3,4)
C
      SB=SPD(IV)
      SPD16=0.0
C
C  Heavy Duty Vehicle classes
C  (HDGV and HDDV do not have Bags, so no Idle SCF calc.)
C
      IF(IV.EQ.4.OR.IV.EQ.7)
     *    SPDGPB=GPBSCO(3,IP,IVB)*SB*SB
     *          +GPBSCO(2,IP,IVB)*SB
     *          +GPBSCO(1,IP,IVB)
C
C  Light Duty Diesel classes
C  (LDDV and LDDT do have Bags, so compute Idle SCF @ 16.0 mph.)
C
      IF(IV.EQ.5.OR.IV.EQ.6) THEN
        SPDGPB=GPBSCO(3,IP,IVB)*(SB*SB-SADJ*SADJ)
     *        +GPBSCO(2,IP,IVB)*(SB-SADJ)
        SPD16=GPBSCO(3,IP,IVB)*(16.0*16.0-SADJ*SADJ)
     *       +GPBSCO(2,IP,IVB)*(16.0-SADJ)
      ENDIF
C
C  Unless HDDV and NOX, exponentiate the 2nd degree polynomial result.
C
      IF(IV.NE.4.OR.IP.NE.3) THEN
        SPDGPB=EXP(SPDGPB)
        SPD16=EXP(SPD16)
      ENDIF
C
      SALHCF(IDX,IP,IV)=SPDGPB
C
      SCIADJ(IDX,IP,IV)=1.0
      IF(IV.EQ.5.OR.IV.EQ.6) SCIADJ(IDX,IP,IV)=SPD16
C
      RETURN
      END
