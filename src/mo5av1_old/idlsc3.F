C Version @(#)$Id$ $Source$ $Date$ 

      SUBROUTINE IDLSC3(IP,MY,IDX)
C
C  IDLSC3 calculates SCFIDL ( Speed Correction Factor at 2.5MPH - Idle)
C  for correcting OFFCO ( CO Offset to Idle ) in Subroutine BEFIDL.
C
C  Called by BIGSAL.
C
C  Calls ISPPTR.
C
C  Input on call:
C
C    parameter list: IP,MY,IDX
C    common blocks:
C    /SCENE1/ SPD
C    /SPEED1/ SK2C,SK4C
C    /SPEED3/ SADJ
C    /SPEED4/ IVA
C
C  Output on return:
C
C    common blocks:
C    /SPEED6/ SCFIDL
C
C  Local variable / array dictionary:
C
C   Name   Type              Description
C  ------  ----  -------------------------------------------------------
C  A        R    coefficient "A" of new scf equation: A/X+B
C  A2       R    SADJ based coefficient "A" of new scf equation: A/X+B
C  B        R    coefficient "B" of new scf equation: A/X+B
C  B2       R    SADJ based coefficient "B" of new scf equation: A/X+B
C  SCADJ2   R    new scf equation computed using opmode adjusted speed value
C  SCUN25   R    uses the 2.5 mph for the calculating of SCFIDL
C
C  Notes:
C
C  IDLSC3 added for MOBILE4.1 v5
C
C  Attempted calling IDLSC3 from Sub. BEFIDL instead of BIGSAL, but when
C  so kept clobbering the index into the Block Data Arrays.
C
C
      COMMON /SPEED1/ SK1(6,18,3),SK2C(2,14,3,3),SK4C(2,14,3,3)
      COMMON /SPEED3/ SADJ,SCUNA1(18,3,8),SCADJ1(3,4),ISPGRP(3,4)
      COMMON /SPEED4/ MAXSP1,LB1STS(3),LBLAST,IVB,IVA,SCUNA0(18,3,8)
      COMMON /SPEED6/ SALHCF(25,3,8),SCFIDL(25,4),SCIADJ(25,3,8)
C
C  The equation form for SCFIDL is (A/X+B) / (A/Y+B), where X = 2.5 MPH( Idle )
C  and Y = operating mode adjusted FTP speed ( SADJ )
C
      ISP=ISPPTR(MY,IP,IVA,3)
      IF(IVA.GT.1) ISP=ISP+1
      IF(MY.LT.1980) ISP=1
      IF(ISP.GT.14) ISP=14
C
      A=SK2C(1,ISP,IP,IVA)
      B=SK2C(2,ISP,IP,IVA)
      A2=A
      B2=B
      IF(SADJ.GT.19.6) THEN
      A2=SK4C(1,ISP,IP,IVA)
      B2=SK4C(2,ISP,IP,IVA)
      ENDIF
C
      SCUN25=(A/2.5)+B
      SCADJ2=(A2/SADJ)+B2
C
C Calculate the BEFIDL variable SCFIDL
C
      SCFIDL(IDX,IVA)=(SCUN25/SCADJ2)
C
      RETURN
      END
