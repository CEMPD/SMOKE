      REAL FUNCTION BIGSC2(IP,IV,MY,SPD,INERR)
C
C  BIGSC2 computes SALHCF for older (MY prior to
C  LB1STS(IV) ) LDGV/T, and for the High Speed ( .GT. 48 MPH ) CF.
C
C  Called by SCFSP.
C
C  Calls ISPPTR and EXP.
C
C  Changes: (Last change first)
C
C  15 Apr 02  AIR Task Bug Fix 388: modification of constants via QUITER.
C  17 Aug 99 @DynTel-ddj 2-698 Changed to a function and removed IDX
C                and ISR from parameter list.
C  12 May 99 @DynTel-HQ 2-698 simplify the code for new speed
C                correction
C  20 Apri 99 @DynTel-HQ 2-698 simplify the code for new speed
C                correction
C
C  Input on call:
C
C    argument list: IP,IV,MY,IDX
C    common blocks:
C    /SPEED1/ SK1
C    /SPEED3/ SADJ
C
C  Output on return:
C
C  Local variable / array dictionary:
C
C   Name   Type              Description
C  ------  ----  -------------------------------------------------------
C  BB       R    slope coefficient "BB" of high speed scf equ.
C  IVA2     I    used to point to the LDGV for motorcycles in SIXT5.
C  MYR      I    model year pointer for SIXT5 array
C  SCF1     R    old scf equation computed using opmode adjusted speed value
C  SP       R    Vehicle speed for high speed equation
C
C  Notes:
C
C  BIGSC2 added for MOBILE4.1 v4.
C  BIGSC2 was modified for MOBILE4.1.
C  BIGSC2 was modified for MOBILE5a in house version for new idle algorithm.
C  15 July 1993 @ ARC-rwk Task 2-277
C  11 Jun 97 DynTel@MLA 2-622  Added IVTYPE.I for Mobile6 and modified the
C            code to use the symbolic vehicle type names and vehicle type
C            vectors that were introduced in Mobile6.
C  16 December 1997 Dyntel@MMS  modified multiple speed to single speed.
C                   SPD array changed to a single variable.
C  20 Apr 98 @DynTel-RJD 2-654 Added an extra dimension to the SALHCF array
C            to distinguish between start and running emissions. No speed
C            correction is needed for start emissions.
C
      IMPLICIT NONE 
      INCLUDE 'SPEED1.I'
      INCLUDE 'SPEED3.I'
      INCLUDE 'SPEED4.I'
C
      INTEGER, INTENT(IN) :: MY
      INTEGER, INTENT(IN) :: IP
      INTEGER, INTENT(IN) :: IV
      INTEGER, INTENT(INOUT) :: INERR
      INTEGER ISP, ICOEF5, ISPPTR
C       
      REAL SCF1, SCF, SPD
C
      BIGSC2=0.0
C
      ISP=ISPPTR(MY,0,IV,1,INERR)
      IF(INERR.GT.0) GOTO 99
C
C  Use MOBILE3's 5th degree polynomial equation.
C
      SCF1=SK1(6,ISP,IP)
      SCF=SK1(6,ISP,IP)

      DO 10 ICOEF5 = 1, 5
         SCF1=SCF1*SADJ+SK1(6-ICOEF5,ISP,IP)
         SCF=SCF*SPD+SK1(6-ICOEF5,ISP,IP)
   10 CONTINUE
C
      IF(IP.NE.3) THEN
         SCF1=EXP(SCF1)
         SCF=EXP(SCF)
      ENDIF
C
      BIGSC2 = SCF/SCF1
C
   99 RETURN
      END
