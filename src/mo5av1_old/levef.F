C Version @(#)$Id$ $Source$ $Date$ 

      SUBROUTINE LEVEF(ILEV,IP,ICAR,AEF,VMTAGE,MY)
C
C  Determine emissions using California LEV standards
C
C  Parameter list: ILEV, IP, ICAR, AEF, VMTAGE, MY
C
C  Called by GETLEV
C
C  Calls NONE
C
C  Output on return: AEF
C
C    common blocks:
C
C    /BASEQ9/ LEVSTD
C    /BASE11/ SLOIM1, SLOIM2
C    /BASE12/ CAPIV1, CAPIV2
C    /LEVBLK/ IMFLG
C    /TAMEQ2/ F50K
C
C  Local variable / array dictionary:
C
C   Name   Type              Description
C  ------  ----  -------------------------------------------------
C  KINK50   I     kinked at 50K miles ef curve flag: 1 = no, 2 = yes
C  ABOV50   R     if GT 0.0 it contains the remainder over 50K miles
C  DET1     R     holder for SLOPE1 or SLOPE1-0.0082
C  DET2     R     holder for SLOPE2 or SLOPE2-0.0082
C
C
C  Notes:
C
C  LEVEF added for MOBILE5, now it processes Low Emitting Vehicles
C
C  The method of giving or not giving I/M credit to the LFV's
C  is unique in MOBILE5, in that the slope changes:
C
C
C  Determine emissions using California LEV standards
C
C  Parameter list: ILEV, IP, ICAR, AEF, VMTAGE, MY
C
C  Called by GETLEV
C
C  Calls NONE
C
C  Output on return: AEF
C
C    common blocks:
C
C    /BASEQ9/ LEVSTD
C    /BASE11/ SLOIM1, SLOIM2
C    /BASE12/ CAPIV1, CAPIV2
C    /LEVBLK/ IMFLG
C    /TAMEQ2/ F50K
C
C  Local variable / array dictionary:
C
C   Name   Type              Description
C  ------  ----  -------------------------------------------------
C  KINK50   I     kinked at 50K miles ef curve flag: 1 = no, 2 = yes
C  ABOV50   R     if GT 0.0 it contains the remainder over 50K miles
C  DET1     R     holder for SLOPE1 or SLOPE1-0.0082
C  DET2     R     holder for SLOPE2 or SLOPE2-0.0082
C
C
C  Notes:
C
C  LEVEF added for MOBILE5, now it processes Low Emitting Vehicles
C
C  The method of giving or not giving I/M credit to the LFV's
C  is unique in MOBILE5, in that the slope changes:
C
C  IMFLG = 1, => not Special I/M, use Tier 1 Slope(s)
C  IMFLG = 2, => with Special I/M, use LEV Slope(s)
C
C
      REAL LEVSTD, LEVIMP, LDTAB
C
      COMMON /BASEQ9/ LEVSTD(3,6,5,3),LEVIMP(8,12,5),LDTAB(2,2)
      COMMON /BASE11/ SLOIM1, SLOIM2
      COMMON /BASE12/ CAPIV1, CAPIV2
      COMMON /IOUCOM/ IOUIMD,IOUGEN,IOUREP,IOUERR,IOUASK
      COMMON /LEVBLK/ LEVFLG,IMFLG,LEVYRS
      COMMON /TAMEQ2/ MYGTAM(3,4),F50K
C
C  VCAPML caps VMTAGE at CAPIV1 or CAPIV2 mileage
C  capiv1 = 100,000 miles   capiv2 = 120,000 miles
C
      VCAPML=VMTAGE
C
C  LEV milage cap with special I/M
C
      IF(MY.GE.LEVYRS.AND.IMFLG.EQ.2) THEN
         IF(ICAR.EQ.1) THEN
            IF(VCAPML.GT.CAPIV1) VCAPML = CAPIV1
         ELSE
            IF(VCAPML.GT.CAPIV2) VCAPML = CAPIV2
         ENDIF
      ENDIF
C
C  redo the kink50 to be sure of it's timing
C
      ABOV50=VCAPML-F50K
      KINK50=1
      IF(ABOV50.GT.0.0) KINK50=2
C
C  slope explanation
C
C  IMFLG = 1 = > not Special I/M
C  Tier 1 slope1 => SLOIM1-0.0082
C  Tier 1 slope2 => SLOIM2-0.0082
C
C  IMFLG = 2 = > Special I/M
C  LEV slope1 => SLOIM1
C  LEV slope2 => SLOIM2
C
      ZPOINT = 0.0
      DET1   = 0.0
      DET2   = 0.0
C
C  Generic CODE gets the ZPOINT, DET1, and DET2
C
      ZPOINT=LEVSTD(1,ILEV,ICAR,IP)
      DET1=LEVSTD(2,ILEV,ICAR,IP)
      IF(KINK50.EQ.2) DET2=LEVSTD(3,ILEV,ICAR,IP)
C
C  If not special I/M => use the 1993 Tier 1 slopes
C  and subtract out the methane ( 0.0082 )
C
      IF(IMFLG.EQ.1) THEN
C
        DET1 = SLOIM1
        DET2 = SLOIM2
C
C  Only HC get the methane correction
        IF(IP.EQ.1) THEN
          DET1 = (SLOIM1 - 0.0082 )
          DET2 = (SLOIM2 - 0.0082 )
        ENDIF
C
      ENDIF
C
C  Calculate emission factors
C
      IF(KINK50.EQ.1) AEF=ZPOINT+DET1*VCAPML
      IF(KINK50.EQ.2) AEF=ZPOINT+DET1*F50K+DET2*ABOV50
C
      RETURN
      END
