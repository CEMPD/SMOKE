C Version "@(#)$Id$ $Source$ $Date$ 

      SUBROUTINE BEFUN(MY,IP,IV,VMTAGE,BER,IBEFSW)
C
C  BEFUN returns the basic emission factor Uncorrected for
C  operating mode, CO temperature offset, methane offset, temperature,
C  inspection / maintenance, tampering offset and reformulated fuels.
C
C  Called by BEF and GETLEV.
C
C  Calls IERPTR
C
C  Input on call or from calls:
C
C    parameter list: MY,IP,IV,VMTAGE,IBEFSW
C    common blocks:
C    /BASEQ1/ ERBZML,ERBDR,ERB50K
C    /BASEQ5/ ERUZML,ERUDR,ERU50K
C    /FLAGS2/ NEWFLG
C    /MAXIMA/ MAXYRS
C    /REGION/ IREJN
C    /TAMEQ2/ F50K
C    /TIER1S/ TR1SL1
C
C  Output on return:
C
C    function: BER
C    common blocks:
C    /BASE11/SLOIM1,SLOIM2
C    /LEVBLK/ LEVFLG,APPLEV
C
C  Local variable / array dictionary:
C
C   Name   Type              Description
C  ------  ----  -----------------------------------------------------
C  ABOV50   R    for a given Year & MY, the vmt in excess of
C                50,000 miles
C  IGER     I    MY pointer into emission slope and intercept
C  IGER50   I    MY pointer into ERB50K
C  ITR1PT   I    MY pointer into TIER1 slope array (TR1SL1).
C  KINK50   I    key to whether or not to use the 50K+ deterioration
C                rates:
C                1 = no:  no 50K kink in the curve or do not reach it
C                2 = yes: have 50K kink & do reach it (mileage exceeds
C                50K)
C  LMY      I    Local model year can be changed from MY depending on
C                the value of NEWFLG.
C  SLOPE1   R    slope of ef curve; if 50K kink, then only for 1st 50K
C                miles
C  SLOPE2   R    slope of kinked ef curve when beyond 50K miles
C  ZPOINT   R    zero mileage point = intercept of ef curve
C
C
C  Notes:
C  BER was created from "BEF" in MOBILE5
C
C  (22-April-1994) @ CSC ljn Request 2-420 added LDGT2 (IV=3) to the LEV's
C        
C  20-July-1994 @CSC pme  request 2-446.  Removed Tier1 mileage cap
C               with special I/M
C  24-August-1994  @CSC-pme request 2-446. Include file BASE12.I was removed. 
C
C October-5-1994 @CSC pme request 2-472. A check for APPLEV variable to 
C see whether or not to apply LEV stds to LDGT2 was included.
C  5-February-1996 @DynTel request 2-613. Warning messages appeared during the
C                  compilation of BEF.FOR was fixed by setting flag IBEFSW
C                  on when BEFUN is called by BEF and off when BEFUN is  
C                  called by GETLEV.  SLOIM1 and SLOIM2 are set in BEFUN only
C                  when flag IBEFSW is on.
C
C
      INCLUDE 'BASEQ1.I'
      INCLUDE 'BASEQ5.I'
      INCLUDE 'BASEQ6.I'
      INCLUDE 'BASE10.I'
      INCLUDE 'BASE11.I'
      INCLUDE 'FLAGS2.I'
      INCLUDE 'IOUCOM.I'
      INCLUDE 'LEVBLK.I'
      INCLUDE 'REGION.I'
      INCLUDE 'TAMEQ2.I'
      INCLUDE 'TIER1S.I'
      INCLUDE 'YEARS4.I'
C
      SLOPE1 = 0.0
      SLOPE2 = 0.0
C
C   redo the kink50 to be sure of it's timing
C
      ABOV50=VMTAGE-F50K
      KINK50=1
      IF(MY.GE.1981.AND.IV.LE.3.AND.ABOV50.GT.0.0) KINK50=2
C
C  Look up base and deterioration rates.  Pointer function IERPTR
C  sets both the arrays to access (block data default or user
C  supplied) and the cell pointers to the coefficients.  KINK50
C  determines if a second slope is needed. Currently, only 1991+
C  LDGT1/2 have 50K+ kinked EFs.
C
C  while if NEWFLG=5, No Clean Air Act (CAA) then set the
C  Model Year (MY) to 1993.( so ALL EF's are Pre-1993)
C  and LEVFLG=1 (FALSE)
C
      LMY=MY
C
      IF (NEWFLG.GE.5.AND.MY.GE.1994) THEN
        LEVFLG=1
        LMY=1993
      END IF
C
C  Set IER LMY (Year Index) to 1999 for LEV's if MY is 1994 or later
C
C     (22-April-1994) @ CSC ljn Request 2-420
C
      IF (LEVFLG.EQ.2.AND.MY.GE.LEVYRS) THEN
        IF (IV.LE.2.OR.(IV.EQ.3.AND.APPLEV.EQ.2)) LMY=1999
      END IF
C
      IGER=IERPTR(LMY,IP,IV,KEYER)
C
C  test for user input
C
      IF(KEYER.EQ.2) GOTO 10
C
      ZPOINT=ERBZML(IGER,IP,IV,IREJN)
      SLOPE1=ERBDR( IGER,IP,IV,IREJN)
C
      IF (KINK50.EQ.1) GOTO 20
      IGER50=1
C
      IF (IV.LE.3) THEN
        IF (LMY.LE.1998) IGER50=LMY-1980
        IF (LMY.GT.1998) IGER50=18
      END IF
C
      SLOPE2=ERB50K(IGER50,IP,IV,IREJN)
      GOTO 20
C
   10 ZPOINT=ERUZML(IGER,IP,IV,IREJN)
      SLOPE1=ERUDR( IGER,IP,IV,IREJN)
      IF(KINK50.EQ.2) SLOPE2=ERU50K(IGER,IP,IV,IREJN)
C
C  save Emission Factor slopes for use in the LEV Calculations
C  ( LEVEF ) used for Nonspecial I/M Cases - Now SLOIM1 & SLOIM1 are
C  set when IBEFSW = 2.
C  (5-February-1996) @ DynTel Request 2-613 If the flag IBEFSW is 2, SLOIM1
C                      and SLOIM2 are set.
C
 20   CONTINUE
C
      IF (IBEFSW.EQ.2) THEN
        SLOIM1 = SLOPE1
        SLOIM2 = SLOPE2
      END IF
C
C  process the "normal" EF's
C
C     (22-April-1994) @ CSC ljn Request 2-420
C
      IF (IV.GE.4) THEN
        IF (KINK50.EQ.1) BER=ZPOINT+SLOPE1*VMTAGE
        IF (KINK50.EQ.2) BER=ZPOINT+SLOPE1*F50K+SLOPE2*ABOV50
        GOTO 99
      END IF
C
C  Tier 1 emission factor for LDGV's MY greater than 1994
C  If special I/M is selected, replace standard slope
C  with special I/M slope (TR1SL1).
C
C     (22-April-1994) @ CSC ljn Request 2-420
C
      IF (IV.LE.3) THEN
C
        IF (IMTIER.EQ.1.AND.LMY.GE.1994) THEN
C
C  get index into tr1sl1
C
          ITR1PT=LMY-1993
          IF (ITR1PT.LT.1) ITR1PT=1
          IF (ITR1PT.GT.6) ITR1PT=6
          SLOPE1=TR1SL1(ITR1PT,IP,IV)
          IF (KINK50.EQ.1) BER=ZPOINT+SLOPE1*VMTAGE
          IF (KINK50.EQ.2) BER=ZPOINT+SLOPE1*F50K+SLOPE1*ABOV50
          GOTO 99
        END IF
C
C  Use standard slopes for non-I/M and standard I/M cases of Tier 1
C  and all non-Tier 1 model years.
C
        IF (IMTIER.EQ.0.OR.LMY.LT.1994) THEN
          IF (KINK50.EQ.1) BER=ZPOINT+SLOPE1*VMTAGE
          IF (KINK50.EQ.2) BER=ZPOINT+SLOPE1*F50K+SLOPE2*ABOV50
        END IF
C
      END IF
C
 99   CONTINUE
      RETURN
      END
