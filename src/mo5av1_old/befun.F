C Version @(#)$Id$ $Source$ $Date$ 

      SUBROUTINE BEFUN(MY,IP,IV,VMTAGE,S1,S2,BER)
C
C  BEFUN returns the basic emission factor Uncorrected for
C  operating mode, CO temperature offset, methane offset, temperature,
C  inspection / maintenance, tampering offset and reformulated fuels.
C
C  Called by BEF, BEFIDL, BEFID3, and GETB92
C
C  Calls IERPTR
C
C  Input on call or from calls:
C
C    parameter list: MY,IP,IV,VMTAGE,S1,S2
C    common blocks:
C    /BASEQ1/ ERBZML,ERBDR,ERB50K
C    /BASEQ5/ ERUZML,ERUDR,ERU50K
C    /BASEQ6/ KEYER,IGER
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
C    /BASE12/ CAPIV1, CAPIV2
C    /LEVBLK/ LEVFLG
C
C  Local variable / array dictionary:
C
C   Name   Type              Description
C  ------  ----  -----------------------------------------------------
C  ABOV50   R    for a given Year & MY, the vmt in excess of
C                50,000 miles
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
C  S1       R    saved variable holding slope of ef curve;
C                if 50K kink, then only for 1st 50K miles
C  SLOPE2   R    slope of kinked ef curve when beyond 50K miles
C  S2       R    saved variable hloding slope of kinked ef curve
C                when beyond 50K miles
C  ZPOINT   R    zero mileage point = intercept of ef curve
C
C
C  Notes:
C  BER was created from "BEF" in MOBILE5
C
C
      INTEGER ALHFLG
C
      COMMON /BASEQ1/ERBZML(25,3,8,2),ERBDR(25,3,8,2),ERB50K(18,3,3,2)
      COMMON /BASEQ5/ERUZML(18,3,8,2),ERUDR(18,3,8,2),ERU50K(18,3,3,2)
      COMMON /BASEQ6/MYGERU(12,2,3,8,2),MAXERU,NUMERU(3,8,2),KEYER,IGER
      COMMON /BASE10/ ZEVTF(25,3),ZEVTOT(3)
      COMMON /BASE12/ CAPIV1, CAPIV2
      COMMON /FLAGS2/ MYMRFG,NEWFLG,IMFLAG,ALHFLG,IMTIER
      COMMON /IOUCOM/ IOUIMD,IOUGEN,IOUREP,IOUERR,IOUASK
      COMMON /LEVBLK/ LEVFLG,IMFLG,LEVYRS
      COMMON /REGION/ FEET(2),IREJN,ALT,INITPR
      COMMON /TAMEQ2/ MYGTAM(3,4),F50K
      COMMON /TIER1S/ TR1SL1(6,3,3)
      COMMON /YEARS4/ IY1941,IY1960,IY1994,IY2020
C
      SLOPE1 = 0.0
      SLOPE2 = 0.0
C
C     VCAPML caps VMTAGE at CAPIV1 OR CAPIV2 mileage
C     capiv1 = 100,000   miles capiv2 = 120,000
C
      VCAPML=VMTAGE
C
C     Tier 1 mileage cap with special I/M
C
      IF(MY.GE.IY1994.AND.IMTIER.EQ.1.AND.IV.LE.2) THEN
         IF(IV.EQ.1) THEN
            IF(VCAPML.GT.CAPIV1) VCAPML = CAPIV1
         ELSE
            IF(VCAPML.GT.CAPIV2) VCAPML = CAPIV2
         ENDIF
      ENDIF
C
C   redo the kink50 to be sure of it's timing
C
      ABOV50=VCAPML-F50K
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
      IF(NEWFLG.GE.5.AND.MY.GE.1994) THEN
         LEVFLG=1
         LMY=1993
      ENDIF
C
C  Set IER LMY (Year Index) to 1999 for LEV's if MY is 1994 or later
C
      IF(LEVFLG.EQ.2.AND.MY.GE.LEVYRS)THEN
         IF(IV.LE.2) LMY=1999
      ENDIF
C
C  Note "IGER" is passed in through Global Common while
C  "IDLGER" is returned (for use in BEFID? HDGV CO idle
C  MY greater than or equal to 1984)  So you get two pointers
C  variables from One call to IERPTR
C
      IDLGER=1
      CALL IERPTR(LMY,IP,IV,IDLGER)
C
C  test for user input
C
      IF(KEYER.EQ.2) GOTO 10
C
      ZPOINT=ERBZML(IGER,IP,IV,IREJN)
      SLOPE1=ERBDR( IGER,IP,IV,IREJN)
C
      IF(KINK50.EQ.1) GOTO 20
      IGER50=1
      IF(IV.LE.3) THEN
         IF(LMY.LE.1998) IGER50=LMY-1980
         IF(LMY.GT.1998) IGER50=18
      ENDIF
      SLOPE2=ERB50K(IGER50,IP,IV,IREJN)
      GOTO 20
C
   10 ZPOINT=ERUZML(IGER,IP,IV,IREJN)
      SLOPE1=ERUDR( IGER,IP,IV,IREJN)
      IF(KINK50.EQ.2) SLOPE2=ERU50K(IGER,IP,IV,IREJN)
C
C  save Emission Factor slopes for use in the LEV Calculations
C  ( LEVEF ) used for Nonspecial I/M Cases - Now S1 & S1 are
C  returned in BUFUN
C
   20 S1 = SLOPE1
      S2 = SLOPE2
C
C  process the "normal" EF's
C
      IF(IV.GE.3) THEN
         IF(KINK50.EQ.1) BER=ZPOINT+SLOPE1*VCAPML
         IF(KINK50.EQ.2) BER=ZPOINT+SLOPE1*F50K+SLOPE2*ABOV50
         GOTO 99
      ENDIF
C
C  Tier 1 emission factor for LDGV's MY greater than 1994
C  If special I/M is selected, replace standard slope
C  with special I/M slope (TR1SL1).
C
      IF(IV.LE.2) THEN
        IF(IMTIER.EQ.1.AND.LMY.GE.1994) THEN
C
C  get index into tr1sl1
C
          ITR1PT=LMY-1993
          IF(ITR1PT.LT.1) ITR1PT=1
          IF(ITR1PT.GT.6) ITR1PT=6
          SLOPE1=TR1SL1(ITR1PT,IP,IV)
          IF(KINK50.EQ.1) BER=ZPOINT+SLOPE1*VCAPML
          IF(KINK50.EQ.2) BER=ZPOINT+SLOPE1*F50K+SLOPE1*ABOV50
          GOTO 99
        ENDIF
C
C  Use standard slopes for non-I/M and standard I/M cases of Tier 1
C  and all non-Tier 1 model years.
C
        IF(IMTIER.EQ.0.OR.LMY.LT.1994) THEN
           IF(KINK50.EQ.1) BER=ZPOINT+SLOPE1*VCAPML
           IF(KINK50.EQ.2) BER=ZPOINT+SLOPE1*F50K+SLOPE2*ABOV50
        ENDIF
      ENDIF
C
   99 RETURN
      END
