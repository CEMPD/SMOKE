C Version "@(#)$Id$ $Source$ $Date$ 

      FUNCTION PCLEFT(MY,ICY,IP,IV)
C
C  PCLEFT determines the basic emission factor multiplicative
C  adjustment ( = percentage left of BEF = PCLEFT ) for the
C  effects of an inspection / maintanance (I/M) program, after
C  checking whether or not I/M applies to the factor being
C  computed by the calling function.
C
C  Called by BEF and GETLEV.
C
C  Calls IMPTR, IMOVLP, and OVLP2
C
C  Input on call:
C
C    parameter list: MY,ICY,IP,IV
C    common blocks:
C    /FLAGS2/ IMFLAG
C    /ICR4VA/ CRD4VA
C    /ICR4VB/ CRD4VB
C    /IMPAR1/ ISTRIN,MODYR,WAIVER,CRIM
C    /IMPAR2/ ILDT
C    /IMPAR5/ CRHDGV,DISCNT
C    /IMPAR6/ IFREQ,INTYP
C    /IM12HC/ CR12HC
C    /IM12CO/ CR12CO
C    /TTCINF/ TTCDAT,ITTCFG,TTCEFF
C
C  Output on return:
C
C    function: PCLEFT
C
C  Local variable / array dictionary:
C
C   Name   Type              Description
C  ------  ----  ------------------------------------------------
C  AGE1ST   I    age of the vehicle at first inspection
C  BCLEFT   R    temporary variable for holding function return
C  BY       I    benefit year for technology 1 or 2 vehicle
C  DISHYB   R    Retest Based Hybrid Discount for Pre-1981 vehicles &
C                TTC credits
C  IBY      I    benefit year for technology 4 plus vehicle
C  IREM     I    remainder = stringency - greatest multiple of
C                10 < stringency
C  ISTRN    I    stringency index into the technology 1 or 2
C                credits array
C  ITECH    I    pointer to technology and vehicle type
C  LICYIM   I    I/M program start year that applies
C  MAXIMP   R    maximum number of I/M programs
C  NNOVLP   I    number of I/M program in effect 0-3
C                0) No program
C                1) 1st program
C                2) 2nd program (no overlap with the 1st)
C                3) biennial, 2nd program (and overlap with
C                   the 1st)
C  PCLFTN   R    percentage correction factor for new I/M
C  PCLFTO   R    percentage correction factor for old I/M
C  PCREDN   R    I/M credit as determined from credit deck
C                and is appropriate to new I/M program
C  PCREDO   R    I/M credit as determined from credit deck
C                and is appropriate to old I/M program
C  REM      R    IREM converted to REAL value
C  WAIV     R    waiver rate used on this call, depends on MY
C
C
C
C  Notes:
C
C  Non-TECH4+ vehicles use 19th value for 20-24 benefit year.
C  TECH4+ vehicles have arrays expanded to MAXYRS-1.
C  PCLEFT was modified to include the IM240 check, and NOx
C  credits.
C  PCLEFT was modified to fix the index into I/M credits.
C  PCLEFT was modified to fix the biennial phase-in for
C  MOBILE5 corrected.  PCLEFT was modified again and turned into
C  a driver routine to decide to return the value of PCLFT1 or
C  PCLFT2 depending or ICY.
C  PCLFT2 was modified to include the IM240 check, and NOx credits.
C  PCLFT2 was modified to fix the index into I/M credits
C  For calendar years less than or equal to 1990 the old PCLEFT
C  code is maintained in PCLFT2.  See PCLFT2 for the new code.
C  PCLFT1 was modified to include the IM240 check, and NOx
C  credits.
C  PCLFT1 was modified to fix the index into I/M credits.
C  PCLFT1 was modified to fix the biennial phase-in
C  2-10-93.
C  June-15-1993 @ ARC-bk Subtask 244 the common block IMPAR5
C  was modified, the HDGV I/M credit was expanded to include
C  an initial value for NOx.
C  September-16-1993 @ ARC-bsg Subtask 311 Regrouped PCLEFT, PCLFT1 
C  and PCLFT2 into one module called PCLEFT.
C
C  May-4-1994 @ CSC-tmm Subtask 2-422 
C    Discovered bug in the way the call to OVLP2 was handled.
C    The function OVLP2 returns an integer value. Here in PCLEFT
C    the compiler assumed that OVLP2 returned a real, because of
C    implicit typing. To correct the problem, a line was added to
C    declare OVLP2 to be an integer. This was the only change needed.
C
C
C  April-1995 @ CSC-zp Req. 2-495 "Add MOB5a-H to MOBILE5b"
C  Added:
C  Retest Based Hybrid Discount for Pre-1981 vehicles & TTC credits.
C  Adjusting the Tech 1, 2 & 3 credits for Retest Based Hybrid Programs.
C  Adjusting the I/M credit for Technician Training & Certification
C  Adjusting the Tech 4 TTC credits for Retest Based Hybrid Programs.
C  Added warning 166 when I/M and TTC credit exeed 100%.
C   Reset credit to 100%
C
C  1 Feb 96 @DynTel-MLA 2-611  Added the IP dimension to array DISCNT. This
C     array may now contain either user-supplied values or default values.
C
C  1 Apr 96 @DynTel-MLA 2-611 (part 2) I/M credits for test-and-repair
C     programs, which were disounted by 50% relative to test-only programs,
C     were not subject to further discounts for waivers unless the waiver
C     rate was greater than 50%. Now, for waiver rates less than 50% the
C     waiver discount is phased in as program effectiveness for test-and-
C     repair programs varies between 50% (the old default) and 100%.
C 
      INCLUDE 'FLAGS2.I'
      INCLUDE 'IM12HC.I'
      INCLUDE 'IM12CO.I'
      INCLUDE 'ICR4VA.I'
      INCLUDE 'ICR4VB.I'
      INCLUDE 'IMPAR1.I'
      INCLUDE 'IMPAR2.I'
      INCLUDE 'IMPAR5.I'
      INCLUDE 'IMPAR6.I'
      INCLUDE 'MAXIMA.I'
      INCLUDE 'TTCINF.I'
C       
      INTEGER BY,AGE1ST
      DIMENSION WAIV(2)
      DIMENSION DISHYB(24,3)
C
C  Retest Based Hybrid Discount for Pre-1981 vehicles & TTC credits.
C
      DATA DISHYB /
     H 0.002,0.013,0.027,0.036,0.035,0.044,0.054,0.054,0.051,0.048,
     H 0.044,0.038,0.034,0.030,0.026,0.022,0.020,0.018,0.017,0.015,
     H 0.015,0.013,0.013,0.012,
     C 0.008,0.033,0.056,0.057,0.047,0.051,0.055,0.051,0.044,0.042,
     C 0.038,0.031,0.027,0.024,0.022,0.019,0.017,0.016,0.015,0.013,
     C 0.013,0.012,0.011,0.010,
     N 0.003,0.010,0.016,0.011,0.010,0.022,0.033,0.035,0.033,0.031,
     N 0.029,0.027,0.025,0.023,0.021,0.019,0.017,0.015,0.012,0.011,
     N 0.010,0.009,0.009,0.008 /
C      
C
C  May-4-1994 @ CSC-tmm Subtask 2-422 
      INTEGER OVLP2
C
C
C  Initialize BCLEFT to "no reduction" value = 1.0.
C  Initialize PCRED to "no credit" value 0.0.
C
      BCLEFT=1.0
      PCREDN=0.0
      PCREDO=0.0
      PCLFTN=1.0
      PCLFTO=1.0
C
C  Discover the input data: has to be one of four cases:
C  "0" No I/M program is in effect
C  "1" The first I/M program is in effect
C  "2" The second I/M program is in effect and no "overlap"
C  "3" New special case of "overlap".  I/M program #2 is in
C      in effect for this model year, vehicle, and evaluation
C      year as well as the first I/M program.
C
      NNOVLP=IMOVLP(MY,IV,ICY)
      IF (NNOVLP.EQ.0) GOTO 99
      MAXIMP=2
      IF (NNOVLP.EQ.1) MAXIMP=1
      LICYIM=OVLP2(MY,IV,ICY)
C
      DO 70 IMPNUM=1,MAXIMP
C
C  No I/M Program Credit is applied when the cars are
C  "zero" years old
C
        IF (ICY.EQ.MY) GOTO 70
C
C  Assign waiver rates
C  Control WAIV (waver rates) for INTYP (I/M program inspection
C  type), where 2 = I/M test and repair (Computerized),
C  & 3 = test and repair (Manual), if INTYP is not a 1 = test
C  only I/M and the WAIV is LE 50%, WAIV is set to zero.
C
        IF (MY.LE.1980) WAIV(IMPNUM)=WAIVER(1,IMPNUM)
        IF (MY.GT.1980) WAIV(IMPNUM)=WAIVER(2,IMPNUM)
C
        IF (WAIV(IMPNUM).LE.0.50.AND.INTYP(IMPNUM).GT.1) THEN
C
          IF (DISCNT(IP, INTYP(IMPNUM)).LT.0.50) THEN
            WAIV(IMPNUM) = 0.0
          ELSE
            WAIV(IMPNUM)=WAIV(IMPNUM)*(2*DISCNT(IP,INTYP(IMPNUM))-1)
          ENDIF
C
        ENDIF
C
        IF (IV.EQ.4) GOTO 50
C
C  Selecting I/M reduction for LDGV or LDGT.  Several parameters
C  must be set.
C
C  Determine technology by model year and vehicle type.
C
        ITECH=IMPTR(MY,IV)
C
C  Branch on technology type: TECH 1 & 2 form 1 group, TECH 4+
C  another.
C
        IF (ITECH.GE.4) GOTO 40
C
C  Find the benefit year:
C  Limit value according to technology type.
C
        BY=ICY-LICYIM
        IF (MY.GT.LICYIM) BY=ICY-MY
        IF (ITECH.LE.3.AND.BY.GT.19) BY=19
C
C  Find the age of the vehicle at first inspection.
C  Limit value according to technology type.
C
        AGE1ST=1
        IF (MY.LT.LICYIM) AGE1ST=LICYIM-MY+1
C
        IF (ITECH.LE.3) THEN
          IF (AGE1ST.GT.19) AGE1ST=19
          IF (AGE1ST+BY.GT.19) BY=20-AGE1ST
        ENDIF
C
C  For now, continue MOBILE3 policy of using technology 2 credits for
C  ITECH=3.
C
        IF (ITECH.EQ.3) ITECH=2
C
C  There are no NOx credits for Tech 1 & 2.
C
        IF (IP.EQ.3) GOTO 70
C
C  Select correct I/M credits for TECH 1 & 2.  The same credits
C  array is used for LDGV and LDGT, but the ITECH mygs differ, so
C  that the same MY may yield a different credit for LDGV than
C  for LDGT.
C
C  Interpolate between 10, 20, 30, 40 & 50% stringency.
C
        IREM=ISTRIN(IMPNUM)-(ISTRIN(IMPNUM)/10)*10
        REM=IREM*.1
        ISTRN=(ISTRIN(IMPNUM)-IREM)/10
C
        NFREQ=IFREQ(IMPNUM)
C
        IF (NFREQ.EQ.1) THEN
C
          IF (IP.EQ.1) THEN
            PCRED=CR12HC(BY,AGE1ST,ISTRN,ITECH)
            IF (ISTRN.LT.5.AND.IREM.GT.0)  PCRED=
     * (CR12HC(BY,AGE1ST,ISTRN+1,ITECH)-PCRED)*REM+PCRED
          ELSE IF (IP.EQ.2) THEN
            PCRED=CR12CO(BY,AGE1ST,ISTRN,ITECH)
            IF (ISTRN.LT.5.AND.IREM.GT.0)  PCRED=
     * (CR12CO(BY,AGE1ST,ISTRN+1,ITECH)-PCRED)*REM+PCRED
          ENDIF
C
        ELSE IF (NFREQ.EQ.2) THEN
C
          IF (IP.EQ.1) THEN
            PCRED=CR12HC(20-BY,21-AGE1ST,ISTRN,ITECH)
            IF (ISTRN.LT.5.AND.IREM.GT.0)  PCRED=
     * (CR12HC(20-BY,21-AGE1ST,ISTRN+1,ITECH)-PCRED)*REM+PCRED
C
          ELSE IF(IP.EQ.2) THEN
            PCRED=CR12CO(20-BY,21-AGE1ST,ISTRN,ITECH)
            IF (ISTRN.LT.5.AND.IREM.GT.0)  PCRED=
     * (CR12CO(20-BY,21-AGE1ST,ISTRN+1,ITECH)-PCRED)*REM+PCRED
C
          ENDIF
C
        ENDIF
C
C  Adjust the Tech 1, 2 & 3 credits for Retest Based Hybrid Programs.
C
        IF (IHYB(IMPNUM).EQ.2) THEN
          PCRED=PCRED-DISHYB(AGE1ST+BY-1,IP)
          IF (PCRED .LT. 0.0) PCRED=0.0
        ENDIF
      
C
        GOTO 60
C
C  Select I/M credits for TECH 4+ vehicles.  For TECH 4+, there
C  are separate credits arrays, one for annual and one for
C  biennial inspection programs.
C
   40   IBY=ICY-MY
C
C  If no CAAA, use 1993 I/M credits for 1993+ MYs
C
        IF (NEWFLG.GE.5.AND.MY.GT.1993) ITECH=IMPTR(1993,IV)
C
C  CRD4VA contains the I/M credits for first I/M description
C  CRD4VB contains the I/M credits for second I/M description
C
        IF (IMPNUM.EQ.1) PCRED=CRD4VA(IBY,IP,ITECH-3)
        IF (IMPNUM.EQ.2) PCRED=CRD4VB(IBY,IP,ITECH-3)
C  
C
C  Adjust the I/M credit for Technician Training & Certification
C  Hybrid I/M PCRED has already been reduced by DISHYB, add DISHYB to
C  PCRED to simulate full credit. Determine ratio of DISHYB to full
C  credit. DISHYB is a delta value.  The absolute value of PCRED has
C  to be taken.
C
        IF (ITTCFG .GT. 1) THEN
          TTCADJ=DISHYB(IBY,IP)/(ABS(PCRED)+DISHYB(IBY,IP))
C
C  Adjust the Tech 4 TTC credits for Retest Based Hybrid Programs.
C
          TTCRED=TTCDAT(IBY,IP,ITECH-3)
          IF (IHYB(IMPNUM).EQ.2) TTCRED=TTCRED*(1.0-TTCADJ)
          IF (TTCRED.LT.0.0) TTCRED=0.0
          TTCRED=TTCRED*TTCEFF
          PCRED=PCRED+TTCRED
C
          IF (PCRED .GT. 1.0) THEN
            PCRED = 1.0
            CALL QUITER(0.,0,166,INERR)
          ENDIF
C
        ENDIF
C  
        GOTO 60
C
C  Assign HDGV I/M credit
C  and there are no TIER1 HDGV I/M credits
C
   50   CONTINUE
        NFREQ=IFREQ(IMPNUM)
        PCRED=CRHDGV(IP,NFREQ)
C
C  Inside IMPNUM loop assign appropriate credit to the
C  old I/M and the new I/M programs.
C
   60   CONTINUE
        IF (IMPNUM.EQ.1) PCREDO=PCRED
        IF (IMPNUM.EQ.2) PCREDN=PCRED
C
   70 CONTINUE
C
C  Calculate the correction factor only if there is one I/M
C  program.
C
      IF (NNOVLP.EQ.1) THEN
        PCLFTO=1.0-(PCREDO*(1.0-WAIV(1))*ENFORC(CRIM(NNOVLP),1)*
     * DISCNT(IP, INTYP(NNOVLP)))
        IF (ICY.EQ.ICYIM(NNOVLP)+1.
     * .AND.IFREQ(NNOVLP).EQ.2) PCLFTO=PCLFTO*0.5+0.5
        BCLEFT=PCLFTO
      ENDIF
C
C  Calculate the correction factor only if there are two I/M
C  programs and the second is in effect but there is no "overlap"
C
      IF (NNOVLP.EQ.2) THEN
        PCLFTN=1.0-(PCREDN*(1.0-WAIV(2))*ENFORC(CRIM(NNOVLP),1)*
     * DISCNT(IP, INTYP(NNOVLP)))
        IF (ICY.EQ.ICYIM(NNOVLP)+1.
     * .AND.IFREQ(NNOVLP).EQ.2) PCLFTN=PCLFTN*0.5+0.5
        BCLEFT=PCLFTN
      ENDIF
C
C  Calculate the correction factor only if there are two I/M
C  programs the second is controlling the first one and the
C  second is biennial.
C
      IF (NNOVLP.EQ.3) THEN
        PCLFTN=1.0-(PCREDN*(1.0-WAIV(2))*ENFORC(CRIM(2),1)*
     * DISCNT(IP, INTYP(2)))
C
        IF (ICY.EQ.ICYIM(2)+1) THEN
C
C  Second program controls, ignore first program
C
          IF (ICYIM(1).EQ.ICYIM(2)) THEN
            BCLEFT=PCLFTN*0.5+0.5
          ELSE
            PCLFTO=1.0-(PCREDO*(1.0-WAIV(1))*ENFORC(CRIM(1),1)*
     * DISCNT(IP,INTYP(1)))
            BCLEFT=PCLFTN+PCLFTO
            BCLEFT=BCLEFT*0.5
          ENDIF
C
        ELSE
          BCLEFT=PCLFTN
        ENDIF
C
      ENDIF
C
   99 PCLEFT=BCLEFT
C
      RETURN
      END
