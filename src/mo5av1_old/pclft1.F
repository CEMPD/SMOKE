C Version @(#)$Id$ $Source$ $Date$ 

      FUNCTION PCLFT1(MY,ICY,IP,IV)
C
C  PCLFT1 determines the basic emission factor multiplicative
C  adjustment ( = percentage left of BEF = PCLFT1 ) for the
C  effects of an inspection / maintanance (I/M) program, after
C  checking whether or not I/M applies to the factor being
C  computed by the calling function.
C
C  Called by PCLEFT
C
C  Calls IMPTR, IMOVLP
C
C  Input on call:
C
C    parameter list: MY,ICY,IP,IV
C    common blocks:
C    /FLAGS2/ IMFLAG
C    /ICR4VA/ CRD4VA
C    /ICR4VB/ CRD4VB
C    /IMPAR1/ ICYIM,ISTRIN,MODYR,WAIVER,CRIM
C    /IMPAR2/ ILDT
C    /IMPAR5/ CRHDGV,DISCNT
C    /IMPAR6/ IFREQ,INTYP
C    /IM12HC/ CR12HC
C    /IM12CO/ CR12CO
C
C  Output on return:
C
C    function: PCLFT1
C
C  Local variable / array dictionary:
C
C   Name   Type              Description
C  ------  ----  ------------------------------------------------
C  AGE1ST   I    age of the vehicle at first inspection
C  BCLEFT   R    temporary variable for holding function return
C  BY       I    benefit year for technology 1 or 2 vehicle
C  IBY      I    benefit year for technology 4 plus vehicle
C  IREM     I    remainder = stringency - greatest multiple of
C                10 < stringency
C  ISTRN    I    stringency index into the technology 1 or 2
C                credits array
C  ITECH    I    pointer to technology and vehicle type
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
C  Notes:
C
C  Non-TECH4+ vehicles use 19th value for 20-24 benefit year.
C  TECH4+ vehicles have arrays expanded to MAXYRS-1.
C  PCLFT1 was modified to include the IM240 check, and NOx
C  credits.
C  PCLFT1 was modified to fix the index into I/M credits.
C  PCLFT1 was modified to fix the biennial phase-in
C  2-10-93.
C
C
      INTEGER ALHFLG,BY,AGE1ST
      DIMENSION WAIV(2)
C
      COMMON /FLAGS2/ MYMRFG,NEWFLG,IMFLAG,ALHFLG,IMTIER
      COMMON /IM12HC/ CR12HC(19,20,5,2)
      COMMON /IM12CO/ CR12CO(19,20,5,2)
      COMMON /ICR4VA/ CRD4VA(24,3,18)
      COMMON /ICR4VB/ CRD4VB(24,3,18)
      COMMON /IMPAR1/ ICYIM(2),ISTRIN(2),MODYR(2,2),WAIVER(2,2),CRIM(2)
      COMMON /IMPAR2/ ILDT(4,2),ITEST(2)
      COMMON /IMPAR5/ MYGIM2(3,3),LBIM4P,CRHDGV(2,2),DISCNT(3)
      COMMON /IMPAR6/ IFREQ(2),INTYP(2),IFILET(2)
      COMMON /MAXIMA/ MAXVEH,MAXLTW,MAXPOL,MAXREG,MAXYRS
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
      IF(NNOVLP.EQ.0) GOTO 99
      MAXIMP=2
      IF(NNOVLP.EQ.1) MAXIMP=1
      DO 70 IMPNUM=1,MAXIMP
C
C  No I/M Program Credit is applied when the cars are
C  "zero" years old
C
        IF(ICY.EQ.MY) GOTO 70
C
C  Assign waiver rates
C  Control WAIV (waver rates) for INTYP (I/M program inspection
C  type), where 2 = I/M test and repair (Computerized),
C  & 3 = test and repair (Manual), if INTYP is not a 1 = test
C  only I/M and the WAIV is LE 50%, WAIV is set to zero.
C
        IF(MY.LE.1980) WAIV(IMPNUM)=WAIVER(1,IMPNUM)
        IF(MY.GT.1980) WAIV(IMPNUM)=WAIVER(2,IMPNUM)
        IF(WAIV(IMPNUM).LE.0.50.AND.INTYP(IMPNUM).GT.1)
     *     WAIV(IMPNUM)=0.0
C
        IF(IV.EQ.4) GOTO 50
C
C  Selecting I/M reduction for LDGV or LDGT.  Several parameters
C  must be set.
C
C  Determine technology by model year and vehicle type.
C
        ITECH=IMPTR(MY,IV)
C
C  Find the benefit year:
C  Limit value according to technology type.
C
        BY=ICY-ICYIM(IMPNUM)
        IF(MY.GT.ICYIM(IMPNUM)) BY=ICY-MY
        IF(ITECH.LE.3.AND.BY.GT.19) BY=19
        IF(ITECH.GE.4.AND.BY.GT.MAXYRS-1) BY=MAXYRS-1
C
C  Find the age of the vehicle at first inspection.
C  Limit value according to technology type.
C
        AGE1ST=1
        IF(MY.LT.ICYIM(IMPNUM)) AGE1ST=ICYIM(IMPNUM)-MY+1
        IF(ITECH.LE.3) THEN
          IF(AGE1ST.GT.19) AGE1ST=19
          IF(AGE1ST+BY.GT.19) BY=20-AGE1ST
        ENDIF
C
C  For now, continue MOBILE3 policy of using technology 2 credits for
C  ITECH=3.
C
        IF(ITECH.EQ.3) ITECH=2
C
C  Branch on technology type: TECH 1 & 2 form 1 group, TECH 4+
C  another.
C
        IF(ITECH.GE.4) GOTO 40
C
C  There are no NOx credits for Tech 1 & 2.
C
        IF(IP.EQ.3) GOTO 70
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
        IF(NFREQ.EQ.1) THEN
          IF(IP.EQ.1) THEN
            PCRED=CR12HC(BY,AGE1ST,ISTRN,ITECH)
            IF(ISTRN.LT.5.AND.IREM.GT.0)  PCRED=
     *        (CR12HC(BY,AGE1ST,ISTRN+1,ITECH)-PCRED)*REM+PCRED
          ELSE IF(IP.EQ.2) THEN
            PCRED=CR12CO(BY,AGE1ST,ISTRN,ITECH)
            IF(ISTRN.LT.5.AND.IREM.GT.0)  PCRED=
     *        (CR12CO(BY,AGE1ST,ISTRN+1,ITECH)-PCRED)*REM+PCRED
          ENDIF
C
        ELSE IF(NFREQ.EQ.2) THEN
          IF(IP.EQ.1) THEN
            PCRED=CR12HC(20-BY,21-AGE1ST,ISTRN,ITECH)
            IF(ISTRN.LT.5.AND.IREM.GT.0)  PCRED=
     *        (CR12HC(20-BY,21-AGE1ST,ISTRN+1,ITECH)-PCRED)*REM+
     *        PCRED
          ELSE IF(IP.EQ.2) THEN
            PCRED=CR12CO(20-BY,21-AGE1ST,ISTRN,ITECH)
            IF(ISTRN.LT.5.AND.IREM.GT.0)  PCRED=
     *        (CR12CO(20-BY,21-AGE1ST,ISTRN+1,ITECH)-PCRED)*REM+
     *        PCRED
          ENDIF
C
        ENDIF
C
        GOTO 60
C
C  Select I/M credits for TECH 4+ vehicles.  For TECH 4+, there
C  are separate credits arrays, one for annual and one for
C  biennial inspection programs.
C
   40   IBY=BY+AGE1ST-1
C
C  If no CAAA, use 1993 I/M credits for 1993+ MYs
C
        IF(NEWFLG.GE.5.AND.MY.GT.1993) ITECH=IMPTR(1993,IV)
C
C  CRD4VA contains the I/M credits for first I/M description
C  CRD4VB contains the I/M credits for second I/M description
C
        IF(IMPNUM.EQ.1) PCRED=CRD4VA(IBY,IP,ITECH-3)
        IF(IMPNUM.EQ.2) PCRED=CRD4VB(IBY,IP,ITECH-3)
C
        GOTO 60
C
C  Assign HDGV I/M credit
C  and there are no TIER1 HDGV I/M credits
C
   50   NFREQ=IFREQ(IMPNUM)
        IF(IP.NE.3) PCRED=CRHDGV(IP,NFREQ)
C
C  Inside IMPNUM loop assign appropriate credit to the
C  old I/M and the new I/M programs.
C
   60   IF(IMPNUM.EQ.1) PCREDO=PCRED
        IF(IMPNUM.EQ.2) PCREDN=PCRED
   70 CONTINUE
C
C  Calaculate the correction factor only if there is one I/M
C  program.
C
      IF(NNOVLP.EQ.1) THEN
        PCLFTO=1.0-(PCREDO*(1.0-WAIV(1))*ENFORC(CRIM(NNOVLP),1)*
     *  DISCNT(INTYP(NNOVLP)))
        IF(ICY.EQ.ICYIM(NNOVLP)+1.
     *   .AND.IFREQ(NNOVLP).EQ.2) PCLFTO=PCLFTO*0.5+0.5
        BCLEFT=PCLFTO
      ENDIF
C
C  Calculate the correction factor only if there are two I/M
C  programs and the second is in effect but there is no "overlap"
C
      IF(NNOVLP.EQ.2) THEN
        PCLFTN=1.0-(PCREDN*(1.0-WAIV(2))*ENFORC(CRIM(NNOVLP),1)*
     *  DISCNT(INTYP(NNOVLP)))
        IF(ICY.EQ.ICYIM(NNOVLP)+1.
     *   .AND.IFREQ(NNOVLP).EQ.2) PCLFTN=PCLFTN*0.5+0.5
        BCLEFT=PCLFTN
      ENDIF
C
C  Calculate the correction factor only if there are two I/M
C  programs the second is controlling the first one and the
C  second is biennial.
C
      IF(NNOVLP.EQ.3) THEN
        PCLFTN=1.0-(PCREDN*(1.0-WAIV(2))*ENFORC(CRIM(2),1)*
     *  DISCNT(INTYP(2)))
        IF(ICY.EQ.ICYIM(2)+1) THEN
C
C  Second program controls, ignore first program
C
          IF(ICYIM(1).EQ.ICYIM(2)) THEN
            BCLEFT=PCLFTN*0.5+0.5
          ELSE
            PCLFTO=1.0-(PCREDO*(1.0-WAIV(1))*ENFORC(CRIM(1),1)*
     *      DISCNT(INTYP(1)))
            BCLEFT=PCLFTN+PCLFTO
            BCLEFT=BCLEFT*0.5
          ENDIF
        ELSE
          BCLEFT=PCLFTN
        ENDIF
      ENDIF
C
   99 PCLFT1=BCLEFT
C
      RETURN
      END
