C Version @(#)$Id$ $Source$ $Date$ 

      FUNCTION PCLFT2(MY,ICY,IP,IV)
C
C  PCLFT2 determines the basic emission factor multiplicative
C  adjustment ( = percentage left of BEF = PCLFT2 ) for the effects
C  of an inspection / maintanance (I/M) program, after checking
C  whether or not I/M applies to the factor being computed by the
C  calling function.
C
C  Called by PCLEFT
C
C  Calls IMPTR, MYRIMT
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
C    function: PCLFT2
C
C  Local variable / array dictionary:
C
C   Name   Type              Description
C  ------  ----  -----------------------------------------------------
C  AGE1ST   I    age of the vehicle at first inspection
C  BY       I    benefit year for technology 1 or 2 vehicle
C  IBY      I    benefit year for technology 4 plus vehicle
C  IREM     I    remainder = stringency - greatest multiple of
C                10 < stringency
C  ISTRN    I    stringency index into technology 1 or 2 credits array
C  REM      R    IREM converted to REAL value
C  WAIV     R    waiver rate used on this call, depends on MY
C  IMPNUM   I    see program header.
C
C
C  Notes:
C
C  Non-TECH4+ vehicles use 19th value for 20-24 benefit year.
C  TECH4+ vehicles have arrays expanded to MAXYRS-1.
C  PCLFT2 was modified to include the IM240 check, and NOx credits.
C  PCLFT2 was modified to fix the index into I/M credits
C  For calendar years less than or equal to 1990 the old PCLEFT
C  code is maintained in PCLFT2.  See PCLFT2 for the new code.
C
C
      INTEGER ALHFLG
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
      INTEGER BY,AGE1ST
C
C  Initialize PCLFT2 to "no reduction" value = 1.0.
C
      PCLFT2=1.0
      PCRED=0.0
C
C  If factor being calculated by BEF is not covered by I/M, then RETURN.
C
      IMPNUM=MYRIMT(MY,IV,ICY)
      IF((IMPNUM.EQ.0).OR.(ICY.EQ.MY)) GOTO 99
C
C  Assign waiver rates
C  Control WAIV (waver rates) for INTYP (I/M program inspection
C  type), where 2 = I/M test and repair (Computerized),
C  & 3 = test and repair (Manual), if INTYP is not a 1 = test
C  only I/M and the WAIV is LE 50%, WAIV is set to zero.
C
      IF(MY.LE.1980) WAIV=WAIVER(1,IMPNUM)
      IF(MY.GT.1980) WAIV=WAIVER(2,IMPNUM)
      IF(WAIV.LE.0.50.AND.INTYP(IMPNUM).GT.1) WAIV=0.0
C
      IF(IV.EQ.4) GOTO 50
C
C  Selecting I/M reduction for LDGV or LDGT.  Several parameters must be set.
C  be set.
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
C  Branch on technology type: TECH 1 & 2 form 1 group, TECH 4+ another.
C
      IF(ITECH.GE.4) GOTO 40
C
C  There are no NOx credits for Tech 1 & 2.
C
      IF(IP.EQ.3) GOTO 99
C
C  Select correct I/M credits for TECH 1 & 2.  The same credits array
C  is used for LDGV and LDGT, but the ITECH mygs differ, so that the
C  same MY may yield a different credit for LDGV than for LDGT.
C
C  Interpolate between 10,20,30,40 & 30% stringency.
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
     *      (CR12HC(BY,AGE1ST,ISTRN+1,ITECH)-PCRED)*REM+PCRED
        ELSE IF(IP.EQ.2) THEN
          PCRED=CR12CO(BY,AGE1ST,ISTRN,ITECH)
          IF(ISTRN.LT.5.AND.IREM.GT.0)  PCRED=
     *      (CR12CO(BY,AGE1ST,ISTRN+1,ITECH)-PCRED)*REM+PCRED
        ENDIF
C
      ELSE IF(NFREQ.EQ.2) THEN
        IF(IP.EQ.1) THEN
          PCRED=CR12HC(20-BY,21-AGE1ST,ISTRN,ITECH)
          IF(ISTRN.LT.5.AND.IREM.GT.0)  PCRED=
     *      (CR12HC(20-BY,21-AGE1ST,ISTRN+1,ITECH)-PCRED)*REM+PCRED
        ELSE IF(IP.EQ.2) THEN
          PCRED=CR12CO(20-BY,21-AGE1ST,ISTRN,ITECH)
          IF(ISTRN.LT.5.AND.IREM.GT.0)  PCRED=
     *      (CR12CO(20-BY,21-AGE1ST,ISTRN+1,ITECH)-PCRED)*REM+PCRED
        ENDIF
C
      ENDIF
C
      GOTO 60
C
C  Select I/M credits for TECH 4+ vehicles.  For TECH 4+, there are
C  separate credits arrays, one for annual and one for biennial
C  inspection programs.
C
   40 IBY=BY+AGE1ST-1
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
   50 NFREQ=IFREQ(IMPNUM)
      IF(IP.NE.3) PCRED=CRHDGV(IP,NFREQ)
C
C  Set I/M benefit adjusted for waivers and enforcement and discount
C
   60 PCLFT2=1.0-(PCRED*(1.0-WAIV)*ENFORC(CRIM(IMPNUM),1)
     *      *DISCNT(INTYP(IMPNUM)))
C
   99 RETURN
      END
