      REAL FUNCTION IMCP81(MY,ICY,IP,IV) 
C 
C  IMCP81 determines the basic emission factor multiplicative 
C  adjustment for the effects of an I/M program. 
C 
C  This subroutine contains the code that was in PCLEFT in earlier 
C  versions of Mobile. In Mobile6 the maximum number of I/M programs 
C  was expanded from 2 to 5 and it was determined that it should not 
C  be legal for I/M programs to overlap. That is, no more than one 
C  I/M program should cover any model year range for any vehicle type 
C  in a given calendar year. The code in PCLEFT that made decisions 
C  about multiple I/M programs was deleted when that code was moved 
C  to this subroutine. 
C 
C  Called by PCLEFT. 
C 
C  Calls FINDIM, IMPTR, IMPSYR, QUITER 
C 
C  Changes: (Last change first) 
C 
C  20 Sep 00 @EPA- BG Added constant IEXHAUST=1 for call to FINDIM 
C 15 Sep 00 AIR Task 08: Updated FINDIM calls to select exhaust I/M. 
C  2 Jun 99 @DynTel-MLA 2-694  Removed EXEMPTAGE test (done now in FINDIM). 
C 10 Feb 99 @DynTel-MLA 2-663  This is a new routine for Mobile6. 
C 
C  Input on call: 
C 
C    argument list: MY,ICY,IP,IV 
C    common blocks: 
C    /IM12HC/ CR12HC 
C    /IM12CO/ CR12CO 
C    /IMPAR1/ ISTRIN,MODYR,WAVPRE81,WAV81PLUS,CRIM,EXEMPTAGE 
C    /IMPAR5/ CRHDGV,DISCNT 
C    /IMPAR6/ IFREQ,INTYP 
C    /VVNAMS/ VVHDGV  (IVTYPE.I) 
C 
C  Output on return: 
C 
C    function: IMCP81 
C 
C  Local variable / array dictionary: 
C 
C   Name   Type              Description 
C  ------  ----  ------------------------------------------------ 
C  AGE1ST   I    age of the vehicle at first inspection 
C  BY       I    benefit year for technology 1 or 2 vehicle 
C  IBY      I    benefit year for technology 4 plus vehicle 
C  IMPGM    I    The I/M program number that applies to the current 
C                caslendar year, model year, and vehicle type. 
C  IMSTART  I    starting year for continuous I/M ocverage. 
C  IREM     I    remainder = stringency - greatest multiple of 
C                10 < stringency 
C  ISTRN    I    stringency index into the technology 1 or 2 
C                credits array 
C  ITECH    I    pointer to technology and vehicle type 
C  LICYIM   I    I/M program start year that applies 
C  NFREQ    I     Program frequency for current I/M program 
C  PCRED    R    I/M credit as determined from credit deck 
C  REM      R    IREM converted to REAL value 
C  WAIV     R    waiver rate used on this call, depends on MY 
C 
C 
C 
C  Notes: 
C 
C     In earlier versions of the code, I/M credits for test-and-repair 
C     programs, which were disounted by 50% relative to test-only programs, 
C     were not subject to further discounts for waivers unless the waiver 
C     rate was greater than 50%. Now, for waiver rates less than 50% the 
C     waiver discount is phased in as program effectiveness for test-and- 
C     repair programs varies between 50% (the old default) and 100%. 
C 
C     The TECHIV+ code was removed for the Mobile6 version of this code. 
C 
C 
      IMPLICIT NONE 
C 
      INCLUDE 'IM12HC.I' 
      INCLUDE 'IM12CO.I' 
      INCLUDE 'IMPAR1.I' 
      INCLUDE 'IMPAR2.I' 
      INCLUDE 'IMPAR5.I' 
      INCLUDE 'IMPAR6.I' 
      INCLUDE 'IVTYPE.I' 
C 
      INTEGER     FINDIM 
      INTEGER     IMPSYR 
      INTEGER     IMPTR 
C 
      INTEGER     ICY 
      INTEGER     IP  
      INTEGER     IV  
      INTEGER     MY  
C 
      INTEGER     AGE1ST 
      INTEGER     BY 
      INTEGER     IMPGM 
      INTEGER     IMSTART 
      INTEGER     INERR 
      INTEGER     IREM 
      INTEGER     ISTRN 
      INTEGER     ITECH 
      INTEGER     NFREQ 
      REAL        PCRED 
      REAL        REM 
      INTEGER, SAVE :: IEXHAUST = 1 
C 
C . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 
C 
      IMCP81 = 0.0 
      PCRED = 0.0 
C 
C     If there was no I/M program operating in this calendar year for the 
C     specified model year and vehicle type, there is no I/M credit. 
C 
      IMPGM = FINDIM(ICY,MY,IV,IEXHAUST) 
      IF (IMPGM.EQ.0) GOTO 70 
C 
C     No I/M Program Credit is applied when the cars are 
C     "zero" years old 
C 
      IF (ICY.EQ.MY) GOTO 70 
C 
      IF (VVHDGV(IV).EQ.1) GOTO 50 
C 
C     Selecting I/M reduction for LDGV or LDGT.  Several parameters 
C     must be set. Determine technology by model year and vehicle type. 
C 
      ITECH=IMPTR(MY,IV) 
C 
C     Branch on technology type: TECH 1 & 2 form 1 group, TECH 4+ 
C     another. (This subroutine should not be called for TECHIV+ 
C     vehicles. We now treat this as a fatal error. INERR is a 
C     dummy argument to QUITER in this case.) 
C 
      IF (ITECH.GE.4) THEN 
        INERR = 0 
        CALL QUITER(FLOAT(MY),IV,171,INERR) 
        PCRED = 0.0 
        GOTO 70 
      ENDIF 
C 
C     Find the benefit year: 
C     Limit value according to technology type. 
C 
      IMSTART = IMPSYR(ICY,MY,IV) 
      IF (ICY.EQ.IMSTART) THEN 
        PCRED = 0.0 
        GOTO 99 
      ELSE 
        BY = ICY-IMSTART 
      ENDIF 
      IF (MY.GT.IMSTART) BY=ICY-MY 
      IF (BY.GT.19) BY=19 
C 
C     Find the age of the vehicle at first inspection. 
C     Limit value according to technology type. 
C 
      AGE1ST=1 
      IF (MY.LT.IMSTART) AGE1ST=IMSTART-MY+1 
      IF (AGE1ST.GT.19) AGE1ST=19 
      IF (AGE1ST+BY.GT.19) BY=20-AGE1ST 
C 
C     For now, continue MOBILE3 policy of using technology 2 credits for 
C     ITECH=3. 
C 
      IF (ITECH.EQ.3) ITECH=2 
C 
C     There are no NOx credits for Tech 1 & 2. 
C 
      IF (IP.EQ.3) GOTO 70 
C 
C     Select correct I/M credits for TECH 1 & 2.  The same credits 
C     array is used for LDGV and LDGT, but the ITECH mygs differ, so 
C     that the same MY may yield a different credit for LDGV than 
C     for LDGT. 
C 
C     Interpolate between 10, 20, 30, 40 & 50% stringency. 
C 
      IREM=ISTRIN(IMPGM)-(ISTRIN(IMPGM)/10)*10 
      REM=IREM*.1 
      ISTRN=(ISTRIN(IMPGM)-IREM)/10 
C 
      NFREQ=IFREQ(IMPGM) 
C 
      IF (NFREQ.EQ.1) THEN 
C 
        IF (IP.EQ.1) THEN 
          PCRED=CR12HC(BY,AGE1ST,ISTRN,ITECH) 
          IF (ISTRN.LT.5.AND.IREM.GT.0)  PCRED= 
     *      (CR12HC(BY,AGE1ST,ISTRN+1,ITECH)-PCRED)*REM+PCRED 
        ELSE IF (IP.EQ.2) THEN 
          PCRED=CR12CO(BY,AGE1ST,ISTRN,ITECH) 
          IF (ISTRN.LT.5.AND.IREM.GT.0)  PCRED= 
     *      (CR12CO(BY,AGE1ST,ISTRN+1,ITECH)-PCRED)*REM+PCRED 
        ENDIF 
C 
      ELSE IF (NFREQ.GE.2) THEN 
C 
        IF (IP.EQ.1) THEN 
           PCRED=CR12HC(20-BY,21-AGE1ST,ISTRN,ITECH) 
          IF (ISTRN.LT.5.AND.IREM.GT.0)  PCRED= 
     *      (CR12HC(20-BY,21-AGE1ST,ISTRN+1,ITECH)-PCRED)*REM+PCRED 
C 
        ELSE IF(IP.EQ.2) THEN 
          PCRED=CR12CO(20-BY,21-AGE1ST,ISTRN,ITECH) 
          IF (ISTRN.LT.5.AND.IREM.GT.0)  PCRED= 
     *      (CR12CO(20-BY,21-AGE1ST,ISTRN+1,ITECH)-PCRED)*REM+PCRED 
C 
        ENDIF 
C 
      ENDIF 
C 
      GOTO 70 
C 
C     Assign HDGV I/M credit 
C     and there are no TIER1 HDGV I/M credits 
C 
   50 CONTINUE 
      NFREQ=IFREQ(IMPGM) 
      IF (NFREQ.GT.2) NFREQ = 2 
      PCRED=CRHDGV(IP,NFREQ) 
C 
   70 CONTINUE 
      IMCP81 = PCRED 
C 
   99 CONTINUE 
      RETURN 
      END 
