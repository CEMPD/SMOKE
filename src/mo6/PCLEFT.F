      REAL FUNCTION PCLEFT(MY,ICY,IP,IV,INERR)
C
C  PCLEFT determines the basic emission factor multiplicative
C  adjustment ( = percentage left of BEF = PCLEFT ) for the
C  effects of an inspection / maintenance (I/M) program.
C
C  Called by BEF and BEF94P.
C
C  Calls IMCP81, IMPSYR
C
C  Changes: (Last change first)
C
C  15 Apr 02  AIR Task Bug Fix 388: modification of constants via QUITER.
C  20 Sep 00 @EPA- BG Added constant IEXHAUST=1 for call to FINDIM
C  15 Sep 00 AIR Task 08: Updated FINDIM calls to select exhaust I/M.
C  11 Feb 99 @DynTel-MLA 2-663  This routine was restructured for Mobile6.
C            Much of the code that used to be in this routine is now in
C            subroutine IMCP81, which is called by this one.
C  17 Dec 98 @DynTel-ZK 2-000 Explicit Typing 
C
C
C  Input on call:
C
C    argument list: MY,ICY,IP,IV
C    common blocks:
C    /IMPAR1/ WAVPRE81,WAV81PLUS,CRIM
C    /IMPAR5/ DISCNT
C    /IMPAR6/ IFREQ,INTYP
C
C  Output on return:
C
C    function: PCLEFT
C
C  Local variable / array dictionary:
C
C   Name   Type              Description
C  ------  ----  ------------------------------------------------
C  IMPGM     I    The I/M program number that applies to the current
C                 calendar year, model year, and vehicle type.
C  IMCREDIT  R    I/M credit returned by IMCP81.
C  IMPGM     I    I/M program number.
C  IMSTART   I    First year of I/ M coverage.
C  WAIVER    R    I/M waiver rate.
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
      INCLUDE 'IMPAR1.I'
      INCLUDE 'IMPAR5.I'
      INCLUDE 'IMPAR6.I'
C
      REAL        ENFORC
      INTEGER     FINDIM
      REAL        IMCP81
      INTEGER     IMPSYR
C
      INTEGER, INTENT(IN) :: ICY
      INTEGER, INTENT(IN) :: IP
      INTEGER, INTENT(IN) :: IV
      INTEGER, INTENT(IN) :: MY
      INTEGER, INTENT(INOUT) :: INERR
C
      REAL        IMCREDIT
      INTEGER     IMPGM
      INTEGER     IMSTART
      REAL        WAIVER
      INTEGER, SAVE :: IEXHAUST=1
C
C . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
C
      PCLEFT = 1.0
      IMCREDIT = 0.0
      WAIVER = 0.0
C
C     If there was no I/M program operating in this calendar year for the
C     specified model year and vehicle type, there is no I/M credit.
C
      IMPGM = FINDIM(ICY,MY,IV,IEXHAUST)
      IF (IMPGM.EQ.0) GOTO 99
C
C     Call IMCP81 to get the I/M credits.
C
      IMCREDIT = IMCP81(MY,ICY,IP,IV,INERR)
      IF(INERR.GT.0) GOTO 99
C
C     Check to see if the credits need to be phased in. No I/M credit
C     is given in the first year of operation of any program. After
C     that, if the inspection frequency is 'n' years, 1/n of the benefits
C     are added each year until the program is operating at full credit.
C
      IMSTART = IMPSYR(ICY,MY,IV)
      IF (ICY.LT.(IMSTART+IFREQ(IMPGM))) THEN
C
        IF (MY.EQ.IMSTART) THEN
          IMCREDIT = 0.0
        ELSE
          IMCREDIT = IMCREDIT *
     &      (ICY-IMSTART)/FLOAT(IFREQ(IMPGM))
        ENDIF
C
      ENDIF
C
C     Assign waiver rates
C     Control WAIVER (waver rates) for INTYP (I/M program inspection
C     type), where 2 = I/M test and repair (Computerized),
C     & 3 = test and repair (Manual), if INTYP is not a 1 = test
C     only I/M and the WAIVER is LE 50%, WAIVER is set to zero.
C
      IF (MY.LE.1980) WAIVER=WAVPRE81(IMPGM)
      IF (MY.GT.1980) WAIVER=WAV81PLUS(IMPGM)
C
      IF (WAIVER.LE.0.50.AND.INTYP(IMPGM).GT.1) THEN
        IF (DISCNT(IP, INTYP(IMPGM)).LT.0.50) THEN
          WAIVER = 0.0
        ELSE
          WAIVER=WAIVER*(2*DISCNT(IP,INTYP(IMPGM))-1)
        ENDIF
      ENDIF
C
   80 PCLEFT = 1.0-(IMCREDIT*(1.0-WAIVER)*ENFORC(CRIM(IMPGM),1)*
     *              DISCNT(IP,INTYP(IMPGM)))
C
   99 CONTINUE
      RETURN
      END
