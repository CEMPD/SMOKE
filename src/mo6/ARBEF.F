      REAL FUNCTION ARBEF(MY,JDX,IP,IV,VMTAGE,IFDS,ISR,IMPGM,ISTD,INERR)
C
C  ARBEF returns the after repair emission factors for 1981-1996 model
C  year light duty gasoline vehicles and trucks. See report M6.IM.001
C  pages 16-21 for documentation of running emissions.
C
C  Called by IMPROG
C
C  Calls EF81_93, EF_POST93
C
C  Changes:
C
C   15 Apr 02  AIR Task Bug Fix 388: modification of constants via QUITER.
C   17 Oct 01 @EPA-BAG bug331 Biennial/ annual I&M program discrepancy
C   25 May 01 @EPA-elg Loaded idle I&M
C   30 Aug 00 @EPA-elg bug244 Technician credit effects are no longer applied
C                       just to the IM240 tests, but to ASM and OBD as well.
C   14 Aug 00 @ EPA-djb Bug 226, Added local variables AIM and LFDS.
C   29 Jul00 @EPA-jk  after repair start emissions conditional for
C            model years 1981 through 1993; included 1981, i.e.,
C            used ".GE.1981" rather than ".GT.1981"
C   18Jul 00 @EPA-ddj combined conditionals for armax and armin
C   07Jul 00 AIR Task 05: Removed unused JDX from function EF_POST93
C            parameter list.
C  09 Jun 00 EPA-BAG  removed JDX from EF81_93 parameter list
C  24 Apr 00 AIR Task 02: Fixed call to EF91-93 and EF_POST93 to remove
C             spurious 2 in parameter lists.
C   8 Sep 99 @Dyntel-bg  2-660 Renamed HNBEF and HN94B to EF81_93 and
C            EF_POST93, respectively
C   2 Aug 99 @Dyntel-bg  2-663 Replaced IF statements over pollutants
C            and start/running emissions with CASE statements; determination
C            of 1981 after repair emissions was being skipped; replaced
C            constants with descriptive parameters
C  17 Mar 99 @Dyntel-RJD 2-694 Extended the scope of ARBEF so that it
C            also returns the after repairs emissions for 1994+ model
C            year LDGV/Ts.
C
C  Input on call or from calls:
C
C    argument list: MY,JDX,VAGE_LSTIM,IP,IV,VMTAGE,IFDS,ISR,IMPGM,ISTD
C    common blocks:
C    /IMPAR2/ ITEST
C    /IMPAR4/ CUTHC, CUTCO, CUTNO, IVCUTIDX
C    /LDGBER/ CRT050K
C    /LDGIMC/ ASEBZML
C    /TTCINF/ ITTCFG
C
C  Output on return:
C
C    function: ARBEF
C    common blocks:
C
C  Local variable / array dictionary:
C
C   Name    Type              Description
C  ------   ----  -----------------------------------------------------
C   AIM       R   Local REAL value of JDX.
C   ARMAX     R   Maximum after repairs emission factor
C   ARMIN     R   Minimum after repairs emission factor
C   CERTRATIO R   Adjusts Tier 0 emissions to account for higher 1994+
C                 emissions standards
C   COPC      R   The after repairs emission factor is equal to COPC
C                 times the normal emission factor for IP=2.
C   HCPC      R   The after repairs emission factor is equal to HCPC
C                 times the normal emission factor for IP=1.
C   LFDS      I   Local value of IFDS.
C   NOPC      R   The after repairs emission factor is equal to NOPC
C                 times the normal emission factor for IP=3.
C
C  Notes:
C
C  04 Aug 98 @DynTel-RJD 2-663 Created ARBEF to calculate the after
C            repair emissions for light duty gasoline vehicles.
C
      IMPLICIT NONE
      INCLUDE 'IMPAR2.I'
      INCLUDE 'IMPAR4.I'
      INCLUDE 'LDGBER.I'
      INCLUDE 'LDGIMC.I'
      INCLUDE 'TTCINF.I'
C
C     Declare external functions.
C
      REAL    EF81_93,EF_POST93
C
C     Declare parameter list.
C
      INTEGER, INTENT(IN) :: MY
      INTEGER, INTENT(IN) :: JDX
      INTEGER, INTENT(IN) :: IP
      INTEGER, INTENT(IN) :: IV
      INTEGER, INTENT(IN) :: IFDS
      INTEGER, INTENT(IN) :: ISR
      INTEGER, INTENT(IN) :: IMPGM
      INTEGER, INTENT(IN) :: ISTD
      INTEGER, INTENT(INOUT) :: INERR
      REAL,    INTENT(IN) :: VMTAGE
C
C     Declare local variables.
C
      INTEGER IG
      INTEGER LFDS
      REAL    AIM
      REAL    HCPC,COPC,NOPC
      REAL    ARMIN,ARMAX,CERTRATIO
C
      INTEGER, PARAMETER :: RUN=1,START=2
      INTEGER, PARAMETER :: NORMAL=1,HIGH=2
      INTEGER, PARAMETER :: HC=1,CO=2,NOX=3
C
C     Initialize local variables
C
      AIM = REAL(JDX)
      LFDS = IFDS
C
C  After repairs running emissions. These are calculated as a
C  percentage of the normal running emissions.
C
      START_OR_RUN : SELECT CASE(ISR)
C
      CASE(RUN) START_OR_RUN
C
        IF (AIM.GE.15.0) AIM=15.0
C
        POLLUTANT : SELECT CASE (IP)
C
        CASE (HC) POLLUTANT
          HCPC = 2.2400- 0.07595*AIM
          HCPC = (0.4990 *CUTHC(JDX,IVCUTIDX(IV),IMPGM)
     *         - 1.011E-4*CUTCO(JDX,IVCUTIDX(IV),IMPGM)+0.398)*HCPC
C
C         The idle I&M tests, ITEST(1) and ITEST(2), have a 1.5 times
C         larger after repair emission than the I&M 240 test. (See page 67
C         of M6.IM.001.)
C
          IF (ITEST(IMPGM).LE.2) HCPC=1.5*HCPC
C
C         If the I&M 240 test has an additional factor for technician training
C         programs. (See page 22-23 of M6.IM.001.)
C
          IF (ITTCFG(IMPGM).EQ.1) HCPC=1.78*HCPC
C
          IF (MY.LE.1993) THEN
            ARBEF=HCPC*EF81_93(MY,IP,IV,VMTAGE,LFDS,ISR,NORMAL,INERR)
          ELSE
            ARBEF=HCPC*EF_POST93(MY,IP,IV,VMTAGE,ISR,NORMAL,ISTD,INERR)
          END IF
C
        CASE(CO) POLLUTANT
C
          COPC = 2.1582- 0.07825*AIM
          COPC = (0.0249*CUTHC(JDX,IVCUTIDX(IV),IMPGM)
     *         +  0.0168*CUTCO(JDX,IVCUTIDX(IV),IMPGM)+0.620)*COPC
C
C         idle I&M tests :
C
          IF (ITEST(IMPGM).LE.2) COPC=1.50*COPC
C
C         I&M 240 technician training :
C
          IF (ITTCFG(IMPGM).EQ.1) COPC=2.74*COPC
C
          IF (MY.LE.1993) THEN
            ARBEF=COPC*EF81_93(MY,IP,IV,VMTAGE,LFDS,ISR,NORMAL,INERR)
          ELSE
            ARBEF=COPC*EF_POST93(MY,IP,IV,VMTAGE,ISR,NORMAL,ISTD,INERR)
          END IF
C
        CASE(NOX) POLLUTANT
C
          NOPC = 1.6410- 0.04348*AIM
          NOPC = (0.2538*CUTNO(JDX,IVCUTIDX(IV),IMPGM)+0.2613)*NOPC
C
C         I&M 240 technician training :
C
          IF (ITTCFG(IMPGM).EQ.1) NOPC=1.39*NOPC
C
          IF (MY.LE.1993) THEN
            ARBEF=NOPC*EF81_93(MY,IP,IV,VMTAGE,LFDS,ISR,NORMAL,INERR)
          ELSE
            ARBEF=NOPC*EF_POST93(MY,IP,IV,VMTAGE,ISR,NORMAL,ISTD,INERR)
          END IF
C
        END SELECT POLLUTANT
C
C After repairs starts emissions. The data used to calculate these
C emissions is stored in BD43.
C
      CASE(START) START_OR_RUN
C
        IG=1
        IF (MY.GE.1981.AND.MY.LE.1993) THEN
          IG=MY-1980
          ARBEF= ASEBZML(IG,IP,IV,LFDS)
        ELSE IF (MY.GT.1993) THEN
          IG=13
          LFDS=3
C
C Calculate the certification ratio
C
          IF (VMTAGE.LT.5.0) THEN
            CERTRATIO = CRT050K(ISTD,IP,IV)/CRT050K(1,IP,1)
          ELSE
            CERTRATIO = CRT050K(ISTD,IP,IV)/CRT050K(1,IP,1)
          END IF
          ARBEF = CERTRATIO*ASEBZML(IG,IP,1,LFDS)
        END IF
C
      END SELECT START_OR_RUN
C
C After repair emissions may become higher than the high basic emission
C factor. In this case the after repair emissions will be set equal to
C the high emission factor.
C
      IF (MY.LE.1993) THEN
        ARMAX=EF81_93(MY,IP,IV,VMTAGE,LFDS,ISR,HIGH,INERR)
        ARMIN=EF81_93(MY,IP,IV,VMTAGE,LFDS,ISR,NORMAL,INERR)
      ELSE
        ARMAX=EF_POST93(MY,IP,IV,VMTAGE,ISR,HIGH,ISTD,INERR)
        ARMIN=EF_POST93(MY,IP,IV,VMTAGE,ISR,NORMAL,ISTD,INERR)
      END IF
C
      IF (ARBEF.GE.ARMAX) ARBEF=ARMAX
C
C After repair emissions may become lower than the normal basic emission
C factor. In this case the after repair emissions will be set equal to
C the normal emission factor.
C
      IF (ARBEF.LE.ARMIN) ARBEF=ARMIN
C
      RETURN
      END
