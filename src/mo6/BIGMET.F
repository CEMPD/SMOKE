      SUBROUTINE BIGMET(IV,MY,IDX,JDX,INERR)
C
C  BIGMET computes the additive methane offset to the composite
C  HC emission factor without a cap for tier 1 and lev.
C
C  Called by BIGCFX.
C
C  Calls EXMPTR, IVPTRT, IEVPTR, EF81_93_AVE
C
C  Changes: (Last change first)
C
C  15 Apr 02 AIR Task Bug Fix 388: modification of constants via QUITER.
C  04 Oct 00 @EPA-djb Bug276 Removed standards ratio from 1994+ MOFFMTH
C            calculation and cleaned up code.
C  07 Jul 00 AIR Task 05: Removed unused MYR from function METH81 parameter list.
C  22 Jun 00 EPA-BAG  removed JDX from EF81_93_AVE and STRT_AVE parameter lists
C  15 Jun 00 @EPA-dj   Initialized MOFFMTH Array to 0.0 rather than 1.0.
C  25-Apr-2000  @EPA-elg Added a multiplicative methane adjustment for 1981 LDGV/T.
C   4 Oct 99 @DynTel-ddj 2-714  Added IFDS index to OFFMTH.
C  28 Sep 99 @DynTel-ddj for EPA-DJB Fixed bag weighting and removed
C            start/run loop.
C   5 Mar 99 @DynTel-DDJ 2-696 Added hourly loop
C  14 Aug 98 @DynTel-zk 2-000  Explicit typing
C  11 Aug 98 @DynTel-MLA 2-000  Changed HDG and HDD names to conform to
C            new naming conventions.
C  03 Aug 98 @DynTel-bag 2-000  Added comments and declared NBAG as an
C            integer variable
C  30 Jul 98 @DynTel-RJD 2-660  Disabled temperature correction for
C            vehicles with separated start and running emissions.
C  23 Jul 98 @DynTel-MLA 2-000  Changed LDGT names to conform to
C            certification naming conventions.
C
C  Input on call:
C
C    Constants: 
C    MAXVEH from include file IVTYPE.I
C
C    argument list: IV,MY,IDX,JDX,ISR
C    common blocks:
C    /CUMCOM/ CUMMIL
C    /NMETH1/ ZDMTH1,ZDMTH2,ZDMTH3
C    /NMETH2/ MTHMYG
C    /REGION/ IREJN
C    /VVNAMS/ VVDSL, VVHEAVY, VVLDGAS, VVHDGV  (IVTYPE.I)
C    /VVTYPS/ VTMC  (IVTYPE.I)
C
C  Output on return:
C
C    common blocks:
C    /OFFSET/ OFFMTH, MOFFMTH
C
C  Local array subscripts
C
C  CRBFNJ(2) - CRBFNJ(IFDS)
C
C  Local variable / array dictionary:
C
C   Name   Type              Description
C  ------  ----  ------------------------------------------------
C  CRBFNJ   R    fraction of carb and fuel injected vehicles
C  IB       I    FTP bags (1-3)
C  IFDS     I    fuel deleivery system index,
C                1=carb, 2=injection
C  IVSCAN   I
C  IVTEMP   I    temporary IV index into ZDMTHi arrays
C  DENOM    R    emission rate denominator used to calculate MOFFMTH ratio
C  OFF      R    by-bag methane offset
C  VMT      R    local mileage version of VMTAGE
C
C  Notes:
C
C  BIGMET was modified for MOBILE4.1 to include methane offset
C  ZML and DR by fuel delivery system and region.
C  BIGMET was modified for MOBILE5v13 to include a cap for
C  methane offset for lev and tier 1 with special I/M.
C
C  20-July-1994 @CSC-pme request 2-446.
C               Cap for  Tier 1 and Lev was removed. 
C  24-August-1994 @CSC-pme Request 2-446.
C               Include file BASE12.I was removed.
C               The third dimension of OFFMTH was removed.
C  11 Jun 97 DynTel@MLA 2-622  Added IVTYPE.I for Mobile6 and modified the
C            code to use the symbolic vehicle type names and vehicle type
C            vectors that were introduced in Mobile6.
C  16 Apr 98 @DynTel-RJD 2-654 Added ISR loop variable to input parameter 
C            list. Expanded OFFMTH array to include separate starts and
C            running emissions. No adjustment for methane is made to starts 
C            emissions (ISR=2).
C  17 Sep 99 @EPA-DJB Took BIGMTH out of the start/run loop and changed
C            the routine to calculate both start and running methane
C            emissions in one call.  BIGMTH now calls METH81 and METH94
C            to calculate start and running methane emissions for 1981
C            and newer light duty cars and trucks.
C
      IMPLICIT NONE
      INCLUDE 'CUMCOM.I'
      INCLUDE 'IVTYPE.I'
      INCLUDE 'NMETH1.I'
      INCLUDE 'NMETH2.I'
      INCLUDE 'OFFSET.I'
      INCLUDE 'REGION.I'
C
C     Declare parameter list
C
      INTEGER, INTENT(IN) :: JDX
      INTEGER, INTENT(IN) :: IV
      INTEGER, INTENT(IN) :: MY
      INTEGER, INTENT(IN) :: IDX
      INTEGER, INTENT(INOUT) :: INERR
C
C     Declare local variables
C
      INTEGER   IVSCAN, IB, IFDS, ISR, IVTEMP
      INTEGER   ISTD
      REAL  VMT, METH81, METH93
      REAL, DIMENSION(3,2) :: OFF
C
C     Declare external functions
C
      REAL  EF81_93_AVE, DENOM, STRT_AVE
C
C  Initialize.
C
      DO ISR = 1, 2
        DO ISTD = 1, 10
          DO IFDS = 1, 3
            OFFMTH(IDX,IV,ISR,ISTD,IFDS)  = 0.0
            MOFFMTH(IDX,IV,ISR,ISTD,IFDS) = 0.0
          END DO
        END DO
      END DO
C
C  Calculate vehicle accumulated mileage
C
      VMT=CUMMIL(JDX,IV)/10000.
C
C  Set model year group pointer
C
      CALL EXMPTR(MY,IV,INERR)
      IF(INERR.GT.0) GOTO 99
C
      IF (MY.LT.1981.OR.VVLDGAS(IV).NE.1) THEN
C
C  Determine LDGV, LDGT12, and LDGT34 offsets.
C
        IF(VVLDGAS(IV).EQ.1) THEN
          DO IB = 1,3
            DO IFDS = 1, 2
               OFF(IB,IFDS) = ZDMTH1(1,IB,IFDS,MTHMYG,IV,IREJN)
     *           + ZDMTH1(2,IB,IFDS,MTHMYG,IV,IREJN)*VMT
            END DO
          END DO
C
          DO IFDS = 1, 2
C
C  Determine runnning emission methane emissions (ISR=1) in g/mi
C  using the combined Bag 2 and Bag 3 methane emission rates
C  using the LA4 VMT bag weighting factors.
C
            OFFMTH(IDX,IV,1,1,IFDS) = 0.521*OFF(2,IFDS) +
     &                                0.479*OFF(3,IFDS)
C
C  Determine the engine start methane emissions (ISR=2)
C  using the difference between Bag 1 and Bag 3
C  and convert to grams by multiplication by the number
C  of miles in the bag.
C
            OFFMTH(IDX,IV,2,1,IFDS) = (OFF(1,IFDS) - OFF(3,IFDS))*3.49
          END DO
          OFFMTH(IDX,IV,1,1,3) = OFFMTH(IDX,IV,1,1,2)
          OFFMTH(IDX,IV,2,1,3) = OFFMTH(IDX,IV,2,1,2)
C
C  Determine HDGV and MC offsets.
C
        ELSE IF(VVHDGV(IV).EQ.1.OR.IV.EQ.VTMC) THEN
C
C       Convert IV to an index for array ZDMTH2, which is numbered from 1 
C       (for HDGV2B) to 9 (for MC) to 10 (for GASBUS).
C       
          IF(IV.EQ.24) THEN
            IVTEMP = 9
          ELSE IF(IV.EQ.25) THEN
            IVTEMP = 10
          ELSE
            IVTEMP = IV-5
          END IF
C
          OFFMTH(IDX,IV,1,1,1) = ZDMTH2(1,MTHMYG,IVTEMP,IREJN) +
     *                        ZDMTH2(2,MTHMYG,IVTEMP,IREJN)*VMT
          DO IFDS = 2, 3
            OFFMTH(IDX,IV,1,1,IFDS) = OFFMTH(IDX,IV,1,1,1)
          END DO
C
C  Determine LDDV, LDDT, HDDV, and Diesel BUS offsets.
C
        ELSE IF (VVDSL(IV).EQ.1) THEN
C
C      Convert IV to an index for array ZDMTH3, which is numbered from 1
C      (for LDDV) to 13 (for SCHOOL BUS).
C       
          IVTEMP = 0
          DO IVSCAN = 1, MAXVEH
             IF (VVDSL(IVSCAN).EQ.1) IVTEMP = IVTEMP + 1
             IF (IVSCAN.EQ.IV) GO TO 17
          END DO
  17      OFFMTH(IDX,IV,1,1,1) = ZDMTH3(MTHMYG,IVTEMP,IREJN)
          DO IFDS = 2, 3
            OFFMTH(IDX,IV,1,1,IFDS) = OFFMTH(IDX,IV,1,1,1)
          END DO
C
C  Caution: Methane values for 1987 and newer HDGV and HDDV are
C           recalculated in routine METHANE.
C           Methane values for all model year MC are recalculated
C           in routine METHANE and converted to start/running.
C           OFFMTH for start for MC here are still in g/mi units.
C
        END IF
C
C  If the vehicle is LDV or LDT and the model year is 1981 or newer,
C  determine OFFMTH as a ratio specifically from coefficients.
C
      ELSE
         DO ISR = 1,2
            IF(MY.LE.1993) THEN
              DO IFDS = 1,3
                IF(ISR.EQ.1) THEN
                  DENOM = EF81_93_AVE(MY,1,IV,VMT,IFDS,INERR)
                ELSEIF(ISR.EQ.2) THEN
                  DENOM = STRT_AVE(MY,1,IV,VMT,IFDS,INERR)
                ENDIF
                  MOFFMTH(IDX,IV,ISR,1,IFDS) =
     &                 METH81(IV,MY,ISR,VMT,IFDS) / DENOM    !Ratio for THC case.
              END DO
            ELSEIF(MY.GE.1994) THEN
              DO IFDS = 1,3
                METH93 = METH81(IV,1993,ISR,VMT,IFDS)
                IF(ISR.EQ.1) THEN
                  DENOM  = EF81_93_AVE(1993,1,IV,VMT,IFDS,INERR)
                ELSEIF(ISR.EQ.2) THEN
                  DENOM  = STRT_AVE(1993,1,IV,VMT,IFDS,INERR)
                ENDIF
              DO ISTD = 1,10
                    MOFFMTH(IDX,IV,ISR,ISTD,IFDS) =
     &                  METH93 / DENOM                       !Ratio for THC case.
                END DO       !ISTD Loop
              END DO         !IFDS Loop
            END IF           !MY Conditional
        END DO               !ISR Loop
      END IF                 !VVLDGAS Conditional
C
  99  RETURN
      END
