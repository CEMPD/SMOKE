      REAL FUNCTION EF_HI_FRCN(MY,IP,IV,VMTAGE,IFDS,ISTD,INERR)
C
C EF_HI_FRCN returns the proportion of high emitters amongst the light duty
C gasoline vehicles and trucks. The effect of tampering is not included.
C
C  Called by HRLOOP and IMBENEFIT
C
C  Calls EF81_93, EF81_93_AVE, QUITER
C
C  Changes: (Last change first)
C
C  15 Apr 02  AIR Task Bug Fix 388: modification of constants via QUITER.
C  17 Oct 01 @EPA-BAG bug331 COMMENTS - Biennial/ annual I&M program discrepancy
C  22 Jun 00 EPA-BAG  removed JDX from EF81_93_AVE parameter list
C  09 Jun 00 EPA-BAG  removed JDX from EF81_93 parameter list
C  08 Jun 00 @EPA-djb Added ISTD to the call to choose method used.
C  24 Apr 00 AIR Task02: Fixed the parameter list to include INERR.
C   7 Jan 00 @DynTel-ddj 1-024, Modified DEN and NUM equations for HC
C            and NOx.  Added condition on EF_HI_FRCN to cap at 1.0 if
C            NUM is greater than DEM.
C  12 Oct 99 @DynTel-BAG Fixed discrepancy in the  mileage dependence for
C            determining high emitter vehicle fractions. There was an
C            inconsistency between what was in the code for NOx and what is
C            in M6.EXH.007. Eqs. 8b-8d (page 9) have mileage ranges whereas
C            the NOx code does not.
C  09 SEP 99 @DynTel-BAG Renamed to EF_HI_FRCN.FOR
C  26 AUG 99 @DynTel-HXQ Deleted the non-used flag IGM
C  16 Jun 99 @DynTel-BG 2-000 Added RUN, NORMAL, and HIGH parameters
C  27 May 99 @DynTel-ddj 2-700 Updated EF_HI_FRCN calculation for NOx.  Made
C            it consistent with M6.EXH.007, March 1999, Eqn 10.
C  02 Mar 99 @Dyntel-RJD 2-694 Updated EF_HI_FRCN to calculate high emitter
C            fraction for 1994+ LDGV/Ts.
C  05 Nov 1998 @DynTel-ZK 2-000 Explicit Typing
C  22 Sep 98 @Dyntel-RJD 2-663 Updated EF81_93 function parameter list.
C  14 Sep 98 @Dyntel-RJD 0-632 Added call to QUITER if denominator in
C            highs formula is in the neighborhood of zero. This is
C            primarily to trap errors in user input data.
C  05 Aug 98 @DynTel-RJD 2-663 Changed integer fuel delivery system
C            parameter to IFDS.
C
C  Input on call or from calls:
C
C    argument list: MY,IP,IV,VMTAGE,IFDS
C    common blocks:
C
C  Output on return:
C
C    function: EF_HI_FRCN
C    common blocks:
C
C  Local variable / array dictionary:
C
C   Name   Type              Description
C  ------  ----  -----------------------------------------------------
C  VMTAGE   R    vehicle mileage in units of 1/10,000 miles
C
C  Notes:
C
C  13 Jul 98 @DynTel-RJD 2-660 Created EF_HI_FRCN to calculate EF_MCDSLPRE81.
C
      IMPLICIT NONE
C
C     declare external functions
C
      REAL    EF81_93
      REAL    EF81_93_AVE
C
C     declare parameter list
C
      INTEGER, INTENT(IN) :: MY
      INTEGER, INTENT(IN) :: IP
      INTEGER, INTENT(IN) :: IV
      INTEGER, INTENT(IN) :: IFDS
      INTEGER, INTENT(IN) :: ISTD
      INTEGER, INTENT(INOUT) :: INERR
      REAL,    INTENT(IN) :: VMTAGE
C
C     declare local variables
C
      INTEGER MY1
      REAL    DEN
      REAL    NUM
      INTEGER, PARAMETER :: RUN=1
      INTEGER, PARAMETER :: NORMAL=1,HIGH=2
      INTEGER, PARAMETER :: HC=1, CO=2
C
      IF ( (MY.LE.1993) .OR.
     *     (MY.GT.1993 .AND. ISTD.EQ.1 ) ) THEN
        MY1=MY
        IF(MY.GT.1993) MY1=1993
        DEN = EF81_93(MY1,IP,IV,VMTAGE,IFDS,RUN,HIGH,INERR)
     &        - EF81_93(MY1,IP,IV,VMTAGE,IFDS,RUN,NORMAL,INERR)
        NUM = EF81_93_AVE(MY1,IP,IV,VMTAGE,IFDS,INERR)
     &        - EF81_93(MY1,IP,IV,VMTAGE,IFDS,RUN,NORMAL,INERR)
C
      ELSE
        SELECT CASE (IP)
          CASE (HC)
C
C Assume the fraction of high emitters is the same as for
C tier 0 PFI vehicles
C
            DEN = EF81_93(1993,HC,1,VMTAGE,3,RUN,HIGH,INERR)
     &            - EF81_93(1993,HC,1,VMTAGE,3,RUN,NORMAL,INERR)
            NUM = EF81_93_AVE(1993,HC,1,VMTAGE,3,INERR)
     &            - EF81_93(1993,HC,1,VMTAGE,3,RUN,NORMAL,INERR)
          CASE (CO)
            DEN = EF81_93(1993,CO,IV,VMTAGE,3,RUN,HIGH,INERR)
     &            - EF81_93(1993,CO,IV,VMTAGE,3,RUN,NORMAL,INERR)
            NUM = EF81_93_AVE(1993,CO,IV,VMTAGE,3,INERR)
     &            - EF81_93(1993,CO,IV,VMTAGE,3,RUN,NORMAL,INERR)
          CASE DEFAULT   ! NOx
            DEN =  1.141 - 0.0294*VMTAGE  ! see eq. 5-10 p.8-12
            NUM = -0.036 + 0.02142*VMTAGE ! M6.EXH.007, 12/99
          END SELECT
      END IF
C
      IF (DEN.LT.0.000001) THEN
        EF_HI_FRCN=0.0
        CALL QUITER(0.,0,152,INERR)
      ELSE IF(NUM .LE. 0.0) THEN
           EF_HI_FRCN = 0.0
        ELSE IF(NUM .GE. DEN) THEN
             EF_HI_FRCN = 1.0
          ELSE
             EF_HI_FRCN = NUM / DEN
      END IF
C
      RETURN
      END
