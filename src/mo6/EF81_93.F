      REAL FUNCTION EF81_93(MY,IP,IV,VMTAGE,IFDS,ISR,INH,INERR)
C
C  EF81_93 returns the normal and high emission factors for 1981-1993
C  model year LDGV/Ts.
C
C  Called by ARBEF, EF_DSPTCH, EF_MCDSLPRE81, EF_HI_FRCN, EF_POST93
C
C  Calls
C        CHK_USER(in USERDAT),IERPTR
C
C  Changes: (Last change first)
C
C  15 Apr 02  AIR Task Bug Fix 388: modification of constants via QUITER.
C  12 Dec 01 EPA-BAG  Added variable initializations for EF81_93, ZPOINT, and SLOPE
C  09 Jun 00 EPA-BAG  Declared parameter list variables with INTENT(IN)
C                     qualifier; removed JDX from parameter list
C  28 Apr 00 EPA-elg  Corrected Kmiles to be 10*VMTAGE.  Erroneous in previous code.
C  25 Apr 00 AIR Task 02: Change code to upper case to avoid problems
C            with some FORTRAN compilers.
C  06 Mar 00 @EPA-eg           hrloop revisions; kmiles=vmyage rather than
C                              kmiles=10*vmtage
C  08 Sep 99 @DynTel-bag 2-660 Renamed from HNBEF.FOR to EF81_93.FOR.
C  26 AUG 99 @DynTel-HXQ Deleted the non-used flag IGM
C  05 Aug 99 @Dyntel NH. Replaced integer constant(vehicle index) in HREBZML
C            with IV,i.e., HREBZML(IGER,NOX,1,IFDS) with HREBZML(IGER,NOX,IV,IFDS)
C            Added comments.
C  07 Jul 99 @Dyntel-NH  Added USE USERDAT reference and code to invoke
C            use of applicable user record.
C  27 May 99 @Dyntel-DDJ 2-700 Replaced IF statements with CASE statements
C            and added constants to hold values of standard MOBILE6
C            variables
C  18 May 99 @Dyntel-HXQ 2-650 delete the code of changing unit of basic
C                        emission from gram/start to gram/mile
C  10 Mar 99 @DynTel-RJD 2-694 Adapted EF81_93 for use in the default
C            calculation for the 1994+ LDGV/T high emitter fraction.
C  22 Sep 98 @DynTel-RJD 2-663 Added new IGM parameter to parameter
C            list. EF81_93 will return a value in g/mile if IGM=1 else
C            it will return a value in g/start.
C  18 Nov 98 @DynTel-ZK  2-000 Explicit typing
C  25 Aug 98 @DynTel-RJD 2-660 Updated IERPTR function parameter list
C  05 Aug 98 @DynTel-RJD 2-663 Changed integer fuel delivery system
C            parameter to IFDS.
C
C
C  Input on call or from calls:
C
C    argument list: MY,IP,IV,VMTAGE,IFDS,ISR,INH
C    common blocks:
C
C    /CEVBMY/ BMYMPD,BMYTPD
C    /LDGBER/ NSEBZML,NSEBDR,HSEBZML,NREBZML,NREBDR,HREBZML
C
C  Output on return:
C
C    function: EF81_93
C    common blocks:
C
C  Local variable / array dictionary:
C
C   Name   Type              Description
C  ------  ----  -----------------------------------------------------
C  SLOPE1   R    slope of ef curve
C  ZPOINT   R    zero mileage point = intercept of ef curve
C  KMILES   R    vehicle miles traveled in kilomiles.
C
C  Notes:
C
C  14 Jul 98 @DynTel-RJD 2-660 Created EF81_93 to calculate the normal
C            (INH=1) and high emission factors for light duty gasoline
C            vehicles.
C  24 Jul 98 @DynTel-BAG 2-000 Deleted unused variable LMY
C  03 Aug 98 @DynTel-RJD 2-663 Moved calculation of after repairs
C            emissions to separate module called ARBEF.
C
      USE USERDAT
      IMPLICIT NONE
      INCLUDE 'CEVBMY.I'
      INCLUDE 'LDGBER.I'
C
      INTEGER, INTENT(IN) :: IFDS
      INTEGER, INTENT(IN) :: MY
      INTEGER, INTENT(IN) :: IP
      INTEGER, INTENT(IN) :: IV
      INTEGER, INTENT(IN) :: ISR
      INTEGER, INTENT(IN) :: INH
      INTEGER, INTENT(INOUT) :: INERR
      REAL,    INTENT(IN) :: VMTAGE
C
      INTEGER IGER,IERPTR,LMY,KEYER
      REAL KMILES, ZPOINT,SLOPE1,USER_ER
      LOGICAL USER_REC
C
      INTEGER, PARAMETER :: NOX=3
      INTEGER, PARAMETER :: RUN=1
      INTEGER, PARAMETER :: NORMAL=1
C
      EF81_93 = 0.0
      ZPOINT  = 0.0
      SLOPE1  = 0.0
      KMILES = VMTAGE*10.0
C
      IF (MY.LE.1993) THEN
        LMY  = MY
        IGER = IERPTR(ISR,LMY,IP,IV,KEYER,INERR)
        IF(INERR.GT.0) GOTO 99
      ELSE
C
C 1981-1993 model year user input should not override default basic
C emission factor coefficients in 1994+ calculations.
C
        IGER  = 13
        KEYER =  0
      END IF
C
C Test for BD40 block data input or user input
      CALL CHK_USER(MY,IP,IV,IFDS,ISR,INH,KMILES,USER_ER,USER_REC)
      IF (USER_REC) THEN
          EF81_93=USER_ER !Use user input data
      ELSE
C
C
C  BD40 block data input
C
      SELECT CASE (IP)
         CASE(NOX)
           IF(INH.EQ.NORMAL.AND.ISR.EQ.RUN) THEN
              ZPOINT = NREBZML(IGER,NOX,IV,IFDS)
              SLOPE1 = NREBDR(IGER,NOX,IV,IFDS)
              EF81_93  = ZPOINT + SLOPE1*KMILES
           ELSE
C
C  High running emissions
C
              IF(ISR.EQ.RUN) THEN
                 EF81_93 = HREBZML(IGER,NOX,IV,IFDS)
              ELSE
C
C  High NOx start emissions (use normal zml and dr)
C
                 ZPOINT = NSEBZML(IGER,NOX,IV,IFDS)
                 SLOPE1 = NSEBDR(IGER,NOX,IV,IFDS)
                 EF81_93  = ZPOINT + SLOPE1*KMILES
              END IF
           END IF
C
         CASE DEFAULT
C
C  DEFAULT CASE is FOR HC and CO.
C
C  Normal running emissions
C
           IF(INH.EQ.NORMAL.AND.ISR.EQ.RUN) THEN
              ZPOINT = NREBZML(IGER,IP,IV,IFDS)
              SLOPE1 = NREBDR(IGER,IP,IV,IFDS)
              EF81_93  = ZPOINT + SLOPE1*KMILES
           ELSE
C
C  Normal start emissions
C
             IF(INH.EQ.NORMAL) THEN
               ZPOINT = NSEBZML(IGER,IP,IV,IFDS)
               SLOPE1 = NSEBDR(IGER,IP,IV,IFDS)
               EF81_93  = ZPOINT + SLOPE1*KMILES
             ELSE
C
C  High running emissions
C
               IF(ISR.EQ.RUN) THEN
                  EF81_93 = HREBZML(IGER,IP,IV,IFDS)
               ELSE
C
C  High start emissions
C
                 EF81_93 = HSEBZML(IGER,IP,IV,IFDS)
               END IF
             END IF
           END IF
      END SELECT
      END IF!not user_rec
C
   99 RETURN
      END
