      REAL FUNCTION EF_HDV(MY,IP,IV,VMTAGE,INERR)
C
C  EF_HDV returns the basic emission factor Uncorrected for
C  operating mode, CO temperature offset, methane offset,
C  temperature, inspection / maintenance, tampering offset 
C  and reformulated fuels, used for Heavy Duty Vehicles.
C
C  Called by EF_CASE.
C
C  Calls IERPTR
C
C  Changes:
C
C  15 APR 02 AIR Task Bug Fix 388: modification of constants via QUITER.
C  26 Jun 01 @EPA-bag BUG325 restructured LOGIC FOR HEAVY DUTY 2007,
C            2005, and tier2 RULEs
C  23 Mar 01 @EPA-bag BUG320 Added LOGIC FOR HEAVY DUTY 2007 RULE
C  19 Jan 01 @EPA-bag Added HDR2007 for no 2007 rule for heavy duty
C            vehicle emissions reductions
C   5 Jun 00 @EPA-ddj Bug125, Added conditional to restrict MDPV
C            calculation to only HDGVs.
C  30 May 00 @EPA-ddj Bug115, Set BERs for HDGV2b from BD04 to
C            2005-2007 BERs when NO Tier 2 Standard.
C  24 May 00 @EPA-ddj E9 Added ZML and DR data for Tier 2 basic
C            emissions and modified algorithm logic.
C  15 Sep 99 @DynTel-HXQ 2-680 Removed kink code
C  24 Aug 99 @DynTel-bg  2-660 Renamed from BEFM5.FOR to EF_HDV.FOR
C  22 Jun 99 @DynTel-MLA 2-694  Removed LEVFLG and all LEV variables.
C  11 Mar 99 @DynTel-rwk 2-694 Removed ILEV from function parameter
C            list as well as reference to P94REF array.
C  17 Feb 1999 @DynTel-rwk 2-668 Pre-81 code, removed MOBILE5 50K
C                                deterioration code.
C  21 January 99 @DynTel-ZK 2-000 Explicit Typing
C  14 Sep 98 @Dyntel-RJD 2-660 Removed '0' argument list in IERPTR
C  11 Aug 98 @Dyntel-RJD 2-663 Removed IFDS dependence from IERPTR.
C
C  Input on call or from calls:
C
C    argument list: MY,IP,IV,VMTAGE
C    common blocks:
C    /BASEQ1/ ERBZML,ERBDR
C    /BASEQ5/ ERUZML,ERUDR
C    /FLAGS2/ NEWFLG
C    /REGION/ IREJN
C    /VVNAMS/ VVLDGT12,VVLDGT34  (IVTYPE.I)
C    /VVTYPS/ VTLDGV  (IVTYPE.I)
C
C  Output on return:
C
C    function: EF_HDV
C    common blocks:
C
C  Local variable / array dictionary:
C
C   Name   Type              Description
C  ------  ----  -----------------------------------------------------
C  IGER     I    MY pointer into emission slope and intercept
C  LMY      I    Local model year can be changed from MY depending on
C                the value of NEWFLG.
C  SLOPE1   R    slope of ef curve; if 50K kink, then only for 1st 50K
C                miles
C  ZPOINT   R    zero mileage point = intercept of ef curve
C
C  Notes:
C
C 30 Sep 96 MLA 2-622  Added IVTYPE.I for Mobile6 and modified the code to
C            use the symbolic vehicle type names and vehicle type vectors
C            that were introduced in Mobile6.
C October-15-1997 @DYNTEL-id Request 2-643.Calculation of BEF for TIER1 
C           and for LEVs was moved to the function EF_DSPTCH. 
C 07 Apr 98 @DynTel-RJD 2-654 EF_HDV.FOR is called EF_DSPTCH in Mobile5.
C           EF_DSPTCH in Mobile6 calls EF_HDV in default cases where there
C           is no seperate start and running emissions data.
C 23 Jul 98 @DynTel-RJD 2-660 Expanded argument list in IGRPTR function
C           to include FDS technolgy. 
C 23 Jul 98 @DynTel-MLA 2-000  Changed LDGT names to conform to
C           certification naming conventions.
C
       IMPLICIT NONE
       INCLUDE 'BASEQ1.I'
       INCLUDE 'BASEQ5.I'
       INCLUDE 'FLAGS2.I'
       INCLUDE 'FLAGS5.I'
       INCLUDE 'REGION.I'
       INCLUDE 'IVTYPE.I'
C
C      Declare external functions
C
       INTEGER IERPTR
C
C      Declare parameter list
C
      INTEGER, INTENT(IN) :: MY
      INTEGER, INTENT(IN) :: IP
      INTEGER, INTENT(IN) :: IV
      INTEGER, INTENT(INOUT) :: INERR
      REAL,    INTENT(IN) :: VMTAGE
C
C      Declare local variables
C
      REAL    Z2B, S2B
      REAL    ZMDPV, SMDPV
      REAL    SLOPE1
      REAL    ZPOINT
      INTEGER LMY
      INTEGER KEYER
      INTEGER IGER
      INTEGER IVI
      REAL    ,SAVE :: PRCNT_MDPV=0.155
      INTEGER, SAVE :: IYES=1
      INTEGER, SAVE :: INO =0
C
C    The following emission parameters include effects of the 2005 heavy duty
C    vehicle rule and Tier2 vehicle standards; the heavy duty gas 2b's have the
C    MDPV (tier2) emission parameters weighted together with the HDG2B emission
C    parameters
       REAL, SAVE, DIMENSION(9,3,2,3)    :: ZML = RESHAPE(
C  NMOG Low Altitude
C        hdg2b  hdg3   hdg4   hdg5   hdg6   hdg7   hdg8a  hdg8b  gas bus
     &(/ 0.118, 0.104, 0.096, 0.113, 0.111, 0.123, 0.131, 0.131, 0.153,  !2005-2007
     &   0.112, 0.104, 0.096, 0.113, 0.111, 0.123, 0.131, 0.131, 0.153,  !2008
     &   0.106, 0.104, 0.096, 0.113, 0.111, 0.123, 0.131, 0.131, 0.153,  !2009+
C  NMOG High Altitude
     &   0.219, 0.193, 0.179, 0.209, 0.207, 0.228, 0.243, 0.243, 0.283,
     &   0.208, 0.193, 0.179, 0.209, 0.207, 0.228, 0.243, 0.243, 0.283,
     &   0.197, 0.193, 0.179, 0.209, 0.207, 0.228, 0.243, 0.243, 0.283,
C  CO Low Altitude
     &   7.782, 8.165, 8.051, 9.400, 9.308,10.267,10.934,15.430,12.750,
     &   7.051, 8.165, 8.051, 9.400, 9.308,10.267,10.934,15.430,12.750,
     &   6.320, 8.165, 8.051, 9.400, 9.308,10.267,10.934,15.430,12.750,
C  CO High Altitude
     &  24.761,25.981,25.620,29.912,29.618,32.668,34.792,49.098,40.569,
     &  22.436,25.981,25.620,29.912,29.618,32.668,34.792,49.098,40.569,
     &  20.112,25.981,25.620,29.912,29.618,32.668,34.792,49.098,40.569,
C  NOx Low Altitude
     &   0.573, 0.594, 0.578, 0.675, 0.669, 0.737, 0.785, 0.785, 0.916,
     &   0.534, 0.594, 0.578, 0.675, 0.669, 0.737, 0.785, 0.785, 0.916,
     &   0.495, 0.594, 0.578, 0.675, 0.669, 0.737, 0.785, 0.785, 0.916,
C  NOx High Altitude
     &   0.469, 0.486, 0.473, 0.552, 0.547, 0.603, 0.642, 0.642, 0.749,
     &   0.437, 0.486, 0.473, 0.552, 0.547, 0.603, 0.642, 0.642, 0.749,
     &   0.405, 0.486, 0.473, 0.552, 0.547, 0.603, 0.642, 0.642, 0.749/)
     &  ,(/ 9,3,2,3 /))
C
       REAL, SAVE, DIMENSION(9,3,2,3)    :: DR = RESHAPE(
C  NMOG Low Altitude
     &(/ 0.008, 0.006, 0.006, 0.007, 0.007, 0.007, 0.008, 0.008, 0.009,
     &   0.008, 0.006, 0.006, 0.007, 0.007, 0.007, 0.008, 0.008, 0.009,
     &   0.007, 0.006, 0.006, 0.007, 0.007, 0.007, 0.008, 0.008, 0.009,
C  NMOG High Altitude
     &   0.015, 0.012, 0.011, 0.012, 0.012, 0.013, 0.014, 0.014, 0.017,
     &   0.014, 0.012, 0.011, 0.012, 0.012, 0.013, 0.014, 0.014, 0.017,
     &   0.013, 0.012, 0.011, 0.012, 0.012, 0.013, 0.014, 0.014, 0.017,
C  CO Low Altitude
     &   0.279, 0.293, 0.289, 0.338, 0.334, 0.369, 0.393, 0.554, 0.458,
     &   0.253, 0.293, 0.289, 0.338, 0.334, 0.369, 0.393, 0.554, 0.458,
     &   0.227, 0.293, 0.289, 0.338, 0.334, 0.369, 0.393, 0.554, 0.458,
C  CO High Altitude
     &   0.889, 0.932, 0.920, 1.074, 1.064, 1.173, 1.250, 1.763, 1.457,
     &   0.806, 0.932, 0.920, 1.074, 1.064, 1.173, 1.250, 1.763, 1.457,
     &   0.722, 0.932, 0.920, 1.074, 1.064, 1.173, 1.250, 1.763, 1.457,
C  NOx Low Altitude
     &   0.008, 0.008, 0.008, 0.009, 0.009, 0.010, 0.011, 0.011, 0.013,
     &   0.007, 0.008, 0.008, 0.009, 0.009, 0.010, 0.011, 0.011, 0.013,
     &   0.007, 0.008, 0.008, 0.009, 0.009, 0.010, 0.011, 0.011, 0.013,
C  NOx High Altitude
     &   0.007, 0.007, 0.006, 0.008, 0.008, 0.008, 0.009, 0.009, 0.010,
     &   0.006, 0.007, 0.006, 0.008, 0.008, 0.008, 0.009, 0.009, 0.010,
     &   0.006, 0.007, 0.006, 0.008, 0.008, 0.008, 0.009, 0.009, 0.010/)
     &  ,(/ 9,3,2,3 /))
C
C  Look up base and deterioration rates.  Pointer function IERPTR
C  sets both the arrays to access (block data default or user
C  supplied) and the cell pointers to the coefficients.  KINK50
C  determines if a second slope is needed. Currently, only 1991+
C  LDGT have 50K+ kinked EFs.
C
C  while if NEWFLG=5, No Clean Air Act (CAA) then set the
C  Model Year (MY) to 1993.( so ALL EF's are Pre-1993)
C
      IF (NEWFLG.GE.5 .AND. MY.GE.1994) THEN
        LMY = 1993
      ELSE
        LMY = MY
      END IF
C
C   Slopes and intercept for the non light duty vehicles.
C
      IGER = IERPTR(1,LMY,IP,IV,KEYER,INERR)
C
C   test for user input
      IF(KEYER/=1) THEN
C   Use user input
        ZPOINT = ERUZML(IGER,IP,IV,IREJN)
        SLOPE1 = ERUDR (IGER,IP,IV,IREJN)
        GOTO 99
      END IF
C
C
      IVI = IV - 5
      IF(IVI .GT. 9) IVI = 9
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C DIESEL VEHICLES
C
      IF(VVHDDV(IV)==IYES) THEN
C       
        SELECT CASE(MY)
C
C         Use BD04 block data
C
          CASE(:2003)
            ZPOINT = ERBZML(IGER,IP,IV,IREJN)
            SLOPE1 = ERBDR (IGER,IP,IV,IREJN)
C
          CASE(2004:2006)
C          Use BD04 for HDDVs.
C          If the 2005 rule (NO TIER2) for heavy duty vehicles is not
c          in effect use the emission line slopes and intercepts which do not
C          take the heavy duty rule into account.
C
           ZPOINT = ERBZML(IGER,IP,IV,IREJN)
           SLOPE1 = ERBDR (IGER,IP,IV,IREJN)
C
          CASE(2007:)
C          If the 2007 rule for heavy duty vehicles is not in effect
C          use the emission line slopes and intercepts which do not
C          take the heavy duty rule into account.
           LMY=2006
           IF(HDR2007==INO) IGER=IERPTR(1,LMY,IP,IV,KEYER,INERR) !NO 2007 rule(2007+)
           ZPOINT = ERBZML(IGER,IP,IV,IREJN)
           SLOPE1 = ERBDR (IGER,IP,IV,IREJN)
C
        END SELECT
      END IF
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C GASOLINE VEHICLES
C
      IF(VVHDGV(IV)==IYES) THEN
C       
        SELECT CASE(MY)
          CASE(:2004)
C
C           Use BD04 block data
C
            ZPOINT = ERBZML(IGER,IP,IV,IREJN)
            SLOPE1 = ERBDR (IGER,IP,IV,IREJN)
C
          CASE(2005:2007)
C
C           Use BD04 for HDGVs when Tier 2 not in effect.
            LMY=2004
            IF(TIER2_FLAG==INO) IGER=IERPTR(1,LMY,IP,IV,KEYER,INERR) !NO 2005 rule(2005-2007)
            ZPOINT = ERBZML(IGER,IP,IV,IREJN)
            SLOPE1 = ERBDR (IGER,IP,IV,IREJN)
C
C           Use the heavy duty gas 2005 rule emission line parameters which are defined
C           above in this subroutine (the Tier2 MDPV's and 2005 rule HDG2B's are weighted
C           together)
            IF(TIER2_FLAG==IYES) THEN
              LMY = MY - 2006
              LMY = MAX(LMY,1)
              IF(LMY .GT. 3) LMY = 3
              ZPOINT = ZML(IVI,LMY,IREJN,IP)
              SLOPE1 = DR (IVI,LMY,IREJN,IP)
            END IF
C
          CASE(2008:)
C
            IF(TIER2_FLAG==INO.AND.HDR2007==INO) THEN
              LMY=2004
              IGER = IERPTR(1,LMY,IP,IV,KEYER,INERR) !NO 2005 rule(2005-2007) & NO 2007 RULE(2008+)
              ZPOINT = ERBZML(IGER,IP,IV,IREJN)
              SLOPE1 = ERBDR (IGER,IP,IV,IREJN)
            ENDIF
C
            IF(TIER2_FLAG==IYES.AND.HDR2007==INO) THEN
C             Use the heavy duty gas 2005 rule emission line parameters which
C             are defined above in this subroutine (the Tier2 MDPV's and 2005
C             rule HDG2B's are weighted together)
              LMY = MY - 2006
              IF(LMY .GT. 3) LMY = 3
              ZPOINT = ZML(IVI,LMY,IREJN,IP)
              SLOPE1 = DR (IVI,LMY,IREJN,IP)
            ENDIF
C
            IF(TIER2_FLAG==IYES.AND.HDR2007==IYES) THEN
C             Use BD04 heavy duty gas 2007 rule emission line parameters (the Tier2
C             MDPV's and 2007 rule HDG2B's are weighted together)
              IF(IP==2) THEN
                LMY = MY - 2006
                IF(LMY .GT. 3) LMY = 3
                ZPOINT = ZML(IVI,LMY,IREJN,IP)
                SLOPE1 = DR (IVI,LMY,IREJN,IP)
              ELSE
                ZPOINT = ERBZML(IGER,IP,IV,IREJN)
                SLOPE1 = ERBDR (IGER,IP,IV,IREJN)
              ENDIF
            ENDIF
C
            IF(TIER2_FLAG==INO.AND.HDR2007==IYES) THEN
C
              IF(IV==6) THEN
C               Use BD04 2004 MDPV emission line parameters without 2005/tier2 rule
C               weighted with HDGV2B emission line parameters with 2007 rule reductions
C
C               HDGV2B's :
                IF(MY==2008) THEN
                  Z2B = 0.0709  !these number are found by subtracting the tier2 MDPV
                  S2B = 0.00674 !portion from the 2009 emission line parameters from BD04.for
                ELSE
                  Z2B = 0.0865
                  S2B = 0.00552
                ENDIF
C               MDPV's :
                LMY=2004
                IGER = IERPTR(1,LMY,IP,IV,KEYER,INERR) !NO 2005/MDPV TIER2 rule
                ZMDPV = ERBZML(IGER,IP,IV,IREJN)
                SMDPV = ERBDR (IGER,IP,IV,IREJN)
C
                ZPOINT = PRCNT_MDPV*ZMDPV + (1-PRCNT_MDPV)*Z2B
                SLOPE1 = PRCNT_MDPV*SMDPV + (1-PRCNT_MDPV)*S2B
              ELSE
C               Use BD04 heavy duty gas 2007 rule emission line parameters
                ZPOINT = ERBZML(IGER,IP,IV,IREJN)
                SLOPE1 = ERBDR (IGER,IP,IV,IREJN)
              ENDIF
C
            ENDIF

        END SELECT
      ENDIF
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
  99  EF_HDV = ZPOINT + SLOPE1*VMTAGE
C
      RETURN
      END
