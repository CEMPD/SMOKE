      SUBROUTINE REFUEL(ICY,INERR) 
C 
C  REFUEL computes the refueling losses HC in grams / mile by relative 
C  (IDX ordered) model year for the 4 vehicle types affected. 
C 
C  Called by EFCALX. 
C 
C  Calls QUITER 
C 
C  Changes: (Last change first) 
C 
C 06 Nov 01 AIR Task 27: Updated all fuel economy calculations via PMFLECON. 
C 06 Mar 01 @EPA-djb Applied CAAFLG to ORVR refueling emissions by 
C           setting RLFLAG.  This replaces NEWFLG completely. 
C 26 Feb 01 @EPA-djb Added HDGV2b to the ORVR calculations (OBIPY & OBPIP). 
C 18 Jul 00 @EPA-bag bug#54 removed TEMFLG 
C 23 Jun 00 @EPA-bag bug153 Removed call to DBREFL; refueling losses are 
C           written to database output in DBEVAP 
C 08 Jun 00 AIR Task 03: Removed nonblock DO constructs. 
C  6 Apr 00 AIR Task 02: Removed warning message for ob vrs regulations (1994) 
C 10 Aug 99 @DynTel-MLA 2-684  Added INERR to the parameter list, added 
C           a call to the DBREFL routine, and added the USE DATABASE line. 
C  2 Jun 99 @DynTel-MLA 2-694  Removed obsolete variable OBDF. 
C 26 Apr 99 @DynTel-bg 2-000 Explicit Typing 
C  1 Mar 99 @DynTel-DDJ 2-696 Added hourly temperatures 
C 
C  Input on call: 
C 
C    Constants:  
C    MAXVEH from include file IVTYPE.I 
C 
C    argument list: ICY 
C 
C    common blocks: 
C    /CITRV1/ RVPUWX 
C    /FLAGS3/ RLFLAG, CAAFLG 
C    /MAXIMA/ MAXYRS 
C    /RLCOM1/ SPILL,OBED,OBES,OBPIP,OBPIY 
C    /RLCOM2/ IS2SY,NPHASE,S2EFF(4) 
C    /TAMOU2/ TOB 
C    /TEMPS/  TEMMIN,TEMMAX,TEMEVP 
C    /VVNAMS/ VVGAS, VVLDGAS  (IVTYPE.I) 
C 
C  Output on return: 
C 
C    common blocks: 
C    /RLCOM3/ RLRATE 
C 
C  Local array subscripts: 
C 
C  S2LEFT(4)       -  S2LEFT ( IVTAM ) 
C 
C  Local variable dictionary: 
C 
C  Name   Type              Description 
C  ------  ----  ------------------------------------------------------- 
C  ICY      I    user supplied calendar year 
C  IDX      I    yearly loop index (1-25) 
C  MY            I    Vehicle model year; MY=ICY+IDX-25 
C  IV       I    vehicle type loop index 
C  IVTAM    I    local vehicle type loop index (gasoline vehicles only) 
C  IH       I    hourly loop index 
C  INERR    I    error count 
C  IALLYR   I    all stations year: year by which phase in of Stage II vrs 
C                controls is completed => 100% of volume of gas passes through 
C                Stage II vrs. 
C  IOBCHK   I    check for ob vrs installation, given ob vrs program exists: 
C                    1 = no ob vrs    2 = ob vrs installed 
C  DFTEMP   R    dispensed fuel temperature (F). 
C  RVP      R    RVP of the fuel in use that calendar year 
C  S2LEFT   R    percentage left of uncontrolled refueling losses after Stage 
C                II controls are applied. 
C  S2VOL    R    percentage of gasoline pumped that has S2 controls applied. 
C  TDFDIF   R    difference in temperature of fuel in vehicle 
C                tank and dispensed fuel 
C  IDXOB    I    Index to the ob vrs phase-in tables 
C  I        I    index for OBPIY, Onboard VRS phase-in years 
C  RFE      I    Road fuel economy 
C 
C  Notes: 
C 
C  RLFLAG impact: 
C 
C    1 or (3 & my/iv  -  maps to RLFLAG = 3 
C    not in OB program) 
C 
C    2 or (4 & my/iv  -  maps to RLFLAG = 4 
C    not in OB program) 
C 
C    3 or 4 and my/iv -  weight uncontrolled & controlled rate by best estimate 
C    in OB program       of ob vrs tamper rate & its complement, respectively, 
C                        to get ob only vrs controlled refueling losses. 
C                        Assume 100% installation. 
C 
C    4                -  compute and assign S2 & ob vrs controlled losses. 
C                        See above. 
C 
C    5                -  assign 0.0 (no refueling losses case). 
C 
C  The factors returned are weighted by their corresponding TF's and summed 
C  across IDX to get the refueling loss for ICY for each IVTAM. 
C 
C  REFUEL was changed in MOBILE4.1 to use 
C  an equation to compute Displacement. 
C 
C  5-May-1994 @ CSC-tmm Request 2-411 
C    Code added to implement 1994 Onboard Refueling Vapor Recovery Regulations 
C 
C  21-Nov-1994 @ CSC-ked Request 0-411 
C    Allow printing of Warning message 154 only once for July runs 
C  09-June-1995 @ CSC-wcs Request # 0-411 ORVR Fix (Allow printing of warning 
C    message 156 only once for July runs.) 
C 
C  08-Feb-1996 @ DynTel Request 2-615 ORVR effects are disabled when NEWFLG  
C                is set to 5 or 6. 
C 
C  11 Jun 97 DynTel@MLA 2-622  Added IVTYPE.I for Mobile6 and modified the 
C            code to use the symbolic vehicle type names and vehicle type 
C            vectors that were introduced in Mobile6. 
C 
C  26 Feb 97 @Dyntel-rjd (after gjr) Request 2-204 Removed IY from parameter  
C            list and changed code described above in ORVR Fix.  This is due 
C            to a change in the way July runs will be calculated.  
C            Interpolation between two consecutive Januaries will no longer  
C            be used.  Thus July runs will call this subroutine only once  
C            per scenario thereby eliminating the problem described in request  
C            0-411. 
C 
C  20 Apr 98 @Dyntel-bag Request 2-661 OUTFMT dependencies have been deleted.  
C 
      USE DATABASE, ONLY: DBFLAG, DBSELEFT 
      IMPLICIT NONE 
C 
      INCLUDE 'CITRV1.I' 
      INCLUDE 'FLAGS2.I' 
      INCLUDE 'FLAGS3.I' 
      INCLUDE 'IVTYPE.I' 
      INCLUDE 'MAXIMA.I' 
      INCLUDE 'RLCOM1.I' 
      INCLUDE 'RLCOM2.I' 
      INCLUDE 'RLCOM3.I' 
      INCLUDE 'TAMOU2.I' 
      INCLUDE 'TEMPS.I' 
C 
C     Declare called functions 
C 
      REAL, EXTERNAL :: PMFLECON 
C 
C     Declare parameter list 
C 
      INTEGER, INTENT(IN) :: ICY 
      INTEGER, INTENT(INOUT) :: INERR 
C 
C     Declare local variables 
C 
      INTEGER IV 
      INTEGER IH 
      INTEGER I 
      INTEGER IDX 
      INTEGER IDXOB 
      INTEGER IALLYR 
      INTEGER MY 
      INTEGER IOBCHK 
      INTEGER IVTAM 
      REAL    S2LEFT(14) 
      REAL    RVP 
      REAL    S2VOL 
      REAL    DFTEMP 
      REAL    TDFDIF 
      REAL    RFE 
C 
C  Initialize 
C 
      S2LEFT = 1.0 
      RLRATE = 0.0 
C 
C  Refueling losses are zero if the no losses case is selected. 
C 
      IF(RLFLAG.EQ.5) GOTO 99 
C 
C   If Clean Air Act is disabled (CAAFLG=1), 
C     disable Onboard Vapor Recovery. 
C 
C     RLFLAG is never 3 or 4, but check anyway. 
C 
      IF (CAAFLG.EQ.1) THEN 
        CALL QUITER(0.,0,156,INERR) 
C 
        IF (RLFLAG.EQ.3) THEN 
          RLFLAG = 1 
        ELSEIF (RLFLAG.EQ.4) THEN 
          RLFLAG = 2 
        ENDIF 
C 
      ELSE 
C 
C  03/25/94 TMM@CSC REQUEST #2-411 
C    Turn on onboard vapor recovery unconditionally 
C      and implement the phase-in scheme in OBPIY and OBPIP 
C 
C    RLFLAG value 1 now maps to 3; 2 maps to 4; 
C 
        IF (RLFLAG.EQ.1) THEN 
          RLFLAG = 3 
        ELSEIF (RLFLAG.EQ.2) THEN 
          RLFLAG = 4 
        ENDIF 
C 
      ENDIF 
C 
C  If Stage II vrs requested, compute the corresponding cf by IVTAM. 
C 
      IF(RLFLAG.EQ.1.OR.RLFLAG.EQ.3) GOTO 20 
C 
      IALLYR=IS2SY+NPHASE 
      IF(ICY.LT.IS2SY) S2VOL=0.0 
      IF(ICY.GE.IS2SY.AND.ICY.LT.IALLYR) S2VOL=(ICY-IS2SY)/(NPHASE*1.0) 
      IF(ICY.GE.IALLYR) S2VOL=1.0 
C 
      S2LEFT(1)=1.0-S2EFF(1)*S2VOL 
      S2LEFT(2)=S2LEFT(1) 
      S2LEFT(3)=S2LEFT(2) 
      S2LEFT(4)=S2LEFT(3) 
      S2LEFT(5)=S2LEFT(4) 
      DO 15 IVTAM=6,14 
        S2LEFT(IVTAM)=1.0-S2EFF(IVTAM)*S2VOL 
   15 CONTINUE 
C 
C  Set and limit the Dispensed Fuel Temperature (F). 
C 
   20 CONTINUE 
      DO 60 IH = 1, 24 
        DFTEMP = TEMHRLY(IH) 
C 
        IF(DFTEMP.LT.20.0) DFTEMP=20.0 
        IF(DFTEMP.GT.95.0) DFTEMP=95.0 
C 
C  Set and limit the difference in temperature of the 
C  fuel in vehicle tank and dispensed fuel. 
C 
        TDFDIF=0.418*DFTEMP-16.6 
C 
        IF(TDFDIF.GT.20.0) TDFDIF=20.0 
C 
C  Set the displacement RVP. 
C 
        RVP=RVPUWX 
C 
C  Calculate and limit the Displacement (g/gal). 
C 
        DISPL=-5.909-0.0949*TDFDIF+0.0884*DFTEMP+0.485*RVP 
C 
        IF(DISPL.LT.1.80) DISPL=1.80 
C 
C  REFUEL is only called once per scenario, computing and saving for later 
C  use all IVTAM * IDX onboard rl offsets.  If RLFLAG = 3, the results are 
C  program applies. 
C 
        IVTAM = 0 
        DO 40 IV=1,MAXVEH 
C 
          IF (VVGAS(IV).EQ.0) CYCLE 
          IVTAM = IVTAM + 1 
C 
          DO 30 IDX=1,MAXYRS 
C 
C  Set the road fuel economy. 
C 
            MY=ICY+IDX-MAXYRS 
C 
            RFE=PMFLECON(MY,IV,INERR) 
C 
C  Set the onboard refueling vapor recovery (ORVR) system flag (IOBCHK), 
C    1 : No ORVR system 
C    2 : With ORVR system 
C    and select the system phase in year index (IDXOB). 
C 
            IDXOB = 0 
            IOBCHK = 1 
            IF (VVLDGAS(IV).EQ.1 .OR. IV.EQ.6) THEN 
              DO 50 I=1,3 
                IF (MY.GE.OBPIY(I,IVTAM)) THEN 
                  IDXOB=I 
                  IOBCHK=2 
                ENDIF 
   50         CONTINUE 
            ENDIF 
C 
C  RLFLAG = 1: Use uncontrolled refueling loss. 
C 
          IF(RLFLAG.EQ.1.OR.(RLFLAG.EQ.3.AND.IOBCHK.EQ.1)) THEN 
            RLRATE(IDX,IVTAM,IH)=(DISPL+SPILL)/RFE 
C 
C  RLFLAG = 2: Calculate Stage II only vrs controlled loss. 
C 
          ELSE IF((RLFLAG.EQ.2.OR.(RLFLAG.EQ.4.AND.IOBCHK.EQ.1))) THEN 
            RLRATE(IDX,IVTAM,IH)=(S2LEFT(IVTAM)*DISPL+SPILL)/RFE 
C 
C  RLFLAG = 3: Calculate an onboard only vrs controlled refueling loss. 
C 
C  Use the evap cannister only (ID = 6) probability of tampering + 
C  100% of the given IV x MY fleet is equipped.  Obtain the onboard loss by 
C  weighting the uncontrolled and controlled refueling loss rates by the 
C  tampering rate and its complement, respectively. 
C 
C  Weight the calculation by the phase-in of onboard refueling vapor 
C  recovery regulations of 1994. 
C 
          ELSE IF(RLFLAG.GE.3.AND.IOBCHK.EQ.2) THEN 
            RLRATE(IDX,IVTAM,IH) = 
     *    (DISPL/RFE)*((1-OBED)*(1-TOB(IDX,IVTAM)) 
     *      *OBPIP(IDXOB,IVTAM) + TOB(IDX,IVTAM) 
     *      *OBPIP(IDXOB,IVTAM)*S2LEFT(IVTAM) 
     *      + (1 - OBPIP(IDXOB,IVTAM))*S2LEFT(IVTAM) ) 
     *      + (SPILL/RFE)*((1 - OBES) 
     *      *OBPIP(IDXOB,IVTAM) + (1 - OBPIP(IDXOB,IVTAM))) 
          ENDIF 
C 
C  RLFLAG = 4: Calculate S2 & ob vrs controlled refueling loss: 
C 
C  Handled above.  Reduces to either RLFLAG=2 or RLFLAG=3 case, 
C  depending on IOBCHK.  Onboard takes precedence over Stage II and 
C  you can't get credit for both. 
C 
C 
   30     CONTINUE 
   40   CONTINUE 
   60 END DO 
C 
   99 RETURN 
      END 
