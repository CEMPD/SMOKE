      SUBROUTINE PMNH3(ICY) 
C 
C  Computes the Ammonia (NH3) emission rates. 
C 
C  Changes : 
C 
C  17 Apr 01 AIR Task 11: Changed code to handle new FVMT dimensions. 
C  01 Feb 01 AIR Task 11: New Module 
C 
C  Called by PMCALX. 
C 
C  Calls P3WPOX 
C 
C  Parameter list: ICY 
C 
C  Common Blocks: See Includes below. 
C 
C  Local variable / array dictionary: 
C 
C   Name        Type              Description 
C  ------       ----  ----------------------------------------------------- 
C  BER           R    basic emission rate 
C  BERCAT        R    BER weighted by CAT 
C  BERROAD       R    BER weighted by ROAD 
C  FCAT          R    Fraction no catalyst (1), 3-way (2), ox (3) 
C  ICAT          I    Loop for catalyst-equipped, 1=without, 2=with 
C  IDX           I    IDX=(MAXYRS+1)-JDX 
C  IH            I    Hour index 
C  IROAD         I    Raod index 
C  IV            I    Vehicle class 
C  JDX           I    Age 
C  MY            I    Model year 
C  P3W           R    Fraction of 3-way catalysts 
C  POX           R    Fraction of oxidation catalysts 
C  SULF          R    fuel sulfur level 
C  ZEVADJ        R    Adjustment for ZEVs 
C 
      IMPLICIT NONE 
C 
      INCLUDE 'IVTYPE.I'  ! MAXVEH, VVGASMC 
      INCLUDE 'MAXIMA.I'  ! MAXYRS 
      INCLUDE 'MYRCAL.I'  ! TF 
      INCLUDE 'PART1.I'   ! D_PPM 
      INCLUDE 'PART6.I'   ! NH3BER, NH3SDR 
      INCLUDE 'PART3.I'   ! EFPM, DBPM, VPM 
      INCLUDE 'SPEED9.I'  ! FVMT, HVMT 
      INCLUDE 'SULFUR.I'  ! S_PPM 
      INCLUDE 'VMXCOM.I'  ! VMTMIX 
C 
C  Declare parameter list. 
C 
      INTEGER ICY 
C 
C  Declare external functions. 
C 
      REAL, EXTERNAL :: PCTLEV 
C 
C  Declare locate variables/arrays 
C 
      INTEGER ICAT 
      INTEGER IDX 
      INTEGER IH 
      INTEGER IROAD 
      INTEGER IV 
      INTEGER JDX 
      INTEGER MY 
      REAL    BER 
      REAL    BERCAT 
      REAL    BERROAD 
      REAL    FCAT(3) 
      REAL    P3W 
      REAL    POX 
      REAL    SULF 
      REAL    ZEVADJ 
C 
C 
C  Loop for each vehicle class 
C 
      DO IV=1,MAXVEH 
C 
C  Skip if VMTMIX is zero 
C 
      IF(VMTMIX(IV).EQ.0.) CYCLE 
C 
C  Assign the sulfur level (gas or diesel). 
C  Sulfur level must be between 30 and 500 ppm. 
C 
        IF(VVGASMC(IV).EQ.1) THEN 
          SULF=S_PPM 
        ELSE 
          SULF=D_PPM 
        ENDIF 
C 
        SULF=MIN(500.,SULF) 
        SULF=MAX( 30.,SULF) 
C 
C  Loop for all ages 
C 
        DO JDX=1,MAXYRS 
C 
C  Set IDX 
C 
          IDX=(MAXYRS+1)-JDX 
C 
C  If travel fraction for this age/vehicle is 0, skip. 
C 
          IF(TF(IDX,IV) .EQ. 0.) CYCLE 
C 
C  Set model year 
C 
          MY=ICY-JDX+1 
C 
C  Determine the ZEV adjustment for LDGV, LDGT1-LDGT4. 
C  (NGVs are assummed to emit NH3 PM like their gasoline counterparts, 
C   so no other adjustments are needed.) 
C 
          ZEVADJ = 1.0 - PCTLEV(MY,IV,11) 
C 
C  Determine the catalyst fractions (1=no, 2=ox, 3=3-way) 
C 
          CALL P3WPOX(MY,IV,P3W,POX) 
C 
          FCAT(1)=1.0-P3W-POX 
          FCAT(2)=POX 
          FCAT(3)=P3W 
C 
C  COMPUTE THE NH3 BERS 
C 
C  Initialize BER weighted by road and catalyst 
C 
          BERROAD=0. 
C 
C  Loop by road type 
C 
          DO IROAD=1,4 
C 
C  Initialize BER weighted by catalyst 
C 
            BERCAT=0. 
C  
C  Loop for catalyst 
C 
            DO ICAT=1,3 
C 
              IF(FCAT(ICAT).EQ.0.) CYCLE 
C 
C  Compute the basic emission rates (BER).  
C 
              BER = NH3BER(ICAT,IROAD,IV) 
     *            + NH3SDR (ICAT,IROAD,IV) 
     *            * SULF 
C 
C  Build the catalyst weighted BER 
C 
              BERCAT = BERCAT + BER * FCAT(ICAT) 
C 
            END DO !ICAT 
C 
C  Adjust for ZEVs 
C 
            BERCAT = BERCAT * ZEVADJ 
C 
C  Loop by hour 
C 
            DO IH=1,MAXIH 
C 
C  Store the values for database.  
C 
              DBPM(7,IROAD,IH,JDX,IV) = BERCAT 
C 
C  Build the road weighted BER 
C 
              BERROAD = BERROAD + BERCAT * FVMT(IROAD,IH,IV) * HVMT(IH) 
C 
            END DO  !IH 
           
          END DO  !IROAD 
C 
          EFPM(7,IV) = EFPM(7,IV) + BERROAD * TF(IDX,IV) 
C 
        END DO  !JDX 
C 
        VPM(7) = VPM(7) + EFPM(7,IV) * VMTMIX(IV) 
C 
      END DO  !IV 
C 
      RETURN 
      END 
