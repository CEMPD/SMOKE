      SUBROUTINE PMDSO4SO2(ICY,INERR) 
C 
C  This routine computes the diesel SO4 and SO2 based on the PART5 model. 
C  It uses the BERs associated sulfur level from PMEMFAC so that the 
C  associated SO4 can be subtracted from the EMFAC-based EFs. 
C 
C Changes: 
C 
C 15 Nov 01 AIR Task 27: Added INERR support due to PMFLECON update. 
C 03 May 01 AIR Task 12: Apply assocated fuel suflur 
C 14 Feb 01 AIR Task 11: New routine for MOBILE6. 
C 14 Feb 01 AIR Task 11: Removed code for indirect sulfur, 
C           buses, hardcoded vehicles. 
C 14 Nov 00 AIR Task 06: Change DO loop maximum to NIYRFE. 
C 10 Nov 00 AIR Task 06: Updated code to MOBILE6 standards. 
C 
C  Input on call: 
C 
C  Called by PMCALX 
C 
C  Calls: None 
C 
C  Parameter list: ICY 
C 
C  Common Blocks: See Includes below. 
C 
C  Local variable dictionary: 
C 
C   Name     Type             Description 
C  ------    ----  ----------------------------------------------------- 
C  DSWGTR      R   Diesel weight ratio of sulfur content in fuel, wt% 
C  FECONM      R   Fuel Economy 
C  IDX         I   IDX=(MAXYRS+1)-JDX 
C  IH          I   Hour 
C  IROAD       I   Road Type 
C  IV          I   Vehicle Class 
C  JDX         I   Age 
C  LCNVRT      R   Local version of DCNVRT for LDDV, and LDDT 
C  MY          I   Model year  
C  SO2         R   Temporary SO2 
C  SO4         R   Temporary SO4 
C  SUMSO2IH    R   Temporary summation of SO2 
C  SUMSO2JDX   R   Temporary summation of SO2 
C  SUMSO4IH    R   Temporary summation of SO4 
C  SUMSO4JDX   R   Temporary summation of SO4 
C 
C 
      IMPLICIT NONE 
C 
      INCLUDE 'FRACT.I'   ! WATER, FDNSTY, DCNVRT 
      INCLUDE 'IVTYPE.I'  ! MAXVEH, VDM6PRT, VVDSL, VTLDDV, VTLDDT12, VTLDDT34 
      INCLUDE 'MAXIMA.I'  ! MAXIH, MAXYRS 
      INCLUDE 'MYRCAL.I'  ! TF 
      INCLUDE 'PART1.I'   ! BERSULF 
      INCLUDE 'PART3.I'   ! VPM, DBPM, EFPM 
      INCLUDE 'SPEED9.I'  ! FVMT, HVMT 
      INCLUDE 'VMXCOM.I'  ! VMTMIX 
C 
C  Declare external functions. 
C 
      REAL, EXTERNAL :: PMFLECON 
C 
C  Declare parameter list. 
C 
      INTEGER ICY 
      INTEGER INERR 
C 
C  Declare local variables/arrays. 
C 
      INTEGER IDX 
      INTEGER IH 
      INTEGER IROAD 
      INTEGER IV 
      INTEGER JDX 
      INTEGER MY 
      REAL    DSWGTR 
      REAL    FECONM 
      REAL    LCNVRT 
      REAL    SO2 
      REAL    SO4 
      REAL    SUMSO2IH 
      REAL    SUMSO2JDX 
      REAL    SUMSO4IH 
      REAL    SUMSO4JDX 
C 
C 
C  Loop for each vehicle class 
C 
      DO IV=1,MAXVEH 
C 
C  If not diesel fueled, skip 
C 
      IF(VVDSL(IV).EQ.0) CYCLE 
C 
C  Skip if VMTMIX is zero 
C 
      IF(VMTMIX(IV).EQ.0.) CYCLE 
C 
C  May-24-1994 @ CSC-bsg Request 417 Changed DCNVRT to "0.015" for  
C  light duty diesel vehicles and trucks.  This will reflect the  
C  Regulatory Impact Analysis which states that 1.5% sulfur in fuel 
C  is converted to sulfate for these types of vehicles 
C       
      IF(IV.EQ.VTLDDV   .OR. 
     *   IV.EQ.VTLDDT12 .OR. 
     *   IV.EQ.VTLDDT34) THEN  
        LCNVRT=0.015 
      ELSE 
        LCNVRT=DCNVRT 
      ENDIF 
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
C  In order for the SO4 to be subtracted from the EMFAC-based EFs, 
C  the sulfur levels must agree. Use the associated fuel sulfur 
C  level from PMEMFAC and convert from ppm to wt%.  
C 
          DSWGTR = BERSULF(JDX,IV)/10000. 
C 
C  Obtain the diesel fuel economy for the model year and vehicle class. 
C  (The diesel fuel economies, in miles/gallon, are from PART5.) 
C 
          FECONM = PMFLECON(MY,IV,INERR) 
C 
C     Compute the grams/mile emission factor of the direct sulfate material (SO4): 
C  
C     13.6078 = (453.592 * 3.)/100. 
C     where: 453.592 = number of grams in a pound 
C                  3 = weight ratio of SO4 to sulfur 
C     and the division by 100 is to correct for DSWGTR the diesel weight ratio 
C     of sulfur content in fuel. 
C 
C     WATER  = weight ratio of seven water molecules to sulfate  
C              (see p. 6-25 of "Regulatory Impact Analysis, 
C              Control of Sulfur and Aromatics Contents of On- 
C              Highway Diesel Fuel") 
C     FDNSTY = fuel density in lbs/gallon, 1=gasoline, 2=diesel 
C     FECONM = fuel economy 
C     DSWGTR = Diesel weight ratio of sulfur content in fuel 
C     DCNVRT  = percent of sulfur exhausted that is directly converted 
C               into sulfuric acid. 
C 
      SO4 = 13.6078*(1.+ WATER)*FDNSTY(2)*DSWGTR*LCNVRT/FECONM 
C 
C     Compute the exhausted SO2: 
C 
C     9.072 = (453.592 * 2.)/100. 
C     where: 453.592 = number of grams in a pound 
C                 2. = weight ratio of SO2 to sulfur 
C     and the division by 100 is to correct for the user supplied or default  
C     weight percent of sulfur, DSWGTU. 
C 
      SO2 = 9.072*FDNSTY(2)*DSWGTR*(1.-LCNVRT)/FECONM 
C 
C  Loop for hour and road type.  Compute and store composite EFs. 
C  PROGRAMMING NOTE: Currently, diesel emissions are assumed to be  
C  independent of hour, road type, and speed. This may change in  
C  the future, so weighting code is kept active here. 
C 
          SUMSO4JDX = 0. 
          SUMSO2JDX = 0. 
C 
          DO IH=1,MAXIH 
C 
            SUMSO4IH = 0. 
            SUMSO2IH = 0. 
C 
            DO IROAD=1,4 
C 
C  Load the arrays. 
C 
              DBPM(1,IROAD,IH,JDX,IV) = SO4 
              DBPM(6,IROAD,IH,JDX,IV) = SO2 
C 
              SUMSO4IH = SUMSO4IH + SO4 * FVMT(IROAD,IH,IV)  !SO4 
              SUMSO2IH = SUMSO2IH + SO2 * FVMT(IROAD,IH,IV)  !SO2 
C 
            END DO !IROAD 
C 
            SUMSO4JDX = SUMSO4JDX + SUMSO4IH * HVMT(IH)  !SO4 
            SUMSO2JDX = SUMSO2JDX + SUMSO2IH * HVMT(IH)  !SO2 
C 
          END DO !IH 
C 
          EFPM(1,IV)=EFPM(1,IV)+SUMSO4JDX*TF(IDX,IV)  !SO4 
          EFPM(6,IV)=EFPM(6,IV)+SUMSO2JDX*TF(IDX,IV)  !SO2 
C 
        END DO !JDX 
C 
        VPM(1)=VPM(1)+EFPM(1,IV)*VMTMIX(IV)  !SO4 
        VPM(6)=VPM(6)+EFPM(6,IV)*VMTMIX(IV)  !SO2 
C 
      END DO !IV 
C 
      RETURN 
      END 
