      REAL FUNCTION HCOXY(MYR,IV,BEF) 
C 
C  HCOXY computes the oxy fuel. Due to the introduction of oxidation and 
C  3-way catalyst on LDTs, special handling of those technologies was needed. 
C 
C  Called by HRLOOP. 
C 
C  Calls P3WPOX and NMOGFID 
C 
C  Changes: (Last change first) 
C 
C  04 Dec 00 @EPA-bag Task X6, removed MY, JDX, and IDX from MYCODE.I; removed 
C            MYCODE.I; model year now passed through parameter lists 
C  16 May 00 @EPA-ddj BUG27 Modified calculation from separate run/start 
C            to FTP. 
C  23 Mar 00 @DynTel-rwk 1-028, Put code back in to calculate oxy effects 
C            if area is RFG. 
C  25 Feb 00 @DynTel-ddj 1-028, Removed PHSFAC array. Added RETURN 
C            if diesel. 
C  18 Feb 00 @DynTel-BAG 1-028, Removed RFG adjustment.  Adjustment now 
C            performed in HRLOOP via call to RFG_CF. 
C  04 Apr 00 AIR Task 02: Removed ASTM class code, replaced with a hardcoded 
C            RFG Region 2 (ASTM class C). This code needs to be fixed when 
C            the RFG Task in installed. 
C  29 Mar 00 @EPA-eg  HRLOOP revisions. 
C  28 Sep 99 @DynTel-ddj for EPA  Added ISTD to HCOXY parameter list 
C            and call to OFFMTH. 
C  17 Aug 99 @DynTel-HXQ 2-667 HC oxygen effect  
C  16 Aug 99 @DynTel-HXQ 2-667 modified this routine from FUEL.FOR for HC  
C 
C  Input on call: 
C 
C    argument list: MYR,IV,BEFADJ 
C    common blocks: 
C     /FLAGS1/ OXYFLG 
C     /OFFSET/ OFFMTH 
C     /OXY1/   SHRMKT,OXYCNT 
C     /OXY2/   MXOXMY,MAXOXY,OXYVAL,OXYBF 
C     /RFORM1/ RFGON 
C     /VVNAMS/ VVDSL, VVHDGV, VVLDGT12, VVLDGT34  (IVTYPE.I) 
C     /VVTYPS/ VTLDGV, VTMC  (IVTYPE.I) 
C 
C  Output on return: 
C 
C    function: FUEL 
C    common blocks: 
C 
C 
C  Local array subscripts: 
C 
C 
C  Local variable dictionary: 
C 
C   Name   Type              Description 
C  ------  ----  ------------------------------------------------------- 
C  AVGFAC   R    percent decrease in emissions from industry average fuel to 
C                Indolene fuel 
C  BENFAC   R    benefit factor for determining oxy adjustment 
C  CAT      R    fraction of vehicles that are catalyst equipped 
C  FUEL1    R    1st oxy fuel adjusment factor: pre-81 LDV, HDGV, MC and 
C                pre-91 LDTs 
C  FUEL2    R    2nd oxy fuel adjustment factor: 1981+ LDV and 86+ LDTs 
C  NOCAT    R    fraction of vehicles that are not catalyst equipped 
C  OFMTH    R    value holder for the methane offset to basic emission rate 
C  OXYCT    R    local variable for the oxygen content 
C  P3W      R    percent of vehicles with 3-way catalyst 
C  POX      R    percent of vehicles with ox catalyst 
C  TOG1     R    value holder for the TOG FID correction factor 
C  TOG2     R    value holder for the TOG FID oxygen correction factor 
C  WT1      R    1st fuel adjustment factor weighting factor 
C  WT2      R    2nd fuel adjustment factor weighting factor 
C  X1       R    lower x-value (IPOXY) for interpolation 
C  X2       R    upper x-value (IPOXY) for interpolation 
C  Y1       R    lower y-value (benefit factor) for interpolation 
C  Y2       R    upper y-value (benefit factor) for interpolation 
C 
C 
C 
      IMPLICIT NONE 
      INCLUDE 'FLAGS1.I' 
      INCLUDE 'IVTYPE.I' 
      INCLUDE 'OFFSET.I' 
      INCLUDE 'OXY1.I' 
      INCLUDE 'OXY2.I' 
      INCLUDE 'RFORM1.I' 
      INCLUDE 'RFORM2.I' 
      INCLUDE 'SOAK.I' 
C 
      REAL, DIMENSION(2), INTENT(IN) :: BEF 
C 
      INTEGER, INTENT(IN)     :: IV 
      INTEGER, INTENT(IN)     :: MYR 
C 
      REAL BEFFTP 
      REAL NOCAT 
      REAL CAT 
      REAL P3W 
      REAL POX 
      REAL WT1 
      REAL FUEL1 
      REAL FUEL2 
      REAL BENFAC 
      REAL Y1 
      REAL Y2 
      REAL WT2 
      REAL X1 
      REAL X2 
C 
      REAL, DIMENSION(2)       :: OXYCT 
C 
      INTEGER IFUEL 
      INTEGER MYMAP 
      INTEGER IGOXY 
      INTEGER IG 
      INTEGER IPOXY 
C 
C  Sep-22-1994 @ CSC-bsg Request 426 Updated PHSFAC. 
C  Continuation code is IV 
C  22-Jul-1996 @ DynTel-wcs Request 2-627 Alter PHSFAC values in order 
C                meet VOC HC Exhaust target values for LDGV & LDGT 
C 
C     Initialize. 
C 
      HCOXY = 1.0 
C 
C  If vehicle is diesel fueled, no adjustment is necessary. 
C 
      IF( VVDSL(IV).NE.1 ) THEN 
C 
C  Return if no reformulated or oxygenated fuel effects are to be calculated. 
C 
        IF(RFGON .OR. OXYFLG.NE.1) THEN 
C 
C  Initialize 
C 
          BENFAC = 0.0 
          BEFFTP = BEF(1) + BEF(2)*(0.43 + 0.57*0.160)/7.5 
C 
          IF(BEFFTP .GT. 0.000001) THEN 
C  Local variable for oxygen content cannot be greater than 2.7% for HC. 
C 
C  Initialize oxy variables. 
C 
            FUEL1 = 1.0 
            FUEL2 = 1.0 
            WT1   = 1.0 
            WT2   = 0.0 
            DO IFUEL = 2, 3 
              IF(OXYCNT(IFUEL-1).GT..027) THEN 
                 OXYCT(IFUEL-1)=.027 
              ELSE 
                 OXYCT(IFUEL-1) = OXYCNT(IFUEL-1) 
              END IF 
            END DO 
C 
C  Sep-22-1994 @ CSC-bsg Request 426  Modularized the calculation of the 
C  3-Way oxidation values, by calling a function to calculate them. 
C  Compute the number of catalyst and non-catalyst vehicles, 
C 
            CALL P3WPOX(MYR,IV,P3W,POX) 
            CAT   = P3W + POX 
            NOCAT = 1.0 - CAT 
C 
C  Compute the benefit factor and the 1st oxy fuel adjustment factor, 
C  as necessary. 
C 
            IF((IV.EQ.VTLDGV .AND. MYR.LT.1981).OR. 
     &        ((VVLDGT12(IV).EQ.1 .OR. VVLDGT34(IV).EQ.1) .AND. 
     &        MYR.LE.1990).OR. (VVHDGV(IV).EQ.1 .OR. IV.EQ.VTMC)) THEN 
              BENFAC = NOCAT*OXYFAC(1,1) + CAT*OXYFAC(1,2) 
              FUEL1 = SHRMKT(1) 
     &            + SHRMKT(2)*(1.0-OXYCT(1)*BENFAC) 
     &            + SHRMKT(3)*(1.0-OXYCT(2)*BENFAC) 
            ENDIF 
C 
C  Compute the 2nd oxy factor, as necessary. 
C 
            IF((IV.NE.VTLDGV .OR. MYR.GE.1981) .AND. 
     &        ((VVLDGT12(IV).NE.1 .AND. VVLDGT34(IV).NE.1) .OR. 
     &        MYR.GE.1986) .AND. 
     &        (VVHDGV(IV).NE.1 .AND. IV.NE.VTMC)) THEN 
C 
C  For 1986-87 LDTs, map the model year. 
C 
              MYMAP = MYR 
              IF((VVLDGT12(IV).EQ.1.OR.VVLDGT34(IV).EQ.1).AND. 
     &          (MYR.EQ.1986 .OR. MYR.EQ.1987)) MYMAP = MYMAP + 1 
C 
C  Look up the bracketing values and set up for interpolation. 
C 
              IGOXY = MYMAP - 1981 + 1 
              IF(IGOXY.GT.MXOXMY) IGOXY = MXOXMY 
C 
C Compute 1st HC fuel adjustment factor  
C 
              DO IPOXY = 1, 2 
                IG = 0 
                befadj: DO 
                  IG = IG + 1 
                  IF(BEFFTP.LE.OXYVAL(IG,IGOXY,IPOXY) .OR. 
     &              IG.EQ.MAXOXY) EXIT befadj 
                END DO befadj 
C 
                IF(IG.EQ.1) THEN 
                  BENFAC=OXYBF(IG,IGOXY,IPOXY) 
                ELSE IF (IG.GE.MAXOXY) THEN 
                  BENFAC = OXYBF(MAXOXY,IGOXY,IPOXY) 
                ELSE 
                  X1 = OXYVAL(IG-1,IGOXY,IPOXY) 
                  Y1 = OXYBF(IG-1,IGOXY,IPOXY) 
                  X2 = OXYVAL(IG,IGOXY,IPOXY) 
                  Y2 = OXYBF(IG,IGOXY,IPOXY) 
C 
C  Interpolate the benefit factor. 
C 
                  BENFAC = (BEFFTP*(Y2-Y1) + X2*Y1 - X1*Y2)/(X2-X1) 
                END IF 
                IF(IPOXY.EQ.1) FUEL2 = SHRMKT(2)*(1.0-OXYCT(1)*BENFAC) 
C 
                IF(IPOXY.EQ.2) THEN 
                  FUEL2 = SHRMKT(1) 
     &                  + FUEL2 
     &                  + SHRMKT(3)*(1.0-OXYCT(2)*BENFAC) 
                ENDIF 
              END DO     ! IPOXY 
C 
C Determine the 1986-90 LDT oxy technology weighting factors. 
C 
              IF((VVLDGT12(IV).EQ.1.OR.VVLDGT34(IV).EQ.1).AND. 
     &          (MYR.GE.1986 .AND. MYR.LE.1990)) THEN 
                WT1 = 1.0 - P3W 
                WT2 = P3W 
              ELSE 
                WT1 = 0.0 
                WT2 = 1.0 
              END IF 
            END IF 
C 
C  Compute the weighted oxy fuel adjustment factor. 
C 
            HCOXY = WT1*FUEL1 + WT2*FUEL2 
C 
          END IF 
        END IF 
      END IF 
C 
      RETURN 
      END 
