      SUBROUTINE ATCOMPLEX(IV,MY,IOXY,CAT) 
C 
C  This subroutine computes the Unconsolidated Complex Model values 
C  for VOC, Acetaldehyde, Benzene, Formaldehyde, 1,3 Butadiene, and MTBE. 
C  (The values match those produced by the EPA Tier0_CM.XLS spreadsheet, 
C   including the bug found for the Acetaldehyde centering values.) 
C  The ratios of these toxics to the VOC Complex values are also computed. 
C 
      IMPLICIT NONE 
C 
      INCLUDE 'ATOX1.I' 
      INCLUDE 'ATOX2.I' 
      INCLUDE 'ATOX3.I'   ! VOCCOMPLEX 
      INCLUDE 'EVAL.I'    ! MEVAL 
      INCLUDE 'IVTYPE.I'  ! VVLDGAS 
C 
      INTEGER IASTM            ! The fuel ASTM Class, 1=Summer, 2=Winter 
      INTEGER IFP              ! Fuel property 
      INTEGER IOXY             ! Oxygenate Type, includes no oxygenate 
      INTEGER ITG              ! Technology Group (1-9 are Normal Emitters 10 is High Emitters) 
      INTEGER IV               ! Vehicle class 
      INTEGER JFP              ! Fuel property, mapped 
      INTEGER JTG              ! Technology Group, mapped 
      INTEGER MAP(7)           ! Fuel property mapping array 
      INTEGER MY               ! Model year 
      REAL    ACETCOMPLEX(10)  ! Acetaldehyde, mg/mi, from Complex Model 
      REAL    BASEEMIS1(10,14) ! Storage array for complex model and fuel property combinations, base fuel 
      REAL    BASEEMIS2(10)    ! Tech group emissions, in mg/mi, base fuel 
      REAL    BENZCOMPLEX(10)  ! Benzene, mg/mi, from Complex Model 
      REAL    BUTACOMPLEX(10)  ! 1,3 Butadiene, mg/mi, from Complex Model 
      REAL    CAT(3)           ! Catalyst fractions 
      REAL    CENTERMTBE       ! Centering value for MTBE     
      REAL    EFFECT(10,14)    ! Storage array for complex model coefficients 
      REAL    FORMCOMPLEX(10)  ! Formaldehyde, mg/mi, from Complex Model 
      REAL    MTBECOMPLEX(10)  ! MTBE, mg/mi, from Complex Model 
      REAL    OFFSET           ! MTBE offset 
      REAL    SUMBASE          ! Summation variable, for base fuel      
      REAL    SUMM             ! Summation for MTBECOMPLEX 
      REAL    SUMTRGT          ! Summation variable, for target fuel 
      REAL    SUMV             ! Summation for VOCCOMPLEX 
      REAL    TRGTEMIS1(10,14) ! Storage array for complex model and fuel property combinations, target fuel 
      REAL    TRGTEMIS2(10)    ! Tech group emissions, in mg/mi, target fuel 
C 
C  Initialize ratio arrays 
C 
      BENZRATIO=0. 
      BUTARATIO=0. 
      ACETRATIO=0. 
      FORMRATIO=0. 
      MTBERATIO=0. 
C 
C  Set up the oxygenate components for the Complex Model. All oxygen 
C  content values are in wt%. 
C 
      ATTRGT( 1)=ATOXY(1,IOXY)                  !OXYGENATE, wt% 
C 
      ATTRGT( 9)=0. 
      IF(IOXY.EQ.1) ATTRGT( 9)=ATOXY(1,IOXY)    !MTBE, wt% 
C 
      ATTRGT(10)=0. 
      IF(IOXY.EQ.2) ATTRGT(10)=ATOXY(1,IOXY)    !ETBE, wt% 
C 
      ATTRGT(11)=0. 
      IF(IOXY.EQ.3) ATTRGT(11)=ATOXY(1,IOXY)    !ETOH, wt% 
C 
      ATTRGT(12)=0. 
      IF(IOXY.EQ.4) ATTRGT(12)=ATOXY(1,IOXY)    !TAME, wt% 
C 
C  Set the ASTM Class based on the Evaluation Month 
C 
      IASTM = 1                !Summer C Class 
      IF(MEVAL.EQ.1) IASTM = 2 !Winter C Class 
C 
C  ********************** VOC COMPLEX ********************************* 
C 
C  Initialize the Effect array 
C 
      EFFECT=0. 
C 
C  Build the Effect array for Fuel Parameters: 
C  1:Oxygen, 2:Sulfur, 3:RVP, 4:E200 
C  5:E300, 6:Aromatic, 7:Olefins. 
C    
C   
      DO IFP=1,7 
        DO ITG=1,10 
          JTG=ITG+1 
          EFFECT(ITG,IFP)=ATVOC(IFP,1)+ATVOC(IFP,JTG) 
        END DO 
      END DO 
C 
C  Build the Effect array for Fuel Parameters: 
C  12:E200SQR, 13:E300SQR, 14:E300ARO. 
C    
C   
      DO IFP=12,14 
        JFP=IFP-4     ! Map IFP to obtain E200sq, E300sq or E300AROM term 
        DO ITG=1,10 
          JTG=ITG+1 
          EFFECT(ITG,IFP)=ATVOC(JFP,1)+ATVOC(JFP,JTG) 
        END DO 
      END DO 
C 
C  Initialize the Interim Base and Target arrays 
C 
      BASEEMIS1=0. 
      TRGTEMIS1=0. 
C 
C  Compute the Base and Target Emissions for: 
C  1:Oxygen, 2:Sulfur, 3:RVP, 4:E200 
C  5:E300, 6:Aromatic, 7:Olefins.  
C 
      DO IFP=1,7 
        DO ITG=1,10 
            BASEEMIS1(ITG,IFP)= 
     &         EFFECT(ITG,IFP)*(ATBASE(IASTM,IFP)-ATCENTER(IFP)) 
            TRGTEMIS1(ITG,IFP)= 
     &         EFFECT(ITG,IFP)*(ATTRGT(      IFP)-ATCENTER(IFP)) 
        END DO 
      END DO 
C 
C  Compute the Base and Target Emissions for: 
C  12:E200SQR, 13:E300SQR 
C 
      DO IFP=12,13 
        JFP=IFP-8     ! Map IFP to obtain E200 or E300 term 
        DO ITG=1,10 
            BASEEMIS1(ITG,IFP)= 
     &         EFFECT(ITG,IFP)*(ATBASE(IASTM,JFP)-ATCENTER(JFP))**2 
            TRGTEMIS1(ITG,IFP)= 
     &         EFFECT(ITG,IFP)*(ATTRGT(      JFP)-ATCENTER(JFP))**2 
        END DO 
      END DO 
C 
C  Compute the Base and Target Emissions for: 
C  14:E300AROM 
C 
      DO ITG=1,10 
          BASEEMIS1(ITG,14)= 
     &       EFFECT(ITG,14)*(ATBASE(IASTM,5)-ATCENTER(5)) !E300 term 
     &                     *(ATBASE(IASTM,6)-ATCENTER(6)) !Aromatic term 
          TRGTEMIS1(ITG,14)= 
     &       EFFECT(ITG,14)*(ATTRGT(      5)-ATCENTER(5)) !E300 term 
     &                     *(ATTRGT(      6)-ATCENTER(6)) !Aromatic term 
      END DO 
C 
C  Initialize the Final Base and Target arrays 
C 
      BASEEMIS2=0. 
      TRGTEMIS2=0. 
C 
C  Compute the Base and Target Emissions, in mg/mi, for the Tech Groups 
C 
      DO ITG=1,10 
C 
        SUMBASE=0. 
        SUMTRGT=0. 
C 
        DO IFP=1,7 
          SUMBASE=SUMBASE+BASEEMIS1(ITG,IFP) 
          SUMTRGT=SUMTRGT+TRGTEMIS1(ITG,IFP) 
        END DO 
C 
        DO IFP=12,14 
          SUMBASE=SUMBASE+BASEEMIS1(ITG,IFP) 
          SUMTRGT=SUMTRGT+TRGTEMIS1(ITG,IFP) 
        END DO 
C 
        BASEEMIS2(ITG)=EXP(SUMBASE) 
        TRGTEMIS2(ITG)=EXP(SUMTRGT) 
C 
      END DO 
C 
C  Compute the Complex Model VOC, in mg/mi 
C 
      VOCCOMPLEX=0. 
C 
      DO ITG=1,10 
        IF(BASEEMIS2(ITG).NE.0.) 
     &    VOCCOMPLEX(ITG)=ATVOCBASE(ITG)*TRGTEMIS2(ITG)/BASEEMIS2(ITG) 
      END DO 
C 
C  ********************** BENZENE COMPLEX ***************************** 
C 
C  Initialize the Effect array 
C 
      EFFECT=0. 
C 
C  Build the Effect array for Fuel Parameters: 
C  1:Sulfur, 2:E200, 3:E300, 4:Aromatic, 5:Benzene, 6:Oxygen. 
C   
      DO IFP=1,6 
        DO ITG=1,10 
          JTG=ITG+1 
          EFFECT(ITG,IFP)=ATBENZ(IFP,1)+ATBENZ(IFP,JTG) 
        END DO 
      END DO 
C 
C  Initialize the Interim Base and Target arrays 
C 
      BASEEMIS1=0. 
      TRGTEMIS1=0. 
C 
C  Map and Compute the Base and Target Emissions for: 
C  1:Sulfur, 2:E200, 3:E300, 4:Aromatic, 5:Benzene, 6:Oxygen. 
C 
      MAP(1)=2 !Sulfur 
      MAP(2)=4 !E200 
      MAP(3)=5 !E300 
      MAP(4)=6 !Aromatic 
      MAP(5)=8 !Benzene 
      MAP(6)=1 !Oxygen 
C 
      DO IFP=1,6 
        JFP=MAP(IFP) 
        DO ITG=1,10 
            BASEEMIS1(ITG,IFP)= 
     &         EFFECT(ITG,IFP)*(ATBASE(IASTM,JFP)-ATCENTER(JFP)) 
            TRGTEMIS1(ITG,IFP)= 
     &         EFFECT(ITG,IFP)*(ATTRGT(      JFP)-ATCENTER(JFP)) 
        END DO 
      END DO 
C 
C  Initialize the Final Base and Target arrays 
C 
      BASEEMIS2=0. 
      TRGTEMIS2=0. 
C 
C  Compute the Base and Target Emissions, in mg/mi, for the Tech Groups 
C 
      DO ITG=1,10 
C 
        SUMBASE=0. 
        SUMTRGT=0. 
C 
        DO IFP=1,6 
          SUMBASE=SUMBASE+BASEEMIS1(ITG,IFP) 
          SUMTRGT=SUMTRGT+TRGTEMIS1(ITG,IFP) 
        END DO 
C 
        BASEEMIS2(ITG)=EXP(SUMBASE) 
        TRGTEMIS2(ITG)=EXP(SUMTRGT) 
C 
      END DO 
C 
C  Compute the Complex Model Benzene, in mg/mi, and the Benzene Ratios. 
C  If model year is 1982+, apply the normal and high emitter Off-cycle 
C  Emissions effects. 
C 
      BENZCOMPLEX=0. 
C 
      DO ITG=1,10 
C 
        IF(BASEEMIS2(ITG).NE.0.) 
     &    BENZCOMPLEX(ITG)=ATBENZBASE(ITG)*TRGTEMIS2(ITG)/BASEEMIS2(ITG) 
C 
        IF(VOCCOMPLEX(ITG).NE.0.) 
     &    BENZRATIO(ITG)=BENZCOMPLEX(ITG)/VOCCOMPLEX(ITG) 
C 
        IF(MY.GE.1982) THEN 
          IF(ITG.LE.9) THEN 
            BENZRATIO(ITG)=BENZRATIO(ITG)*1.315  !Normals 
          ELSE 
            BENZRATIO(ITG)=BENZRATIO(ITG)*1.126  !Highs 
          END IF 
        END IF 
C 
      END DO 
C 
C  ********************** BUTADIENE COMPLEX *************************** 
C 
C  Initialize the Effect array 
C 
      EFFECT=0. 
C 
C  Build the Effect array for Fuel Parameters: 
C  1:Oxygen, 2:Sulfur, 3:E200, 4:E300, 5=Aromatic, 6:Olefins. 
C   
      DO IFP=1,6 
        DO ITG=1,10 
          JTG=ITG+1 
          EFFECT(ITG,IFP)=ATBUTA(IFP,1)+ATBUTA(IFP,JTG) 
        END DO 
      END DO 
C 
C  Initialize the Interim Base and Target arrays 
C 
      BASEEMIS1=0. 
      TRGTEMIS1=0. 
C 
C  Map and Compute the Base and Target Emissions for: 
C  1:Oxygen, 2:Sulfur, 3:E200, 4:E300, 5=Aromatic, 6:Olefins. 
C 
      MAP(1)=1 !Oxygen 
      MAP(2)=2 !Sulfur 
      MAP(3)=4 !E200 
      MAP(4)=5 !E300 
      MAP(5)=6 !Aromatic 
      MAP(6)=7 !Olefins 
C 
      DO IFP=1,6 
        JFP=MAP(IFP) 
        DO ITG=1,10 
            BASEEMIS1(ITG,IFP)= 
     &         EFFECT(ITG,IFP)*(ATBASE(IASTM,JFP)-ATCENTER(JFP)) 
            TRGTEMIS1(ITG,IFP)= 
     &         EFFECT(ITG,IFP)*(ATTRGT(      JFP)-ATCENTER(JFP)) 
        END DO 
      END DO 
C 
C  Initialize the Final Base and Target arrays 
C 
      BASEEMIS2=0. 
      TRGTEMIS2=0. 
C 
C  Compute the Base and Target Emissions, in mg/mi, for the Tech Groups 
C 
      DO ITG=1,10 
C 
        SUMBASE=0. 
        SUMTRGT=0. 
C 
        DO IFP=1,6 
          SUMBASE=SUMBASE+BASEEMIS1(ITG,IFP) 
          SUMTRGT=SUMTRGT+TRGTEMIS1(ITG,IFP) 
        END DO 
C 
        BASEEMIS2(ITG)=EXP(SUMBASE) 
        TRGTEMIS2(ITG)=EXP(SUMTRGT) 
C 
      END DO 
C 
C  Compute the Complex Model Butadiene, in mg/mi, and the Butdiene Ratios. 
C  If model year is 1982+, apply the normal and high emitter Off-cycle 
C  Emissions effects. 
C 
      BUTACOMPLEX=0. 
C 
      DO ITG=1,10 
C 
        IF(BASEEMIS2(ITG).NE.0.) 
     &    BUTACOMPLEX(ITG)=ATBUTABASE(ITG)*TRGTEMIS2(ITG)/BASEEMIS2(ITG) 
C 
        IF(VOCCOMPLEX(ITG).NE.0.) 
     &    BUTARATIO(ITG)=BUTACOMPLEX(ITG)/VOCCOMPLEX(ITG) 
C 
        IF(MY.GE.1982) THEN 
          IF(ITG.LE.9) THEN 
            BUTARATIO(ITG)=BUTARATIO(ITG)*1.037  !Normals 
          ELSE 
            BUTARATIO(ITG)=BUTARATIO(ITG)*0.708  !Highs 
          END IF 
        END IF 
C 
      END DO 
C 
C  ********************** ACETALDEHYDE COMPLEX ************************ 
C 
C  Initialize the Effect array 
C 
      EFFECT=0. 
C 
C  Build the Effect array for Fuel Parameters: 
C  1:Sulfur, 2:RVP 3:E300, 4:Aromatic, 5:MTBE, 6:ETBE, 7:Ethanol. 
C   
      DO IFP=1,7 
        DO ITG=1,10 
          JTG=ITG+1 
          EFFECT(ITG,IFP)=ATACET(IFP,1)+ATACET(IFP,JTG) 
        END DO 
      END DO 
C 
C  Initialize the Interim Base and Target arrays 
C 
      BASEEMIS1=0. 
      TRGTEMIS1=0. 
C 
C  Map and Compute the Base and Target Emissions for: 
C  1:Sulfur, 2:RVP 3:E300, 4:Aromatic, 5:MTBE, 6:ETBE, 7:Ethanol. 
C  Programming Note: special use of acetaldehyde centering values. 
C 
      MAP(1)=2  !Sulfur 
      MAP(2)=3  !RVP 
      MAP(3)=5  !E300 
      MAP(4)=6  !Aromatic 
      MAP(5)=9  !MTBE 
      MAP(6)=10 !ETBE 
      MAP(7)=11 !ETOH 
C 
      DO IFP=1,7 
        JFP=MAP(IFP) 
        DO ITG=1,10 
            BASEEMIS1(ITG,IFP)= 
     &         EFFECT(ITG,IFP)*(ATBASE(IASTM,JFP)-ATCENTERA(JFP)) 
            TRGTEMIS1(ITG,IFP)= 
     &         EFFECT(ITG,IFP)*(ATTRGT(      JFP)-ATCENTERA(JFP)) 
        END DO 
      END DO 
C 
C  Initialize the Final Base and Target arrays 
C 
      BASEEMIS2=0. 
      TRGTEMIS2=0. 
C 
C  Compute the Base and Target Emissions, in mg/mi, for the Tech Groups 
C 
      DO ITG=1,10 
C 
        SUMBASE=0. 
        SUMTRGT=0. 
C 
        DO IFP=1,7 
          SUMBASE=SUMBASE+BASEEMIS1(ITG,IFP) 
          SUMTRGT=SUMTRGT+TRGTEMIS1(ITG,IFP) 
        END DO 
C 
        BASEEMIS2(ITG)=EXP(SUMBASE) 
        TRGTEMIS2(ITG)=EXP(SUMTRGT) 
C 
      END DO 
C 
C  Compute the Complex Model Acetaldehyde, in mg/mi, and the Acetaldehyde Ratios. 
C  If model year is 1982+, apply the normal and high emitter Off-cycle 
C  Emissions effects. 
C 
      ACETCOMPLEX=0. 
C 
      DO ITG=1,10 
C 
        IF(BASEEMIS2(ITG).NE.0.) 
     &    ACETCOMPLEX(ITG)=ATACETBASE(ITG)*TRGTEMIS2(ITG)/BASEEMIS2(ITG) 
C 
        IF(VOCCOMPLEX(ITG).NE.0.) 
     &    ACETRATIO(ITG)=ACETCOMPLEX(ITG)/VOCCOMPLEX(ITG) 
C 
        IF(MY.GE.1982) THEN 
          IF(ITG.LE.9) THEN 
            ACETRATIO(ITG)=ACETRATIO(ITG)*1.020  !Normals 
          ELSE 
            ACETRATIO(ITG)=ACETRATIO(ITG)*0.919  !Highs 
          END IF 
        END IF 
C 
      END DO 
C 
C  ********************** FORMALDEHYDE COMPLEX ************************* 
C 
C  Initialize the Effect array 
C 
      EFFECT=0. 
C 
C  Build the Effect array for Fuel Parameters: 
C  1:Sulfur, 2:RVP 3:E300, 4:Aromatic, 5:Olefins, 6:MTBE, 7:Ethanol. 
C   
      DO IFP=1,7 
        DO ITG=1,10 
          JTG=ITG+1 
          EFFECT(ITG,IFP)=ATFORM(IFP,1)+ATFORM(IFP,JTG) 
        END DO 
      END DO 
C 
C  Initialize the Interim Base and Target arrays 
C 
      BASEEMIS1=0. 
      TRGTEMIS1=0. 
C 
C  Map and Compute the Base and Target Emissions for: 
C  1:Sulfur, 2:RVP 3:E300, 4:Aromatic, 5:Olefins, 6:MTBE, 7:Ethanol. 
C 
      MAP(1)=2  !Sulfur 
      MAP(2)=3  !RVP 
      MAP(3)=5  !E300 
      MAP(4)=6  !Aromatic 
      MAP(5)=7  !Olefins 
      MAP(6)=9  !MTBE 
      MAP(7)=11 !ETOH 
C 
      DO IFP=1,7 
        JFP=MAP(IFP) 
        DO ITG=1,10 
            BASEEMIS1(ITG,IFP)= 
     &         EFFECT(ITG,IFP)*(ATBASE(IASTM,JFP)-ATCENTER(JFP)) 
            TRGTEMIS1(ITG,IFP)= 
     &         EFFECT(ITG,IFP)*(ATTRGT(      JFP)-ATCENTER(JFP)) 
        END DO 
      END DO 
C 
C  Initialize the Final Base and Target arrays 
C 
      BASEEMIS2=0. 
      TRGTEMIS2=0. 
C 
C  Compute the Base and Target Emissions, in mg/mi, for the Tech Groups 
C 
      DO ITG=1,10 
C 
        SUMBASE=0. 
        SUMTRGT=0. 
C 
        DO IFP=1,7 
          SUMBASE=SUMBASE+BASEEMIS1(ITG,IFP) 
          SUMTRGT=SUMTRGT+TRGTEMIS1(ITG,IFP) 
        END DO 
C 
        BASEEMIS2(ITG)=EXP(SUMBASE) 
        TRGTEMIS2(ITG)=EXP(SUMTRGT) 
C 
      END DO 
C 
C  Compute the Complex Model Formaldehyde, in mg/mi, and the Formaldehyde Ratios. 
C  If model year is 1982+, apply the normal and high emitter Off-cycle 
C  Emissions effects. 
C 
      FORMCOMPLEX=0. 
C 
      DO ITG=1,10 
C 
        IF(BASEEMIS2(ITG).NE.0.) 
     &    FORMCOMPLEX(ITG)=ATFORMBASE(ITG)*TRGTEMIS2(ITG)/BASEEMIS2(ITG) 
C 
        IF(VOCCOMPLEX(ITG).NE.0.) 
     &    FORMRATIO(ITG)=FORMCOMPLEX(ITG)/VOCCOMPLEX(ITG) 
C 
        IF(MY.GE.1982) THEN 
          IF(ITG.LE.9) THEN 
            FORMRATIO(ITG)=FORMRATIO(ITG)*1.163  !Normals 
          ELSE 
            FORMRATIO(ITG)=FORMRATIO(ITG)*0.894  !Highs 
          END IF 
        END IF 
C 
      END DO 
C 
C  ********************** MTBE COMPLEX ******************************** 
C 
C  Compute the CenterMTBE term and assign the MTBE equation offset 
C  (The 11./2. is the Complex Models coversion factor from wt% to vol%) 
C 
      CENTERMTBE= ATTRGT(9)*(11./2.) - ATCENTER(14) 
C 
      OFFSET=-0.01 
C 
C  Compute the MTBE, in mg/mi, for Tech Group 1 
C 
      MTBECOMPLEX(1)=EXP(ATMTBE(1,1) 
     &                  +ATMTBE(2,1)*CENTERMTBE 
     &                  +ATMTBE(3,1)*CENTERMTBE**2 
     &                  +ATMTBE(4,1)*CENTERMTBE**2) 
     &              +OFFSET 
C 
C  Compute the MTBE, in mg/mi, for Tech Groups 2-9 
C 
      DO ITG=2,9 
        MTBECOMPLEX(ITG)=EXP(ATMTBE(1,ITG) 
     &                      +ATMTBE(2,ITG)*CENTERMTBE 
     &                      +ATMTBE(3,ITG)*CENTERMTBE**2) 
     &                  +OFFSET 
      END DO 
C 
C  Compute the MTBE, in mg/mi, for Tech Group 10 
C 
      MTBECOMPLEX(10)=EXP(ATMTBE(1,10) 
     &                   +ATMTBE(2,10)*CENTERMTBE 
     &                   +ATMTBE(3,10)*CENTERMTBE**2 
     &                   +ATMTBE(4,10)*CENTERMTBE) 
     &               +OFFSET 
CDFK 
CDFK      WRITE(99,1) IV,MY,CENTERMTBE,MTBECOMPLEX,VOCCOMPLEX 
CDFK   1  FORMAT('MTBE1',I3,I5,11F10.6,/23X,10F10.6) 
C 
C  Following the Complex Model method EXACTLY, convert values to VOC rates. 
C  (Be sure that the base MTBECOMPLEX values are not negative.) 
C  Determine special values for Tech Groups 5, 8 and 9 via weighting 
C  of Tech Groups 1-4 and 6-7. 
C 
      SUMM=0. 
      SUMV=0. 
C 
      DO ITG=1,10 
C 
      MTBECOMPLEX(ITG)=MAX(0.,MTBECOMPLEX(ITG)) 
C 
      IF(CMHC(ITG).GT.0.)  
     * MTBECOMPLEX(ITG)=MTBECOMPLEX(ITG) 
     *                 *VOCCOMPLEX(ITG) 
     *                 /(CMHC(ITG)*1000.)  ! Convert grams to milligrams 
C 
      IF(ITG.LE.4. OR. ITG.EQ.6 .OR. ITG.EQ.7) THEN 
        SUMM=SUMM+MTBECOMPLEX(ITG) 
        SUMV=SUMV+VOCCOMPLEX(ITG) 
      END IF 
C 
      END DO 
C 
      MTBECOMPLEX(5)=VOCCOMPLEX(5)*SUMM/SUMV 
CDFK      MTBECOMPLEX(8)=0. ! Per Complex Model 
      MTBECOMPLEX(8)=MTBECOMPLEX(2) ! Per EPA assumption 
      MTBECOMPLEX(9)=VOCCOMPLEX(9)*SUMM/SUMV 
CDFK 
CDFK      WRITE(99,2) IV,MY,CENTERMTBE,MTBECOMPLEX,VOCCOMPLEX 
CDFK   2  FORMAT('MTBE2',I3,I5,11F10.6,/23X,10F10.6) 
C 
C  Compute the MTBE Ratios. If model year is 1982+, apply the normal 
C  and high emitter Off-cycle Emissions effects. 
C 
      DO ITG=1,10 
C 
        IF(VOCCOMPLEX(ITG).NE.0.) 
     &    MTBERATIO(ITG)=MTBECOMPLEX(ITG)/VOCCOMPLEX(ITG) 
C 
        IF(MY.GE.1982) THEN 
          IF(ITG.LE.9) THEN 
            MTBERATIO(ITG)=MTBERATIO(ITG)*0.825  !Normals 
          ELSE 
            MTBERATIO(ITG)=MTBERATIO(ITG)*0.965  !Highs 
          END IF 
        END IF 
C 
      END DO 
C 
C  ********************** TECH GROUP WEIGHTING ******************************** 
C 
C  Weight the 10 AT Tech Group Conversions to form a single Normal and High value. 
C  Then weight the single normal and high rates by the oxygenate market share so that 
C  the two values are compatible with the MOBILE6 VOC and TOG values. 
C   
      CALL ATNH(IV,MY,ATOXY(2,IOXY),CAT) 
C 
      RETURN 
      END 
