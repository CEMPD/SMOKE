      BLOCK DATA BD48 
C 
C  BLOCK DATA Subprogram 48: /ATOX1/,/ATOX3/,/ATOX4/,ATOX5/ 
C 
C  BLOCK DATA 43 initializes the Air Toxics commons. 
C 
C  Changes: (Last change first) 
C 
C  08 Mar 02 AIR Task 18: Moved ADDITIONAL HAPS related arrays to ATWORK. 
C            Updated assocated initializations and comments. 
C 28 Sep 01: AIR TASK 25: Changed ATSSFLAG to SSFLAG to accomodate PM output. 
C  24 Sep 01 AIR Task 24: Comments 
C  13 Sep 01 AIR Task 18: Added HAPS, and removed initialization 
C            of large arrays. 
C  06 Sep 01 AIR Task 15: New Block Data. 
C 
C  Common block array subscripts: 
C 
C  ATOXY(3,5) -  ATOXY(I,IEVAP) 
C  TOXCHK(6)  -  TOXCHK(ICHK) 
C 
C  ACETRATIO(10)  - ACETRATIO (ITG) 
C  BENZRATIO(10)  - BENZRATIO (ITG) 
C  BUTARATIO(10)  - BUTARATIO (ITG) 
C  FORMRATIO(10)  - FORMRATIO (ITG) 
C  GROUP(10)      - GROUP     (ITG) 
C  MTBERATIO(10)  - MTBERATIO (ITG) 
C  VOCCOMPLEX(10) - VOCCOMPLEX(ITG) 
C 
C  ATCONV(2,2,6)         - ATCONV(INH,IVT,IAT) 
C  ATEMFAC(2,10,5,3,2)   - ATEMFAC(ISR,ISTD,IFDS,INH,IVT) 
C  ATHEF(24,5,2)         - ATHEF(IH,IROAD,IVT)  
C  ATHEV(24,5,5)         - ATHEV(IH,IROAD,IEVAP)  
C  ATNHWT(2)             - ATNHWT(INH) 
C 
C  Common block dictionary: 
C 
C   Name    Type                  Description 
C  ------   ----  ----------------------------------------------------- 
C  /ATOX1/ 
C  AROM       R    percentage of aromatic compounds in gasoline on a Volume Basis (required input) 
C  BENZ       R    percentage of benzene in gasoline on a Volume Basis (required input) 
C  E200       R    percentage of vapor of gasoline fuel at 200F (required input) 
C  E300       R    percentage of vapor of gasoline fuel at 300F (required input) 
C  OLEF       R    percentage of olefin in gasoline on a Volume Basis (required input) 
C  TOXCHK     I    checklist of required air toxics run inputs (1=not input, 2=input) for: 
C                  1: E200, 2: E300, 3: Aromatics, 4: Olefins, 5: Benzene, 6: Oxygenate 
C  TOXFLG     L    switch to run the air toxics subroutine 
C  ATOXY      R     Oxygenate content, wt%, for 1:MTBE, 2:ETBE, 3:ETHANOL, 4:TAME, 5:NoOxy 
C  MAXIAT     I    Number of AT pollutant types 
C  SSFLAG     L    switch to print out the air toxics and/or PM spreadsheet 
C 
C  /ATOX3/ 
C  ACETRATIO  R    Ratio of Air Toxic Acetaldehyde emissions to MOBILE6 VOC/TOG emissions 
C  BENZRATIO  R    Ratio of Air Toxic Benzene emissions to MOBILE6 VOC/TOG emissions 
C  BUTARATIO  R    Ratio of Air Toxic 1,3Butadiene emissions to MOBILE6 VOC/TOG emissions 
C  FORMRATIO  R    Ratio of Air Toxic Formaldehyde emissions to MOBILE6 VOC/TOG emissions 
C  GROUP      R    Technology group fractions  
C  MTBERATIO  R    Ratio of Air Toxic MTBE emissions to MOBILE6 VOC/TOG emissions 
C  VOCCOMPLEX R    VOC emissions from the Complex Model 
C   
C  All /ATOX3/ Common Block arrays have the following ITG index meaning: 
C                    1: PFI , 3way   , No Air, EGR 
C                    2: PFI , 3way   , No Air, No EGR 
C                    3: TBI , 3way   , No Air, EGR 
C                    4: PFI , 3way+Ox, Air   , EGR 
C                    5: PFI , 3way   , Air   , EGR 
C                    6: TBI , 3way   , Air   , EGR 
C                    7: TBI , 3way+Ox, Air   , EGR 
C                    8: TBI , 3way   , No Air, No EGR 
C                    9: CARB, 3way+Ox, Air   , EGR 
C                   10: All High Emitters 
C 
C  /ATOX4/ 
C  ATCONV    R     Final ratios of Air Toxic to MOBILE6 VOC/TOG emissions  
C  ATEMFAC   R     The basic MOBILE6 exhaust VOC/TOG emission factors corrected for tampering and methane 
C  ATHEF     R     The final, corrected MOBILE6 exhaust VOC/TOG emission factors 
C  ATHEV     R     The final, corrected MOBILE6 evap emission factors 
C  ATNHWT    R     Air Toxic normal and high weighting factors 
C 
C  /ATOX5/ 
C  ATHAPS    L     Air Toxic HAPS toggle flag 
C  NATH      I     Number of valid Air Toxic Additional HAPS entries 
C  NMAP      I     Mapping index for the pollutant number       
C 
C 
      IMPLICIT NONE 
C 
      INCLUDE 'ATOX1.I' 
      INCLUDE 'ATOX3.I' 
      INCLUDE 'ATOX4.I' 
      INCLUDE 'ATOX5.I' 
C 
C  /ATOX1/ 
C 
      DATA AROM    /0./ 
      DATA ATOXY   /15*0./ 
      DATA BENZ    /0./ 
      DATA E200    /0./ 
      DATA E300    /0./ 
      DATA MAXIAT  /6/ 
      DATA OLEF    /0./ 
      DATA SSFLAG/.FALSE./ 
      DATA TOXCHK  /6*1/ 
      DATA TOXFLG  /.FALSE./ 
C 
C  /ATOX3/ 
C 
      DATA ACETRATIO  /10*0./ 
      DATA BENZRATIO  /10*0./ 
      DATA BUTARATIO  /10*0./ 
      DATA FORMRATIO  /10*0./ 
      DATA GROUP      /10*0./ 
      DATA MTBERATIO  /10*0./ 
      DATA VOCCOMPLEX /10*0./ 
C 
C  /ATOX4/ 
C 
      DATA ATCONV  /24*0./ 
      DATA ATEMFAC /240*0./ 
      DATA ATHEF   /240*0./ 
      DATA ATHEV   /600*0./ 
      DATA ATNHWT  /2*0./ 
C 
C  /ATOX5/ 
C 
      DATA ATHAPS    /.FALSE./ 
      DATA NATH      /0/ 
      DATA NMAP      /0/ 
C 
      END 
       
