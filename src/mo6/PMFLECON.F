      REAL FUNCTION PMFLECON(MY,IV) 
C 
C  PMFLECON gets the gas and diesel fuel economies,  
C  from MOBILE6 or PART5. 
C 
C Changes: 
C 
C 13 Feb 01 Air TASK 11: New MOBILE6 function (from PART5). 
C 10 Nov 00 AIR Task 06: Updated code to MOBILE6 standards. 
C 
C  Called by DB_AGGR, DBEXHRS, PMGAS, PMDB_AGGR, 
C            PMDSO4SO2, PMGSO4SO2, PMLEADED, PMOUTDB 
C 
C    common blocks: 
C     /DIEFE / DIESFE,IYRFE,NIYRFE 
C     /RLCOM1/ ROADFE 
C     /RLCOM4/ ROADFEMC 
C     /VTNAMS/ VTMC,VTGASBUS 
C     /VVTYPS/ VVGASMC,VDM6PRT 
C 
C  Output on return: 
C     PMFLECON 
C 
C  Local variable dictionary: 
C 
C   Name   Type              Description 
C  ------  ----  ------------------------------------------------------------- 
C  IV       I    Vehicle class 
C  IVPART   I    MOBILE6 to PART5 vehicle class translation 
C  IYR      I    Model year pointer 
C  MDX      I    Model year index mapping to Model year(MY) for ROADFE and ROADFEMC 
C  MY       I    Model year 
C 
C 
      IMPLICIT NONE 
C 
      INCLUDE 'FRACT.I' 
      INCLUDE 'IVTYPE.I' 
      INCLUDE 'RLCOM1.I' 
      INCLUDE 'RLCOM4.I' 
C 
C  Declare parameter list. 
C 
      INTEGER IV 
      INTEGER MY 
C 
C  Declare local variables/arrays. 
C 
      INTEGER IVPART 
      INTEGER IYR 
      INTEGER MDX 
C 
C  Determine fuel economy for gas 
C 
      IF(VVGASMC(IV).EQ.1) THEN 
C 
C  For each Gas Vehicle Type select the corresponding MDX to MY  
C  for ROADFE/ROADFEMC 
C 
      IF(IV .EQ. VTMC) THEN 
        MDX=MY-1976+1 
        IF(MDX .LT.  1) MDX=1 
        IF(MDX .GT. 12) MDX=12 
        PMFLECON = ROADFEMC(MDX) 
C 
      ELSE 
        MDX=MY-1969+1 
        IF(MDX .LT.  1) MDX=1 
        IF(MDX .GT. 32) MDX=32 
          IF(IV.EQ.VTGASBUS) THEN 
            PMFLECON=ROADFE(MDX,14) 
          ELSE 
            PMFLECON=ROADFE(MDX,IV) 
          ENDIF 
      ENDIF 
C 
C  Determine fuel economy for diesel. Diesel values are from PART5. 
C 
      ELSE 
C 
C  Map the MOBILE6 vehicle classes into PART5 diesel space. 
C 
      IVPART = VDM6PRT(IV) - 5 
C 
      DO IYR=1,NIYRFE 
        IF (MY.LE.IYRFE(IYR)) THEN 
          MDX=IYR 
          GO TO 10 
        ENDIF 
      END DO 
C 
   10 PMFLECON = DIESFE(MDX,IVPART) 
C 
      END IF 
C 
      RETURN 
      END 
