      SUBROUTINE PMCALX(ICY,INERR,RC) 
C 
C  Directs to calculation of EMFAC-PM, RCP, SOF, SO4, NH3, SO2. 
C 
C  Changes : 
C 
C  29 Aug 01 AIR Task 17: Added call to FUELPAR to be sure the in-use 
C            fuel parameters are current. 
C  24 Aug 01 AIR Task 17: Removed duplicate INERR update. 
C  23 May 01 AIR Task 17: Expanded PRTCHK to 3 (ZML/DR filenames). 
C  10 May 01 AIR Task 11/5: Added PM Idle. 
C  17 Apr 01 AIR Task 11: Lahey FORTRANv4.0-5.5 have several BLOCK DATA and 
C            MODULE bugs. To circumvent these problems, the data arrays 
C            for Speed Correction Factors, Zero Mile Levels and Deterioration 
C            Rates (for PM and NH3) are now read in. 
C  16 Feb 01 AIR Task 11: New Module 
C 
C  Called by: EFCALX 
C 
C  Calls: PMBRAKETIRE, PMDRCPSOF, PMDSO4SO2, PMEMFAC PMGASPM, 
C         PMGSO4SO2, PMLEAD, PMNH3, PMOUTDB, PMOUTDSC, PMPARSIZ, 
C         PMSO4SO2, QUITER 
C 
C  Input on call: 
C 
C    Common blocks: 
C    /PART1/      PRTCHK 
C    /PART3/      DBPM,EFPM,VPM 
C 
C  Output on return: 
C 
C    Common blocks: 
C    /PART3/      DBPM,EFPM,VPM 
C 
C  Local variable / array dictionary: 
C 
C   Name        Type              Description 
C  ------       ----  ----------------------------------------------------- 
C  ICHK          I    Loop for checking the required particulate input 
C 
      IMPLICIT NONE 
C 
      INCLUDE 'FLAGS3.I' 
      INCLUDE 'FLAGS5.I' 
      INCLUDE 'IDLE.I' 
      INCLUDE 'PART1.I' 
      INCLUDE 'PART3.I' 
      INCLUDE 'RFORM1.I' 
      INCLUDE 'SULFUR.I' 
C 
C  Declare parameter list 
C 
      INTEGER ICY 
      INTEGER INERR 
      INTEGER RC 
C 
C  Declare local variable 
C 
      INTEGER ICHK 
      INTEGER IT2 
      INTEGER IYR 
C 
C 
C  Initialize PM arrays 
C 
      DBPM = 0. 
      EFPM = 0. 
      VPM  = 0. 
C 
C  Check that all required particulate data has been entered. 
C 
      DO ICHK=1,3 
        IF(PRTCHK(ICHK).EQ.1) THEN 
          CALL QUITER(0.0,ICHK,184,INERR) 
        END IF 
      END DO 
C 
      IF(INERR.GT.0) RETURN 
C 
C  Read in the Ammonia (NH3) Basic Emission Rates and Sulfur Deterioration Rates. 
C  (This data could not be initialized in BLOCK DATA or MODULES due to bugs in the  
C  Lahey compilers v4.0-5.5) 
C  If any errors are encountered, abort immediately. 
C 
      CALL GETNH3BER(INERR,RC) 
      IF(RC.NE.0.OR.INERR.GT.0) RETURN 
C 
      CALL GETNH3SDR(INERR,RC) 
      IF(RC.NE.0.OR.INERR.GT.0) RETURN 
C 
C  Set the gasoline sulfur content. This has to be done here because 
C  the PM calculations can be done independently of the criteria pollutants. 
C 
      IF(CAAFLG.EQ.1 .AND. ICY.GE.CAAMY) THEN 
        IT2 = 1 
      ELSE 
        IT2 = TIER2_FLAG + 1 
      ENDIF 
C 
      IF(ICY .GT. 2015) THEN 
         S_PPM = S_FUEL(1,FUEL_P,16,IT2)  ! Post-2015, use matrix 
      ELSE IF(ICY.GE.2000) THEN 
         IYR = ICY - 1999 
         S_PPM = S_FUEL(1,FUEL_P,IYR,IT2) ! Post-1999, use matrix 
      END IF 
C 
C  Compute the particle size distribution 
C 
      CALL PMPARSIZ 
C 
C  Compute the EPA-based PM emission factors. 
C  Note: EMFAC PM = SOF + RCP + SO4. 
C 
      CALL PMEMFAC(ICY,INERR) 
C 
C  Compute the ammonia emission factors. 
C 
      CALL PMNH3(ICY)  
C 
C  Compute the brake and tire wear emission factors, 
C  adjusted for particle size. 
C 
      CALL PMBRAKETIRE 
C 
C  Compute the lead emission factors, adjusted for 
C  particle size. 
C 
      CALL PMLEAD(ICY,INERR) 
C 
C  Compute the gasoline SO4 and SO2 emission factors. 
C 
      CALL PMGSO4SO2(ICY,INERR) 
C 
C  Compute the diesel SO4 and SO2 emission factors. 
C 
      CALL PMDSO4SO2(ICY) 
C 
C  Compute GASPM (which is EMFAC PM - SO4 or more simply RCP+SOF),  
C  adjusted for particle size. 
C 
      CALL PMGASPM(ICY) 
C 
C  Compute diesel RCP and SOF, adjusted for particle size. 
C 
      CALL PMDRCPSOF(ICY) 
C 
C  Adjust the gas and diesel SO4 and SO2 emission factors 
C  for in-use fuel sulfur content. 
C 
      CALL PMSO4SO2(ICY) 
C 
C  Compute idle PM emissions, if requested 
C 
      IF(IDLEPMFLG) CALL PMIDLE(ICY) 
C 
C  Output the descriptive results. 
C 
      CALL PMOUTDSC(ICY) 
C 
C  Output the database results. 
C 
      CALL PMOUTDB(ICY) 
C 
      RETURN 
      END 
