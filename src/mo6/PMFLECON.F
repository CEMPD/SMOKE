      REAL FUNCTION PMFLECON(MY,IV,INERR) 
C 
C  PMFLECON sets the gas and diesel fuel economies,  
C  either from MOBILE6 or the user-supplied data. 
C 
C Changes: 
C 
C 06 Nov 01 AIR Task 27: Replace code to process updated fuel 
C           economies associated with CO2 work. Added INERR 
C           support to catch MPH<=0.0. 
C 13 Feb 01 AIR Task 11: New MOBILE6 function (from PART5). 
C 10 Nov 00 AIR Task 06: Updated code to MOBILE6 standards. 
C 
C  Called by DB_AGGR, DBEXHRS, PMGAS, PMDB_AGGR, 
C            PMDSO4SO2, PMGSO4SO2, PMLEADED, PMOUTDB 
C 
C    common blocks: 
C     /RLCOM1/ ROADFE 
C 
C  Output on return: 
C     PMFLECON 
C 
C  Local variable dictionary: 
C 
C   Name   Type              Description 
C  ------  ----  ------------------------------------------------------------- 
C  IM6U     I    MOBILE6 (1) or User-Supplied (2) fuel economy switch 
C  IV       I    Vehicle class 
C  MDX      I    Model year index mapping to Model year (MY) for ROADFE 
C  MY       I    Model year 
C 
C 
      IMPLICIT NONE 
C 
      INCLUDE 'IOUCOM.I' 
      INCLUDE 'RLCOM1.I' 
C 
C  Declare parameter list. 
C 
      INTEGER INERR 
      INTEGER IV 
      INTEGER MY 
C 
C  Declare local variables/arrays. 
C 
      INTEGER IM6U 
      INTEGER MDX 
C 
C  Set the MOBILE6 or User pointer 
C 
      IM6U=1                ! MOBILE6 
      IF(MPGFLAG) IM6U=2    ! User-supplied 
C 
C  Set the model year pointer 
C 
      MDX=MY-1951 
      IF(MDX.LT. 1) MDX=1 
      IF(MDX.GT.99) MDX=99 
C 
C  Determine the fuel economy  
C 
      PMFLECON = ROADFE(MDX,IV,IM6U) 
C 
C  Apply speed and temperature corrections (FUTURE MODIFICATIONS) 
C 
C 
C  Trap error that the requested/corrected fuel economy is zero or less. 
C 
      IF(PMFLECON.LE.0.) THEN 
      WRITE(IOUERR,10) PMFLECON,IV,MY 
   10 FORMAT(' Error: An invalid fuel economy value of ',F5.1/ 
     *       '        was computed by PMFLECON for vehicle'/ 
     *       '        class ',I2,' and model year ',I4) 
      INERR=INERR+1 
      END IF 
C 
      RETURN 
      END 
