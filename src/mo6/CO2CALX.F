      SUBROUTINE CO2CALX(ICY,INERR) 
C 
C  This subroutine returns the CO2 emission factors (in g/mi) 
C  for calendar year ICY. 
C 
C Changes:
C 
C 24 Sep 03 BAG@EPA commented out CO and HC subtractions from CO2 
C 04 Jan 01 AIR Task 27: Removed unused variables/code 
C 06 Nov 01 AIR Task 27: New MOBILE6 routine. 
C 
C  Called by CO2CALX. 
C 
C  Calls CO2HRLOOP 
C 
C  Input on call: 
C    common blocks: /CO2COM/ CO2 
C 
C 
C  Output on return: CO2 
C 
C     Local Variable Dictionary: 
C 
C  Variable    Type      Description 
C ----------------------------------------------------------------------- 
C   IDX        I         Model year indicator 
C   IV         I         Vehicle class 
C   MY         I         Model year 
C   JDX        I         Age 
C 
C 
      IMPLICIT NONE 
C 
      INCLUDE 'CO2COM.I'  ! DBCO2, CO2HEF, MPG, FMPG 
      INCLUDE 'IVTYPE.I'  ! MAXVEH, VVGASMC 
      INCLUDE 'MAXIMA.I'  ! MAXIH, MAXYRS 
      INCLUDE 'MYRCAL.I'  ! TF 
      INCLUDE 'RESUL1.I'  ! EFFTP 
      INCLUDE 'RESUL3.I'  ! VFTP 
      INCLUDE 'VMXCOM.I'  ! VMTMIX 
      INCLUDE 'SPEED9.I'  ! HVMT 
C 
C  Declare parameter list. 
C 
      INTEGER ICY 
      INTEGER INERR 
C 
C  Externals 
C 
      REAL, EXTERNAL :: PMFLECON 
C 
C  Declare local variables/arrays. 
C 
      INTEGER IDX 
      INTEGER IH 
      INTEGER IROAD 
      INTEGER IV 
      INTEGER JDX 
      INTEGER MY 
      REAL    CO2EF 
      REAL    RFE 
      REAL    SUMIDX 
      REAL    SUMIH 
      REAL    SUMMPG(28) 
C 
C  Initialize the database CO2 array and "All Vehicles" to zero. 
C 
      DBCO2=0.0 
      VFTP(4)=0.0 
C 
C  Initialize the harmonic mean MPG arrays and variables. 
C 
      SUMMPG = 0. 
      MPG = 0. 
      FMPG = 0. 
C 
C  Loop for each vehicle class 
C 
      DO IV=1,MAXVEH 
C 
C  Initialize 
C 
      EFFTP(4,IV)=0. 
C 
C  Skip if VMTMIX is zero 
C 
      IF(VMTMIX(IV).EQ.0.) CYCLE 
C 
C  Loop for all ages 
C 
      DO IDX=1,MAXYRS 
C 
C  If travel fraction for this age/vehicle is 0, skip. 
C 
          IF(TF(IDX,IV) .EQ. 0.) CYCLE 
C 
C  Compute JDX (age) 
C 
      JDX=(MAXYRS+1)-IDX 
C 
C  Determine model year 
C 
      MY=ICY-JDX+1 
C 
C  Compute the THC and CO for all hours weighted together for this IDX. 
C  (THC does not include evap or refueling) 
C 
      CALL CO2HRLOOP(IV,MY,IDX,INERR) 
C 
      IF(INERR.GT.0) RETURN 
C 
C  Compute the Road Fuel Economy 
C 
          RFE=PMFLECON(MY,IV,INERR) 
C 
          IF(INERR.GT.0) RETURN 
C 
C  Summate the inverse Road Fuel Economy weighted by Travel Fraction to 
C  set up for computing the harmonic mean MPG below. 
C 
      SUMMPG(IV) = SUMMPG(IV) + (1./RFE) * TF(IDX,IV) 
C 
C  Compute the gas and diesel CO2. Since CO2HRLOOP computes the 
C  THC and CO weighted over the entire day, CO2 is hour (and road) 
C  independent. 
C 
            IF(VVGASMC(IV).EQ.1) THEN 
              CO2EF = 8868.132/RFE      ! Gas 
CCCC     *              - 3.1722 *CO2HEF(1) ! THC adjustment 
CCCC     *              - 1.57143*CO2HEF(2) ! CO adjustment 
            ELSE 
              CO2EF = 10175.824/RFE     ! Diesel 
CCCC     *              - 3.1722 *CO2HEF(1) ! THC adjustment 
CCCC     *              - 1.57143*CO2HEF(2) ! CO adjustment 
            END IF 
C 
C  Initialize 
C 
          SUMIDX=0. 
C 
C  Loop for each hour 
C 
          DO IH=1,MAXIH 
C 
          SUMIH=0. 
C 
C  Loop for each road type 
C 
          DO IROAD=1,4 
C 
C  Store in database and EFFTP arrays. 
C 
           DBCO2(IROAD,IH,IDX,IV) = CO2EF 
C 
           SUMIH = SUMIH + CO2EF * FVMT(IROAD,IH,IV) 
C 
        END DO  !IROAD 
C 
        SUMIDX = SUMIDX + SUMIH * HVMT(IH) 
C 
        END DO  !IH 
C 
        EFFTP(4,IV) = EFFTP(4,IV) + SUMIDX * TF(IDX,IV) 
C 
       END DO ! IDX 
C 
       VFTP(4) = VFTP(4) + EFFTP(4,IV) * VMTMIX(IV) 
C 
       FMPG = FMPG + SUMMPG(IV) * VMTMIX(IV) 
C 
      END DO ! IV 
C 
C  Compute the harmonic mean MPGs 
C 
      DO IV=1,MAXVEH 
      IF(VMTMIX(IV).GT.0..AND.SUMMPG(IV).GT.0.) MPG(IV)=1./SUMMPG(IV) 
      END DO 
      IF(FMPG.GT.0.) FMPG = 1./FMPG 
C 
C  Output the database results. 
C 
      CALL CO2OUTDB(ICY,INERR) 
C 
      RETURN 
      END 
