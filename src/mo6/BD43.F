      BLOCK DATA BD43 
C 
C  BLOCK DATA Subprogram 43: /PART1/,/PART2/,/PART3/ 
C 
C  BLOCK DATA 43 initializes the particulate commons. 
C 
C  Changes: (Last change first) 
C 
C  25 May 01 AIR Task 17: Extended PRTCHK to 3 necessary inputs. 
C  09 May 01 AIR Task 11/5: Updated to handle idle. 
C  03 May 01 AIR Task 12: Added BERSULF initialization. 
C  27 Feb 01 AIR Task 11: Updated to hold fuel sulfur levels. 
C  09 Feb 01 AIR Task 11: Updated to hold tampering information. 
C  31 Jan 01 AIR Task 11: New Block Data. 
C 
C  Common block array subscripts: 
C 
C  BERSULF(25)  -  BERSULF(MY,IV) 
C  EFPM(11,28)  -  EFPM (IPPM,IV) 
C  FRACA(4)     -  FRACA (I) 
C  IVEHWL(28)   -  IVEHWL(IV) 
C  PRTCHK(3)    -  PRTCHK (ICHK) 
C  VPM(28)      -  VPM(IV) 
C  FUEL_S(3,29) -  FUEL_S(fuel type, calendar year)  
C 
C  Common block dictionary: 
C 
C   Name    Type                  Description 
C  ------   ----  ----------------------------------------------------- 
C  /PART1/ 
C  BERSULF   R    Fuel sulfur level associated with the basic emission factor (BER). 
C  D_PPM     R    diesel sulfur content, in ppm (required input) 
C  IVEHWL    I    Average Number of Vehicle Wheels for 28 Vehicle Classes 
C  MAXIPPM   I    Maximum number of particulate pollutants, as with IPPM as follows: 
C                 1=SO4, 2=OCARBON, 3=ECARBON, 4=GASPM, 5=Lead, 6=SO2, 7=MH3, 8=Brake, 9=Tire 
C                 10=EMFAC_PM, 11=Idle 
C  PARTFLG   L    switch to run the particulate subroutine 
C  PRTCHK    I    checklist of required particulate run inputs (1=not input, 2=input) 
C                 (1=particle size, 2=diesel sulfur level, 3=ZML/DR files) 
C  PSZCUT    R    particle size cutoff, in micrometers (required input) 
C 
C  /PART2/ 
C  BSWPTR   R    Base gasoline fuel sulfur content, in wt% 
C  CTLFRC   R    Catalyst fraction for a given IV and MY  
C  FECONM   R    Fuel economy, mpg 
C  FRACA    R    Fraction of vehicles with: 
C                 1: Ox No Air 
C                 2: Ox Air 
C                 3: 3w No Air 
C                 4: 3w Air 
C  GSWPTR   R    In-use gasoline fuel sulfur content, in wt% 
C  RMIS     R    Misfueling rate 
C  SULFCF   R    Fuel Sulfur Correction Factor 
C  TAMFRC   R    PM Tampering fraction 
 
C /PART3/ 
C   DPPM    R    PM Database Emissions 
C   EFPM    R    PM Emissions weighted for each veh. class 
C   VPM     R    PM Emissions weighted over all vehicle classes 
C 
C 
      INCLUDE 'PART1.I' 
      INCLUDE 'PART2.I' 
      INCLUDE 'PART3.I' 
C 
C  /PART1/ 
C 
      DATA BERSULF /700*0./ 
      DATA D_PPM   /0./ 
      DATA FUEL_S  /87*0./ 
      DATA IVEHWL  / 
     * 6* 4, !LDGV,LDGT1-LDGT4,HDGV2B 
     * 5* 6, !HDGV3-HDGV7 
     * 2*18, !HDGV8A-HDGV8B 
     * 3* 4, !LDDV,LDDT12,HDDV2B 
     * 5* 6, !HDDV3-HDDV7 
     * 2*18, !HDDV8A-HDDV8B 
     *    2, !MC 
     * 3* 6, !HDGB,HDDBT,HDDBS 
     *    4/ !LDDT34 
      DATA MAXIPPM /11/ 
      DATA PARTFLG /.FALSE./ 
      DATA PRTCHK  /3*1/ 
      DATA PSZCUT  /0./ 
C 
C  /PART2/ 
C 
      DATA BSWGTR /0./ 
      DATA CTLFRC /0./ 
      DATA FECONM /0./ 
      DATA FRACA  /4*0./ 
      DATA GSWGTR /0./ 
      DATA RMIS   /0./ 
      DATA SULFCF /0./ 
      DATA TAMFRC /0./ 
C 
C  /PART3/ 
C 
      DATA DBPM /672000*0./ 
      DATA EFPM /   308*0./ 
      DATA VPM  /    11*0./ 
C 
      END 
