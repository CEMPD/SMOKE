      BLOCK DATA BD06 
C 
C 
C 
C  BLOCK DATA Subprogram 06:  /VTNAMS/,/CVNAMS/,/VVTYPS/ 
C 
C  BD06 defines and initializes the symbolic names for individual vehicle  
C  types, composite vehicle types and the vectors that define various  
C  classes of vehicles. 
C 
C  Changes: (Last change first) 
C 
C  13 Feb 01 AIR Task 11:       Added VDM6M5B array for diesel particulate 
C                               and gasoline particulate mapping 
C  20 jUN 00 @AIR-dK 1-03       Added VM6M5B array for temperature 
C                               correction 
C   2 Feb 00 @DynTel-ddj 1-032  Changed vehicle type 28 to 0 for 
C            VVBUS. 
C  24 Jan 00 @DynTel-MLA 1-011  Split LDDT into 2 new types. 
C  22 Oct. 98 @DynTel-HXQ 2-680  delete diesel bus - vtcomml. 
C  14 Sep 98 @DynTel-RD 2-660  Added light duty gasoline trucks to VVSTARTS. 
C  11 Aug 98 @DynTel-MLA 2-000  Changed HDG and HDD names to conform to 
C            new naming conventions. 
C 
C  Common block array subscripts: 
C 
C     VVBUS    (MAXVEH) 
C     VVCNG    (MAXVEH) 
C     VVDSL    (MAXVEH) 
C     VVGAS    (MAXVEH) 
C     VVGASMC  (MAXVEH) 
C     VVHDDV   (MAXVEH) 
C     VVHDGV   (MAXVEH) 
C     VVHEAVY  (MAXVEH) 
C     VVLDDSL  (MAXVEH) 
C     VVLDDT   (MAXVEH) 
C     VVLDGAS  (MAXVEH) 
C     VVLDGT12 (MAXVEH) 
C     VVLDGT34 (MAXVEH) 
C     VVVEH    (MAXVEH) 
C 
C  Notes: 
C 
C     Names that begin with 'VT' are individual vehicle types. Names that 
C     begin with 'VV' are arrays that define related groups of vehicles. 
C 
C  Variable Dictionary: 
C 
C   Name    Type                  Description 
C  ------   ----  ----------------------------------------------------- 
C 
C  /VTNAMS/ 
C 
C  VTHDDV2B  I    Heavy Duty Diesel Vehicles - 2B 
C  VTHDGV2B  I    Heavy Duty Gasoline Vehicles - 2B 
C  VTHDDV3   I    Heavy Duty Diesel Vehicles - 3 
C  VTHDGV3   I    Heavy Duty Gasoline Vehicles - 3 
C  VTHDDV4   I    Heavy Duty Diesel Vehicles - 4 
C  VTHDGV4   I    Heavy Duty Gasoline Vehicles - 4 
C  VTHDDV5   I    Heavy Duty Diesel Vehicles - 5 
C  VTHDGV5   I    Heavy Duty Gasoline Vehicles - 5 
C  VTHDDV6   I    Heavy Duty Diesel Vehicles - 6 
C  VTHDGV6   I    Heavy Duty Gasoline Vehicles - 6 
C  VTHDDV7   I    Heavy Duty Diesel Vehicles - 7 
C  VTHDGV7   I    Heavy Duty Gasoline Vehicles - 7 
C  VTHDDV8A  I    Heavy Duty Diesel Vehicles - 8A 
C  VTHDGV8A  I    Heavy Duty Gasoline Vehicles - 8A 
C  VTHDDV8B  I    Heavy Duty Diesel Vehicles - 8B 
C  VTHDGV8B  I    Heavy Duty Gasoline Vehicles - 8B 
C  VTLDDT12  I    Light Duty Diesel Trucks 1 and 2 
C  VTLDDT34  I    Light Duty Diesel Trucks 3 and 4 
C  VTLDDV    I    Light Duty Diesel Vehicles 
C  VTLDGT1   I    Light Duty Gasoline Trucks 1 (previously 1A) 
C  VTLDGT2   I    Light Duty Gasoline Trucks 2 (previously 1B) 
C  VTLDGT3   I    Light Duty Gasoline Trucks 3 (previously 2A) 
C  VTLDGT4   I    Light Duty Gasoline Trucks 4 (previously 1B) 
C  VTLDGV    I    Light Duty Gasoline Vehicles 
C  VTMC      I    Motorcycles 
C 
C  VTGASBUS  I    Gasoline-powered buses 
C  VTURBAN   I    Diesel bus - intra-urban 
C  VTSCHOOL  I    Diesel school buses 
C 
C  VTCNG1    I    Compressed natural gas vehicles - category 1 
C  VTCNG2    I    Compressed natural gas vehicles - category 2 
C  VTCNG3    I    Compressed natural gas vehicles - category 3 
C  VTCNG4    I    Compressed natural gas vehicles - category 4 
C  VTCNG5    I    Compressed natural gas vehicles - category 5 
C  VTCNG6    I    Compressed natural gas vehicles - category 6 
C  VTCNG7    I    Compressed natural gas vehicles - category 7 
C  VTCNG8    I    Compressed natural gas vehicles - category 8 
C 
C  /VVTYPS/ 
C 
C  VVBUS     I    All buses - gasoline, school, urban, and commercial. 
C  VVDSL     I    All Diesel vehicles: LDDV, LDDT12, LDDT34 and HDDV. 
C  VVGAS     I    Gas vehicles, excluding MC: LDGV, LDGT, HDGV. 
C  VVGASMC   I    VVGAS plus MC: LDGV, LDGT, HDGV, MC. 
C  VVHDDV    I    Heavy Duty Diesel Vehicles: HDDV2B, HDDV3, ..., HDDV8B. 
C  VVHDGV    I    Heavy Duty Gasoline Vehicles: HDGV2B, HDGV3, ..., HDGV8B. 
C  VVHEAVY   I    Heavy Duty Vehicles: HDGV, HDDV. 
C  VVLDDSL   I    Light-duty Diesel vehicles: LDDV, LDDT12 and LDDT34. 
C  VVLDDT    I    Light-duty Diesel Trucks: LDDT12, LDDT34 
C  VVLDGAS   I    Light-duty gasoline vehicles: LDGV, LDGT1-4. 
C  VVLDGT12  I    Light-duty gasoline trucks <6000 lbs: LDGT1, LDGT2. 
C  VVLDGT34  I    Light-duty gasoline trucks >6000 lbs: LDGT3, LDGT4. 
C  VVVEH     I    Non-trucks: LDGV, LDDV, MC. 
C  VVSTARTS  I    Vehicles with seperated start and running emissions 
C  VM6M5B    I    Maps the gasoline MOBILE6 vehicle classes into the 
C                 MOBILE5b/PART5 classes for TCF and gasoline particulate. 
C  VDM6PRT   I    Maps the diesel MOBILE6 vehicle classes into the PART5 
C                 classes for diesel particulate. 
C   
C  /CVNAMS/ 
C  ILDGV     I    All Light-duty Gasoline Vehicles 
C  ILDGT12   I    All Light-Duty Gasoline Trucks 1 
C  ILDGT34   I    All Light-Duty Gasoline Trucks 2 
C  IHDGV     I    All Heavy-duty Gasoline Vehicles 
C  ILDDV     I    All Light-duty Diesel Vehicles 
C  ILDDT     I    All Light-duty Diesel Trucks 
C  IHDDV     I    All Heavy-duty Diesel Vehicles 
C  IMC       I    All Motorcycles 
C  ILDGT     I    All Light-Duty Gasoline Trucks 
C 
C  Note: 
C       
C  10 Mar 98 DynTel@RJD 2-622 Added integer labels for composite 
C            vehicle types. 
C  12 Mar 98 DynTel@RJD 2-622 Added LDG trucks to the list of 
C            composite vehicle types. 
C  14 July 98 Dyntel@MMS 2-666 added 1, for running, in array VVSTARTS 
C            for VI = 14,15,24 
C  23 Jul 98 @DynTel-MLA 2-000  Changed LDGT names to conform to 
C            certification naming conventions. 
C 
C  /CVNAMS/ 
C 
C  ILDGV     I         All Light-duty Gasoline Vehicles 
C  ILDGT12   I    All Light-Duty Gasoline Trucks (<6000 lbs. GVW) 
C  ILDGT34   I    All Light-Duty Gasoline Trucks (>6000 lbs. GVW) 
C  IHDGV     I    All Heavy-duty Gasoline Vehicles 
C  ILDDV     I    All Light-duty Diesel Vehicles 
C  ILDDT     I    All Light-duty Diesel Trucks  
C  IHDDV     I    All Heavy-duty Diesel Vehicles 
C  IMC       I    All Motorcycles 
C  ILDGT     I    All Light-Duty Gasoline Trucks 
C 
C .................................................................... 
C 
      IMPLICIT NONE 
      INCLUDE 'IVTYPE.I' 
C 
C     Note: If the numbers assigned to different vehicle types are 
C     changed, this defines a permutation on the 'VV' vectors,  
C     initialized immediately below, and on the vectors IVPTRA, IVPTRB, 
C     IVPTRL, and IVPTRT initialized in BD16 (block data 16). 
C 
      DATA VTLDGV,VTLDGT1,VTLDGT2,VTLDGT3,VTLDGT4    /1, 2, 3, 4, 5/ 
      DATA VTHDGV2B,VTHDGV3,VTHDGV4,VTHDGV5,VTHDGV6  /6, 7, 8, 9,10/ 
      DATA VTHDGV7,VTHDGV8A,VTHDGV8B,VTLDDV,VTLDDT12 /11,12,13,14,15/ 
      DATA VTHDDV2B,VTHDDV3,VTHDDV4,VTHDDV5,VTHDDV6  /16,17,18,19,20/ 
      DATA VTHDDV7,VTHDDV8A,VTHDDV8B,VTMC            /21,22,23,24/ 
      DATA VTGASBUS,VTURBAN,VTSCHOOL,VTLDDT34        /25,26,27,28/ 
C 
      DATA ILDGV,ILDGT12,ILDGT34,IHDGV,ILDDV,ILDDT   /1,2,4,6,14,15/ 
      DATA IHDDV,IMC,ILDGT                           /16,24,29/ 
C 
      DATA VVDSL   /0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,0, 
     &             0,1,1,1/ 
C 
      DATA VVGAS   /1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0, 
     &             1,0,0,0/ 
C 
      DATA VVGASMC /1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,1, 
     &             1,0,0,0/ 
C 
      DATA VVHDDV  /0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,0, 
     &             0,1,1,0/ 
C 
      DATA VVHDGV  /0,0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0, 
     &             1,0,0,0/ 
C 
      DATA VVHEAVY /0,0,0,0,0,1,1,1,1,1,1,1,1,0,0,1,1,1,1,1,1,1,1,0, 
     &             1,1,1,0/ 
C 
      DATA VVLDDSL /0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0, 
     &             0,0,0,1/ 
C 
      DATA VVLDDT  /0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0, 
     &             0,0,0,1/ 
C 
      DATA VVLDGAS /1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, 
     &             0,0,0,0/ 
C 
      DATA VVLDGT12 /0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, 
     &             0,0,0,0/ 
C 
      DATA VVLDGT34 /0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, 
     &             0,0,0,0/ 
C     
      DATA  VVSTARTS /1,1,1,1,1,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0, 
     &             1,0,0,0,1/ 
C 
      DATA  VM6M5B /1,2,2,3,3,4,4,4,4,4,4,4,4,0,0,0,0,0,0,0,0,0,0, 
     &             5,4,0,0,0/ 
C 
      DATA VDM6PRT /0,0,0,0,0,0,0,0,0,0,0,0,0,6,7,8,9,9,9,10,10,10,11,0, 
     &             0,12,12,7/ 
C 
C  Currently, there are no references in the code to the two following 
C  vehicle type categories. 
C 
CCC   DATA VVVEH   /1,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,1, 
CCC  &             0,0,0,0/ 
C 
CCC   DATA VVBUS   /0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, 
CCC  &             1,1,1,0/ 
C 
      END 
