C 
C  IVTYPE defines the symbolic names for individual vehicle types, 
C  composite vehicle types and the vectors that define various classes  
C  of vehicles. 
C 
C  Changes: (Last change first) 
C 
C  13 Feb 01 AIR Task 03: Added VDM6PRT mapping array 
C  08 May 00 AIR Task 03: Added VM6M5B mapping array 
C  24 Jan 00 @DynTel-MLA 1-011  Split LDDT into 2 new types; changed 
C            MAXVEH from 27 to 28. 
C  22 Oct. 98 @DynTel-HXQ 2-680 delete the vtcomml bus. 
C 
C  11 Aug 98 @DynTel-MLA 2-000  Changed HDG and HDD names to conform to 
C            new naming conventions. 
C 
C  Notes: 
C 
C     Names that begin with 'VT' are individual vehicle types. Names that 
C     begin with 'VV' are categories that include more than one vehicle 
C     type. For example, VVGAS includes VTLDGV (LDGV), VVLDGT12 (LDGT1/2), 
C     VVLDGT34 (LDGT3/4, VVHDGV (all HDGV types), and VTMC (MC). Names 
C     that begin with I are composite vehicle types. 
C 
C  10 Mar 98 DynTel@RJD 2-622 Added integer labels for composite 
C            vehicle types. Names that begin with I are composite vehicle  
C            types. 
C  12 Mar 98 DynTel@RJD 2-622 Added LDG trucks to the list of composite 
C            vehicle types.  
C 
C  The data is initialized in BD06 
C 
C  Variable Dictionary: 
C 
C   Name    Type                  Description 
C  ------   ----  ----------------------------------------------------- 
C 
C  /VTNAMS/ 
C 
C  VTHDDV2B  I    Heavy-duty Diesel Vehicles -2B 
C  VTHDGV2B  I    Heavy-duty Gasoline Vehicles - 2B 
C  VTHDDV3   I    Heavy-duty Diesel Vehicles - 3 
C  VTHDGV3   I    Heavy-duty Gasoline Vehicles - 3 
C  VTHDDV4   I    Heavy-duty Diesel Vehicles - 4 
C  VTHDGV4   I    Heavy-duty Gasoline Vehicles - 4 
C  VTHDDV5   I    Heavy-duty Diesel Vehicles - 5 
C  VTHDGV5   I    Heavy-duty Gasoline Vehicles - 5 
C  VTHDDV6   I    Heavy-duty Diesel Vehicles - 6 
C  VTHDGV6   I    Heavy-duty Gasoline Vehicles - 6 
C  VTHDDV7   I    Heavy-duty Diesel Vehicles - 7 
C  VTHDGV7   I    Heavy-duty Gasoline Vehicles - 7 
C  VTHDDV8A  I    Heavy-duty Diesel Vehicles - 8A 
C  VTHDGV8B  I    Heavy-duty Gasoline Vehicles - 8B 
C  VTHDDV8A  I    Heavy-duty Diesel Vehicles - 8A 
C  VTHDGV8B  I    Heavy-duty Gasoline Vehicles - 8B 
C  VTLDDT12  I    Light-duty Diesel Trucks 1 and 2 
C  VTLDDT34  I    Light-duty Diesel Trucks 3 and 4 
C  VTLDDV    I    Light-duty Diesel Vehicles 
C  VTLDGT1   I    Light-Duty Gasoline Trucks 1A 
C  VTLDGT2   I    Light-Duty Gasoline Trucks 1B 
C  VTLDGT3   I    Light-Duty Gasoline Trucks 2A 
C  VTLDGT4   I    Light-Duty Gasoline Trucks 2B 
C  VTLDGV    I    Light-duty Gasoline Vehicles 
C  VTMC      I    Motorcycles 
C 
C  /CVNAMS/ 
C 
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
C  /VVTYPS/ 
C 
C  VVDSL     I    All Diesel vehicles: LDDV, LDDT, HDDV. 
C  VVGAS     I    Gas vehicles, excluding MC: LDGV, LDGT1-4, HDGV. 
C  VVGASMC   I    VVGAS plus MC: LDGV, LDGT1-4, HDGV, MC. 
C  VVHDDV    I    Heavy Duty Diesel Vehicles: HDDV2B, HDDV3, ..., HDDV8B. 
C  VVHDGV    I    Heavy Duty Gasoline Vehicles: HDGV2B, HDGV3, ..., HDGV8B. 
C  VVHEAVY   I    Heavy Duty Vehicles: HDGV and HDDV. 
C  VVLDDSL   I    Light Duty Diesels: LDDV and LDDT. 
C  VVLDDT    I    Light Duty Diesel Trucks: LDDT12 and LDDT34. 
C  VVLDGAS   I    Light Duty Gasoline vehicles: LDGV, LDGT1-4. 
C  VVLDGT12  I    Light Duty Gasoline Trucks <6000 lbs: LDGT1, LDGT2. 
C  VVLDGT34  I    Light Duty Gasoline Trucks >6000 lbs: LDGT3, LDGT4. 
C  VVVEH     I    Non-trucks: LDGV, LDDV, MC. 
C..................................................................... 
C 
      INTEGER VTHDDV2B 
      INTEGER VTHDGV2B 
      INTEGER VTHDDV3 
      INTEGER VTHDGV3 
      INTEGER VTHDDV4 
      INTEGER VTHDGV4 
      INTEGER VTHDDV5 
      INTEGER VTHDGV5 
      INTEGER VTHDDV6 
      INTEGER VTHDGV6 
      INTEGER VTHDDV7 
      INTEGER VTHDGV7 
      INTEGER VTHDDV8A 
      INTEGER VTHDGV8A 
      INTEGER VTHDDV8B 
      INTEGER VTHDGV8B 
      INTEGER VTLDDT12 
      INTEGER VTLDDT34 
      INTEGER VTLDDV 
      INTEGER VTLDGT1 
      INTEGER VTLDGT2 
      INTEGER VTLDGT3 
      INTEGER VTLDGT4 
      INTEGER VTLDGV 
      INTEGER VTMC 
C 
      INTEGER VTGASBUS 
      INTEGER VTSCHOOL 
      INTEGER VTURBAN 
C 
      COMMON /VTNAMS/ VTLDGV,VTLDGT1,VTLDGT2,VTLDGT3,VTLDGT4, 
     &          VTHDGV2B,VTHDGV3,VTHDGV4,VTHDGV5,VTHDGV6,VTHDGV7, 
     &          VTHDGV8A,VTHDGV8B,VTLDDV,VTLDDT12,VTHDDV2B,VTHDDV3, 
     &          VTHDDV4,VTHDDV5,VTHDDV6,VTHDDV7,VTHDDV8A,VTHDDV8B, 
     &          VTMC,VTGASBUS,VTSCHOOL,VTURBAN,VTLDDT34 
C 
      INTEGER  ILDGV 
      INTEGER  ILDGT12 
      INTEGER  ILDGT34 
      INTEGER  IHDGV 
      INTEGER  ILDDV 
      INTEGER  ILDDT 
      INTEGER  IHDDV 
      INTEGER  IMC 
      INTEGER  ILDGT 
C 
      COMMON /CVNAMS/ ILDGV,ILDGT12,ILDGT34,IHDGV,ILDDV,ILDDT, 
     &            IHDDV,IMC,ILDGT     
C 
      INTEGER  MAXVEH 
      PARAMETER (MAXVEH = 28) 
C 
      INTEGER  VVDSL(MAXVEH) 
      INTEGER  VVGAS(MAXVEH) 
      INTEGER  VVGASMC(MAXVEH) 
      INTEGER  VVHDDV(MAXVEH) 
      INTEGER  VVHDGV(MAXVEH) 
      INTEGER  VVHEAVY(MAXVEH) 
      INTEGER  VVLDDSL(MAXVEH) 
      INTEGER  VVLDDT(MAXVEH) 
      INTEGER  VVLDGAS(MAXVEH) 
      INTEGER  VVLDGT12(MAXVEH) 
      INTEGER  VVLDGT34(MAXVEH) 
      INTEGER  VVSTARTS(MAXVEH) 
      INTEGER  VM6M5B(MAXVEH) 
      INTEGER  VDM6PRT(MAXVEH) 
C  There are no current references in the code to the next three vehicle  
C  type categories. 
C 
      COMMON /VVTYPS/ VVDSL, VVGAS, VVGASMC, VVHDDV, VVHDGV, 
     &         VVHEAVY, VVLDDSL, VVLDGAS, VVLDGT12, VVLDGT34, 
     &         VVSTARTS,VVLDDT,VM6M5B,VDM6PRT 
