      BLOCK DATA BD35 
C 
C  BLOCK DATA Subprogram 35: /BYMYC1/,/BYMYC2/,/BYMYC3/,/BYMYC4/,/BYMYC5/ 
C 
C  BLOCK DATA 35 defines the By Model Year Table generation common blocks. 
C 
C  Changes: (Last change first) 
C 
C 16 Jun 00 @EPA-BAG  Removed IVPICK. This variable is obsolete in Mobile6. 
C   2 Feb 00 @DynTel-ddj 1-032  Changed LDDT to LDDT12 and created 
C            LDDT34 by copying LDDT12 data to place holder for 
C            vehicle type 28. 
C   5 Aug 99 @DynTel-ddj 2-699  Removed array BYTAM. 
C  22 Feb 99 @DynTel-ddj 2-696 Added hourly index to BYTAM. 
C  06 Nov 98 @Dyntel-BAG 2-000 light duty truck naming; replaced composite 
C            Light Duty Gas Trucks 1/2 and Light Duty Gas Trucks 3/4 with 
C            Light Duty Gas Trucks 1 and Light Duty Gas Trucks 3, 
C            respectively. 
C  22 Oct. 98 @DynTel-HXQ 2-680 replace vtcomml bus with school bus in BYVTN. 
C 
C  Common block array subscripts: 
C 
C  BYBEF4(3,25,28,2) -  BYBEF4 ( IP, JDX, IV, ISR ) 
C  BYCUM(25,28)      -  BYCUM  (     JDX, IV ) 
C  BYEVAP(25,28)     -  BYEVAP (     JDX, IV ) 
C  BYFER(3,25,28)    -  BYFER  ( IP, JDX, IV ) 
C  BYIMCR(3,25,14)   -  BYIMCR ( IPIM, JDX, IVIM ) 
C  BYIMRE(3,25,14)   -  BYIMRE ( IPIM, JDX, IVIM ) 
C  BYREFL(25,14)     -  BYREFL (     JDX, IVTAM ) 
C  BYRSTL(25,28)     -  BYRSTL (     JDX, IV ) 
C  BYRUNL(25,14)     -  BYRUNL (     JDX, IVTAM ) 
C  BYVTN(28)         -  BYVTN  (          IV ) 
C  IVPICK(28)        -  IVPICK (          IV ) 
C  PN(10)            -  PN     ( ICH ) 
C  PSEC(3)           -  PSEC   ( IP ) 
C  REDSUM(3,14)      -  REDSUM ( IPIM, IVIM ) 
C 
C  Common block dictionary: 
C 
C   Name   Type               Description 
C  ------  ----  -------------------------------------------- 
C  /BYMYC1/: 
C 
C  IVPICK   I    By Model Year (BYMY) Table vehicle class case selection vector. 
C                For each IV, 1 = skip the Table(s), 2 = print the Table(s), 
C                and 3 = print both BYMY and I/M tables (IVIM only). 
C  BYCUM    R    BYMY version of CUMMIL: same except zeroed out where the 
C                corresponding TF is zero.  BYCUM is used as a "clean" 
C                proxy for CUMMIL in the BYMY tables.  CUMMIL could be used 
C                directly if it was recomputed every scenario TF was redone, 
C                which is not true when the user enters AMAR. 
C 
C  /BYMYC2/: 
C 
C  BYBEF4   R    by model year basic emission rate for given IDX x IV x IV, 
C                corrected for RVP, temperature, operating mode and, if 
C                applicable, methane and / or I/M: BEF4 = BEF - OMTTAM 
C  BYFER    R    by model year final emission rate for IDX x IV x IP = COMPEF 
C 
C  /BYMYC3/: 
C 
C  BYEVAP   R    by model year evaporative HC component 
C  BYRUNL   R    by model year running loss HC component = RNLOSS 
C  BYREFL   R    by model year refueling loss HC component = RLRATE 
C  BYRSTL   R    by model year resting loss HC component = RSTLOS 
C 
C  /BYMYC4/: 
C 
C  BYVTN   C*26  full BYMY Table Vehicle Type names 
C  YFM     C*17  column titles for model year, travel fraction & cumulative 
C                mileage 
C  PN      C*7   pollutant names for column headings of BYMY tables 
C  U7      C*7   underline of column heading title in BYMY table 
C  PSEC    C*3   pollutant names for section headings of BYMY table 
C  VB      C*1   vertical bar 
C 
C 
C  Notes: 
C 
C  BD35 was added for MOBILE4 In-House Version 02 (By Model Year Tables). 
C  BD 35 provide the parameters and by model year results save areas. 
C 
C  9 Dec 96 DynTel@MLA 2-622  Expanded the arrays to accomodate the new  
C           LDGT, HDGV, and new BUS vehicle categories. 
C 
C 26 Feb 97 DynTel@MLA  2-622 Deleted IMPICK. For the Mobile6 transition 
C           code, each element in the IVPICK array can take on values 1 (Do 
C           not print BYMY tables for the vehicle type), or 2 (print BYMY  
C           tables). IVIM vehicles can also use the value 3 (print BYMY  
C           and I/M tables for a particular vehicle type). 
C 20 Apr 98 DynTel@BAG 2-661 Updated the names of the light duty gas truck. 
C           The new names are LDGT1(old name LDGT1A), LDGT2(old name LDGT1B),  
C           LDGT3(old name LDGT2A), and LDGT4(old name LDGT2B). 
C 23 Apr 98 @DynTel-RJD 2-654 Added an extra dimension to the BYBEF array 
C           to distinguish between start and running emissions. 
C 30 Jun 98 DynTel@BAG 2-000 Deleted I/M credit block data, BYMYC5.I, and 
C           variables associated with it, BYIMCR, BYIMRE, and REDSUM 
C  4 Aug 98 DynTel@MLA 2-000 Corrected the LDGT names. 
C 
      IMPLICIT NONE 
      INCLUDE 'BYMYC1.I' 
      INCLUDE 'BYMYC2.I' 
      INCLUDE 'BYMYC3.I' 
      INCLUDE 'BYMYC4.I' 
C 
C  /BYMYC1/: 
C 
      DATA BYCUM/700*0.0/ 
C 
C  /BYMYC2/: 
C 
      DATA BYBEF4/4200*0.0/,BYFER/2100*0.0/ 
C 
C  /BYMYC3/: 
C 
      DATA BYEVAP,BYRUNL,BYREFL,BYRSTL/2100*0.0/ 
C 
C 
C  /BYMYC4/: 
C 
C  The following are the character strings for the vehicle type names. The 
C  continuation codes are the old Mobile5 vehicle type numbers (1=LDGV, 
C  2=LDGT1, 3=LDGT2, ... The new code B is used for buses.) 
C 
      DATA BYVTN/ 
     1            'Light Duty Gas Vehicles   ', 
     2            'Light Duty Gas Trucks 1   ', 
     2            'Light Duty Gas Trucks 2   ', 
     3            'Light Duty Gas Trucks 3   ', 
     3            'Light Duty Gas Trucks 4   ', 
CCC  4            '2B Heavy Duty Gas Veh.    ', 
     4            'Heavy Duty Gas Vehicles   ', 
     4            '3 Heavy Duty Gas Veh.     ', 
     4            '4 Heavy Duty Gas Veh.     ', 
     4            '5 Heavy Duty Gas Veh.     ', 
     4            '6 Heavy Duty Gas Veh.     ', 
     4            '7 Heavy Duty Gas Veh.     ', 
     4            '8A Heavy Duty Gas Veh.    ', 
     4            '8B Heavy Duty Gas Veh.    ', 
     5            'Light Duty Diesel Vehicles', 
     6            'Light Duty Diesel Trucks12', 
CCC  7            '2B Heavy Duty Diesel Veh. ', 
     7            'Heavy Duty Diesel Vehicles', 
     7            '3 Heavy Duty Diesel Veh.  ', 
     7            '4 Heavy Duty Diesel Veh.  ', 
     7            '5 Heavy Duty Diesel Veh.  ', 
     7            '6 Heavy Duty Diesel Veh.  ', 
     7            '7 Heavy Duty Diesel Veh.  ', 
     7            '8A Heavy Duty Diesel Veh. ', 
     7            '8B Heavy Duty Diesel Veh. ', 
     8            'Motorcycles               ', 
     B            'Gasoline Bus              ', 
     B            'Urban Bus                 ', 
     B            'School Bus                ', 
     6            'Light Duty Diesel Trucks34'/ 
C    
      DATA YFM/'Year   TF   Miles'/ 
      DATA PN/' BEF4  ','Tamper ','SALHCF ','Evapor ','Refuel ', 
     *        'Runnin ','Restin ','  FER  ','Credit ','Reduc  '/ 
      DATA U7/' ______'/ 
      DATA PSEC/'HC ','CO ','NOX'/ 
      DATA VB/'|'/ 
C 
      END 
