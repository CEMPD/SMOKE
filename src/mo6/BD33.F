      BLOCK DATA BD33 
C 
C  BLOCK DATA Subprogram 33: /REGION/,/SCENE2/ 
C 
C  BLOCK DATA 33 defines common blocks related to the first scenario record. 
C 
C  Changes: (Last change first) 
C 
C  13 Apr 00 AIR Task 02: Set default ALTITUDE (IREJN) to Low (1) 
C  22 Jun 99 @DynTel-MLA 2-694  Removed Lev variables unused in Mobile6. 
C  26 May 99 @DynTel-HQ 2-698 delete code for initializing spd, pccn 
C                             pchc, pccc 
C  15 Apr 99 @DynTel-HQ 2-698 operating mode data are initialized 
C 
C  Common block array subscripts: 
C 
C  Common block dictionary: 
C 
C   Name   Type              Description 
C  ------  ----  -------------------------------------------------- 
C  /REGION/: 
C 
C  IREJN    I    user supplied region code: 1 = low altitude, 2 = high 
C                1 = low altitude region 
C                2 = high altitude region 
C  INITPR   I    switch: 1 = 1st scenario's output pass, print run parameters 
C                        0 = 2nd or later, do not 'CALL OUTHD3' 
C 
C 
C  /SCENE2/: 
C 
C  FCC      R    fraction of vmt in cold start mode, vehicle catalyst equipped 
C  FCN      R    fraction of vmt in cold start mode, non-catalyst equipped 
C  FHC      R    fraction of vmt in hot start mode, catalyst equipped 
C  FHN      R    fraction of vmt in hot start mode, non-catalyst equipped 
C 
C  Notes: 
C 
C  3 Dec 96 DynTel@MLA 2-622  Expanded the arrays to accomodate the new 
C           Heavy Duty Gasoline, Heavy Duty Diesel, LDGT, and BUS vehicle  
C           categories. 
C  16 December 1997 Dyntel@MMS  2-645 modified speed from a multiple speed 
C                   to a single speed. SPD array reduced to a single variable. 
C 
C  23 Jul 98 @DynTel-MLA 2-000  Changed LDGT names to conform to 
C            certification naming conventions. 
C 
C 24-Aug-98 @Dyntel-nh Request 0-632 Removed dead code related to special I/M 
C 
      IMPLICIT NONE 
      INCLUDE 'REGION.I' 
      INCLUDE 'SCENE2.I' 
C 
C  Only INITPR and IREJN need to be initialized. 
C 
      DATA INITPR /1/ 
      DATA IREJN /1/ 
C 
C  The rest are read in or set during execution, so dummy zeroes are used here. 
C 
      DATA FCC,FCN,FHC,FHN/4*0.0/ 
C 
      END 
