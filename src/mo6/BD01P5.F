      BLOCK DATA BD01P5 
C 
C  BLOCK DATA Subprogram 01: /DFSBLK/ 
C 
C  BLOCK DATA 01 initializes 
C     1) the discretionary fuel switching rates 
C     2) the fraction of lead burned that is exhausted 
C 
C Changes: 
C 
C 12 Feb 01: AIR Task 11: Changed last model year from 2020 to NEYR. 
C 10 Nov 00 AIR Task 06: Updated code to MOBILE6 standards. 
C 
C  Common block array subscripts: 
C 
C  DFS(4,5) - DFS (IYPT, IV) 
C  MYLEAD(4,5) - MYLEAD (IYPT, IV) 
C 
C  Common block dictionary: 
C 
C   Name   Type                Description 
C  ------  ----  ---------------------------------------------- 
C  /DFSBLK/: 
C 
C  DFS      R    Fraction of vehicles designed to run on leaded 
C                gasoline that use unleaded gasoline, i.e. 
C                discretionary fuel switching rates, 1-DFS is the 
C                fraction of vehicles designed to run on leaded 
C                gasoline that use leaded gasoline 
C  MYLEAD   I    Model year groups for discretionary fuel switching 
C  ABURN    R    Fraction of lead burned that is exhausted 
C                (denoted as a in particulate document) 
C  MAXLED   I    Number of model year groups for base lead factors 
C 
C 
C  Notes: 
C  Nov-29-93 @ ARC-bsg Subtask 340 
C  Changed DFS rates and lead start years for LDGT2 
C  this fixes errors in SAI code delivered 10-22-93 that was not 
C  using what the comments said they were using 
C  Jul-1-94 @ CSC-bsg Request 447 Investigate graphical trends in the 
C  particulate model by comparing IMFLAG=1 to IMFLAG=2 for various  
C  calendar years.  Discovered a mistake with HDGV DFS rate. 
C 
      IMPLICIT NONE 
C 
      INCLUDE 'DFSBLK.I' 
      INCLUDE 'PARAME.I' 
C 
      DATA ABURN/0.75,0.40,0.44/ 
      DATA MAXLED/4/ 
C 
C  last year by vehicle type for base lead fueling percents. 
C 
      DATA MYLEAD/ 
C        ---- LDGV ----       ---- LDGT1 ---       --- LDGT2 ---- 
     1 1970,1974,1989,NEYR,   1970,1974,1989,NEYR,1970,1978,1989,NEYR, 
C        ---- HDGV ----       ---- MC ---- 
     2 1989,NEYR,   0,   0,   1977,1989,NEYR,   0/ 
C 
      DATA DFS/ 
C  LDGV 
     * 0.113,0.084,0.276,1.000, 
C  LDGT1 
     * 0.113,0.084,0.276,1.000, 
C  LDGT2 
     * 0.113,0.084,0.276,1.000, 
C  HDGV 
     * 0.000,1.000,1.000,1.000, 
C  MC 
     * 0.802,0.954,1.000,1.000/ 
C 
      END 
