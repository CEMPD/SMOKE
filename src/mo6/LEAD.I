C  LEAD.I 
C 
C Changes: 
C 
C 14 Nov 00 AIR Task 06: Added NIYRPB to LEDCOM. 
C 10 Nov 00 AIR Task 06: Updated code to MOBILE6 standards. 
C 
C  Common block dictionary: 
C 
C   Name    Type             Description 
C  ------   ----   ------------------------------------------------------- 
C  /LEDCOM/:       LEAD COMMON BLOCK 
C 
C  PB        R     Lead content of gasoline (in grams/gallon) by calendar year 
C                  PB = PB(IYPTR,LNL) 
C                        LNL=1: leaded gasoline 
C                        LNL=2: unleaded gasoline 
C                  Note: Under 40 CFR 80, CY1992+ lead content of leaded 
C                        fuel shall be from 0.5 to 0.1 g/gal. 
C  IYRPB     I     Model year groups for PB 
C                  IYRPB = IYRPB(IYPTR) 
C  NIYRPB    I     Number of cells in array IYRPB 
C 
C  /PSBLK/: 
C 
C  PSDIST    R     The fraction of particles less than or equal to the 
C                  particle size cutoff input by the user, PSZCUT. 
C                  PSDIST=PSDIST(ipts,I, ISRCE) 
C                  There are six categories of sources: 
C 
C                  ISRCE = 1 fraction for leaded gasoline 
C                  ISRCE = 2 fraction for unleaded gasoline, Catalyst 
C                  ISRCE = 3 fraction for unleaded gasoline, no Catalyst 
C                  ISRCE = 4 fraction for diesel fuel 
C                  ISRCE = 5 fraction for break wear 
C                  ISRCE = 6 fraction for tire wear 
C                  ISRCE = 7 fraction for unpaved road dust 
C                  ISRCE = 8 Paved Base emission factor (g/VMT) 
C 
C                  ipts = the number of particle sizes for which there is 
C                         data about the fraction of total particles less 
C                         than or equal to those sizes (MXDATA=5) 
C                  I = 1  corresponds to the particle size 
C                  I = 2  corresponds to the fraction of particles less 
C                         than or equal to that size 
C 
C  /TAMBLK/:       Tampering common block 
C 
C  TAMRAT    R     TAMRAT = TAMRAT(IMFLAG,ZML/DET1/2,IG,CAR/TRUCK,CAT/RMIS) 
C                  where: 
C                           IMFLAG = 1 for no IM program 
C                           IMFLAG = 2 for IM program 
C                           1=ZML 2=DET1 PRE 50,000 MILES 3=DET2 50-130K 
C                           IG=1 PRE 81 IG=2 1981-1983 IG=3 1984+ 
C                           1=CAR 2=TRUCK 
C                           1=CAT (Catalyst Tampering) vs 2=RMIS (misfueling) 
C                  TAMRAT is the fraction of vehicles that have 
C                         their catalysts removed if fifth dimension=1 and 
C                  TAMRAT is the fraction of vehicles that have 
C                         misfueled their vehicles if fifth dimension=2. 
C 
C  /RBCOM2/: 
C 
C  IPBPTR    I     Index for fuel lead content by calendar year 
C 
C 
C  The data are initialized in BD02P5. 
C 
C 
      INTEGER IPBPTR 
      INTEGER IYRPB 
      INTEGER NIYRPB 
      REAL    PB 
      REAL    PSDIST 
      REAL    TAMRAT 
C 
      COMMON /PSBLK/  PSDIST(5,2,8) 
      COMMON /LEDCOM/ PB(22,2), IYRPB(22), NIYRPB 
      COMMON /TAMBLK/ TAMRAT(2,3,3,2,2) 
      COMMON /RBCOM2/ IPBPTR 
C 
