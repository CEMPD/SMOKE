      BLOCK DATA BD02P5 
C 
C  BLOCK DATA Subprogram 02: /PSBLK/,/LEDCOM/,/TAMBLK/,/RLCOM2/, 
C                            /PSIZE/ 
C 
C  BLOCK DATA 02 initializes: 
C     1) The fraction of particles less than or equal to the 
C        particle size cutoff, PSZCUT, that is input by the user. 
C        This fraction depends on the source category. 
C     2) Tampering rates for misfueling and and catalyst removal 
C     3) The lead content of gasoline (in grams/gallon) by calendar year 
C 
C Changes: 
C 
C 12 Feb 01 AIR Task 11: Removed dust and some particle size variables. 
C           Changed last model year from 2020 to NEYR. 
C 04 Jan 01 AIR Task 08: Updated the k-factor for paved roads at PM2.5 
C           from 3.3 to 1.8, per Table 13.2-1.1. Updated the k-factors for 
C           unpaved roads per Table 13.2.2-2 and installed coefficients 
C           for the new equation (Equation 2 in document C13S02-1UNPAV.PDF). 
C           Corrected the documentation. Initialized MOIST to default 0.2%. 
C 14 Nov 00 AIR Task 06: Added initialization of COMMONs RLCOM2, PSIZE 
C           and PSDUST. Added NIYRPB to LEDCOM. 
C 10 Nov 00 AIR Task 06: Updated code to MOBILE6 standards. 
C 
C  Common block array subscripts: 
C 
C  PB(22,2)            - PB(IYPTR,LNL) 
C  IYRPB(22)           - IYRPB(IYPTR) 
C  PSDIST(5,2,8)       - PSDIST(IPTS,I,ISRCE) 
C  TAMRAT(2,3,3,2,2)   - TAMRAT(IMFLAG,ZML/DET1/2,IG,CAR/TRUCK,CAT/RMIS) 
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
C  PSDIST    R     For sources 1-6: The fraction of particles less than or  
C                  equal to the particle size cutoff input by the user, PSZCUT. 
C                  For sources 7-8: The fugitive dust k-factor associated with 
C                  the given particle size.                   
C                  PSDIST=PSDIST(ipts,I, ISRCE) 
C                  There are eight categories of sources: 
C 
C                  ISRCE = 1 fraction for leaded gasoline 
C                  ISRCE = 2 fraction for unleaded gasoline, Catalyst 
C                  ISRCE = 3 fraction for unleaded gasoline, no Catalyst 
C                  ISRCE = 4 fraction for diesel fuel 
C                  ISRCE = 5 fraction for break wear 
C                  ISRCE = 6 fraction for tire wear 
C                  ISRCE = 7 Unpaved Road Base emission factor (lb/VMT) 
C                  ISRCE = 8 Paved Road Base emission factor (g/VMT) 
C 
C                  ipts = the number of particle sizes for which there is 
C                         data about the fraction of total particles less 
C                         than or equal to those sizes (MXDATA=5) 
C                  I = 1  corresponds to the particle size 
C                  I = 2  corresponds to the fraction of particles less 
C                         than or equal to that size 
C 
C  /TAMBLK/:       Tampering common block 
C                  TAMRAT = TAMRAT(IMFLAG,ZML/DET1/2,IG,CAR/TRUCK,CAT/RMIS) 
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
C  /RLCOM2/: 
C 
C  IPBPTR    I     Index for fuel lead content by calendar year 
C 
C  /PSIZE/: 
C 
C  PSBRK     R    fraction of PM for break wear 
C  PSDIE     R    fraction of PM for diesel fuel 
C  PSL       R    fraction of PM for leaded gasoline 
C  PSNL      R    fraction of PM for unleaded gasoline, no Catalyst 
C  PSNLCT    R    fraction of PM for unleaded gasoline, Catalyst 
C  PSTIRE    R    fraction of PM for tire wear 
C                       
C  Notes: Modified Sept'93 to conform to new ap42. 
C  4-11-94 @ CSC-bsg Request 419 Change Tampering rates to be like MOBILE5a 
C  5-12-94 @ CSC-bsg Request 428 Change Lead content of fuel to updated #'s 
C  Jun-29-94 @ CSC-bsg Request 447 Investigate graphical trends in the 
C  particulate model by comparing IMFLAG=1 to IMFLAG=2 for various  
C  calendar years.  Discovered a mistake with the Tampering rates block data. 
C 
C 
      IMPLICIT NONE 
C 
      INCLUDE 'LEAD.I' 
      INCLUDE 'PARAME.I' 
      INCLUDE 'PSIZE.I' 
C 
C     Common Block PSBLK -- particle size, fractions, and k-factor distirbutions 
C 
      DATA PSDIST/ 
C     ISRCE = 1 = LEADED (fractions) 
     *     0.20, 2.00, 10.00, 0.00, 0.00, 
     *     0.23, 0.43,  0.64, 0.00, 0.00, 
C     ISRCE = 2 = NONLEADED, CATALYST (fractions) 
     *     0.20, 2.00, 10.00, 0.00, 0.00, 
     *     0.87, 0.89,  0.97, 0.00, 0.00, 
C     ISRCE = 3 = NONLEADED, NO CATALYST (fractions) 
     *     0.20, 2.00, 10.00, 0.00, 0.00, 
     *     0.42, 0.66,  0.90, 0.00, 0.00, 
C     ISRCE = 4 = DIESEL (fractions) 
     *     0.20, 1.00, 2.00, 2.50, 10.00, 
     *     0.73, 0.86, 0.90, 0.92,  1.00, 
C     ISRCE = 5 = BREAK (fractions) 
     *     0.43, 1.10, 4.70, 7.00, 10.00, 
     *     0.09, 0.16, 0.82, 0.90,  0.98, 
C     ISRCE = 6 = TIRE (fractions) 
     *     0.10, 10.00, 0.00, 0.00, 0.00, 
     *     0.01,  1.00, 0.00, 0.00, 0.00, 
C     ISRCE = 7 = UNPAVED ROAD DUST (k-factor, lb/VMT) 
     *     2.50, 10.00, 0.00, 0.00, 0.00, 
     *     0.38,  2.6 , 0.00, 0.00, 0.00, 
C     ISRCE = 8 = PAVED ROAD DUST (K-factor, g/VMT) 
     *     2.50, 10.00, 0.00, 0.00, 0.00, 
     *     1.8 ,  7.3 , 0.00, 0.00, 0.00/ 
C  
C  Apr-8-1994 @ CSC-bsg Request 419 made PART5 like MOBILE5a with respect 
C  to the Tampering rates of catalyst removal and misfueling. 
C  Jun-29-94 @ CSC-bsg Request 447 Investigate graphical trends in the 
C  particulate model by comparing IMFLAG=1 to IMFLAG=2 for various  
C  calendar years.  Discovered a mistake with the Tampering rates block data. 
C 
C    TAMBLK -- Tampering common block 
 
      DATA TAMRAT/ 
C CATS CARS 
     1  0.0084,  0.0067, 0.0193, 0.0084, 0.9731, 1.0, 
     2 -0.0019, -0.0152, 0.0038, 0.0060, 3.6154, 1.0, 
     3 -0.0008,  0.0,    0.0014, 0.0003, 1.04,   1.0, 
C CATS TRUCKS 
     1  0.0785,  0.1037, 0.0444, 0.0203, 0.981,  1.0, 
     2 -0.0018, -0.0062, 0.0258, 0.0109, 0.5767, 1.0, 
     3  0.0024,  0.0006, 0.0035, 0.0043, 4.0011, 5.3523, 
C RMIS CARS 
     1 -0.0259,  0.0907, 0.0342, 0.0051, 0.7653, 1.0, 
     2  0.0110,  0.0108, 0.0034, 0.0017, 2.6554, 0.08, 
     3  0.0,     0.0,    0.0025, 0.0036, 1.0,    1.0, 
C RMIS TRUCKS 
     1  0.1519,  0.0491, 0.0297, 0.0388, 1.6721, 1.0, 
     2  0.0113, -0.0379, 0.0203, 0.0229, 0.7318, 0.7812, 
     3  0.0,     0.0075, 0.0057, 0.0086, 1.0,    1.0/ 
C 
C  Lead content of gasoline (in grams/gallon) by calendar year 
C  Note: Under 40 CFR 80, CY1992+ lead content of leaded fuel shall be 
C        from 0.5 to 0.1 g/gal. 
C  1: leaded gasoline    2: unleaded gasoline 
C 
      DATA PB/ 
     1 2.07,1.82,2.02,2.03,1.76,1.76,1.33,1.01,1.02,0.83,0.84, 
     1 0.59,0.31,0.15,0.15,0.002, 0.0004, 0.0002, 0.0002,3*0.0, 
     2 .014,.014,.014,.014,.010,.016,.028,.009,.005,.003,.006, 
     2 .002,.002,.001,.001,.002,.0004,.0003,.0002, 0.0, 0.0, 0.0/ 
C 
      DATA IYRPB/ 
     *1974,1975,1976,1977,1978,1979,1980,1981,1982,1983,1984, 
     *1985,1986,1987,1988,1989,1990,1991,NEYR,   0,   0,   0/ 
C 
      DATA NIYRPB /22/ 
C 
C  COMMONS PSIZE, RLCOM2 
C 
      DATA PSL,PSNL,PSNLCT,PSDIE,PSBRK,PSTIRE /6*0.0/ 
      DATA IPBPTR /0/ 
C 
      END 
