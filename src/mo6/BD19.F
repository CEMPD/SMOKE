      BLOCK DATA BD19 
C 
C  BLOCK DATA Subprogram 19: /OPMOD2/,/OPMOD3/ 
C 
C  BLOCK DATA 19 defines the last years of the model year groups used for 
C  indexing BFRCOF and the mileage related parameters for the operating mode 
C  component of the operating mode and temperature correction factor algorithm. 
C 
C  Changes: (Last change first) 
C 
C 
C  Common block array subscripts: 
C 
C  MYGBFR(8,3,8)  -  MYGBFR ( IG, IP, IVLM ) 
C 
C  Common block dictionary: 
C 
C   Name   Type              Description 
C  ------  ----  --------------------------------------------- 
C  /OPMOD2/: 
C 
C  MYSTEP   I    flag = 2 if MY = 1981+ & IV = 1   = LDGV, 
C                      or if MY = 1984+ & IV = 2,3 = LDGT12, 
C                      or if MY = 1984+ & IV = 4,5 = LDGT34, = 1 otherwise 
C  MAXBFR   I    maximum number of model year groups allowed for MYGBFR 
C  MYGBFR   I    last year of model year groups for BFRCOF 
C 
C  /OPMOD3/: 
C 
C  BFMILE   R    cumulative mileage in ICY for IV from MY scaled by 10**-4 
C  KINKBF   I    flag: 1 = mileage 50,000 or less, 2 = mileage over 50,000 
C  BFGT50   R    mileage greater than 50,000 (BFMILE=60,000 -> BFGT50=10,000) 
C 
C  Local array subscripts: 
C 
C  IBFY01(8,3) - IBFY01( IG, IP ) 
C  IBFY02(8,3) - IBFY02( IG, IP ) 
C  IBFY03(8,3) - IBFY03( IG, IP ) 
C  IBFY04(8,3) - IBFY04( IG, IP ) 
C  IBFY05(8,3) - IBFY05( IG, IP ) 
C  IBFY14(8,3) - IBFY14( IG, IP ) 
C  IBFY15(8,3) - IBFY15( IG, IP ) 
C  IBFY24(8,3) - IBFY24( IG, IP ) 
C 
C  Local array dictionary: 
C 
C   Name   Type                Description 
C  ------  ----  -------------------------------------------------- 
C  IBFY<IV>  I   IV's BFRCOF myg upper bounds (ex: IBFY24 = MC myg ub's), 
C                equivalenced to the corresponding section of MYGBFR 
C 
C  Notes: 
C 
C 20 Nov 96 DynTel@MLA 2-622  Expanded array MYGBFR to accomodate the new  
C           LDGT vehicle categories. 
C 
C 23 Jul 98 @DynTel-MLA 2-000  Changed LDGT names to conform to 
C           certification naming conventions. 
C 
C 
      IMPLICIT NONE 
      INCLUDE 'OPMOD2.I' 
      INCLUDE 'OPMOD3.I' 
      INCLUDE 'PARAME.I' 
C 
      INTEGER 
     *   IBFY01(8,3), 
     *   IBFY02(8,3), 
     *   IBFY03(8,3), 
     *   IBFY04(8,3), 
     *   IBFY05(8,3), 
     *   IBFY14(8,3), 
     *   IBFY15(8,3), 
     *   IBFY24(8,3) 
C 
      EQUIVALENCE 
     *   (MYGBFR(1,1,1),IBFY01(1,1)), 
     *   (MYGBFR(1,1,2),IBFY02(1,1)), 
     *   (MYGBFR(1,1,3),IBFY03(1,1)), 
     *   (MYGBFR(1,1,4),IBFY04(1,1)), 
     *   (MYGBFR(1,1,5),IBFY05(1,1)), 
     *   (MYGBFR(1,1,6),IBFY14(1,1)), 
     *   (MYGBFR(1,1,7),IBFY15(1,1)), 
     *   (MYGBFR(1,1,8),IBFY24(1,1)) 
C 
C  /OPMOD2/: 
C 
      DATA MAXBFR/8/,MYSTEP/1/ 
C 
C  Last year of model year groups for normalized bag fractions.  Continuation 
C  code is IP. 
C  LDGV 
      DATA IBFY01/ 
     H 1967,1969,1971,1974,1979,1980,2*0, 
     C 1967,1969,1971,1974,1979,1980,2*0, 
     N 1967,1972,1974,1976,1979,1980,2*0/ 
C  LDGT1 
      DATA IBFY02/ 
     H 1967,1969,1971,1974,1978,1980,1983,0000, 
     C 1967,1969,1971,1974,1978,1980,1983,0000, 
     N 1967,1972,1974,1978,1980,1983,2*0/ 
C  LDGT2 
      DATA IBFY03/ 
     H 1967,1969,1971,1974,1978,1980,1983,0000, 
     C 1967,1969,1971,1974,1978,1980,1983,0000, 
     N 1967,1972,1974,1978,1980,1983,2*0/ 
C  LDGT3 
      DATA IBFY04/ 
     H 1969,1973,1978,1980,1983,3*0, 
     C 1969,1973,1978,1980,1983,3*0, 
     N 1969,1973,1978,1980,1983,3*0/ 
C  LDGT4 
      DATA IBFY05/ 
     H 1969,1973,1978,1980,1983,3*0, 
     C 1969,1973,1978,1980,1983,3*0, 
     N 1969,1973,1978,1980,1983,3*0/ 
C  LDDV 
      DATA IBFY14/ 
     H 1974,1976,1977,1978,1979,NEYR,2*0, 
     C 1974,1976,1977,1978,1979,NEYR,2*0, 
     N 1974,1976,1977,1978,1979,1980,1982,NEYR/ 
C  LDDT 
      DATA IBFY15/ 
     H 1978,1979,1982,NEYR,4*0, 
     C 1978,1979,1982,NEYR,4*0, 
     N 1978,1979,1984,NEYR,4*0/ 
C  MC 
      DATA IBFY24/ 
     H 1977,1979,NEYR,5*0, 
     C 1977,1979,NEYR,5*0, 
     N 1977,NEYR,6*0/ 
C 
C /OPMOD3/: 
C 
      DATA BFMILE/0.0/,KINKBF/1/,BFGT50/0.0/ 
C 
      END 
