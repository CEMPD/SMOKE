C Version "@(#)$Id$ $Source$ $Date$ 

      BLOCK DATA BD19
C
C  BLOCK DATA Subprogram 19: /OPMOD2/,/OPMOD3/
C
C  BLOCK DATA 19 defines the last years of the model year groups used for
C  indexing BFRCOF and the mileage related parameters for the operating mode
C  component of the operating mode and temperature correction factor algorithm.
C
C  Common block array subscripts:
C
C  MYGBFR(8,3,6)  -  MYGBFR ( IG, IP, IVLM )
C
C  Common block dictionary:
C
C   Name   Type              Description
C  ------  ----  ---------------------------------------------
C  /OPMOD2/:
C
C  MYSTEP   I    flag = 2 if MY = 1981+ & IV = 1 = LDGV,
C                      or if MY = 1984+ & IV = 2 = LDGT1,
C                      or if MY = 1984+ & IV = 3 = LDGT2,  = 1 otherwise
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
C  IBFYR1(8,3) - IBFYR1( IG, IP )
C  IBFYR2(8,3) - IBFYR2( IG, IP )
C  IBFYR3(8,3) - IBFYR3( IG, IP )
C  IBFYR5(8,3) - IBFYR5( IG, IP )
C  IBFYR6(8,3) - IBFYR6( IG, IP )
C  IBFYR8(8,3) - IBFYR8( IG, IP )
C
C  Local array dictionary:
C
C   Name   Type                Description
C  ------  ----  --------------------------------------------------
C  IBFYR<IV> I   IV's BFRCOF myg upper bounds (ex: IBFYR8 = MC myg ub's),
C                equivalenced to the corresponding section of MYGBFR
C
C  Notes:
C
C  None.
C
C
      INCLUDE 'OPMOD2.I'
      INCLUDE 'OPMOD3.I'
      INCLUDE 'PARAME.I'
C
      DIMENSION
     *   IBFYR1(8,3),
     *   IBFYR2(8,3),
     *   IBFYR3(8,3),
     *   IBFYR5(8,3),
     *   IBFYR6(8,3),
     *   IBFYR8(8,3)
C
      EQUIVALENCE
     *   (MYGBFR(1,1,1),IBFYR1(1,1)),
     *   (MYGBFR(1,1,2),IBFYR2(1,1)),
     *   (MYGBFR(1,1,3),IBFYR3(1,1)),
     *   (MYGBFR(1,1,4),IBFYR5(1,1)),
     *   (MYGBFR(1,1,5),IBFYR6(1,1)),
     *   (MYGBFR(1,1,6),IBFYR8(1,1))
C
C  /OPMOD2/:
C
      DATA MAXBFR/8/,MYSTEP/1/
C
C  Last year of model year groups for normalized bag fractions.  Continuation
C  code is IP.
C
      DATA IBFYR1/
     H 1967,1969,1971,1974,1979,1980,2*0,
     C 1967,1969,1971,1974,1979,1980,2*0,
     N 1967,1972,1974,1976,1979,1980,2*0/
C
      DATA IBFYR2/
     H 1967,1969,1971,1974,1978,1980,1983,0000,
     C 1967,1969,1971,1974,1978,1980,1983,0000,
     N 1967,1972,1974,1978,1980,1983,2*0/
C
      DATA IBFYR3/
     H 1969,1973,1978,1980,1983,3*0,
     C 1969,1973,1978,1980,1983,3*0,
     N 1969,1973,1978,1980,1983,3*0/
C
      DATA IBFYR5/
     H 1974,1976,1977,1978,1979,NEYR,2*0,
     C 1974,1976,1977,1978,1979,NEYR,2*0,
     N 1974,1976,1977,1978,1979,1980,1982,NEYR/
C
      DATA IBFYR6/
     H 1978,1979,1982,NEYR,4*0,
     C 1978,1979,1982,NEYR,4*0,
     N 1978,1979,1984,NEYR,4*0/
C
      DATA IBFYR8/
     H 1977,1979,NEYR,5*0,
     C 1977,1979,NEYR,5*0,
     N 1977,NEYR,6*0/
C
C /OPMOD3/:
C
      DATA BFMILE/0.0/,KINKBF/1/,BFGT50/0.0/
C
      END
