C Version @(#)$Id$ $Source$ $Date$ 

      BLOCK DATA BD15
C
C  BLOCK DATA Subprogram 15: /CRACOM/
C
C  BLOCK DATA 15 initializes the evaporative crankcase HC emissions.
C
C  Common block array subscripts:
C
C  CRAVAL(3,8,2)  -  CRAVAL ( IG, IV, IR )
C  MYGCRA(3,8,2)  -  MYGCRA ( IG, IV, IR )
C
C  Common block dictionary:
C
C   Name   Type              Description
C  ------  ----  -------------------------------------------
C  /CRACOM/:
C
C  MAXCRA   I    maximum number of crankcase HC model year groups allowed
C  CRAVAL   R    crankcase HC emissions
C  MYGCRA   I    last year of crankcase HC model year groups
C
C  Local array subscripts:
C
C  CRAVLO(3,8) - CRAVLO ( IG, IV )
C  CRAVHI(3,8) - CRAVHI ( IG, IV )
C  ICRALO(3,8) - ICRALO ( IG, IV )
C  ICRAHI(3,8) - ICRAHI ( IG, IV )
C
C  Local array dictionary:
C
C   Name   Type                Description
C  ------  ----  -----------------------------------------
C  CRAVLO   R    crankcase HC emissions in high altitude region
C  CRAVHI   R    crankcase HC emissions in high altitude region
C  ICRALO   I    last year of crankcase HC model year groups
C  ICRAHI   I    last year of crankcase HC model year groups
C
C  Notes:
C
C  None.
C
C
      COMMON /CRACOM/ MAXCRA,CRAVAL(3,8,2),MYGCRA(3,8,2)
C
      DIMENSION
     *   CRAVLO(3,8),
     *   CRAVHI(3,8),
     *   ICRALO(3,8),
     *   ICRAHI(3,8)
C
      EQUIVALENCE
     *   (CRAVAL(1,1,1),CRAVLO(1,1)),
     *   (CRAVAL(1,1,2),CRAVHI(1,1)),
     *   (MYGCRA(1,1,1),ICRALO(1,1)),
     *   (MYGCRA(1,1,2),ICRAHI(1,1))
C
      DATA MAXCRA/3/
C
C  Value of crankcase HC emissions
C
C  Low altitude region
C
      DATA CRAVLO/
     *    4.10,  0.80,  0.00,
     *    4.10,  0.80,  0.00,
     *    5.70,  0.00,  0.00,
     *    5.70,  0.00,  0.00,
     *    0.00,  0.00,  0.00,
     *    0.00,  0.00,  0.00,
     *    0.00,  0.00,  0.00,
     *    0.31,  0.00,  0.00/
C
C  High altitude region
C
      DATA CRAVHI/
     *    5.29,  1.03,  0.00,
     *    5.29,  1.03,  0.00,
     *    7.35,  0.00,  0.00,
     *    7.35,  0.00,  0.00,
     *    0.00,  0.00,  0.00,
     *    0.00,  0.00,  0.00,
     *    0.00,  0.00,  0.00,
     *    0.40,  0.00,  0.00/
C
C  Last year of crankcase HC model year groups
C
C  Low altitude region
C
      DATA ICRALO/
     *   1962,1967,2020,
     *   1962,1967,2020,
     *   1967,2020,0000,
     *   1967,2020,0000,
     *   2020,0000,0000,
     *   2020,0000,0000,
     *   2020,0000,0000,
     *   1977,2020,0000/
C
C  High altitude region
C
      DATA ICRAHI/
     *   1962,1967,2020,
     *   1962,1967,2020,
     *   1967,2020,0000,
     *   1967,2020,0000,
     *   2020,0000,0000,
     *   2020,0000,0000,
     *   2020,0000,0000,
     *   1977,2020,0000/
C
      END
