C Version @(#)$Id$ $Source$ $Date$ 

      BLOCK DATA BD14
C
C  BLOCK DATA Subprogram 14: /MPDCOM/,/TPDCOM/
C
C  BLOCK DATA 14 initializes the average miles per day and trips per day values.
C
C  Common block array subscripts:
C
C  VALMPD(2,8)  -  VALMPD ( IG, IV )
C  MYGMPD(2,8)  -  MYGMPD ( IG, IV )
C  TPDCO(2,5)   -  TPDCO  ( IMDL, IV )
C
C  Common block dictionary:
C
C   Name   Type              Description
C  ------  ----  ------------------------------------------
C  /MPDCOM/:
C
C  MAXMPD   I    maximum number of MPD model year groups allowed
C  VALMPD   R    average miles per day (MPD)
C  MYGMPD   I    last year of MPD model year groups
C
C  /TPDCOM/:
C
C  TPDCO    R    the coefficients (intercept and slope) of the TPD equation
C
C  Notes:
C
C  None.
C
C
      COMMON /MPDCOM/ MAXMPD,VALMPD(2,8),MYGMPD(2,8)
      COMMON /TPDCOM/ TPDCO(2,8)
C
      DATA MAXMPD/2/
C
C  Continuation line code is IV.
C
C  Average miles per day
C
      DATA VALMPD/
     1    31.1,  0.0,
     2    26.3,  0.0,
     3    33.7,  0.0,
     4    36.7,  0.0,
     5    31.1,  0.0,
     6    29.8,  0.0,
     7   138.3,  0.0,
     8     8.3,  0.0/
C
C  Last year of MPD model year groups
C
      DATA MYGMPD/
     1   2020,0000,
     2   2020,0000,
     3   2020,0000,
     4   2020,0000,
     5   2020,0000,
     6   2020,0000,
     7   2020,0000,
     8   2020,0000/
C
C  /TPDCOM/:
C
C  The intercept and slope of the trips per day by model year equation
C
      DATA TPDCO/
     1   4.7187, 0.058508,
     2   4.7187, 0.058508,
     3   4.7187, 0.058508,
     4   6.88  , 0.0     ,
     5   3.05  , 0.0     ,
     6   3.05  , 0.0     ,
     7   6.65  , 0.0     ,
     8   1.35  , 0.0     /
C
      END
