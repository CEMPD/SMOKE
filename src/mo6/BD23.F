      BLOCK DATA BD23 
C 
C  BLOCK DATA Subprogram 23: /GSFCOM/ 
C 
C  BLOCK DATA 23 initializes the gasoline vs diesel sales fractions. 
C 
C  Changes: (Last change first) 
C 
C  26 Apr 00 @EPA-BAG changed diesel fractions for light duty diesel 
C                     trucks #3 to match diesel sales fractions of 
C                     light duty diesel trucks #4 
C  07 Apr 00 @EPA-djb Replaced arrays with final values and 
C            restructured the array sizes. 
C 
C  Common block array subscripts: 
C 
C  DSFRAC(25,14)   -  DSFRAC ( model year grouping , vehicle class ) 
C 
C  Common block dictionary: 
C 
C   Name   Type              Description 
C  ------  ----  ------------------------------------------------ 
C  /GSFCOM/: 
C 
C  MAXDSF   I    maximum number of sales fractions for any IV case 
C  DSFRAC   R    diesel sales as fraction of total model year sales 
C                (gasoline and diesel). 
C 
C  Local array subscripts: 
C 
C  DSF1(25,6) -  DSF ( model year grouping , vehicle class groups 1 - 6 ) 
C  DSF2(25,6) -  DSF ( model year grouping , vehicle class groups 7 - 12 ) 
C  DSF3(25,2) -  DSF ( model year grouping , vehicle class groups 13 & 14 ) 
C 
C  Model year groupings are 1 = Pre-1973 (1972 and earlier model years) 
C                      2 - 24 = 1973 through 1995 model years 
C                          25 = 1996 and later model years 
C 
C  Vehicle class groupings are: 
C 
C     1 : LDV 
C     2 : LDT1 
C     3 : LDT2 
C     4 : LDT3 
C     5 : LDT4 
C     6 : HDV2B 
C     7 : HDV3 
C     8 : HDV4 
C     9 : HDV5 
C    10 : HDV6 
C    11 : HDV7 
C    12 : HDV8A 
C    13 : HDV8B 
C    14 : HDBS (school bus) 
C 
C  Local array dictionary: 
C 
C   Name   Type              Description 
C  ------  ----  ----------------------------------------------- 
C  DSF1     R    diesel sales fractions vehicle class groups 1 - 6 
C  DSF2     R    diesel sales fractions vehicle class groups 7 - 12 
C  DSF3     R    diesel sales fractions vehicle class groups 13 & 14 
C 
C  Notes: 
C 
C 
      IMPLICIT NONE 
      INCLUDE 'DSFCOM.I' 
C 
      REAL DSF1(25,6),DSF2(25,6),DSF3(25,2) 
C 
      EQUIVALENCE 
     *   (DSFRAC(1, 1),DSF1(1,1)), 
     *   (DSFRAC(1, 7),DSF2(1,1)), 
     *   (DSFRAC(1,13),DSF3(1,1)) 
C 
      DATA MAXDSF/25/ 
C 
C  diesel sales as a fraction of total (g + d) model year sales 
C 
      DATA DSF1/ 
C  LDV 
     *0.0067,0.0067,0.0067,0.0155,0.0137,0.0093,0.0114,0.0269, 
     *0.0390,0.0706,0.0510,0.0241,0.0162,0.0097,0.0032,0.0027, 
     *0.0001,0.0004,0.0004,0.0013,0.0006,0.0003,0.0001,0.0006,0.0009, 
C  LDT1 
     *0.1170,0.1170,0.1170,0.1038,0.0187,0.0000,0.0259,0.0316, 
     *0.0439,0.0616,0.0656,0.0223,0.0120,0.0048,0.0033,0.0007, 
     *0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000, 
C  LDT2 
     *0.1170,0.1170,0.1170,0.1038,0.0187,0.0000,0.0259,0.0316, 
     *0.0439,0.0616,0.0656,0.0223,0.0120,0.0048,0.0033,0.0007, 
     *0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000, 
C  LDT3 
     *0.0001,0.0001,0.0001,0.0000,0.0000,0.0000,0.0001,0.0011, 
     *0.0006,0.0013,0.0256,0.0209,0.0169,0.0135,0.0124,0.0082, 
     *0.0072,0.0083,0.0096,0.0129,0.0115,0.0145,0.0111,0.0115,0.0126, 
C  LDT4 
     *0.0001,0.0001,0.0001,0.0000,0.0000,0.0000,0.0001,0.0011, 
     *0.0006,0.0013,0.0256,0.0209,0.0169,0.0135,0.0124,0.0082, 
     *0.0072,0.0083,0.0096,0.0129,0.0115,0.0145,0.0111,0.0115,0.0126, 
C  HDV2B 
     *0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000, 
     *0.0000,0.0138,0.2859,0.2918,0.3004,0.2743,0.2726,0.1958, 
     *0.1756,0.2058,0.2384,0.2963,0.2784,0.3263,0.2515,0.2578,0.1998/ 
C 
      DATA DSF2/ 
C  HDV3 
     *0.0000,0.0000,0.0248,0.0028,0.0010,0.0003,0.0001,0.0000, 
     *0.0000,0.0079,0.4277,0.5032,0.5139,0.6145,0.7842,0.7789, 
     *0.7488,0.7940,0.8477,0.8280,0.8068,0.8105,0.7910,0.7715,0.6774, 
C  HDV4 
     *0.0112,0.0090,0.0004,0.0078,0.0259,0.0000,0.0000,0.0000, 
     *0.0003,0.0414,0.0341,0.0738,0.1570,0.1968,0.2207,0.3178, 
     *0.5647,0.7158,0.7275,0.7316,0.7901,0.8331,0.8048,0.8473,0.8606, 
C  HDV5 
     *0.0016,0.0013,0.0001,0.0011,0.0037,0.0000,0.0000,0.0000, 
     *0.0000,0.0060,0.0049,0.0111,0.0255,0.0333,0.0383,0.0615, 
     *0.1543,0.2616,0.2730,0.2771,0.3462,0.4125,0.3670,0.4384,0.4647, 
C  HDV6 
     *0.0297,0.0274,0.0288,0.0365,0.0476,0.0808,0.1390,0.1679, 
     *0.3094,0.4413,0.3690,0.3569,0.4310,0.4525,0.4705,0.4734, 
     *0.4216,0.4537,0.5617,0.5788,0.5289,0.5767,0.5246,0.6078,0.6300, 
C  HDV7 
     *0.0966,0.0897,0.0940,0.1170,0.1489,0.2353,0.3610,0.4140, 
     *0.6107,0.7344,0.6717,0.6602,0.7261,0.7431,0.7567,0.7588, 
     *0.7184,0.7440,0.8177,0.8279,0.7972,0.8266,0.7943,0.8443,0.8563, 
C  HDV8A 
     *0.9720,0.9812,0.9819,0.9936,0.9920,0.9949,0.9964,0.9965, 
     *0.9974,0.9982,0.9978,0.9969,0.9976,0.9979,0.9980,0.9978, 
     *0.9969,0.9979,0.9982,0.9984,0.9977,0.9989,0.9987,0.9989,0.9992/ 
C 
      DATA DSF3/ 
C  HDV8B 
     *1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,1.0000, 
     *1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,1.0000, 
     *1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,1.0000, 
C  HDBS 
     *0.0000,0.0000,0.0000,0.0087,0.0086,0.0240,0.0291,0.0460, 
     *0.0594,0.2639,0.3260,0.3238,0.3845,0.5155,0.6733,0.7345, 
     *0.7502,0.7710,0.8760,0.9105,0.9900,0.8795,0.8525,0.8857,0.9585/ 
C 
      END 
 
