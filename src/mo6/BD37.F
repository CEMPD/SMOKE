      BLOCK DATA BD37 
C 
C  BLOCK DATA Subprogram 37: /FID/ 
C 
C  BLOCK DATA 37 initializes VOC & TOG FID model yr groups & correction factors. 
C 
C  Changes: (Last change first) 
C 
C  29 MAR 00 @EPA-elg   Fixed typo bug in TOGGAS Block Data. 
C                       The incorrect TOG effect was being applied to 
C                       the 1987 model year. 
C 
C  Common block array subscripts: 
C 
C  MYFID(16, 7)            -  MYFID ( IG, IVFID ) 
C  TOGCF(16, 7, 3)         -  TOGCF (  IGFID, IVFID, IFUEL  ) 
C  VOCCF(16, 7, 3)         -  VOCCF (  IGFID, IVFID, IFUEL  ) 
C 
C  Common block dictionary: 
C 
C   Name   Type               Description 
C  ------  ----  -------------------------------------------- 
C  MAXFID   I    Number of VOC and TOG FID model year groups. 
C  MYFID    I    VOC and TOG FID model year group years. 
C  TOGCF    R    TOG FID correction factors by my, vehicle, and fuel type. 
C  VOCCF    R    VOC FID correction factors by my, vehicle, and fuel type. 
C 
C  Local array subscripts: 
C 
C  TOGGAS(16,7)  -  TOGGAS ( IGFID, IVFID ) 
C  TOGETH(16,7)  -  TOGGAS ( IGFID, IVFID ) 
C  TOGHOL(16,7)  -  TOGGAS ( IGFID, IVFID ) 
C  VOCGAS(16,7)  -  VOCGAS ( IGFID, IVFID ) 
C  VOCETH(16,7)  -  VOCGAS ( IGFID, IVFID ) 
C  VOCHOL(16,7)  -  VOCGAS ( IGFID, IVFID ) 
C 
C  Local array dictionary: 
C 
C   Name   Type              Description 
C  ------  ----  ----------------------------------------------- 
C  TOGGAS   R    TOG FID correction factor for gasoline fueled vehicles 
C  TOGETH   R    TOG FID correction factor for ether blend fueled vehicles 
C  TOGHOL   R    TOG FID correction factor for alcohol blend fueled vehicles 
C  VOCGAS   R    VOC FID correction factor for gasoline fueled vehicles 
C  VOCETH   R    VOC FID correction factor for ether blend fueled vehicles 
C  VOCHOL   R    VOC FID correction factor for alcohol blend fueled vehicles 
C 
C  Notes: 
C 
C  BD37 was added for MOBILE4.1. 
C 
C 11 Jun 97 DynTel@MLA 2-622  Expanded the arrays to accomodate the new 
C           GAS and Diesel BUS vehicle categories. 
C 
      IMPLICIT NONE 
      INCLUDE 'FID.I' 
      INCLUDE 'PARAME.I' 
C 
      REAL TOGGAS(16,7),TOGETH(16,7),TOGHOL(16,7) 
      REAL VOCGAS(16,7),VOCETH(16,7),VOCHOL(16,7) 
C 
      EQUIVALENCE 
     *     (TOGCF(1,1,1),TOGGAS(1,1)), 
     *     (TOGCF(1,1,2),TOGETH(1,1)), 
     *     (TOGCF(1,1,3),TOGHOL(1,1)), 
     *     (VOCCF(1,1,1),VOCGAS(1,1)), 
     *     (VOCCF(1,1,2),VOCETH(1,1)), 
     *     (VOCCF(1,1,3),VOCHOL(1,1)) 
C 
C  Number of FID my groups 
C 
      DATA MAXFID/16/ 
C 
C  TOG and VOC FID model year groups: 
C  (Continuation character represents: 
C   1=LDGV, 2=LDGT, 3=HDGV, 4=LDDV/T, 5=HDDV, 6=GAS BUS, 7=DIESEL BUS) 
C 
      DATA MYFID/ 
     1    1974,   1975,   1977,   1979,   1980,   1981,   1982,   1983, 
     1    1984,   1985,   1986,   1987,   1988,   1989,   1990,   NEYR, 
     2    1974,   1975,   1977,   1979,   1980,   1981,   1982,   1983, 
     2    1984,   1985,   1986,   1987,   1988,   1989,   1990,   NEYR, 
     3    1982,   1986,   1987,   1989,   NEYR,   11*0, 
     4    NEYR,   15*0, 
     5    NEYR,   15*0, 
     6    1982,   1986,   1987,   1989,   NEYR,   11*0, 
     7    NEYR,   15*0/ 
C 
C  TOG FID correction factors: 
C  (Continuation character represents: 
C   1=LDGV, 2=LDGT, 3=HDGV, 4=LDDV/T, 5=HDDV, 6=GAS BUS, 7=DIESEL BUS) 
C 
      DATA TOGGAS/ 
C  Gasoline 
     1  1.0352, 1.0226, 1.0219, 1.0211, 1.0206, 1.0211, 1.0211, 1.0218, 
     1  1.0213, 1.0202, 1.0196, 1.0190, 1.0175, 1.0172, 1.0167, 1.0163, 
     2  1.0352, 1.0226, 1.0219, 1.0211, 1.0206, 1.0200, 1.0204, 1.0202, 
     2  1.0204, 1.0197, 1.0178, 1.0179, 1.0171, 1.0176, 1.0181, 1.0176, 
     3  1.0478, 1.0465, 1.0399, 1.0372, 1.0332, 11*0.0, 
     4  1.0516, 15*0.0, 
     5  1.0349, 15*0.0, 
     6  1.0478, 1.0465, 1.0399, 1.0372, 1.0332, 11*0.0, 
     7  1.0349, 15*0.0/ 
      DATA TOGETH/ 
C  Ether Blend 
     1  1.0352, 1.0226, 1.0219, 1.0211, 1.0206, 1.0211, 1.0211, 1.0218, 
     1  1.0213, 1.0202, 1.0196, 1.1090, 1.0175, 1.0172, 1.0167, 1.0163, 
     2  1.0352, 1.0226, 1.0219, 1.0211, 1.0206, 1.0200, 1.0204, 1.0202, 
     2  1.0204, 1.0197, 1.0178, 1.0179, 1.0171, 1.0176, 1.0181, 1.0176, 
     3  1.0478, 1.0465, 1.0399, 1.0372, 1.0332, 11*0.0, 
     4  1.0516, 15*0.0, 
     5  1.0349, 15*0.0, 
     6  1.0478, 1.0465, 1.0399, 1.0372, 1.0332, 11*0.0, 
     7  1.0349, 15*0.0/ 
      DATA TOGHOL/ 
C  Alcohol Blend 
     1  1.0352, 1.0226, 1.0219, 1.0211, 1.0206, 1.0211, 1.0211, 1.0218, 
     1  1.0213, 1.0202, 1.0196, 1.1090, 1.0175, 1.0172, 1.0167, 1.0163, 
     2  1.0352, 1.0226, 1.0219, 1.0211, 1.0206, 1.0200, 1.0204, 1.0202, 
     2  1.0204, 1.0197, 1.0178, 1.0179, 1.0171, 1.0176, 1.0181, 1.0176, 
     3  1.0478, 1.0465, 1.0399, 1.0372, 1.0332, 11*0.0, 
     4  1.0516, 15*0.0, 
     5  1.0349, 15*0.0, 
     6  1.0478, 1.0465, 1.0399, 1.0372, 1.0332, 11*0.0, 
     7  1.0349, 15*0.0/ 
C 
C  VOC FID correction factors: 
C  (Continuation character represents: 
C   1=LDGV, 2=LDGT, 3=HDGV, 4=LDDV/T, 5=HDDV, 6=GAS BUS, 7=DIESEL BUS) 
C 
      DATA VOCGAS/ 
C  Gasoline 
     1  1.0239, 0.9975, 0.9959, 0.9942, 0.9905, 0.9710, 0.9709, 0.9690, 
     1  0.9680, 0.9705, 0.9719, 0.9734, 0.9768, 0.9777, 0.9787, 0.9797, 
     2  1.0239, 0.9975, 0.9959, 0.9942, 0.9905, 0.9913, 0.9915, 0.9886, 
     2  0.9833, 0.9829, 0.9821, 0.9770, 0.9748, 0.9736, 0.9757, 0.9766, 
     3  1.0312, 1.0281, 1.0126, 1.0064, 0.9971, 11*0.0, 
     4  1.0497, 15*0.0, 
     5  1.0332, 15*0.0, 
     6  1.0312, 1.0281, 1.0126, 1.0064, 0.9971, 11*0.0, 
     7  1.0332, 15*0.0/ 
      DATA VOCETH/ 
C  Ether Blend 
     1  1.0239, 0.9975, 0.9959, 0.9942, 0.9905, 0.9710, 0.9709, 0.9690, 
     1  0.9680, 0.9705, 0.9719, 0.9734, 0.9768, 0.9777, 0.9787, 0.9797, 
     2  1.0239, 0.9975, 0.9959, 0.9942, 0.9905, 0.9913, 0.9915, 0.9886, 
     2  0.9833, 0.9829, 0.9821, 0.9770, 0.9748, 0.9736, 0.9757, 0.9766, 
     3  1.0312, 1.0281, 1.0126, 1.0064, 0.9971, 11*0.0, 
     4  1.0497, 15*0.0, 
     5  1.0332, 15*0.0, 
     6  1.0312, 1.0281, 1.0126, 1.0064, 0.9971, 11*0.0, 
     7  1.0332, 15*0.0/ 
      DATA VOCHOL/ 
C  Alcohol Blend 
     1  1.0239, 0.9975, 0.9959, 0.9942, 0.9905, 0.9710, 0.9709, 0.9690, 
     1  0.9680, 0.9705, 0.9719, 0.9734, 0.9768, 0.9777, 0.9787, 0.9797, 
     2  1.0239, 0.9975, 0.9959, 0.9942, 0.9905, 0.9913, 0.9915, 0.9886, 
     2  0.9833, 0.9829, 0.9821, 0.9770, 0.9748, 0.9736, 0.9757, 0.9766, 
     3  1.0312, 1.0281, 1.0126, 1.0064, 0.9971, 11*0.0, 
     4  1.0497, 15*0.0, 
     5  1.0332, 15*0.0, 
     6  1.0312, 1.0281, 1.0126, 1.0064, 0.9971, 11*0.0, 
     7  1.0332, 15*0.0/ 
C 
      END 
