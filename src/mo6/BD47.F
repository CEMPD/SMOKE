      BLOCK DATA BD47 
C 
C  Block Data 47 initializes the values used in the Complex Model and 
C  the Acrolein calculations. 
C 
C  Changes (last listed first) 
C 
C  14 Jan 02: AIR Task 15-22: Reformatted to get around LAHEY compiler bugs. 
C  24 Sep 01: AIR Task 15-22: Comments 
C  18 Jun 01: AIR Task 15-22: New MOBILE6 Block Data 
C       
C 
C  Common block array subscripts: 
C 
C      ATACET(7,11)     - ATACET() 
C      ATACETBASE(10)   - ATACETBASE(ITG) 
C      ATACROLEIN(28,2) - ATACROLEIN(IV,INH) 
C      ATBASE(3,13)     - ATBASE(IASTM,IFP) 
C      ATBENZ(6,11)     - ATBENZ(IFP,ITG) 
C      ATBENZBASE(10)   - ATBENZBASE(ITG) 
C      ATBUTA(6,11)     - ATBUTA(IFP,ITG) 
C      ATBUTABASE(10)   - ATBUTABASE(ITG) 
C      ATCENTER(14)     - ATCENTER(IFP) 
C      ATCENTERA(14)    - ATCENTERA(IFP) 
C      ATFORM(7,11)     - ATFORM(IFP,ITG) 
C      ATFORMBASE(10)   - ATFORMBASE(ITG) 
C      ATMTBE(4,10)     - ATMTBE(IFP,ITG) 
C      ATTRGT(13)       - ATTRGT(IMAP) 
C      ATVOC(10,11)     - ATVOC(IFP,ITG) 
C      ATVOCBASE(10)    - ATVOCBASE(ITG) 
C      CMHC(10)         - CMHC(ITG) 
C      VOLTOWT(4)       - VOLTOWT(IOXY) 
C 
C 
C  Common block dictionary: 
C 
C   Name       Type                  Description 
C  ------      ----  ----------------------------------------------------- 
C  /ATOX2/ 
C  ATACET       R    Complex model acetaldehyde fuel property coefficients 
C  ATACETBASE   R    Complex model acetaldehyde base values by technology group 
C  ATACROLEIN   R    Acrolein ratios by vehicle class and no cat/cat 
C  ATBASE       R    Complex model base fuel properties 
C  ATBENZ       R    Complex model benzene fuel property coefficients 
C  ATBENZBASE   R    Complex model benzene base values by technology group 
C  ATBUTA       R    Complex model 1,3butadiene fuel property coefficients 
C  ATBUTABASE   R    Complex model 1,3butadiene base values by technology group 
C  ATCENTER     R    Complex model fuel property centering values 
C  ATCENTERA    R    Complex model fuel property centering values for acetaldehyde only (may be a bug in the Complex Model) 
C  ATFORM       R    Complex model formaldehyde fuel property coefficients 
C  ATFORMBASE   R    Complex model formaldehyde base values by technology group 
C  ATMTBE       R    Complex model MTBE fuel property coefficients 
C  ATTRGT       R    Complex model target fuel property values 
C  ATVOC        R    Complex model VOC fuel property coefficients 
C  ATVOCBASE    R    Complex model VOC base values by technology group 
C  ATWTOXY      R    Air Toxic weighted oxygenate 
C  CMHC         R    Complex Model HC values 
C  VOLTOWT      R    Vol% to Wt% conversion factors  
C 
C 
      IMPLICIT NONE 
C 
      INCLUDE 'ATOX2.I' 
C 
      DATA ATCENTER/ !For all toxics except Acetaldehyde 
     *   1.774834,    !Oxygen, wt% 
     * 204.577894,    !Sulfur, ppm 
     *   8.6114785,   !RVP, psi 
     *  46.7257717,   !E200, % 
     *  85.898961979, !E300, % 
     *  28.2610891,   !Aromatic, vol% 
     *   7.3187162,   !Olefin, vol% 
     *   1.0666825,   !Benzene, vol% 
     *   0.9472399,   !MTBE, wt% 
     *   0.0232037,   !ETBE, wt% 
     *   0.3143523,   !Ethanol, wt% 
     *   0.0164433,   !TAME, wt% 
     *   0.0     0,   !Other 
     *  10.1731114/   !MTBE Fuel Concentration Centering 
C 
      DATA ATCENTERA/ !For Acetaldehyde only 
     *   1.2846497,   !Oxygen, wt% 
     * 202.0532096,   !Sulfur, ppm 
     *   8.6054497,   !RVP, psi 
     *  46.7151013,   !E200, % 
     *  85.8610736,   !E300, % 
     *  28.216531,    !Aromatic, vol% 
     *   7.260841,    !Olefin, vol% 
     *   1.0680169,   !Benzene, vol% 
     *   0.9472399,   !MTBE, wt% 
     *   0.0232037,   !ETBE, wt% 
     *   0.3143523,   !Ethanol, wt% 
     *   0.0164433,   !TAME, wt% 
     *   0.0     0,   !Other 
     *  10.1731114/   !MTBE Fuel Concentration Centering 
C 
      DATA ATBASE/ 
C       ASTM C  ASTM C  ASTM B 
C       SUMMER  WINTER    ALL 
     *   0.0 ,   0.0 ,   0.0 , !Oxygen, wt% 
     * 339.  , 338.  , 339.  , !Sulfur, ppm 
     *   8.7 ,   8.7 ,   8.7 , !RVP, psi 
     *  41.  ,  50.  ,  41.  , !E200, % 
     *  83.  ,  83.  ,  83.  , !E300, % 
     *  32.  ,  26.4 ,  32.  , !Aromatic, vol% 
     *   9.2 ,  11.9 ,   9.2 , !Olefins, vol% 
     *   1.53,   1.64,   1.53, !Benzene, vol% 
     *   0.0 ,   0.0 ,   0.0 , !MTBE, wt% 
     *   0.0 ,   0.0 ,   0.0 , !ETBE, wt% 
     *   0.0 ,   0.0 ,   0.0 , !Ethanol, wt% 
     *   0.0 ,   0.0 ,   0.0 , !TAME, wt% 
     *   0.0 ,   0.0 ,   0.0 / !Other 
C 
      DATA ATVOC/                                            !Tech 
C                                                            Group 
     * -0.003626,  0.000515,  0.035158, -0.003598, -0.008239, !All   Oxygen,    Sulfur,       RVP,      E200,      E300, 
     *  0.003552, -0.002858,  0.000106,  0.000408, -0.000287, !All Aromatic,   Olefins,    E200Sq,    E300Sq,  300*Arom, 
     *  0.0     ,  0.0     ,  0.0     , -0.003231,  0.0     , ! 1    Oxygen,    Sulfur,       RVP,      E200,      E300, 
     *  0.0     ,  0.0     ,  0.0     ,  0.0     ,  0.0     , ! 1  Aromatic,   Olefins,    E200Sq,    E300Sq,  300*Arom, 
     *  0.0     , -0.000241, -0.024777,  0.0     ,  0.0     , ! 2    Oxygen,    Sulfur,       RVP,      E200,      E300, 
     *  0.0     ,  0.0     ,  0.0     ,  0.0     ,  0.0     , ! 2  Aromatic,   Olefins,    E200Sq,    E300Sq,  300*Arom, 
     *  0.0     ,  0.0     ,  0.0     ,  0.0     ,  0.0     , ! 3    Oxygen,    Sulfur,       RVP,      E200,      E300, 
     *  0.0     ,  0.0     ,  0.0     ,  0.0     , -0.000332, ! 3  Aromatic,   Olefins,    E200Sq,    E300Sq,  300*Arom, 
     *  0.0     ,  0.000321,  0.0     ,  0.0     ,  0.0     , ! 4    Oxygen,    Sulfur,       RVP,      E200,      E300, 
     * -0.005971,  0.0     ,  0.0     ,  0.0     ,  0.0     , ! 4  Aromatic,   Olefins,    E200Sq,    E300Sq,  300*Arom, 
     *  0.0     ,  0.000361,  0.0     ,  0.0     ,  0.0     , ! 5    Oxygen,    Sulfur,       RVP,      E200,      E300, 
     *  0.0     ,  0.0     ,  0.0     ,  0.0     ,  0.0     , ! 5  Aromatic,   Olefins,    E200Sq,    E300Sq,  300*Arom, 
     *  0.0     , -0.000306,  0.0     ,  0.0     ,  0.0     , ! 6    Oxygen,    Sulfur,       RVP,      E200,      E300, 
     *  0.0     ,  0.0     ,  0.0     ,  0.0     ,  0.0     , ! 6  Aromatic,   Olefins,    E200Sq,    E300Sq,  300*Arom, 
     *  0.0     ,  0.0     ,  0.0     ,  0.0     ,  0.0     , ! 7    Oxygen,    Sulfur,       RVP,      E200,      E300, 
     *  0.0     ,  0.0     ,  0.0     ,  0.0     ,  0.0     , ! 7  Aromatic,   Olefins,    E200Sq,    E300Sq,  300*Arom, 
     *  0.0     ,  0.0     ,  0.0     ,  0.0     ,  0.0     , ! 8    Oxygen,    Sulfur,       RVP,      E200,      E300, 
     *  0.0     ,  0.0     ,  0.0     ,  0.0     ,  0.0     , ! 8  Aromatic,   Olefins,    E200Sq,    E300Sq,  300*Arom, 
     *  0.0     ,  0.0     ,  0.0     ,  0.0     ,  0.0     , ! 9    Oxygen,    Sulfur,       RVP,      E200,      E300, 
     *  0.0     ,  0.0     ,  0.0     ,  0.0     ,  0.0     , ! 9  Aromatic,   Olefins,    E200Sq,    E300Sq,  300*Arom, 
     *  0.0     , -0.000569,  0.008137,  0.0     ,  0.007892, !10    Oxygen,    Sulfur,       RVP,      E200,      E300, 
     *  0.0     ,  0.0     ,  0.0     ,  0.0     ,  0.0     / !10  Aromatic,   Olefins,    E200Sq,    E300Sq,  300*Arom 
C  
      DATA ATVOCBASE/ 
C          1      2      3      4      5 : Tech Group 
     * 493.0, 404.0, 408.0, 771.0, 317.0, 
C          6      7      8      9     10 : Tech Group 
     * 354.0, 689.0, 404.0, 457.0,3075.0/ 
C 
      DATA ATBENZ/                  !Tech 
C                                   Group 
     * 0.000337, 0.0     , 0.0     , !All   Sulfur,     E200,     E300, 
     * 0.025880, 0.222318, 0.0     , !All Aromatic,  Benzene,   Oxygen, 
     * 0.000717,-0.009481, 0.0     , ! 1    Sulfur,     E200,     E300, 
     * 0.0     , 0.0     , 0.0     , ! 1  Aromatic,  Benzene,   Oxygen, 
     * 0.0     , 0.0     , 0.0     , ! 2    Sulfur,     E200,     E300, 
     * 0.0     , 0.0     , 0.0     , ! 2  Aromatic,  Benzene,   Oxygen, 
     * 0.000850,-0.005776, 0.0     , ! 3    Sulfur,     E200,     E300, 
     * 0.0     , 0.0     , 0.0     , ! 3  Aromatic,  Benzene,   Oxygen, 
     * 0.0     , 0.0     , 0.0     , ! 4    Sulfur,     E200,     E300, 
     * 0.0     , 0.0     , 0.0     , ! 4  Aromatic,  Benzene,   Oxygen, 
     * 0.0     , 0.0     , 0.0     , ! 5    Sulfur,     E200,     E300, 
     * 0.022714, 0.0     , 0.0     , ! 5  Aromatic,  Benzene,   Oxygen, 
     * 0.0     , 0.0     , 0.0     , ! 6    Sulfur,     E200,     E300, 
     * 0.0     , 0.0     , 0.0     , ! 6  Aromatic,  Benzene,   Oxygen, 
     * 0.0     , 0.0     , 0.0     , ! 7    Sulfur,     E200,     E300, 
     * 0.0     , 0.0     , 0.0     , ! 7  Aromatic,  Benzene,   Oxygen, 
     * 0.0     , 0.0     , 0.0     , ! 8    Sulfur,     E200,     E300, 
     * 0.0     , 0.0     , 0.0     , ! 8  Aromatic,  Benzene,   Oxygen, 
     *-0.002283, 0.0     , 0.0     , ! 9    Sulfur,     E200,     E300, 
     * 0.0     , 0.0     , 0.0     , ! 9  Aromatic,  Benzene,   Oxygen, 
     * 0.0     , 0.0     , 0.011251, !10    Sulfur,     E200,     E300, 
     *-0.013998, 0.0     ,-0.096047/ !10  Aromatic,  Benzene,   Oxygen 
C  
      DATA ATBENZBASE/ 
C           1       2       3       4       5 : Tech Group 
     * 27.298, 22.388, 22.585, 42.715, 17.577, 
C           6       7       8       9      10 : Tech Group 
     * 19.639, 38.198, 22.388, 25.334,190.650/ 
C 
      DATA ATBUTA/                  !Tech 
C                                   Group 
     * 0.0     , 0.0     ,-0.007311, !All Oxygen,   Sulfur,     E200, 
     *-0.016785,-0.004005, 0.028238, !All   E300, Aromatic,  Olefins, 
     * 0.0     , 0.000506, 0.0     , ! 1  Oxygen,   Sulfur,     E200, 
     * 0.0     , 0.0     , 0.0     , ! 1    E300, Aromatic,  Olefins, 
     * 0.0     , 0.0     , 0.0     , ! 2  Oxygen,   Sulfur,     E200, 
     * 0.0     , 0.0     , 0.0     , ! 2    E300, Aromatic,  Olefins, 
     * 0.0     , 0.000544, 0.0     , ! 3  Oxygen,   Sulfur,     E200, 
     * 0.010535, 0.0     , 0.0     , ! 3    E300, Aromatic,  Olefins, 
     * 0.0     ,-0.000406, 0.0     , ! 4  Oxygen,   Sulfur,     E200, 
     * 0.0     , 0.0     , 0.0     , ! 4    E300, Aromatic,  Olefins, 
     * 0.0     , 0.0     , 0.0     , ! 5  Oxygen,   Sulfur,     E200, 
     * 0.0     , 0.0     , 0.0     , ! 5    E300, Aromatic,  Olefins, 
     * 0.0     , 0.0     , 0.013097, ! 6  Oxygen,   Sulfur,     E200, 
     * 0.0     , 0.0     , 0.0     , ! 6    E300, Aromatic,  Olefins, 
     * 0.0     , 0.0     , 0.0     , ! 7  Oxygen,   Sulfur,     E200, 
     * 0.0     , 0.0     , 0.0     , ! 7    E300, Aromatic,  Olefins, 
     * 0.0     , 0.0     , 0.0     , ! 8  Oxygen,   Sulfur,     E200, 
     * 0.0     , 0.0     , 0.0     , ! 8    E300, Aromatic,  Olefins, 
     * 0.0     , 0.0     , 0.0     , ! 9  Oxygen,   Sulfur,     E200, 
     * 0.0     , 0.0     , 0.0     , ! 9    E300, Aromatic,  Olefins, 
     *-0.060771, 0.0     , 0.0     , !10  Oxygen,   Sulfur,     E200, 
     * 0.008727, 0.0     , 0.015458/ !10    E300, Aromatic,  Olefins 
C 
      DATA ATBUTABASE/ 
C          1      2      3      4       5 : Tech Group 
     * 2.671, 2.190, 2.210, 4.179, 1.720 , 
C          6      7      8      9      10 : Tech Group 
     * 1.921, 3.737, 2.190, 2.479,43.9725/ 
C 
      DATA ATACET/                               !Tech 
C                                                Group 
     * 0.0002627, 0.0      ,-0.012157,-0.0055480, !All Sulfur,       RVP,     E300,  Aromatic, 
     * 0.0      , 0.3164665, 0.2493259,           !All   MTBE,      ETBE,     ETOH, 
     * 0.0      , 0.0      , 0.0      , 0.0     , ! 1  Sulfur,       RVP,     E300,  Aromatic, 
     *-0.0364640, 0.0      , 0.0      ,           ! 1    MTBE,      ETBE,     ETOH, 
     * 0.0      , 0.0      , 0.0      , 0.0     , ! 2  Sulfur,       RVP,     E300,  Aromatic, 
     * 0.0      , 0.0      , 0.0      ,           ! 2    MTBE,      ETBE,     ETOH, 
     * 0.0      , 0.0      , 0.0      , 0.0     , ! 3  Sulfur,       RVP,     E300,  Aromatic, 
     * 0.0      , 0.0      , 0.0      ,           ! 3    MTBE,      ETBE,     ETOH, 
     * 0.0      , 0.2423032, 0.0      , 0.0     , ! 4  Sulfur,       RVP,     E300,  Aromatic, 
     * 0.0      , 0.0      , 0.0      ,           ! 4    MTBE,      ETBE,     ETOH, 
     * 0.0      , 0.0      , 0.0      , 0.0     , ! 5  Sulfur,       RVP,     E300,  Aromatic, 
     * 0.0      , 0.0      , 0.0      ,           ! 5    MTBE,      ETBE,     ETOH, 
     * 0.0      , 0.0      , 0.0      , 0.0     , ! 6  Sulfur,       RVP,     E300,  Aromatic, 
     * 0.0      , 0.0      , 0.0      ,           ! 6    MTBE,      ETBE,     ETOH, 
     * 0.0      , 0.0      , 0.0      , 0.0     , ! 7  Sulfur,       RVP,     E300,  Aromatic, 
     * 0.0      , 0.0      , 0.0      ,           ! 7    MTBE,      ETBE,     ETOH, 
     * 0.0      , 0.0      , 0.0      , 0.0     , ! 8  Sulfur,       RVP,     E300,  Aromatic, 
     * 0.0      , 0.0      , 0.0      ,           ! 8    MTBE,      ETBE,     ETOH, 
     * 0.0      , 0.0      , 0.0      , 0.0     , ! 9  Sulfur,       RVP,     E300,  Aromatic, 
     * 0.0      , 0.0      , 0.0      ,           ! 9    MTBE,      ETBE,     ETOH, 
     * 0.0      , 0.0      , 0.0      , 0.0     , !10  Sulfur,       RVP,     E300,  Aromatic, 
     *-0.0559804, 0.0      , 0.0      /           !10    MTBE,      ETBE,     ETOH 
C 
      DATA ATACETBASE/ 
C          1      2      3      4      5 : Tech Group 
     * 2.424, 1.988, 2.006, 3.793, 1.561, 
C          6      7      8      9     10 : Tech Group 
     * 1.744, 3.392, 1.988, 2.250,15.006/ 
C 
      DATA ATFORM/                               !Tech 
C                                                Group 
     * 0.0      , 0.0      ,-0.010226,-0.0071660, !All  Sulfur,       RVP,     E300,  Aromatic, 
     * 0.0      , 0.0462131, 0.0      ,           !All Olefins,      MTBE,     ETOH, 
     * 0.0      , 0.0      , 0.0      , 0.0     , ! 1   Sulfur,       RVP,     E300,  Aromatic, 
     * 0.0      , 0.0      , 0.0      ,           ! 1  Olefins,      MTBE,     ETOH, 
     * 0.0      , 0.0      , 0.0      , 0.0     , ! 2   Sulfur,       RVP,     E300,  Aromatic, 
     * 0.0      , 0.0      , 0.0      ,           ! 2  Olefins,      MTBE,     ETOH, 
     * 0.0      , 0.0      , 0.0      , 0.0     , ! 3   Sulfur,       RVP,     E300,  Aromatic, 
     * 0.0      , 0.0      , 0.0      ,           ! 3  Olefins,      MTBE,     ETOH, 
     * 0.0      , 0.0      , 0.0      , 0.0     , ! 4   Sulfur,       RVP,     E300,  Aromatic, 
     * 0.0      , 0.0      , 0.0      ,           ! 4  Olefins,      MTBE,     ETOH, 
     * 0.0      , 0.0      , 0.0      , 0.0     , ! 5   Sulfur,       RVP,     E300,  Aromatic, 
     * 0.0      , 0.0      , 0.0      ,           ! 5  Olefins,      MTBE,     ETOH, 
     * 0.0      , 0.0      , 0.0      , 0.0     , ! 6   Sulfur,       RVP,     E300,  Aromatic, 
     * 0.0      , 0.0      , 0.0      ,           ! 6  Olefins,      MTBE,     ETOH, 
     * 0.0      , 0.0      , 0.0      , 0.0     , ! 7   Sulfur,       RVP,     E300,  Aromatic, 
     * 0.0      , 0.0      , 0.0      ,           ! 7  Olefins,      MTBE,     ETOH, 
     * 0.0      , 0.0      , 0.0      , 0.0     , ! 8   Sulfur,       RVP,     E300,  Aromatic, 
     * 0.0      , 0.0      , 0.0      ,           ! 8  Olefins,      MTBE,     ETOH, 
     * 0.0      , 0.0      , 0.0      , 0.0     , ! 9   Sulfur,       RVP,     E300,  Aromatic, 
     * 0.0      , 0.0      , 0.0      ,           ! 9  Olefins,      MTBE,     ETOH, 
     * 0.0      , 0.0      , 0.0      , 0.0     , !10   Sulfur,       RVP,     E300,  Aromatic, 
     *-0.0313520, 0.0      , 0.0      /           !10  Olefins,      MTBE,     ETOH 
C 
      DATA ATFORMBASE/ 
C          1      2      3      4       5 : Tech Group 
     * 5.913, 4.849, 4.892, 9.252, 3.807 , 
C          6      7      8      9      10 : Tech Group 
     * 4.254, 8.274, 4.849, 5.488,29.5815/ 
C 
      DATA ATMTBE/ 
C                                               Tech 
C      Intercept,  MTBEVP, MTBEVP^2,  MTBEVP^3, Group 
     * 1.468692, 0.256401,-0.045650, 0.0089880, ! 1 
     * 2.163555, 0.256401,-0.045650, 0.0      , ! 2 
     * 2.169000, 0.256401,-0.045650, 0.0      , ! 3 
     * 2.245980, 0.256401,-0.045650, 0.0      , ! 4 
     * 2.853995, 0.256401,-0.045650, 0.0      , ! 5 
     * 2.885609, 0.256401,-0.045650, 0.0      , ! 6 
     * 1.903058, 0.256401,-0.045650, 0.0      , ! 7 
     * 1.865054, 0.256401,-0.045650, 0.0      , ! 8 
     * 2.620530, 0.256401,-0.045650, 0.0      , ! 9 
     * 3.955618, 0.256401,-0.045650, 0.2175550/ !10 
C 
      DATA ATACROLEIN/ 
C 
C     No Catalyst 
C       LDGV  LDGT1  LDGT2  LDGT3  LDGT4 HDGV2B  HDGV3  
     * .0006, .0006, .0006, .0006, .0006, .0045, .0045, 
C      HDGV4  HDGV5  HDGV6  HDGV7 HDGV8A HDGV8B   LDDV  
     * .0045, .0045, .0045, .0045, .0045, .0045, .0035, 
C     LDDT12 HDDV2B  HDDV3  HDDV4  HDDV5  HDDV6  HDDV7  
     * .0035, .0035, .0035, .0035, .0035, .0035, .0035, 
C     HDDV8A HDDV8B     MC GASBUS  URBAN SCHOOL LDDT34 
     * .0035, .0035, .0006, .0045, .0035, .0035, .0035, 
C     With Catalyst 
C       LDGV  LDGT1  LDGT2  LDGT3  LDGT4 HDGV2B  HDGV3  
     * .0006, .0006, .0006, .0006, .0006, .0005, .0005, 
C      HDGV4  HDGV5  HDGV6  HDGV7 HDGV8A HDGV8B   LDDV  
     * .0005, .0005, .0005, .0005, .0005, .0005, .0035, 
C     LDDT12 HDDV2B  HDDV3  HDDV4  HDDV5  HDDV6  HDDV7  
     * .0035, .0035, .0035, .0035, .0035, .0035, .0035, 
C     HDDV8A HDDV8B     MC GASBUS  URBAN SCHOOL LDDT34 
     * .0035, .0035, .0006, .0005, .0035, .0035, .0035/ 
C 
      DATA VOLTOWT/                 ! Factors to change vol% to wt% 
C 
C        MTBE    ETBE    ETOH    TAME 
     * 0.1786, 0.1533, 0.3488, 0.1636/ 
C 
      DATA CMHC/ ! Complex Model HC (g/mi) for MTBE conversion to VOC 
C         1     2     3     4     5     6     7     8     9    10 Tech Group 
     * 0.23, 0.40, 0.36, 0.39, 0.00, 0.48, 0.34, 0.00, 0.00, 2.11/ 
C 
      END 
