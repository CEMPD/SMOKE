C Version "@(#)$Id$ $Source$ $Date$ 

      BLOCK DATA BD26
C
C  BLOCK DATA Subprogram 26: /HDDMAR/,/HDDREG/
C
C  BLOCK DATA 26 initializes the HDDV fleet characterization data by weight
C  subclass (C2B,L,M,H array name prefixes) & by relative model yr (JDX & KDX):
C  (1) annual mileage accumulation rates
C  (2) total new registrations
C
C  Common block array subscripts:
C
C  C2BDVM(25) - C2BDVM ( JDX )
C   LHDVM(25) -  LHDVM ( JDX )
C   MHDVM(25) -  MHDVM ( JDX )
C   HHDVM(25) -  HHDVM ( JDX )
C
C  C2BDRG(26) - C2BDRG ( KDX )
C   LHDRG(21) -  LHDRG ( KDX )
C   MHDRG(26) -  MHDRG ( KDX )
C   HHDRG(26) -  HHDRG ( KDX )
C    KDX=1 for pre-1981, 2 for 1981, ... , 21 for 2000+.
C
C  Common block dictionary:
C
C   Name   Type           Description
C  ------  ----  ---------------------------------------------
C  /HDDMAR/:
C
C  C2BDVM   R    class 2b HDDV mileage accumulation rates
C  LHDVM    R    light HDDV mileage accumulation rates
C  MHDVM    R    medium HDDV mileage accumulation rates
C  HHDVM    R    heavy HDDV mileage accumulation rates
C
C  /HDDREG/:
C
C  C2BDRG   R    class 2b HDDV new registrations
C  LHDRG    R    light HDDV new registrations
C  MHDRG    R    medium HDDV new registrations
C  HHDRG    R    heavy HDDV new registrations
C
C  Notes:
C
C  None.
C
C
      INCLUDE 'HDDMAR.I'
      INCLUDE 'HDDREG.I'
C
C  /HDDMAR/:
C
C  Annual average Class 2B truck mileage accumulation (divided by 100000)
C
      DATA C2BDVM/
     * .22517, .20009, .17779, .15798, .14038,
     * .12474, .11084, .09849, .08752, .07777,
     * .06910, .06140, .05456, .04848, .04308,
     * .03828, .03402, .03023, .02686, .02387,
     * .02121, .01884, .01675, .01488, .01322/
C
C  Annual average Light HDDV mileage accumulation (divided by 100000)
C
      DATA LHDVM/
     *.24239, .23273, .22346, .21455, .20600,
     *.19779, .18990, .18233, .17507, .16809,
     *.16139, .15496, .14878, .14285, .13716,
     *.13169, .12644, .12140, .11656, .11192,
     *.10746, .10317, .09906, .09511, .09132/
C
C  Annual average medium HDDV mileage accumulation (divided by 100000)
C
      DATA MHDVM/
     *.26081, .25204, .24357, .23538, .22746,
     *.21982, .21243, .20528, .19838, .19171,
     *.18527, .17904, .17302, .16720, .16158,
     *.15614, .15089, .14582, .14092, .13618,
     *.13160, .12718, .12290, .11877, .11478/
C
C  Annual average heavy HDDV mileage accumulation (divided by 100000)
C
      DATA HHDVM/
     *.62176, .58663, .55348, .52220, .49269,
     *.46485, .43858, .41380, .39042, .36836,
     *.34754, .32790, .30937, .29189, .27540,
     *.25983, .24515, .23130, .21823, .20590,
     *.19426, .18328, .17293, .16315, .15393/
C
C  /HDDREG/:
C
C  Annual Class 2B truck new registrations (1980- through 2000+)
C
      DATA C2BDRG/
     * .001, .001, .047, .119, .200, .309, .397, .495, .599, .710,
     * .831, .956,1.086,1.214,1.341,1.464,1.585,1.702,1.816,1.925,
     *2.030/
C
C  Annual light HDDV new registrations (1980- through 2000+)
C
      DATA LHDRG/
     * .003, .004, .009, .015, .024, .041, .057, .074, .091, .107,
     * .122, .136, .150, .164, .176, .188, .200, .211, .221, .232,
     * .241/
C
C  Annual medium HDDV new registrations (1980- through 2000+)
C
      DATA MHDRG/
     * .536, .570, .590, .609, .648, .711, .761, .811, .858, .902,
     * .944, .988,1.033,1.079,1.127,1.175,1.225,1.277,1.328,1.379,
     *1.430/
C
C  Annual heavy HDDV new registrations (1980- through 2000+)
C
      DATA HHDRG/
     * .990,1.033,1.049,1.067,1.129,1.189,1.228,1.260,1.287,1.308,
     *1.324,1.340,1.356,1.372,1.387,1.406,1.425,1.447,1.470,1.495,
     *1.522/
C
      END
