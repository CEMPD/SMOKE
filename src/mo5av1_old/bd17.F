C Version @(#)$Id$ $Source$ $Date$ 

      BLOCK DATA BD17
C
C  BLOCK DATA Subprogram 17: /TEMPC1/,/TEMPC2/,TEMPC3/,/TEMPC4/,/TEMPC5/
C                            /TEMPC6/,/TEMPC7/,/TEMPS/
C
C  BLOCK DATA 17 initializes the temperature correction factor (tcf)
C  coefficents and their model year groups (myg), the computed tcf parameters,
C  the tcf vector and the cold start CO additive tcf offset parameters.
C  The tcf coefficients and offsets have cases with the Fuel Delivery Systems
C  (FDS) combined and with them differentiated.
C  BD 17 also initializes the following temperatures: ambient, minimum,
C  maximum, exhaust, and evaporative
C
C  In this cold CO standard version, BD 17 also initializes all cold CO
C  common block parameters and data.
C
C  Common block array subscripts:
C
C  CHFRAC(3)      -  CHFRAC ( IB )
C  CMUL(3,2)      -  CMUL   ( IVLDG, ICOLD )
C  CCSTD(3,2)     -  CCSTD  ( IVLDG, ICOLD )
C  COOTFC(3)      -  COOTFC ( IGO )
C  COOTFD(3)      -  COOTFD ( IFDS )
C  CSUB(3,2)      -  CSUB   ( IVLDG, ICOLD )
C  MDLOHI(3)      -  MDLOHI ( IP )
C  MYCOL(2)       -  MYCOL  ( ICOLD )
C  MYGCOO(3)      -  MYGCOO ( IVLDG )
C  MYGTFD(3)      -  MYGTFD ( IVLDG )
C  MYGTCF(9,3,5)  -  MYGTCF ( IG, IP, IVTEMP )
C  PIVPCT(3,3,2)  -  PIVPCT ( MY-MYCOL+1, IVLDG, ICOLD )
C  TCF(3)         -  TCF    ( IB )
C  TDIFF(3)       -  TDIFF ( IP )
C  TEMEXH(3)      -  TEMEXH ( IP )
C  TEMEVP(6)      -  TEMEVP ( IHRD )
C  TTFC(2,3,9,3)  -  TTFC ( MDLOHI, IB, IT, IP )
C  TTFD(3,3,3)    -  TTFD ( IB, IP, IFDS )
C  TT4(2,6,3)     -  TT4 ( MDLOHI, IT, IP )
C  TT8(2,3,4,3)   -  TT8 ( MDLOHI, IB, IT, IP )
C
C  Common block dictionary:
C
C   Name   Type              Description
C  ------  ----  ---------------------------------------------
C  /TEMPC1/:
C
C  TTFC     R    low/high tcf coefficients for IVLDG with FDSs Combined
C  TTFD     R    low tcf's for IVLDG with FDSs Differentiated
C  TT4      R    low/high tcf coefficients for HDGV
C  TT8      R    low/high tcf coefficients for MC
C
C  /TEMPC2/:
C
C  MAXTCF   I    maximum number of tcf model year groups allowed
C  MYGTCF   I    the last year of the combined FDS IVTEMP tcf model year groups
C  MYGTFD   I    the first year of differentiating IVLDG tcf by FDS
C  ISHIFT   I    IVLDG MYG shift to change array indexed from TTFC to COOTFC
C
C  /TEMPC3/:
C
C  MDLOHI   I    TCF case: 0 = mid (no correction), 1 = low, 2 = high
C  TDIFF    R    difference between entered ambient or calculated
C                exhaust temperature & FTP temp
C  CHFRAC   R    computed operating mode fractions for given MY & IV
C  TCF      R    temperature correction factor for each bag
C
C  /TEMPC4/:
C
C  MYGCOO   I    first my for using the CO additive offset, if IVLDG&Bag1&LowT
C  LOWCO    I    switch flags whether or not case of CO additive offset for
C                IVLDG cold start/low temperature condition: 1 = no & 2 = yes
C  COOTFC   R    low temp CO offsets (1980/1/2) for IVLDG with FDS Combined
C  COOTFD   R    low temp CO offset for IVLDG with FDS Differentiated
C
C  /TEMPC5/:
C
C  NCC      I    number of cold CO standards, 1 or 2 (2 = read in a second
C                (MYCOL,CCSTD) pair, with MYCOL(1).LT.MYCOL(2) )
C  MYCOL    I    first MY for which cold CO standard applies (ex: 1992)
C  CCSTD    R    cold CO emissions standard(s); must be met by LDGV/T produced
C                in corresponding MYCOL or later (ex: 3.4 for 2000+)
C
C  /TEMPC6/:
C
C  CSUB     R    correction to be subtracted from the additive low temp bag 1
C                CO offset to adjust for a "lenient" strategy cold CO standard.
C  CMUL     R    correction to be multiplied times both the additive low temp
C                bag 1 CO offset & the multiplicative low temp bags 2 & 3 cfs,
C                to adjust for a "stringent" cold CO standard.
C  PIVPCT   R    percentage of IVLDG classes covered by cold CO standard ICOLD
C                in model years MYCOL, MYCOL+1, and MYCOL+2 and later
C
C  /TEMPC7/:
C
C  ICOLD    I    defined in Subscript Dictionary: here indexes CSUB/MUL.
C  ADDCSW   R    switch (1 = off & 2 = on) for subtracting CSUB correction
C                from the additive low temperature bag 1 CO offset
C  MULCSW   R    switch (1 = off & 2 = on) for multiplying CMUL cf times both
C                bag 1 low CO offset & bags 2 & 3 low CO (cf - 1.0)
C
C  /TEMPS/:
C
C  AMBT     R    user supplied ambient temperature
C  TEMMIN   R    user supplied SCNAME's minimum ambient temperature,
C                used as fuel tank minimum temperature in diurnal index,
C                exhaust temperature, and evaporative temperature
C                calculations
C  TEMMAX   R    user supplied SCNAME's maximum ambient temperature,
C                used as fuel tank maximum temperature in diurnal index,
C                exhaust temperature, and evaporative temperature
C                calculations
C  TEMEXH   R    exhaust temperatures for HC, CO, and NOx
C  TEMEVP   R    evaporative temperatures for diurnal, hot soak,
C                and running loss
C
C  Local array subscripts:
C
C  TFCHC(2,3,9)   -  TFCHC ( MDLOHI, IB, IT )
C  TFCCO(2,3,9)   -  TFCCO ( MDLOHI, IB, IT )
C  TFCNOX(2,3,9)  -  TFCNOX ( MDLOHI, IB, IT )
C
C  Local array dictionary:
C
C   Name   Type                Description
C  ------  ----  --------------------------------------------------
C  For the local arrays equivalenced to TTFC:
C
C  TFC<IP>  R    low/high tcf coefficients for pollutant IP by bag and model
C                year group for vehicle classes IVLDG, with FDSs Combined.
C
C  Notes:
C
C  BD 17 was modified for MOBILE4.1 to include COLD CO vlaues.
C  BD 17 was modified for MOBILE4.1 to remove TCFF (oxy fuels support.)
C
      INTEGER ADDCSW
C
      COMMON /TEMPC1/TTFC(2,3,9,3),TTFD(3,3,3),TT4(2,6,3),TT8(2,3,4,3)
      COMMON /TEMPC2/ MAXTCF,MYGTCF(9,3,5),MYGTFD(3),ISHIFT
      COMMON /TEMPC3/ MDLOHI(3),TDIFF(3),CHFRAC(3),TCF(3)
      COMMON /TEMPC4/ MYGCOO(3),LOWCO,COOTFC(3),COOTFD(3)
      COMMON /TEMPC5/ NCC,MYCOL(2),CCSTD(3,2)
      COMMON /TEMPC6/ CSUB(3,2),CMUL(3,2),PIVPCT(3,3,2)
      COMMON /TEMPC7/ ICOLD,ADDCSW,MULCSW
      COMMON /TEMPS/ AMBT,TEMMIN,TEMMAX,TEMEXH(3),TEMEVP(6)
C
      DIMENSION TFCHC(2,3,9),TFCCO(2,3,9),TFCNOX(2,3,9)
C
      EQUIVALENCE
     *   (TTFC(1,1,1,1),TFCHC(1,1,1)),
     *   (TTFC(1,1,1,2),TFCCO(1,1,1)),
     *   (TTFC(1,1,1,3),TFCNOX(1,1,1))
C
C  LDGV, LDGT1 and LDGT2 (IVLDG) carbureted
C
C            Bag 1                  Bag 2                 Bag 3
C  AmbT:  Low        High       Low        High       Low        High
C  HC
      DATA TFCHC/
     H-.02062300,-.01438100,-.00240320, .00132190,-.00100810, .00347990,
     H-.02446200,-.01255200,-.00320170, .00426670,-.00086884, .00758430,
     H-.02125500,-.01088800,-.00052755,-.00047925, .00093659, .00766660,
     H-.02142700,-.00661070,-.00039442, .00262880, .00497310, .01232000,
     H-.02351700,-.01409500,-.00880570, .02617900,-.00162220, .02429700,
     H-.02682,    .00000000,-.0075815,  .00000000,-.0051660,  .00000000,
     H-.0327753,  .00000000,-.0083176,  .00000000,-.0090264,  .00000000,
     H-.032082,   .00000000,-.008513,   .00000000,-.0090264,  .00000000,
     H .00000000, .00000000, .00000000, .00000000, .00000000, .00000000/
C  CO  (Note: Bag 1/LowT MYG's 6 - 8 0.0 coefs -> neutral mult adj to TAMBAG)
      DATA TFCCO/
     C-.01348700,-.01469100, .00157840, .00374620, .00110970, .01101400,
     C-.02112600,-.03876700,-.00152890, .00846850, .00157490, .02517900,
     C-.02084300,-.02116500,-.00599510, .02360300, .00182530, .02848300,
     C-.01909100,-.01314600,-.00042373, .02471700, .00579820, .02584800,
     C-.02483500,-.01961200,-.00883360, .04853700,-.00115530, .03143900,
     C .0,        .00000000,-.012478,   .00000000,-.0074106,  .00000000,
     C .0,        .00000000,-.014584,   .00000000,-.011371,   .00000000,
     C .0,        .00000000,-.014584,   .00000000,-.011371,   .00000000,
     C .00000000, .00000000, .00000000, .00000000, .00000000, .00000000/
C  NOX
      DATA TFCNOX/
     N-.00016897, .00388410,-.00892450,-.00873250,-.00725800,-.01083900,
     N-.00025074,-.00103890,-.00597910,-.00924660,-.00626900,-.01010800,
     N .00388550,-.01830100,-.00241560,-.01092500,-.00211880,-.01804200,
     N-.45504E-4,-.00714200,-.00125750,-.00879100,-.00053153,-.00754700,
     N-.00760440,-.02615300,-.00680450,-.01860300,-.00541980,-.02087800,
     N-.0019,     .00000000,-.0061656,  .00000000,-.0049643,  .00000000,
     N-.0045479,  .00000000,-.0074823,  .00000000,-.0090882,  .00000000,
     N-.0047657,  .00000000,-.006989,   .00000000,-.0090882,  .00000000,
     N .00000000, .00000000, .00000000, .00000000, .00000000, .00000000/
C
C  IVLDG fuel injected (Bag 1 / LowT coef: see Note above on TFCCO)
C
      DATA TTFD/
C  CARB
     H-.039143,             -.0060254,            -.0047641,
     C .0,                  -.011402,             -.0090777,
     N-.0043258,            -.010664,             -.010515,
C  TBI
     H-.033227,             -.013174,             -.010889,
     C .0,                  -.022124,             -.0090777,
     N-.0043258,            -.0082652,            -.010515,
C  PFI
     H-.021306,             -.0070187,            -.005731,
     C .0,                  -.013518,             -.0090777,
     N-.0043258,            -.0053098,            -.0061644/
C
C  HDGV
C
      DATA TT4/
C  HC
     H-.00589030, .00134580,
     H-.00738700, .00523170,
     H-.00497590, .00546510,
     H-.00285490, .01008200,
     H-.00741070, .02054600,
     H-.00928590, .00000000,
C  CO
     C-.00205760, .0081720,
     C-.00455410, .02026800,
     C-.00428990, .02412700,
     C-.00130850, .02206100,
     C-.00771170, .02701900,
     C-.00601950, .00000000,
C  NOX
     N-.00643150,-.00839860,
     N-.00554560,-.00868800,
     N-.00139690,-.01807900,
     N-.00046352,-.00748890,
     N-.00575240,-.02159300,
     N-.00197330, .00000000/
C
C  MC
C
      DATA TT8/
C  HC
     H-.02062300,-.01438100,-.00240320, .00132190,-.00100810, .00347990,
     H-.02446200,-.01255200,-.00320170, .00426670,-.00086884, .00758430,
     H-.02125500,-.01088800,-.00052755,-.00047925, .00093659, .00766660,
     H .00000000, .00000000, .00000000, .00000000, .00000000, .00000000,
C  CO
     C-.01348700,-.01469100, .00157840, .00374620, .00110970, .01101400,
     C-.02112600,-.03876700,-.00152890, .00846850, .00157490, .02517900,
     C-.02084300,-.02116500,-.00599510, .02360300, .00182530, .02848300,
     C .00000000, .00000000, .00000000, .00000000, .00000000, .00000000,
C  NOX
     N-.00016897, .00388410,-.00892450,-.00873250,-.00725800,-.01083900,
     N-.00025074,-.00103890,-.00597910,-.00924660,-.00626900,-.01010800,
     N .00000000, .00000000, .00000000, .00000000, .00000000, .00000000,
     N .00000000, .00000000, .00000000, .00000000, .00000000, .00000000/
C
C  /TEMPC2/:
C
      DATA MAXTCF/9/,MYGTFD/1983,1986,1986/,ISHIFT/5/
C
C  Note that the last MYG ub in MYGTCF for IVLDG is the year before the
C  MYGTFD FDS differentiation start year.  ITCPTR is not called from that
C  start year onward, hence there is no need to cover that range in MYGTCF.
C
      DATA MYGTCF/
C  MYG#:   1    2    3    4    5    6    7    8    9
C  LDGV
     H   1967,1969,1971,1974,1979,1980,1981,1982,0000,
     C   1967,1969,1971,1974,1979,1980,1981,1982,0000,
     N   1967,1972,1974,1976,1979,1980,1981,1982,0000,
C  LDGT1
     H   1967,1969,1971,1974,1980,1983,1984,1985,0000,
     C   1967,1969,1971,1974,1980,1983,1984,1985,0000,
     N   1967,1972,1974,1978,1980,1983,1984,1985,0000,
C  LDGT2
     H   1969,1973,1978,0000,1980,1983,1984,1985,0000,
     C   1969,1973,1978,0000,1980,1983,1984,1985,0000,
     N   1969,1973,1978,0000,1980,1983,1984,1985,0000,
C  HDGV
     H   1969,1973,1978,1983,1984,2020,0000,0000,0000,
     C   1969,1973,1978,1983,1984,2020,0000,0000,0000,
     N   1969,1973,1978,1983,1984,2020,0000,0000,0000,
C  MC
     H   1977,1979,2020,0000,0000,0000,0000,0000,0000,
     C   1977,1979,2020,0000,0000,0000,0000,0000,0000,
     N   1977,2020,0000,0000,0000,0000,0000,0000,0000/
C
C  /TEMPC3/:
C
      DATA MDLOHI/3*0/,TDIFF/3*0.0/,CHFRAC/3*0.0/,TCF/3*1.0/
C
C  /TEMPC4/:
C
      DATA MYGCOO/1980,1984,1984/, LOWCO/1/
C
C                    MY1      MY2       MY3
      DATA COOTFC/-1.2448, -1.30945, -1.28402/
C
C                   CARB      TBI       PFI
      DATA COOTFD/-1.2262, -1.2262,  -0.65727/
C
C  /TEMPC5/:
C
C  From the Cold CO CAA: LDGV=10.0 g/mi 1994+, LDGT1/2=12.0 g/mi 1994+
C  (Note: code allows for two standards, if ever needed again.)
C
      DATA NCC/1/,MYCOL/1994,1994/,CCSTD/10.0,2*12.0,10.0,2*12.0/
C
C  /TEMPC6/:
C
      DATA CSUB/6*0.0/,CMUL/6*1.0/
C
C  COLD CO Phase-in Percentages for LDGV/T.
C
C         1994   1995   1996+
      DATA PIVPCT/
     1    40.0,  80.0, 100.0,
     2    40.0,  80.0, 100.0,
     3    40.0,  80.0, 100.0,
     1    40.0,  80.0, 100.0,
     2    40.0,  80.0, 100.0,
     3    40.0,  80.0, 100.0/
C
C  /TEMPC7/:
C
C  These are dummy values.  Actual values are set during execution.
C
      DATA ICOLD/1/,ADDCSW/1/,MULCSW/1/
C
C  /TEMPS/:
C
      DATA AMBT,TEMMIN,TEMMAX,TEMEXH,TEMEVP/12*0.0/
C
      END
