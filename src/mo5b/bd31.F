C Version "@(#)$Id$ $Source$ $Date$ 

      BLOCK DATA BD31
C
C  BLOCK DATA Subprogram 31: /RVPEX1/,/RVPEX2/,/RVPNAT/
C
C  BLOCK DATA 31 initializes exhaust emission factor RVP correction factor (cf)
C  calculation variables.
C
C  Common block array subscripts:
C
C  MYGRVP(4,4)      -  MYGRVP ( IG, IVTAM )
C  RVPCF(25,3,4)    -  RVPCF  ( IDX, IP, IVTAM )
C  RVPDIU(2)        -  RVPDIU ( LIMIT )
C  RVPLIM(2)        -  RVPLIM ( LIMIT )
C  RG23(3,3,3,2)    -  RG23   ( IRQ3, IB, IP, IG )
C  R9G1(3,2)        -  R9G1   ( IRQ1, IP )    (IP=1,2)
C  TLL(2)           -  TLL    ( IRQ2 )
C
C  Common block dictionary:
C
C   Name   Type              Dictionary
C  ------  ----  ------------------------------------------------
C  /RVPEX1/:
C
C  R9G1     R    RVP cf normalized at RVP=9.0 equation coeff used by RVP IG=1
C  RG23     R    RVP & high T combined cf equation coeff used by RVP IG=2,3
C                The 1st (a1) coeff is also used to set the ub (FTP T) of the
C                low T IG=2,3 interpolation interval.
C  TLL      R    RVP cf temperature range lb & length for MYG 1 & low T MYG 2&3
C
C  /RVPEX2/:
C
C  MYGRVP   I    upper bound years for RVP CF MYGs
C  RVPCF    R    multiplicative RVP cf's, applied to corresponding exhaust ef's
C                in EFCALX & HCCALX
C
C  /RVPNAT/:
C
C  RVPCER   R    certification fuel control program RVP, fixed wrt fuel type
C  RVPLIM   R    RVP value limits for constraining RVPBAS and RVPIUS
C                1 = lower bound    2 = upper bound
C  RVPDIU   R    RVP value limits for constraining RVPICY, again 1=lb & 2=ub
C  PFUL     R    Percent of fuel tank filled
C
C  Notes:
C
C  Variables of the oxygenated and certification fuels support from MOBILE4
C  were removed in MOBILE4.1.
C
      INCLUDE 'PARAME.I'
      INCLUDE 'RVPEX1.I'
      INCLUDE 'RVPEX2.I'
      INCLUDE 'RVPNAT.I'
C
C  /RVPEX1/: These are the hardcoded terms of the RVP cf equations.
C
C  The continuation line code for R9G1 is IP.
C
      DATA R9G1/
C           Numerator        Denominator
C      Intercept     Slope
     1  0.56222,   0.012512,  0.67483,
     2  7.1656,    0.33413,  10.17277/
C
C  The continuation line code for RG23 is IB.
C
      DATA RG23/
C                 a1           a2           a3
C  IG=2 (= 1980-82 for LDGV)
C  HC
     1         0.091402,   0.0044270,   0.0029466,
     2         0.042060,   0.0048358,   0.0,
     3         0.093179,   0.0074688,   0.0047276,
C  CO
     1         0.091345,   0.0062182,   0.0,
     2         0.139680,   0.0149430,   0.0,
     3         0.163220,   0.0149230,   0.0,
C  NOx
     1         0.0,        0.0,         0.0,
     2        -0.040024,   0.0,         0.0,
     3         0.0,        0.0,         0.0,
C
C  IG=3 (= 1983+ for LDGV)
C  HC
     1         0.023202,   0.0,         0.0,
     2         0.153730,   0.0086550,   0.0,
     3         0.132630,   0.0083730,   0.0056009,
C  CO
     1         0.040748,   0.0035170,   0.0,
     2         0.262140,   0.0149660,   0.0056416,
     3         0.232180,   0.0206950,   0.0082344,
C  NOx
     1         0.014219,   0.0,         0.0,
     2         0.027491,   0.0037789,   0.0,
     3         0.0,        0.0,         0.0/
C
      DATA TLL/
C  LB of Temperature Interval      Interval Length
     *        45.0,                      30.0/
C
C  /RVPEX2/: RVPCF is computed by RVPEXH.
C
C  The MYGRVP continuation code is IVTAM.
C
      DATA MYGRVP/
C                   MYG0   MYG1   MYG2   MYG3
     1              1970,  1979,  1982,  NEYR,
     2              1970,  1980,  1985,  NEYR,
     3              1978,  1980,  1985,  NEYR,
     4              1984,  0000,  NEYR,  0000/
C
      DATA  RVPCF/300*1.0/
C
C  /RVPNAT/: All variables but RVPLIM were required to be read into One-time
C            Data.  MOBILE4 hardcoded all the parameters, removing them
C            from user control.  MOBILE4 added RVPDIU derived from RVPLIM,
C            so initialized here to the same limits.
C
      DATA RVPLIM/6.5,15.2/,RVPDIU/5.0,15.2/
      DATA PFUL/54.57/
C
      END
