C Version "@(#)$Id$ $Source$ $Date$ 

      BLOCK DATA BD21
C
C  BLOCK DATA Subprogram 21: /ALUHAC/,/ALUHLT/
C
C  BLOCK DATA 21 initializes the correction factors for air conditioning
C  usage (a/c or AC), extra load (XL), and trailer towing (TW).
C
C  Common block array subscripts:
C
C  ACEQIP(8,3,2)    -  ACEQIP ( IGEPTR, IVP, IR )
C  ACCF(2,3,2,3,2)  -  ACCF ( IAC, IP, IGUPTR, IVP, IR )
C  XLCF(3,6,3,2)    -  XLCF ( IP, IGLPTR, IVP, IR )
C  TWCF(3,6,3,2)    -  TWCF ( IP, IGTPTR, IVP, IR )
C  IAC = AC correction equation coefficient number (2 terms per equation)
C
C  Common block dictionary:
C
C   Name   Type              Description
C  ------  ----  ----------------------------------------------
C  /ALUHAC/:
C
C  ACEQIP   R    fractions of LDGV & LDGT1/2 vehicles equipped with a/c
C  ACCF     R    a/c correction factors
C
C  /ALUHLT/:
C
C  XLCF     R    extra load correction factors
C  TWCF     R    trailer towing correction factors
C
C  Local array subscripts:
C
C  TWCFLO(3,6,3)  -  TWCFLO ( IP, IGLPTR, IVP )
C  TWCFHI(3,6,3)  -  TWCFHI ( IP, IGLPTR, IVP )
C  XLCFLO(3,6,3)  -  XLCFLO ( IP, IGLPTR, IVP )
C  XLCFHI(3,6,3)  -  XLCFHI ( IP, IGLPTR, IVP )
C
C  Local array dictionary:
C
C   Name   Type              Description
C  ------  ----  -----------------------------------------------
C  TWCFLO   R    trailer towing correction factors for low altitude region
C  TWCFHI   R    trailer towing corrrection factors for high altitude region
C  XLCFLO   R    extra load correction factors for low altitude region
C  XLCFHI   R    extra load correction factors for high altitude region
C
C  Notes:
C
C  None.
C
C
      INCLUDE 'ALUHAC.I'
      INCLUDE 'ALUHLT.I'
C
      DIMENSION TWCFLO(3,6,3),TWCFHI(3,6,3)
      DIMENSION XLCFLO(3,6,3),XLCFHI(3,6,3)
C
      EQUIVALENCE
     *   (TWCF(1,1,1,1),TWCFLO(1,1,1)),
     *   (TWCF(1,1,1,2),TWCFHI(1,1,1)),
     *   (XLCF(1,1,1,1),XLCFLO(1,1,1)),
     *   (XLCF(1,1,1,2),XLCFHI(1,1,1))
C
C  Fraction of vehicles equipped with a/c
      DATA ACEQIP/
C  Low altitude region
     *   .07,  .14,  .24,  .37,  .51,  .61,  .72,  0.0,
     *   .32,  .52,  .39,  0.0,  0.0,  0.0,  0.0,  0.0,
     *   .32,  .52,  .39,  0.0,  0.0,  0.0,  0.0,  0.0,
C  High altitude region
     *   .07,  .14,  .24,  .37,  .51,  .61,  .72,  0.0,
     *   .32,  .52,  .39,  0.0,  0.0,  0.0,  0.0,  0.0,
     *   .32,  .52,  .39,  0.0,  0.0,  0.0,  0.0,  0.0/
C
C  A/c correction factors
      DATA ACCF/
C  Low altitude region
     * 1.0230E+0, 0.3344E-2, 1.2020E+0, 0.1808E-2, 1.2990E+0, 0.5643E-4,
     * 1.0000E+0, 0.3512E-2, 1.1300E+0, 0.1528E-2, 1.2210E+0, 0.4262E-3,
     * 1.0230E+0, 0.3344E-2, 1.2020E+0, 0.1808E-2, 1.2990E+0, 0.5643E-4,
     * 1.0000E+0, 0.3512E-2, 1.1300E+0, 0.1528E-2, 1.2210E+0, 0.4262E-3,
     * 1.0230E+0, 0.3344E-2, 1.2020E+0, 0.1808E-2, 1.2990E+0, 0.5643E-4,
     * 1.0000E+0, 0.3512E-2, 1.1300E+0, 0.1528E-2, 1.2210E+0, 0.4262E-3,
C  High altitude region
     * 1.0230E+0, 0.3344E-2, 1.2020E+0, 0.1808E-2, 1.2990E+0, 0.5643E-4,
     * 1.0000E+0, 0.3512E-2, 1.1300E+0, 0.1528E-2, 1.2210E+0, 0.4262E-3,
     * 1.0230E+0, 0.3344E-2, 1.2020E+0, 0.1808E-2, 1.2990E+0, 0.5643E-4,
     * 1.0000E+0, 0.3512E-2, 1.1300E+0, 0.1528E-2, 1.2210E+0, 0.4262E-3,
     * 1.0230E+0, 0.3344E-2, 1.2020E+0, 0.1808E-2, 1.2990E+0, 0.5643E-4,
     * 1.0000E+0, 0.3512E-2, 1.1300E+0, 0.1528E-2, 1.2210E+0, 0.4262E-3/
C
C  Extra load correction factors
C        HC      CO      NOX
C  Low altitude region
      DATA XLCFLO/
C  LDGV
     *   1.0786, 1.2765, 0.9535,
     *   1.0495, 1.1384, 1.0313,
     *   1.0852, 1.2478, 1.0313,
     *   1.0556, 1.1347, 1.0313,
     *   1.0556, 1.1347, 1.0753,
     *   1.0455, 1.3058, 1.0719,
C  LDGT1
     *   1.0786, 1.2765, 0.9535,
     *   1.0495, 1.1384, 1.0313,
     *   1.0852, 1.2478, 1.0313,
     *   1.0556, 1.1347, 1.0313,
     *   1.0556, 1.1347, 1.0753,
     *   1.0455, 1.3058, 1.0719,
C  LDGT2
     *   1.0786, 1.2765, 0.9535,
     *   1.0495, 1.1384, 1.0313,
     *   1.0556, 1.1347, 1.0753,
     *   1.0455, 1.3058, 1.0719,
     *   0.0000, 0.0000, 0.0000,
     *   0.0000, 0.0000, 0.0000/
C
C  High altitude region
      DATA XLCFHI/
C  LDGV
     *   1.0786, 1.2765, 0.9535,
     *   1.0495, 1.1384, 1.0313,
     *   1.0852, 1.2478, 1.0313,
     *   1.0556, 1.1347, 1.0313,
     *   1.0556, 1.1347, 1.0753,
     *   1.0455, 1.3058, 1.0719,
C  LDGT1
     *   1.0786, 1.2765, 0.9535,
     *   1.0495, 1.1384, 1.0313,
     *   1.0852, 1.2478, 1.0313,
     *   1.0556, 1.1347, 1.0313,
     *   1.0556, 1.1347, 1.0753,
     *   1.0455, 1.3058, 1.0719,
C  LDGT2
     *   1.0786, 1.2765, 0.9535,
     *   1.0495, 1.1384, 1.0313,
     *   1.0556, 1.1347, 1.0753,
     *   1.0455, 1.3058, 1.0719,
     *   0.0000, 0.0000, 0.0000,
     *   0.0000, 0.0000, 0.0000/
C
C  Trailer towing correction factors
C        HC      CO      NOX
C  Low altitude region
      DATA TWCFLO/
C  LDGV
     *   1.2614, 1.9327, 1.1184,
     *   1.2762, 1.8940, 1.1384,
     *   1.4598, 2.4753, 1.1384,
     *   1.7288, 2.1414, 1.1384,
     *   1.7288, 2.1414, 1.2170,
     *   1.5909, 3.9722, 1.3875,
C  LDGT1
     *   1.2614, 1.9327, 1.1184,
     *   1.2762, 1.8940, 1.1384,
     *   1.4598, 2.4753, 1.1384,
     *   1.7288, 2.1414, 1.1384,
     *   1.7288, 2.1414, 1.2170,
     *   1.5909, 3.9722, 1.3875,
C  LDGT2
     *   1.2614, 1.9327, 1.1184,
     *   1.2762, 1.8940, 1.1384,
     *   1.7288, 2.1414, 1.2170,
     *   1.5909, 3.9722, 1.3875,
     *   0.0000, 0.0000, 0.0000,
     *   0.0000, 0.0000, 0.0000/
C
C  High altitude region
      DATA TWCFHI/
C  LDGV
     *   1.2614, 1.9327, 1.1184,
     *   1.2762, 1.8940, 1.1384,
     *   1.4598, 2.4753, 1.1384,
     *   1.7288, 2.1414, 1.1384,
     *   1.7288, 2.1414, 1.2170,
     *   1.5909, 3.9722, 1.3875,
C  LDGT1
     *   1.2614, 1.9327, 1.1184,
     *   1.2762, 1.8940, 1.1384,
     *   1.4598, 2.4753, 1.1384,
     *   1.7288, 2.1414, 1.1384,
     *   1.7288, 2.1414, 1.2170,
     *   1.5909, 3.9722, 1.3875,
C  LDGT2
     *   1.2614, 1.9327, 1.1184,
     *   1.2762, 1.8940, 1.1384,
     *   1.7288, 2.1414, 1.2170,
     *   1.5909, 3.9722, 1.3875,
     *   0.0000, 0.0000, 0.0000,
     *   0.0000, 0.0000, 0.0000/
C
      END
