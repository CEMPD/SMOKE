C Version "@(#)$Id$ $Source$ $Date$ 

      BLOCK DATA BD22
C
C  BLOCK DATA Subprogram 22: /ALUPAC/,/ALUPLT/
C
C  BLOCK DATA 22 defines the last years of the model year groups used for
C  the a/c, XL, & TW correction factors.
C
C  Common block array subscripts:
C
C  MYGACE(8,3,2)  -  MYGACE ( IG, IVP, IR )
C  MYGACC(2,3,2)  -  MYGACC ( IG, IVP, IR )
C  MYGXLC(6,3,2)  -  MYGXLC ( IG, IVP, IR )
C  MYGTWC(6,3,2)  -  MYGTWC ( IG, IVP, IR )
C
C  Common block dictionary:
C
C   Name   Type             Description
C  ------  ----  ----------------------------------------
C  /ALUPAC/:
C
C  MYGACE   I    the last year of the model year groups for ACEQIP
C                (fraction of vehicles equipped with a/c)
C  MYGACC   I    the last year of the model year groups for ACCF
C                (a/c correction factors)
C
C  /ALUPLT/:
C
C  MYGXLC   I    the last year of the model year groups for XLCF
C                (extra load correction factors)
C  MYGTWC   I    the last year of the model year groups for TWCF
C                (trailer towing correction factors)
C
C  Notes:
C
C  None.
C
C
      INCLUDE 'ALUPAC.I'
      INCLUDE 'ALUPLT.I'
      INCLUDE 'PARAME.I'
C
      DATA MYGACE/
C  Low altitude region
     *   1961,1964,1966,1968,1971,1976,NEYR,0000,
     *   1976,1977,NEYR,0000,0000,0000,0000,0000,
     *   1976,1977,NEYR,0000,0000,0000,0000,0000,
C  High altitude region
     *   1961,1964,1966,1968,1971,1976,NEYR,0000,
     *   1976,1977,NEYR,0000,0000,0000,0000,0000,
     *   1976,1977,NEYR,0000,0000,0000,0000,0000/
C
      DATA MYGACC/
C  Low altitude region
     *   1974,NEYR,
     *   1974,NEYR,
     *   1978,NEYR,
C  High altitude region
     *   1974,NEYR,
     *   1974,NEYR,
     *   1978,NEYR/
C
      DATA MYGXLC/
C  Low altitude region
     *   1967,1969,1971,1972,1974,NEYR,
     *   1967,1969,1971,1972,1974,NEYR,
     *   1969,1973,1978,NEYR,0000,0000,
C  High altitude region
     *   1967,1969,1971,1972,1974,NEYR,
     *   1967,1969,1971,1972,1974,NEYR,
     *   1969,1973,1978,NEYR,0000,0000/
C
      DATA MYGTWC/
C  Low altitude region
     *   1967,1969,1971,1972,1974,NEYR,
     *   1967,1969,1971,1972,1974,NEYR,
     *   1969,1973,1978,NEYR,0000,0000,
C  High altitude region
     *   1967,1969,1971,1972,1974,NEYR,
     *   1967,1969,1971,1972,1974,NEYR,
     *   1969,1973,1978,NEYR,0000,0000/
C
      END
