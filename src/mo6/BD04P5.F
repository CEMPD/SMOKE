      BLOCK DATA BD04P5 
C 
C   BLOCK DATA Subprogram 04:  /CARBON/,/DIESEL/,/SULFR/ 
C 
C   Block Data 04 initializes 
C      1) the fraction of total carbon that is the soluable organic 
C         fraction (SOF) 
C      2) The total exhaust PM10 emission factors for diesel vehicles 
C         in grams/BHP-hr 
C      3) The weight percent of sulfur in gasoline and diesel fuel, the 
C         density of gasoline and diesel fuel (lbs/gallon), the fraction 
C         of the sulfur in the fuel that gets converted directly into 
C         sulfate matter (sulfuric acid), the exstimated weight ratio of 
C         water molecules in sulfuric acid to SO4, the fraction of SO2 
C         exhausted that is estimated converts to sulfate material in the 
C         atmosphere (either ammonium sulfate or ammonium bisulfate), the 
C         weight ratio of the combination of ammonium sulfate and ammonium 
C         bisulfate to SO2. 
C      4) The technology fractions for LDGV, LDGT1, LDGT2, HDGV, and MC for 
C         use in the direct sulfate equations.  
C 
C Changes: 
C 
C 06 Nov 01 AIR Task 27: Moved all fuel economies to BD28. 
C 12 Feb 01 AIR Task 11: New BLOCK DATA for MOBILE6. 
C 12 Feb 01 AIR Task 11: Removed FRACA, changed year 2020 to NEYR. 
C 14 Nov 00 AIR Task 06: Added NDIEYR to common DIESEL. Corrected docs. 
C                        Added NIYRFE to common DIEFE. 
C 10 Nov 00 AIR Task 06: Updated code to MOBILE6 standards. 
C 
C  Common Block Array Subscripts: 
C 
C         EFDPM(9,10)     --  EFDPM (IYEAR,IV) 
C         IDIEYR(9,10)    --  IDIEYR (IYEAR,IV) 
C         FDNSTY(2)       --  FDNSTY (IGASDSL) 
C         SOFPNT(12)      --  SOFPNT (IVPART) 
C 
C  Common Block Dictionary: 
C 
C 
C   Name   Type                Description 
C  ------  ----  ---------------------------------------------- 
C  /CARBON/: 
C 
C   SOFPNT  R    The fraction of the carbon particulate that is 
C                considered the soluable organic fraction (SOF) 
CC 
C  /DIESEL/ : 
C 
C   IDIEYR  I    The model years corresponding to the exhaust PM10 emission 
C                factors for diesel vehicles, EFDPM. 
C   EFDPM   R    The exhaust PM10 emissions in grams/BHP-hr for heavy-duty 
C                diesel vehicles, and in grams/mile for light duty diesels. 
C   NDIEYR  I    Number of years in IDIEYR. 
C 
C  /SULFR/ : 
C 
C   GSWGTR  R    The sulfur weight percent of gasoline fuel 
C   DSWGTR  R    The sulfur weight percent of diesel fuel where 1=1990 and 
C                earlier, 2=1991-1992 and 3=1993+. 
C   DCNVRT  R    The fraction of primary sulfate, i.e. the fraction of 
C                sulfur that is directly converted into SO4 + bound water 
C   EXPECV  R    The correction factor for RFG for exhaust particulate  
C                emissions for catalyst vehicles. 
C   ICNVRT  R    The fraction of primary sulfate, i.e. the fraction of 
C                sulfur that is directly converted into SO4 + bound water 
C                (sulfuric acid) 
C   FDNSTY  R    Fuel density, lb/gal, 1=gas, 2=diesel 
C   WATER   R    weight ratio of seven water molecules to sulfate 
C                (see p. 6-25 of "Regulatory Impact Analysis, 
C                Control of Sulfur and Aromatics Contents of On- 
C                Highway Diesel Fuel") 
C   AMNWGT  R    The extimated weight ratio of the combination of ammonium 
C                bisulfate and ammonium sulfate in the atmosphere to SO4 
C                (see p 6-25). 
C                 
C 
      IMPLICIT NONE 
C 
      INCLUDE 'FRACT.I' 
      INCLUDE 'PARAME.I' 
C 
C 
C  Notes: 
C  Oct-19-1993 @ ARC-bsg Subtask 325  Making the other busses EF by 
C  creating FE for use in Direct SO4, Exhaust SO2 an Indirect SO4 
C  the same FE number is used for all model years. 
C  Feb-7-1994 @ CSC-bsg Request 390 New Sulfate equations 
C  Mar-9-1994 @ CSC-bsg Request 401  Add user input variables for 
C  sulfur weight percent for gasoline and diesel fuel or default values 
C  Mar-15-1994 @ CSC-bsg Request 407 Change the fuel density for both 
C  gasoline and diesel fuels. 
C  Mar-22-1994 @ CSC-bsg  Request 409 Changed SULFLG to RFGFLG          
C  removed user input of sulfur weight percent 
C  The correction factor for exhaust particulate emission factors for 
C  catalyst vehicles 
C  Jun-7-1994 @ CSC-bsg Request 443 changed bus category exhaust particulate 
C  default no trap rates for model years 1992 and later due to standard 
C  changes.  Changed other heavy duty rates for 1994 and later due to a better 
C  prediction of new technology 
C  Jun-22-1994 @ CSC-bsg Request 429 Changed block data LDDV/T emission values 
C  increased array size of IDIEYR and EFDPM to include new LDDV/T standards. 
C  Jul-21-1994 @ CSC-bsg Changed block data LDDV/T emission values and IDIEYR 
C  values to take pre-control standards into account 
C 
C 
      DATA EXPECV/1.0/ 
C 
C     The density of gasoline (1) and  diesel fuel (2) in pounds/gallon: 
C 
      DATA FDNSTY/6.09, 7.11/ 
C 
C     The fraction of primary sulfate from diesel vehicles, i.e. the 
C     fraction of sulfur that is directly converted into SO4 + bound water: 
C 
      DATA DCNVRT/0.02/ 
C 
C     The weight ratio of seven water molecules to sulfate (see p. 6-25 
C     of "Regulatory Impact Analysis, Control of Sulfur and Aromatics 
C     Contents of On-Highway Diesel Fuel") -- used to calculate sulfuric 
C     acid weight in grams. 
C 
      DATA WATER/1.2857/ 
C 
C     The fraction of the non-sulfate particulate that is considered 
C     the soluable organic fraction (SOF): 
C 
C     NOTE:  The data for the vehicle classes LHDDE, MHDDE, and HHDDE was 
C            taken from "Regulatory Impact Analysis, Control of Sulfur 
C            and Aromatics Contents of On-Highway Diesel Fuel", EPA, 
C            June 1990, page 3-9.  The fraction for buses was assumed to 
C            be the same as MHDDV.  For LDDV and LDDT, the data came 
C            from Table 3-9.  The emission factors for LDDT was calculated 
C            as the average of LDDT1 and LDDT2, as well as the SOF percent 
C            of (total PM10 - SO4). 
C 
      DATA SOFPNT/0.,0.,0.,0.,0.,0.22,0.59,0.51,0.51,0.44,0.24,0.44/ 
C 
C  Jun-7-1994 @ CSC-bsg Request 443 changed bus category exhaust particulate 
C  default no trap rates for model years 1992 and later due to standard 
C  changes.  Changed other heavy duty rates for 1994 and later due to a better 
C  prediction of new technology 
C  Jun-22-1994 @ CSC-bsg Request 429 Changed block data LDDV/T emission values 
C  increased array size of IDIEYR and EFDPM to include new LDDV/T standards. 
C  Jul-21-1994 @ CSC-bsg Changed block data LDDV/T emission values and IDIEYR 
C  values to take pre-control standards into account 
C 
      DATA IDIEYR/ 
C    LDDVs 
     * 1900, 1981, 1982, 1985, 1987, 1988, 1991, 1994, 1996, 
C    LDDT 
     * 1900, 1981, 1982, 1985, 1987, 1988, 1991, 1994, 1997, 
C    2BHDDV 
     * 1900, 1987, 1988, 1991, 1994, NEYR, NEYR, NEYR, NEYR, 
C    LHDDV 
     * 1900, 1987, 1988, 1991, 1994, NEYR, NEYR, NEYR, NEYR, 
C    MHDDV 
     * 1900, 1987, 1988, 1991, 1994, NEYR, NEYR, NEYR, NEYR, 
C    HHDDV 
     * 1900, 1987, 1988, 1991, 1994, NEYR, NEYR, NEYR, NEYR, 
C    HBUS 
     * 1900, 1987, 1988, 1991, 1992, 1993, 1994, NEYR, NEYR, 
C    HBUST1 
     * 1900, 1987, 1988, 1991, 1992, 1993, 1994, NEYR, NEYR, 
C    HBUST2 
     * 1900, 1987, 1988, 1991, 1992, 1993, 1994, NEYR, NEYR, 
C    HBUST3 
     * 1900, 1987, 1988, 1991, 1992, 1993, 1994, NEYR, NEYR/ 
C 
      DATA NDIEYR /9/ 
C 
C  Jun-7-1994 @ CSC-bsg Request 443 changed bus category exhaust particulate 
C  default no trap rates for model years 1992 and later due to standard 
C  changes.  Changed other heavy duty rates for 1994 and later due to a better 
C  prediction of new technology 
C  Jun-22-1994 @ CSC-bsg Request 429 Changed block data LDDV/T emission values 
C  increased array size of IDIEYR and EFDPM to include new LDDV/T standards. 
C  Jul-21-1994 @ CSC-bsg Changed block data LDDV/T emission values and IDIEYR 
C  values to take pre-control standards into account 
C 
      DATA EFDPM/ 
C    LDDV 
     * .7,.259,.256,.255,.134,.132,.131,.128,.1, 
C    LDDT 
     * .7,.309,.354,.358,.334,.291,.294,.13,.109, 
C    2BHDDV 
     * .5156,.5156,.5140,.2873,.1011,0.,0.,0.,0., 
C    LHDDV 
     * .5156,.5156,.5140,.2873,.1011,0.,0.,0.,0., 
C    MHDDV 
     * .6946,.6946,.4790,.2747,.0948,0.,0.,0.,0., 
C    HHDDV 
     * .6444,.6444,.4360,.2709,.0836,0.,0.,0.,0., 
C    HBUS 
     * .6931,.6931,.4790,.2772,.1716,.1457,.0591,0.,0., 
C    HBUST1 
     * .6931,.6931,.4790,.2772,.1716,.1457,.0591,0.,0., 
C    HBUST2 
     * .6931,.6931,.4790,.2772,.1716,.1457,.0591,0.,0., 
C    HBUST3 
     * .6931,.6931,.4790,.2772,.1716,.1457,.0591,0.,0./ 
C 
      END 
