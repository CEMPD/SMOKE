      BLOCK DATA BD04P5 
C 
C   BLOCK DATA Subprogram 04:  /CARBON/,/DIEFE/,/DIESEL/,/SULFR/ 
C 
C   Block Data 04 initializes 
C      1) the fraction of total carbon that is the soluable organic 
C         fraction (SOF) 
C      2) The diesel fueled vehicle fuel economies based on model year 
C      3) The total exhaust PM10 emission factors for diesel vehicles 
C         in grams/BHP-hr 
C      4) The weight percent of sulfur in gasoline and diesel fuel, the 
C         density of gasoline and diesel fuel (lbs/gallon), the fraction 
C         of the sulfur in the fuel that gets converted directly into 
C         sulfate matter (sulfuric acid), the exstimated weight ratio of 
C         water molecules in sulfuric acid to SO4, the fraction of SO2 
C         exhausted that is estimated converts to sulfate material in the 
C         atmosphere (either ammonium sulfate or ammonium bisulfate), the 
C         weight ratio of the combination of ammonium sulfate and ammonium 
C         bisulfate to SO2. 
C      5) The technology fractions for LDGV, LDGT1, LDGT2, HDGV, and MC for 
C         use in the direct sulfate equations.  
C 
C Changes: 
C 
C 12 Feb 01 AIR Task 11: New BLOCK DATA for MOBILE6. 
C 12 Feb 01 AIR Task 11: Removed FRACA, changed year 2020 to NEYR. 
C 14 Nov 00 AIR Task 06: Added NDIEYR to common DIESEL. Corrected docs. 
C                        Added NIYRFE to common DIEFE. 
C 10 Nov 00 AIR Task 06: Updated code to MOBILE6 standards. 
C 
C  Common Block Array Subscripts: 
C 
C         DIESFE(32,10)   --  DIESFE (IYEAR,IV) 
C         EFDPM(9,10)     --  EFDPM (IYEAR,IV) 
C         IDIEYR(9,10)    --  IDIEYR (IYEAR,IV) 
C         IYRFE(32)       --  IYRFE (IYEAR) 
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
C 
C  /DIEFE/: 
C 
C   DIESFE  R    Fuel economy (miles/gal) for diesel fueled vehicles: 
C   IYRFE   I    The years corresponding to the fuel economy array DIESFE 
C                for the diesel fueled vehicles 
C   NIYRFE  I    Number of years in IYRFE 
C 
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
      REAL DLDV 
      REAL DLDT 
      REAL DBHD 
      REAL DLHD 
      REAL DMHD 
      REAL DHHD 
      REAL DBUS 
      REAL DBUS1 
      REAL DBUS2 
      REAL DBUS3 
C 
      DIMENSION DLDV(32),DLDT(32),DBHD(32),DLHD(32),DMHD(32), 
     *                   DHHD(32),DBUS(32), 
     *                   DBUS1(32),DBUS2(32),DBUS3(32) 
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
      EQUIVALENCE (DIESFE(1,1),DLDV(1)), 
     *            (DIESFE(1,2),DLDT(1)), 
     *            (DIESFE(1,3),DBHD(1)), 
     *            (DIESFE(1,4),DLHD(1)), 
     *            (DIESFE(1,5),DMHD(1)), 
     *            (DIESFE(1,6),DHHD(1)), 
     *            (DIESFE(1,7),DBUS(1)), 
     *            (DIESFE(1,8),DBUS1(1)), 
     *            (DIESFE(1,9),DBUS2(1)), 
     *            (DIESFE(1,10),DBUS3(1)) 
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
C    The years corresponding to the fuel economy array DIESFE for 
C    the diesel fueled vehicles: 
C 
      DATA IYRFE/ 
     *  1975, 1990, 1991, 1992, 1993, 1994, 1995, 1996, 1997, 1998, 
     *  1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007, 2008, 
     *  2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018, 
     *  2019, NEYR/ 
C 
      DATA NIYRFE /32/ 
C 
C 
C   Fuel economy (miles/gal) for diesel fueled vehicles: 
C 
C   NOTE 1:  This data comes from Table 6-A-1 of "Regulatory Impact 
C            Analysis, Control of Sulfur and Aromatics Contents of On- 
C            Highway Diesel Fuel", EPA, June 1990 
C 
C   NOTE 2:  The fuel economy of the vehicle class 2BHDDV was assumed 
C            to be the same as LHDDV since this data was not provided 
C            in Table 6-A-1 of the reference used. 
C 
C   LDDV 
      DATA DLDV/ 
     * 21.85,28.26,28.76,29.32,29.87,30.40,30.92,31.41,31.86,32.29, 
     * 32.71,33.11,33.44,33.72,33.95,34.15,34.31,34.44,34.57,34.66, 
     * 34.74,34.80,34.84,34.88,34.91,34.92,34.94,34.95,34.95,34.96, 
     * 34.96,34.96/ 
C   LDDT 
      DATA DLDT/ 
     * 20.00,23.57,23.76,23.96,24.13,24.28,24.44,24.64,24.84,25.08, 
     * 25.31,25.54,25.74,25.90,26.04,26.15,26.25,26.33,26.40,26.46, 
     * 26.52,26.56,26.60,26.63,26.65,26.67,26.69,26.70,26.71,26.72, 
     * 26.72,26.72/ 
C   2BHDDV 
      DATA DBHD/ 
     *  8.23,14.04,14.19,14.35,14.51,14.67,14.85,15.01,15.19,15.38, 
     * 15.56,15.75,15.92,16.08,16.20,16.37,16.47,16.56,16.63,16.69, 
     * 16.75,16.79,16.82,16.84,16.87,16.89,16.90,16.90,16.91,16.92, 
     * 16.92,16.92/ 
C   LHDDV 
      DATA DLHD/ 
     *  8.23,14.04,14.19,14.35,14.51,14.67,14.85,15.01,15.19,15.38, 
     * 15.56,15.75,15.92,16.08,16.20,16.37,16.47,16.56,16.63,16.69, 
     * 16.75,16.79,16.82,16.84,16.87,16.89,16.90,16.90,16.91,16.92, 
     * 16.92,16.92/ 
C   MHDDV 
      DATA DMHD/ 
     *  6.07, 7.28, 7.37, 7.46, 7.54, 7.63, 7.71, 7.78, 7.86, 7.93, 
     *  8.00, 8.07, 8.13, 8.17, 8.20, 8.22, 8.24, 8.25, 8.26, 8.27, 
     *  8.28, 8.29, 8.29, 8.29, 8.30, 8.30, 8.30, 8.30, 8.31, 8.31, 
     *  8.31, 8.31/ 
C   HHDDV 
      DATA DHHD/ 
     *  4.35, 5.93, 6.02, 6.12, 6.21, 6.31, 6.41, 6.51, 6.61, 6.69, 
     *  6.76, 6.82, 6.90, 6.96, 7.01, 7.05, 7.08, 7.10, 7.11, 7.12, 
     *  7.13, 7.14, 7.15, 7.15, 7.15, 7.16, 7.16, 7.16, 7.16, 7.16, 
     *  7.16, 7.16/ 
C   BUSES 
      DATA DBUS/ 
     *  4.78, 6.06, 6.16, 6.26, 6.39, 6.51, 6.62, 6.73, 6.84, 6.97, 
     *  7.08, 7.20, 7.25, 7.29, 7.34, 7.40, 7.44, 7.49, 7.54, 7.59, 
     *  7.63, 7.67, 7.72, 7.76, 7.81, 7.85, 7.89, 7.94, 7.97, 8.02, 
     *  8.06, 8.06/ 
C 
C  Oct-19-1993 @ ARC-bsg Subtask 325  Making the other busses EF by 
C  creating FE for use in Direct SO4, Exhaust SO2 an Indirect SO4 
C  the same FE number is used for all model years. 
C 
C   Transit buses 
      DATA DBUS1/32*4.04/ 
C 
C   CBD buses 
      DATA DBUS2/32*3.33/ 
C 
C   Heavy Urban Buses 
      DATA DBUS3/32*1.65/ 
C 
      END 
