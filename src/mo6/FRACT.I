C  FRACT.I 
C 
C Changes: 
C 
C 15 Feb 01 AIR Task 11: New INCLUDE for MOBILE6. 
C 14 Feb 01 AIR Task 11: Removed unused variables. 
C 14 Nov 00 AIR Task 06: Added NDIEYR to common DIESEL. Corrected docs. 
C                        Added NIYRFE to common DIEFE. 
C 10 Nov 00 AIR Task 06: Updated code to MOBILE6 standards. 
C 
C  Common Block Dictionary: 
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
c   NDIEYR  i    Number of years in IDIEYR. 
C 
C  /SULFR/ : 
C 
C   DCNVRT  R    The fraction of primary sulfate, i.e. the fraction of 
C                sulfur that is directly converted into SO4 + bound water 
C   EXPECV  R    The correction factor for RFG for exhaust particulate  
C                emissions for catalyst vehicles. 
C   FDNSTY  R    The density of gasoline (1) and diesel (2) fuel in pounds/gallon 
C   WATER   R    weight ratio of seven water molecules to sulfate 
C                (see p. 6-25 of "Regulatory Impact Analysis, 
C                Control of Sulfur and Aromatics Contents of On- 
C                Highway Diesel Fuel") 
C 
C 
C  The data are initialized in BD04P5. 
C 
C 
      INTEGER IDIEYR 
      INTEGER IYRFE 
      INTEGER NDIEYR 
      INTEGER NIYRFE 
      REAL    DCNVRT 
      REAL    DIESFE 
      REAL    EFDPM 
      REAL    EXPECV 
      REAL    FDNSTY 
      REAL    SOFPNT 
      REAL    WATER 
C 
      COMMON /CARBON/ SOFPNT(12) 
      COMMON /DIEFE/  DIESFE(32,10), IYRFE(32), NIYRFE 
      COMMON /DIESEL/ IDIEYR(9,10), EFDPM(9,10), NDIEYR 
      COMMON /SULFR/  DCNVRT, FDNSTY(2), WATER, EXPECV 
C 
