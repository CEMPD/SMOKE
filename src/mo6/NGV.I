C 
C  NGV.I contains the variables that hold the VALUES OF the NGV user 
C  supplied Emission Factors. 
C 
C  Changes: (Last change first) 
C 
C   1 Jun 00 @EPA-ddj E11, Added NGVFRAC. 
C  24 Apr 00 @EPA-ELG increased size of NGVPHIN array 
C  24 Feb 00 @EPA-elg, New Information File 
C 
C 
C  Common Block Dictionary: 
C 
C   Name    Type                  Description 
C  ------   ----  ----------------------------------------------------- 
C  NATGAS    R    Final NGV Emission Factors passed to HRLOOP. 
C  NGVDET    R    Deterioration rate (10,000K miles) of NGV vehicles. 
C  NGVFLG1   I    Toggle for Alternate NGV Phase-in Fractions 0-DEFAULT 1-Alternate. 
C  NGVFLG2   I    Toggle for Alternate NGV Emission Factors 0-DEFAULT 1-Alternate. 
C  NGVHIGH   R    High Emitter Emission Level for NGVs. 
C  NGVPHIN   R    NGV Phase-in Fractions BY model year and vehicles TYPE. 
C  NGV_START I    Earliest NGV Vehicle that can be modeled.  Default = 1994 MY 
C  NGVZML    R    Zero mile Emission Level of NGV vehicles. 
C 
C 
C 
C 
C  Notes: 
C 
C 
C 
      INTEGER    NGVFLG1 
      INTEGER    NGVFLG2 
      INTEGER    NGV_START 
C 
      REAL,  DIMENSION(28,3,2)      :: NGVZML 
      REAL                          :: NGVFRAC 
      REAL,  DIMENSION(28,3,2)      :: NGVDET 
      REAL,  DIMENSION(28,3,2)      :: NGVHIGH 
      REAL,  DIMENSION(57,28)       :: NGVPHIN 
      REAL,  DIMENSION(3,5)         :: NGV_NZML 
      REAL,  DIMENSION(3,5)         :: NGV_NDET 
      REAL,  DIMENSION(3,5)         :: NGV_HIGH 
      REAL,  DIMENSION(3,3)         :: NGV_HDZML 
      REAL,  DIMENSION(3,3)         :: NGV_HDDET 
      REAL,  DIMENSION(3,2)         :: NGV_P 
      REAL,  DIMENSION(28,3,25,2,2) :: NATGAS 
C 
C 
C     NGVDAT2 Contains the User Supplied NGV Data if any. 
C 
      COMMON /NGVDAT1/ NGVZML, NGVDET, NGVHIGH, NGVPHIN, 
     *                 NGVFLG1, NGVFLG2, NGV_START, 
     *                 NATGAS 
C 
C     NGVDAT2 Contains the Default NGV Data. 
C 
      COMMON /NGVDAT2/ NGV_NZML, NGV_NDET, NGV_HIGH, 
     *                 NGV_HDZML, NGV_HDDET, NGV_P 
 
C 
C     NGVDAT3 Contains MOBILE6 generated NGV data. 
C 
      COMMON /NGVDAT3/ NGVFRAC 
