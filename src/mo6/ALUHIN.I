C
C  Changes: (Last change first)
C
C 02 Jul 02 AIR Task 36: Added hourly relative humidity variables/arrays
C           and BARO variables
C 21 Sep 00 AIR Task 07 Bug Fix 256: Added GCLD_CVR for reinitialization
C           via INISCN.
C 13 Jan 00 @DynTel ddj 1-013, Removed variables for towing and
C           trailers and wetbulb and drybulb temperatures.  Added
C           variables for sunrise, sunset and peak sun.
C 12 Dec 96 DynTel@MLA Expanded arrays TRAILR and XLOAD to accomodate the 
C           new LDGT vehicle categories.
C
C The initialization of the data can be found in BD27.
C
      INTEGER GSUNRISE
      INTEGER GSUNSET
      INTEGER GPEAKSUN
      LOGICAL RH_DIST, BARO_FLAG
C
      REAL AC, ABSHUM, GLAC, GLAHUM, GCLD_CVR, RH_HRLY, RHMIN, RHMAX,
     *     GBARO, BARO
C
      COMMON /ALUHIN/ AC, ABSHUM, GLAC, GSUNRISE, GSUNSET, GPEAKSUN(2)
     &     ,GLAHUM,GCLD_CVR, RH_HRLY(24), RHMAX, RHMIN, BARO, GBARO,
     &      RH_DIST, BARO_FLAG
