      MODULE DDDATA 
C 
C 
C  Changes: (Last change first) 
C 
C 13 Feb 01 @EPA-elg    Added bus category to each array. 
C 13 Feb 01 @EPA-elg    Changed Pull Ahead Emission Delta to reflect 
C                       Changes in the underlying emission factors. 
C 09 Feb 01 @EPA-djb    Updated IRBLD, NEFFCT, RBLD & FLEET values. 
C                       Added IRBSTRT and IRBEND variables. 
C 01 May 00 @EPA- BG    Removed use of DATA statement initializations 
C 28 Apr 00 @EPA-elg    Corrected 8A Vehicle FLEET block data. 
C 06 Dec 99 @DynTel-bg  1-004  defined logical variable FIRST_DEFEAT 
C 19 Oct 99 @DynTel-NH  1-004  New for the HDDV defeat task 
C 
C  Dimension 18 reflects the MY's and 4 reflects the vehicle classes. 
C  In RBLD - 7 reflects the MY's and 4 reflects the vehicle classes. 
C 
C  Values for NEFFCT are the increase in NOx emission levels 
C  as the result of the heavy-duty NOx defeat device. 
C 
C  Light duty includes all diesel except classes 6,7,8a and 8b. 
C 
C..First entry is for model year 2004.  Last entry is 
C  model year 1987 and earlier.  Negatives reflect Pull ahead effects. 
C  Value of Zero is assigned if no Pull Ahead or effect.  There are no 
C  effects for model years 2004+ or 1987 and earlier. 
C 
C..URBAN is the fraction of the time that the defeat device operates 
C  when an equipped vehicle is operating in Urban type driving. 
C 
C..ARTERL is the fraction of the time that the defeat device operates 
C  when an equipped vehicle is operating in Arterial type driving. 
C 
C..EXPRES is the fraction of the time that the defeat device operates 
C  when an equipped vehicle is operating in Expressway type driving. 
C 
C..FLEET is the fraction of the fleet that is equipped with a defeat 
C  device. 
C 
C..First year is 1999+ and last is 1993 and earlier. 
C  Rebuilds only in 1994 through 1998 MYs for Plan A. 
C..7.00 g/bhp-hr is the rebuild emission target for HDD8a and HDD8b. 
C..6.00 g/bhp-hr is the rebuild emission target for MDDV. 
C  Rebuild effect is the difference between the value and 7.00 or the 
C  value and 6.00. 
C 
C..IRBLD is the age of a vehicle when the first rebuild is done. 
C  parameters are available for the four vehicle classes. 
C  These parameters are currently implicit in the results;  thus, 
C  the array IRBLD is not currently used in Defeat. 
C 
C..IRBSTRT is the first model year where defeat device effects are 
C  rebuilt. Rebuilds are affected by IRBLD parameters. 
C 
C..IRBEND is the last model year where defeat device effects are rebuilt. 
C 
C 
      IMPLICIT NONE 
      SAVE 
C 
      INTEGER, PARAMETER :: NCY = 18 
      INTEGER, PARAMETER :: NVCLS = 5 
C 
      REAL                       :: RBLEV=0.90 
      INTEGER                    :: IBUILD=0 
      INTEGER                    :: IPULL=0 
      LOGICAL                    :: DEFEAT_FLAG 
      LOGICAL                    :: FIRST_DEFEAT 
C 
      INTEGER, DIMENSION(NVCLS)  :: IRBLD=(/ 99, 12, 5, 5, 99 /) 
      INTEGER                    :: IRBSTRT=1994 
      INTEGER                    :: IRBEND=1998 
C 
      REAL, DIMENSION(18,5) :: NEFFCT=RESHAPE((/ 
C  LIGHT DIESELS 
     *0.0000, -1.1400, -0.2850, 0.0000, 0.000, 
     *0.0000,  0.0000,  0.0000, 0.0000, 
     *0.0000,  0.0000,  0.0000, 0.0000, 
     *0.0000,  0.0000,  0.0000, 0.0000, 0.000, 
C  MEDIUM DIESELS (CLS 6 & 7) 
     * 0.0000, -1.5700,-0.3925, 0.0000, 0.0000, 
     * 1.3533,  3.3623, 2.7163, 2.7408, 
     * 2.9136,  2.4900, 2.5700, 0.0680, 
     * 0.0680,  0.7890, 1.5000, 1.5000, 0.0000, 
C  8A HEAVY-DUTIES 
     * 0.0000, -1.5600,-0.3900, 2.3200, 2.3200, 
     * 2.3200,  8.2377, 7.1143, 6.9407, 
     * 7.1260,  6.9118, 7.4985, 7.9025, 
     * 8.8713,  8.3125, 8.9200, 8.9200, 0.0000, 
C  8B HEAVY-DUTIES 
     * 0.0000, -1.5600,-0.3900, 2.3200, 2.3200, 
     * 2.3200,  8.2617, 7.2067, 6.9501, 
     * 7.1022,  6.8980, 7.3978, 7.7489, 
     * 8.8179,  7.8961, 8.9200, 8.9200, 0.0000, 
C  DIESEL BUSES 
     *0.0000, -1.1400, -0.2850, 0.0000, 0.000, 
     *0.0000,  0.0000,  0.0000, 0.0000, 
     *0.0000,  0.0000,  0.0000, 0.0000, 
     *0.0000,  0.0000,  0.0000, 0.0000, 0.000/),(/18,5/)) 
C 
C 
      REAL, DIMENSION(18,5) :: RBLD=RESHAPE((/ 
C  LIGHT DIESEL 
     *0.0000, -1.6300, -0.4075, 0.0000, 0.000, 
     *0.0000,  0.0000,  0.0000, 0.0000, 
     *0.0000,  0.0000,  0.0000, 0.0000, 
     *0.0000,  0.0000,  0.0000, 0.0000, 0.000, 
C  MEDIUM DIESELS (CLS 6 & 7) 
     * 0.0000, -1.8500,-0.4625, 0.0000, 0.0000, 
     * 1.3533,  2.3152, 1.4226, 1.4251, 
     * 1.4424,  1.4000, 2.5700, 0.0680, 
     * 0.0680,  0.7890, 1.5000, 1.5000, 0.0000, 
C  8A HEAVY-DUTIES 
     * 0.0000,-1.8400, -0.4600, 2.3200, 2.3200, 
     * 2.3200,  3.8118, 2.8624, 2.8451, 
     * 2.8636,  2.8422, 7.4985, 7.9025, 
     * 8.8713,  8.3125, 8.9200, 8.9200, 0.0000, 
C  8B HEAVY-DUTIES 
     * 0.0000, -1.8400,-0.4600, 2.3200, 2.3200, 
     * 2.3200,  3.8142, 2.8717, 2.8460, 
     * 2.8612,  2.8408, 7.3978, 7.7489, 
     * 8.8179,  7.8961, 8.9200, 8.9200, 0.0000, 
C  BUSES 
     *0.0000, -1.6300, -0.4075, 0.0000, 0.000, 
     *0.0000,  0.0000,  0.0000, 0.0000, 
     *0.0000,  0.0000,  0.0000, 0.0000, 
     *0.0000,  0.0000,  0.0000, 0.0000, 0.000/),(/18,5/)) 
C 
C 
      REAL, DIMENSION(18,5) :: URBAN=RESHAPE((/ 
C  LIGHT DIESEL 
     *0.0000,  1.0000,  1.0000, 0.0000, 0.000, 
     *0.0000,  0.0000,  0.0000, 0.0000, 
     *0.0000,  0.0000,  0.0000, 0.0000, 
     *0.0000,  0.0000,  0.0000, 0.0000, 0.000, 
C  MEDIUM DIESELS (CLS 6 & 7) 
     * 0.0000,  1.0000,  1.0000, 0.0009, 0.0009, 
     * 0.0009,  0.0009, 0.0042, 0.0049, 
     * 0.0692,  0.1200, 0.1200, 0.0000, 
     * 0.0000,  0.0000, 0.0000, 0.0000,  0.0000, 
C  8A HEAVY-DUTIES 
     * 0.0000,  1.0000,  1.0000, 0.0210, 0.0210, 
     * 0.0210,  0.0210, 0.0083, 0.0107, 
     * 0.0468,  0.0341, 0.0271, 0.0188, 
     * 0.0295,  0.0612, 0.0000, 0.0000,  0.0000, 
C  8B HEAVY-DUTIES 
     * 0.0000,  1.0000,  1.0000, 0.0267, 0.0267, 
     * 0.0267,  0.0267, 0.0116, 0.0115, 
     * 0.0487,  0.0362, 0.0299, 0.0216, 
     * 0.0363,  0.0753, 0.0000, 0.0000,  0.0000, 
C  DIESEL BUSES 
     *0.0000,  1.0000,  1.0000, 0.0000, 0.000, 
     *0.0000,  0.0000,  0.0000, 0.0000, 
     *0.0000,  0.0000,  0.0000, 0.0000, 
     *0.0000,  0.0000,  0.0000, 0.0000, 0.000/),(/18,5/)) 
C 
C 
      REAL, DIMENSION(18,5) :: ARTERL=RESHAPE((/ 
C  LIGHT DIESEL 
     * 0.0000,  1.0000,  1.0000, 0.0000, 0.000, 
     * 0.0000,  0.0000,  0.0000, 0.0000, 
     * 0.0000,  0.0000,  0.0000, 0.0000, 
     * 0.0000,  0.0000,  0.0000, 0.0000, 0.000, 
C  MEDIUM DIESELS (CLS 6 & 7) 
     * 0.0000,  1.0000,  1.0000, 0.2000, 0.2000, 
     * 0.2000,  0.2000, 0.0385, 0.0247, 
     * 0.3516,  0.6100, 0.6100, 0.0000, 
     * 0.0000,  0.0000, 0.0000, 0.0000,  0.0000, 
C  8A HEAVY-DUTIES 
     * 0.0000,  1.0000,  1.0000, 0.3651, 0.3651, 
     * 0.3651,  0.3651, 0.3592, 0.3399, 
     * 0.4022,  0.3704, 0.4080, 0.4392, 
     * 0.5095,  0.5093, 0.6600, 0.6600,  0.0000, 
C  8B HEAVY-DUTIES 
     * 0.0000,  1.0000,  1.0000, 0.3949, 0.3949, 
     * 0.3949,  0.3949, 0.3881, 0.3511, 
     * 0.4130,  0.3825, 0.4103, 0.4384, 
     * 0.5188,  0.4810, 0.6600, 0.6600,  0.0000, 
C  DIESEL BUSES 
     * 0.0000,  1.0000,  1.0000, 0.0000, 0.000, 
     * 0.0000,  0.0000,  0.0000, 0.0000, 
     * 0.0000,  0.0000,  0.0000, 0.0000, 
     * 0.0000,  0.0000,  0.0000, 0.0000, 0.000/),(/18,5/)) 
C 
      REAL, DIMENSION(18,5) :: EXPRES=RESHAPE((/ 
C  LIGHT DIESEL 
     * 0.0000,  1.0000,  1.0000, 0.0000, 0.000, 
     * 0.0000,  0.0000,  0.0000, 0.0000, 
     * 0.0000,  0.0000,  0.0000, 0.0000, 
     * 0.0000,  0.0000,  0.0000, 0.0000, 0.000, 
C  MEDIUM DIESELS (CLS 6 & 7) 
     * 0.0000,  1.0000,  1.0000, 0.8220, 0.8220, 
     * 0.8220,  0.8220, 0.5026, 0.5396, 
     * 0.9249,  0.9800, 0.9800, 0.0000, 
     * 0.0000,  0.0000, 0.0000, 0.0000,  0.0000, 
C  8A HEAVY-DUTIES 
     * 0.0000,  1.0000,  1.0000, 0.9225, 0.9225, 
     * 0.9225,  0.9225, 0.9154, 0.9153, 
     * 0.9439,  0.9357, 0.9510, 0.9568, 
     * 0.9023,  0.8925, 0.9800, 0.9800,  0.0000, 
C  8B HEAVY-DUTIES 
     * 0.0000,  1.0000,  1.0000, 0.9377, 0.9377, 
     * 0.9377,  0.9377, 0.9373, 0.9285, 
     * 0.9533,  0.9478, 0.9625, 0.9697, 
     * 0.9064,  0.8756, 0.9800, 0.9800,  0.0000, 
C  DIESEL BUSES 
     * 0.0000,  1.0000,  1.0000, 0.0000, 0.000, 
     * 0.0000,  0.0000,  0.0000, 0.0000, 
     * 0.0000,  0.0000,  0.0000, 0.0000, 
     * 0.0000,  0.0000,  0.0000, 0.0000, 0.000/),(/18,5/)) 
C 
C 
      REAL, DIMENSION(18,5) :: FLEET=RESHAPE((/ 
C  LIGHT DIESELS 
     * 0.0000,  0.0000,  0.0000, 0.0000, 0.0000, 
     * 0.0000,  0.0000,  0.0000, 0.0000, 
     * 0.0000,  0.0000,  0.0000, 0.0000, 
     * 0.0000,  0.0000,  0.0000, 0.0000, 0.0000, 
C  MEDIUM DIESELS (CLS 6 & 7) 
     * 0.0000,  0.4000,  0.4000, 1.0000, 1.0000, 
     * 1.0000,  1.0000,  0.4950, 0.2470, 
     * 0.0592,  0.0120,  0.0060, 0.0000, 
     * 0.0000,  0.0000,  0.0000, 0.0000, 0.0000, 
C  8A HEAVY-DUTIES 
     * 0.0000,  1.0000,  1.0000, 1.0000, 1.0000, 
     * 1.0000,  1.0000,  1.0000, 1.0000, 
     * 1.0000,  1.0000,  0.7590, 0.5780, 
     * 0.2900,  0.2110,  0.1860, 0.1300, 0.0000, 
C  8B HEAVY-DUTIES 
     * 0.0000,  1.0000,  1.0000, 1.0000, 1.0000, 
     * 1.0000,  1.0000,  1.0000, 1.0000, 
     * 1.0000,  1.0000,  1.0000, 0.9940, 
     * 0.4680,  0.3400,  0.1510, 0.0880, 0.0000, 
C  DIESEL BUSES 
     * 0.0000,  0.7500,  0.7500, 0.0000, 0.0000, 
     * 0.0000,  0.0000,  0.0000, 0.0000, 
     * 0.0000,  0.0000,  0.0000, 0.0000, 
     * 0.0000,  0.0000,  0.0000, 0.0000, 0.0000/),(/18,5/)) 
C 
      END MODULE DDDATA 
