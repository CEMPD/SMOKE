      MODULE DATABASE 
C 
C     This module defines the data structures that are used by the  
C     code in Mobile that produces database output. 
C 
C  Changes : (Last change first) 
C 
C  13 Dec 01 AIR Task 27: Changed CO2 to be not selected by default. 
C  06 NOV 01 AIR Task 27: Added CO2 support to DBSELPOL and DBPOLNAMES. 
C            DBVEHNAMES now includes ALL VEH. 
C  01 Aug 01 AIR Task 15: Added DBSELTOX for Air Toxics. 
C  10 May 01 AIR Task 11/5: Extended DBSELPART to handle idle PM. 
C  06 Mar 01 AIR Task 11: Added particulate elements to DBSELEFT. 
C  30 Jan 01 AIR Task 11: Added particulate and air toxic elements. 
C  24 Oct 00 @EPA-BG  New labels for database age, model year and hour selection 
C  19 Oct 00 @EPA-djb Task W2, Added DBAGGR and DBAGGRCOLS array. 
C  20 Oct 00 AIR Task 03: Added MYR to daily and hourly output 
C  05 MAY 00 @EPA - HMC E22 ADDED DBHOURCOLS, DBDAYCOLS, DBHEAD, AND DBDAILY 
C  24 Jan 00 @DynTel-MLA 1-011  Expanded array DBSELVEH to 28 places. 
C 
C  Notes: 
C 
C 
      IMPLICIT NONE 
C 
      INTEGER                 ::  DBFLAG  ! DB output requested? 
      LOGICAL                 ::  DBAGGR  ! AGGREGATED OUTPUT REQUESTED? 
      LOGICAL                 ::  DBDAILY ! DAILY OUTPUT REQUESTED? 
      LOGICAL                 ::  DBHEAD  ! FIRST RECORD OF FIELDNAMES? 
C 
      INTEGER, DIMENSION(2)   ::  DBMODELYR 
      INTEGER, DIMENSION(2)   ::  DBSELAGE = (/0,24/) 
      INTEGER, DIMENSION(2)   ::  DBSELHR  = (/1,24/) 
      INTEGER, DIMENSION(10)  ::  DBSELEFT = (/2,2,2,2,2,2,2,2,2,2/) 
      INTEGER, DIMENSION(5)   ::  DBSELFAC = (/2,2,2,2,2/) 
      INTEGER, DIMENSION(4)   ::  DBSELPOL = (/2,2,2,1/)  ! CO2 is not selected by default 
      INTEGER, DIMENSION(28)  ::  DBSELVEH = (/2,2,2,2,2,2,2,2, 
     &               2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2/) 
      INTEGER, DIMENSION(11)  ::  DBSELPART= (/1,1,1,1,1,1,1,1,1,1,1/) 
      INTEGER, DIMENSION(6)   ::  DBSELTOX = (/1,1,1,1,1,1/) 
C 
      CHARACTER, PARAMETER        :: TABCHAR = CHAR(9) 
      LOGICAL                     :: DBOPEN = .FALSE. 
      INTEGER, DIMENSION(7), SAVE :: DBUNIT = (/51,52,53,54,55,56,57/) 
      CHARACTER (LEN=80)          :: DBBATCH 
      CHARACTER (LEN=80)          :: DBINPUT 
C 
      INTEGER                     :: DB_MDLYR = 0 !database output selected by model year; 
                                                  !default value is all models years 
      INTEGER                     :: DBNFILE 
      INTEGER                     :: DBNRUN 
      INTEGER                     :: DBNSCEN 
C 
      CHARACTER(LEN=10), DIMENSION(7), PARAMETER ::  DBTABLES = 
     &  (/'EMISSIONS ','INPUTS    ','DATAFILES ','MESSAGES  ', 
     &    'RUNDATA   ','SCENARIOS ','ETYPES    '/) 
C 
      CHARACTER(LEN=80), DIMENSION(7)            ::  DBFNAMES  
C     INTEGER, DIMENSION(7)                      ::  DBOPENCOUNT = 
C    &   (/0,0,0,0,0,0,0/) 
      CHARACTER(LEN=7), DIMENSION(7)             ::  DBACTION = 
     &   (/'REPLACE','REPLACE','REPLACE','REPLACE','REPLACE', 
     &     'REPLACE','REPLACE'/) 
C 
      CHARACTER(LEN=11), DIMENSION(7), PARAMETER :: DBDEFAULTNAMES = 
     &   (/'table01.txt','table02.txt','table03.txt','table04.txt', 
     &     'table05.txt','table06.txt','table07.txt'/) 
C 
      CHARACTER(LEN=16), DIMENSION(9), PARAMETER ::  DBEFTNAMES = 
     &   (/'Exhaust Running ','Exhaust Start   ','Hot Soak        ', 
     &     'Diurnal         ','Resting Loss    ','Running Loss    ', 
     &     'Crankcase       ','Refueling Loss  ','Exhaust Idle    '/) 
      CHARACTER(LEN=4), DIMENSION(5), PARAMETER  ::  DBFACNAMES = 
     &   (/'Fwy.','Art.','Loc.','Ramp','None'/) 
      CHARACTER(LEN=4), DIMENSION(4), PARAMETER  ::  DBPOLNAMES = 
     &   (/'HC  ','CO  ','NOX ','CO2 '/) 
      CHARACTER(LEN=8), DIMENSION(29), PARAMETER ::  DBVEHNAMES = 
     &   (/'LDGV    ','LDGT1   ','LDGT2   ','LDGT3   ','LDGT4   ', 
     &     'HDGV2B  ','HDGV3   ','HDGV4   ','HDGV5   ','HDGV6   ', 
     &     'HDGV7   ','HDGV8A  ','HDGV8B  ','LDDV    ','LDDT12  ', 
     &     'HDDV2B  ','HDDV3   ','HDDV4   ','HDDV5   ','HDDV6   ', 
     &     'HDDV7   ','HDDV8A  ','HDDV8B  ','MC      ','GAS BUS ', 
     &     'URB BUS ','COM BUS ','LDDT34  ','ALL VEH '/) 
C  OPTIONAL COLUMN HEADERS FOR HOURLY DATABASE OUTPUT 
      CHARACTER(LEN=8), DIMENSION(22)            ::  DBHOURCOLS = 
     &   (/'FILE    ','RUN     ','SCEN    ','POL     ','VTYPE   ', 
     &     'ETYPE   ','FTYPE   ','AGE     ','HOUR    ','GM_MILE ', 
     &     'GM_HOUR ','STARTS  ','ENDS    ','MILES   ','MPG     ', 
     &     'HRVMT   ','FACVMT  ','REGDIST ','VCOUNT  ','AMBTEMP ', 
     &     'DIURTEMP','MYR     '/) 
C  OPTIONAL COLUMN HEADERS FOR DAILY DATABASE OUTPUT 
      CHARACTER(LEN=8), DIMENSION(18)            ::  DBDAYCOLS = 
     &   (/'FILE    ','RUN     ','SCEN    ','POL     ','VTYPE   ', 
     &     'ETYPE   ','FTYPE   ','AGE     ','GM_MILE ', 
     &     'GM_DAY  ','STARTS  ','ENDS    ','MILES   ','MPG     ', 
     &     'FACVMT  ','REGDIST ','VCOUNT  ','MYR     '/) 
C 
C  OPTIONAL COLUMN HEADERS FOR COMPOSITE DATABASE OUTPUT 
      CHARACTER(LEN=8), DIMENSION(13)            ::  DBAGGRCOLS = 
     &   (/'FILE    ','RUN     ','SCEN    ','CAL_YEAR','POL     ', 
     &     'VTYPE   ','GM_MILE ','GM_DAY  ','STARTS  ','ENDS    ', 
     &     'MILES   ','MPG     ','VMT     '/) 
C 
      INTEGER                                    ::  DBLINENUM 
      INTEGER                                    ::  DBSTRNUM 
C 
C   DBRUNFLAGS describes whether user supplied data for or specified: 
C     1. No tampering. 
C     2. Disable 1990 CAAA. 
C     3. EFs do not include refueling loss. 
C     4. Stage II parameters. 
C     5. Alternate tampering rates. 
C     6. Alternate registration distributions. 
C     7. Alternate mileage accumulation rates. 
C     8. Alternate basic emission rates. (81-93, 94+) 
C     9. I/M program parameters. 
C    10. ATP program parameters. 
C    11. Pressure/purge functional test parameters. 
C    12. New LEV implementation shcedule. 
C    13. Alternate start EFs. 
C    14. Alternate running EFs. 
 
      CHARACTER (LEN=1),DIMENSION(14)     ::  DBRUNFILES 
      CHARACTER (LEN=80),DIMENSION(14)    ::  DBRUNFLAGS 
 
      CHARACTER (LEN=1),DIMENSION(8)      ::  DBSCFLAGS 
 
      CHARACTER (LEN=150)                 ::  DBSCTITLE 
 
      END MODULE DATABASE 
