C 
C This file contains variables that control the expanded output options -- 
C printing emission tables for the vehicle types that are new in Mobile6 
C and the production of output in the spreadsheet and descriptive formats. 
C 
C  Changes: (Last change first) 
C 
C  15 May 01 @EPA-BAG Bug120, Removed SSFLAG and SSONLY. 
C  10 Apr 00 - AIR Task 02: Added EXPEXH and NODESC 
C  20 Apr 98 - DynTel@BAG  Added SSONLY variable for spreadsheet output. 
C  3 Apr 97 - DynTel@MLA  Expanded array OPTABL to accomodate the new BUS  
C             and CNG vehicle categories. 
C 
C 
C OPTABL  -  I  -  This array defines whether expanded output tables will  
C                  be printed for LDGT, HDGV, HDDV, BUS, or CNG. 
C EXPEXH  -  L  -  Should descriptive output include expanded exhaust? 
C NODESC  -  L  -  Should descriptive output be skipped? 
C 
C The data is initialized in BD07. 
C 
      LOGICAL   EXPEXH, NODESC 
      INTEGER   OPTABL(4) 
      COMMON   /OPFMT/ OPTABL, EXPEXH, NODESC 
 
