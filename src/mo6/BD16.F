      BLOCK DATA BD16 
C 
C  BLOCK DATA Subprogram 16: /IOUCOM/,/IVPCOM/,/MAXIMA/,/QUITXQ/,/YEARS4/ 
C 
C  BLOCK DATA 16 initializes miscellaneous integer-valued parameters. 
C 
C  Changes : 
C 
C  31 Aug 01 AIR Task 22-24: Added IOUATSS 
C  11 Jul 01 AIR Task 15: added I/O for Air Toxic output. 
C  10 Apr 11 AIR Task 11: Redirected IOUERR to separate file. 
C  09 Mar 01 AIR Task 11: Added I/O for PM 
C  02 Oct 00 @EPA-djb Bug230 Deleted some I/O variables. 
C  18 Sep 00 @EPA-elg  Fixed Bug 230.  IOU variables are initialized consistently 
C                      with DRIVER. 
C  08 May 00 AIR Task 03: Added IVPTRR array to handle BIGRVP calls. 
C   8 Jun 00 @EPA-ddj, Bug120 Removed code associated with spreadsheet 
C            output. 
C  11 Apr 00 AIR Task 02: Replaced all occurences IOUASK with IOUIN and IOUOUT 
C            so that Lahey Standard for keyboard/monitor I/O is followed. 
C   2 Feb 00 @DynTel-ddj 1-032  Removed commented out FTP1-3 coding. 
C            Added IVPTRC for mapping TIER2_50K array. 
C  27 Jan 00 @DynT3l-MLA 1-011 Added index 22 to IVPTRB for LDDT34. 
C  29 Nov 99 @DynTel-bag 1-006 removed variables for output of hourly evaporative 
C            emissions; i.e., those variables in HOURLY.I 
C  11 Aug 99 @DynTel-MLA 2-684 Added comment about reserved unit numbers 
C            for the database output tables. 
C  18 Jun 99 @DynTel-ddj 2-000 Added FTP1, FTP2, and FTP3 to allow writing 
C            out of FTP emission #'s; these are currently commented out 
C  06 Apr 99 @DynTel-bag 2-000 MAXIH_EVAP replaced with MAXIH 
C  12 Mar 99 @DynTel-bag 2-664  IWEEK is now set by a labelled input 
C            and is passed through WEEKLY.I(BD16.FOR). It is initialized 
C            (default value) to weekday (IWEEK = 1) for each scenario in 
C            INISCN.FOR . For temporary output of hourly emissions added 
C            HOURLY.I with the variables HOURLY_OUT. 
C    4 Mar 99 @DynTel-bag 2-664  MAXIH_EVAP now 24. (It was 18.) 
C   11 Feb 99 @DynTel-MLA 2-663  Added unit IOUXIM to read external I/M 
C             descriptions. 
C   30 Oct 98 @DynTel-bag 2-677 Added MAXIH_EVAP evaporative hourly loop. 
C   22 Oct. 98 @DynTel-HXQ 2-680 delete the last number of IVPTRB. 
C   16 Sep 98 @DynTel-MLA 2-659 Added IVPTRS for hourly loop 
C   10 Sep 98 @DynTel-BAG 2-661 WMOUT now 80 characters long and added 
C             logical variable READOUTFILE 
C 
C   03 Sep 98 @DynTel-RJD 0-632 Added new MCYRS variable that specifies the 
C            maximum model year index for motorcycles included in basic 
C            emission factor calculations. 
C 
C  Common block array subscripts: 
C 
C  IVPTRL(28)  -  IVPTRL ( IVLM ) 
C  IVPTRT(28)  -  IVPTRT ( IVTEMP ) 
C  IVPTRA(28)  -  IVPTRA ( IVA ) 
C  IVPTRB(28)  -  IVPTRB ( IVB ) 
C  IVPTRC(28)  -  IVPTRC ( IV  ) 
C  IVPTRS(28)  -  IVPTRS ( IVST ) 
C 
C  Common block dictionary: 
C 
C  /IOUCOM/: a special format is utilized for the i/o unit number (#) 
C          COMMON block: 
C 
C  READOUTFILE L determines whether or not the user supplied filename 
C                and filename extension of output files (.TRUE. if 
C                user supplied output filename with extension; .FALSE. 
C                if WMOUT does not contain a filename extension); 
C                see DRIVER.FOR for assignment of READOUTFILE 
C  WMOUT   C*80  name of the user supplied output file 
C 
C  Unit#   Usage                     Description 
C  ------  ------  ---------------------------------------------------------- 
C  IOUALT  input   user supplied alternate (external) data file 
C  IOUBAT  input   user supplied batch input data 
C  IOUGEN  input   user supplied general input data 
C  IOUIMD  input   user supplied I/M credits (when IMFLAG=2) 
C  IOUIN   input   input device (keyboard) 
C  IOUREP  output  formatted report file 
C  IOUERR  output  error/warning messages file 
C  IOUOUT  output  prompting messages device (monitor) 
C  IOUXIM  input   user supplied I/M parameters from an external file 
C  IOUNULL output  Dummy IO unit - output written here disappears! 
C  NOTE: Unit numbers 51-57 are reserved for the database ouptut tables. 
C 
C   Name   Type               Description 
C  ------  ----  -------------------------------------------- 
C  /IVPCOM/: 
C  IVPTRL   I    IV subgroup pointer into MYGSP1 & BFRCOF (see IVLM subscript 
C                definition) 
C  IVPTRT   I    IV subgroup pointer into MYGTCF (see IVTEMP subscript 
C                definition) 
C  IVPTRA   I    IV subgroup pointer into ALHRET (see IVA 
C                subscript definition) 
C  IVPTRB   I    IV subgroup pointer into GPBSCO (see IVB subscript 
C                definition) 
C  IVPTRC   I    IV subgroup pointer into TIER2_50K array. 
C  IVPTRS   I    Maps the vehicles that have Tier2 emissions 
C 
C  /QUITXQ/: 
C 
C 
C  N1QUIT   I    switch: 1 = first QUITER call has been made, 0 = it has not 
C                On first call, extra vertical spacing is inserted between the 
C                run title and the QUITER message. 
C 
C  /MAXIMA/: 
C 
C  MAXPOL     I    maximum number of pollutants covered by MOBILE4 
C                ( = 3 = HC, CO & NOx ) 
C  MAXREG     I    maximum number of regions allowed in MOBILE4 
C                ( = 3 = low, high & mid ) 
C  MAXYRS     I    maximum number of model years in ICY's fleet (1 thru 20+) 
C  MAXIH      I    maximum number of time intervals in a day 
C  MAXEMIT    I    maximum number of emitter categories (high/normal) 
C 
C  /YEARS4/: 
C 
C  IYSTART  I    4-digit representation of start year 
C  IY1960   I    4-digit representation of year 1960 
C  IY1994   I    4-digit representation of year 1994 
C  IYEND    I    4-digit representation of end year 
C 
C  Notes: 
C 
C 15 Nov 96 DynTel@MLA 2-622  Expanded the arrays to accomodate the new  
C           Heavy Duty Gasoline, Heavy Duty Diesel, and LDGT vehicle 
C           categories. 
C 
C  7 Apr 97 DynTel@MLA 2-622  Expanded arrays IVPTRA, IVPTRB, IVPTRL, and  
C           IVPTRT to accomodate the new BUS vehicle types. Constants 
C           MAXLTW and MAXVEH were removed from this subroutine. MAXLTW is 
C           obsolete in Mobile6. MAXVEH is now defined in BD06, along with  
C           the vehicle type arrays and symbolic names for vehicle types. 
C 
C   25 Mar 98 @DynTel-MLA 2-659 Added IOUALT for Mobile6. 
C   24 Apr 98 @DynTel-BAG 2-661 Added spreadsheet output file unit number 
C             variable (IOUSSH=15) for Mobile6. 
C 
C 
      IMPLICIT NONE 
C 
      INCLUDE 'IOUCOM.I' 
      INCLUDE 'IVPCOM.I' 
      INCLUDE 'MAXIMA.I' 
      INCLUDE 'PARAME.I' 
      INCLUDE 'PARAMS.I' 
      INCLUDE 'QUITXQ.I' 
      INCLUDE 'WEEKLY.I' 
      INCLUDE 'YEARS4.I' 
C 
      DATA WMOUT/ ' '/ 
C 
C     NOTE: Unit numbers 51-57 are reserved for the database output 
C           tables and are initialized in MODULE DATABASE. 
      DATA IOUGEN,IOUIMD,IOUIN,IOUALT,IOUBAT,IOUXIM/1, 4, 5, 7, 19, 3/ 
      DATA IOUREP,IOUERR,IOUOUT,IOUNULL/2, 13, 6, -1/ 
      DATA IOUPM,IOUAT,IOUATSS/8,9,10/ 
C 
      DATA IVPTRL/1, 2, 3, 4, 5, 0, 0, 0, 0, 0, 0, 0, 0, 6, 7, 0, 0, 
     &            0, 0, 0, 0, 0, 0, 8, 0, 0, 0, 9/ 
C 
      DATA IVPTRT/1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13, 0, 0, 0, 0, 
     &            0, 0, 0, 0, 0, 0,14,15, 0, 0, 0/ 
C 
      DATA IVPTRA/1, 2, 3, 4, 5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
     &            0, 0, 0, 0, 0, 0, 6, 0, 0, 0, 0/ 
C 
      DATA IVPTRB/0, 0, 0, 0, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12, 
     &           13,14,15,16,17,18, 0,19,20,21,22/ 
C 
      DATA IVPTRC/1, 2, 3, 4, 5, 0, 0, 0, 0, 0, 0, 0, 0, 1, 3, 0, 0, 
     &            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 5/ 
C 
      DATA IVPTRS/1, 2, 3, 4, 5, 0, 0, 0, 0, 0, 0, 0, 0, 6, 7, 0, 0, 
     &            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 8/ 
C 
      DATA IVPTRR/1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13, 0, 0, 0, 0, 
     &            0, 0, 0, 0, 0, 0, 0,14, 0, 0, 0/ 
     
C 
      DATA MAXPOL,MAXREG,MAXYRS,MCYRS/3,2,25,12/ 
      DATA MAXIH,MAXEMIT/24,2/ 
C 
      DATA N1QUIT/1/ 
C 
      DATA IWEEK/ 1 / 
C 
      DATA IYSTART,IY1960,IY1994,IYEND/NSYR,1960,1994,NEYR/ 
C 
      END 
