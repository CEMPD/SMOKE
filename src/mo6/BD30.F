      BLOCK DATA BD30 
C 
C  BLOCK DATA Subprogram 30: /CITPAR/,/CITRV1/,/CITRV2/, /INJECT/, 
C                            /RFORM1/, /RFORM2/ 
C 
C  BLOCK DATA 30 defines the parameters for reading in, storing and using 
C  the second required scenario record (local data) parameters, plus fuel 
C  injection rates and model year groups. 
C 
C  Changes: (Last change first) 
C 
C  09 May 01 @EPA-BAG Added fractions and additional model years for post 2004 fuel injection vehicle percentages 
C 06 Oct 00 @EPA-BAG  bug155 Removed INDUSTRY AVERAGE FUEL parameters IAFCNT 
C           and GLIAFCNT 
C  08 Aug 00 @EPA-eg  bug216; added GLRFG 
C  21 Jul 00 @EPA-RAG running loss; removed pass/fail fractions, PFRATE's, 
C            for purge/pressure evap. tests 
C  26 Jun 00 @EPA-elg Changed RFGRVP array to (2,18) 
C  12 Jun 00 @EPA-djb bug#99  Removed RVPX and added RVPNOX and OXYRVP. 
C  17 May 00 AIR Task 03: Removed hourly dependency of RVPHS. 
C  01 Jun 00 @EPA-lap Bug 86, corrected fuel injection rates for HDG in 
C            DATA PFI, were all zero. 
C  16 May 00 @EPA-ddj RFG BUG, Removed RFG_S array. 
C  02 May 00 @EPA-ddj D-28,  Added data to arrays RFGRVP and GLRFGRVP 
C            to account for additional RFG configurations. 
C  28 Mar 00 @DynTel-rwk 1-028, Added global RFG RVP variable for 
C            reinitialization 
C  25 Feb 00 @DynTel-ddj 1-028, Changed DATA for RFGSYR and RFGRVP to 
C            reflect change in dimension of arrays. 
C  04 Apr 00 AIR Task 02: Removed ASTM class info. 
C  04 Apr 00 @EPA-djb   Updated TBI arrays to match MOBILE5b for 
C            pre-1990. 
C  04 Dec 98 @DynTel-bag 2-664  Removed uncontrolled diurnal arrays UDI, 
C            and IUDI. Added hourly arrays to hot soak RVP, RVPHS(14). 
c            Temporarily, exhaust RVP is a hourly temperature weighted average. 
C            (See the calculations in LOCAL.) Deleted RVPLOS. 
C  11 Aug 98 @DynTel-MLA 2-000  Changed HDG and HDD names to conform to 
C            new naming conventions. 
C 
C  Common block array subscripts: 
C 
C  PFI(13,14)    -  PFI ( IG5, IVTAM ) 
C  RFGSYR(2)     -  RFGSYR( JPHASE ) 
C  RFGRVP(2,5)   -  RFGRVP( REGION, IYR ) 
C  GLRFGRVP(2,5) -  GLRFGRVP( REGION, IYR ) 
C  SCNAME        -  SCNAME ( ICH ) 
C  TBI(13,14)    -  TBI ( IG5, IVTAM ) 
C 
C  Common block dictionary: 
C 
C   Name   Type              Dictionary 
C  ------  ----  ------------------------------------------------ 
C 
C  /CITPAR/: 
C 
C  SCNAME   C    scenario identification field: echoed in numeric formats 
C 
C  /CITRV2/: 
C 
C  RVPHS    R    RVP in untampered/uncontrolled hot soak 
C  RVP090   R    RVP in RVPCER cases prior to ICERSY, in denominator of Diurnal 
C                Index equation & to normalize RVPCER. 
C  RVP115   R    Uncontrolled RVP. 
C  RVPNOX   R    Weathered 8.7 psi RVP value. 
C  OXYRVP   R    Weathered RVP value without oxy commingling effects. 
C 
C  /INJECT/: 
C 
C  TBI      I    throttle body fuel injected fraction of IVTAM built in IG5 
C  PFI      I    port fuel injected fraction of veh class IVTAM built in IG5 
C  MAXFIY   I    first year of highest myg for indexing TBI & PFI. 
C  MINFIY   I    last year of lowest myg for indexing TBI & PFI. 
C 
C  /RFORM1/: 
C 
C  JPHASE   I    corresponds to the reformulated fuel phase 
C                1 = Phase I 
C                2 = Phase II 
C  RFGFLG   I    flag indicating a reformulated gas area 
C                1 = no reform program 
C                2 = apply reform program 
C                3 =  apply reform program phase I only 
C  RFGON    L    flag indicating a reformulated gas area with a calender 
C                year of evaluation that is >= reform phase year 
C  WINFLG   I    Corresponds to the season that was input by the user. 
C                (This will override the value of MEVAL for the 
C                 reformulated fuel calculations.) 
C                  0 = No override desired 
C                  1 = Summer 
C                  2 = Winter 
C  SEAFLG   I    Corresponds to the season that is either supplied by 
C                the user in WINFLG, or calculated from the contents 
C                of MEVAL. If WINFLG is 1 or 2, it has priority. 
C                  1 = Summer 
C                  2 = Winter 
C 
C  /RFORM2/: 
C 
C  GLRFGRVP   R    reformulated RVP 
C  RFGSYR   I    start year for phase of reformulated fuel effects 
C  RFGRVP   R    reformulated RVP 
C  RFG_CF   R    Array holding the residual effects of RFG after 
C                sulfur, RVP and Oxy fuel effects have been applied. 
C 
C 
C  Notes: 
C 
C  RVP and gas fraction values for oxygenated fuels support from MOBILE4 
C  were removed in MOBILE4.1.  Certification fuels switch also removed. 
C 
C  TMM@ARC 12/15/93 Request #2-366. 
C  Added comment about reset of RFGSYR in RDCITY. 
C  Oct-20-1994 @ CSC-bsg Request 426  Summer RVP values for Phase 2 
C  for Reformulated Gasoline have changed from the values found in MOBILE5a. 
C 
C  3 Dec 96 DynTel@MLA 2-622 Expanded the arrays to accomodate the new 
C           Heavy Duty Gasoline, LDGT, and GAS BUS vehicle categories. 
C  Oct-21-97 DynTel@gjr 2-645 Removed RVPIUS and IUSESY. 
C  March-2-98 DynTel@ID 2-643 Updated LDGV and LDGT technologies. Commented 
C  data for CARB and OPLP technologies (will be added eventually). 
C 18 May 98 DynTel@MLA 2-659 Changed SCNAME from REAL to CHAR data type. 
C 
C 23 Jul 98 @DynTel-MLA 2-000  Changed LDGT names to conform to 
C           certification naming conventions. 
C 
      IMPLICIT NONE 
      INCLUDE 'CITPAR.I' 
      INCLUDE 'CITRV1.I' 
      INCLUDE 'CITRV2.I' 
      INCLUDE 'INJECT.I' 
      INCLUDE 'RFORM1.I' 
      INCLUDE 'RFORM2.I' 
C 
C 
C  The /CIT???/ variables / arrays get their values by being: 
C 
C  Assigned in LOCAL/RDCITY as follows: 
C 
C  supplied by user: SCNAME,RVPBAS 
C 
C  calculated by code: RVPICY,RVPHS 
C 
C  The only constants are: 
C  the controlled and uncontrolled RVP standard values of /CITRV2/: 
C 
      DATA RVP090/9.0/,RVP115/11.5/ 
C 
C  Dummy zeroes are assigned the remaining variables and arrays. 
C 
      DATA SCNAME/'                '/ 
      DATA RVPBAS,RVPICY,RVPUWX/3*0.0/ 
      DATA RVPHS,RVPNOX,OXYRVP/3*0.0/ 
      DATA RFG_CF /540*1.0/ 
      DATA GLRFG /1/ 
C 
C  The reformulated fuel variables. 
C 
      DATA JPHASE,RFGFLG,WINFLG,SEAFLG,RFGON 
     *     /1, 1, 1, 1, .FALSE. / 
C 
C  Start Year by phase; RVP by phase. 
C 
C  RFGSYR(3) is reset to 2000 in RDCITY because it might have been 
C  modified there;  12/15/93 TMM@ARC Request #2-366. 
C  Oct-20-1994 @ CSC-bsg Request 426  Summer RVP values for Phase 2 
C  for Reformulated Gasoline have changed from the values found in MOBILE5a. 
C 
C 
      DATA RFGSYR, GLRFGSYR/ 1995, 2000, 1995, 2000/ 
C 
C  RFGRVP is the RVP corresponding to the RFG region, Phase and calendar year. 
C  The first two numbers correspond to Phase 1 (1995-1999) for regions 
C  South (7.1) and North (8.0).  The remaining values are for Phase 2 starting in 
C  calendar year 2000 through 2015+. 
C 
      DATA RFGRVP  / 7.1, 8.0, 6*6.7, 28*6.8/ 
ccc      DATA GLRFGRVP/ 7.1, 8.0, 6*6.7, 2*6.8/ 
C 
C  /INJECT/: no provision exists for user modification of the fuel injection 
C            fractions or their model year groups. 
C 
      DATA MAXFIY/1995/,MINFIY/1980/ 
C 
C  Both TBI and PFI have the old Mobile5 IVTAM continuation code: 
C    1 = LDGV 
C    2 = LDGT12 
C    3 = LDGT34 
C    4 = HDGV2B, HDGV3, HDGV4, ..., HDGV8B 
C    B = GAS BUS (new in Mobile6) 
C  Each line has 8 (model year group) values 
C 
C  Observe: 
C 
C    1) LDGT12 = LDGT34 values 
C    2) All HDGV values = 0.0/0, since for now HDGV's are assumed to 
C       be 100% carbureted. 
C    3) For tampering, the fractions are only used to weight hot soak rates. 
C    4) Before 1981, it is assumed that no vehicle of any type had 
C       fuel injection and therefore hot soak rates are for those years 
C       based solely on evaporative cannister only tampering. 
C 
C    pre-81 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992  
C      1993 1994 1995 2004+ 
C 
      DATA TBI/ 
     1 .000,.029,.106,.183,.282,.208,.284,.369,.407,.275,.188,.208,.096, 
     1 .109,.039,.012,.000, 
     2 .000,.000,.000,.000,.000,.047,.136,.280,.452,.369,.394,.495,.319, 
     2 .301,.272,.254,.000, 
     2 .000,.000,.000,.000,.000,.047,.136,.280,.452,.369,.394,.495,.319, 
     2 .301,.272,.254,.000, 
     3 .000,.000,.000,.000,.000,.047,.136,.280,.452,.369,.394,.495,.319, 
     3 .301,.272,.254,.000, 
     3 .000,.000,.000,.000,.000,.047,.136,.280,.452,.369,.394,.495,.319, 
     3 .301,.272,.254,.000, 
     4 .000,.000,.000,.000,.000,.000,.000,.000,.000,.000,.000,.000,.000, 
     4 .000,.000,.000,.000, 
     4 .000,.000,.000,.000,.000,.000,.000,.000,.000,.000,.000,.000,.000, 
     4 .000,.000,.000,.000, 
     4 .000,.000,.000,.000,.000,.000,.000,.000,.000,.000,.000,.000,.000, 
     4 .000,.000,.000,.000, 
     4 .000,.000,.000,.000,.000,.000,.000,.000,.000,.000,.000,.000,.000, 
     4 .000,.000,.000,.000, 
     4 .000,.000,.000,.000,.000,.000,.000,.000,.000,.000,.000,.000,.000, 
     4 .000,.000,.000,.000, 
     4 .000,.000,.000,.000,.000,.000,.000,.000,.000,.000,.000,.000,.000, 
     4 .000,.000,.000,.000, 
     4 .000,.000,.000,.000,.000,.000,.000,.000,.000,.000,.000,.000,.000, 
     4 .000,.000,.000,.000, 
     4 .000,.000,.000,.000,.000,.000,.000,.000,.000,.000,.000,.000,.000, 
     4 .000,.000,.000,.000, 
     B .000,.000,.000,.000,.000,.000,.000,.000,.000,.000,.000,.000,.000, 
     B .000,.000,.000,.000/ 
C 
      DATA PFI/ 
     1 .000,.061,.062,.088,.110,.307,.392,.372,.492,.597,.793,.790,.902, 
     1 .891,.961,.988,1.00, 
     2 .000,.000,.000,.002,.022,.066,.238,.327,.416,.540,.552,.484,.660, 
     2 .688,.728,.746,1.00, 
     2 .000,.000,.000,.002,.022,.066,.238,.327,.416,.540,.552,.484,.660, 
     2 .688,.728,.746,1.00, 
     3 .000,.000,.000,.002,.022,.066,.238,.327,.416,.540,.552,.484,.660, 
     3 .688,.728,.746,1.00, 
     3 .000,.000,.000,.002,.022,.066,.238,.327,.416,.540,.552,.484,.660, 
     3 .688,.728,.746,1.00, 
     4 .000,.000,.000,.000,.000,.000,.000,.892,.953,1.00,1.00,1.00,1.00, 
     4 1.00,1.00,1.00,1.00, 
     4 .000,.000,.000,.000,.000,.000,.000,.892,.953,1.00,1.00,1.00,1.00, 
     4 1.00,1.00,1.00,1.00, 
     4 .000,.000,.000,.000,.000,.000,.000,.350,.354,.691,.700,1.00,1.00, 
     4 1.00,1.00,1.00,1.00, 
     4 .000,.000,.000,.000,.000,.000,.000,.350,.354,.691,.700,1.00,1.00, 
     4 1.00,1.00,1.00,1.00, 
     4 .000,.000,.000,.000,.000,.000,.000,.350,.354,.691,.700,1.00,1.00, 
     4 1.00,1.00,1.00,1.00, 
     4 .000,.000,.000,.000,.000,.000,.000,.350,.354,.691,.700,1.00,1.00, 
     4 1.00,1.00,1.00,1.00, 
     4 .000,.000,.000,.000,.000,.000,.000,.350,.354,.691,.700,1.00,1.00, 
     4 1.00,1.00,1.00,1.00, 
     4 .000,.000,.000,.000,.000,.000,.000,.350,.354,.691,.700,1.00,1.00, 
     4 1.00,1.00,1.00,1.00, 
     B .000,.000,.000,.000,.000,.000,.000,.350,.354,.691,.700,1.00,1.00, 
     B 1.00,1.00,1.00,1.00/ 
C 
       END 
