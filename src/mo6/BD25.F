      BLOCK DATA BD25 
C 
C  BLOCK DATA Subprogram 25: /ATPAR1/,/ATPAR2/,/ATPAR3/,/ATPAR4/,/EGSCAL/, 
C                            /EMECAL/,/EVPDAT/,/LOOKUP/,/MYCODE/,/MYUB1/ , 
C                            /MYUB2/ ,/OMTCOM/,/SIZCAL/,/SIZPAR/,/TAMEQ1/, 
C                            /TAMEQ2/,/TAMEQ3/,/TAMEQ4/,/TAMOUT/,/TAMOU2/, 
C                            /TAMPB1/ 
C 
C  BLOCK DATA 25 defines and initializes the common blocks created by the 
C  tampering code. 
C 
C  Changes: (Last change first) 
C 
C  07 Mar 01 @EPA-bag Bug317, fixed typo in HD technology fractions (TGS arrays) 
C            for HDG2b & 3, model years 1996-1999(IG1=4), EGR with 3-way 
C            catalyst group (IQG=5) (was 0.37 now 0.27); typo is from m6.flt.008 
C            in table 4-2 
C  04 Dec 00 @EPA-bag Task X6, removed MY, JDX, and IDX from MYCODE.I 
C  28 Jun 00 @EPA-ddj Bug158  Incorporated HD technology fractions as 
C            described in M6.FLT.008 for 1990 and later model years. 
C  20 Jun 00 AIR Task 03: Removed OMTTAM array 
C  20 Jun 00 @AIR-DK  T-03   replaced OMTCF with TCFSTART, TCFRUN, and TCF12 
C  12 Jun 00 @EPA-ddj Bug141 Corrected TGS typos for HDGV2b. 
C  23 May 00 @EPA-elg   Removed TEMP_ST array 
C  26 Apr 00 @EPA-elg  Added TEX and TAMBAS Arrays. 
C  31 Mar 99 @DynTel-DDJ 2-696  OMTCF array expanded for hourly intervals and added 
C            DIST_SUM 
C  6 Nov 98 @DynTel-MLA 2-663  A new dimension was added to ICYIM to allow 
C           the array to hold program end dates as well as start dates. All 
C           references to ICYIM in existing code were updated accordingly. 
C  11 Aug 98 @DynTel-MLA 2-000  Changed HDG and HDD names to conform to 
C            new naming conventions. 
C 
C  Common block array subscripts: 
C 
C  AER(11,11,2,2)    -  AER    ( IC11, IC11, IPG, IAY ) 
C  AIR(2)            -  AIR    (  IAY ) 
C  ATR(5,2)          -  ATR    ( IDBAG , IAY ) 
C  BTR(9,2)          -  BTR    ( ID, IAY ) 
C  CAN(2)            -  CAN    (  IAY ) 
C  CAP(2)            -  CAP    (  IAY ) 
C  CAT(2)            -  CAT    (  IAY ) 
C  CCEMI(7,14)       -  CCEMI  ( IG4, IVTAM ) 
C  CSAE(12)          -  CSAE   ( IC ) 
C  CSIZE(12,2)       -  CSIZE  ( IC , IAY ) 
C  DIST_SUM(24)      -  DIST_SUM (IH) 
C  DISTYP(8)         -  DISTYP  ( IDUSER )  (see ID defn in (1,999) ) 
C  EFF(2,3,2,15)     -  EFF    (  ATPPGM, ICOL , IATP , IROW ) 
C  EFFSAV(2,3,2,15)  -  EFFSAV (  ATPPGM, ICOL , IATP , IROW ) 
C  EGR(2)            -  EGR    (  IAY ) 
C  EGS(7,2)          -  EGS    ( IEG, IPG ) 
C  EMEGR(4,3,3)      -  EMEGR  ( IG2, IB, IP ) 
C  EMI(7,3,3)        -  EMI    ( IEG, IB, IP ) 
C  EVPAER(2,2)       -  EVPAER ( IHG, IAY ) 
C  EVPTGS(3,2,14)    -  EVPTGS ( IG3, IHG, IVTAM ) 
C  GCONLY(25,14)     -  GCONLY ( IDX, IVTAM ) 
C  INDXCS(5,4)       -  INDXCS ( NT, IS ) 
C  LVTFLG(14)        -  LVTFLG ( IVTAM ) 
C  MYCCEI(7,14)      -  MYCCEI ( IG4, IVTAM ) 
C  MYEGR(4,14)       -  MYEGR  ( IG2, IVTAM ) 
C  MYGTAM(3,14)      -  MYGTAM ( IVTAM ) 
C  MYTGS(20,6,14)    -  MYTGS  ( IG1, IQG, IVTAM ) 
C  MYVTS(3,2,14)     -  MYVTS  ( IG3, IHG, IVTAM ) 
C  ECK(2)            -  ECK    (  NCKIAY ) 
C  NTERMS(4)         -  NTERMS ( IS ) 
C  PCV(2)            -  PCV    (  IAY ) 
C  PBBTR(3,2,25,14)  -  PBBTR  ( IDPB, IAY, IDX, IVTAM ) 
C  TAMBAG(3,3,25,14) -  TAMBAG ( IB, IP, IDX, IVTAM ) 
C  TAMA50(9,14,3,2)  -  TAMA50 ( ID, IVTAM, IGD, IM ) 
C  TAMDR(9,14,3,2)   -  TAMDR  ( ID, IVTAM, IGD, IM ) 
C  TAMZML(9,14,3,2)  -  TAMZML ( ID, IVTAM, IGD, IM ) 
C  TCC(25,14)        -  TCC    ( IDX, IVTAM ) 
C  TCFSTART(69)      -  TCFSTART(ST_INDX) 
C  TDU(25,14)        -  TDU    ( IDX, IVTAM ) 
C  TGS(20,6,14)      -  TGS    ( IG1, IQG, IVTAM ) 
C  THS(2,25,14)      -  THS    ( IFDS, IDX, IVTAM ) 
C  TNK(2)            -  TNK    ( IAY ) 
C  TOB(25,14)        -  TOB    ( IDX, IVTAM ) 
C   
C  Common block dictionary: 
C 
C   Name   Type              Description 
C  ------  ----  ------------------------------------------------------- 
C  /ATPAR1/: 
C 
C  ATPFLG = 2 => user enters EVPAER, AER & all of this block: 
C 
C  LAPSY    I    ATP start year 
C  LAP1ST   I    1st (oldest) model year of the ATP 
C  LAPLST   I    last (youngest) model year of the ATP 
C  LVTFLG   I    which vehicle types covered by the anti-tampering program. 
C                for each of 14 (LDGV, LDGT, HDGV, GAS BUS) flags: 
C                  1 = not covered 
C                  2 = covered 
C                all "1"s = same result as ATPFLG = 1 = no program - the ATP 
C                effectiveness matrices are not applied. 
C                 
C  /ATPAR2/: 
C 
C  ATPFLG = 2 => user enters all of this block: 
C 
C  ATPPGM   I    ATP inspection program type: value must be either 
C                1 (test only) or 2 (test and repair) 
C  ATPFQT   I    ATP inspection program frequency: value must be either 
C                1 (annual) or 2 (biennial) 
C  CRATP    R    compliance rate for anti-tampering program, value must 
C                be from 0 to 100%, inclusive 
C  DISTYP   I    flag indicating which ATP inspection type is conducted: 
C                value must be either 1 (no) or 2 (yes): 
C 
C                DISTYP(1) = Air pump system disablement 
C                DISTYP(2) = Catalyst removal 
C                DISTYP(3) = Fuel inlet restrictor disablement (Requiring 
C                            catalyst replacement) 
C                DISTYP(4) = Tailpipe lead deposit test (Requiring 
C                            catalyst replacement) 
C                DISTYP(5) = EGR disablement 
C                DISTYP(6) = Evaporative system disablement 
C                DISTYP(7) = PCV system disablement 
C                DISTYP(8) = Missing gas cap 
C 
C  /EGSCAL/: 
C 
C  AER      R    anti-tampering program effectiveness rates: ATPFLG = 2 
C                   11 X 11 matrix 1 - previous (up thru LAPSY), HC/CO 
C                   11 X 11 matrix 2 - previous, NOX 
C                   11 X 11 matrix 3 - subsequent (after LAPSY), HC/CO 
C                   11 X 11 matrix 4 - subsequent, NOX 
C                row IC of each 11 X 11 is multiplied times vector CSIZE to 
C                get CSAE(IC). 
C  TGS      R    technology (pollution control equipment) group size 
C  EGS      R    tampering effects group size 
C 
C  /EMECAL/: 
C 
C  EM3W     R    emission impact rates for 3-way catalysts 
C  EMOX     R    emission impact rates for oxidation catalysts 
C  EMEGR    R    emission offsets for EGR failure 
C  EMI      R    emission impact offset assuming 100% tampering rate 
C 
C  /EVPDAT/: 
C 
C  EVPAER   R    evap/PCV ATP effectiveness rates: ATPFLG = 2 => user supplies 
C                one 2 X 2 matrix using 4F5.2 format: 
C                    (1,1) - previous evap    (1,2) - subsequent evap 
C                    (2,1) - previous CC      (2,2) - subsequent CC 
C  EVPTGS   R    evap/PCV technology group sizes 
C  CCEMI    R    PCV ("CC") emission impact rates 
C 
C  Abbreviations - 
C      ATP   -  Anti-Tampering Program 
C      ATER  -  Anti-Tampering Effectiveness Rate 
C      PREV  -  PREVious   to the start of an ATP (IAY case = 1) 
C      SUBS  -  SUBSequent to the start of an ATP (IAY case = 2) 
C      fP&S  -  for PREV and SUBS cases 
C 
C  /ATPAR3/: 
C 
C  AIR      R    ATER for air pump disablement fP&S 
C  CAT      R    ATER for catalyst removal fP&S 
C  ECK      R    ATER for fuel neck/inlet disablement for PREV CAT (NCKIAY = 1) 
C                PREV (NCKIAY = 2), and SUBS (NCKIAY = 3) 
C  TNK      R    ATER for tank/tailpipe lead deposit tests fP&S 
C  EGR      R    ATER for EGR disablement fP&S 
C 
C  /ATPAR4/: 
C 
C  CAN      R    ATER for evap canister disablement fP&S 
C  PCV      R    ATER for PCV disablement fP&S 
C  CAP      R    ATER for missing gas cap fP&S 
C  EFFSAV   R    Saves the default values for EFF 
C  EFF      R    array containing basic ATERs, as follows: 
C 
C        ATPPGM: (1) (2) (1) (2) (1) (2) (1=TEST ONLY 2=TEST AND REPAIR) 
C               PREVIOUS        SUBSEQUENT 
C                  CAT   PREVIOUS 
C               REMOVED  ------- ------- 
C ID=1 
C 
C   ATPFQT=1    --- ---  .80 .40 .80 .80  AIR PUMP INSP. ANNUAL 
C   ATPFQT=2    --- ---  .70 .35 .70 .70                 BIENNIAL 
C 
C ID=2 
C 
C   ----------- --- ---  .95 .48 .95 .95  CATALYST INSP. 
C 
C ID=3 
C 
C   DISTYP(2)=1 .00 .00  .00 .00 .00 .00  FUEL INLET INSP. W/O CAT INSP 
C 
C   DISTYP(2)=2 .59 .30  .33 .17 .70 .70  FUEL INLET INSP. WITH I/M. 
C & IMFLAG   =2 
C & DISTYP(4)=1 
C 
C   DISTYP(2)=2 .33 .17  .17 .08 .30 .30  FUEL INLET INSP. W/O  I/M 
C & IMFLAG   =1 
C & DISTYP(4)=1 
C 
C   DISTYP(2)=2 .75 .38  .75 .38 .85 .85        WITH PLUMBTESMO CHECK 
C & DISTYP(4)=2 
C 
C ID=4 
C 
C   DISTYP(2)=1 --- ---  .00 .00 .00 .00  PLUMBTESMO CHECK WITH CAT INSP 
C   DISTYP(2)=2 --- ---  .75 .38 .85 .85                   W/O  CAT INSP 
C 
C ID=5 
C 
C   ----------- --- ---  .00 .00 .00 .00  EGR INSP 
C 
C ID=6 
C 
C   ATPFQT=1    --- ---  .70 .35 .70 .35  EVAP SYS INSP  ANNUAL 
C   ATPFQT=2    --- ---  .57 .29 .57 .29                 BIENNIAL 
C 
C ID=7 
C 
C   ATPFQT=1    --- ---  .70 .35 .70 .35  PCV  SYS INSP  ANNUAL 
C   ATPFQT=2    --- ---  .56 .28 .56 .28                 BIENNIAL 
C 
C ID=8 
C 
C   ----------- --- ---  .50 .50 .50 .50  GAS CAP INSP 
C 
C 
C  /LOOKUP/: 
C 
C  IVTAM, IQG, IPG, JPGD, IHG & IGCSF are defined in the subscript dictionary. 
C 
C  /MYCODE/: 
C 
C  MY       I    model year being evaluated on current IDX loop pass 
C  LDXSY    I    pointer: anti-tampering program start year's location in the 
C                calendar year's 25 model year window; JDX ordering 
C  LMYRVT   I    flag: 1 = IDX not in ATP my range &/or vt not covered by ATP 
C                      2 = IDX is in ATP my range & vt is covered by ATP 
C  IMDXSY   I    pointer: inspection / maintenance program start year's 
C                location in the calendar year's 20 model year window; 
C                JDX ordering 
C  IMKINK   I    flag: 1 = tampering rate curves straight; 
C                      2 = tampering rate curves have I/M kink 
C  IDX, JDX and IAY are defined in the subscript dictionary. 
C 
C 
C  /MYUB1/: 
C 
C  MYTGS    I    model year range upper bound years for ITAMPT indexing TGS 
C  MYEGR    I    model year range upper bound years for ITAMPT indexing EMEGR 
C  MYVTS    I    model year range upper bound years for ITAMPT indexing EVPTGS 
C 
C  /MYUB2/: 
C 
C  MYCCEI   I    model year range upper bound years for ITAMPT indexing CCEMI 
C 
C  /OMTCOM/: 
C 
C  TCFSTART R    Start emission temperature correction factor, calculated in BIGTCF. 
C  TCFRUN   R    Running emission temperature correction factor, calculated in BIGTCF. 
C  TCF12    R    Start emission temperature correction factor, 12-hour soak, 
C                calculated in BIGTCF. 
C 
C  /SIZCAL/: 
C 
C  ATR      R    exhaust ef BTRs adjusted for equipment type & size 
C  TGSUSE   R    technology group size from TGS used on EFFGRP IQG loop pass. 
C  CSIZE    R    effects category size 
C  CSAE     R    CSIZE modified by AER 
C 
C  /SIZPAR/: 
C 
C  NTERMS   I    number of overlap category sizes to be subtracted from 
C                ATR(IS,IAY) to obtain CSIZE(IS+7,IAY). 
C  INDXCS   I    pointer to the category sizes to be subtracted from 
C                ATR(IS,IAY) to obtain CSIZE(IS+7,IAY) 
C 
C  /TAMEQ1/: 
C 
C  TAMZML   R    basic tampering rate equation's zero mile level (intercept) 
C  TAMDR    R    basic tampering rate equation's deterioration rate (slope) 
C  TAMA50   R    adjustment to TAMDR for mileage over 50,000 miles 
C 
C  TAMFLG = 2 => user replaces all of TAMZML & TAMDR with values read in by 
C                ONESEC (1-time Data Section). 
C 
C  /TAMEQ2/: 
C 
C  MYGTAM   I    ub of myg of each tampering rate case for each vehicle class 
C  F50K     R    50,000 miles vmt, scaled by division by 10**4 (= 5.0) 
C                in BEF = the LDGV 81+ HC/CO curves deflection point 
C                in DISCAL = the IVTAM IGD IP curves deflection point 
C 
C  /TAMEQ3/: 
C 
C  CUMCAP   R    cumulative mileage cap (/10**4) for tampering rate calculation 
C                 (i.e., at 13.0, tampering rates all top off at 130,000 miles) 
C  CUMALL    R    mileage accumulated segments for straight line subcases 
C                (no kink from I/M startup): 
C                  (1) MY through ICY for DISATP cases 1, 2 & 3: 
C                      all "previous", all "subsequent", overall p+s mileage 
C                  (2) MY through LAPSY for DISATP case 3: "previous" mileage 
C 
C  Kinked rate lines mileage accumulated sections: 
C 
C  CUMNIM   R    from MY through ICYIM(1/2) 
C  CUMIM    R    after ICYIM(1/2) through last year in segment (ICY or LAPSY) 
C  GT50KA   R    mileage on the CUMALL segment in excess of 50K & under CUMCAP 
C  GT50KN   R    mileage in excess of 50,000 on the CUMNIM segment 
C  GT50KI   R    mileage in excess of 50,000 on the CUMNIM+CUMIM segments, 
C                except if CUMNIM GE 50,000, GT50KI = CUMIM 
C 
C  /TAMEQ4/: 
C 
C  BTR      R    base tampering rate for standalone disablement cases 
C 
C  /TAMOUT/: 
C 
C  TAMBAG   R    tampering effects on the exhaust bag data.  Calculated in 
C                subroutine TAMPER.  Passed to BIGCFX to get OMTTAM. 
C  TAMBAS   R    Base tampering effects. 
C  TEX      R    ATP tampering effects.  Difference between TAMBAG and TAMBAS. 
C  THS      R    percentage of given IV & MY subject to tampered HS loss 
C  TDU      R    percentage of given IV & MY subject to tampered DU loss 
C  TCC      R    tampered CC offset = affected fleet %age * tampered CC loss 
C 
C  /TAMOU2/: 
C 
C  TOB     R    Probability (Evap Cannister is tampered | IVTAM is EC equipped) 
C               for each IDX x IVTAM case.  Unlike THS & TDU, TOB is not wgted 
C               by PCTEQP = the Probability (IVTAM is EC equipped), because it 
C               is used to figure the tampering impact on OnBoard Evaporative 
C               Cannisters, where it is assumed that if OB is in effect in a MY 
C               then 100% of the vehicle types covered are OB EC equipped. 
C  TDU      R    percentage of given IV & MY subject to tampered gas caps 
C 
C  /TAMPB1/: 
C 
C  LPOD     I    Lead Phase Out Date: last year of leaded gasoline, currently 
C                hard coded at 1995 
C  PBBTR    R    misfueling tampering categories rates, computed as if ICY=LPOD 
C 
C  Local array subscripts: 
C 
C  MYTGS1(20,6,3)  -  MYTGS1 ( IG1, IQG, IVTAM ) 
C  MYTGS2(20,6,3)  -  MYTGS2 ( IG1, IQG, IVTAM ) 
C  MYTGS3(20,6,3)  -  MYTGS3 ( IG1, IQG, IVTAM ) 
C  MYTGS4(20,6,3)  -  MYTGS4 ( IG1, IQG, IVTAM ) 
C  MYTGS5(20,6,2)  -  MYTGS5 ( IG1, IQG, IVTAM ) 
C  TA5IM(9,4,3)    -  TA5IM  ( ID, IVTAM, IGD ) 
C  TA5NON(9,4,3)   -  TA5NON ( ID, IVTAM, IGD ) 
C  TDRIM(9,4,3)    -  TDRIM  ( ID, IVTAM, IGD ) 
C  TDRNON(9,4,3)   -  TDRNON ( ID, IVTAM, IGD ) 
C  TGS01(20,6,3)   -  TGS01  ( IG1, IQG, IVTAM) 
C  TGS02(20,6,3)   -  TGS02  ( IG1, IQG, IVTAM ) 
C  TGS03(20,6,3)   -  TGS03  ( IG1, IQG, IVTAM ) 
C  TGS04(20,6,3)   -  TGS04  ( IG1, IQG, IVTAM ) 
C  TGS05(20,6,2)   -  TGS05  ( IG1, IQG, IVTAM ) 
C  TZLIM(9,4,3)    -  TZLIM  ( ID, IVTAM, IGD ) 
C  TZLNON(9,4,3)   -  TZLNON ( ID, IVTAM, IGD ) 
C 
C  Local array dictionary: 
C 
C   Name   Type              Description 
C  ------  ----  ------------------------------------------------------- 
C  MYTGS1   I    my range ub years for LDGV/T1A/T1B for ITAMPT indexing TGS 
C  MYTGS2   I    my range ub years for LDGT34,HDGV2B ITAMPT indexing TGS 
C  MYTGS3   I    my range ub years for 3/4/HDGV5 for ITAMPT indexing TGS 
C  MYTGS4   I    my range ub years for 6/7/HDGV8A for ITAMPT indexing TGS 
C  MYTGS5   I    my range ub years for HDGV8B, G. BUS for ITAMPT indexing TGS 
C  TA5IMx   R    TAMA50s for the I/M case 
C  TA5NIx   R    TAMA50s for the Non-I/M case 
C  TDRIMx   R    TAMDRs for the I/M case 
C  TDRNIx   R    TAMDRs for the Non-I/M case 
C  TGSHDV   R    technology group sizes for HDGV 
C  TGSLDV   R    technology group sizes for LDGV/T1/T2 
C  TZLIMx   R    TAMZMLs for the I/M case 
C  TZLNIx   R    TAMZMLs for the Non-I/M case 
C 
C  Notes: 
C 
C  The following common blocks defined and used elsewhere in MOBILE4 contain 
C  data used by the tampering effects section: 
C 
C  /CUMCOM/: CUMMIL - cumulative mileage array (JDX,IV) 
C  /FLAGS2/: IMFLAG  - I/M program flag 
C  /IMPAR1/: ICYIM  - I/M program(s) start year and end years 
C            MODYR -  first and last model year covered by I/M program(s) 
C            ILDT   - vehicle types covered by I/M program(s) 
C 
C  See the block data sections for these commons for more information on their 
C  variables / arrays. 
C 
C  FORTRAN 90 limits the number of continuation lines to 39.  Use equivalenced 
C  local arrays to initialize common block arrays if this limit will be 
C  exceeded otherwise.  This is true throughout all the block datas. 
C 
C  MOBILE4 classifies tampering by fuel delivery system (TBI vs PFI vs CAR). 
C  MOBILE4 has 8 disablement categories (ID = 8), adding a second evaporative 
C  case (gas cap + evap cannisters) as case #8. 
C 
C  MOBILE4 replaced EVPEMI, MYVEI, IH and TAMEVP with CCEMI, MYCCEI, THS, 
C  TDU and TCC, in support of the Evaporative Study II change requirements. 
C 
C  OPENLP deleted by removal of MOBILE4 oxy. fuels support in MOBILE 4.1. 
C 
C 
C  BD 25 was modified for MOBILE4.1 v4 & v5 to expand the tampering rate MYGs 
C  from 2 to 3, change the over 50K cases to adjustments to the under 50K rates, 
C  and update all tampering rates. 
C 
C  23 Feb 96 DynTel@MLA 2-611 Added array EFFSAV to save default values 
C    for EFF. 
C 
C 27 Nov 97 DynTel@MLA 2-622  Expanded the arrays to accomodate the new  
C           Heavy Duty Gasoline, LDGT, and GAS BUS vehicle categories. 
C 
C 13 Apr 98 DynTel@ID 2-643 Updated TGS and MYTGS data using Technology 
C spreadsheet from Mobile 4.1 for 1981-1989 years for Air Pump (no EGR data 
C was available) and data supplied by Pechan/EEA Projections for 1990-1995 years. 
C The data for 1980 and earlier remained without changes.  
C 
C 23 Jul 98 @DynTel-MLA 2-000  Changed LDGT names to conform to 
C           certification naming conventions. 
C 
      IMPLICIT NONE 
      INCLUDE 'ATPAR1.I' 
      INCLUDE 'ATPAR2.I' 
      INCLUDE 'ATPAR3.I' 
      INCLUDE 'ATPAR4.I' 
      INCLUDE 'EGSCAL.I' 
      INCLUDE 'EMECAL.I' 
      INCLUDE 'EVPDAT.I' 
      INCLUDE 'LOOKUP.I' 
      INCLUDE 'MYCODE.I' 
      INCLUDE 'MYUB1.I' 
      INCLUDE 'MYUB2.I' 
      INCLUDE 'OMTCOM.I' 
      INCLUDE 'PARAME.I' 
      INCLUDE 'SIZCAL.I' 
      INCLUDE 'SIZPAR.I' 
      INCLUDE 'TAMEQ1.I' 
      INCLUDE 'TAMEQ2.I' 
      INCLUDE 'TAMEQ3.I' 
      INCLUDE 'TAMEQ4.I' 
      INCLUDE 'TAMOUT.I' 
      INCLUDE 'TAMOU2.I' 
      INCLUDE 'TAMPB1.I' 
C 
      INTEGER MYTGS1(20,6,3),MYTGS2(20,6,3),MYTGS3(20,6,3), 
     *          MYTGS4(20,6,3),MYTGS5(20,6,2) 
      REAL TGS01(20,6,3),TGS02(20,6,3),TGS03(20,6,3), 
     *          TGS04(20,6,3),TGS05(20,6,2) 
      REAL TZLNI1(9,14),TZLNI2(9,14),TZLNI3(9,14) 
      REAL TZLIM1(9,14),TZLIM2(9,14),TZLIM3(9,14) 
      REAL TDRNI1(9,14),TDRNI2(9,14),TDRNI3(9,14) 
      REAL TDRIM1(9,14),TDRIM2(9,14),TDRIM3(9,14) 
      REAL TA5NI1(9,14),TA5NI2(9,14),TA5NI3(9,14) 
      REAL TA5IM1(9,14),TA5IM2(9,14),TA5IM3(9,14) 
C 
      EQUIVALENCE (MYTGS(1,1,1), MYTGS1(1,1,1)) 
      EQUIVALENCE (MYTGS(1,1,4), MYTGS2(1,1,1)) 
      EQUIVALENCE (MYTGS(1,1,7), MYTGS3(1,1,1)) 
      EQUIVALENCE (MYTGS(1,1,10),MYTGS4(1,1,1)) 
      EQUIVALENCE (MYTGS(1,1,13),MYTGS5(1,1,1)) 
      EQUIVALENCE (TAMA50(1,1,1,1),TA5NI1(1,1)) 
      EQUIVALENCE (TAMA50(1,1,2,1),TA5NI2(1,1)) 
      EQUIVALENCE (TAMA50(1,1,3,1),TA5NI3(1,1)) 
      EQUIVALENCE (TAMA50(1,1,1,2),TA5IM1(1,1)) 
      EQUIVALENCE (TAMA50(1,1,2,2),TA5IM2(1,1)) 
      EQUIVALENCE (TAMA50(1,1,3,2),TA5IM3(1,1)) 
      EQUIVALENCE (TAMDR(1,1,1,1),TDRNI1(1,1)) 
      EQUIVALENCE (TAMDR(1,1,2,1),TDRNI2(1,1)) 
      EQUIVALENCE (TAMDR(1,1,3,1),TDRNI3(1,1)) 
      EQUIVALENCE (TAMDR(1,1,1,2),TDRIM1(1,1)) 
      EQUIVALENCE (TAMDR(1,1,2,2),TDRIM2(1,1)) 
      EQUIVALENCE (TAMDR(1,1,3,2),TDRIM3(1,1)) 
      EQUIVALENCE (TAMZML(1,1,1,1),TZLNI1(1,1)) 
      EQUIVALENCE (TAMZML(1,1,2,1),TZLNI2(1,1)) 
      EQUIVALENCE (TAMZML(1,1,3,1),TZLNI3(1,1)) 
      EQUIVALENCE (TAMZML(1,1,1,2),TZLIM1(1,1)) 
      EQUIVALENCE (TAMZML(1,1,2,2),TZLIM2(1,1)) 
      EQUIVALENCE (TAMZML(1,1,3,2),TZLIM3(1,1)) 
      EQUIVALENCE (TGS(1,1,1), TGS01(1,1,1)) 
      EQUIVALENCE (TGS(1,1,4), TGS02(1,1,1)) 
      EQUIVALENCE (TGS(1,1,7), TGS03(1,1,1)) 
      EQUIVALENCE (TGS(1,1,10),TGS04(1,1,1)) 
      EQUIVALENCE (TGS(1,1,13),TGS05(1,1,1)) 
C 
C  /ATPAR1/: values can be replaced by user entered values. 
C 
      DATA LAPSY,LAP1ST,LAPLST/3*NEYR/,LVTFLG/14*1/ 
C 
C  /ATPAR2/: values can be replaced by user entered values. 
C 
      DATA ATPPGM,ATPFQT,DISTYP/10*1/,CRATP/100./ 
C 
C  /ATPAR3/: 
C 
      DATA AIR,CAT,ECK,TNK,EGR/11*0./ 
C 
C  /ATPAR4/: 
C 
      DATA CAN,PCV,CAP/6*0./ 
C 
C  Array EFF: continuation line code corresponds to ID: 
C 
C          IATP=1    IATP=1    IATP=1   IATP=2    IATP=2    IATP=2 
C          ICOL=1    ICOL=2    ICOL=3   ICOL=1    ICOL=2    ICOL=3 
      DATA EFF/ 
     1    .00,.00,  .80,.40,  .80,.80, .00,.00,  .80,.40,  .80,.80, 
     1    .00,.00,  .70,.35,  .70,.70, .00,.00,  .70,.35,  .70,.70, 
     2    .00,.00,  .95,.48,  .95,.95, .00,.00,  .95,.48,  .95,.95, 
     3    .00,.00,  .00,.00,  .00,.00, .00,.00,  .00,.00,  .00,.00, 
     3    .59,.30,  .33,.17,  .70,.70, .59,.30,  .33,.17,  .70,.70, 
     3    .33,.17,  .17,.08,  .30,.30, .33,.17,  .17,.08,  .30,.30, 
     3    .75,.38,  .75,.38,  .85,.85, .75,.38,  .75,.38,  .85,.85, 
     4    .00,.00,  .00,.00,  .00,.00, .00,.00,  .00,.00,  .00,.00, 
     4    .00,.00,  .75,.38,  .85,.85, .00,.00,  .75,.38,  .85,.85, 
     5    .00,.00,  .00,.00,  .00,.00, .00,.00,  .00,.00,  .00,.00, 
     6    .00,.00,  .70,.35,  .70,.35, .00,.00,  .70,.35,  .70,.35, 
     6    .00,.00,  .57,.29,  .57,.29, .00,.00,  .57,.29,  .57,.29, 
     7    .00,.00,  .70,.35,  .70,.35, .00,.00,  .70,.35,  .70,.35, 
     7    .00,.00,  .56,.28,  .56,.28, .00,.00,  .56,.28,  .56,.28, 
     8    .00,.00,  .50,.25,  .50,.25, .00,.00,  .50,.25,  .50,.25/ 
C 
C  Array EFFSAV. This duplicates array EFF. These values never change. 
C  (Some values in EFF may be replaced by user inputs.) 
C 
      DATA EFFSAV/ 
     1    .00,.00,  .80,.40,  .80,.80, .00,.00,  .80,.40,  .80,.80, 
     1    .00,.00,  .70,.35,  .70,.70, .00,.00,  .70,.35,  .70,.70, 
     2    .00,.00,  .95,.48,  .95,.95, .00,.00,  .95,.48,  .95,.95, 
     3    .00,.00,  .00,.00,  .00,.00, .00,.00,  .00,.00,  .00,.00, 
     3    .59,.30,  .33,.17,  .70,.70, .59,.30,  .33,.17,  .70,.70, 
     3    .33,.17,  .17,.08,  .30,.30, .33,.17,  .17,.08,  .30,.30, 
     3    .75,.38,  .75,.38,  .85,.85, .75,.38,  .75,.38,  .85,.85, 
     4    .00,.00,  .00,.00,  .00,.00, .00,.00,  .00,.00,  .00,.00, 
     4    .00,.00,  .75,.38,  .85,.85, .00,.00,  .75,.38,  .85,.85, 
     5    .00,.00,  .00,.00,  .00,.00, .00,.00,  .00,.00,  .00,.00, 
     6    .00,.00,  .70,.35,  .70,.35, .00,.00,  .70,.35,  .70,.35, 
     6    .00,.00,  .57,.29,  .57,.29, .00,.00,  .57,.29,  .57,.29, 
     7    .00,.00,  .70,.35,  .70,.35, .00,.00,  .70,.35,  .70,.35, 
     7    .00,.00,  .56,.28,  .56,.28, .00,.00,  .56,.28,  .56,.28, 
     8    .00,.00,  .50,.25,  .50,.25, .00,.00,  .50,.25,  .50,.25/ 
C 
C  /EGSCAL/: TGS never changes, EGS is always computed 
C            during execution. 
C 
C  AER is initialized to 4 11 X 11 matrices with 1.'s down the major diagonal 
C  and 0.'s for all other cells.  The continuation line code is: 
C    1 = previous (before ATP start year), HC/CO 
C    2 = previous, NOX 
C    3 = subsequent (at or after ATP start year), HC/CO 
C    4 = subsequent, NOX 
C 
      DATA AER/ 
     1  1.,11*0., 1.,11*0., 1.,11*0., 1.,11*0., 1.,11*0., 1.,11*0., 
     1  1.,11*0., 1.,11*0., 1.,11*0., 1.,11*0., 1., 
     2  1.,11*0., 1.,11*0., 1.,11*0., 1.,11*0., 1.,11*0., 1.,11*0., 
     2  1.,11*0., 1.,11*0., 1.,11*0., 1.,11*0., 1., 
     3  1.,11*0., 1.,11*0., 1.,11*0., 1.,11*0., 1.,11*0., 1.,11*0., 
     3  1.,11*0., 1.,11*0., 1.,11*0., 1.,11*0., 1., 
     4  1.,11*0., 1.,11*0., 1.,11*0., 1.,11*0., 1.,11*0., 1.,11*0., 
     4  1.,11*0., 1.,11*0., 1.,11*0., 1.,11*0., 1./ 
C 
C  TGS has up to 10 percentages for each vehicle type / technology 
C  combination.  The remaining cells are dummied out with zeroes. 
C  For each IVTAM, the continuation line technology code is: 
C    1 = air pump only 
C    2 = air pump / catalyst 
C    3 = no air pump / catalyst  
C    4 = EGR only 
C    5 = EGR / 3-way catalyst 
C    6 = no EGR / 3-way catalyst 
C 
      DATA TGS01/ 
C  LDGV 
     1  .0000,.0500,.1000,.3000,.1500,.1000,.0500,.0000,.0000,.0000, 
     1  .0000,.0000,.0000,.0000,.0000,.0000,.0000,.0000,.0000,.0000, 
     2  .0000,.3000,.2000,.2500,.6500,.7000,.7800,.6500,.5200,.4400, 
     2  .3600,.1900,.1400,.2720,.2272,.1786,.0920,.0752,.0586,.0000,  
     3  .0000,.5000,.5500,.6500,.3000,.2200,.3500,.4800,.5600,.6400,  
     3  .8100,.8600,.7280,.7728,.8214,.9080,.9248,.9414,.0000,.0000, 
     4  .0000,.8000,.9000,.0500,.0000,.0041,.0021,.0029,.0000,.0000, 
     4  .0000,.0000,.0000,.0000,.0000,.0000,.0000,.0000,.0000,.0000, 
     5  .0000,.0700,.8500,.9300,.9000,.7538,.7494,.6344,.7825,.8152, 
     5  .9082,.0000,.0000,.0000,.0000,.0000,.0000,.0000,.0000,.0000, 
     6  .0000,.0100,.0300,.0700,.1000,.2421,.2485,.3627,.2175,.1848,  
     6  .0918,.0000,.0000,.0000,.0000,.0000,.0000,.0000,.0000,.0000, 
C  LDGT1 
     1  .0000,.0500,.1000,.3000,.1000,.0300,.0000,.0300,.0000,.0000, 
     1  .0000,.0000,.0000,.0000,.0000,.0000,.0000,.0000,.0000,.0000, 
     2  .0000,.3000,.2000,.4000,.9600,.9300,.8000,.6200,.4800,.2700,  
     2  .3200,.5317,.3360,.2461,.2394,.1841,.1587,.0000,.0000,.0000, 
     3  .0000,.4000,.5000,.5500,.4000,.0100,.0400,.0700,.2000,.3800, 
     3  .5200,.7300,.6800,.4683,.6640,.7539,.7606,.8159,.8413,.0000, 
     4  .0000,.8000,.9000,1.000,.9500,.9000,.8000,.7000,.6000,.5000,  
     4  .0500,.0314,.0105,.0135,.0000,.0000,.0000,.0000,.0000,.0000, 
     5  .0000,.0500,.1000,.2000,.3000,.4000,.5000,.9500,.6947,.7022,  
     5  .6437,.7630,.9167,.8904,.0000,.0000,.0000,.0000,.0000,.0000, 
     6  .0000,.2739,.2873,.3428,.2370,.0833,.1096,.0000,.0000,.0000, 
     6  .0000,.0000,.0000,.0000,.0000,.0000,.0000,.0000,.0000,.0000,   
C  LDGT2 
     1  .0000,.0500,.1000,.3000,.1000,.0300,.0000,.0300,.0000,.0000, 
     1  .0000,.0000,.0000,.0000,.0000,.0000,.0000,.0000,.0000,.0000, 
     2  .0000,.3000,.2000,.4000,.9600,.9300,.8000,.6200,.4800,.2700,  
     2  .3200,.5317,.3360,.2461,.2394,.1841,.1587,.0000,.0000,.0000, 
     3  .0000,.4000,.5000,.5500,.4000,.0100,.0400,.0700,.2000,.3800, 
     3  .5200,.7300,.6800,.4683,.6640,.7539,.7606,.8159,.8413,.0000, 
     4  .0000,.8000,.9000,1.000,.9500,.9000,.8000,.7000,.6000,.5000,  
     4  .0500,.0314,.0105,.0135,.0000,.0000,.0000,.0000,.0000,.0000, 
     5  .0000,.0500,.1000,.2000,.3000,.4000,.5000,.9500,.6947,.7022,  
     5  .6437,.7630,.9167,.8904,.0000,.0000,.0000,.0000,.0000,.0000, 
     6  .0000,.2739,.2873,.3428,.2370,.0833,.1096,.0000,.0000,.0000, 
     6  .0000,.0000,.0000,.0000,.0000,.0000,.0000,.0000,.0000,.0000/ 
C   
      DATA TGS02/ 
C  LDGT3 
     1  .0000,.0300,.0000,.0300,.0000,.0000,.0000,.0000,.0000,.0000, 
     1  .0000,.0000,.0000,.0000,.0000,.0000,.0000,.0000,.0000,.0000, 
     2  .0000,.5000,.9600,.9300,.8000,.6200,.4800,.2700,.3200,.5317,  
     2  .3360,.2461,.2394,.1841,.1587,.0000,.0000,.0000,.0000,.0000, 
     3  .0000,.5000,.0100,.0400,.0700,.2000,.3800,.5200,.7300,.6800,  
     3  .4683,.6640,.7539,.7606,.8159,.8413,.0000,.0000,.0000,.0000, 
     4  .0000,.3000,1.000,.9500,.9000,.8000,.7000,.6000,.5000,.0500,  
     4  .0314,.0105,.0135,.0000,.0000,.0000,.0000,.0000,.0000,.0000, 
     5  .0000,.0500,.1000,.2000,.3000,.4000,.5000,.9500,.6947,.7022,  
     5  .6437,.7630,.9167,.8904,.0000,.0000,.0000,.0000,.0000,.0000, 
     6  .0000,.2739,.2873,.3428,.2370,.0833,.1096,.0000,.0000,.0000,  
     6  .0000,.0000,.0000,.0000,.0000,.0000,.0000,.0000,.0000,.0000, 
C  LDGT4 
     1  .0000,.0300,.0000,.0300,.0000,.0000,.0000,.0000,.0000,.0000, 
     1  .0000,.0000,.0000,.0000,.0000,.0000,.0000,.0000,.0000,.0000, 
     2  .0000,.5000,.9600,.9300,.8000,.6200,.4800,.2700,.3200,.5317,  
     2  .3360,.2461,.2394,.1841,.1587,.0000,.0000,.0000,.0000,.0000, 
     3  .0000,.5000,.0100,.0400,.0700,.2000,.3800,.5200,.7300,.6800,  
     3  .4683,.6640,.7539,.7606,.8159,.8413,.0000,.0000,.0000,.0000, 
     4  .0000,.3000,1.000,.9500,.9000,.8000,.7000,.6000,.5000,.0500,  
     4  .0314,.0105,.0135,.0000,.0000,.0000,.0000,.0000,.0000,.0000, 
     5  .0000,.0500,.1000,.2000,.3000,.4000,.5000,.9500,.6947,.7022,  
     5  .6437,.7630,.9167,.8904,.0000,.0000,.0000,.0000,.0000,.0000, 
     6  .0000,.2739,.2873,.3428,.2370,.0833,.1096,.0000,.0000,.0000,  
     6  .0000,.0000,.0000,.0000,.0000,.0000,.0000,.0000,.0000,.0000, 
C  HDGV2B 
     1  .0000,.0700,.0000,.0000,.0000,.0000,.0000,.0000,.0000,.0000, 
     1  .0000,.0000,.0000,.0000,.0000,.0000,.0000,.0000,.0000,.0000, 
     2  .0000,.3000,.5000,.5500,.6500,1.000,1.000,.0000,.0000,.0000, 
     2  .0000,.0000,.0000,.0000,.0000,.0000,.0000,.0000,.0000,.0000, 
     3  .0000,.0500,.0000,.1000,.5000,.4500,.3500,.0000,.0000,.0000, 
     3  .0000,.0000,.0000,.0000,.0000,.0000,.0000,.0000,.0000,.0000, 
     4  .0000,.9500,1.000,.8500,1.000,.7300,.7000,.0000,.0000,.0000, 
     4  .0000,.0000,.0000,.0000,.0000,.0000,.0000,.0000,.0000,.0000, 
     5  .0000,.1500,.0000,.2700,.0000,.7500,.7500,.0000,.0000,.0000, 
     5  .0000,.0000,.0000,.0000,.0000,.0000,.0000,.0000,.0000,.0000, 
     6  .0000,.0000,.0000,.3000,.2500,.2500,.0000,.0000,.0000,.0000, 
     6  .0000,.0000,.0000,.0000,.0000,.0000,.0000,.0000,.0000,.0000/ 
      DATA TGS03/ 
C  HDGV3 
     1  .0000,.0700,.0000,.0000,.0000,.0000,.0000,.0000,.0000,.0000, 
     1  .0000,.0000,.0000,.0000,.0000,.0000,.0000,.0000,.0000,.0000, 
     2  .0000,.3000,.5000,.5500,.6500,1.000,1.000,.0000,.0000,.0000, 
     2  .0000,.0000,.0000,.0000,.0000,.0000,.0000,.0000,.0000,.0000, 
     3  .0000,.0500,.0000,.1000,.5000,.4500,.3500,.0000,.0000,.0000, 
     3  .0000,.0000,.0000,.0000,.0000,.0000,.0000,.0000,.0000,.0000, 
     4  .0000,.9500,1.000,.8500,1.000,.7300,.7000,.0000,.0000,.0000, 
     4  .0000,.0000,.0000,.0000,.0000,.0000,.0000,.0000,.0000,.0000, 
     5  .0000,.1500,.0000,.2700,.0000,.7500,.7500,.0000,.0000,.0000, 
     5  .0000,.0000,.0000,.0000,.0000,.0000,.0000,.0000,.0000,.0000, 
     6  .0000,.0000,.0000,.3000,.2500,.2500,.0000,.0000,.0000,.0000, 
     6  .0000,.0000,.0000,.0000,.0000,.0000,.0000,.0000,.0000,.0000, 
C  HDGV4 
     1  .0000,.0700,.7000,.6000,.2000,.0000,.0000,.0000,.0000,.0000, 
     1  .0000,.0000,.0000,.0000,.0000,.0000,.0000,.0000,.0000,.0000, 
     2  .0000,.3000,.3000,.4000,.8000,1.000,1.000,.0000,.0000,.0000, 
     2  .0000,.0000,.0000,.0000,.0000,.0000,.0000,.0000,.0000,.0000, 
     3  .0000,.0500,.0000,.1000,.0000,.0000,.0000,.0000,.0000,.0000, 
     3  .0000,.0000,.0000,.0000,.0000,.0000,.0000,.0000,.0000,.0000, 
     4  .0000,.9500,1.000,.8500,1.000,1.000,1.000,.0000,.0000,.0000, 
     4  .0000,.0000,.0000,.0000,.0000,.0000,.0000,.0000,.0000,.0000, 
     5  .0000,.1500,.0000,.0000,.0000,1.000,1.000,.0000,.0000,.0000, 
     5  .0000,.0000,.0000,.0000,.0000,.0000,.0000,.0000,.0000,.0000, 
     6  .0000,.0000,.0000,.0000,.0000,.0000,.0000,.0000,.0000,.0000, 
     6  .0000,.0000,.0000,.0000,.0000,.0000,.0000,.0000,.0000,.0000, 
C  HDGV5 
     1  .0000,.0700,.7000,.6000,.2000,.0000,.0000,.0000,.0000,.0000, 
     1  .0000,.0000,.0000,.0000,.0000,.0000,.0000,.0000,.0000,.0000, 
     2  .0000,.3000,.3000,.4000,.8000,1.000,1.000,.0000,.0000,.0000, 
     2  .0000,.0000,.0000,.0000,.0000,.0000,.0000,.0000,.0000,.0000, 
     3  .0000,.0500,.0000,.1000,.0000,.0000,.0000,.0000,.0000,.0000, 
     3  .0000,.0000,.0000,.0000,.0000,.0000,.0000,.0000,.0000,.0000, 
     4  .0000,.9500,1.000,.8500,1.000,1.000,1.000,.0000,.0000,.0000, 
     4  .0000,.0000,.0000,.0000,.0000,.0000,.0000,.0000,.0000,.0000, 
     5  .0000,.1500,.0000,.0000,.0000,1.000,1.000,.0000,.0000,.0000, 
     5  .0000,.0000,.0000,.0000,.0000,.0000,.0000,.0000,.0000,.0000, 
     6  .0000,.0000,.0000,.0000,.0000,.0000,.0000,.0000,.0000,.0000, 
     6  .0000,.0000,.0000,.0000,.0000,.0000,.0000,.0000,.0000,.0000/ 
      DATA TGS04/ 
C  HDGV6 
     1  .0000,.0700,.7000,.6000,.2000,.0000,.0000,.0000,.0000,.0000, 
     1  .0000,.0000,.0000,.0000,.0000,.0000,.0000,.0000,.0000,.0000, 
     2  .0000,.3000,.3000,.4000,.8000,1.000,1.000,.0000,.0000,.0000, 
     2  .0000,.0000,.0000,.0000,.0000,.0000,.0000,.0000,.0000,.0000, 
     3  .0000,.0500,.0000,.1000,.0000,.0000,.0000,.0000,.0000,.0000, 
     3  .0000,.0000,.0000,.0000,.0000,.0000,.0000,.0000,.0000,.0000, 
     4  .0000,.9500,1.000,.8500,1.000,1.000,1.000,.0000,.0000,.0000, 
     4  .0000,.0000,.0000,.0000,.0000,.0000,.0000,.0000,.0000,.0000, 
     5  .0000,.1500,.0000,.0000,.0000,1.000,1.000,.0000,.0000,.0000, 
     5  .0000,.0000,.0000,.0000,.0000,.0000,.0000,.0000,.0000,.0000, 
     6  .0000,.0000,.0000,.0000,.0000,.0000,.0000,.0000,.0000,.0000, 
     6  .0000,.0000,.0000,.0000,.0000,.0000,.0000,.0000,.0000,.0000, 
C  HDGV7 
     1  .0000,.0700,.7000,.6000,.2000,.0000,.0000,.0000,.0000,.0000, 
     1  .0000,.0000,.0000,.0000,.0000,.0000,.0000,.0000,.0000,.0000, 
     2  .0000,.3000,.3000,.4000,.8000,1.000,1.000,.0000,.0000,.0000, 
     2  .0000,.0000,.0000,.0000,.0000,.0000,.0000,.0000,.0000,.0000, 
     3  .0000,.0500,.0000,.1000,.0000,.0000,.0000,.0000,.0000,.0000, 
     3  .0000,.0000,.0000,.0000,.0000,.0000,.0000,.0000,.0000,.0000, 
     4  .0000,.9500,1.000,.8500,1.000,1.000,1.000,.0000,.0000,.0000, 
     4  .0000,.0000,.0000,.0000,.0000,.0000,.0000,.0000,.0000,.0000, 
     5  .0000,.1500,.0000,.0000,.0000,1.000,1.000,.0000,.0000,.0000, 
     5  .0000,.0000,.0000,.0000,.0000,.0000,.0000,.0000,.0000,.0000, 
     6  .0000,.0000,.0000,.0000,.0000,.0000,.0000,.0000,.0000,.0000, 
     6  .0000,.0000,.0000,.0000,.0000,.0000,.0000,.0000,.0000,.0000, 
C  HDGV8A 
     1  .0000,.0700,.7000,.6000,.2000,.0000,.0000,.0000,.0000,.0000, 
     1  .0000,.0000,.0000,.0000,.0000,.0000,.0000,.0000,.0000,.0000, 
     2  .0000,.3000,.3000,.4000,.8000,1.000,1.000,.0000,.0000,.0000, 
     2  .0000,.0000,.0000,.0000,.0000,.0000,.0000,.0000,.0000,.0000, 
     3  .0000,.0500,.0000,.1000,.0000,.0000,.0000,.0000,.0000,.0000, 
     3  .0000,.0000,.0000,.0000,.0000,.0000,.0000,.0000,.0000,.0000, 
     4  .0000,.9500,1.000,.8500,1.000,1.000,1.000,.0000,.0000,.0000, 
     4  .0000,.0000,.0000,.0000,.0000,.0000,.0000,.0000,.0000,.0000, 
     5  .0000,.1500,.0000,.0000,.0000,1.000,1.000,.0000,.0000,.0000, 
     5  .0000,.0000,.0000,.0000,.0000,.0000,.0000,.0000,.0000,.0000, 
     6  .0000,.0000,.0000,.0000,.0000,.0000,.0000,.0000,.0000,.0000, 
     6  .0000,.0000,.0000,.0000,.0000,.0000,.0000,.0000,.0000,.0000/ 
      DATA TGS05/ 
C  HDGV8B 
     1  .0000,.0700,.7000,.6000,.2000,.0000,.0000,.0000,.0000,.0000, 
     1  .0000,.0000,.0000,.0000,.0000,.0000,.0000,.0000,.0000,.0000, 
     2  .0000,.3000,.3000,.4000,.8000,1.000,1.000,.0000,.0000,.0000, 
     2  .0000,.0000,.0000,.0000,.0000,.0000,.0000,.0000,.0000,.0000, 
     3  .0000,.0500,.0000,.1000,.0000,.0000,.0000,.0000,.0000,.0000, 
     3  .0000,.0000,.0000,.0000,.0000,.0000,.0000,.0000,.0000,.0000, 
     4  .0000,.9500,1.000,.8500,1.000,1.000,1.000,.0000,.0000,.0000, 
     4  .0000,.0000,.0000,.0000,.0000,.0000,.0000,.0000,.0000,.0000, 
     5  .0000,.1500,.0000,.0000,.0000,1.000,1.000,.0000,.0000,.0000, 
     5  .0000,.0000,.0000,.0000,.0000,.0000,.0000,.0000,.0000,.0000, 
     6  .0000,.0000,.0000,.0000,.0000,.0000,.0000,.0000,.0000,.0000, 
     6  .0000,.0000,.0000,.0000,.0000,.0000,.0000,.0000,.0000,.0000, 
C  GAS BUS 
     1  .0000,.0700,.7000,.6000,.2000,.0000,.0000,.0000,.0000,.0000, 
     1  .0000,.0000,.0000,.0000,.0000,.0000,.0000,.0000,.0000,.0000, 
     2  .0000,.3000,.3000,.4000,.8000,1.000,1.000,.0000,.0000,.0000, 
     2  .0000,.0000,.0000,.0000,.0000,.0000,.0000,.0000,.0000,.0000, 
     3  .0000,.0500,.0000,.1000,.0000,.0000,.0000,.0000,.0000,.0000, 
     3  .0000,.0000,.0000,.0000,.0000,.0000,.0000,.0000,.0000,.0000, 
     4  .0000,.9500,1.000,.8500,1.000,1.000,1.000,.0000,.0000,.0000, 
     4  .0000,.0000,.0000,.0000,.0000,.0000,.0000,.0000,.0000,.0000, 
     5  .0000,.1500,.0000,.0000,.0000,1.000,1.000,.0000,.0000,.0000, 
     5  .0000,.0000,.0000,.0000,.0000,.0000,.0000,.0000,.0000,.0000, 
     6  .0000,.0000,.0000,.0000,.0000,.0000,.0000,.0000,.0000,.0000, 
     6  .0000,.0000,.0000,.0000,.0000,.0000,.0000,.0000,.0000,.0000/ 
C 
      DATA EGS/14*0.0/ 
C 
C  /EMECAL/: EM3W, EMOX, EMEGR and EMI(1/6/7,IB,IP) never change. 
C            EMI(2-5,IB,IP) values are computed during execution. 
C 
C  Each of the 4 arrays is ordered 1st by pollutant and then by bag. 
C  The continuation line code is thus: 
C    H = HC 
C    C = CO 
C    N = NOX 
C 
C  Within each line the values are blocked by bag 1,2,3 in that order. 
C  Within each block the values are ordered by disablement effects group 
C  for EM3W, EMOX AMD EMI or by model year group for EMEGR. 
C 
C             Bag 1               Bag 2               Bag 3 
C        AirP  CatR  Misf    AirP  CatR  Misf    AirP  CatR  Misf 
C        ----  ----  ----    ----  ----  ----    ----  ----  ---- 
      DATA EM3W/ 
     H   1.36, 1.80, 1.42,  0.765, 2.25, 1.56,  0.611, 1.81, 1.21, 
     C 31.796,16.32, 8.08, 18.207,14.11, 6.60,  18.25,11.07, 5.37, 
     N   0.00, 1.49, 0.64,   0.00, 1.36, 0.45,   0.00, 1.83, 0.74/ 
C 
      DATA EMOX/ 
     H   1.80, 2.31, 2.30,   1.37, 3.40, 2.57,   1.04, 2.95, 2.40, 
     C  34.67,41.40,46.50,  33.90,28.97,13.13,  21.28,16.06,16.62, 
     N   0.00, 0.00, 0.00,   0.00, 0.00, 0.00,   0.00, 0.00, 0.00/ 
C 
C  Column headings above only apply to EM3W and EMOX. 
C 
      DATA EMEGR/ 
     H  4*0.00,  4*0.00,  4*0.00, 
     C  4*0.00,  4*0.00,  4*0.00, 
     N  1.40,3.82,4.11,1.36, 0.96,2.63,2.68,1.19, 1.54,4.21,4.53,1.21/ 
C 
      DATA EMI/ 
     H  1.80,6*0.00,       1.37,6*0.00,       1.04,6*0.00, 
     C  34.67,6*0.00,      33.90,6*0.00,      21.28,6*0.00, 
     N  5*0.00,3.02,2.12,  5*0.00,3.46,1.85,  5*0.00,3.55,2.16/ 
C 
C  /EVPDAT/: EVPTGS and CCEMI never change.  The user will be allowed to 
C            replace the EVPAER rates with values from the One-time Data 
C            section. 
C 
C  All the evap and PCV effectiveness rates start at the neutral value, 1. 
C 
      DATA EVPAER/6*1./ 
C 
C  EVPTGS continuation line code is: 
C    1 = LDGV 
C    2 = LDGT12 
C    3 = LDGT34 
C    4 = HDGV2B, HDGV3, HDGV4, ..., HDGV8B 
C    B = GAS BUS 
C  Each line is divided into evap and PCV blocks.  Each block is filled with 
C  model year group values and dummy zeroes. 
C 
      DATA EVPTGS/ 
     1  .00,1.00, .00,   .00, .00,1.00, 
     2  .00,1.00, .00,   .00, .00,1.00, 
     2  .00,1.00, .00,   .00, .00,1.00, 
     3  .00, .05,1.00,   .00,1.00, .00, 
     3  .00, .05,1.00,   .00,1.00, .00, 
     4  .00,1.00, .00,   .00,1.00, .00, 
     4  .00,1.00, .00,   .00,1.00, .00, 
     4  .00,1.00, .00,   .00,1.00, .00, 
     4  .00,1.00, .00,   .00,1.00, .00, 
     4  .00,1.00, .00,   .00,1.00, .00, 
     4  .00,1.00, .00,   .00,1.00, .00, 
     4  .00,1.00, .00,   .00,1.00, .00, 
     4  .00,1.00, .00,   .00,1.00, .00, 
     B  .00,1.00, .00,   .00,1.00, .00/ 
C 
C  The CCEMI continuation code is the old Mobile5 IVTAM number 
C  (= EVPTGS code above). 
C 
      DATA CCEMI/ 
     1  0.00, 1.28, 1.27, 1.26, 1.24, 1.22, 1.21, 
     2  0.00, 1.28, 1.27, 1.24, 1.23, 0.00, 0.00, 
     2  0.00, 1.28, 1.27, 1.24, 1.23, 0.00, 0.00, 
     3  0.00, 1.28, 1.27, 1.24, 1.23, 0.00, 0.00, 
     3  0.00, 1.28, 1.27, 1.24, 1.23, 0.00, 0.00, 
     4  0.00, 1.34, 1.29, 0.00, 0.00, 0.00, 0.00, 
     4  0.00, 1.34, 1.29, 0.00, 0.00, 0.00, 0.00, 
     4  0.00, 1.34, 1.29, 0.00, 0.00, 0.00, 0.00, 
     4  0.00, 1.34, 1.29, 0.00, 0.00, 0.00, 0.00, 
     4  0.00, 1.34, 1.29, 0.00, 0.00, 0.00, 0.00, 
     4  0.00, 1.34, 1.29, 0.00, 0.00, 0.00, 0.00, 
     4  0.00, 1.34, 1.29, 0.00, 0.00, 0.00, 0.00, 
     4  0.00, 1.34, 1.29, 0.00, 0.00, 0.00, 0.00, 
     B  0.00, 1.34, 1.29, 0.00, 0.00, 0.00, 0.00/ 
C 
C  /LOOKUP/ & /MYCODE/ are all indices computed during execution. 
C 
      DATA IVTAM,IQG,IPG,JPGD,IHG,IGCSF/6*1/ 
      DATA LDXSY,LMYRVT,IAY,IMDXSY,IMKINK/5*1/ 
C 
C  /MYUB/: the model year upper bound pointers never change. 
C 
C  MYTGS continuation line code corresponds to the IQG code: 
C    1 = air pump only 
C    2 = air pump / catalyst 
C    3 = no air pump / catalyst  
C    4 = EGR only 
C    5 = EGR / 3-way catalyst 
C    6 = no EGR / 3-way catalyst 
C 
      DATA MYTGS1/ 
C  LDGV 
     1  1967,1971,1972,1974,1975,1977,1979,NEYR,0000,0000, 
     1  0000,0000,0000,0000,0000,0000,0000,0000,0000,0000, 
     2  1974,1976,1977,1979,1980,1982,1983,1984,1985,1986,    
     2  1987,1988,1989,1990,1991,1992,1993,1994,NEYR,0000, 
     3  1974,1975,1976,1979,1982,1983,1984,1985,1986,1987, 
     3  1988,1989,1990,1991,1992,1993,1994,NEYR,0000,0000, 
     4  1972,1973,1980,1983,1989,1990,1991,1992,NEYR,0000, 
     4  0000,0000,0000,0000,0000,0000,0000,0000,0000,0000, 
     5  1979,1980,1983,1986,1989,1990,1991,1992,1993,1994, 
     5  NEYR,0000,0000,0000,0000,0000,0000,0000,0000,0000, 
     6  1981,1982,1983,1986,1989,1990,1991,1992,1993,1994, 
     6  NEYR,0000,0000,0000,0000,0000,0000,0000,0000,0000, 
C  LDGT1 
     1  1967,1971,1972,1974,1980,1981,1982,1983,NEYR,0000, 
     1  0000,0000,0000,0000,0000,0000,0000,0000,0000,0000, 
     2  1974,1976,1978,1980,1982,1984,1985,1986,1987,1988, 
     2  1989,1990,1991,1992,1993,1994,NEYR,0000,0000,0000, 
     3  1974,1975,1976,1978,1980,1981,1983,1984,1985,1986, 
     3  1987,1988,1989,1990,1991,1992,1993,1994,NEYR,0000, 
     4  1972,1973,1978,1980,1981,1982,1983,1984,1985,1986, 
     4  1989,1990,1991,1992,NEYR,0000,0000,0000,0000,0000, 
     5  1980,1981,1982,1983,1984,1985,1986,1989,1990,1991, 
     5  1992,1993,1994,NEYR,0000,0000,0000,0000,0000,0000, 
     6  1989,1990,1991,1992,1993,1994,NEYR,0000,0000,0000, 
     6  0000,0000,0000,0000,0000,0000,0000,0000,0000,0000, 
C  LDGT2 
     1  1967,1971,1972,1974,1980,1981,1982,1983,NEYR,0000, 
     1  0000,0000,0000,0000,0000,0000,0000,0000,0000,0000, 
     2  1974,1976,1978,1980,1982,1984,1985,1986,1987,1988, 
     2  1989,1990,1991,1992,1993,1994,NEYR,0000,0000,0000, 
     3  1974,1975,1976,1978,1980,1981,1983,1984,1985,1986, 
     3  1987,1988,1989,1990,1991,1992,1993,1994,NEYR,0000, 
     4  1972,1973,1978,1980,1981,1982,1983,1984,1985,1986, 
     4  1989,1990,1991,1992,NEYR,0000,0000,0000,0000,0000, 
     5  1980,1981,1982,1983,1984,1985,1986,1989,1990,1991, 
     5  1992,1993,1994,NEYR,0000,0000,0000,0000,0000,0000, 
     6  1989,1990,1991,1992,1993,1994,NEYR,0000,0000,0000, 
     6  0000,0000,0000,0000,0000,0000,0000,0000,0000,0000/ 
C 
      DATA MYTGS2/ 
C  LDGT3 
     1  1980,1981,1982,1983,NEYR,0000,0000,0000,0000,0000, 
     1  0000,0000,0000,0000,0000,0000,0000,0000,0000,0000, 
     2  1978,1980,1982,1984,1985,1986,1987,1988,1989,1990, 
     2  1991,1992,1993,1994,NEYR,0000,0000,0000,0000,0000, 
     3  1978,1980,1981,1983,1984,1985,1986,1987,1988,1989, 
     3  1990,1991,1992,1993,1994,NEYR,0000,0000,0000,0000, 
     4  1972,1978,1980,1981,1982,1983,1984,1985,1986,1989, 
     4  1990,1991,1992,NEYR,0000,0000,0000,0000,0000,0000, 
     5  1980,1981,1982,1983,1984,1985,1986,1989,1990,1991, 
     5  1992,1993,1994,NEYR,0000,0000,0000,0000,0000,0000, 
     6  1989,1990,1991,1992,1993,1994,NEYR,0000,0000,0000, 
     6  0000,0000,0000,0000,0000,0000,0000,0000,0000,0000, 
C  LDGT4 
     1  1980,1981,1982,1983,NEYR,0000,0000,0000,0000,0000, 
     1  0000,0000,0000,0000,0000,0000,0000,0000,0000,0000, 
     2  1978,1980,1982,1984,1985,1986,1987,1988,1989,1990, 
     2  1991,1992,1993,1994,NEYR,0000,0000,0000,0000,0000, 
     3  1978,1980,1981,1983,1984,1985,1986,1987,1988,1989, 
     3  1990,1991,1992,1993,1994,NEYR,0000,0000,0000,0000, 
     4  1972,1978,1980,1981,1982,1983,1984,1985,1986,1989, 
     4  1990,1991,1992,NEYR,0000,0000,0000,0000,0000,0000, 
     5  1980,1981,1982,1983,1984,1985,1986,1989,1990,1991, 
     5  1992,1993,1994,NEYR,0000,0000,0000,0000,0000,0000, 
     6  1989,1990,1991,1992,1993,1994,NEYR,0000,0000,0000, 
     6  0000,0000,0000,0000,0000,0000,0000,0000,0000,0000, 
C  HDGV2B 
     1  1986,1989,1995,1999,2004,2009,NEYR,0000,0000,0000, 
     1  0000,0000,0000,0000,0000,0000,0000,0000,0000,0000, 
     2  1986,1989,1995,1999,2004,2009,NEYR,0000,0000,0000, 
     2  0000,0000,0000,0000,0000,0000,0000,0000,0000,0000, 
     3  1981,1986,1987,1989,1995,1999,2004,2009,NEYR,0000, 
     3  0000,0000,0000,0000,0000,0000,0000,0000,0000,0000, 
     4  1984,1985,1986,1989,1995,1999,2004,2009,NEYR,0000, 
     4  0000,0000,0000,0000,0000,0000,0000,0000,0000,0000, 
     5  1986,1989,1995,1999,2004,2009,NEYR,0000,0000,0000, 
     5  0000,0000,0000,0000,0000,0000,0000,0000,0000,0000, 
     6  1989,1995,1999,2004,2009,NEYR,0000,0000,0000,0000, 
     6  0000,0000,0000,0000,0000,0000,0000,0000,0000,0000/ 
      DATA MYTGS3/ 
C  HDGV3 
     1  1986,1989,1995,1999,2004,2009,NEYR,0000,0000,0000, 
     1  0000,0000,0000,0000,0000,0000,0000,0000,0000,0000, 
     2  1986,1989,1995,1999,2004,2009,NEYR,0000,0000,0000, 
     2  0000,0000,0000,0000,0000,0000,0000,0000,0000,0000, 
     3  1981,1986,1987,1989,1995,1999,2004,2009,NEYR,0000, 
     3  0000,0000,0000,0000,0000,0000,0000,0000,0000,0000, 
     4  1984,1985,1986,1989,1995,1999,2004,2009,NEYR,0000, 
     4  0000,0000,0000,0000,0000,0000,0000,0000,0000,0000, 
     5  1986,1989,1995,1999,2004,2009,NEYR,0000,0000,0000, 
     5  0000,0000,0000,0000,0000,0000,0000,0000,0000,0000, 
     6  1989,1995,1999,2004,2009,NEYR,0000,0000,0000,0000, 
     6  0000,0000,0000,0000,0000,0000,0000,0000,0000,0000, 
C  HDGV4 
     1  1986,1989,1995,1999,2004,2009,NEYR,0000,0000,0000, 
     1  0000,0000,0000,0000,0000,0000,0000,0000,0000,0000, 
     2  1986,1989,1995,1999,2004,2009,NEYR,0000,0000,0000, 
     2  0000,0000,0000,0000,0000,0000,0000,0000,0000,0000, 
     3  1981,1986,1987,1989,1995,1999,2004,2009,NEYR,0000, 
     3  0000,0000,0000,0000,0000,0000,0000,0000,0000,0000, 
     4  1984,1985,1986,1989,1995,1999,2004,2009,NEYR,0000, 
     4  0000,0000,0000,0000,0000,0000,0000,0000,0000,0000, 
     5  1986,1989,1995,1999,2004,2009,NEYR,0000,0000,0000, 
     5  0000,0000,0000,0000,0000,0000,0000,0000,0000,0000, 
     6  1989,1995,1999,2004,2009,NEYR,0000,0000,0000,0000, 
     6  0000,0000,0000,0000,0000,0000,0000,0000,0000,0000, 
C  HDGV5 
     1  1986,1989,1995,1999,2004,2009,NEYR,0000,0000,0000, 
     1  0000,0000,0000,0000,0000,0000,0000,0000,0000,0000, 
     2  1986,1989,1995,1999,2004,2009,NEYR,0000,0000,0000, 
     2  0000,0000,0000,0000,0000,0000,0000,0000,0000,0000, 
     3  1981,1986,1987,1989,1995,1999,2004,2009,NEYR,0000, 
     3  0000,0000,0000,0000,0000,0000,0000,0000,0000,0000, 
     4  1984,1985,1986,1989,1995,1999,2004,2009,NEYR,0000, 
     4  0000,0000,0000,0000,0000,0000,0000,0000,0000,0000, 
     5  1986,1989,1995,1999,2004,2009,NEYR,0000,0000,0000, 
     5  0000,0000,0000,0000,0000,0000,0000,0000,0000,0000, 
     6  1989,1995,1999,2004,2009,NEYR,0000,0000,0000,0000, 
     6  0000,0000,0000,0000,0000,0000,0000,0000,0000,0000/ 
      DATA MYTGS4/ 
C  HDGV6 
     1  1986,1989,1995,1999,2004,2009,NEYR,0000,0000,0000, 
     1  0000,0000,0000,0000,0000,0000,0000,0000,0000,0000, 
     2  1986,1989,1995,1999,2004,2009,NEYR,0000,0000,0000, 
     2  0000,0000,0000,0000,0000,0000,0000,0000,0000,0000, 
     3  1981,1986,1987,1989,1995,1999,2004,2009,NEYR,0000, 
     3  0000,0000,0000,0000,0000,0000,0000,0000,0000,0000, 
     4  1984,1985,1986,1989,1995,1999,2004,2009,NEYR,0000, 
     4  0000,0000,0000,0000,0000,0000,0000,0000,0000,0000, 
     5  1986,1989,1995,1999,2004,2009,NEYR,0000,0000,0000, 
     5  0000,0000,0000,0000,0000,0000,0000,0000,0000,0000, 
     6  1989,1995,1999,2004,2009,NEYR,0000,0000,0000,0000, 
     6  0000,0000,0000,0000,0000,0000,0000,0000,0000,0000, 
C  HDGV7 
     1  1986,1989,1995,1999,2004,2009,NEYR,0000,0000,0000, 
     1  0000,0000,0000,0000,0000,0000,0000,0000,0000,0000, 
     2  1986,1989,1995,1999,2004,2009,NEYR,0000,0000,0000, 
     2  0000,0000,0000,0000,0000,0000,0000,0000,0000,0000, 
     3  1981,1986,1987,1989,1995,1999,2004,2009,NEYR,0000, 
     3  0000,0000,0000,0000,0000,0000,0000,0000,0000,0000, 
     4  1984,1985,1986,1989,1995,1999,2004,2009,NEYR,0000, 
     4  0000,0000,0000,0000,0000,0000,0000,0000,0000,0000, 
     5  1986,1989,1995,1999,2004,2009,NEYR,0000,0000,0000, 
     5  0000,0000,0000,0000,0000,0000,0000,0000,0000,0000, 
     6  1989,1995,1999,2004,2009,NEYR,0000,0000,0000,0000, 
     6  0000,0000,0000,0000,0000,0000,0000,0000,0000,0000, 
C  HDGV8A 
     1  1986,1989,1995,1999,2004,2009,NEYR,0000,0000,0000, 
     1  0000,0000,0000,0000,0000,0000,0000,0000,0000,0000, 
     2  1986,1989,1995,1999,2004,2009,NEYR,0000,0000,0000, 
     2  0000,0000,0000,0000,0000,0000,0000,0000,0000,0000, 
     3  1981,1986,1987,1989,1995,1999,2004,2009,NEYR,0000, 
     3  0000,0000,0000,0000,0000,0000,0000,0000,0000,0000, 
     4  1984,1985,1986,1989,1995,1999,2004,2009,NEYR,0000, 
     4  0000,0000,0000,0000,0000,0000,0000,0000,0000,0000, 
     5  1986,1989,1995,1999,2004,2009,NEYR,0000,0000,0000, 
     5  0000,0000,0000,0000,0000,0000,0000,0000,0000,0000, 
     6  1989,1995,1999,2004,2009,NEYR,0000,0000,0000,0000, 
     6  0000,0000,0000,0000,0000,0000,0000,0000,0000,0000/ 
      DATA MYTGS5/ 
C  HDGV8B 
     1  1986,1989,1995,1999,2004,2009,NEYR,0000,0000,0000, 
     1  0000,0000,0000,0000,0000,0000,0000,0000,0000,0000, 
     2  1986,1989,1995,1999,2004,2009,NEYR,0000,0000,0000, 
     2  0000,0000,0000,0000,0000,0000,0000,0000,0000,0000, 
     3  1981,1986,1987,1989,1995,1999,2004,2009,NEYR,0000, 
     3  0000,0000,0000,0000,0000,0000,0000,0000,0000,0000, 
     4  1984,1985,1986,1989,1995,1999,2004,2009,NEYR,0000, 
     4  0000,0000,0000,0000,0000,0000,0000,0000,0000,0000, 
     5  1986,1989,1995,1999,2004,2009,NEYR,0000,0000,0000, 
     5  0000,0000,0000,0000,0000,0000,0000,0000,0000,0000, 
     6  1989,1995,1999,2004,2009,NEYR,0000,0000,0000,0000, 
     6  0000,0000,0000,0000,0000,0000,0000,0000,0000,0000, 
C  GAS BUS 
     1  1986,1989,1995,1999,2004,2009,NEYR,0000,0000,0000, 
     1  0000,0000,0000,0000,0000,0000,0000,0000,0000,0000, 
     2  1986,1989,1995,1999,2004,2009,NEYR,0000,0000,0000, 
     2  0000,0000,0000,0000,0000,0000,0000,0000,0000,0000, 
     3  1981,1986,1987,1989,1995,1999,2004,2009,NEYR,0000, 
     3  0000,0000,0000,0000,0000,0000,0000,0000,0000,0000, 
     4  1984,1985,1986,1989,1995,1999,2004,2009,NEYR,0000, 
     4  0000,0000,0000,0000,0000,0000,0000,0000,0000,0000, 
     5  1986,1989,1995,1999,2004,2009,NEYR,0000,0000,0000, 
     5  0000,0000,0000,0000,0000,0000,0000,0000,0000,0000, 
     6  1989,1995,1999,2004,2009,NEYR,0000,0000,0000,0000, 
     6  0000,0000,0000,0000,0000,0000,0000,0000,0000,0000/ 
C 
C  MYEGR keys only on vehicle type, even though EMEGR varies with pollutant 
C  and bag number.  Also, cell 2 for LDGT34 is a dummy, since there is no 
C  NOX 3.1 standard my range for LDGT34.  Otherwise, the 4 ranges for each 
C  IVTAM correspond to the no standard, 3.1, 2.0 and 1.0 cases. 
C  The continuation line code is the old Mobile5 IVTAM code: 
C    1 = LDGV 
C    2 = LDGT12 
C    3 = LDGT34 
C    4 = HDGV2B, HDGV3, HDGV4, ..., HDGV8B 
C 
      DATA MYEGR/ 
     1  1974,1976,1980,NEYR, 
     2  1974,1978,1987,NEYR, 
     2  1974,1978,1987,NEYR, 
     3  1978,1978,1987,NEYR, 
     3  1978,1978,1987,NEYR, 
     4  1989,NEYR,0000,0000, 
     4  1989,NEYR,0000,0000, 
     4  1989,NEYR,0000,0000, 
     4  1989,NEYR,0000,0000, 
     4  1989,NEYR,0000,0000, 
     4  1989,NEYR,0000,0000, 
     4  1989,NEYR,0000,0000, 
     4  1989,NEYR,0000,0000, 
     B  1989,NEYR,0000,0000/ 
C 
C  MYVTS continuation line code is also the old Mobile5 IVTAM code. Each  
C  line is divided into the evap and PCV blocks. 
C 
      DATA MYVTS/ 
     1  1970,NEYR,0000,   1962,1967,NEYR, 
     2  1970,NEYR,0000,   1962,1967,NEYR, 
     2  1970,NEYR,0000,   1962,1967,NEYR, 
     3  1970,1978,NEYR,   1967,NEYR,0000, 
     3  1970,1978,NEYR,   1967,NEYR,0000, 
     4  1984,NEYR,0000,   1967,NEYR,0000, 
     4  1984,NEYR,0000,   1967,NEYR,0000, 
     4  1984,NEYR,0000,   1967,NEYR,0000, 
     4  1984,NEYR,0000,   1967,NEYR,0000, 
     4  1984,NEYR,0000,   1967,NEYR,0000, 
     4  1984,NEYR,0000,   1967,NEYR,0000, 
     4  1984,NEYR,0000,   1967,NEYR,0000, 
     4  1984,NEYR,0000,   1967,NEYR,0000, 
     B  1984,NEYR,0000,   1967,NEYR,0000/ 
C 
C  MYCCEI continuation code is the old Mobile5 IVTAM code  
C  (= EVPTGS code above). 
C 
      DATA MYCCEI/ 
     1  1963,1970,1974,1977,1979,1980,NEYR, 
     2  1963,1977,1979,1980,NEYR,0000,0000, 
     2  1963,1977,1979,1980,NEYR,0000,0000, 
     3  1963,1977,1979,1980,NEYR,0000,0000, 
     3  1963,1977,1979,1980,NEYR,0000,0000, 
     4  1963,1980,NEYR,4*0, 
     4  1963,1980,NEYR,4*0, 
     4  1963,1980,NEYR,4*0, 
     4  1963,1980,NEYR,4*0, 
     4  1963,1980,NEYR,4*0, 
     4  1963,1980,NEYR,4*0, 
     4  1963,1980,NEYR,4*0, 
     4  1963,1980,NEYR,4*0, 
     B  1963,1980,NEYR,4*0/ 
C 
C  /OMTCF/: all values are computed during execution by BIGCFX. 
C 
      DATA TCFSTART,TCFRUN,TCF12/71*1.0/, DIST_SUM/24*0.0/ 
C 
C  /SIZCAL/: values are all computed during execution. 
C 
      DATA ATR,TGSUSE,CSIZE,CSAE/47*0.0/ 
C 
C  /SIZPAR/: values never change. 
C 
      DATA NTERMS/5,5,3,3/ 
      DATA INDXCS/1,2,3,4,5,1,4,5,6,7,2,4,6,0,0,3,5,7,0,0/ 
C 
C  /TAMEQ1/: TAMZML and TAMDR have hardcoded values that can be replaced with 
C            user supplied values in the One-time Data input section. 
C 
C  Continuation line code for TAMZML, TAMDR, and TAMA50 is the old Mobile5 
C  IVTAM number (1=LDGV, 2=LDGT12, 3=LDGT34, 4=HDGV, B=GAS BUS). 
C 
C  Each line contains the 9 disablement intercepts (TAMZML) and deterioration 
C  rates (TAMDR) for the given vehicle type, myg and I/M case. 
C 
      DATA TZLNI1/ 
C         AIRS   CATS   NCKS   TNKS   EGRS   EVAP   PCVS   CAPS   MISF 
C  Non-I/M & MYG 1 
     1   .0072, .0084, .0027,-.0286, .0412, .0192, .0107, .0334,-.0259, 
     2   .0320, .0785, .1519, .0000,-.2435,-.0195, .0000,-.1849, .1519, 
     2   .0320, .0785, .1519, .0000,-.2435,-.0195, .0000,-.1849, .1519, 
     3   .0320, .0785, .1519, .0000,-.2435,-.0195, .0000,-.1849, .1519, 
     3   .0320, .0785, .1519, .0000,-.2435,-.0195, .0000,-.1849, .1519, 
     4   .0320, .0785, .1519, .0000,-.2435,-.0195, .0000,-.1849, .1519, 
     4   .0320, .0785, .1519, .0000,-.2435,-.0195, .0000,-.1849, .1519, 
     4   .0320, .0785, .1519, .0000,-.2435,-.0195, .0000,-.1849, .1519, 
     4   .0320, .0785, .1519, .0000,-.2435,-.0195, .0000,-.1849, .1519, 
     4   .0320, .0785, .1519, .0000,-.2435,-.0195, .0000,-.1849, .1519, 
     4   .0320, .0785, .1519, .0000,-.2435,-.0195, .0000,-.1849, .1519, 
     4   .0320, .0785, .1519, .0000,-.2435,-.0195, .0000,-.1849, .1519, 
     4   .0320, .0785, .1519, .0000,-.2435,-.0195, .0000,-.1849, .1519, 
     B   .0320, .0785, .1519, .0000,-.2435,-.0195, .0000,-.1849, .1519/ 
      DATA TZLNI2/ 
C  Non-I/M & MYG 2 
     1   .0033,-.0019,-.0110, .0220, .0225, .0078,-.0036, .0110, .0110, 
     2  -.0389,-.0018, .0017, .0097,-.0204, .0009, .0145,-.0590, .0113, 
     2  -.0389,-.0018, .0017, .0097,-.0204, .0009, .0145,-.0590, .0113, 
     3  -.0389,-.0018, .0017, .0097,-.0204, .0009, .0145,-.0590, .0113, 
     3  -.0389,-.0018, .0017, .0097,-.0204, .0009, .0145,-.0590, .0113, 
     4  -.0389,-.0018, .0017, .0097,-.0204, .0009, .0145,-.0590, .0113, 
     4  -.0389,-.0018, .0017, .0097,-.0204, .0009, .0145,-.0590, .0113, 
     4  -.0389,-.0018, .0017, .0097,-.0204, .0009, .0145,-.0590, .0113, 
     4  -.0389,-.0018, .0017, .0097,-.0204, .0009, .0145,-.0590, .0113, 
     4  -.0389,-.0018, .0017, .0097,-.0204, .0009, .0145,-.0590, .0113, 
     4  -.0389,-.0018, .0017, .0097,-.0204, .0009, .0145,-.0590, .0113, 
     4  -.0389,-.0018, .0017, .0097,-.0204, .0009, .0145,-.0590, .0113, 
     4  -.0389,-.0018, .0017, .0097,-.0204, .0009, .0145,-.0590, .0113, 
     B  -.0389,-.0018, .0017, .0097,-.0204, .0009, .0145,-.0590, .0113/ 
      DATA TZLNI3/ 
C  Non-I/M & MYG 3 
     1   .0006,-.0008,-.0009, .0009, .0061,-.0013,-.0002,-.0012, .0000, 
     2  -.0039, .0024, .0075,-.0075,-.0019, .0100, .0011, .0092, .0000, 
     2  -.0039, .0024, .0075,-.0075,-.0019, .0100, .0011, .0092, .0000, 
     3  -.0039, .0024, .0075,-.0075,-.0019, .0100, .0011, .0092, .0000, 
     3  -.0039, .0024, .0075,-.0075,-.0019, .0100, .0011, .0092, .0000, 
     4  -.0039, .0024, .0075,-.0075,-.0019, .0100, .0011, .0092, .0000, 
     4  -.0039, .0024, .0075,-.0075,-.0019, .0100, .0011, .0092, .0000, 
     4  -.0039, .0024, .0075,-.0075,-.0019, .0100, .0011, .0092, .0000, 
     4  -.0039, .0024, .0075,-.0075,-.0019, .0100, .0011, .0092, .0000, 
     4  -.0039, .0024, .0075,-.0075,-.0019, .0100, .0011, .0092, .0000, 
     4  -.0039, .0024, .0075,-.0075,-.0019, .0100, .0011, .0092, .0000, 
     4  -.0039, .0024, .0075,-.0075,-.0019, .0100, .0011, .0092, .0000, 
     4  -.0039, .0024, .0075,-.0075,-.0019, .0100, .0011, .0092, .0000, 
     B  -.0039, .0024, .0075,-.0075,-.0019, .0100, .0011, .0092, .0000/ 
C 
      DATA TZLIM1/ 
C         AIRS   CATS   NCKS   TNKS   EGRS   EVAP   PCVS   CAPS   MISF 
C  I/M & MYG 1 
     1   .0668, .0067, .0379, .0528, .0744, .0099, .0048,-.0116, .0907, 
     2   .0615, .1037, .1376,-.0885, .0808, .0695,-.0156, .0471, .0491, 
     2   .0615, .1037, .1376,-.0885, .0808, .0695,-.0156, .0471, .0491, 
     3   .0615, .1037, .1376,-.0885, .0808, .0695,-.0156, .0471, .0491, 
     3   .0615, .1037, .1376,-.0885, .0808, .0695,-.0156, .0471, .0491, 
     4   .0615, .1037, .1376,-.0885, .0808, .0695,-.0156, .0471, .0491, 
     4   .0615, .1037, .1376,-.0885, .0808, .0695,-.0156, .0471, .0491, 
     4   .0615, .1037, .1376,-.0885, .0808, .0695,-.0156, .0471, .0491, 
     4   .0615, .1037, .1376,-.0885, .0808, .0695,-.0156, .0471, .0491, 
     4   .0615, .1037, .1376,-.0885, .0808, .0695,-.0156, .0471, .0491, 
     4   .0615, .1037, .1376,-.0885, .0808, .0695,-.0156, .0471, .0491, 
     4   .0615, .1037, .1376,-.0885, .0808, .0695,-.0156, .0471, .0491, 
     4   .0615, .1037, .1376,-.0885, .0808, .0695,-.0156, .0471, .0491, 
     B   .0615, .1037, .1376,-.0885, .0808, .0695,-.0156, .0471, .0491/ 
      DATA TZLIM2/ 
C  I/M & MYG 2 
     1   .0066,-.0152,-.0005, .0112, .0203,-.0158,-.0028,-.0173, .0108, 
     2  -.0318,-.0062, .0201,-.0580, .0538,-.0126,-.0261,-.0314,-.0379, 
     2  -.0318,-.0062, .0201,-.0580, .0538,-.0126,-.0261,-.0314,-.0379, 
     3  -.0318,-.0062, .0201,-.0580, .0538,-.0126,-.0261,-.0314,-.0379, 
     3  -.0318,-.0062, .0201,-.0580, .0538,-.0126,-.0261,-.0314,-.0379, 
     4  -.0318,-.0062, .0201,-.0580, .0538,-.0126,-.0261,-.0314,-.0379, 
     4  -.0318,-.0062, .0201,-.0580, .0538,-.0126,-.0261,-.0314,-.0379, 
     4  -.0318,-.0062, .0201,-.0580, .0538,-.0126,-.0261,-.0314,-.0379, 
     4  -.0318,-.0062, .0201,-.0580, .0538,-.0126,-.0261,-.0314,-.0379, 
     4  -.0318,-.0062, .0201,-.0580, .0538,-.0126,-.0261,-.0314,-.0379, 
     4  -.0318,-.0062, .0201,-.0580, .0538,-.0126,-.0261,-.0314,-.0379, 
     4  -.0318,-.0062, .0201,-.0580, .0538,-.0126,-.0261,-.0314,-.0379, 
     4  -.0318,-.0062, .0201,-.0580, .0538,-.0126,-.0261,-.0314,-.0379, 
     B  -.0318,-.0062, .0201,-.0580, .0538,-.0126,-.0261,-.0314,-.0379/ 
      DATA TZLIM3/ 
C  I/M & MYG 3 
     1   .0053, .0000, .0000, .0000, .0012,-.0013,-.0018,-.0012, .0000, 
     2  -.0262, .0006,-.0122, .0197,-.0021, .0100, .0000, .0092, .0075, 
     2  -.0262, .0006,-.0122, .0197,-.0021, .0100, .0000, .0092, .0075, 
     3  -.0262, .0006,-.0122, .0197,-.0021, .0100, .0000, .0092, .0075, 
     3  -.0262, .0006,-.0122, .0197,-.0021, .0100, .0000, .0092, .0075, 
     4  -.0262, .0006,-.0122, .0197,-.0021, .0100, .0000, .0092, .0075, 
     4  -.0262, .0006,-.0122, .0197,-.0021, .0100, .0000, .0092, .0075, 
     4  -.0262, .0006,-.0122, .0197,-.0021, .0100, .0000, .0092, .0075, 
     4  -.0262, .0006,-.0122, .0197,-.0021, .0100, .0000, .0092, .0075, 
     4  -.0262, .0006,-.0122, .0197,-.0021, .0100, .0000, .0092, .0075, 
     4  -.0262, .0006,-.0122, .0197,-.0021, .0100, .0000, .0092, .0075, 
     4  -.0262, .0006,-.0122, .0197,-.0021, .0100, .0000, .0092, .0075, 
     4  -.0262, .0006,-.0122, .0197,-.0021, .0100, .0000, .0092, .0075, 
     B  -.0262, .0006,-.0122, .0197,-.0021, .0100, .0000, .0092, .0075/ 
C 
      DATA TDRNI1/ 
C        AIRS   CATS   NCKS    TNKS   EGRS   EVAP   PCVS   CAPS   MISF 
C  Non-I/M & MYG 1 & vmt <= 50K 
     1   .0261, .0193, .0242, .0100, .0347, .0097, .0022, .0096, .0342, 
     2   .0767, .0444, .0297, .0000, .1198, .0438, .0073, .0839, .0297, 
     2   .0767, .0444, .0297, .0000, .1198, .0438, .0073, .0839, .0297, 
     3   .0767, .0444, .0297, .0000, .1198, .0438, .0073, .0839, .0297, 
     3   .0767, .0444, .0297, .0000, .1198, .0438, .0073, .0839, .0297, 
     4   .0767, .0444, .0297, .0000, .1198, .0438, .0073, .0839, .0297, 
     4   .0767, .0444, .0297, .0000, .1198, .0438, .0073, .0839, .0297, 
     4   .0767, .0444, .0297, .0000, .1198, .0438, .0073, .0839, .0297, 
     4   .0767, .0444, .0297, .0000, .1198, .0438, .0073, .0839, .0297, 
     4   .0767, .0444, .0297, .0000, .1198, .0438, .0073, .0839, .0297, 
     4   .0767, .0444, .0297, .0000, .1198, .0438, .0073, .0839, .0297, 
     4   .0767, .0444, .0297, .0000, .1198, .0438, .0073, .0839, .0297, 
     4   .0767, .0444, .0297, .0000, .1198, .0438, .0073, .0839, .0297, 
     B   .0767, .0444, .0297, .0000, .1198, .0438, .0073, .0839, .0297/ 
      DATA TDRNI2/ 
C  Non-I/M & MYG 2 & vmt <= 50K 
     1   .0048, .0038, .0067,-.0033, .0040, .0013, .0022, .0046, .0034, 
     2   .0299, .0258, .0188, .0015, .0258, .0057, .0003, .0270, .0203, 
     2   .0299, .0258, .0188, .0015, .0258, .0057, .0003, .0270, .0203, 
     3   .0299, .0258, .0188, .0015, .0258, .0057, .0003, .0270, .0203, 
     3   .0299, .0258, .0188, .0015, .0258, .0057, .0003, .0270, .0203, 
     4   .0299, .0258, .0188, .0015, .0258, .0057, .0003, .0270, .0203, 
     4   .0299, .0258, .0188, .0015, .0258, .0057, .0003, .0270, .0203, 
     4   .0299, .0258, .0188, .0015, .0258, .0057, .0003, .0270, .0203, 
     4   .0299, .0258, .0188, .0015, .0258, .0057, .0003, .0270, .0203, 
     4   .0299, .0258, .0188, .0015, .0258, .0057, .0003, .0270, .0203, 
     4   .0299, .0258, .0188, .0015, .0258, .0057, .0003, .0270, .0203, 
     4   .0299, .0258, .0188, .0015, .0258, .0057, .0003, .0270, .0203, 
     4   .0299, .0258, .0188, .0015, .0258, .0057, .0003, .0270, .0203, 
     B   .0299, .0258, .0188, .0015, .0258, .0057, .0003, .0270, .0203/ 
      DATA TDRNI3/ 
C  Non-I/M & MYG 3 & vmt <= 50K 
     1   .0028, .0014, .0006, .0019, .0017, .0026, .0013, .0055, .0025, 
     2   .0088, .0035, .0006, .0051, .0053, .0011, .0009, .0047, .0057, 
     2   .0088, .0035, .0006, .0051, .0053, .0011, .0009, .0047, .0057, 
     3   .0088, .0035, .0006, .0051, .0053, .0011, .0009, .0047, .0057, 
     3   .0088, .0035, .0006, .0051, .0053, .0011, .0009, .0047, .0057, 
     4   .0088, .0035, .0006, .0051, .0053, .0011, .0009, .0047, .0057, 
     4   .0088, .0035, .0006, .0051, .0053, .0011, .0009, .0047, .0057, 
     4   .0088, .0035, .0006, .0051, .0053, .0011, .0009, .0047, .0057, 
     4   .0088, .0035, .0006, .0051, .0053, .0011, .0009, .0047, .0057, 
     4   .0088, .0035, .0006, .0051, .0053, .0011, .0009, .0047, .0057, 
     4   .0088, .0035, .0006, .0051, .0053, .0011, .0009, .0047, .0057, 
     4   .0088, .0035, .0006, .0051, .0053, .0011, .0009, .0047, .0057, 
     4   .0088, .0035, .0006, .0051, .0053, .0011, .0009, .0047, .0057, 
     B   .0088, .0035, .0006, .0051, .0053, .0011, .0009, .0047, .0057/ 
C 
      DATA TDRIM1/ 
C        AIRS   CATS   NCKS    TNKS   EGRS   EVAP   PCVS   CAPS   MISF 
C  I/M & MYG 1 & vmt <= 50K 
     1   .0055, .0084, .0089,-.0038, .0105, .0085, .0023, .0156, .0051, 
     2   .0304, .0203, .0210, .0178, .0114, .0093, .0040, .0150, .0388, 
     2   .0304, .0203, .0210, .0178, .0114, .0093, .0040, .0150, .0388, 
     3   .0304, .0203, .0210, .0178, .0114, .0093, .0040, .0150, .0388, 
     3   .0304, .0203, .0210, .0178, .0114, .0093, .0040, .0150, .0388, 
     4   .0304, .0203, .0210, .0178, .0114, .0093, .0040, .0150, .0388, 
     4   .0304, .0203, .0210, .0178, .0114, .0093, .0040, .0150, .0388, 
     4   .0304, .0203, .0210, .0178, .0114, .0093, .0040, .0150, .0388, 
     4   .0304, .0203, .0210, .0178, .0114, .0093, .0040, .0150, .0388, 
     4   .0304, .0203, .0210, .0178, .0114, .0093, .0040, .0150, .0388, 
     4   .0304, .0203, .0210, .0178, .0114, .0093, .0040, .0150, .0388, 
     4   .0304, .0203, .0210, .0178, .0114, .0093, .0040, .0150, .0388, 
     4   .0304, .0203, .0210, .0178, .0114, .0093, .0040, .0150, .0388, 
     B   .0304, .0203, .0210, .0178, .0114, .0093, .0040, .0150, .0388/ 
      DATA TDRIM2/ 
C  I/M & MYG 2 & vmt <= 50K 
     1   .0029, .0060, .0009, .0008, .0014, .0085, .0021, .0122, .0017, 
     2   .0148, .0109, .0031, .0198, .0000, .0128, .0079, .0184, .0229, 
     2   .0148, .0109, .0031, .0198, .0000, .0128, .0079, .0184, .0229, 
     3   .0148, .0109, .0031, .0198, .0000, .0128, .0079, .0184, .0229, 
     3   .0148, .0109, .0031, .0198, .0000, .0128, .0079, .0184, .0229, 
     4   .0148, .0109, .0031, .0198, .0000, .0128, .0079, .0184, .0229, 
     4   .0148, .0109, .0031, .0198, .0000, .0128, .0079, .0184, .0229, 
     4   .0148, .0109, .0031, .0198, .0000, .0128, .0079, .0184, .0229, 
     4   .0148, .0109, .0031, .0198, .0000, .0128, .0079, .0184, .0229, 
     4   .0148, .0109, .0031, .0198, .0000, .0128, .0079, .0184, .0229, 
     4   .0148, .0109, .0031, .0198, .0000, .0128, .0079, .0184, .0229, 
     4   .0148, .0109, .0031, .0198, .0000, .0128, .0079, .0184, .0229, 
     4   .0148, .0109, .0031, .0198, .0000, .0128, .0079, .0184, .0229, 
     B   .0148, .0109, .0031, .0198, .0000, .0128, .0079, .0184, .0229/ 
      DATA TDRIM3/ 
C  I/M & MYG 3 & vmt <= 50K 
     1   .0011, .0003, .0002, .0034, .0006, .0026, .0014, .0055, .0036, 
     2   .0164, .0043, .0081, .0005, .0027, .0011, .0000, .0047, .0086, 
     2   .0164, .0043, .0081, .0005, .0027, .0011, .0000, .0047, .0086, 
     3   .0164, .0043, .0081, .0005, .0027, .0011, .0000, .0047, .0086, 
     3   .0164, .0043, .0081, .0005, .0027, .0011, .0000, .0047, .0086, 
     4   .0164, .0043, .0081, .0005, .0027, .0011, .0000, .0047, .0086, 
     4   .0164, .0043, .0081, .0005, .0027, .0011, .0000, .0047, .0086, 
     4   .0164, .0043, .0081, .0005, .0027, .0011, .0000, .0047, .0086, 
     4   .0164, .0043, .0081, .0005, .0027, .0011, .0000, .0047, .0086, 
     4   .0164, .0043, .0081, .0005, .0027, .0011, .0000, .0047, .0086, 
     4   .0164, .0043, .0081, .0005, .0027, .0011, .0000, .0047, .0086, 
     4   .0164, .0043, .0081, .0005, .0027, .0011, .0000, .0047, .0086, 
     4   .0164, .0043, .0081, .0005, .0027, .0011, .0000, .0047, .0086, 
     B   .0164, .0043, .0081, .0005, .0027, .0011, .0000, .0047, .0086/ 
C 
      DATA TA5NI1/ 
C        AIRS   CATS   NCKS    TNKS   EGRS   EVAP   PCVS   CAPS   MISF 
C  Non-I/M & MYG 1 & vmt > 50K 
     1  1.5794, .9731, 1.0913,-.0211, .3245,2.0404, 1.00, 2.3463, .7653, 
     2   .2207, .9810, 1.6721, .0000, .0071, .3830, 1.00, .1806, 1.6721, 
     2   .2207, .9810, 1.6721, .0000, .0071, .3830, 1.00, .1806, 1.6721, 
     3   .2207, .9810, 1.6721, .0000, .0071, .3830, 1.00, .1806, 1.6721, 
     3   .2207, .9810, 1.6721, .0000, .0071, .3830, 1.00, .1806, 1.6721, 
     4   .2207, .9810, 1.6721, .0000, .0071, .3830, 1.00, .1806, 1.6721, 
     4   .2207, .9810, 1.6721, .0000, .0071, .3830, 1.00, .1806, 1.6721, 
     4   .2207, .9810, 1.6721, .0000, .0071, .3830, 1.00, .1806, 1.6721, 
     4   .2207, .9810, 1.6721, .0000, .0071, .3830, 1.00, .1806, 1.6721, 
     4   .2207, .9810, 1.6721, .0000, .0071, .3830, 1.00, .1806, 1.6721, 
     4   .2207, .9810, 1.6721, .0000, .0071, .3830, 1.00, .1806, 1.6721, 
     4   .2207, .9810, 1.6721, .0000, .0071, .3830, 1.00, .1806, 1.6721, 
     4   .2207, .9810, 1.6721, .0000, .0071, .3830, 1.00, .1806, 1.6721, 
     B   .2207, .9810, 1.6721, .0000, .0071, .3830, 1.00, .1806, 1.6721/ 
      DATA TA5NI2/ 
C  Non-I/M & MYG 2 & vmt > 50K 
     1  7.4853, 3.6154, .9121,-.8487,2.3522,1.0129, .2863,2.8154,2.6554, 
     2  1.0969, .5767, .5832, 2.5794, .5597,3.8992,1.4240, .5891, .7318, 
     2  1.0969, .5767, .5832, 2.5794, .5597,3.8992,1.4240, .5891, .7318, 
     3  1.0969, .5767, .5832, 2.5794, .5597,3.8992,1.4240, .5891, .7318, 
     3  1.0969, .5767, .5832, 2.5794, .5597,3.8992,1.4240, .5891, .7318, 
     4  1.0969, .5767, .5832, 2.5794, .5597,3.8992,1.4240, .5891, .7318, 
     4  1.0969, .5767, .5832, 2.5794, .5597,3.8992,1.4240, .5891, .7318, 
     4  1.0969, .5767, .5832, 2.5794, .5597,3.8992,1.4240, .5891, .7318, 
     4  1.0969, .5767, .5832, 2.5794, .5597,3.8992,1.4240, .5891, .7318, 
     4  1.0969, .5767, .5832, 2.5794, .5597,3.8992,1.4240, .5891, .7318, 
     4  1.0969, .5767, .5832, 2.5794, .5597,3.8992,1.4240, .5891, .7318, 
     4  1.0969, .5767, .5832, 2.5794, .5597,3.8992,1.4240, .5891, .7318, 
     4  1.0969, .5767, .5832, 2.5794, .5597,3.8992,1.4240, .5891, .7318, 
     B  1.0969, .5767, .5832, 2.5794, .5597,3.8992,1.4240, .5891, .7318/ 
      DATA TA5NI3/ 
C  Non-I/M & MYG 3 & vmt > 50K 
     1  3.7188, 1.0400, 1.00, 1.000, 1.3764, .3976, .3986, .2305, 1.000, 
     2  1.5120, 4.0011, 1.00, 1.00, 2.0764,  10.39, 1.00, 3.2234, 1.000, 
     2  1.5120, 4.0011, 1.00, 1.00, 2.0764,  10.39, 1.00, 3.2234, 1.000, 
     3  1.5120, 4.0011, 1.00, 1.00, 2.0764,  10.39, 1.00, 3.2234, 1.000, 
     3  1.5120, 4.0011, 1.00, 1.00, 2.0764,  10.39, 1.00, 3.2234, 1.000, 
     4  1.5120, 4.0011, 1.00, 1.00, 2.0764,  10.39, 1.00, 3.2234, 1.000, 
     4  1.5120, 4.0011, 1.00, 1.00, 2.0764,  10.39, 1.00, 3.2234, 1.000, 
     4  1.5120, 4.0011, 1.00, 1.00, 2.0764,  10.39, 1.00, 3.2234, 1.000, 
     4  1.5120, 4.0011, 1.00, 1.00, 2.0764,  10.39, 1.00, 3.2234, 1.000, 
     4  1.5120, 4.0011, 1.00, 1.00, 2.0764,  10.39, 1.00, 3.2234, 1.000, 
     4  1.5120, 4.0011, 1.00, 1.00, 2.0764,  10.39, 1.00, 3.2234, 1.000, 
     4  1.5120, 4.0011, 1.00, 1.00, 2.0764,  10.39, 1.00, 3.2234, 1.000, 
     4  1.5120, 4.0011, 1.00, 1.00, 2.0764,  10.39, 1.00, 3.2234, 1.000, 
     B  1.5120, 4.0011, 1.00, 1.00, 2.0764,  10.39, 1.00, 3.2234, 1.000/ 
C 
      DATA TA5IM1/ 
C        AIRS   CATS   NCKS    TNKS   EGRS   EVAP   PCVS   CAPS   MISF 
C  I/M & MYG 1 & vmt > 50K 
     1  5.6555, 1.000, 1.000, 1.000, 1.000, 1.4620, 1.000, .8646, 1.000, 
     2   .9643, 1.000, 1.000, 1.000, 1.000, 1.000, 3.1065, 1.000, 1.000, 
     2   .9643, 1.000, 1.000, 1.000, 1.000, 1.000, 3.1065, 1.000, 1.000, 
     3   .9643, 1.000, 1.000, 1.000, 1.000, 1.000, 3.1065, 1.000, 1.000, 
     3   .9643, 1.000, 1.000, 1.000, 1.000, 1.000, 3.1065, 1.000, 1.000, 
     4   .9643, 1.000, 1.000, 1.000, 1.000, 1.000, 3.1065, 1.000, 1.000, 
     4   .9643, 1.000, 1.000, 1.000, 1.000, 1.000, 3.1065, 1.000, 1.000, 
     4   .9643, 1.000, 1.000, 1.000, 1.000, 1.000, 3.1065, 1.000, 1.000, 
     4   .9643, 1.000, 1.000, 1.000, 1.000, 1.000, 3.1065, 1.000, 1.000, 
     4   .9643, 1.000, 1.000, 1.000, 1.000, 1.000, 3.1065, 1.000, 1.000, 
     4   .9643, 1.000, 1.000, 1.000, 1.000, 1.000, 3.1065, 1.000, 1.000, 
     4   .9643, 1.000, 1.000, 1.000, 1.000, 1.000, 3.1065, 1.000, 1.000, 
     4   .9643, 1.000, 1.000, 1.000, 1.000, 1.000, 3.1065, 1.000, 1.000, 
     B   .9643, 1.000, 1.000, 1.000, 1.000, 1.000, 3.1065, 1.000, 1.000/ 
      DATA TA5IM2/ 
C  I/M & MYG 2 & vmt > 50K 
     1  5.4872, 1.000, .4492,-.3720, 1.1930, .0310, 1.000, .1665, .0800, 
     2  1.0000, 1.000, 1.000, .7472, .0000, .0119, 1.0000, .0344, .7812, 
     2  1.0000, 1.000, 1.000, .7472, .0000, .0119, 1.0000, .0344, .7812, 
     3  1.0000, 1.000, 1.000, .7472, .0000, .0119, 1.0000, .0344, .7812, 
     3  1.0000, 1.000, 1.000, .7472, .0000, .0119, 1.0000, .0344, .7812, 
     4  1.0000, 1.000, 1.000, .7472, .0000, .0119, 1.0000, .0344, .7812, 
     4  1.0000, 1.000, 1.000, .7472, .0000, .0119, 1.0000, .0344, .7812, 
     4  1.0000, 1.000, 1.000, .7472, .0000, .0119, 1.0000, .0344, .7812, 
     4  1.0000, 1.000, 1.000, .7472, .0000, .0119, 1.0000, .0344, .7812, 
     4  1.0000, 1.000, 1.000, .7472, .0000, .0119, 1.0000, .0344, .7812, 
     4  1.0000, 1.000, 1.000, .7472, .0000, .0119, 1.0000, .0344, .7812, 
     4  1.0000, 1.000, 1.000, .7472, .0000, .0119, 1.0000, .0344, .7812, 
     4  1.0000, 1.000, 1.000, .7472, .0000, .0119, 1.0000, .0344, .7812, 
     B  1.0000, 1.000, 1.000, .7472, .0000, .0119, 1.0000, .0344, .7812/ 
      DATA TA5IM3/ 
C  I/M & MYG 3 & vmt > 50K 
     1   .5021, 1.000, 1.000, 1.000, 1.5360, .3976,2.5328, .2305, 1.000, 
     2  1.000, 5.3523, 1.000, 1.000, 2.0765, 10.39, .0000,3.2234, 1.000, 
     2  1.000, 5.3523, 1.000, 1.000, 2.0765, 10.39, .0000,3.2234, 1.000, 
     3  1.000, 5.3523, 1.000, 1.000, 2.0765, 10.39, .0000,3.2234, 1.000, 
     3  1.000, 5.3523, 1.000, 1.000, 2.0765, 10.39, .0000,3.2234, 1.000, 
     4  1.000, 5.3523, 1.000, 1.000, 2.0765, 10.39, .0000,3.2234, 1.000, 
     4  1.000, 5.3523, 1.000, 1.000, 2.0765, 10.39, .0000,3.2234, 1.000, 
     4  1.000, 5.3523, 1.000, 1.000, 2.0765, 10.39, .0000,3.2234, 1.000, 
     4  1.000, 5.3523, 1.000, 1.000, 2.0765, 10.39, .0000,3.2234, 1.000, 
     4  1.000, 5.3523, 1.000, 1.000, 2.0765, 10.39, .0000,3.2234, 1.000, 
     4  1.000, 5.3523, 1.000, 1.000, 2.0765, 10.39, .0000,3.2234, 1.000, 
     4  1.000, 5.3523, 1.000, 1.000, 2.0765, 10.39, .0000,3.2234, 1.000, 
     4  1.000, 5.3523, 1.000, 1.000, 2.0765, 10.39, .0000,3.2234, 1.000, 
     B  1.000, 5.3523, 1.000, 1.000, 2.0765, 10.39, .0000,3.2234, 1.000/ 
C 
C  /TAMEQ2/: tampering rate equation parameters used to select and apply 
C            tampering rate slope and intercept. 
C 
C  Continuation line code for MYGTAM is the old Mobile5 IVTAM code (1=LDGV, 
C  2=LDGT12, 3=LDGT34, 4=HDGV). The code for the new GAS BUS category is B. 
C 
      DATA MYGTAM/ 
C                          MYG 
C                    1      2      3 
     1             1980,  1983,  NEYR, 
     2             1980,  1983,  NEYR, 
     2             1980,  1983,  NEYR, 
     3             1980,  1983,  NEYR, 
     3             1980,  1983,  NEYR, 
     4             1980,  1983,  NEYR, 
     4             1980,  1983,  NEYR, 
     4             1980,  1983,  NEYR, 
     4             1980,  1983,  NEYR, 
     4             1980,  1983,  NEYR, 
     4             1980,  1983,  NEYR, 
     4             1980,  1983,  NEYR, 
     4             1980,  1983,  NEYR, 
     B             1980,  1983,  NEYR/ 
      DATA F50K/5.0/ 
C 
C  /TAMEQ3/: contains the cumulative mileage segments calculated in DISIMP. 
C 
      DATA CUMCAP/13.0/,CUMALL,CUMNIM,CUMIM,GT50KA,GT50KN,GT50KI/6*0.0/ 
C 
C  /TAMEQ4/: contains BTR, whose rates are calculated during execution. 
C 
      DATA BTR/18*0./ 
C 
C  /TAMOUT/: all values are computed during execution - TAMBAG by BAGEME 
C            and THS, TDU & TCC by EVPEME. TEX and TAMBAS are computed by 
C            TEX1. 
C 
C 
      DATA TAMBAG/3150*0.0/,THS,TDU,TCC/1400*0.0/ 
      DATA TAMBAS/3150*0.0/ 
C 
      DATA TEX/4200*0.0/ 
C 
C  /TAMOU2/: all values are computed during execution by EVPEME. 
C 
      DATA TOB,GCONLY/700*0.0/ 
C 
C  /TAMPB1/: LPOD is fixed at fixed-by-law 1995, PBBTR is computed in SAVEPB. 
C 
      DATA LPOD/1995/, PBBTR/2100*0.0/ 
C 
      END 
