      BLOCK DATA BD27 
C 
C  BLOCK DATA Subprogram 27: /ALUHIN/,/ALUOUT/,/BASEQ3/,/CEVBMY/, 
C                            /CUMCOM/,/FLAGS1/,/FLAGS2/,/FLAGS3/,/FLAGS4/, 
C                            /IMPAR1/,/IMPAR2/,/IM240P/,/MYRCAL/,/OFFSET/, 
C                            /REGISF/,/RESUL1/,/RESUL2/,/RESUL3/,/RESUL4/ 
C 
C  BLOCK DATA 27's common blocks did not need to be initialized, so were 
C  assigned dummy zeroes, to facilitate debugging, to prevent control statement 
C  failure when compiled under a FORTRAN that does not employ the left to right 
C  rule, etc. 
C 
C  Changes: (Last change first) 
C 
C  09 May 01 @EPA-elg Deactivated POST96IM Flag by setting it to 1. 
C  16 Apr 01 @EPA-elg Added initialization of POST96IM Flag to zero. 
C  19 Feb 01 AIR Task 11: Added IEFPM, particulates. 
C  20 Nov 00 @EPA-ELG Added the No Clean Air Act Flag Initialization to zero. 
C 27 Nov 00 @EPA-bag bug#301 Removed FVMTFLG, HVMTFLG, SVMTFLG 
C  06 Oct 00 @EPA-BG  Bug 155, Removed IDLFLG 
C  25 Sep 00 AIR Task 07: Bug Fixes 256 and 260: Added initialization of 
C            of global AC parameters. 
C  15 Sep 00 AIR Task 08: Added EVIMFLAG 
C  29 Aug 00 @EPA-elg bug241 Removed PCTREM array. 
C  07 Aug 00 @EPA-djb Bug 218, added PCTREM initialization to 1. 
C  27 Jul 00 @EPA-djb Removed HOUREF and added PCTREM.  Initialized MOFFMTH to 0. 
C  18 Jul 00 @EPA-bag bug#54 removed TEMFLG 
C  08 May 00 AIR Task 03: Added the array OFFSET. 
C  15 Jun 00 @EPA-dj   Initialized MOFFMTH Array to 0.0 rather than 1.0. 
C  25 Apr 00 @EPA-elg  Added MOFFMTH Array. 
C  10 Apr 00 AIR Task 02: Removed PRTFLG (which is now controlled via DBSELPOL) 
C  17 Apr 00 @EPA-djb Initialize USRGSF. 
C  30 Jan 00 @DynTel-BG  1-005  Removed DAILY_DIURNAL.I 
C  13 Jan 00 @DynTel-ddj 1-013  Removed ALHFLG, drybulb and wetbulb 
C            temperatures and towing and trailer parameters.  Also 
C            set default for absolute humidity (ABSHUM) to 75 and 
C            fraction of AC usage to 1.0. 
C  12 Nov 99 @DynTel-BG  1-006  Cleaned up evaporative emissions arrays. Initialized HCFLAG to 1 
C            rather than 0. 
C  14 Oct 99 @DynTel-ddj 2-714  Removed array OFFCO.  OFFCO is now a 
C            local variable in MAINTCF. 
C   4 Oct 99 @DynTel-ddj 2-714  Added IFDS indes to OFFMTH. 
C  28 Sep 99 @DynTel-ddj for EPA  Added ISTD index to OFFMTH. 
C  28 Jul 99 @Dyntel-bag 2-663  Moved I&M common block data to a separate 
C            I&M block data file, BD_IM1.FOR. 
C  18 May 99 @Dyntel-HXQ 2-650  Delete array BMYTPD 
C  30 Mar 99 @Dyntel-HXQ 2-698  Initialize speed correction related 
C            coefficient 
C  29 Mar 99 @DynTel-ddj 2-696 Added hourly indices to OFFCO and ALHRET. 
C  19 Mar 99 @DynTel-MLA 2-701 Added /RESUL4/ 
C   3 Mar 99 @DynTel-bag 2-664 Increased hourly indices to 24. Removed RNLFLG. 
C  12 Jan 99 @DynTel-bag 2-664 Added hourly running loss array, RNLOSSHR, and 
C            modified arrays for diurnal and resting losses, added hot soak 
C            variables HSK(14) and EVAPHRLY from HOTSK1.I . 
C  20 Nov 1998 @DynTel-MLA 2-663 Removed RDTTCF (Not used anywhere in code) 
C  30 Oct 98 @DynTel-bag 2-677 Added hourly diurnal and resting loss arrays 
C  29 Sep 98 @DynTel-HQ 2-678 Conversion Factor 
C  24 Sep 98 @DynTel-MLA 2-670  Added HSHOUR ( IH, IVGAS ) for the hourly 
C            hot soak emissions. 
C 
C  Common block array subscripts: 
C 
C  ALHRET(25,6,3,24)    -  ALHRET ( JDX, IG, IP, IH ) 
C  BMYMPD(25,28)        -  BMYMPD ( IG, IV ) 
C  BMYTPD(25,28)        -  BMYTPD ( IG, IV ) 
C  BMTPDF(25,28)        -  BMTPDF ( IG, IV ) 
C  CUMMIL(25,28)        -  CUMMIL ( JDX, IV ) 
C  EFEVAP(29)           -  EFEVAP ( IV9 ) 
C  EFEXH(29)            -  EFEXH  ( IV9 ) 
C  EFFTP(3,29)          -  EFFTP  ( IP, IV9 ) 
C  EFIDLE(3,29)         -  EFIDLE ( IP, IV9 ) 
C  REFL_GM_MI(29)       -  REFL_GM_MI ( IV9 ) 
C  EFRUNL(29)           -  EFRUNL ( IV9 ) 
C  G_MI_EVP(6,29)       -  G_MI_EVP (IEVP, IV9 ) 
C  GSFVCT(28)           -  GSFVCT ( IV ) 
C  MEVMYR(25,28)        -  MEVMYR ( JDX, IV ) 
C  MOFFMTH(25,28,2,10,3)-  MOFFMTH ( JDX, IV, ISR, ISTD, IFDS ) 
C  OFFMTH(25,28,2,10,3) -  OFFMTH ( JDX, IV, ISR, ISTD, IFDS ) 
C  RLGGAL(29)           -  RLGGAL ( IV9 ) 
C  RSTGPH(29)           -  RSTGPH ( IV9 ) 
C  RUNEF(3,29)          -  RUNEF ( IP, IV ) 
C  STARTEF(3,29)        -  STARTEF ( IP, IV ) 
C  TF(25,28)            -  TF     ( IDX, IV ) 
C  TFMYM(25,28)         -  TFMYM  ( JDX, IV ) 
C  TOTEVAP(29)          -  TOTEVAP ( IV ) 
C  USRGSF(25,14)        -  USRGSF ( IDX, IVTEMP ) 
C  VFTP(3)              -  VFTP   ( IP ) 
C  VIDLE(3)             -  VIDLE  ( IP ) 
C  XMYM(25,28)          -  XMYM   ( JDX, IV ) 
C 
C  Common block dictionary: 
C 
C   Name   Type              Description 
C  ------  ----  -------------------------------------------------- 
C  /ALUHIN/: 
C 
C  AC       R    a/c usage factor 
C  ABSHUM   R    absolute humidity level 
C  GLAC     R    global a/c usage factor 
C  GSUNRISE R    global sunrise 
C  GSUNSET  R    global sunset 
C  GPEAKSUN R    global peaksun (1) begin, (2) end 
C  GCLD_CVR R    global cloud cover 
C  GLAHUM   R    global absolute humidity 
C 
C  /ALUOUT/: 
C 
C  ALHRET   R    combined correction factor derived from user supplied data 
C                on a/c, extra load & trailer towing usage and on ambient 
C                humidity & temperature. 
C 
C  /BASEQ3/: 
C 
C  NEWVEH   I    new emission rate's vehicle type code 
C  NEWMYF   I    new emission rate's first my of affected range 
C  NEWMYL   I    new emission rate's last my of affected range 
C  ZMLNEW   R    new emission rate's zero mile level (intercept) 
C  DRNEW    R    new emission rate's deterioration rate (slope) 
C  KINKS    I    flag signaling whether or not any kinked er's have been read 
C                for 1981+ LDGV HC/CO:  1 = no  2 = yes 
C 
C  /CEVBMY/: 
C 
C  BMYMPD   R    by model year miles per day (MPD), calculated in CEVMPD 
C  BMYTPD   R    by model year TPD, calculated in CEVTPD 
C  BMTPDF   R    by model year TPD fractions, calculated in CEVTPD 
C 
C  /CUMCOM/: 
C 
C  CUMMIL   R    cumulative mileage accumulation rates 
C 
C  /FLAGS1/: 
C 
C  USR94P   I    decides Control Section flags read format & prompt status: 
C                1 = use Mobile default data for implementation schedules 
C                    and EFs for 94+ LDG vehicles 
C                2 = user has supplied data for implementation schedules or 
C                    EFs for 94+ LDG vehicles. 
C  TAMFLG   I    1 = use default tampering rates 
C                2 = user supplies tampering rates for 8 categories 
C  VMFLAG   I    1 = use default vmt mix 
C                2 = user supplies a new vmt mix for each scenario 
C                3 = user supplies 1 vmt mix to apply to all scenarios 
C  OXYFLG   I    Oxy fuel adjustments flag: 
C                1=no adjustments 
C                2=apply adjustments. 
C  SPCOFLG  I    1 = use default speed correction coefficient 
C                2 = user supplies speed correction coefficient 
C  SCFFLG   I    1 = use default SCF 
C                2 = user supplies SCF 
C  SCFUSER  I    1 = use default SCF 
C                2 = user supplies SCF 
C  /FLAGS2/: 
C 
C  MYMRFG   I    1 = default registration distributions/mileage accrual rates 
C                2 = user supplies mileage accumulation rates 
C                3 = user supplies registration distributions 
C                4 = user supplies both 
C  NEWFLG   I    1 = use all default exhaust emission rates 
C                2 = replace up to 99 of the default rates; caution - 
C                    (a) for IV = 4 or 7 (HDV), use gpbhp units instead of gpm 
C                    (b) each IV x IP x IR case has a limit on new efs - see UG 
C  IMFLAG   I    1 = no I/M program 
C                2 = specify and apply I/M program 
C                3 = specify and apply two I/M programs 
C                4 = specify and apply a single I/M program (same as IMFLAG=2) 
C                    and beginning in 1994, emission credits based on assuming 
C                    a specific I/M program for Tier 1 vehicles override 
C                    the user-specified program for LDGVs & LDGTs               
C                5 = specify and apply twoe I/M program (same as IMFLAG=3) 
C                    and beginning in 1994, emission credits based on assuming 
C                    a specific I/M program for Tier 1 vehicles override 
C                    the user-specified program for LDGVs & LDGTs               
C                6 = indicate the fact that the information about I/M program 
C                    will be entered in the special record just prior to first 
C                    I/M descriptive record in the One-Time Data Section.   
C                    This information include the following flags:  
C                     NUMIMP-the number of the I/M descriptive records;  
C                     TI1FLG-Yes or No  for "Tier 1 " program; 
C                     ITTCFG Yes or No for "TTC" program, 
C                     RSDFLG- Yes or No for "RSD" program. 
C 
C                     If ITTCFG=3, another record should be in the 
C                     One-Time Data Section with the following fields:  
C                     TTCEFF- percent of effectiveness for TTC program. 
C  EVIMFLAG I    1 = no evap I/M program 
C                2+= specify and apply one or more evap I/M program 
C   
C 
C  /FLAGS3/: 
C 
C  ATPFLG   I    1 = no anti-tampering (A-T) program 
C                2 = specify A-T program to apply to all scenarios 
C                3 = Enter Pressure Check record 
C                4 = Enter Purge Check record 
C                5 = Enter ATP record & Pressure Check record 
C                6 = Enter ATP record & Purge Check Record 
C                7 = Enter Pressure Check record and Purge Check record 
C                8 = Enter ATP, Pressure & Purge Check records. 
C 
C  CAAFLG   I    0 = DEFAULT Clean Air ACT is in effect 
C                1 = Clean Air ACT is NOT in effect 
C 
C  POST96IM I    0 = DEFAULT Exhaust I/M programs not allowed for 96+ model years. 
C                1 = Exhaust I/M programs allowed for 96+ model years. 
C                *** POST96IM *** is now always set to 1. Deactivated. 
C 
C  CAAMY    I    Year that the Clean Air ACT is disabled - DEFAULT 1993 
C 
C  RLFLAG   I    1 = use the uncontrolled rl rate (assume no recovery 
C                    systems installed) 
C 
C                2 = compute and assign Stage II (S2) only vrs controlled 
C 
C                3 = weight uncontrolled & controlled rate by best estimate 
C                    of ob vrs tamper rate & its complement, respectively, 
C                    to get ob only vrs controlled refueling losses. 
C                    Assume 100% installation. 
C                        (If my/iv not in OB program, same as 1) 
C 
C                4 = compute and assign S2 or ob vrs controlled losses. 
C                        (If my/iv not in OB program, same as 2) 
C                        (If my/iv is in OB program, same as 3) 
C 
C                5 = assign 0.0 (no losses case: results = release MOBILE3). 
C 
C                Note: Due to the Onboard Refueling Vapor Recovery Regulation 
C                      of 1994, an incoming value of 1 is converted to 3 
C                      and an incoming value of 2 is converted to 4. This 
C                      results in there always being an OB program in force 
C                      except when RLFLAG = 5 or NEWFLG = 5 (CAA effects 
C                      disabled). When RLFLAG = 5, there are no refueling 
C                      losses. When NEWFLG = 5, there are uncontrolled losses 
C                      (RLFLAG is converted back to 1). 
C 
C  /FLAGS4/: 
C 
C  NMHFLG   I    1 = print total HC emission factors 
C                2 = print non-methane HC emission factors 
C                3 = print VOC emission factors 
C                4 = print TOG emission factors 
C                5 = print NMOG emission factors 
C  HCFLAG   I    1 = do not print HC components 
C                2 = print HC components (exhaust, evaporative, 
C                    refueling loss) 
C                3 = print HC components in grams/mile and alternate units 
C                4 = print HC components (exhaust, evaporative, 
C                    refueling loss) but refueling loss is not included in 
C                    total evaporative(HC) and the evap + exhaust 'composite' sums 
C  COLDFG   I    1 = no cold CO standard 
C                2 = cold CO standard CCSTD applies starting in MYCOL; 
C                    uses "lenient" standard algorithm: all cold CO standard 
C                    reductions occur in Bag 1. 
C                3 = cold CO standard CCSTD applies starting in MYCOL; 
C                    uses "stringent" standard algorithm: cold CO standard 
C                    causes proportional reductions in all Bags. 
C 
C  /MYRCAL/: 
C 
C  XMYM     R    January 1 model year mileage accumulation rates used to 
C                compute CUMMIL 
C  MEVMYR   R    January 1 registration distributions 
C  TF       R    travel weighting fractions 
C  TFMYM    R    January 1 mileage accumulation rates used to compute TF 
C 
C  /REGISF/: 
C 
C  GSFVCT   R    fleet gasoline/diesel sales fractions used to split LDV, LDT 
C  USRGSF   R    user-entered gas/diesel sales fractions for LDV & LDT 
C 
C  /RESUL1/: 
C 
C  EFFTP         R    composite (3 bag cycle weighted) FTP emission factor for each 
C                    IV9 class = exhaust (+ evap/cc + refuel loss) 
C  EFEXH         R    exhaust component of IV9 HC ef 
C  EFEVAP        R    evaporative + crankcase component of IV9 HC ef 
C  REFL_GM_MI    R    refueling loss component of IV9 HC ef 
C  EFRUNL        R    running loss component of IV9 HC ef 
C  G_MI_EVP      R    evaporative HC emissions of IV9 in grams/mile 
C 
C  /RESUL2/: 
C 
C  RLGGAL   R    refueling loss component of IV9 HC ef (in grams per gallon) 
C  RSTGPH   R    resting loss component of IV9 HC ef (in grams per hour) 
C  EFIDLE   R    hot stabilized idle emission factor for each IV9 class 
C 
C  /RESUL3/: 
C 
C  Prefix in v/a in /RESUL3/: 
C 
C  'V'  =  vmtmix weighted (ALLVEH IP, each IV's ef weighted by its VMTMIX) 
C 
C  VFTP     R    For each IP, FTP emission factor for each IV class weighted by 
C                corresponding normalized vmt mix fraction and summed together 
C  VEXH     R    vmt mix weighted local (city) exhaust HC 
C  VEVAP    R    vmt mix weighted local (city) evaporative HC 
C  VLOSS    R    vmt mix weighted local (city) refueling loss HC 
C  VRUNLS   R    vmt mix weighted local (city) running loss HC 
C  VRSTLS   R    vmt mix weighted local (city) resting loss HC 
C  VIDLE    R    For each IP, idle emission factor for each IV class weighted by 
C                corresponding normalized vmt mix fraction and summed together 
C 
C  /RESUL4/: 
C 
C  RUNEF    R    Separated running emission factors 
C  STARTEF  R    Separated start emission factors 
C  TOTEVAP  R    Sum of evaporative emission factors 
C 
C 
C 
C  Notes: 
C 
C  BD 27 was modified for MOBILE4 In-House Version 02 to change /RESUL1/'s 
C  EFLOSS to REFL_GM_MI (consistent w/ EFRUNL for running loss). 
C 
C  COLDFG is set to 2 (additive COLD CO) for MOBILE4.1 v D01. 
C  Added G_MI_EVP, RSTGPH and VRSTLS (resting losses) for MOBILE4.1 D03. 
C  Added comments for flags and expanded OFFMTH array for MOBILE5v13. 
C 
C  5-May-1994 @ CSC-tmm request 2-411 
C       Updated comments for RLFLAG to reflect the result of this request. 
C  11-August-1994 @ CSC-pme request 2-455. BD27 was modified for  
C  MOBILE5b to add comments for Prompt flag possible values. 
C  24-August-1994 @CSC-pme request 2-446. Third dimension of OFFMTH was removed 
C  and initialization instruction was updated. A comment to explain variable PPTFLG 
C  and its initial value were added. 
C 
C  3 Dec 96 DynTel@MLA 2-622  Expanded the arrays to accomodate the new  
C           Heavy Duty Gasoline, Heavy Duty Diesel, and LDGT vehicle  
C           categories. 
C  December 16 1997 Dyntel@MMS 2-645 modified multispeed to a single speed.  
C                   SPDFLG renamed RNLFLG.RNLFLG is equivalent of  SPDFLG  
C                   by 1 for : RNLFLG from 2 to 4. 
C  3 Mar 98 DynTel@RJD 2-204  Changed the name of the array storing the Model 
C           Year Registration fractions from JANMYR to MEVMYR to reflect the 
C           fact that January is no longer the only month of evaluation that 
C           the Mobile model can handle. 
C 10 Mar 98 DynTel@RJD 2-622 Added integer labels for emission factors to 
C           facilitate the efficient calculation of emission factors for 
C           composite vehicles in EFCVT. 
C 21 Apr 98 DynTel@RJD 2-654 Added an extra dimension to the OFFMTH array 
C           to distinguish between start and running emissions.       
C 4  Jun 98 @DynTel - mla 2-659 Removed LOCFLG since the LAP record is  
C           obsolete in Mobile6 
C 18 Jun 98 @DynTel - bag 2-000 Removed OUTFMT because it has been  
C           superceded by spreadsheet flags 
C  1 Jul 98 @DynTel - MLA 2-000 Removed PROJID and its include file PROJEC. 
C           (PROJID is obsolete in Mobile6.) 
C 
C 24-Aug-98 @Dyntel-nh Request 0-632 Removed dead code related to special I/M 
C    
      IMPLICIT NONE 
      INCLUDE 'ALUHIN.I' 
      INCLUDE 'ALUOUT.I' 
      INCLUDE 'BASEQ3.I' 
      INCLUDE 'CEVBMY.I' 
      INCLUDE 'CUMCOM.I' 
      INCLUDE 'DAILY_EVAP.I' 
      INCLUDE 'FLAGS1.I' 
      INCLUDE 'FLAGS2.I' 
      INCLUDE 'FLAGS3.I' 
      INCLUDE 'FLAGS4.I' 
      INCLUDE 'MYRCAL.I' 
      INCLUDE 'OFFSET.I' 
      INCLUDE 'PARAME.I' 
      INCLUDE 'REGISF.I' 
      INCLUDE 'RESUL1.I' 
      INCLUDE 'RESUL2.I' 
      INCLUDE 'RESUL3.I' 
      INCLUDE 'RESUL4.I' 
C 
      DATA AC, ABSHUM/1.0, 75.0/ 
      DATA GLAC, GLAHUM/ 1.0, 75.0/ 
      DATA GSUNRISE, GSUNSET, GPEAKSUN/ 1, 16, 5, 11/ 
      DATA GCLD_CVR/0.0/ 
      DATA ALHRET/10800*0.0/ 
      DATA NEWVEH,NEWMYF,NEWMYL,ZMLNEW,DRNEW/3*0,2*0.0/,KINKS/1/ 
      DATA BMYMPD,BMTPDF/1400*0.0/ 
      DATA CUMMIL/700*0.0/ 
      DATA TAMFLG,VMFLAG,SPCOFLG,SCFFLG,OXYFLG,DSFLAG,CONVFLG, 
     *     USR94P,SCFUSER/7*0,2*1/ 
      DATA MYMRFG,NEWFLG,IMFLAG,EVIMFLAG/4*0/ 
      DATA ATPFLG,RLFLAG,CAAFLG/3*0/ 
      DATA POST96IM/1/ 
      DATA CAAMY/1993/ 
      DATA NMHFLG,HCFLAG,COLDFG/0,1,2/ 
      DATA XMYM,MEVMYR,TF,TFMYM/2800*0.0/ 
      DATA OFFSET,OFFMTH/0.0,42000*0.0/ 
      DATA MOFFMTH/42000*0.0/ 
      DATA GSFVCT/28*0.0/,USRGSF/350*0.0/ 
      DATA EFFTP,EFEXH,EFEVAP,REFL_GM_MI,EFRUNL,G_MI_EVP/377*0.0/ 
      DATA IEFFTP,IEFEXH,IEFEVAP,IEFREFL,IEFRUNL /1,2,3,4,5/ 
      DATA IRLGGAL,IRSTGPH,IEFIDLE,IGMIEVP,IEFPM /6,7,8,9,10/ 
      DATA VFTP,VEXH,VEVAP,VLOSS,VRUNLS,VRSTLS,VIDLE/11*0.0/ 
      DATA EVAP_GM_MI,EVAP_GM_HR,EVAP_GM_TRP,EVAP_GM_GA/14*0.0/ 
      DATA RUNEF, STARTEF, TOTEVAP / 203*0.0 / 
      DATA IRUNEF, ISTARTEF, ITEVAP /22,23,24 / 
C 
      END 
