C Version @(#)$Id$ $Source$ $Date$ 

      BLOCK DATA BD27
C
C  BLOCK DATA Subprogram 27: /ALUHIN/,/ALUOUT/,/BASEQ3/,/CEVCOM/,/CEVBMY/,
C                            /CUMCOM/,/FLAGS1/,/FLAGS2/,/FLAGS3/,/FLAGS4/,
C                            /IMPAR1/,/IMPAR2/,/IM240P/,/MYRCAL/,/OFFSET/,
C                            /PROJEC/,/REGISF/,/RESUL1/,/RESUL2/,/RESUL3/
C
C  BLOCK DATA 27's common blocks did not need to be initialized, so were
C  assigned dummy zeroes, to facilitate debugging, to prevent control statement
C  failure when compiled under a FORTRAN that does not employ the left to right
C  rule, etc.
C
C  Common block array subscripts:
C
C  ALHRET(25,4,3)  -  ALHRET ( JDX, IG, IP )
C  BMYMPD(25,8)    -  BMYMPD ( IG, IV )
C  BMYTPD(25,8)    -  BMYTPD ( IG, IV )
C  BMTPDF(25,8)    -  BMTPDF ( IG, IV )
C  CUMMIL(25,8)    -  CUMMIL ( JDX, IV )
C  DSIZE(25,4)     -  DSIZE  ( IDX, IVGAS )
C  EFEVAP(9)       -  EFEVAP ( IV9 )
C  EFEXH(9)        -  EFEXH  ( IV9 )
C  EFFTP(3,9)      -  EFFTP  ( IP, IV9 )
C  EFIDLE(3,9)     -  EFIDLE ( IP, IV9 )
C  EFREFL(9)       -  EFREFL ( IV9 )
C  EFRUNL(9)       -  EFRUNL ( IV9 )
C  EFRSTL(9)       -  EFRSTL ( IV9 )
C  GSFVCT(8)       -  GSFVCT ( IV )
C  ICYIM(2)        -  ICYIM  ( I/M Program )
C  ILDT(4,2)       -  ILDT   ( IVIM, IM PROGRAM )
C  ISTRIN(2)       -  ISTRIN ( I/M Program )
C  ITEST(2)        -  ITEST  ( I/M Program )
C  JANMYR(25,8)    -  JANMYR ( JDX, IV )
C  MODYR(2,2)      -  MODYR  ( 1st/ Last, I/M Program )
C  OFFCO(25,3)     -  OFFCO  ( JDX, IVP )
C  OFFMTH(25,8,2)  -  OFFMTH ( JDX, IV, ICAP )
C  PROJID(20)      -  PROJID ( ICH )
C  RLGGAL(9)       -  RLGGAL ( IV9 )
C  RSTGPH(9)       -  RSTGPH ( IV9 )
C  TF(25,8)        -  TF     ( IDX, IV )
C  TFMYM(25,8)     -  TFMYM  ( JDX, IV )
C  TRAILR(3)       -  TRAILR ( IP )
C  USRGSF(25,2)    -  USRGSF ( IDX, IVTEMP )
C  USRMPD(8)       -  USRMPD ( IV )
C  USRTPD(8)       -  USRTPD ( IV )
C  VFTP(3)         -  VFTP   ( IP )
C  VIDLE(3)        -  VIDLE  ( IP )
C  WAIVER(2,2)     -  WAIVER ( IGW, I/M Program )
C  XLOAD(3)        -  XLOAD  ( IP )
C  XMYM(25,8)      -  XMYM   ( JDX, IV )
C
C  Common block dictionary:
C
C   Name   Type              Description
C  ------  ----  --------------------------------------------------
C  /ALUHIN/:
C
C  AC       R    a/c usage factor
C  XLOAD    R    extra load usage fractions
C  TRAILR   R    trailer towing usage fractions
C  ABSHUM   R    absolute humidity level
C  DB       R    dry bulb temperature
C  WB       R    wet bulb temperature
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
C  /CEVCOM/:
C
C  USRTPD   R    user supplied average trips per day
C  USRMPD   R    user supplied average miles per day
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
C  PROMPT   I    decides Control Section flags read format & prompt status:
C                1 = use vertical flag format; do not prompt
C                2 = use vertical flag format; prompt for input
C                3 = use horizontal flag format; do not prompt
C                4 = use horizontal flag format; prompt for input
C                Internally, VHMODE is assigned the flag format function.
C                PROMPT is then reduced to the prompt function:
C                1 = do not prompt for input
C                2 = prompt for input
C  TAMFLG   I    1 = use default tampering rates
C                2 = user supplies tampering rates for 8 categories
C  SPDFLG   I    1 = user supplies one speed for all vehicle types
C                2 = user supplies 8 speeds, 1 for each vehicle type
C  VMFLAG   I    1 = use default vmt mix
C                2 = user supplies a new vmt mix for each scenario
C                3 = user supplies 1 vmt mix to apply to all scenarios
C  OXYFLG   I    Oxy fuel adjustments flag:
C                1=no adjustments
C                2=apply adjustments.
C
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
C  ALHFLG   I    1 = no special correction factors
C                2 = user supplies 6 values: 1 a/c use factor, 3 extra load
C                    usage fractions, 1 trailer towing fraction & 1 humidity
C                    level
C                3 = user supplies 10 values: the 6 above, except 3 instead of
C                    1 towing fractions plus 1 db & 1 wb temperatures
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
C  TPDFLG   I    1 = use default calculated trips per day rates
C                2 = user supplied data from which trips per day are calculated
C  RLFLAG   I    1 = refueling loss is always 0.0, same as release MOBILE3
C  TEMFLG   I    1 = use default calculated hot soak, exhaust,
C                   and running loss temperatures
C                2 = use ambient temperature for hot soak, exhaust,
C                   and running loss temperatures
C  LOCFLG   I    1 = user enters a local area parameters record after the
C                    the required parameters record of each scenario
C                2 = user enters 1 local area parameters record for all
C                    scenarios, placing it as the last record in one-time data.
C                2 = if MY < 1989, then use 100% uncontrolled rate,
C                    else (MY >= 1989) use 100% controlled rate (no tampering)
C                3 = if MY < 1989, then use 100% uncontrolled rate,
C                    else (MY >= 1989) use best estimate of tampering rate to
C                    weight uncontrolled & controlled (VRS installed) rates
C                4 = use 100% uncontrolled refueling evap HC loss (no VRS's)
C  OUTFMT    I    various output format selections:
C                1 = new for MOBILE4 long (212 column) numerical format
C                2 = new for MOBILE4 short (130 column) numerical format
C                3 = 112 column descriptive format
C                4 = 80 column descriptive format
C                5 = BYMY format
C                6 = new for MOBILE5v8 spreadsheet version format
C
C  /FLAGS4/:
C
C  PRTFLG   I    1 = print only the HC factors
C                2 = print only the CO factors
C                3 = print only the NOX factors
C                4 = print the factors for all three: HC, CO and NOX
C  IDLFLG   I    1 = do not print idle emission factors
C                2 = print idle emission factors
C  NMHFLG   I    1 = print total HC emission factors
C                2 = print non-methane HC emission factors
C                3 = print VOC emission factors
C                4 = print TOG emission factors
C                5 = print NMOG emission factors
C  HCFLAG   I    1 = do not print HC components
C                2 = print HC components (exhaust, evaporative,
C                    refueling loss)
C                3 = print HC components and evaporative break down
C                4 = HCFLAG=3 with refueling zeroed out
C                5 = HCFLAG=3 with emission factors in grams
C  COLDFG   I    1 = no cold CO standard
C                2 = cold CO standard CCSTD applies starting in MYCOL;
C                    uses "lenient" standard algorithm: all cold CO standard
C                    reductions occur in Bag 1.
C                3 = cold CO standard CCSTD applies starting in MYCOL;
C                    uses "stringent" standard algorithm: cold CO standard
C                    causes proportional reductions in all Bags.
C
C  /IMPAR1/:
C
C  ICYIM    I    I/M program start year
C  ISTRIN   I    Technology I/II stringency level of I/M program
C  MODYR    I    first (earliest) and last (latest) year included in the I/M program
C  WAIVER   R    waiver rate - percent of reviewed fleet that is exempted from
C                I/M testing
C  CRIM     R    enforcement level - percent of fleet reviewed for I/M testing
C                WAIVER is based on CRIM: the percent of the total fleet
C                tested = (1.0 - WAIVER(IGW)) * (CRIM * fleetsize),
C                IGW = 1 for pre-1981 and IGW = 2 for 1981+ MY vehicles.
C
C  /IMPAR2/:
C
C  ILDT     I    inclusion vector - vehicle classes to be covered by I/M:
C                    1 = omit    2 = included
C  ITEST    I    technology IV+ test type:
C                1 = idle test
C                2 = 2500 rpm idle test
C                3 = loaded idle test
C                4 = IM240 test
C
C  /IM240P/:
C
C  DSIZE    R    non-ATP tampering rates captured in EVPEME
C
C  /MYRCAL/:
C
C  XMYM     R    January 1 model year mileage accumulation rates used to
C                compute CUMMIL
C  JANMYR   R    January 1 registration distributions
C  TF       R    travel weighting fractions
C  TFMYM    R    January 1 mileage accumulation rates used to compute TF
C
C  /OFFSET/:
C
C  OFFCO    R    CO offset to basic emission rate
C  OFFMTH   R    methane offset to basic emission rate
C
C  /PROJEC/:
C
C  PROJID   R    user supplied project (i.e., run) title
C
C  /REGISF/:
C
C  GSFVCT   R    fleet gasoline/diesel sales fractions used to split LDV, LDT
C  USRGSF   R    user-entered gas/diesel sales fractions for LDV & LDT
C
C  /RESUL1/:
C
C  EFFTP    R    composite (3 bag cycle weighted) FTP emission factor for each
C                IV9 class = exhaust (+ evap/cc + refuel loss)
C  EFEXH    R    exhaust component of IV9 HC ef
C  EFEVAP   R    evaporative + crankcase component of IV9 HC ef
C  EFREFL   R    refueling loss component of IV9 HC ef
C  EFRUNL   R    running loss component of IV9 HC ef
C  EFRSTL   R    resting loss component of IV9 HC ef
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
C
C  Notes:
C
C  BD 27 was modified for MOBILE4 In-House Version 02 to change /RESUL1/'s
C  EFLOSS to EFREFL (consistent w/ EFRUNL for running loss).
C
C  COLDFG is set to 2 (additive COLD CO) for MOBILE4.1 v D01.
C  Added EFRSTL, RSTGPH and VRSTLS (resting losses) for MOBILE4.1 D03.
C  Added comments for flags and expaned OFFMTH array for MOBILE5v13.
C
C
      CHARACTER*4 PROJID
      INTEGER PROMPT,TAMFLG,SPDFLG,VMFLAG,OXYFLG,DSFLAG
      INTEGER ALHFLG
      INTEGER ATPFLG,TPDFLG,RLFLAG,TEMFLG,OUTFMT
      INTEGER PRTFLG,HCFLAG,COLDFG
      REAL JANMYR
C
      COMMON /ALUHIN/ AC,XLOAD(3),TRAILR(3),ABSHUM,DB,WB
      COMMON /ALUOUT/ ALHRET(25,4,3)
      COMMON /BASEQ3/ NEWVEH,NEWMYF,NEWMYL,ZMLNEW,DRNEW,KINKS
      COMMON /CEVCOM/ USRTPD(8),USRMPD(8)
      COMMON /CEVBMY/ BMYMPD(25,8),BMYTPD(25,8),BMTPDF(25,8)
      COMMON /CUMCOM/ CUMMIL(25,8)
      COMMON /FLAGS1/ PROMPT,TAMFLG,SPDFLG,VMFLAG,OXYFLG,DSFLAG
      COMMON /FLAGS2/ MYMRFG,NEWFLG,IMFLAG,ALHFLG,IMTIER
      COMMON /FLAGS3/ ATPFLG,TPDFLG,RLFLAG,LOCFLG,TEMFLG,OUTFMT
      COMMON /FLAGS4/ PRTFLG,IDLFLG,NMHFLG,HCFLAG,COLDFG
      COMMON /IMPAR1/ ICYIM(2),ISTRIN(2),MODYR(2,2),WAIVER(2,2),CRIM(2)
      COMMON /IMPAR2/ ILDT(4,2),ITEST(2)
      COMMON /IM240P/ DSIZE(25,4)
      COMMON /MYRCAL/ XMYM(25,8),JANMYR(25,8),TF(25,8),TFMYM(25,8)
      COMMON /OFFSET/ OFFCO(25,3),OFFMTH(25,8,2)
      COMMON /PROJEC/ PROJID(20)
      COMMON /REGISF/ GSFVCT(8),USRGSF(25,2)
      COMMON /RESUL1/ EFFTP(3,9),EFEXH(9),EFEVAP(9),EFREFL(9),EFRUNL(9)
      COMMON /RESUL2/ RLGGAL(9),RSTGPH(9),EFRSTL(9),EFIDLE(3,9)
      COMMON /RESUL3/ VFTP(3),VEXH,VEVAP,VLOSS,VRUNLS,VRSTLS,VIDLE(3)
C
      DATA AC,XLOAD,TRAILR,ABSHUM,DB,WB/10*0.0/
      DATA ALHRET/300*0.0/
      DATA NEWVEH,NEWMYF,NEWMYL,ZMLNEW,DRNEW/3*0,2*0.0/,KINKS/1/
      DATA USRTPD,USRMPD/16*0.0/
      DATA BMYMPD,BMYTPD,BMTPDF/600*0.0/
      DATA CUMMIL/200*0.0/
      DATA PROMPT,TAMFLG,SPDFLG,VMFLAG,OXYFLG,DSFLAG/6*0/
      DATA MYMRFG,NEWFLG,IMFLAG,ALHFLG,IMTIER/5*0/
      DATA ATPFLG,TPDFLG,RLFLAG,LOCFLG,TEMFLG,OUTFMT/0,3,4*0/
      DATA PRTFLG,IDLFLG,NMHFLG,HCFLAG,COLDFG/4*0,2/
      DATA ICYIM,MODYR,ISTRIN,WAIVER,CRIM/6*2020,0,0,6*0.0/
      DATA ILDT,ITEST/10*0/
      DATA DSIZE/100*0.0/
      DATA XMYM,JANMYR,TF,TFMYM/800*0.0/
      DATA OFFCO,OFFMTH/475*0.0/
      DATA PROJID/20*'    '/
      DATA GSFVCT,USRGSF/58*0.0/
      DATA EFFTP,EFEXH,EFEVAP,EFREFL,EFRUNL,EFRSTL/72*0.0/
      DATA RLGGAL,RSTGPH,EFIDLE/45*0.0/
      DATA VFTP,VEXH,VEVAP,VLOSS,VRUNLS,VRSTLS,VIDLE/11*0.0/
C
      END
