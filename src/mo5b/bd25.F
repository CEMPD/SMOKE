C Version "@(#)$Id$ $Source$ $Date$ 

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
C  Common block array subscripts:
C
C  AER(11,11,2,2)    -  AER    ( IC11, IC11, IPG, IAY )
C  AIR(2)            -  AIR    (  IAY )
C  ATR(5,2)          -  ATR    ( IDBAG , IAY )
C  BTR(9,2)          -  BTR    ( ID, IAY )
C  CAN(2)            -  CAN    (  IAY )
C  CAP(2)            -  CAP    (  IAY )
C  CAT(2)            -  CAT    (  IAY )
C  CCEMI(7,4)        -  CCEMI  ( IG4, IVTAM )
C  CSAE(12)          -  CSAE   ( IC )
C  CSIZE(12,2)       -  CSIZE  ( IC , IAY )
C  DISTYP(8)         -  DISTYP  ( IDUSER )  (see ID defn in (1,999) )
C  EFF(2,3,2,15)     -  EFF    (  ATPPGM, ICOL , IATP , IROW )
C  EFFSAV(2,3,2,15)  -  EFFSAV (  ATPPGM, ICOL , IATP , IROW )
C  EGR(2)            -  EGR    (  IAY )
C  EGS(7,2)          -  EGS    ( IEG, IPG )
C  EMEGR(4,3,3)      -  EMEGR  ( IG2, IB, IP )
C  EMI(7,3,3)        -  EMI    ( IEG, IB, IP )
C  EVPAER(2,2)       -  EVPAER ( IHG, IAY )
C  EVPTGS(3,2,4)     -  EVPTGS ( IG3, IHG, IVTAM )
C  GCONLY(25,4)      -  TOB    ( IDX, IVTAM )
C  INDXCS(5,4)       -  INDXCS ( NT, IS )
C  LVTFLG(4)         -  LVTFLG ( IVTAM )
C  MYCCEI(7,4)       -  MYCCEI ( IG4, IVTAM )
C  MYEGR(4,4)        -  MYEGR  ( IG2, IVTAM )
C  MYGTAM(4)         -  MYGTAM ( IVTAM )
C  MYTGS(11,6,4)     -  MYTGS  ( IG1, IQG, IVTAM )
C  MYVTS(3,2,4)      -  MYVTS  ( IG3, IHG, IVTAM )
C  ECK(2)            -  ECK    (  NCKIAY )
C  NTERMS(4)         -  NTERMS ( IS )
C  OMTCF(25,3,8)     -  OMTCF  ( IDX, IP, IV )
C  OMTTAM(25,3,4)    -  OMTTAM ( IDX, IP, IVTAM )
C  PCV(2)            -  PCV    (  IAY )
C  PBBTR(3,2,25,4)   -  PBBTR  ( IDPB, IAY, IDX, IVTAM )
C  TAMBAG(3,3,25,4)  -  TAMBAG ( IB, IP, IDX, IVTAM )
C  TAMA50(9,4,3,2)   -  TAMA50 ( ID, IVTAM, IGD, IM )
C  TAMDR(9,4,3,2)    -  TAMDR  ( ID, IVTAM, IGD, IM )
C  TAMZML(9,4,3,2)   -  TAMZML ( ID, IVTAM, IGD, IM )
C  TCC(25,4)         -  TCC    ( IDX, IVTAM )
C  TDU(25,4)         -  TDU    ( IDX, IVTAM )
C  TGS(11,6,4)       -  TGS    ( IG1, IQG, IVTAM )
C  THS(2,25,4)       -  THS    ( IFDS, IDX, IVTAM )
C  TNK(2)            -  TNK    ( IAY )
C  TOB(25,4)         -  TOB    ( IDX, IVTAM )
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
C                for each of the 4 (LDGV, LDGT1, LDGT2, HDGV) flags:
C                  1 = not covered
C                  2 = covered
C                all "1"s = same result as ATPFLG = 1 = no program - the ATP
C                effectiveness matrices are not applied.
C
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
C  OMTCF    R    operating mode & temperature correction factor to the
C                composite emission factor.  Calculated in BIGCFX.  Passed to
C                function BEF to be multiplied times the CEF.
C  OMTTAM   R    operating mode and temperature corrected tampering effect on
C                bag emissions.  Calculated in BIGCFX.  Passed to BEF to be
C                added to the (corrected for OMTCF & I/M) CEF.
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
C  MYTGSH(11,6)    -  MYTGSH ( IG1, IQG )
C  MYTGSL(11,6,3)  -  MYTGSL ( IG1, IQG, IVTAM )
C  TA5IM(9,4,3)    -  TA5IM  ( ID, IVTAM, IGD )
C  TA5NON(9,4,3)   -  TA5NON ( ID, IVTAM, IGD )
C  TDRIM(9,4,3)    -  TDRIM  ( ID, IVTAM, IGD )
C  TDRNON(9,4,3)   -  TDRNON ( ID, IVTAM, IGD )
C  TGSHDV(11,6)    -  TGSHDV ( IG1, IQG )
C  TGSLDV(11,6,3)  -  TGSLDV ( IG1, IQG, IVTAM )
C  TZLIM(9,4,3)    -  TZLIM  ( ID, IVTAM, IGD )
C  TZLNON(9,4,3)   -  TZLNON ( ID, IVTAM, IGD )
C
C  Local array dictionary:
C
C   Name   Type              Description
C  ------  ----  -------------------------------------------------------
C  MYTGSH   I    my range ub years for HDGV for ITAMPT indexing TGS
C  MYTGSL   I    my range ub years for LDGV/T1/T2 for ITAMPT indexing TGS
C  TA5IM    R    TAMA50s for the I/M case
C  TA5NON   R    TAMA50s for the Non-I/M case
C  TDRIM    R    TAMDRs for the I/M case
C  TDRNON   R    TAMDRs for the Non-I/M case
C  TGSHDV   R    technology group sizes for HDGV
C  TGSLDV   R    technology group sizes for LDGV/T1/T2
C  TZLIM    R    TAMZMLs for the I/M case
C  TZLNON   R    TAMZMLs for the Non-I/M case
C
C  Notes:
C
C  The following common blocks defined and used elsewhere in MOBILE4 contain
C  data used by the tampering effects section:
C
C  /CUMCOM/: CUMMIL - cumulative mileage array (JDX,IV)
C  /FLAGS2/: IMFLAG  - I/M program flag
C  /IMPAR1/: ICYIM  - I/M program(s) start year
C            MODYR -  first and last model year covered by I/M program(s)
C            ILDT   - vehicle types covered by I/M program(s)
C
C  See the block data sections for these commons for more information on their
C  variables / arrays.
C
C  FORTRAN limits the number of continuation lines to 19.  Use equivalenced
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
C  23 Feb 96 @DynTel-MLA 2-611 Added array EFFSAV to save default values
C    for EFF.
C
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
      DIMENSION MYTGSL(11,6,3),MYTGSH(11,6)
      DIMENSION TGSLDV(11,6,3),TGSHDV(11,6)
      DIMENSION TZLNON(9,4,3),TZLIM(9,4,3),TDRNON(9,4,3),TDRIM(9,4,3)
      DIMENSION TA5NON(9,4,3),TA5IM(9,4,3)
C
      EQUIVALENCE (MYTGS(1,1,1),MYTGSL(1,1,1))
      EQUIVALENCE (MYTGS(1,1,4),MYTGSH(1,1))
      EQUIVALENCE (TGS(1,1,1),TGSLDV(1,1,1))
      EQUIVALENCE (TAMA50(1,1,1,1),TA5NON(1,1,1))
      EQUIVALENCE (TAMA50(1,1,1,2),TA5IM(1,1,1))
      EQUIVALENCE (TAMDR(1,1,1,1),TDRNON(1,1,1))
      EQUIVALENCE (TAMDR(1,1,1,2),TDRIM(1,1,1))
      EQUIVALENCE (TAMZML(1,1,1,1),TZLNON(1,1,1))
      EQUIVALENCE (TAMZML(1,1,1,2),TZLIM(1,1,1))
      EQUIVALENCE (TGS(1,1,4),TGSHDV(1,1))
C
C  /ATPAR1/: values can be replaced by user entered values.
C
      DATA LAPSY,LAP1ST,LAPLST/3*NEYR/,LVTFLG/4*1/
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
C    3 = catalyst only
C    4 = EGR only
C    5 = EGR / 3-way catalyst
C    6 = no EGR / 3-way catalyst
C
      DATA TGSLDV/
C  LDGV
     1  .00, .05, .10, .30, .15, .10, .05, .00, .00, .00, .00,
     2  .00, .30, .20, .25, .65, .75, .70, .60, .40, .30, .00,
     3  .00, .50, .55, .65, .30, .25, .30, .40, .60, .70, .00,
     4  .00, .80, .90, .05, .00, .00, .00, .00, .00, .00, .00,
     5  .00, .07, .85, .93, .90, .00, .00, .00, .00, .00, .00,
     6  .00, .01, .03, .07, .10, .00, .00, .00, .00, .00, .00,
C  LDGT1
     1  .00, .05, .10, .30, .10, .00, .00, .00, .00, .00, .00,
     2  .00, .30, .20, .40, .50, .60, .75, .55, .50, .00, .00,
     3  .00, .40, .50, .55, .40, .50, .40, .25, .45, .50, .00,
     4  .00, .80, .90,1.00, .95, .90, .80, .70, .60, .50, .05,
     5  .00, .05, .10, .20, .30, .40, .50, .95, .00, .00, .00,
     6  .00, .00, .00, .00, .00, .00, .00, .00, .00, .00, .00,
C  LDGT2
     1  .00, .00, .00, .00, .00, .00, .00, .00, .00, .00, .00,
     2  .00, .50, .60, .75, .55, .50, .00, .00, .00, .00, .00,
     3  .00, .50, .40, .25, .45, .50, .00, .00, .00, .00, .00,
     4  .00, .30,1.00, .90, .70, .60, .50, .15, .00, .00, .00,
     5  .00, .10, .30, .40, .50, .85, .00, .00, .00, .00, .00,
     6  .00, .00, .00, .00, .00, .00, .00, .00, .00, .00, .00/
C  HDGV
      DATA TGSHDV/
     1  .00, .07, .00, .00, .00, .00, .00, .00, .00, .00, .00,
     2  .00, .30, .00, .00, .00, .00, .00, .00, .00, .00, .00,
     3  .00, .05, .00, .10, .25, .00, .00, .00, .00, .00, .00,
     4  .00, .95,1.00, .85, .75, .00, .00, .00, .00, .00, .00,
     5  .00, .15, .25, .00, .00, .00, .00, .00, .00, .00, .00,
     6  .00, .00, .00, .00, .00, .00, .00, .00, .00, .00, .00/
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
C    2 = LDGT1
C    3 = LDGT2
C    4 = HDGV
C  Each line is divided into evap and PCV blocks.  Each block is filled with
C  model year group values and dummy zeroes.
C
      DATA EVPTGS/
     1  .00,1.00, .00,   .00, .00,1.00,
     2  .00,1.00, .00,   .00, .00,1.00,
     3  .00, .05,1.00,   .00,1.00, .00,
     4  .00,1.00, .00,   .00,1.00, .00/
C
C  The CCEMI continuation code is IVTAM (= EVPTGS code above).
C
      DATA CCEMI/
     1  0.00, 1.28, 1.27, 1.26, 1.24, 1.22, 1.21,
     2  0.00, 1.28, 1.27, 1.24, 1.23, 0.00, 0.00,
     3  0.00, 1.28, 1.27, 1.24, 1.23, 0.00, 0.00,
     4  0.00, 1.34, 1.29, 0.00, 0.00, 0.00, 0.00/
C
C  /LOOKUP/ & /MYCODE/ are all indices computed during execution.
C
      DATA IVTAM,IQG,IPG,JPGD,IHG,IGCSF/6*1/
      DATA MY,IDX,JDX,LDXSY,LMYRVT,IAY,IMDXSY,IMKINK/NEYR,7*1/
C
C  /MYUB/: the model year upper bound pointers never change.
C
C  MYTGS continuation line code corresponds to the IQG code:
C    1 = air pump only
C    2 = air pump / catalyst
C    3 = catalyst only
C    4 = EGR only
C    5 = EGR / 3-way catalyst
C    6 = no EGR / 3-way catalyst
C
      DATA MYTGSL/
C  LDGV
     1  1967,1971,1972,1974,1975,1977,1979,NEYR,0000,0000,0000,
     2  1974,1976,1977,1979,1980,1981,1982,1985,1986,NEYR,0000,
     3  1974,1975,1976,1979,1980,1981,1982,1985,1986,NEYR,0000,
     4  1972,1973,1980,1983,NEYR,0000,0000,0000,0000,0000,0000,
     5  1979,1980,1983,1986,NEYR,0000,0000,0000,0000,0000,0000,
     6  1981,1982,1983,1986,NEYR,0000,0000,0000,0000,0000,0000,
C  LDGT1
     1  1967,1971,1972,1974,1980,NEYR,0000,0000,0000,0000,0000,
     2  1974,1976,1978,1980,1981,1983,1985,1987,NEYR,0000,0000,
     3  1974,1975,1976,1978,1980,1981,1983,1985,1987,NEYR,0000,
     4  1972,1973,1978,1980,1981,1982,1983,1984,1985,1986,NEYR,
     5  1980,1981,1982,1983,1984,1985,1986,NEYR,0000,0000,0000,
     6  NEYR,0000,0000,0000,0000,0000,0000,0000,0000,0000,0000,
C  LDGT2
     1  NEYR,0000,0000,0000,0000,0000,0000,0000,0000,0000,0000,
     2  1978,1981,1983,1985,1987,NEYR,0000,0000,0000,0000,0000,
     3  1978,1981,1983,1985,1987,NEYR,0000,0000,0000,0000,0000,
     4  1972,1978,1982,1983,1984,1985,1986,NEYR,0000,0000,0000,
     5  1982,1983,1984,1985,1986,NEYR,0000,0000,0000,0000,0000,
     6  NEYR,0000,0000,0000,0000,0000,0000,0000,0000,0000,0000/
C  HDGV
      DATA MYTGSH/
     1  1986,NEYR,0000,0000,0000,0000,0000,0000,0000,0000,0000,
     2  1986,NEYR,0000,0000,0000,0000,0000,0000,0000,0000,0000,
     3  1981,1986,1987,1989,NEYR,0000,0000,0000,0000,0000,0000,
     4  1984,1985,1986,1989,NEYR,0000,0000,0000,0000,0000,0000,
     5  1986,1989,NEYR,0000,0000,0000,0000,0000,0000,0000,0000,
     6  NEYR,0000,0000,0000,0000,0000,0000,0000,0000,0000,0000/
C
C  MYEGR keys only on vehicle type, even though EMEGR varies with pollutant
C  and bag number.  Also, cell 2 for LDGT2 is a dummy, since there is no
C  NOX 3.1 standard my range for LDGT2.  Otherwise, the 4 ranges for each
C  IVTAM correspond to the no standard, 3.1, 2.0 and 1.0 cases.
C  The continuation line code is the IVTAM code:
C    1 = LDGV
C    2 = LDGT1
C    3 = LDGT2
C    4 = HDGV
C
      DATA MYEGR/
     1  1974,1976,1980,NEYR,
     2  1974,1978,1987,NEYR,
     3  1978,1978,1987,NEYR,
     4  1989,NEYR,0000,0000/
C
C  MYVTS continuation line code is also the IVTAM code.  Each line is
C  divided into the evap and PCV blocks.
C
      DATA MYVTS/
     1  1970,NEYR,0000,   1962,1967,NEYR,
     2  1970,NEYR,0000,   1962,1967,NEYR,
     3  1970,1978,NEYR,   1967,NEYR,0000,
     4  1984,NEYR,0000,   1967,NEYR,0000/
C
C  MYCCEI continuation code is IVTAM (= EVPTGS code above).
C
      DATA MYCCEI/
     1  1963,1970,1974,1977,1979,1980,NEYR,
     2  1963,1977,1979,1980,NEYR,0000,0000,
     3  1963,1977,1979,1980,NEYR,0000,0000,
     4  1963,1980,NEYR,4*0/
C
C  /OMTCF/: all values are computed during execution by BIGCFX.
C
      DATA OMTCF/600*1.0/,OMTTAM/300*0.0/
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
C  Continuation line code for TAMZML, TAMDR, and TAMA50 is IVTAM.
C
C  Each line contains the 9 disablement intercepts (TAMZML) and deterioration
C  rates (TAMDR) for the given vehicle type, myg and I/M case.
C
      DATA TZLNON/
C         AIRS   CATS   NCKS   TNKS   EGRS   EVAP   PCVS   CAPS   MISF
C  Non-I/M & MYG 1
     1   .0072, .0084, .0027,-.0286, .0412, .0192, .0107, .0334,-.0259,
     2   .0320, .0785, .1519, .0000,-.2435,-.0195, .0000,-.1849, .1519,
     3   .0320, .0785, .1519, .0000,-.2435,-.0195, .0000,-.1849, .1519,
     4   .0320, .0785, .1519, .0000,-.2435,-.0195, .0000,-.1849, .1519,
C  Non-I/M & MYG 2
     1   .0033,-.0019,-.0110, .0220, .0225, .0078,-.0036, .0110, .0110,
     2  -.0389,-.0018, .0017, .0097,-.0204, .0009, .0145,-.0590, .0113,
     3  -.0389,-.0018, .0017, .0097,-.0204, .0009, .0145,-.0590, .0113,
     4  -.0389,-.0018, .0017, .0097,-.0204, .0009, .0145,-.0590, .0113,
C  Non-I/M & MYG 3
     1   .0006,-.0008,-.0009, .0009, .0061,-.0013,-.0002,-.0012, .0000,
     2  -.0039, .0024, .0075,-.0075,-.0019, .0100, .0011, .0092, .0000,
     3  -.0039, .0024, .0075,-.0075,-.0019, .0100, .0011, .0092, .0000,
     4  -.0039, .0024, .0075,-.0075,-.0019, .0100, .0011, .0092, .0000/
C
      DATA TZLIM/
C         AIRS   CATS   NCKS   TNKS   EGRS   EVAP   PCVS   CAPS   MISF
C  I/M & MYG 1
     1   .0668, .0067, .0379, .0528, .0744, .0099, .0048,-.0116, .0907,
     2   .0615, .1037, .1376,-.0885, .0808, .0695,-.0156, .0471, .0491,
     3   .0615, .1037, .1376,-.0885, .0808, .0695,-.0156, .0471, .0491,
     4   .0615, .1037, .1376,-.0885, .0808, .0695,-.0156, .0471, .0491,
C  I/M & MYG 2
     1   .0066,-.0152,-.0005, .0112, .0203,-.0158,-.0028,-.0173, .0108,
     2  -.0318,-.0062, .0201,-.0580, .0538,-.0126,-.0261,-.0314,-.0379,
     3  -.0318,-.0062, .0201,-.0580, .0538,-.0126,-.0261,-.0314,-.0379,
     4  -.0318,-.0062, .0201,-.0580, .0538,-.0126,-.0261,-.0314,-.0379,
C  I/M & MYG 3
     1   .0053, .0000, .0000, .0000, .0012,-.0013,-.0018,-.0012, .0000,
     2  -.0262, .0006,-.0122, .0197,-.0021, .0100, .0000, .0092, .0075,
     3  -.0262, .0006,-.0122, .0197,-.0021, .0100, .0000, .0092, .0075,
     4  -.0262, .0006,-.0122, .0197,-.0021, .0100, .0000, .0092, .0075/
C
      DATA TDRNON/
C        AIRS   CATS   NCKS    TNKS   EGRS   EVAP   PCVS   CAPS   MISF
C  Non-I/M & MYG 1 & vmt <= 50K
     1   .0261, .0193, .0242, .0100, .0347, .0097, .0022, .0096, .0342,
     2   .0767, .0444, .0297, .0000, .1198, .0438, .0073, .0839, .0297,
     3   .0767, .0444, .0297, .0000, .1198, .0438, .0073, .0839, .0297,
     4   .0767, .0444, .0297, .0000, .1198, .0438, .0073, .0839, .0297,
C  Non-I/M & MYG 2 & vmt <= 50K
     1   .0048, .0038, .0067,-.0033, .0040, .0013, .0022, .0046, .0034,
     2   .0299, .0258, .0188, .0015, .0258, .0057, .0003, .0270, .0203,
     3   .0299, .0258, .0188, .0015, .0258, .0057, .0003, .0270, .0203,
     4   .0299, .0258, .0188, .0015, .0258, .0057, .0003, .0270, .0203,
C  Non-I/M & MYG 3 & vmt <= 50K
     1   .0028, .0014, .0006, .0019, .0017, .0026, .0013, .0055, .0025,
     2   .0088, .0035, .0006, .0051, .0053, .0011, .0009, .0047, .0057,
     3   .0088, .0035, .0006, .0051, .0053, .0011, .0009, .0047, .0057,
     4   .0088, .0035, .0006, .0051, .0053, .0011, .0009, .0047, .0057/
C
      DATA TDRIM/
C        AIRS   CATS   NCKS    TNKS   EGRS   EVAP   PCVS   CAPS   MISF
C  I/M & MYG 1 & vmt <= 50K
     1   .0055, .0084, .0089,-.0038, .0105, .0085, .0023, .0156, .0051,
     2   .0304, .0203, .0210, .0178, .0114, .0093, .0040, .0150, .0388,
     3   .0304, .0203, .0210, .0178, .0114, .0093, .0040, .0150, .0388,
     4   .0304, .0203, .0210, .0178, .0114, .0093, .0040, .0150, .0388,
C  I/M & MYG 2 & vmt <= 50K
     1   .0029, .0060, .0009, .0008, .0014, .0085, .0021, .0122, .0017,
     2   .0148, .0109, .0031, .0198, .0000, .0128, .0079, .0184, .0229,
     3   .0148, .0109, .0031, .0198, .0000, .0128, .0079, .0184, .0229,
     4   .0148, .0109, .0031, .0198, .0000, .0128, .0079, .0184, .0229,
C  I/M & MYG 3 & vmt <= 50K
     1   .0011, .0003, .0002, .0034, .0006, .0026, .0014, .0055, .0036,
     2   .0164, .0043, .0081, .0005, .0027, .0011, .0000, .0047, .0086,
     3   .0164, .0043, .0081, .0005, .0027, .0011, .0000, .0047, .0086,
     4   .0164, .0043, .0081, .0005, .0027, .0011, .0000, .0047, .0086/
C
      DATA TA5NON/
C        AIRS   CATS   NCKS    TNKS   EGRS   EVAP   PCVS   CAPS   MISF
C  Non-I/M & MYG 1 & vmt > 50K
     1  1.5794, .9731, 1.0913,-.0211, .3245,2.0404, 1.00, 2.3463, .7653,
     2   .2207, .9810, 1.6721, .0000, .0071, .3830, 1.00, .1806, 1.6721,
     3   .2207, .9810, 1.6721, .0000, .0071, .3830, 1.00, .1806, 1.6721,
     4   .2207, .9810, 1.6721, .0000, .0071, .3830, 1.00, .1806, 1.6721,
C  Non-I/M & MYG 2 & vmt > 50K
     1  7.4853, 3.6154, .9121,-.8487,2.3522,1.0129, .2863,2.8154,2.6554,
     2  1.0969, .5767, .5832, 2.5794, .5597,3.8992,1.4240, .5891, .7318,
     3  1.0969, .5767, .5832, 2.5794, .5597,3.8992,1.4240, .5891, .7318,
     4  1.0969, .5767, .5832, 2.5794, .5597,3.8992,1.4240, .5891, .7318,
C  Non-I/M & MYG 3 & vmt > 50K
     1  3.7188, 1.0400, 1.00, 1.000, 1.3764, .3976, .3986, .2305, 1.000,
     2  1.5120, 4.0011, 1.00, 1.00, 2.0764,  10.39, 1.00, 3.2234, 1.000,
     3  1.5120, 4.0011, 1.00, 1.00, 2.0764,  10.39, 1.00, 3.2234, 1.000,
     4  1.5120, 4.0011, 1.00, 1.00, 2.0764,  10.39, 1.00, 3.2234, 1.000/
C
      DATA TA5IM/
C        AIRS   CATS   NCKS    TNKS   EGRS   EVAP   PCVS   CAPS   MISF
C  I/M & MYG 1 & vmt > 50K
     1  5.6555, 1.000, 1.000, 1.000, 1.000, 1.4620, 1.000, .8646, 1.000,
     2   .9643, 1.000, 1.000, 1.000, 1.000, 1.000, 3.1065, 1.000, 1.000,
     3   .9643, 1.000, 1.000, 1.000, 1.000, 1.000, 3.1065, 1.000, 1.000,
     4   .9643, 1.000, 1.000, 1.000, 1.000, 1.000, 3.1065, 1.000, 1.000,
C  I/M & MYG 2 & vmt > 50K
     1  5.4872, 1.000, .4492,-.3720, 1.1930, .0310, 1.000, .1665, .0800,
     2  1.0000, 1.000, 1.000, .7472, .0000, .0119, 1.0000, .0344, .7812,
     3  1.0000, 1.000, 1.000, .7472, .0000, .0119, 1.0000, .0344, .7812,
     4  1.0000, 1.000, 1.000, .7472, .0000, .0119, 1.0000, .0344, .7812,
C  I/M & MYG 3 & vmt > 50K
     1   .5021, 1.000, 1.000, 1.000, 1.5360, .3976,2.5328, .2305, 1.000,
     2  1.000, 5.3523, 1.000, 1.000, 2.0765, 10.39, .0000,3.2234, 1.000,
     3  1.000, 5.3523, 1.000, 1.000, 2.0765, 10.39, .0000,3.2234, 1.000,
     4  1.000, 5.3523, 1.000, 1.000, 2.0765, 10.39, .0000,3.2234, 1.000/
C
C  /TAMEQ2/: tampering rate equation parameters used to select and apply
C            tampering rate slope and intercept.
C
C  Continuation line code for MYGTAM is IVTAM.
C
      DATA MYGTAM/
C                          MYG
C                    1      2      3
     1             1980,  1983,  NEYR,
     2             1980,  1983,  NEYR,
     3             1980,  1983,  NEYR,
     4             1980,  1983,  NEYR/
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
C            and THS, TDU & TCC by EVPEME.
C
      DATA TAMBAG/900*0.0/,THS,TDU,TCC/400*0.0/
C
C  /TAMOU2/: all values are computed during execution by EVPEME.
C
      DATA TOB,GCONLY/200*0.0/
C
C  /TAMPB1/: LPOD is fixed at fixed-by-law 1995, PBBTR is computed in SAVEPB.
C
      DATA LPOD/1995/, PBBTR/600*0.0/
C
      END
