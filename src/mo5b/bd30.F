C Version "@(#)$Id$ $Source$ $Date$ 

      BLOCK DATA BD30
C
C  BLOCK DATA Subprogram 30: /CITCIN/,/CITPAR/,/CITRV1/,/CITRV2/, /INJECT/,
C                            /RFORM1/, /RFORM2/
C
C  BLOCK DATA 30 defines the parameters for reading in, storing and using
C  the second required scenario record (local data) parameters, plus fuel
C  injection rates and model year groups.
C
C  Common block array subscripts:
C
C  PFI(13,4)  -  PFI ( IG5, IVTAM )
C  PFRATE(3,2)- PFRATE  (KEYEQ , IPFCU)
C  SCNAME(4)  -  SCNAME ( ICH )
C  TBI(13,4)  -  TBI ( IG5, IVTAM )
C  UDI(5)     -  UDI ( IUDI )
C
C  Common block dictionary:
C
C   Name   Type              Dictionary
C  ------  ----  ------------------------------------------------
C  /CITCIN/:
C
C   UDI     R    Uncontrolled Diurnal rate, g/test, at:
C                1 = standard FTP conditions
C                2 = local area conditions
C                3 = 8AM-11AM conditions
C                4 = 10AM-3PM conditions
C                5 = 8AM-2PM conditions
C  PFRATE   R    The rates for Pass, Failed Purge, and Failed Pressure
C                for calculated (IPFCU=1) or user-supplied (IPFCU=2)
C
C  /CITPAR/:
C
C  SCNAME   R    scenario identification field: echoed in numeric formats
C
C  /CITRV1/:
C
C  RVPBAS   R    SCNAME's base RVP
C  RVPIUS   R    SCNAME's in-use RVP
C  RVPICY   R    RVP to be used for not cert controlled my for SCNAME:
C                    = RVPBAS if ICY .LT. IUSESY
C                    = RVPIUS if ICY .GE. IUSESY
C  IUSESY   I    start year of in-use fuel volatility control program
C  RVPUWX   R    un-weathered RVPICY
C
C  /CITRV2/:
C
C  RVPHS    R    RVP in untampered/uncontrolled hot soak
C  RVPX     R    The RVP for use in calculating the RVP exhaust (and depending
C                on the case, temperature) correction factor.  RVPX initially
C                equals the corresponding evaporative RVPHS.
C  RVP090   R    RVP in RVPCER cases prior to ICERSY, in denominator of Diurnal
C                Index equation & to normalize RVPCER.  ASTM Class A RVP UB.
C  RVP115   R    Uncontrolled RVP.  Also, ASTM Class C RVP UB.
C
C  /INJECT/:
C
C  TBI      I    throttle body fuel injected fraction of IVGAS built in IG5
C  PFI      I    port fuel injected fraction of veh class IVGAS built in IG5
C  MAXFIY   I    first year of highest myg for indexing TBI & PFI.
C  MINFIY   I    last year of lowest myg for indexing TBI & PFI.
C
C  /RFORM1/:
C
C  IASTM    I    corresponds to the character variable for ASTM class
C                1 = A, 2 = B, 3 = C, 4 = D, 5 = E.
C  IPHASE   I    corresponds to the reformulated fuel phase
C                1 = Phase I (simple model)
C                2 = Phase I (complex model)
C                3 = Phase II (complex model)
C  RFGFLG   I    flag indicating a reformulated gas area
C                1 = no reform program
C                2 = apply reform program
C                3 =  apply reform program phase I only
C  RFGON    L    flag indicating a reformulated gas area with a calender
C                year of evaluation that is >= reform phase year
C  IAFCNT   I    base fuel indicator
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
C  RFGSYR   I    start year for phase of reformulated fuel effects
C  RFGRVP   R    reformulated RVP
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
C
      INCLUDE 'CITCIN.I'
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
C  Assigned in LOCAL/RDCITY/CALUDI as follows:
C
C  supplied by user: SCNAME(4),RVPBAS,RVPIUS,IUSESY
C
C  calculated by code: RVPICY,RVPHS,RVPX,PFRATE
C
C  The only constants are:
C  the controlled and uncontrolled RVP standard values of /CITRV2/:
C
      DATA RVP090/9.0/,RVP115/11.5/,RVPLOS/0.0/
C
C  Dummy zeroes are assigned the remaining variables and arrays.
C
      DATA UDI,PFRATE,IPFCU/11*0.0,1/
      DATA SCNAME/4*'    '/
      DATA RVPBAS,RVPIUS,RVPICY,IUSESY,RVPUWX/5*0.0/
      DATA RVPHS,RVPX/2*0.0/
C
C  The reformulated fuel variables.
C
      DATA IASTM,IPHASE,RFGFLG,IAFCNT,WINFLG,SEAFLG,RFGON
     *     /3*1,999,1,1,.FALSE./
C
C  Start Year by phase; RVP by ASTM class and phase.
C
C  RFGSYR(3) is reset to 2000 in RDCITY because it might have been
C  modified there;  12/15/93 TMM@ARC Request #2-366.
C  Oct-20-1994 @ CSC-bsg Request 426  Summer RVP values for Phase 2
C  for Reformulated Gasoline have changed from the values found in MOBILE5a.
C
C
      DATA RFGSYR/ 1995, 1997, 2000/
      DATA RFGRVP/ 2*7.1,3*8.0, 2*7.1,3*8.0, 5*6.7/
C
C  /INJECT/: no provision exists for user modification of the fuel injection
C            fractions or their model year groups.
C
      DATA MAXFIY/1992/,MINFIY/1980/
C
C  Both TBI and PFI have the IVTAM continuation code:
C    1 = LDGV
C    2 = LDGT1
C    3 = LDGT2
C    4 = HDGV
C  Each line has 8 (model year group) values
C
C  Observe:
C
C    1) LDGT1 = LDGT2 values
C    2) All HDGV values = 0.0/0, since for now HDGV's are assumed to
C       be 100% carbureted.
C    3) For tampering, the fractions are only used to weight hot soak rates.
C    4) Before 1981, it is assumed that no vehicle of any type had
C       fuel injection and therefore hot soak rates are for those years
C       based solely on evaporative cannister only tampering.
C
C    pre-81 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992+
C
      DATA TBI/
     1 .000,.029,.106,.183,.282,.208,.284,.369,.407,.275,.219,.202,.202,
     2 .000,.000,.000,.000,.000,.047,.136,.280,.452,.369,.367,.407,.407,
     3 .000,.000,.000,.000,.000,.047,.136,.280,.452,.369,.367,.407,.407,
     4 .000,.000,.000,.000,.000,.000,.000,.000,.000,.000,.000,.000,.000/
C
      DATA PFI/
     1 .000,.061,.062,.088,.110,.307,.392,.372,.492,.597,.761,.795,.795,
     2 .000,.000,.000,.002,.022,.066,.238,.327,.416,.540,.601,.577,.577,
     3 .000,.000,.000,.002,.022,.066,.238,.327,.416,.540,.601,.577,.577,
     4 .000,.000,.000,.000,.000,.000,.000,.000,.000,.000,.000,.000,.000/
C
      END

