C Version @(#)$Id$ $Source$ $Date$ 

      FUNCTION CCEVRT(MY,IDX,IV,IY)
C
C  CCEVRT is a REAL-valued function that returns a hydrocarbon combined
C  evaporative (hot soak + diurnal) and crankcase emission factor.  The
C  hot soak value is weighted by trips per day and the evaporative sum is
C  weighted by miles per day.
C
C  Called by HCCALX.
C
C  Calls CRANKC, DIURNL, HOTSOK, IEVPTR, and PCTLEV.
C
C  Input on call:
C
C    parameter list: MY,IDX,IV
C    common blocks:
C    /CEVBMY/ BMYMPD,BMYTPD,BMTPDF
C    /FLAGS4/ HCFLAG
C    /INJECT/ TBI,PFI
C    /LEVBLK/ LEVFLG, LEVYRS
C    /MAXIMA/ MAXYRS
C    /TAMOUT/ TCC
C    /VMXCOM/ TFNORM
C
C  Output on return:
C
C    common blocks:
C    /EVAPAR/ IFDS,ISTDC,ISTDT,ISTDP,FINJ  (to IEVPTR, HOTSOK, DIURNAL)
C    /EVAPGR/ EVP,PARD
C    function: CCEVRT
C
C  Local variable dictionary:
C
C   Name   Type              Description
C  ------  ----  -------------------------------------------------------
C  CC       R    crankcase HC effect = CRANKC + tampering cc offset, if any
C  DI       R    diurnal   HC effect = DIURNL + tampering di offset, if any
C  HS       R    hot soak  HC effect = HOTSOK + tampering hs offset, if any
C  IG5      I    pointer to fuel injected sales percentage for given my
C  IVGAS    I    see (1,999)'s subscript dictionary
C  MPD      R    miles per day
C  Z        R    discount evap emission by ZEV fraction
C
C
C  Notes:
C
C  Hot soak, diurnal and crankcase without tampering and with tampering HC are
C  figured separately and then added together, according to Evap Study 2.
C  The additions are done before miles per day and trips per day are applied
C  as weights.
C
C  CCEVRT was modified in MOBILE4.1 to transform hot soaks, diurnals
C  and partial diurnals to grams per mile (g/mi), when IV is LDGV/T
C  and MPD is positive.    Also, oxygenated and cert. fuels support
C  from MOBILE4 was removed in MOBILE4.1.
C
C
      INTEGER ATPFLG,TPDFLG,RLFLAG,TEMFLG,OUTFMT
      INTEGER PRTFLG,HCFLAG,COLDFG
      REAL MPD
C
      COMMON /CEVBMY/ BMYMPD(25,8),BMYTPD(25,8),BMTPDF(25,8)
      COMMON /CITCIN/ UDI(5),IUDI,PFRATE(3,2),IPFCU
      COMMON /EVAPAR/ IFDS,ISTDC,ISTDT,ISTDP,FINJ(2)
      COMMON /EVAPGR/ EVP(4),GREVP(4,9),VGREVP(4),PARD(3),GRPD(3,9)
      COMMON /FLAGS3/ ATPFLG,TPDFLG,RLFLAG,LOCFLG,TEMFLG,OUTFMT
      COMMON /FLAGS4/ PRTFLG,IDLFLG,NMHFLG,HCFLAG,COLDFG
      COMMON /INJECT/ TBI(13,4),PFI(13,4),MAXFIY,MINFIY
      COMMON /LEVBLK/ LEVFLG,IMFLG,LEVYRS
      COMMON /MAXIMA/ MAXVEH,MAXLTW,MAXPOL,MAXREG,MAXYRS
      COMMON /REGION/ FEET(2),IREJN,ALT,INITPR
      COMMON /TAMOUT/ TAMBAG(3,3,25,4),THS(2,25,4),TDU(25,4),TCC(25,4)
      COMMON /TPDCOM/ TPDCO(2,8)
      COMMON /VMXCOM/ REGMIX(8),TFNORM(8),VMTMIX(8),VCOUNT(39,8)
c
c  LMOS Modification 04/08/93
c  Add LMOS common block
c
      COMMON /LMOS/ EFDNL(9),EFEVND(9),VDNL,DNLGPM
c
C
C  Set ZEV percentage.  As was done with BEFs the ZEV is taken
C  out first and then the ZEV adjusted number is corrected for
C  FUEL adjustment
C
      Z=1.0
      IF(LEVFLG.EQ.2.AND.IV.LE.2.AND.MY.GE.LEVYRS) THEN
        Z=1.0-PCTLEV(MY,IV,8)
      ENDIF
C
C  Set parameters that are identical for hotsoak and diurnal calculations.
C
      IVGAS=IV
      IF(IV.EQ.8) IVGAS=5
C
C  Determine basic miles per day according to vehicle class.
C
      JDX=MAXYRS-IDX+1
C
      IF(IVGAS.LE.3) MPD=BMYMPD(JDX,IV)
      IF(IVGAS.GE.4) MPD=TFNORM(IV)*100000./365.0
C
C  Warn if MPD is zero.
C
      IF(MPD.EQ.0.0) CALL QUITER(FLOAT(MY),IV,58,INERR)
C
C  Look up HS/DU model year group pointers:
C
C  (1) 1 emission standard code for each fuel delivery system (IFDS) -
C          1 = carbureted (ISTDC)
C          2 = TBI (ISTDT)
C          3 = PFI (ISTDP)
C
C  (2) 1 for the percentage of vehicle type's fleet that is TBI (FINJ(1))
C      and PBI (FINJ(2)). (IG5 dimension of TBI and PFI
C      (IG5 dimension of TBI and PFI).
C
C  Carb
      IFDS=1
      ISTDC=IEVPTR(MY,IVGAS,1)
C
      IF(IVGAS.GE.4.OR.MY.LT.1981) GOTO 10
C
C  TBI
      IFDS=2
      ISTDT=IEVPTR(MY,IVGAS,1)
C
      IG5=IEVPTR(MY,IVGAS,2)
      FINJ(1)=TBI(IG5,IVGAS)
C
C  PFI
      IFDS=3
      ISTDP=IEVPTR(MY,IVGAS,1)
C
      IG5=IEVPTR(MY,IVGAS,2)
      FINJ(2)=PFI(IG5,IVGAS)
C
C  Call functions to calculate hot soak and diurnal and lookup crankcase
C  rates.  The first two are returned with the tampering offset already
C  added in, given LDGV/T or HDGV. Hot soak and crank case are zero when
C  MPD is zero. Diurnals occur regardless of MPD.
C
   10 HS=0.0
      EVP(1)=0.0
      IF(MPD.GT.0.0) HS=HOTSOK(MY,IDX,IVGAS,IY)*Z
C
C     VOC/TOG adjustment for oxy fuels
C
      CALL FUEL2(2,HS)
      EVP(1)=EVP(1)*Z
      CALL FUEL2(2,EVP(1))
C
      DU=DIURNL(MY,IDX,IVGAS,IY)*Z
      PARD(1)=PARD(1)*Z
      PARD(2)=PARD(2)*Z
      PARD(3)=PARD(3)*Z
      CALL FUEL2(1,PARD(1))
      CALL FUEL2(1,PARD(2))
      CALL FUEL2(1,PARD(3))
      CALL FUEL2(1,DU)
      EVP(2)=EVP(2)*Z
      EVP(3)=EVP(3)*Z
      CALL FUEL2(1,EVP(2))
      CALL FUEL2(1,EVP(3))
C
      CC=0.0
      IF(MPD.GT.0.0) CC=CRANKC(MY,IV)
C
C  GV/T => tampering crankcase HC effects are to be added.
C
      IF(IVGAS.LE.4.AND.MPD.GT.0.0) CC=CC+TCC(IDX,IV)
      EVP(4)=CC*Z
C
C  CCEVRT is also zero if MPD is zero.
C
      CCEVRT=0.0
      IF(MPD.EQ.0.0) GOTO 90
C
C  Combining the 3 factors requires converting hotsoak from grams per trip to
C  grams per day (gpd), adding in diurnal (which by definition only occurs once
C  per day), converting the sum from gpd to grams per mile (gpm), and, finally,
C  adding in crankcase (which is already in gpm).
C
c
c   LMOS Modification 04/08/93
c   Separate out diurnals from CCEVRT.
c
c      CCEVRT=(HS*BMYTPD(JDX,IV)+DU)/MPD+CC
      CCEVRT=(HS*BMYTPD(JDX,IV))/MPD+CC
      DNLGPM=DU/MPD
c
C
C  For HCFLAG=3,5 (grams output), generate a composite grams/trip
C  fleet hot soak by applying the trips per day fractions.
C
      IF(HCFLAG.EQ.3.OR.HCFLAG.EQ.5) EVP(1)=EVP(1)*BMTPDF(JDX,IV)
C
C  For HCFLAG=4 (grams/mile output), convert the fleet hot soak,
C  diurnal, and multiple diurnal hc values to g/mi.
C
      IF(HCFLAG.NE.4) GOTO 90
C
      EVP(1)=EVP(1)*BMYTPD(JDX,IV)
      DO 20 IEVP=1,3
      EVP(IEVP)=EVP(IEVP)/MPD
   20 CONTINUE
C
C  Convert the LDGV/T partial diurnals to g/mi also.
C
      IF(IV.GT.3) GOTO 90
      DO 30 IDU=2,4
      IPARTL=IDU-1
      PARD(IPARTL)=PARD(IPARTL)/MPD
   30 CONTINUE
C
   90 RETURN
      END
