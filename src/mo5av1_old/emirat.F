C Version @(#)$Id$ $Source$ $Date$ 

      SUBROUTINE EMIRAT
C
C  EMIRAT computes / assigns emission impact rates for the tampering effects
C  groups.
C
C  Called by TAMPER.
C
C  Calls ITAMPT, OUTPOL.
C
C  Input on call:
C
C    common blocks:
C    /EGSCAL/ TGS
C    /EMECAL/ EM3W,EMOX,EMEGR
C    /FLAGS4/ IDLFLG
C    /LOOKUP/ IVTAM
C    /MAXIMA/ MAXPOL
C    /TAMID1/ EMI3W,EMIOX
C
C  Output on return:
C
C    common blocks:
C    /EMECAL/ EMI
C    /LOOKUP/ IQG
C    /TAMID1/ EMIDL
C
C  Local variable dictionary:
C
C   Name   Type              Description
C  ------  ----  -------------------------------------------------------
C  POX      R    percentage of vehicles equipped with oxidation catalysts
C  POXN     R    POX normalized to ox cat equipped / all cat equipped
C  P3W      R    percentage of vehicles equipped with 3-way catalysts
C  P3WN     R    P3W normalized to 3-way cat equipped / all cat equipped
C  TK       R    percentage of fleet equipped with catalysts
C
C  Notes:
C
C  OPENLP was removed in MOBILE4.1 as part of deleting old oxyfuels support.
C
C
      INTEGER PRTFLG,HCFLAG,COLDFG
C
      COMMON /EGSCAL/ AER(11,11,2,2),TGS(11,6,4),EGS(7,2)
      COMMON /EMECAL/ EM3W(3,3,3),EMOX(3,3,3),EMEGR(4,3,3),EMI(7,3,3)
      COMMON /FLAGS4/ PRTFLG,IDLFLG,NMHFLG,HCFLAG,COLDFG
      COMMON /MAXIMA/ MAXVEH,MAXLTW,MAXPOL,MAXREG,MAXYRS
      COMMON /TAMID1/ EMI3W(3,3,3),EMIOX(3,3,3),EMIDL(3,4,3)
      COMMON /TAMID2/ MYIDLE(3,3),IEK
      COMMON /LOOKUP/ IVTAM,IQG,IPG,JPGD,IHG,IGCSF
C
C  Compute the percentage of vehicles equipped with 3-way and oxidation
C  catalysts.
C
      P3W=0.0
      DO 10 IQG=5,6
      IG1=ITAMPT(1)
      P3W=P3W+TGS(IG1,IQG,IVTAM)
   10 CONTINUE
C
      POX=0.0
      DO 20 IQG=2,3
      IG1=ITAMPT(1)
      POX=POX+TGS(IG1,IQG,IVTAM)
   20 CONTINUE
      POX=POX-P3W
C
C  Normalize the percentages in order to use them as weights.
C
      TK=P3W+POX
C
C  Avoid divide by zero runtime error.
C
      IF(TK.EQ.0.0) TK=1.0
      P3WN=P3W/TK
      POXN=POX/TK
C
C  Determine pollutants
C
      CALL OUTPOL(IP1,IP2)
C
C  Compute / assign emissions impact for each bag for each pollutant.
C
      DO 50 IP=IP1,IP2
      DO 40 IB=1,3
C
C  Air pump disablement (air pump only equipped) emission impact is the
C  oxidation catalyst emission for the given bag and pollutant.  These
C  numbers are constant with respect to model year and vehicle type and
C  therefore are hardcoded in BLOCK DATA into EMI(1,IB,IP).  Note that
C  EMI(1,IB,IP) = EMOX(1,IB,IP) for all IB,IP.  If EMI for air pump
C  disablement becomes variable in MOBILE4+, the code will be inserted here.
C
C  Next 3 effects groups use the catalyst type weighted emissions.
C
      DO 30 IEG=2,4
      IEK=IEG-1
      EMI(IEG,IB,IP)=P3WN*EM3W(IEK,IB,IP)+POXN*EMOX(IEK,IB,IP)
   30 CONTINUE
C
C  Get EGR only MY group pointer and look up corresponding emission impact.
C
      IG2=ITAMPT(2)
      EMI(5,IB,IP)=EMEGR(IG2,IB,IP)
C
C  EGR disabled / catalyst removed and EGR disabled / misfueled emissions are
C  a constant function of bag and pollutant.  As with group 1, the values are
C  hardcoded in BLOCK DATA into EMI(6,IB,IP) and EMI(7,IB,IP), respectively.
C  Future versions of the model may vary these emissions on some other
C  parameter(s) such as model year and vehicle type, in which case the
C  computation / assignment of values would go here.
C
   40 CONTINUE
   50 CONTINUE
C
C  If idle ef requested do emission impact rates for idle tampering effects.
C
      IF(IDLFLG.EQ.1) GOTO 80
      DO 70 IP=1,MAXPOL
      DO 60 IEG=2,4
      IEK=IEG-1
C
C  Get model year group pointer for idle tampering rates
C
      IGTIDL=ITAMPT(5)
      EMIDL(IGTIDL,IEG,IP)=P3WN*EMI3W(IGTIDL,IEK,IP)+
     *      POXN*EMIOX(IGTIDL,IEK,IP)
      IF(IEG.EQ.2) EMIDL(IGTIDL,1,IP)=EMIDL(IGTIDL,IEG,IP)
   60 CONTINUE
   70 CONTINUE
C
   80 RETURN
      END
