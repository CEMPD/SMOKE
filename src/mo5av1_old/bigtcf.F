C Version @(#)$Id$ $Source$ $Date$ 

      SUBROUTINE BIGTCF(ICY,IP,IVTEMP,MY,IDX)
C
C  BIGTCF computes the multiplicative temperature correction factor and, in
C  some cases, the RVP cf and/or the additive temperature CO offset.
C
C  Called by BIGCFX.
C
C  Calls ICCOPT, IGRVPT, ITCPTR, IEVPTR, PIV and BIGRVP.
C
C  Input on call:
C
C    parameter list: IP,IVTEMP,MY,IDX
C    common blocks:
C    /CITRV2/ RVPX
C    /FLAGS4/ COLDFG
C    /INJECT/ PFI
C    /TEMPC1/ TTFC,TTFD,TT4,TT8
C    /TEMPC2/ MYGTFD,ISHIFT
C    /TEMPC3/ MDLOHI,TDIFF,CHFRAC
C    /TEMPC4/ MYGCOO,LOWCO,COOTFC,COOTFD
C    /TEMPC5/ NCC,MYCOL
C    /TEMPC6/ CSUB,CMUL
C    /TEMPS/  TEMEXH
C
C  Output on return:
C
C    common blocks:
C    /OFFSET/ OFFCO
C    /TEMPC3/ TCF
C
C  Local array subscripts:
C
C  FRFDS(3)  -  FRFDS ( IFDS )
C
C  Local variable / array dictionary:
C
C   Name   Type              Description
C  ------  ----  -------------------------------------------------------
C  COOFDS   R    sum of CO offsets by FDS weighted by FDS sales fractions
C  FRFDS    R    fraction of class IVTEMP given MY that is IFDS equipped
C  HDGTCF   R    HDGV T CF - HDGV testing does not use the LDG operating modes,
C                so multiply the 1 HDGTCF x each of the RVP CFs in TCF(IB).
C  TC1      R    LDGV/T tcf before checked for cold CO std adjustment:
C                no std applies = use as is, std applies = weight in std impact
C  TC2      R    LDGV/T tcf checked & adjusted for cold CO std adjustment
C  TFCB     R    tcf by bag for fuel delivery systems combined case
C  TFDS     R    FDS weighted low tcf bag equation coefficient
C
C  Notes:
C
C  For low tcf, LDGV/T use the additive CO offset starting in MYGCOO(IVTEMP) &
C  are differentiated by fuel delivery system (carbureted vs fuel injected)
C  starting in MYGTFD(IVTEMP).  The CO offset start year precedes the fds start
C  year in all 3 cases.  Hence, BIGTCF checks the LOWCO flag on both the LDG
C  MYGTFD keyed branches.
C
C  BIGTCF was added for MOBILE4.
C  BIGTCF was modified for MOBILE4, adding BIGRVP.
C  Oxygenated and cert. fuels support from MOBILE4 was removed in MOBILE4.1.
C
      INTEGER PRTFLG,HCFLAG,COLDFG
      INTEGER ADDCSW
C
      COMMON /CITRV2/ RVPHS,RVPX,RVP090,RVP115,RVPLOS
      COMMON /FLAGS4/ PRTFLG,IDLFLG,NMHFLG,HCFLAG,COLDFG
      COMMON /INJECT/ TBI(13,4),PFI(13,4),MAXFIY,MINFIY
      COMMON /OFFSET/ OFFCO(25,3),OFFMTH(25,8,2)
      COMMON /TEMPC1/TTFC(2,3,9,3),TTFD(3,3,3),TT4(2,6,3),TT8(2,3,4,3)
      COMMON /TEMPC2/ MAXTCF,MYGTCF(9,3,5),MYGTFD(3),ISHIFT
      COMMON /TEMPC3/ MDLOHI(3),TDIFF(3),CHFRAC(3),TCF(3)
      COMMON /TEMPC4/ MYGCOO(3),LOWCO,COOTFC(3),COOTFD(3)
      COMMON /TEMPC5/ NCC,MYCOL(2),CCSTD(3,2)
      COMMON /TEMPC6/ CSUB(3,2),CMUL(3,2),PIVPCT(3,3,2)
      COMMON /TEMPC7/ ICOLD,ADDCSW,MULCSW
      COMMON /TEMPS/ AMBT,TEMMIN,TEMMAX,TEMEXH(3),TEMEVP(6)
C
      DIMENSION TFCB(3),TFDS(3)
C
      DIMENSION FRFDS(3)
C
C  Initialize TCF to the neutral (no correction) value.
C
      DO 10 IB=1,3
      TCF(IB)=1.0
   10 CONTINUE
C
C  First cases evaluated are the new combined high temperature and RVP cf
C  algorithm and bag dependent low temperature RVP cf cases.  Check these first
C  because their year ranges are mostly likely to be evaluated.
C
      IF(IVTEMP.GT.4) GOTO 30
      IG=IGRVPT(MY,IVTEMP)
      IF(IG.LE.1.OR.TEMEXH(IP).LE.45.0) GOTO 30
      CALL BIGRVP(ICY,IP,IG)
C
C  High temperature and RVP combined cf case is done.
C
      IF(MDLOHI(IP).EQ.2) RETURN
C
C  All RVP cf calculations are done.  If T is neutral, no T cf occurs.
C  For the IG >= 2 case returned from BIGRVP, pass on the RVP cf in TCF.
C
   30 IF(MDLOHI(IP).EQ.0) RETURN
C
C  At this point, TCF(IB) = 1.0, if IV = MC or IG <= 1 or (IG > 1 & T = low)
C                           RVP cf, if IG >= 2 & T < 75.0
C  Hence, the formulas below for the multiplicative tcf cases are all of the
C  form "TCF(IB)=f(---)*TCF(IB)".  The result is either a tcf or the product
C  by bag tcf * RVP cf.
C
      GOTO(40,40,40,70,80),IVTEMP
C
C  Low for all MY & high for RVP MYGs 0 & 1 TCF for LDGV, LDGT1, and LDGT2.
C
   40 IF(IP.EQ.2.AND.MDLOHI(IP).EQ.1.AND.MY.GE.MYGCOO(IVTEMP)) LOWCO=2
      ADDCSW=1
      MULCSW=1
      IF(IP.EQ.2.AND.MDLOHI(IP).EQ.1.AND.COLDFG.GT.1) CALL ICCOPT(MY)
      IF(MY.GE.MYGTFD(IVTEMP)) GOTO 50
C
C  IVLDG multiplicative low & high tcf with FDS combined.
C
      IT=ITCPTR(MY,IP,IVTEMP)
      DO 45 IB=LOWCO,3
      TFCB(IB)=TTFC(MDLOHI(IP),IB,IT,IP)
      TC1=EXP(TFCB(IB)*TDIFF(IP))
      IF(MULCSW.EQ.1) GOTO 43
      TC2=((TC1-1.0)*CMUL(IVTEMP,ICOLD))+1.0
      TC1=TC1*(1.0-PIV(MY,IVTEMP))+TC2*PIV(MY,IVTEMP)
   43 TCF(IB)=TC1*TCF(IB)
   45 CONTINUE
      IF(LOWCO.EQ.1) RETURN
C
C  IVLDG additive low temperature CO offset with FDS combined.
C
      IGO=IT-ISHIFT
      TC1=COOTFC(IGO)*TDIFF(IP)*CHFRAC(1)
      IF(ADDCSW.EQ.1.AND.MULCSW.EQ.1) GOTO 47
      IF(ADDCSW.EQ.2) TC2=TC1-CSUB(IVTEMP,ICOLD)
      IF(MULCSW.EQ.2) TC2=TC1*CMUL(IVTEMP,ICOLD)
      TC1=TC1*(1.0-PIV(MY,IVTEMP))+TC2*PIV(MY,IVTEMP)
   47 OFFCO(IDX,IVTEMP)=TC1
      RETURN
C
C  IVLDG multiplicative low tcf with FDS differentiated.
C
   50 IG5=IEVPTR(MY,IVTEMP,2)
      FRFDS(2)=TBI(IG5,IVTEMP)
      FRFDS(3)=PFI(IG5,IVTEMP)
      FRFDS(1)=1.0-FRFDS(2)-FRFDS(3)
      DO 60 IB=LOWCO,3
      TFDS(IB)=0.0
      DO 55 IFDS=1,3
      TFDS(IB)=TFDS(IB)+FRFDS(IFDS)*TTFD(IB,IP,IFDS)
   55 CONTINUE
      TC1=EXP(TFDS(IB)*TDIFF(IP))
      IF(MULCSW.EQ.1) GOTO 57
      TC2=((TC1-1.0)*CMUL(IVTEMP,ICOLD))+1.0
      TC1=TC1*(1.0-PIV(MY,IVTEMP))+TC2*PIV(MY,IVTEMP)
   57 TCF(IB)=TC1*TCF(IB)
   60 CONTINUE
      IF(LOWCO.EQ.1) RETURN
C
C  IVLDG additive low temperature CO offset with FDS differentiated.
C
      COOFDS=0.0
      DO 65 IFDS=1,3
      COOFDS=COOFDS+FRFDS(IFDS)*COOTFD(IFDS)
   65 CONTINUE
      TC1=COOFDS*TDIFF(IP)*CHFRAC(1)
      IF(ADDCSW.EQ.1.AND.MULCSW.EQ.1) GOTO 67
      IF(ADDCSW.EQ.2) TC2=TC1-CSUB(IVTEMP,ICOLD)
      IF(MULCSW.EQ.2) TC2=TC1*CMUL(IVTEMP,ICOLD)
      TC1=TC1*(1.0-PIV(MY,IVTEMP))+TC2*PIV(MY,IVTEMP)
   67 OFFCO(IDX,IVTEMP)=TC1
      RETURN
C
C  Low for all MY & high for RVP MYG 0 (pre-1985) TCF for HDGV.
C
   70 IT=ITCPTR(MY,IP,IVTEMP)
      HDGTCF=EXP(TT4(MDLOHI(IP),IT,IP)*TDIFF(IP))
      DO 75 IB=1,3
      TCF(IB)=HDGTCF*TCF(IB)
   75 CONTINUE
      RETURN
C
C  Lo/hi TCF for MC.
C
   80 IT=ITCPTR(MY,IP,IVTEMP)
      DO 85 IB=1,3
      TCF(IB)=EXP(TT8(MDLOHI(IP),IB,IT,IP)*TDIFF(IP))
   85 CONTINUE
C
      RETURN
      END
