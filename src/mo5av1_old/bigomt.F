C Version @(#)$Id$ $Source$ $Date$ 

      SUBROUTINE BIGOMT(IP,IV,IVLM,MY,IDX)
C
C  Compute bef's multiplicative correction for operating mode and temperature.
C
C  Also computes OMTCFF.
C  Called by BIGCFX.
C
C  Calls IBFPTR.
C
C  Input on call:
C
C    parameter list: IP,IV,IVLM,MY,IDX
C    common blocks:
C    /OPMOD1/ BFRCOF,BFSTEP,BFR50K
C    /OPMOD2/ MYSTEP
C    /OPMOD3/ BFMILE,KINKBF,BFGT50
C    /TAMEQ2/ F50K
C    /TEMPC3/ CHFRAC,TCF
C    /TEMPC4/ LOWCO
C
C  Output on return:
C
C    common blocks:
C    /OMTCOM/ OMTCF, OMTCFF
C
C  Local variable / array dictionary:
C
C   Name   Type              Description
C  ------  ----  -------------------------------------------------------
C  COMT     R    working location for operating mode and temperature c.f. to
C  DBF      R    BF deterioration rate for the first 50K miles
C  D50      R    BF deterioration rate for beyond 50K miles
C  FTPBF    R    FTP BF normalizing factor calculated at JDX x IV mileage
C                be multiplied times the base emission factor
C  OMT      R    unweighted product of operating mode vmt fraction times
C                temperature correction factor
C  ZBF      R    BF zero mile level
C
C  Notes:
C
C  Added OMTCFF [ OMTCF @ FTP ] to correct OMTCF in Subroutine
C  BEFIDL for MOBILE4.1.
C  BIGOMT was modiefied in MOBILE5 to include kinks for LDGT1/2.
C
      COMMON /OMTCOM/ OMTCF(25,3,8),OMTTAM(25,3,4),OMTCFF(25,3,8)
      COMMON /OPMOD1/ BFRCOF(8,8,3,6),BFSTEP(8,12,3,3),BFR50K(4,12,3,3)
      COMMON /OPMOD2/ MYSTEP,MAXBFR,MYGBFR(8,3,6)
      COMMON /OPMOD3/ BFMILE,KINKBF,BFGT50
      COMMON /TAMEQ2/ MYGTAM(3,4),F50K
      COMMON /TEMPC3/ MDLOHI(3),TDIFF(3),CHFRAC(3),TCF(3)
      COMMON /TEMPC4/ MYGCOO(3),LOWCO,COOTFC(3),COOTFD(3)
      COMMON /TEMPS/ AMBT,TEMMIN,TEMMAX,TEMEXH(3),TEMEVP(6)
C
C  HDV: a single RVP & T cf - assume 100% hot stabilized opmode = TCF's bag 2
C
      IF(IV.EQ.4.OR.IV.EQ.7) OMTCF(IDX,IP,IV)=TCF(2)
C
C  fill in the not used IV's
      IF(IV.EQ.4.OR.IV.EQ.7) OMTCFF(IDX,IP,IV)=1.0
C
      IF(IV.EQ.4.OR.IV.EQ.7) RETURN
C
C  LDV/T and MC: bag loop generates the operating mode and temperature cf
C  using the appropriate RIPSTWXN bag term coefficients.
C
C  Note: multiplicative tcf for correcting tampering offset is neutral (1.0)
C  when LOWCO = 2 ( = cold start LDGV/T 1980/4+ CO -> use additive CO offset).
C  So if cold start & low temp CO case, correct only for operating mode & RVP.
C
      COMT=0.0
      IGLM=IBFPTR(MY,IP,IVLM)
      DO 40 IB=1,3
      IF(IB.EQ.2) COMT1=0.0
      OMT=CHFRAC(IB)
      OMT=OMT*TCF(IB)
      IF(MYSTEP.EQ.2) GOTO 10
C
C  Not 1981+ LDGV and not 1984+ LDGT1/2 - access "old"
C  (from MOBILE2) array using myg index lookup.
C
      ZBF=BFRCOF(2*IB-1,IGLM,IP,IVLM)
      DBF=BFRCOF(2*IB,  IGLM,IP,IVLM)
      GOTO 20
C
C  1981+ LDGV or 1984+ LDGT1/2 - access its own
C  bf array using calculated my index.
C
   10 ZBF=BFSTEP(2*IB-1,IGLM,IP,IVLM)
      DBF=BFSTEP(2*IB,  IGLM,IP,IVLM)
      IF(KINKBF.EQ.2) GOTO 30
C
   20 COMT=COMT+OMT*(ZBF+DBF*BFMILE)
C   COMT1 used to calculate OMTCFF
      IF(IB.EQ.2) COMT1=(ZBF+DBF*BFMILE)
      GOTO 40
C
C  1981+ LDGV and 1984+ LDGT1/2 with mileage > 50K - get >50K
C  det rate using same MY index.
C
   30 D50=BFR50K(IB,IGLM,IP,IVLM)
      COMT=COMT+OMT*(ZBF+DBF*F50K+D50*BFGT50)
C  COMT1 used to calculate OMTCFF
      IF(IB.EQ.2) COMT1=(ZBF+DBF*F50K+D50*BFGT50)
C
   40 CONTINUE
C
C  Normalize by dividing through by the formula evaluated at FTP coefficients.
C
C  Case 1: not 1981+ LDGV and not 1984+ LDGT1/2
C
      IF(MYSTEP.EQ.1) FTPBF=
     *  BFRCOF(7,IGLM,IP,IVLM)+(BFRCOF(8,IGLM,IP,IVLM)*BFMILE)
      IF(MYSTEP.EQ.1) GOTO 50
C
C  Case 2: Subcase 1: 1981+ LDGV and 1984+ LDGT1/2 <=50K miles (KINKBF=1).
C          Subcase 2: 1981+ LDGV and 1984+ LDGT1/2 > 50K miles. (KINKBF=2)
C
      IF(KINKBF.EQ.1) FTPBF=BFSTEP(7,IGLM,IP,IVLM)+
     *  (BFSTEP(8,IGLM,IP,IVLM)*BFMILE)
      IF(KINKBF.EQ.2) FTPBF=BFSTEP(7,IGLM,IP,IVLM)+
     *  (BFSTEP(8,IGLM,IP,IVLM)*F50K)+(BFR50K(4,IGLM,IP,IVLM)*BFGT50)
C
C  Stock Operating Mode Temperature Correction Factor
   50 OMTCF(IDX,IP,IV)=COMT/FTPBF
C
C  OMTCFF is OMTCF @ FTP [ 100% stable Mode /0/0/0/, 75 deg F
C  and 9.0 RVP ] for correcting BEFIDL in Subroutine BEFIDL
C  OMTCFF uses Bag TWO(2) ZBF and DBF each of the Three(3) IP's
C
      OMTCFF(IDX,IP,IV)=COMT1/FTPBF
C
      RETURN
      END
