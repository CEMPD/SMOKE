C Version @(#)$Id$ $Source$ $Date$ 

      SUBROUTINE HDDMYM(ICY)
C
C  HDDMYM computes annual weighted calendar year mileage for HDDV.
C
C  Called by REGMOD.
C
C  Input on call:
C
C    parameter list: ICY
C    common blocks:
C    /HDDMAR/ LHDVM,MHDVM,HHDVM
C    /HDDREG/ LHDRG,MHDRG,HHDRG
C    /MAXIMA/ MAXYRS
C
C  Output on return:
C
C    common block: /MYRSAV/ AMAR,NEWCUM
C
C  Notes:
C
C  None.
C
C
      REAL JULMYR
C
      COMMON /HDDMAR/ C2BDVM(25),LHDVM(25),MHDVM(25),HHDVM(25)
      COMMON /HDDREG/ C2BDRG(21),LHDRG(21),MHDRG(21),HHDRG(21)
      COMMON /MAXIMA/ MAXVEH,MAXLTW,MAXPOL,MAXREG,MAXYRS
      COMMON /MYRSAV/ AMAR(25,8),JULMYR(25,8),NEWCUM
C
      REAL LHDVM,MHDVM,HHDVM,LHDRG,MHDRG,HHDRG
C
C  Determine pointer KDX.  KDX=1 for 1980-, 2 for 1981,
C  3 for 1982, ..., 21 for 2000+.
C
      KDX=ICY-1979
      IF(KDX.LT.1) KDX=1
      IF(KDX.GT.21) KDX=21
C
C  Compute registration weighted annual mileage distribution for HDDV.
C
      DO 10 JDX=1,MAXYRS
      AMAR(JDX,7)=(C2BDRG(KDX)*C2BDVM(JDX)+
     *   LHDRG(KDX)*LHDVM(JDX)+MHDRG(KDX)*MHDVM(JDX)+
     *   HHDRG(KDX)*HHDVM(JDX))/(C2BDRG(KDX)+LHDRG(KDX)+
     *   MHDRG(KDX)+HHDRG(KDX))
   10 CONTINUE
C
C  Set flag so that model year mileage is recalculated in EFCALX.
C
      NEWCUM=1
      RETURN
      END
