C Version @(#)$Id$ $Source$ $Date$ 

      SUBROUTINE CEVMPD
C
C  CEVMPD calculates miles per day by the type and relative (to ICY) age of
C  the vehicle.  The formula uses the REGMOD computed TFMYM (January 1
C  mileage accumulation rates).  As with TYMYM and BMYTPD, one set
C  of rates is calculated for use on all scenarios of the run.
C
C  Called by REGMOD.
C
C  Input on call:
C
C    common blocks:
C    /MAXIMA/ MAXVEH,MAXYRS
C    /MYRCAL/ TFMYM
C
C  Output on return:
C
C    common blocks:
C    /CEVBMY/ BMYMPD
C
C  Notes:
C
C  None.
C
C
      REAL JANMYR
C
      COMMON /CEVBMY/ BMYMPD(25,8),BMYTPD(25,8),BMTPDF(25,8)
      COMMON /MAXIMA/ MAXVEH,MAXLTW,MAXPOL,MAXREG,MAXYRS
      COMMON /MYRCAL/ XMYM(25,8),JANMYR(25,8),TF(25,8),TFMYM(25,8)
      COMMON /TPDCOM/ TPDCO(2,8)
C
      DO 10 IV=1,MAXVEH
      DO 10 JDX=1,MAXYRS
      BMYMPD(JDX,IV)=(TFMYM(JDX,IV)*100000)/365.
   10 CONTINUE
C
      RETURN
      END
