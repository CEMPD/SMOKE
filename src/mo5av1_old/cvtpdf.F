C Version @(#)$Id$ $Source$ $Date$ 

      SUBROUTINE CVTPDF
C
C  CVTPDF calculates the normalized trips per day fractions used in
C  composite hot soak calculations.
C
C  Called by REGMOD.
C
C  Input on call:
C
C    common blocks:
C    /CEVBMY/ BMYTPD
C    /MAXIMA/ MAXVEH,MAXYRS
C    /MYRCAL/ JANMYR
C
C  Output on return:
C
C    common blocks:
C    /CEVBMY/ BMTPDF
C
C  Local array subscripts
C
C  TPDNRM(8)  -  TPDNRM( IV )
C
C  Local variable / array dictionary:
C
C   Name   Type              Description
C  ------  ----  -------------------------------------------------------
C  TPDNRM   R    Summation variable for trips per day times registration
C                fraction.
C
C  Notes:
C
C  CVTPDF was added in MOBILE4.1 to calculate the normalized
C  trips per day fractions.
C
      REAL JANMYR
C
      COMMON /CEVBMY/ BMYMPD(25,8),BMYTPD(25,8),BMTPDF(25,8)
      COMMON /MAXIMA/ MAXVEH,MAXLTW,MAXPOL,MAXREG,MAXYRS
      COMMON /MYRCAL/ XMYM(25,8),JANMYR(25,8),TF(25,8),TFMYM(25,8)
C
      DIMENSION TPDNRM(8)
C
C  Compute the normalized trips per day fractions.
C
      DO 30 IV=1,MAXVEH
      TPDNRM(IV)=0.0
      DO 20 JDX=1,MAXYRS
      BMTPDF(JDX,IV)=BMYTPD(JDX,IV)*JANMYR(JDX,IV)
      TPDNRM(IV)=TPDNRM(IV)+BMTPDF(JDX,IV)
   20 CONTINUE
      IF(TPDNRM(IV).LE.0.0) TPDNRM(IV)=1.0
   30 CONTINUE
C
      DO 40 IV=1,MAXVEH
      DO 40 JDX=1,MAXYRS
      BMTPDF(JDX,IV)=BMTPDF(JDX,IV)/TPDNRM(IV)
   40 CONTINUE
C
      RETURN
      END
