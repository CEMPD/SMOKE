C Version "@(#)$Id$ $Source$ $Date$ 

      SUBROUTINE CEVTPD
C
C  CEVTPD calculates trips per day by the type and relative (to ICY) age of
C  the vehicle.  The formula uses the intercept and slope data of the stored
C  tpd on type/age regression model.  Since relative age is used, one set
C  of rates is calculated for use on all scenarios of the run.
C
C  Called by ONESEC.
C
C  Input on call:
C
C    common blocks:
C    /MAXIMA/ MAXVEH,MAXYRS
C    /TPDCOM/ TPDCO
C
C  Output on return:
C
C    common blocks:
C    /CEVBMY/ BMYTPD
C
C  Notes:
C
C  None.
C
      INCLUDE 'CEVBMY.I'
      INCLUDE 'MAXIMA.I'
      INCLUDE 'TPDCOM.I'
C
      DO 10 IV=1,MAXVEH
      DO 10 JDX=1,MAXYRS
      BMYTPD(JDX,IV)=TPDCO(1,IV)-(TPDCO(2,IV)*JDX)
   10 CONTINUE
C
      RETURN
      END
