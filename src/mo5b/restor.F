C Version "@(#)$Id$ $Source$ $Date$ 

      SUBROUTINE RESTOR
C
C  RESTOR restores the BLOCK DATA "default" values to the common block
C  variables and arrays that can be modified through user one-time data option
C  selections.  In subroutine MOBILE4 these must be restored, since the user
C  may elect different options on successive calls.  Hence, RESTOR is called
C  at the outset of each call to MOBILE after the first.
C
C  Called by MOBILE.
C
C  Calls none.
C
C  Input on call:
C
C     common blocks:
C     /SAVE04/ VA007,VA008
C     /SAVE05/ VA009,VA010
C     /SAVE06/ VA011,VA012,IVA013,VA014
C
C  Output on return:
C
C     common blocks:
C     /ETPPAR/ ETP,ETPPIV,NETPYR
C     /IMPAR4/ CUTHC,CUTCO,CUTNO
C     /MAXIMA/ MAXVEH,MAXPOL,MAXREG,MAXYRS
C     /MYRSAV/ AMAR,JULMYR
C     /TAMEQ1/ TAMZML,TAMDR,TAMA50
C     /VMXCOM/ VMTMIX
C
C  Local variable / array dictionary:
C
C   Name   Type              Description
C  ------  ----  -------------------------------------------------------
C  I<d>     I    array subscripts for the save areas
C
C  Notes:
C
C  RESTOR was added for Version 01's conversion of MOBILE4 to a subroutine.
C
C  In Version 03, "IF(I2.LE.2)" removed before TAMZML line (dim now 3).
C  TAMA50 is left alone, since for the release version, users will not be
C  given the option to enter their own TAMA50.
C
C 
      INCLUDE 'SAVE01.I'
      INCLUDE 'SAVE04.I'
      INCLUDE 'SAVE05.I'
      INCLUDE 'SAVE06.I'
C
      INCLUDE 'ETPPAR.I'
      INCLUDE 'IMPAR4.I'
      INCLUDE 'MAXIMA.I'
      INCLUDE 'MYRSAV.I'
      INCLUDE 'TAMEQ1.I'
      INCLUDE 'VMXCOM.I'
C
C  Subroutine GETMRN
C
      DO 40 I1=1,MAXVEH
      DO 40 I2=1,MAXYRS
      AMAR(I2,I1)=VA007(I2,I1)
      JULMYR(I2,I1)=VA008(I2,I1)
  40  CONTINUE
C
C  Subroutine GETTAM
C
      DO 50 I1=1,2
      DO 50 I2=1,3
      DO 50 I3=1,4
      DO 50 I4=1,9
      TAMZML(I4,I3,I2,I1)=VA009(I4,I3,I2,I1)
      TAMDR(I4,I3,I2,I1)=VA010(I4,I3,I2,I1)
  50  CONTINUE
C
C  Subroutine GETVMT
C
      DO 60 I1=1,MAXVEH
      VMTMIX(I1)=VA011(I1)
  60  CONTINUE
C
C  Subroutine GETETP
C
      DO 70 I1=1,3
      DO 70 I2=1,7
      ETP(I2,I1)=VA012(I2,I1)
  70  CONTINUE
C
      DO 80 I1=1,NETPYR
      IETPYR(I1)=IVA013(I1)
      ETPPIV(I1)=VA014(I1)
  80  CONTINUE
C
C  Subroutine GETIMC
C
      CUTHC(1)=220.0
      CUTHC(2)=220.0
      CUTCO(1)=1.2
      CUTCO(2)=1.2
      CUTNO(1)=999.0
      CUTNO(2)=999.0
      RETURN
      END
