C Version "@(#)$Id$ $Source$ $Date$ 

      SUBROUTINE STOERQ(INERR)
C
C  STOERQ attempts to store the user-supplied specifications for new FTP basic
C  emission rate equations and when successful (i.e., when room found in the
C  user supplied rate arrays).
C
C  Called by GETERQ.
C
C  Calls QUITER.
C
C  Input on call:
C
C    parameter list: INERR
C    common blocks:
C    /BASEQ4/ NEWPAR,BERNEW,NEWCT
C    /BASEQ6/ MAXERU
C
C  Output on return:
C
C    parameter list: INERR
C    common blocks:
C    /BASEQ4/ NEWFIT
C    /BASEQ5/ ERUZML,ERUDR,ERU50K
C    /BASEQ6/ MYGERU,NUMERU
C
C  Notes:
C
C  None.
C
      INCLUDE 'BASEQ4.I'
      INCLUDE 'BASEQ5.I'
      INCLUDE 'BASEQ6.I'
C
      DO 20 IFAC=1,NEWCT
      NEWFIT(IFAC)=0
      IR=NEWPAR(1,IFAC)
      IV=NEWPAR(2,IFAC)
      IP=NEWPAR(3,IFAC)
C
C  Check if there's room to store equation IFAC's parameters.
C
      IERNEW=NUMERU(IP,IV,IR)+1
      IF(IERNEW.GT.MAXERU) CALL QUITER(0.,IFAC,51,INERR)
      IF(IERNEW.GT.MAXERU) GOTO 20
C
C  There is room.
C
      NUMERU(IP,IV,IR)=IERNEW
      MYGERU(IERNEW,1,IP,IV,IR)=NEWPAR(4,IFAC)
      MYGERU(IERNEW,2,IP,IV,IR)=NEWPAR(5,IFAC)
      ERUZML(IERNEW,IP,IV,IR)=BERNEW(1,IFAC)
      ERUDR(IERNEW,IP,IV,IR)=BERNEW(2,IFAC)
      IF(BERNEW(3,IFAC).GE.0.0) ERU50K(IERNEW,IP,IV,IR)=BERNEW(3,IFAC)
      NEWFIT(IFAC)=1
C
   20 CONTINUE
C
      RETURN
      END
