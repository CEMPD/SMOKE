C Version @(#)$Id$ $Source$ $Date$ 

      SUBROUTINE STOERQ(INERR)
C
C  STOERQ attempts to store the user-supplied specifications for new FTP basic
C  emission rate equations and when successful (i.e., when room found in the
C  user supplied rate arrays), zeroes out all idle coefficients for the
C  specified pollutant, vehicle and region combination.
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
C    /FLAGS4/ IDLFLG
C    /IDLEQ2/ MAXIDL
C
C  Output on return:
C
C    parameter list: INERR
C    common blocks:
C    /BASEQ4/ NEWFIT
C    /BASEQ5/ ERUZML,ERUDR,ERU50K
C    /BASEQ6/ MYGERU,NUMERU
C    /BASEQ7/ IPMOD
C    /IDLEQ1/ ZMLIDL,DRIDL
C
C  Notes:
C
C  None.
C
C
      INTEGER PRTFLG,HCFLAG,COLDFG
C
      COMMON /BASEQ4/NEWPAR(5,100),BERNEW(3,100),NEWFIT(100),NEWCT,MAXCT
      COMMON /BASEQ5/ ERUZML(18,3,8,2),ERUDR(18,3,8,2),ERU50K(18,3,3,2)
      COMMON /BASEQ6/ MYGERU(12,2,3,8,2),MAXERU,NUMERU(3,8,2),KEYER,IGER
      COMMON /BASEQ7/ IPMOD(3)
      COMMON /FLAGS4/ PRTFLG,IDLFLG,NMHFLG,HCFLAG,COLDFG
      COMMON /IDLEQ1/ ZMLIDL(15,3,8,2),DRIDL(15,3,8,2)
      COMMON /IDLEQ2/ MYGIDL(15,3,8,2),MAXIDL
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
C  Set flag to skip VIDLE calculation for IP affected by accepted new equation.
C  Warn user and then zero out all idle coefficients for the accepted new
C  (IP,IV,IR) equation case.  Skip if idle ef not requested or already zeroed
C  out this case.
C
      IF(IDLFLG.EQ.1.OR.NUMERU(IP,IV,IR).GT.1) GOTO 20
      IF(IPMOD(IP).EQ.0) IPMOD(IP)=1
      CALL QUITER(0.,IFAC,59,INERR)
      DO 10 IEQU=1,MAXIDL
      ZMLIDL(IEQU,IP,IV,IR)=0.0
      DRIDL(IEQU,IP,IV,IR)=0.0
   10 CONTINUE
C
   20 CONTINUE
C
      RETURN
      END
