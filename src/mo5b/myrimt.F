C Version "@(#)$Id$ $Source$ $Date$ 

      FUNCTION MYRIMT(MYR,IVTAM,ICY)
C
C  MYRIMT discovers the controlling I/M Program number 1 or 2. (none=0)
C  for a given MYR, Vehicle, and CY for I/M
C
C  Called by TAMPER, PCLEFT
C
C  Input on call:
C
C    Parameter list: MYR,IVTAM,ICY
C    common blocks:
C    /FLAGS2/ IMFLAG
C    /IMPAR1/ ICYIM,MODYR
C    /IMPAR2/ ILDT
C
C  Output on return:
C
C    Function: MYRIMT
C
C
C  Notes:
C
C  This function was developed for TAMPER modifications made to
C  MOBILE5.3 in order to calculate IM240 effects on catalyst removal
C
C
      INCLUDE 'FLAGS2.I'
      INCLUDE 'IMPAR1.I'
      INCLUDE 'IMPAR2.I'
C
      MYRIMT = 0
C
C  The following are possible exclusions no I/M program, vehicle class is
C  is not covered or scenario is before I/M start year.
C
      IF(IVTAM.GT.4) RETURN
      IF(IMFLAG.EQ.3) GOTO 60
      IF(IMFLAG.EQ.2) GOTO 40
      RETURN
   40 IF((ILDT(IVTAM,1).EQ.1).OR.(ICY.LE.ICYIM(1)).OR.
     *   (MYR.GT.MODYR(2,1)).OR.(MYR.LT.MODYR(1,1))) RETURN
         MYRIMT = 1
         RETURN
   60 IF((ILDT(IVTAM,2).EQ.1).OR.(ICY.LE.ICYIM(2)).OR.
     *   (MYR.GT.MODYR(2,2)).OR.(MYR.LT.MODYR(1,2))) GOTO 40
         MYRIMT = 2
C
C  I/M Program #2 is controlling if there is an I/M program 1 and
C  an I/M program 2 in conflict then I/M program 2 wins.
C
      RETURN
      END
