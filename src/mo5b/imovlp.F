C Version "@(#)$Id$ $Source$ $Date$ 

      FUNCTION IMOVLP(MYR,IVTAM,ICY)
C
C  IMOVLP discovers the controlling I/M Program number 1 or 2
C  for a given MYR, Vehicle, and CY and checks for overlaps.
C
C  Called by PCLEFT
C
C  Input on call:
C
C    Parameter list: MYR,IVTAM,ICY
C    common blocks:
C    /FLAGS2/ IMFLAG
C    /IMPAR1/ ICYIM,MODYR
C    /IMPAR2/ ILDT
C    /IMPAR6/ IFREQ
C
C  Output on return:
C
C    Function: IMOVLP
C
C
C  Notes:
C
C  This function was developed for PCLEFT modifications made to
C  MOBILE5 RFP in order to calculate correctly the phase-in of
C  Biennial I/M credits.
C
C
      INCLUDE 'FLAGS2.I'
      INCLUDE 'IMPAR1.I'
      INCLUDE 'IMPAR2.I'
      INCLUDE 'IMPAR6.I'
C
      IMOVLP=0
C
C  The following are possible exclusions: no I/M program, vehicle
C  class is not covered or scenario is before I/M start year.
C  Returns:  0 - No program,  1 - 1st program,  2 - 2nd program,
C  3 - 2nd with overlap and 2nd program is biennial.
C
      IF(IVTAM.GT.4) RETURN
      IF(IMFLAG.EQ.3) THEN
C
C  Inititalize both, there are two I/M programs in overlap
C
        IMNOT1=1
        IMNOT2=1
C
C  Set the second I/M program to not an I/M program
C
        IF((ILDT(IVTAM,2).EQ.1).OR.(ICY.LE.ICYIM(2)).OR.
     *    (MYR.GT.MODYR(2,2)).OR.(MYR.LT.MODYR(1,2))) IMNOT2=2
C
C  Set the first I/M program to not an I/M program
C
        IF((ILDT(IVTAM,1).EQ.1).OR.(ICY.LE.ICYIM(1)).OR.
     *    (MYR.GT.MODYR(2,1)).OR.(MYR.LT.MODYR(1,1))) IMNOT1=2
C
C  Check for I/M overlap, set to overlap value "3"
C
        IF(IMNOT1.EQ.1.AND.IMNOT2.EQ.1) IMOVLP=3
C
C  Biennial second I/M program, check and reset "3" to "2"
C  if not biennial.
C
        IF(IFREQ(2).EQ.1.AND.IMOVLP.EQ.3) IMOVLP=2
C
C  No first I/M program, second I/M program exists set to "2"
C
        IF(IMNOT1.EQ.2.AND.IMNOT2.EQ.1) IMOVLP=2
C
C  No second I/M program, first I/M program exists set to "1"
C
        IF(IMNOT1.EQ.1.AND.IMNOT2.EQ.2) IMOVLP=1
C
C  endif for, if imflag=3
C
      ENDIF
C
C  If there is only one I/M descriptive record (IMFLAG=2)
C  determine if this model year/vehicle is currently
C  covered by that I/M program.
C
      IF(IMFLAG.EQ.2) THEN
        IF((ILDT(IVTAM,1).EQ.1).OR.(ICY.LE.ICYIM(1)).OR.
     *    (MYR.GT.MODYR(2,1)).OR.(MYR.LT.MODYR(1,1))) RETURN
        IMOVLP=1
      ENDIF
      RETURN
C
      END
