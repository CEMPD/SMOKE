C Version "@(#)$Id$ $Source$ $Date$ 

      FUNCTION OVLP2(MYR,IVTAM,ICY)
C
C  OVLP2 discovers the controlling I/M Program number 1 or 2
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
C    Function: OVLP2
C
C
C  Notes:
C
C  This function was developed for PCLEFT modifications made to
C  MOBILE5 RFP in order to calculate correctly the phase-in of
C  Biennial I/M credits.
C  September-16-1993 @ ARC-bsg Subtask 311 
C  OVLP2 is the second function developed to check for I/M program
C  overlap
C
C
      INTEGER LOVLP2,OVLP2
      LOGICAL IMPRG1,IMPRG2
C
      INCLUDE 'FLAGS2.I'
      INCLUDE 'IMPAR1.I'
      INCLUDE 'IMPAR2.I'
      INCLUDE 'IMPAR6.I'
C
      LOVLP2=0
C
C  The following are possible exclusions: no I/M program, vehicle
C  class is not covered or scenario is before I/M start year.
C  Returns:  0 - No program,  
C  ICYIM(1) - 1st program,
C  ICYIM(2) - 2nd program
C  Both programs apply so use earlier srart year
C
      IF(IVTAM.GT.4) GOTO 99
      IF(IMFLAG.EQ.3) THEN
C
C  Inititalize both, there are two I/M programs in overlap
C
        IMPRG1=.TRUE.
        IMPRG2=.TRUE.
C
C  Set the second I/M program to not an I/M program
C
        IF((ILDT(IVTAM,2).EQ.1).OR.(ICY.LE.ICYIM(2)).OR.
     *    (MYR.GT.MODYR(2,2)).OR.(MYR.LT.MODYR(1,2))) IMPRG2=.FALSE.
C
C  Set the first I/M program to not an I/M program
C
        IF((ILDT(IVTAM,1).EQ.1).OR.(ICY.LE.ICYIM(1)).OR.
     *    (MYR.GT.MODYR(2,1)).OR.(MYR.LT.MODYR(1,1))) IMPRG1=.FALSE.
C
C  No first I/M program, second I/M program use second start year
C
        IF(.NOT.IMPRG1.AND.IMPRG2) THEN
             LOVLP2=ICYIM(2)
C
C  No second I/M program, first I/M program use first start year
C
          ELSE IF(IMPRG1.AND..NOT.IMPRG2) THEN
              LOVLP2=ICYIM(1)
C
C  Both so choose earlier start year
C
          ELSE IF(ICYIM(1).LE.ICYIM(2)) THEN 
              LOVLP2=ICYIM(1)
          ELSE 
              LOVLP2=ICYIM(2)
        ENDIF
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
     *    (MYR.GT.MODYR(2,1)).OR.(MYR.LT.MODYR(1,1))) GOTO 99
        LOVLP2=ICYIM(1)
      ENDIF
   99 OVLP2=LOVLP2
      RETURN
C
      END
