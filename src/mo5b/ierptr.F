C Version "@(#)$Id$ $Source$ $Date$ 

      FUNCTION IERPTR(MY,IP,IV,KEYER)
C
C  IERPTR selects the array pair and the model year group pointers into the
C  intercept and slope arrays for the FTP emission rate equations.
C
C  Called by BEF.
C
C  Calls QUITER.
C
C  Input on call:
C
C    parameter list: MY,IP,IV
C    common blocks:
C    /BASEQ2/ MYGERB,MAXERB
C    /BASEQ6/ MYGERU,NUMERU
C    /FLAGS2/ NEWFLG
C    /REGION/ IREJN
C
C  Output on return:
C    parameter list: KEYER
C    function: IERPTR
C
C  Local variable / array dictionary:
C
C   Name   Type              Description
C  ------  ----  -------------------------------------------------------
C  NERU     I    number of user supplied ER equations for the IPxIVxIREJN
C                case passed in on the call to IERPTR
C
C  Notes:
C
C  STOREQ looks for the match in first the user supplied (if any) and then the
C  hardcoded (initialized in BLOCK DATA) equation myg arrays.
C
C  The user is responsible for insuring entered ER records do not have
C  overlapping mygs.  IERPTR returns the index to the first myg match.
C
C
      INCLUDE 'BASEQ2.I'
      INCLUDE 'BASEQ6.I'
      INCLUDE 'FLAGS2.I'
      INCLUDE 'REGION.I'
C
      NERU=NUMERU(IP,IV,IREJN)
C
C     use mobile5 default BER's
      IF(NEWFLG.EQ.1.OR.NERU.EQ.0) GOTO 20
C
C
C  Might be a match in the user supplied myg's.
C
      KEYER=2
      DO 10 IG=1,NERU
         IF(MY.GE.MYGERU(IG,1,IP,IV,IREJN).AND.
     *      MY.LE.MYGERU(IG,2,IP,IV,IREJN)) GOTO 90
   10 CONTINUE
C
C  Look in the Mobile5 Block Data hardcoded myg's.
C
   20 KEYER=1
      DO 30 IG=1,MAXERB
         IF(MY.LE.MYGERB(IG,IP,IV,IREJN)) GOTO 90
   30 CONTINUE
C
C  Set default.
C
      IG=1
      CALL QUITER(0.,11,61,MY)
C
   90 IERPTR=IG
C
      RETURN
      END
