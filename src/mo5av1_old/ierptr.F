C Version @(#)$Id$ $Source$ $Date$ 

      SUBROUTINE IERPTR(MY,IP,IV,IDLGER)
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
C    /BASEQ8/ HDVMY
C    /FLAGS2/ NEWFLG
C    /REGION/ IREJN
C
C  Output on return:
C    parameter list: IDLGER
C    common block: /BASEQ6/ KEYER,IGER
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
      INTEGER ALHFLG
C
      COMMON /BASEQ2/ MYGERB(25,3,8,2),MAXERB
      COMMON /BASEQ6/ MYGERU(12,2,3,8,2),MAXERU,NUMERU(3,8,2),KEYER,IGER
      COMMON /FLAGS2/ MYMRFG,NEWFLG,IMFLAG,ALHFLG,IMTIER
      COMMON /REGION/ FEET(2),IREJN,ALT,INITPR
      COMMON /BASEQ8/ IHDVMY(23,2),HDVZML(23,2),HDVBDR(23,2)
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
      DO 10 IGER=1,NERU
         IF(MY.GE.MYGERU(IGER,1,IP,IV,IREJN).AND.
     *      MY.LE.MYGERU(IGER,2,IP,IV,IREJN)) GOTO 99
   10 CONTINUE
C
C  Look in the Mobile5 Block Data hardcoded myg's.
C
   20 KEYER=1
      DO 30 IGER=1,MAXERB
         IF(MY.LE.MYGERB(IGER,IP,IV,IREJN)) GOTO 40
   30 CONTINUE
C
C  Set default.
C
      IGER=1
      CALL QUITER(0.,11,61,MY)
      GOTO 99
C
C  Idle Pointer
C
   40 IDLGER=IGER
C
      IF(IV.EQ.4.AND.IP.EQ.2) THEN
         DO 50 J=1,MAXERB
            IF(MY.LE.IHDVMY(J,IREJN)) THEN
               IDLGER=J
               GOTO 99
            ENDIF
   50    CONTINUE
      ENDIF
C
   99 RETURN
C
      END
