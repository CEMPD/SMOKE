C Version "@(#)$Id$ $Source$ $Date$ 

      FUNCTION MYRATP(MY,IV,ICY)
C
C  MYRATP discovers the controlling I/M Program number 1 or 2.
C  (none=0) for a given MY, Vehicle, and CY for I/M
C
C  Called by TAMPER
C
C  Input on call:
C
C    Parameter list: MY,IV,ICY
C    common blocks:
C    /FLAGS2/ IMFLAG
C    /IMPAR1/ ICYIM,MODYR
C    /IMPAR2/ ILDT
C
C  Output on return:
C
C    Function: MYATP
C
C  Local variable dictionary:
C
C   Name   Type              Description
C  ------  ----  ---------------------------------------------------
C  IMNOT1   I    Flag to indicate the I/M program number one covers
C                this time the function is called 1=not covered
C                2=covered
C  IMNOT2   I    Flag to indicate the I/M program number two covers
C                this time the function is called 1=not covered
C                2=covered
C  LOCIND   I    Local variable used to hold function return value
C
C
C  Notes:
C
C  This function was developed for recalculating the I/M deterrence
C  effect, April 13, 1993 @ ARC-bsg Subtask 2-212
C
C
      INCLUDE 'FLAGS2.I'
      INCLUDE 'IMPAR1.I'
      INCLUDE 'IMPAR2.I'
C
      LOCIND = 0
      IMNOT1 = 1
      IMNOT2 = 1
C
C  The following are possible exclusions:
C  no I/M program, vehicle class or model year
C  is not covered or scenario is before I/M start year.
C
      IF(IV.GT.4) GOTO 98
C
C  There are two I/M descriptions to check and compare
C
      IF(IMFLAG.EQ.3) THEN
        IF((ILDT(IV,1).EQ.1).OR.(ICY.LE.ICYIM(1)).OR.
     *    (MY.GT.MODYR(2,1)).OR.(MY.LT.MODYR(1,1))) THEN
            IMNOT1 = 1
        ELSE
            IMNOT1 = 2
        ENDIF
        IF((ILDT(IV,2).EQ.1).OR.(ICY.LE.ICYIM(2)).OR.
     *    (MY.GT.MODYR(2,2)).OR.(MY.LT.MODYR(1,2))) THEN
            IMNOT2 = 1
        ELSE
            IMNOT2 = 2
        ENDIF
C
C  No first I/M program, second exist set to "2"
C
        IF(IMNOT1.EQ.1.AND.IMNOT2.EQ.2) LOCIND=2
C
C  No second I/M program, first exist set to "1"
C
        IF(IMNOT1.EQ.2.AND.IMNOT2.EQ.1) LOCIND=1
C
C  First I/M program and second exists but choose the I/M
C  program which started first.
C
        IF(IMNOT1.EQ.2.AND.IMNOT2.EQ.2) THEN
          IF(ICYIM(1).LT.ICYIM(2)) THEN
            LOCIND = 1
          ELSE
            LOCIND = 2
          ENDIF
        ENDIF
      ENDIF
C
C  There is only one I/M description to check
C
      IF(IMFLAG.EQ.2) THEN
            IF((ILDT(IV,1).EQ.1).OR.(ICY.LE.ICYIM(1)).OR.
     *      (MY.GT.MODYR(2,1)).OR.(MY.LT.MODYR(1,1))) THEN
              LOCIND = 0
            ELSE
              LOCIND = 1
            ENDIF
       ENDIF
C
C  I/M Program #2 is controlling if there is an I/M program 1 and
C  an I/M program 2 in conflict then I/M program 2 wins unless
C  the first I/M program has an earlier start date.
C
   98 MYRATP=LOCIND
      RETURN
      END
