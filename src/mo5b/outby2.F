C Version "@(#)$Id$ $Source$ $Date$ 

      SUBROUTINE OUTBY2(ICY)
C
C  OUTBY2 writes scenario emission factor results by model year for selected
C  vehicle types on unit IOUREP using output format (OUTFMT) 5, when only CO
C  emission factors are requested (PRTFLG=2).
C
C  Called by OUTDT5.
C
C  Calls OBYIMP, and MYRIMT.
C
C  Input on call:
C
C    parameter list: ICY
C    common blocks:
C    /BYMYC1/ IVPICK,IMPICK
C    /IOUCOM/ IOUREP
C    /MAXIMA/ MAXVEH,MAXYRS
C
C  Local variable / array dictionary:
C
C   Name   Type              Description
C  ------  ----  -------------------------------------------------------
C  IMVEH    I    Flag to indicate the presence of an I/M program for that
C                particular IVTAM 1=none 2=present
C  IP,IV,JDX,MY
C           I    The usual definitions of these variables applies.
C
C
C  Notes:
C
C  (1) OUTBY2 was added for MOBILE4 In-House Version 02 (By Model Year Tables).
C  OUTBY1-4 were modified in MOBILE5.3 to make use of the MYRIMT function.
C  OUTBY2 is similiar to OUTBY3, they both had a JDX loop added to calculate
C  IMVEH.
C
C
      INCLUDE 'BYMYC1.I'
      INCLUDE 'BYMYC4.I'
      INCLUDE 'IOUCOM.I'
      INCLUDE 'MAXIMA.I'
C
      IP=2
      DO 20 IV=1,MAXVEH
        IF(IVPICK(IV).EQ.1) GOTO 20
           IMVEH=0
           DO 10 JDX=1,MAXYRS
             MY=ICY-JDX+1
             IF(IMVEH.EQ.0) IMVEH=MYRIMT(MY,IV,ICY)
   10      CONTINUE
C
           CALL OBYHDR(IOUREP)
C
C          Always print the base emission factor components table.
C
           CALL OBYTCN(ICY,IP,IV)
C
C          If By Model Year I/M Table picked, print it.
C
           IF(IV.GT.4) GOTO 20
              IF(IMPICK.EQ.2.AND.IMVEH.NE.0) CALL OBYIMP(ICY,IV)
C
   20 CONTINUE
C
      RETURN
      END
