C Version @(#)$Id$ $Source$ $Date$ 

      SUBROUTINE OUTZEV(IV,ICY)
C
C  OUTZEV writes the Zero Emission Vehicle (ZEV) By Model Year
C  Travel fraction TF on unit IOUREP
C
C
C  Called by OUTBY1, OUTBY2, OUTBY3, and OUTBY4.
C
C  Input on call:
C
C    parameter list: IV,ICY
C    common blocks:
C    /BASE10/ ZEVTF(25,3), ZEVTOT(25)
C    /IOUCOM/ IOUREP
C    /LEVBLK/ LEVYRS
C
C  Notes:
C
C   OUTZEV was added for MOBILE5v3  (By Model Year Tables).
C
C
      COMMON /BASE10/ ZEVTF(25,3),ZEVTOT(3)
      COMMON /IOUCOM/ IOUIMD,IOUGEN,IOUREP,IOUERR,IOUASK
      COMMON /LEVBLK/ LEVFLG,IMFLG,LEVYRS
      COMMON /MAXIMA/ MAXVEH,MAXLTW,MAXPOL,MAXREG,MAXYRS
C
      CHARACTER*4 DIVCHA(25)
C
      DIMENSION TTF(25)
      DIMENSION IYR(25)
C
C     ZEV start in the fleet in 1998, but ZEV's run for 27 years
C     from 1994 to 2020
C
      IF(ICY.LT.1994) GOTO 99
C
C
      IF(ICY.GE.1994) THEN
C
        DATA DIVCHA /
     *  '____', '____','____','____','____','____','____',
     *  '____', '____','____','____','____','____','____',
     *  '____', '____','____','____','____','____','____',
     *  '____', '____','____','____' /
C
C       find the end of array for writes
        IEND=(ICY-LEVYRS)+1
        IF(IEND.GT.25) IEND=25
C
C       fill in IYR array - icy's from 1994 to 2020
        DO 10 I=1,MAXYRS
           IYR(I)=(ICY+1)-I
   10   CONTINUE
C
C       initilize TTF
C
        DO 11 I=1, MAXYRS
           TTF(I)=0.0
  11    CONTINUE
C
C      flip the values in the LDGV TF array for printing
C
        IF(IV.EQ.1) THEN
           DO 15 I=1, MAXYRS
             J=(MAXYRS+1)-I
             TTF(J)=ZEVTF(I,IV)
  15       CONTINUE
        ENDIF
C
C       flip the values in the LDGT1 TF array for printing
C
        IF(IV.EQ.2) THEN
           DO 20 I=1, MAXYRS
             J=(MAXYRS+1)-I
             TTF(J)=ZEVTF(I,IV)
  20       CONTINUE
        ENDIF
C
        WRITE(IOUREP,125) LEVYRS
  125   FORMAT(' ','Travel Fraction attributable to Zero Emitting '
     *,'Vehicles (ZEV)  Model Years ',I4,' thur 2020 ',/,
     *' ','             ' )
C
        IF(IEND.LE.20) THEN
C
C          this format is for LE 20 years
           WRITE(IOUREP,130) (IYR(I),I=1,IEND)
  130      FORMAT(' ',' MY ',1X,23(I4,2X) )
CC
           WRITE(IOUREP,135) (DIVCHA(I),I=1,IEND)
  135      FORMAT('+','____',1X,23(A4,2X) )
CC
           WRITE(IOUREP,140) (TTF(J),J=1,IEND)
 140       FORMAT(' ',' TF ',23(F5.4,1X),/,
     *' ','                              ')
           GOTO 99
CC
        ENDIF
C
C       this format is for the full 25 years of MAXYRS
        WRITE(IOUREP,230) (IYR(I),I=1,IEND)
  230   FORMAT(' ',' MY ',1X,25(I4,1X) )
CC
        WRITE(IOUREP,235) (DIVCHA(I),I=1,IEND)
  235   FORMAT('+','____',1X,25(A4,1X) )
CC
        WRITE(IOUREP,240) (TTF(J),J=1,IEND)
  240   FORMAT(' ',' TF ',25(F5.4),/,
     *' ','                              ')
      ENDIF
CC
 99   CONTINUE
      RETURN
      END
