C Version "@(#)$Id$ $Source$ $Date$ 

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
C    /LEVBLK/ LEVYRS,APPLEV
C    /YEARS4/ IYEND
C
C  Notes:
C
C   OUTZEV was added for MOBILE5v3  (By Model Year Tables).
C
C October-5-1994 @CSC pme request 2-472. A check for APPLEV variable to 
C see whether or not to apply LEV stds to LDGT2 was included.
C
C
      INCLUDE 'BASE10.I'
      INCLUDE 'IOUCOM.I'
      INCLUDE 'LEVBLK.I'
      INCLUDE 'MAXIMA.I'
      INCLUDE 'YEARS4.I'
C
      CHARACTER*4 DIVCHA(25)
C
      DIMENSION TTF(25)
      DIMENSION IYR(25)
C
      DATA DIVCHA /
     *  '____', '____','____','____','____','____','____',
     *  '____', '____','____','____','____','____','____',
     *  '____', '____','____','____','____','____','____',
     *  '____', '____','____','____' /
C
C     ZEV start in the fleet in 1998, but ZEV's run for 57 years
C     from 1994 to 2050
C
      IF(ICY.LT.1994.OR.IV.EQ.3.AND.APPLEV.EQ.2) GOTO 99
C
C       find the end of array for writes
C
        IEND=(ICY-LEVYRS)+1
        IF(IEND.GT.25) IEND=25
C
C       fill in IYR array - icy's from 1994 to 2050
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
C      flip the values in the TF array for printing
C
           DO 15 I=1, MAXYRS
             J=(MAXYRS+1)-I
             TTF(J)=ZEVTF(I,IV)
  15       CONTINUE
C
          WRITE(IOUREP,125) LEVYRS,IYEND-1
  125     FORMAT(' ','Travel Fraction attributable to Zero Emitting ',
     *    'Vehicles (ZEV)  Model Years ',I4,' thru ',I4,/,
     *    ' ','             ' )
C
        IF(IEND.LE.20) THEN
C
C          this format is for LE 20 years
            WRITE(IOUREP,130) (IYR(I),I=1,IEND)
  130       FORMAT(' ',' MY ',1X,23(I4,2X) )
C
            WRITE(IOUREP,135) (DIVCHA(I),I=1,IEND)
  135       FORMAT('+','____',1X,23(A4,2X) )
C
            WRITE(IOUREP,140) (TTF(J),J=1,IEND)
 140        FORMAT(' ',' TF ',23(F5.4,1X),/,
     *             ' ','                              ')
          ELSE
C
C       this format is for the full 25 years of MAXYRS
          WRITE(IOUREP,230) (IYR(I),I=1,IEND)
  230     FORMAT(' ',' MY ',1X,25(I4,1X) )
C
          WRITE(IOUREP,235) (DIVCHA(I),I=1,IEND)
  235     FORMAT('+','____',1X,25(A4,1X) )
C
          WRITE(IOUREP,240) (TTF(J),J=1,IEND)
  240     FORMAT(' ',' TF ',25(F5.4),/,
     *           ' ','                              ')
       ENDIF
C
 99   CONTINUE
      RETURN
      END
