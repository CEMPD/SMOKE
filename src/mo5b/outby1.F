C Version "@(#)$Id$ $Source$ $Date$ 

      SUBROUTINE OUTBY1(ICY)
C
C  OUTBY1 writes scenario emission factor results by model year for selected
C  vehicle types on unit IOUREP using output format (OUTFMT) 5, when only HC
C  emission factors are requested (PRTFLG=1).
C
C  Called by OUTDT5.
C
C  Calls OBYHDR, MYRIMT, OBYIMP, and OUTZEV.
C
C  Input on call:
C
C    parameter list: ICY
C    common blocks:
C    /BYMYC1/ IVPICK,IMPICK,BYCUM
C    /BYMYC2/ BYBEF4,BYTAM,BYFER
C    /BYMYC3/ BYEVAP,BYRUNL,BYREFL,BYRSTL
C    /BYMYC4/ BYVTN,YFM,PN,U7,PSEC,VB
C    /EVAL/   MEVAL,MONTH
C    /IOUCOM/ IOUREP
C    /LEVBLK/ LEVFLG
C    /MAXIMA/ MAXVEH,MAXYRS
C    /MYRCAL/ TF
C    /RESUL1/ EFFTP,EFEVAP,EFREFL,EFRUNL
C    /RESUL2/ EFRSTL
C    /SPEED6/ SALHCF
C
C  Local variable / array dictionary:
C
C   Name   Type              Description
C  ------  ----  -------------------------------------------------------
C  IMVEH    I    Flag to indicate the presence of an I/M program for that
C                particular IVTAM 1=none 2=present
C  IP,IV,JDX,IDX,MY
C           I    The usual definitions of these variables applies.
C
C
C  Notes:
C
C  (1) OUTBY1 was added for MOBILE4 In-House Version 02 (By Model Year Tables).
C  OUTBY1-4 were modified in MOBILE5.3 to make use of the MYRIMT function.
C
C  (22-April-1994) @ CSC ljn Request 2-420 added LDGT2 (IV=3) to the LEV's
C
      INCLUDE 'BYMYC1.I'
      INCLUDE 'BYMYC2.I'
      INCLUDE 'BYMYC3.I'
      INCLUDE 'BYMYC4.I'
      INCLUDE 'EVAL.I'
      INCLUDE 'IOUCOM.I'
      INCLUDE 'LEVBLK.I'
      INCLUDE 'MAXIMA.I'
      INCLUDE 'MYRCAL.I'
      INCLUDE 'RESUL1.I'
      INCLUDE 'RESUL2.I'
      INCLUDE 'SPEED6.I'
C
      IP=1
      DO 90 IV=1,MAXVEH
      IF(IVPICK(IV).EQ.1) GOTO 90
C
      CALL OBYHDR(IOUREP)
C
C  Components printed and hence BYMY table format vary with the vehicle class.
C
      GOTO(10,10,10,10,30,30,30,60),IV
C
C  LDGV, LDGT1, LDGT2 & HDGV include tampering & all non-exhaust HC components.
C
C
   10 WRITE(IOUREP,209) BYVTN(IV), MONTH(MEVAL), ICY
  209 FORMAT(
     *'-',T54,A26,/,
     *' ',T59,A4,' 1 ', I4,/,
     *' ','                 ')
C
C  (22-April-1994) @ CSC ljn Request 2-420
C  write out the Zero Emitting Vehicles BYMY Travel Fraction
C
      IF(LEVFLG.EQ.2.AND.IV.LE.3) CALL OUTZEV(IV,ICY)
C
C    Finish writing
C
      WRITE(IOUREP,210) VB,PSEC(IP),YFM,VB,
     *    (PN(ICH),ICH=1,8),(U7,ICH=1,8)
  210 FORMAT(
     *'0',T2,'Model',T21,A1,T46,A3/
     *' ',T2,A17,T21,A1,1X,8A7/
     *'+',T2,'____',T7,'_____',T13,'_______',T22,8A7)
C
C
      IMVEH=0
      DO 20 JDX=1,MAXYRS
      IDX=MAXYRS-JDX+1
      MY=ICY-JDX+1
      IF(IMVEH.EQ.0) IMVEH=MYRIMT(MY,IV,ICY)
C
      WRITE(IOUREP,215)
     5    MY,TF(IDX,IV),BYCUM(JDX,IV),BYBEF4(IP,JDX,IV),
     5    BYTAM(IP,JDX,IV),SALHCF(IDX,IP,IV),BYEVAP(JDX,IV),
     5    BYREFL(JDX,IV),BYRUNL(JDX,IV),BYRSTL(JDX,IV),BYFER(IP,JDX,IV)
C
  215 FORMAT(' ',T2,I4,T7,F5.4,F8.0,T22,8F7.3)
C
   20 CONTINUE
C
C  Write column sums for TF weighted components, to check for identity with
C  the same sums written out in the OUTFMT = 3 table.
C
      WRITE(IOUREP,220) (U7,ICH=1,5),
     1      EFEVAP(IV),EFREFL(IV),EFRUNL(IV),EFRSTL(IV),EFFTP(IP,IV)
C
  220 FORMAT('+',T43,5A7/' ',T43,5F7.3)
C
C  If By Model Year I/M Table picked, print it.
C
      IF(IMPICK.EQ.2.AND.IMVEH.NE.0) CALL OBYIMP(ICY,IV)
      GOTO 90
C
C  LDDV, LDDT, & HDDV exclude tampering & all non-exhaust HC components.
C
   30 WRITE(IOUREP,230) BYVTN(IV),MONTH(MEVAL),ICY,VB,PSEC(IP),YFM,
     *    VB,PN(1),PN(3),PN(8),(U7,ICH=1,3)
  230 FORMAT(
     1  '-',T10,A26/
     2  ' ',T16,A4,' 1, ', I4/
     3  '0',T2,'Model',T21,A1,T32,A3/
     4  ' ',T2,A17,T21,A1,1X,3A7/
     4  '+',T2,'____',T7,'_____',T13,'_______',T22,3A7)
C
      DO 50 JDX=1,MAXYRS
      IDX=MAXYRS-JDX+1
      MY=ICY-JDX+1
C
      WRITE(IOUREP,240)
     5    MY,TF(IDX,IV),BYCUM(JDX,IV),
     5    BYBEF4(IP,JDX,IV),SALHCF(IDX,IP,IV),BYFER(IP,JDX,IV)
C
  240 FORMAT(' ',T2,I4,T7,F5.4,F8.0,T22,3F7.3)
C
   50 CONTINUE
C
C  Write column sums for TF weighted components, to check for identity with
C  the same sums written out in the OUTFMT = 3 table.
C
      WRITE(IOUREP,250) U7,EFFTP(IP,IV)
  250 FORMAT('+',T36,A7/' ',T36,F7.3)
      GOTO 90
C
C  Motorcycles exclude tampering plus running & refueling loss HC.
C
   60 WRITE(IOUREP,260) BYVTN(IV),MONTH(MEVAL),ICY,VB,PSEC(IP),YFM,
     *    VB,PN(1),PN(3),PN(4),PN(7),PN(8),
     *    (U7,ICH=1,5)
  260 FORMAT(
     1  '-',T13,A26/
     2  ' ',T20,A4,' 1, ', I4/
     3  '0',T2,'Model',T21,A1,T36,A3/
     4  ' ',T2,A17,T21,A1,1X,5A7/
     4  '+',T2,'____',T7,'_____',T13,'_______',T22,5A7)
C
      DO 80 JDX=1,MAXYRS
      IDX=MAXYRS-JDX+1
      MY=ICY-JDX+1
C
      WRITE(IOUREP,270)
     5    MY,TF(IDX,IV),BYCUM(JDX,IV),
     5    BYBEF4(IP,JDX,IV),SALHCF(IDX,IP,IV),BYEVAP(JDX,IV),
     5    BYRSTL(JDX,IV),BYFER(IP,JDX,IV)
C
  270 FORMAT(' ',T2,I4,T7,F5.4,F8.0,T22,5F7.3)
C
   80 CONTINUE
C
C  Write column sums for TF weighted components, to check for identity with
C  the same sums written out in the OUTFMT = 3 table.
C
      WRITE(IOUREP,280) (U7,ICH=1,3),
     *                  EFEVAP(IV),EFRSTL(IV),EFFTP(IP,IV)
  280 FORMAT('+',T36,3A7/' ',T36,3F7.3)
C
   90 CONTINUE
C
      RETURN
      END
