C Version "@(#)$Id$ $Source$ $Date$ 

      SUBROUTINE OUTBY4(ICY)
C
C  OUTBY4 writes scenario emission factor results by model year for selected
C  vehicle types on unit IOUREP using output format (OUTFMT) 5, when all 3
C  pollutants are requested (PRTFLG=4).
C
C  Called by OUTDT5.
C
C  Calls MYRIMT, and OBYIMP.
C
C  Input on call:
C
C    parameter list: ICY
C    common blocks:
C    /BYMYC1/ IVPICK,IMPICK,BYCUM
C    /BYMYC2/ BYBEF4,BYTAM,BYFER
C    /BYMYC3/ BYEVAP,BYRUNL,BYREFL,BYRSTL
C    /BYMYC4/ BYVTN,YFM,PN,U7,PSEC,VB
C     /EVAL/ MEVAL,MONTH(12)
C    /IOUCOM/ IOUREP
C    /MAXIMA/ MAXVEH,MAXYRS
C    /MYRCAL/ TF
C    /LEVBLK/ LEVFLG
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
C  (1) OUTBY4 was added for MOBILE4 In-House Version 02 (By Model Year Tables).
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
      INCLUDE 'MAXIMA.I'
      INCLUDE 'MYRCAL.I'
      INCLUDE 'LEVBLK.I'
      INCLUDE 'RESUL1.I'
      INCLUDE 'RESUL2.I'
      INCLUDE 'SPEED6.I'
C
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
     1'-',T54,A26,/,
     1' ',T59,A4,' 1 ', I4,/,
     2' ','                 ')
C
C  (22-April-1994) @ CSC ljn Request 2-420
C  write out the Zero Emitting Vehicles BYMY Travel Fraction
C
      IF(LEVFLG.EQ.2.AND.IV.LE.3) CALL OUTZEV(IV,ICY)
C
C
C    Finish writing
C
      WRITE(IOUREP,211) (VB,PSEC(IP),IP=1,3),YFM,VB,
     *  (PN(ICH),ICH=1,8),(VB,PN(1),PN(2),PN(3),PN(8),ICH=2,3),
     *  (U7,ICH=1,16)
  211 FORMAT(
     3  '0',T2,'Model',T20,A1,T48,A3,T76,A1,T91,A3,T105,A1,T118,A3/
     4' ',T2,A17,T20,A1,8A7,T76,A1,1X,4A7,T105,A1,4A7/
     5'+',T2,'____',T7,'_____',T13,'_______',T20,8A7,T77,4A7,T105,4A7)
C
      IMVEH=0
      DO 20 JDX=1,MAXYRS
      IDX=MAXYRS-JDX+1
      MY=ICY-JDX+1
      IF(IMVEH.EQ.0) IMVEH=MYRIMT(MY,IV,ICY)
C
      WRITE(IOUREP,215)
     5    MY,TF(IDX,IV),BYCUM(JDX,IV),BYBEF4(1,JDX,IV),
     5    BYTAM(1,JDX,IV),SALHCF(IDX,1,IV),BYEVAP(JDX,IV),
     5    BYREFL(JDX,IV),BYRUNL(JDX,IV),BYRSTL(JDX,IV),BYFER(1,JDX,IV),
     5    (BYBEF4(IP,JDX,IV),BYTAM(IP,JDX,IV),SALHCF(IDX,IP,IV),
     5    BYFER(IP,JDX,IV),IP=2,3)
C
  215 FORMAT(' ',T2,I4,T7,F5.4,F8.0,T20,8F7.3,T77,4F7.3,T105,4F7.3)
C
   20 CONTINUE
C
C  Write column sums for TF weighted components, to check for identity with
C  the same sums written out in the OUTFMT = 3 table.
C
      WRITE(IOUREP,220) (U7,ICH=1,7),
     1      EFEVAP(IV),EFREFL(IV),EFRUNL(IV),EFRSTL(IV),
     1     (EFFTP(IP,IV),IP=1,3)
C
  220 FORMAT('+',T41,5A7,T98,A7,T126,A7/
     *       ' ',T41,5F7.3,T98,F7.3,T126,F7.3)
C
C  If By Model Year I/M Table picked, print it.
C
      IF(IMPICK.EQ.2.AND.IMVEH.NE.0) CALL OBYIMP(ICY,IV)
      GOTO 90
C
C  LDDV, LDDT, & HDDV exclude tampering & all non-exhaust HC components.
C
   30 WRITE(IOUREP,230) BYVTN(IV),MONTH(MEVAL),ICY,(VB,PSEC(IP),IP=1,3),
     *YFM,
     *    (VB,PN(1),PN(3),PN(8),ICH=1,3),(U7,ICH=1,9)
  230 FORMAT(
     1  '-',T33,A26/
     2  ' ',T38,A4,' 1, ', I4/
     3  '0',T2,'Model',T21,A1,T32,A3,T44,A1,T55,A3,T67,A1,T78,A3/
     4  ' ',T2,A17,T21,A1,1X,3A7,T44,A1,1X,3A7,T67,A1,1X,3A7/
     4  '+',T2,'____',T7,'_____',T13,'_______',T22,3A7,T45,3A7,T68,3A7)
C
      DO 50 JDX=1,MAXYRS
      IDX=MAXYRS-JDX+1
      MY=ICY-JDX+1
C
      WRITE(IOUREP,240)
     5    MY,TF(IDX,IV),BYCUM(JDX,IV),
     5    (BYBEF4(IP,JDX,IV),SALHCF(IDX,IP,IV),BYFER(IP,JDX,IV),IP=1,3)
C
  240 FORMAT(' ',T2,I4,T7,F5.4,F8.0,T22,3F7.3,T45,3F7.3,T68,3F7.3)
C
   50 CONTINUE
C
C  Write column sums for TF weighted components, to check for identity with
C  the same sums written out in the OUTFMT = 3 table.
C
      WRITE(IOUREP,250) (U7,ICH=1,3),(EFFTP(IP,IV),IP=1,3)
  250 FORMAT('+',T36,A7,T59,A7,T82,A7/' ',T36,F7.3,T59,F7.3,T82,F7.3)
      GOTO 90
C
C  Motorcycles exclude tampering plus running & refueling loss HC.
C
   60 WRITE(IOUREP,260) BYVTN(IV),MONTH(MEVAL),ICY,(VB,PSEC(IP),IP=1,3),
     *    YFM,VB,PN(1),PN(3),PN(4),PN(7),PN(8),
     *    (VB,PN(1),PN(3),PN(8),ICH=2,3),
     *    (U7,ICH=1,11)
  260 FORMAT(
     1  '-',T39,A26/
     2  ' ',T46,A4,' 1, ', I4/
     3  '0',T2,'Model',T21,A1,T39,A3,T58,A1,T69,A3,T81,A1,T92,A3/
     4  ' ',T2,A17,T21,A1,1X,5A7,T58,A1,1X,3A7,T81,A1,1X,3A7/
     4  '+',T2,'____',T7,'_____',T13,'_______',T22,5A7,T59,3A7,T82,3A7)
C
      DO 80 JDX=1,MAXYRS
      IDX=MAXYRS-JDX+1
      MY=ICY-JDX+1
C
      WRITE(IOUREP,270)
     5    MY,TF(IDX,IV),BYCUM(JDX,IV),BYBEF4(1,JDX,IV),
     5    SALHCF(IDX,1,IV),BYEVAP(JDX,IV),BYRSTL(JDX,IV),
     5    BYFER(1,JDX,IV),
     5    (BYBEF4(IP,JDX,IV),SALHCF(IDX,IP,IV),BYFER(IP,JDX,IV),IP=2,3)
C
  270 FORMAT(' ',T2,I4,T7,F5.4,F8.0,T22,5F7.3,T59,3F7.3,T82,3F7.3)
C
   80 CONTINUE
C
C  Write column sums for TF weighted components, to check for identity with
C  the same sums written out in the OUTFMT = 3 table.
C
      WRITE(IOUREP,280) (U7,ICH=1,5),EFEVAP(IV),EFRSTL(IV),
     *                  (EFFTP(IP,IV),IP=1,3)
  280 FORMAT('+',T36,3A7,T73,A7,T96,A7/
     *       ' ',T36,3F7.3,T73,F7.3,T96,F7.3)
C
   90 CONTINUE
C
      RETURN
      END
