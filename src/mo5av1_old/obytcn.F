C Version @(#)$Id$ $Source$ $Date$ 

      SUBROUTINE OBYTCN(ICY,IP,IV)
C
C  OUTBY4 writes the Terse by model year emission factors tables for selected
C  vehicle type on unit IOUREP using output format (OUTFMT) 5, when either CO
C  or NOX emission factors are requested (PRTFLG=2 or 3).
C
C  Called by OUTBY2 and OUTBY3.
C
C  Input on call:
C
C    parameter list: ICY,IP,IV
C    common blocks:
C    /BYMYC1/ BYCUM
C    /BYMYC2/ BYBEF4,BYTAM,BYFER
C    /BYMYC4/ BYVTN,YFM,PN,U7,PSEC,VB
C    /EVAL/   MEVAL,MONTH
C    /IOUCOM/ IOUREP
C    /LEVBLK/ LEVFLG
C    /MAXIMA/ MAXYRS
C    /MYRCAL/ TF
C    /RESUL1/ EFFTP
C    /SPEED6/ SALHCF
C
C  Notes:
C
C  (1) OBYTCN was added for MOBILE4 In-House Version 02 (By Model Year Tables).
C
C
      CHARACTER*26 BYVTN
      CHARACTER*17 YFM
      CHARACTER*7 PN,U7
      CHARACTER*4 MONTH
      CHARACTER*3 PSEC
      CHARACTER*1 VB
C
      COMMON /BYMYC1/ IVPICK(8),IMPICK,BYCUM(25,8)
      COMMON /BYMYC2/ BYBEF4(3,25,8),BYTAM(3,25,4),BYFER(3,25,8)
      COMMON /BYMYC4/ BYVTN(8),YFM,PN(10),U7,PSEC(3),VB
      COMMON /EVAL/ MEVAL,MONTH(12)
      COMMON /IOUCOM/ IOUIMD,IOUGEN,IOUREP,IOUERR,IOUASK
      COMMON /LEVBLK/ LEVFLG,IMFLG,LEVYRS
      COMMON /MAXIMA/ MAXVEH,MAXLTW,MAXPOL,MAXREG,MAXYRS
      COMMON /MYRCAL/ XMYM(25,8),JANMYR(25,8),TF(25,8),TFMYM(25,8)
      COMMON /RESUL1/ EFFTP(3,9),EFEXH(9),EFEVAP(9),EFREFL(9),EFRUNL(9)
      COMMON /SPEED6/ SALHCF(25,3,8),SCFIDL(25,4),SCIADJ(25,3,8)
C
C  Components printed and hence BYMY table format vary with the vehicle class.
C
      IF(IV.GE.5) GOTO 20
C
C  LDGV, LDGT1, LDGT2 & HDGV include tampering offsets.
C
C
      WRITE(IOUREP,180) BYVTN(IV), MONTH(MEVAL), ICY
  180 FORMAT(
     *'-',T54,A26,/,
     *' ',T59,A4,' 1 ', I4,/,
     *' ','                 ')
C
C    write out the Zero Emitting Vehicles BYMY Travel Fraction
C
      IF(LEVFLG.EQ.2.AND.IV.LE.2) CALL OUTZEV(IV,ICY)
C
C    Finish writing
C
      WRITE(IOUREP,200) VB,PSEC(IP),YFM,VB,
     *    PN(1),PN(2),PN(3),PN(8),(U7,ICH=1,4)
  200 FORMAT(
     *    '0',T2,'Model',T21,A1,T35,A3/
     *    ' ',T2,A17,T21,A1,1X,4A7/
     *    '+',T2,'____',T7,'_____',T13,'_______',T22,4A7)
C
C
      DO 10 JDX=1,MAXYRS
      IDX=MAXYRS-JDX+1
      MY=ICY-JDX+1
C
      WRITE(IOUREP,210)
     5    MY,TF(IDX,IV),BYCUM(JDX,IV),BYBEF4(IP,JDX,IV),
     5    BYTAM(IP,JDX,IV),SALHCF(IDX,IP,IV),BYFER(IP,JDX,IV)
C
  210 FORMAT(' ',T2,I4,T7,F5.4,F8.0,T22,4F7.3)
C
   10 CONTINUE
C
C  Write column sums for TF weighted components, to check for identity with
C  the same sums written out in the OUTFMT = 3 table.
C
      WRITE(IOUREP,220) U7,EFFTP(IP,IV)
C
  220 FORMAT('+',T43,A7/' ',T43,F7.3)
      RETURN
C
C  LDDV, LDDT, HDDV and MC do not have tampering offsets.
C
   20 WRITE(IOUREP,230) BYVTN(IV),MONTH(MEVAL),ICY,VB,PSEC(IP),YFM,
     *    VB,PN(1),PN(3),PN(8),(U7,ICH=1,3)
  230 FORMAT(
     1  '-',T9,A26/
     2  ' ',T16,A4,' 1, ', I4/
     3  '0',T2,'Model',T21,A1,T32,A3/
     4  ' ',T2,A17,T21,A1,1X,3A7/
     4  '+',T2,'____',T7,'_____',T13,'_______',T22,3A7)
C
      DO 30 JDX=1,MAXYRS
      IDX=MAXYRS-JDX+1
      MY=ICY-JDX+1
C
      WRITE(IOUREP,240)
     5    MY,TF(IDX,IV),BYCUM(JDX,IV),
     5    BYBEF4(IP,JDX,IV),SALHCF(IDX,IP,IV),BYFER(IP,JDX,IV)
C
  240 FORMAT(' ',T2,I4,T7,F5.4,F8.0,T22,3F7.3)
C
   30 CONTINUE
C
C  Write column sums for TF weighted components, to check for identity with
C  the same sums written out in the OUTFMT = 3 table.
C
      WRITE(IOUREP,250) U7,EFFTP(IP,IV)
  250 FORMAT('+',T36,A7,/' ',T36,F7.3)
C
      RETURN
      END
