C Version @(#)$Id$ $Source$ $Date$ 

      SUBROUTINE OBYIMP(ICY,IV)
C
C  OBYIMP prints the by model year I/M credits and reductions table for the
C  selected evaluation year, pollutant (HC or CO or both?), and vehicle class
C  (LDGV or LDGT1 or LDGT2 or HDGV), as determined by the argument list.
C
C  Called by OUTBY1, OUTBY2, and OUTBY4.
C
C  Calls OUTPOL.
C
C  Input on call:
C
C    parameter list: ICY,IP,IV
C    common blocks:
C    /BYMYC1/ BYCUM
C    /BYMYC4/ YFM,PN,U7,PSEC,VB
C    /BYMYC5/ BYIMCR,BYIMRE,REDSUM
C    /IOUCOM/ IOUREP
C    /MAXIMA/ MAXPOL,MAXYRS
C    /MYRCAL/ TF
C
C  Local variable / array dictionary:
C
C   Name   Type              Description
C  ------  ----  -------------------------------------------------------
C  IMPROG  C*6   I/M section header with 6 column displacement available
C  IMTABL  C*6   IMPROG with (both HC & CO) or without (HC or CO) displacement
C  NCOLIM   I    number of columns in the I/M section of the optional table
C  NPOLL    I    number of pollutants (1 or 2) covered in the I/M section
C
C  Notes:
C
C  (1) OBYIMP was added for MOBILE4 In-House Version 02 (By Model Year Tables).
C
C
      CHARACTER*26 BYVTN
      CHARACTER*17 YFM
      CHARACTER*7 PN,U7
      CHARACTER*3 PSEC
      CHARACTER*1 VB
C
      COMMON /BYMYC1/ IVPICK(8),IMPICK,BYCUM(25,8)
      COMMON /BYMYC4/ BYVTN(8),YFM,PN(10),U7,PSEC(3),VB
      COMMON /BYMYC5/ BYIMCR(3,25,4),BYIMRE(3,25,4),REDSUM(3,4)
      COMMON /IOUCOM/ IOUIMD,IOUGEN,IOUREP,IOUERR,IOUASK
      COMMON /MAXIMA/ MAXVEH,MAXLTW,MAXPOL,MAXREG,MAXYRS
      COMMON /MYRCAL/ XMYM(25,8),JANMYR(25,8),TF(25,8),TFMYM(25,8)
C
      CHARACTER*8 IMPROG(4),IMTABL(3)
      DATA IMPROG/'        ','   I/M P','rogram  ','        '/
C
C  OBYIMP is designed to handle HC or CO or both, using pointers IP1, IP2,
C  and NCOLIM with implicit DO loops in the WRITE statements.
C
      CALL OUTPOL(IP1,IP2)
      NPOLL=IP2-IP1+1
      NCOLIM=NPOLL*2
      IF(NPOLL.GT.2) NPOLL=2
C
      DO 10 ICH=1,3
      IMTABL(ICH)=IMPROG(ICH-NPOLL+2)
   10 CONTINUE
C
      WRITE(IOUREP,200) VB,IMTABL,VB,(PSEC(IP),IP=IP1,IP2)
      WRITE(IOUREP,210) YFM,VB,(PN(9),PN(10),IP=IP1,IP2)
      WRITE(IOUREP,220) (U7,ICH=1,NCOLIM)
  200 FORMAT('-'//' ',T21,A1,T23,3A8/
     1       ' ',T2,'Model',T21,A1,T28,A3,T42,A3,T56,A3)
  210 FORMAT(' ',T2,A17,T21,A1,1X,6A7)
  220 FORMAT('+',T2,'____',T7,'_____',T13,'_______',T22,6A7)
C
      DO 20 JDX=1,MAXYRS
      IDX=MAXYRS-JDX+1
      MY=ICY-JDX+1
C
      WRITE(IOUREP,230)
     1    MY,TF(IDX,IV),BYCUM(JDX,IV),
     2    (BYIMCR(IP,JDX,IV),BYIMRE(IP,JDX,IV),IP=IP1,IP2)
  230 FORMAT(' ',T2,I4,T7,F5.4,F8.0,T22,6F7.3)
C
   20 CONTINUE
C
C  Write column sums for TF weighted I/M reduction, to check for identity with
C  the difference between without & with I/M in the OUTFMT = 3 tables.
C
      WRITE(IOUREP,240) (U7,IP=IP1,IP2)
      WRITE(IOUREP,250) (REDSUM(IP,IV),IP=IP1,IP2)
  240 FORMAT('+',T29,A7,T43,A7,T57,A7)
  250 FORMAT(' ',T29,F7.3,T43,F7.3,T57,F7.3)
C
      RETURN
      END
