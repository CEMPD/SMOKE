C Version "@(#)$Id$ $Source$ $Date$ 

      SUBROUTINE BIGCFX(ICY)
C
C  BIGCFX computes correction factors applied to composite emission
C  factors, e.g., speed, optional (humidity, trailer towing, a/c), temperature,
C  by-bag RVP, operating mode, methane and tampering.
C
C  Called by EFCALX.
C
C  Calls BIGALH, BIGIN1, BIGIN2, BIGMET, BIGOMT, BIGSAL, BIGTAM and BIGTCF.
C
C  Input on call:
C
C    parameter list: ICY
C    common blocks:
C    /FLAGS2/ ALHFLG
C    /FLAGS4/ PRTFLG,NMHFLG
C    /MAXIMA/ MAXVEH,MAXPOL,MAXYRS
C    /MYRCAL/ TF
C
C  Output on return:
C
C    common blocks:
C    /TEMPC4/ LOWCO
C
C  Notes:
C
C  BIGCFX was modified for MOBILE4 In-House Version 02 to skip calculation of
C  correction factors when because travel fraction for the IV x IDX case is 0,
C  the CFs will not be used.  SALHCF for such cases will print 0.0 in the
C  BYMY table.
C
C
      INCLUDE 'FLAGS2.I'
      INCLUDE 'FLAGS4.I'
      INCLUDE 'MAXIMA.I'
      INCLUDE 'MYRCAL.I'
      INCLUDE 'TEMPC4.I'
C
C  Initialize variables to be used in upcoming loops.
C
      CALL BIGIN1
      IF(ALHFLG.GT.1) CALL BIGALH(ICY)
C
C  Enter the correction factor loops.
C
      DO 12 IDX=1,MAXYRS
      JDX=(MAXYRS+1)-IDX
      MY=ICY-(MAXYRS-IDX)
C
C  Skip calculations if travel fraction is 0.0, zeroing emission factor.
C
      DO 11 IV=1,MAXVEH
C  Do not put TF branch check here, as EF results change in benchmark.
      CALL BIGIN2(IV,MY,JDX,IDX,IVLM,IVTEMP,ISP)
      IF(TF(IDX,IV).EQ.0.0) GOTO 11
C
      DO 10 IP=1,MAXPOL
      IF(PRTFLG.NE.IP.AND.PRTFLG.NE.4) GOTO 10
C
C  Reset CO offset switch and calculate temperature correction factor (TCF).
C  There is no TCF in MOBILE4 for diesel IV (IVTEMP=0).  When the ambient
C  temperature is the FTP standard 75 degrees (MDLOHI=0), TCF represents the
C  RVP CF for RVP MYG cases 2 and 3 for LDGV/T & HDGV, otherwise it remains the
C  initialized default neutral value (1.0) for all IV (1-8).
C
      LOWCO=1
      IF(IVTEMP.GE.1) CALL BIGTCF(ICY,IP,IVTEMP,MY,IDX)
C
C  Compute speed and, if LDGV/T or MC, optional variables correction factor.
C
      CALL BIGSAL(IP,IV,MY,IDX,ISP)
C
C  Check for and calculate if needed the methane offset to composite HC.
C
      IF(IP.EQ.1) CALL BIGMET(IV,MY,IDX,JDX)
C
C  If LDGV/T or HDGV, correct tampering offset for opmode and temperature.
C
      IF(IV.LE.4) CALL BIGTAM(IP,IV,IDX)
C
C  Compute bef's multiplicative correction for operating mode and temperature.
C
      CALL BIGOMT(IP,IV,IVLM,MY,IDX)
C
   10 CONTINUE
   11 CONTINUE
   12 CONTINUE
C
      RETURN
      END
