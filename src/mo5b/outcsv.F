C Version "@(#)$Id$ $Source$ $Date$ 

      SUBROUTINE OUTCSV(ICY)
C
C  OUTCSV makes the output into a CSV format and in one single record.
C
C  Called by OUTPUT.
C
C  Calls OUTPOL
C
C  Input on call:
C
C    common blocks:
C    /ATPAR1/ LAPSY,LAP1ST,LAPLST,LVTFLG
C    /ATPAR2/ ATPPGM,ATPFQT,CRATP,DISTYP
C    /BASE10/ ZEVTOT
C    /CITPAR/ SCNAME
C    /CITRV1/ RVPBAS,RVPIUS,IUSESY
C    /ETPPAR/ IETPYR,ETPPIV
C    /EVAL/   MEVAL
C    /FLAGS1/ TAMFLG,VMFLAG,OXYFLG,DSFLAG
C    /FLAGS2/ MYMRFG,NEWFLG,IMFLAG,ALHFLG,IMTIER
C    /FLAGS3/ ATPFLG,OUTFMT
C    /FLAGS4/ PRTFLG,NMHFLG,HCFLAG
C    /IMPAR1/ ICYIM,ISTRIN,MODYR,WAIVER,CRIM
C    /IMPAR2/ ILDT,ITEST
C    /IMPAR4/ CUTHC,CUTCO,CUTNO
C    /IMPAR6/ IFREQ,INTYP,IFILET
C    /IOUCOM/ IOUREP
C    /LEVBLK/ LEVFLG,IMFLG,LEVYRS
C    /MAXIMA/ MAXVEH
C    /OXY1/   SHRMKT,OXYCNT,IGASHW
C    /REGION/ IREJN
C    /RESUL1/ EFFTP,EFEXH
C    /RESUL3/ VFTP
C    /RFORM1/ RFGFLG
C    /RLCOM2/ IOBMY,IVOB,IS2SY,NPHASE,S2EFF
C    /RUNLS1/ TLVMTU
C    /PRGCH1/ LPGSY,LPG1ST,LPGLST,LVPRG
C    /PRGCH2/ PRGPGM,PRGFQT,CRPRG,PRGFLG
C    /PROJEC/ PROJID
C    /PRSCH1/ LPRSY,LPR1ST,LPRLST,LVPRS
C    /PRSCH2/ PRSPGM,PRSFQT,CRPRS,PRSFLG
C    /SCENE1/ SPD,PCCN,PCHC,PCCC
C    /STRING/ NOYES
C    /TEMPS/  TEMMIN,TEMMAX,TEMEXH
C    /YEARS4/ IY1994
C    /VMXCOM/ VMTMIX
C
C  Local array subscripts
C
C  NAMEVP (  5 )  IEVAP
C  NAMNMH (  5 )  NMHFLG
C  NAMREG (  2 )  IREJN
C  NAMVEH ( 10 )  IV 1-9 and I 10
C  POLBL  (  3 )  IP
C
C
C  Local array / variable dictionary:
C
C   Name   Type              Description
C  ------  ----  ------------------------------------------------
C  F200    C*23  Project identification, and evaluation month
C  F201    C*9   Calendar year format
C  F202    C*9   HC composite EF type format, RFGFLG (yes/no),
C                Region (Low/High) and Tier 1(I/M) (yes/no),
C                use F214 when not defined.
C  F203    C*9   Vehicle classes format
C  F204    C*11  Composite EF value, and VMT and speed value
C                format when defined, see F207 when not defined.
C                Also ZEV% when LEVFLG=2, LEVYRS is less than or
C                equal ICY, and IV=1,2; use F214 when not defined.
C  F205    C*10  non-HC exhaust EF type format
C  F206    C*9   HC component EF type format
C  F208    C*85  Ambient temperature is a function of IP, also
C                the three temp, Amb, Min, and Max are in one
C                field.  RVP period 1 and 2 and the start year
C                for period 2 are in one field.
C                Operating mode (3) are in one field
C  F209    C*6   IMFLAG=1 (when I/M is not applicable)
C  F210    C*108 I/M Program #1 Descriptive Parameter Record
C  F211    C*239 I/M Program 1&2 Descriptive Parameter Records
C  F212    C*49  Anti-Tampering Descriptive Parameter Record
C  F213    C*40  Pressure/Purge Descriptive Parameter Records
C  F214    C*5   ATPFLG=1, PRSFLG=1, PRGFLG=1, RFGFLG=1,
C                SEAFLG=2, OXYFLG=1, IMTIER=0, IMFLG=1, LEVFLG=1,
C                NEWFLG=1, RLFLAG=1, SPDFLG=1
C  F215    C*20  Oxy record
C  F216    C*9   Scenario
C  F217    C*17  LEV Parameter Record use F214 when not defined.
C  F218    C*5   ALHFLG, DSFLAG, MYMRFG, MYMRFG, TAMFLG, and VMFLAG
C                in the last field for flags
C  F219    C*9   NEWFLG in the last field for flags
C  F220    C*24  ETP record
C  F221    C*28  Stage 2 Vapor Recovery System record
C  F222    C*18  Onboard Vapor Recovery System record
C  F223    C*17  Trip Length Distribution record
C  I       I     Vehicle class in horizontal appearance for
C                (OUTFMT=3) ranges from 1 to 10, where 1=LDGV
C                2=LDGT1, 3=LDGT2, 4=LDGT, 5=HDGV, 6=LDDV,
C                7=LDDT, 8=HDDV, 9=MC, and 10=total
C  IEVAP   I     HC component in vertical appearance
C                1=Exhuast HC, 2=Evaporative HC, 3=Refueling,
C                4=Running Losses, and 5=Resting Losses
C  INTLBL  C*7   Interpolated month of evalution is chosen
C                between interpolated month or regular month
C  IP      I     Pollutant values range from 1 to 3
C  IP1     I     Pollutant lower bound
C  IP2     I     Pollutant upper bound
C  IV      I     Vehicle class values range from 1 to 8, and
C                1-8 are as defined elsewhere 9=LDGT
C  MTHLBL  C*7   Month of evaluation
C  NAMEVP  C*12  Type of HC EF component
C  NAMNMH  C*4   Type of HC emission factor
C  NAMREG  C*4   region names
C  NAMVEH  C*5   Vehicle classes
C  ONE     R     All vehicle miles traveled fractions sum to one
C  REGLBL  C*7   Regular month of evalution
C  POLBL   C*11  Pollutant exhaust types
C  Q       C*1   Double quote can be used as an open or close
C  QW      C*3   Double quote can only be used as a separator
C
C
C  Notes:
C
C  OUTFMT 7 contains the same set of fields as MOBILE5's OUTFMT 3,
C  except that
C
C  (1) The fields are all on 1 record, for a given scenario,
C      vehicle and pollutant.
C
C  (2) Field print order is not comparable.
C
C  (2) The emission factors are printed with 3, instead of 2,
C      decimal places.
C
C  (3) The project identification label (new for MOBILE4) is
C      appended as the first field.
C
C  (4) The scenario identification label (new for MOBILE4) is
C      appended as the last field.
C
C  24-May-1993 @ ARC-bsg Subtask 222, Comma Separated Variables
C  Cadillac version preceded by header
C
C  29 July 1994 @ CSC-tmm Request #2-457 Replaced WINFLG with SEAFLG 
C
      INCLUDE 'ATPAR1.I'
      INCLUDE 'ATPAR2.I'
      INCLUDE 'BASE10.I'
      INCLUDE 'CITPAR.I'
      INCLUDE 'CITRV1.I'
      INCLUDE 'ETPPAR.I'
      INCLUDE 'EVAL.I'
      INCLUDE 'FLAGS1.I'
      INCLUDE 'FLAGS2.I'
      INCLUDE 'FLAGS3.I'
      INCLUDE 'FLAGS4.I'
      INCLUDE 'IMPAR1.I'
      INCLUDE 'IMPAR2.I'
      INCLUDE 'IMPAR4.I'
      INCLUDE 'IMPAR6.I'
      INCLUDE 'IOUCOM.I'
      INCLUDE 'LEVBLK.I'
      INCLUDE 'MAXIMA.I'
      INCLUDE 'OXY1.I'
      INCLUDE 'PRGCH1.I'
      INCLUDE 'PRGCH2.I'
      INCLUDE 'PROJEC.I'
      INCLUDE 'PRSCH1.I'
      INCLUDE 'PRSCH2.I'
      INCLUDE 'REGION.I'
      INCLUDE 'RESUL1.I'
      INCLUDE 'RESUL2.I'
      INCLUDE 'RESUL3.I'
      INCLUDE 'RFORM1.I'
      INCLUDE 'RLCOM2.I'
      INCLUDE 'RUNLS1.I'
      INCLUDE 'SCENE1.I'
      INCLUDE 'STRING.I'
      INCLUDE 'TEMPS.I'
      INCLUDE 'YEARS4.I'
      INCLUDE 'VMXCOM.I'
C
      CHARACTER*7 REGLBL,INTLBL,MTHLBL
      CHARACTER*4 NAMNMH(5),HCLBL,NAMREG(2)
      CHARACTER*5 NAMVEH(10)
      CHARACTER*12 NAMEVP(5)
      CHARACTER*1 Q
      REAL ONE
      CHARACTER*3 QW
      CHARACTER*23 F200
      CHARACTER*11 POLBL(3)
      CHARACTER*9 F201, F202, F203, F216, F219
      CHARACTER*11 F204
      CHARACTER*10 F205, F206
      CHARACTER*6  F209
      CHARACTER*85 F208
      CHARACTER*108 F210
      CHARACTER*239 F211
      CHARACTER*49 F212
      CHARACTER*40 F213
      CHARACTER*5 F214, F218
      CHARACTER*20 F215
      CHARACTER*17 F217, F223
      CHARACTER*24 F220
      CHARACTER*29 F221
      CHARACTER*18 F222
      DATA F200 / '(A1, 20A4, A3, A7, A3\)'/
      DATA F201 / '(I4, A3\)' /
      DATA F202 / '(A4, A3\)' /
      DATA F203 / '(A5, A3\)' /
      DATA F204 / '(F8.3, A3\)' /
      DATA F205 / '(A11, A3\)' /
      DATA F206 / '(A12, A3\)' /
      DATA F208 / '(F5.1, 2(1X, F4.0), '' (F)'', A3, 2(F4.1, 1X),
     *I4, A3, F5.1, 2('' \ '', F5.1), A3\)' /
      DATA F209 / '(2A3\)' /
      DATA F210 / '(I4, 1X, I2, 2(1X, I4), 2(1X, F3.2), 1X, F8.3,
     *2(1X, I1), 1X, 4I1, 1X, I1, 1X, 2I1, 3(1X, F6.2), 2A3\)' /
      DATA F211 / '(I4, 1X, I2, 2(1X, I4), 2(1X, F3.2), 1X, F8.3,
     *2(1X, I1), 1X, 4I1, 1X, I1, 1X, 2I1, 3(1X, F6.2), A3,
     *I4, 1X, I2, 2(1X, I4), 2(1X, F3.2), 1X, F8.3,
     *2(1X, I1), 1X, 4I1, 1X, I1, 1X, 2I1, 3(1X, F6.2), A3\)' /
      DATA F212 / '(3(I4, 1X), 4I1, 1X, 2I1, 1X, F8.3, 1X, 8I1, A3\)' /
      DATA F213 / '(3(I4, 1X), 4I1, 1X, 2I1, 1X, F8.3, A3\)' /
      DATA F214 / '(A3\)' /
      DATA F215 / '(4F6.3, 1X, A4, A3\)' /
      DATA F216 / '(4A4, A1)' /
      DATA F217 / '(I4, 1X, I1, A3\)' /
      DATA F218 / '(I1\)' /
      DATA F219 / '(I1, A3\)' /
      DATA F220 / '(4(I4,1X,F5.1,1X), A3\)' /
      DATA F221 / '(I4, 1X, I1, 2(1X,F5.3), A3\)' /
      DATA F222 / '(I4, 1X, 4I1, A3\)' /
      DATA F223 / '(6(F5.3,1X), A3\)' /
C       0123456789012345678901234567890123456789012345678901
      DATA Q / '"' /
      DATA QW / '","' /
      DATA ONE / 1. /
      DATA NAMREG / 'Low ', 'High' /
      DATA NAMVEH / 'LDGV ', 'LDGT1', 'LDGT2', 'HDGV ', 'LDDV ',
     *              'LDDT ', 'HDDV ', 'MC   ', 'LDGT ', 'All  ' /
      DATA REGLBL / 'January' /
      DATA INTLBL / 'July   ' /
      DATA NAMNMH / 'THC ', 'NMHC', 'VOC ', 'TOG ', 'NMOG' /
      DATA NAMEVP / 'HC Exhaust  ', 'Evap      ', 'Refueling   ',
     *              'Running     ', 'Resting   '/
      DATA POLBL / '           ', 'CO Exhaust ', 'NOx Exhaust' /
C
C  Set the HC label to the appropriate function of the NMHC FLAG
C  value.
C
      HCLBL=NAMNMH(NMHFLG)
C
C  Set the month label to either January or July
C
      IF(MEVAL.EQ.1) MTHLBL=REGLBL
      IF(MEVAL.EQ.7) MTHLBL=INTLBL
C
C  Define upper and lower bounds of pollutant values
C
      CALL OUTPOL(IP1,IP2)
C
C  Iterate for pollutant vlues
C
      DO 25 IP=IP1,IP2
        DO 10 I=1,10
          IV=I
          IF(I.EQ.4) IV=9
          IF(I.GT.4) IV=IV-1
          WRITE(IOUREP,F200) Q,PROJID,QW,MTHLBL,QW
          WRITE(IOUREP,F201) ICY,QW
C
C  Choose between HC type and non-HC exhaust type
C
          IF(IP.EQ.1) WRITE(IOUREP,F202) HCLBL,QW
          IF(IP.NE.1) WRITE(IOUREP,F205) POLBL(IP),QW
C
C  Choose between 9 vehicle names and total
C
          IF(I.NE.10) WRITE(IOUREP,F203) NAMVEH(IV),QW
          IF(I.EQ.10) WRITE(IOUREP,F203) NAMVEH(10),QW
C
C  Choose between 9 vehicle types and total for EF
C
          IF(I.NE.10) WRITE(IOUREP,F204) EFFTP(IP,IV),QW
          IF(I.EQ.10) WRITE(IOUREP,F204) VFTP(IP),QW
C
C  Choose appropriate VMT, leave blank or 1
C
          IF(I.NE.4.AND.I.NE.10)
     *      WRITE(IOUREP,F204) VMTMIX(IV),QW
          IF(I.EQ.4) WRITE(IOUREP,F214) QW
          IF(I.EQ.10) WRITE(IOUREP,F204) ONE,QW
C
C  Choose appropriate speed or leave blank
C
          IF(I.EQ.4.OR.I.EQ.10) THEN
            WRITE(IOUREP,F214) QW
          ELSE
            WRITE(IOUREP,F204) SPD(IV),QW
          ENDIF
          WRITE(IOUREP,F208) TEMEXH(IP),TEMMIN,TEMMAX,QW,
     *    RVPBAS,RVPIUS,IUSESY,QW,PCCN,PCHC,PCCC,QW
          IF(IMFLAG.EQ.1) WRITE(IOUREP,F209) QW,QW
          IF(IMFLAG.EQ.2) WRITE(IOUREP,F210) ICYIM(1),ISTRIN(1),
     *      MODYR(1,1),MODYR(2,1),WAIVER(1,1),WAIVER(2,1),CRIM(1),
     *      INTYP(1),IFREQ(1),(ILDT(IVIM,1),IVIM=1,4),ITEST(1),
     *      IFILET(1),IFILET(2),CUTHC(1),CUTCO(1),CUTNO(1),QW,QW
          IF(IMFLAG.EQ.3) WRITE(IOUREP,F211) ICYIM(1),ISTRIN(1),
     *      MODYR(1,1),MODYR(2,1),WAIVER(1,1),WAIVER(2,1),CRIM(1),
     *      INTYP(1),IFREQ(1),(ILDT(IVIM,1),IVIM=1,4),ITEST(1),
     *      IFILET(1),IFILET(2),CUTHC(1),CUTCO(1),CUTNO(1),QW,
     *      ICYIM(2),ISTRIN(2),
     *      MODYR(1,2),MODYR(2,2),WAIVER(1,2),WAIVER(2,2),CRIM(2),
     *      INTYP(2),IFREQ(2),(ILDT(IVIM,2),IVIM=1,4),ITEST(2),
     *      IFILET(1),IFILET(2),CUTHC(2),CUTCO(2),CUTNO(2),QW
C
C  Choose when to print ATP record
C
          IF(ATPFLG.EQ.1) WRITE(IOUREP,F214) QW
          IF(ATPFLG.EQ.2) WRITE(IOUREP,F212)
     *      LAPSY,LAP1ST,LAPLST,LVTFLG,ATPPGM,ATPFQT,CRATP,DISTYP,QW
          IF(PRSFLG.EQ.1) WRITE(IOUREP,F214) QW
          IF(PRSFLG.EQ.2) WRITE(IOUREP,F213)
     *      LPRSY,LPR1ST,LPRLST,LVPRS,PRSPGM,PRSFQT,CRPRS,QW
          IF(PRGFLG.EQ.1) WRITE(IOUREP,F214) QW
          IF(PRGFLG.EQ.2) WRITE(IOUREP,F213)
     *      LPGSY,LPG1ST,LPGLST,LVPRG,PRGPGM,PRGFQT,CRPRG,QW
          WRITE(IOUREP,F202) NOYES(RFGFLG),QW
C
C  Only write Oxy record when OXYFLG=2, provided that RFGFLG=1,
C  or when OXYFLG=2 and when RFGFLG=2 and SEAFLG=2
C
          IF(OXYFLG.EQ.2.AND.RFGFLG.EQ.2.AND.SEAFLG.EQ.2)
     *      WRITE(IOUREP,F215) SHRMKT(2),SHRMKT(3),
     *      OXYCNT,NOYES(IGASHW),QW
          IF(OXYFLG.EQ.2.AND.RFGFLG.EQ.2.AND.SEAFLG.EQ.1)
     *      WRITE(IOUREP,F214) QW
          IF(OXYFLG.EQ.2.AND.RFGFLG.EQ.1)
     *      WRITE(IOUREP,F215) SHRMKT(2),SHRMKT(3),
     *      OXYCNT,NOYES(IGASHW),QW
          IF(OXYFLG.EQ.1)
     *      WRITE(IOUREP,F214) QW
          WRITE(IOUREP,F202) NAMREG(IREJN),QW
C
C  Only write Tier 1(I/M), LEV (I/M) or ZEV% when IV=1,2, and
C  ICY is greater than or equal to 1994.  And only write
C  Tier 1(I/M) when there is an I/M.  And only write LEV (I/M)
C  and ZEV% when there is an LEV (ICY is greater than or equal
C  to LEVYRS and LEVFLG=2).
C
          IF((I.EQ.1.OR.I.EQ.2).AND.IMFLAG.NE.1.AND.
     *    ICY.GE.IY1994) THEN
            WRITE(IOUREP,F202) NOYES(IMTIER+1),QW
          ELSE
            WRITE(IOUREP,F214) QW
          ENDIF
C
          IF((I.EQ.1.OR.I.EQ.2).AND.LEVFLG.NE.1.AND.
     *    ICY.GE.LEVYRS) THEN
            WRITE(IOUREP,F217) LEVYRS,IMFLG,QW
          ELSE
            WRITE(IOUREP,F214) QW
          ENDIF
C
          IF((I.EQ.1.OR.I.EQ.2).AND.LEVFLG.NE.1.AND.
     *    ICY.GE.LEVYRS) THEN
            WRITE(IOUREP,F204) ZEVTOT(I),QW
          ELSE
            WRITE(IOUREP,F214) QW
          ENDIF
          IF(ALHFLG.EQ.1) WRITE(IOUREP,F218) ALHFLG
          IF(ALHFLG.NE.1) WRITE(IOUREP,F218) 2
          WRITE(IOUREP,F218) DSFLAG
C
C  Write a 1 for input mileage accumulation rates if 1 or 3
C
          IF(MYMRFG.EQ.1.OR.MYMRFG.EQ.3) WRITE(IOUREP,F218) 1
          IF(MYMRFG.EQ.2.OR.MYMRFG.EQ.4) WRITE(IOUREP,F218) 2
C
C  Write a 1 for input registration distributions if 1 or 2
C
          IF(MYMRFG.EQ.1.OR.MYMRFG.EQ.2) WRITE(IOUREP,F218) 1
          IF(MYMRFG.EQ.3.OR.MYMRFG.EQ.4) WRITE(IOUREP,F218) 2
          IF(TAMFLG.EQ.2) WRITE(IOUREP,F218) 2
          IF(TAMFLG.NE.2) WRITE(IOUREP,F218) 1
          IF(VMFLAG.EQ.1) WRITE(IOUREP,F218) 1
          IF(VMFLAG.NE.1) WRITE(IOUREP,F218) 2
          IF(NEWFLG.EQ.2.OR.NEWFLG.EQ.4.OR.NEWFLG.EQ.6)
     *      WRITE(IOUREP,F219) 2,QW
          IF(NEWFLG.EQ.1.OR.NEWFLG.EQ.3.OR.NEWFLG.EQ.5)
     *      WRITE(IOUREP,F219) 1,QW
          IF(NEWFLG.EQ.3.OR.NEWFLG.EQ.4) WRITE(IOUREP,F220)
     *      IETPYR(1),ETPPIV(1),IETPYR(2),ETPPIV(2),
     *      IETPYR(3),ETPPIV(3),IETPYR(4),ETPPIV(4),QW
          IF(NEWFLG.NE.3.AND.NEWFLG.NE.4) WRITE(IOUREP,F214) QW
C
C  Field for Stage II VRS input record, RLFLAG=2,4
C
          IF(RLFLAG.EQ.2.OR.RLFLAG.EQ.4) WRITE(IOUREP,F221)
     *      IS2SY,NPHASE,S2EFF(1),S2EFF(4),QW
          IF(RLFLAG.NE.2.AND.RLFLAG.NE.4) WRITE(IOUREP,F214) QW
C
C  Field for Stage II VRS input record, RLFLAG=3,4
C
          IF(RLFLAG.EQ.3.OR.RLFLAG.EQ.4) WRITE(IOUREP,F222)
     *      IOBMY,IVOB,QW
          IF(RLFLAG.NE.3.AND.RLFLAG.NE.4) WRITE(IOUREP,F214) QW
          IF(SPDFLG.EQ.3.OR.SPDFLG.EQ.4) WRITE(IOUREP,F223)
     *      TLVMTU,QW
          IF(SPDFLG.NE.3.AND.SPDFLG.NE.4) WRITE(IOUREP,F214) QW
          WRITE(IOUREP,F216) SCNAME,Q
   10   CONTINUE
C
C  Only print HC components if user has requested it.
C
        IF(HCFLAG.NE.1.AND.IP.EQ.1) THEN
          DO 20 IEVAP=1,5
            DO 15 I=1,10
              IV=I
              IF(I.EQ.4) IV=9
              IF(I.GT.4) IV=IV-1
              WRITE(IOUREP,F200) Q,PROJID,QW,MTHLBL,QW
              WRITE(IOUREP,F201) ICY,QW
              WRITE(IOUREP,F206) NAMEVP(IEVAP),QW
C
C  Choose between 9 vehicle names and total
C
              IF(I.NE.10) WRITE(IOUREP,F203) NAMVEH(IV),QW
              IF(I.EQ.10) WRITE(IOUREP,F203) NAMVEH(10),QW
              IF(IEVAP.EQ.1) THEN
                IF(I.NE.10) WRITE(IOUREP,F204) EFEXH(IV),QW
                IF(I.EQ.10) WRITE(IOUREP,F204) VEXH,QW
              ENDIF
C
              IF(IEVAP.EQ.2) THEN
                IF(I.LE.5.OR.I.EQ.9) WRITE(IOUREP,F204) EFEVAP(IV),QW
                IF(I.EQ.10) WRITE(IOUREP,F204) VEVAP,QW
C
C  Evap HC emissions for LDDV, LDDT and HDDV are left blank
C
                IF(I.GE.6.AND.I.LE.8) WRITE(IOUREP,F214) QW
              ENDIF
C
              IF(IEVAP.EQ.3) THEN
                IF(I.LE.5) WRITE(IOUREP,F204) EFREFL(IV),QW
                IF(I.EQ.10) WRITE(IOUREP,F204) VLOSS,QW
C
C  Refuel EF for LDDV, LDDT, HDDV, and MC are left blank
C
                IF(I.GE.6.AND.I.LE.9) WRITE(IOUREP,F214) QW
              ENDIF
C
              IF(IEVAP.EQ.4) THEN
                IF(I.LE.5) WRITE(IOUREP,F204) EFRUNL(IV),QW
                IF(I.EQ.10) WRITE(IOUREP,F204) VRUNLS,QW
C
C  Running losses for LDDV, LDDT, HDDV, and MC are left blank
C
                IF(I.GE.6.AND.I.LE.9) WRITE(IOUREP,F214) QW
              ENDIF
C
              IF(IEVAP.EQ.5) THEN
                IF(I.LE.5.OR.I.EQ.9) WRITE(IOUREP,F204) EFRSTL(IV),QW
                IF(I.EQ.10) WRITE(IOUREP,F204) VRSTLS,QW
C
C  Resting losses for LDDV, LDDT and HDDV are left blank
C
                IF(I.GE.6.AND.I.LE.8) WRITE(IOUREP,F214) QW
              ENDIF
C
C  Choose appropriate VMT, leave blank or 1
C
              IF(I.NE.4.AND.I.NE.10)
     *          WRITE(IOUREP,F204) VMTMIX(IV),QW
              IF(I.EQ.4) WRITE(IOUREP,F214) QW
              IF(I.EQ.10) WRITE(IOUREP,F204) ONE,QW
C
C  Choose appropriate speed or leave blank
C
              IF(I.EQ.4.OR.I.EQ.10) THEN
                WRITE(IOUREP,F214) QW
              ELSE
                WRITE(IOUREP,F204) SPD(IV),QW
              ENDIF
              WRITE(IOUREP,F208) TEMEXH(IP),TEMMIN,TEMMAX,QW,
     *        RVPBAS,RVPIUS,IUSESY,QW,PCCN,PCHC,PCCC,QW
              IF(IMFLAG.EQ.1) WRITE(IOUREP,F209) QW,QW
              IF(IMFLAG.EQ.2) WRITE(IOUREP,F210) ICYIM(1),ISTRIN(1),
     *          MODYR(1,1),MODYR(2,1),WAIVER(1,1),WAIVER(2,1),CRIM(1),
     *          INTYP(1),IFREQ(1),(ILDT(IVIM,1),IVIM=1,4),ITEST(1),
     *          IFILET(1),IFILET(2),CUTHC(1),CUTCO(1),CUTNO(1),QW,QW
              IF(IMFLAG.EQ.3) WRITE(IOUREP,F211) ICYIM(1),ISTRIN(1),
     *          MODYR(1,1),MODYR(2,1),WAIVER(1,1),WAIVER(2,1),CRIM(1),
     *          INTYP(1),IFREQ(1),(ILDT(IVIM,1),IVIM=1,4),ITEST(1),
     *          IFILET(1),IFILET(2),CUTHC(1),CUTCO(1),CUTNO(1),QW,
     *          ICYIM(2),ISTRIN(2),
     *          MODYR(1,2),MODYR(2,2),WAIVER(1,2),WAIVER(2,2),CRIM(2),
     *          INTYP(2),IFREQ(2),(ILDT(IVIM,2),IVIM=1,4),ITEST(2),
     *          IFILET(1),IFILET(2),CUTHC(2),CUTCO(2),CUTNO(2),QW
C
C  Choose when to print ATP record
C
              IF(ATPFLG.EQ.1) WRITE(IOUREP,F214) QW
              IF(ATPFLG.EQ.2) WRITE(IOUREP,F212)
     *          LAPSY,LAP1ST,LAPLST,LVTFLG,ATPPGM,ATPFQT,CRATP,DISTYP,
     *          QW
              IF(PRSFLG.EQ.1) WRITE(IOUREP,F214) QW
              IF(PRSFLG.EQ.2) WRITE(IOUREP,F213)
     *          LPRSY,LPR1ST,LPRLST,LVPRS,PRSPGM,PRSFQT,CRPRS,QW
              IF(PRGFLG.EQ.1) WRITE(IOUREP,F214) QW
              IF(PRGFLG.EQ.2) WRITE(IOUREP,F213)
     *          LPGSY,LPG1ST,LPGLST,LVPRG,PRGPGM,PRGFQT,CRPRG,QW
              WRITE(IOUREP,F202) NOYES(RFGFLG),QW
C
C  Only write Oxy record when OXYFLG=2, provided that RFGFLG=1,
C  or when OXYFLG=2 and when RFGFLG=2 and SEAFLG=2
C
              IF(OXYFLG.EQ.2.AND.RFGFLG.EQ.2.AND.SEAFLG.EQ.2)
     *          WRITE(IOUREP,F215) SHRMKT(2),SHRMKT(3),
     *          OXYCNT,NOYES(IGASHW),QW
              IF(OXYFLG.EQ.2.AND.RFGFLG.EQ.2.AND.SEAFLG.EQ.1)
     *          WRITE(IOUREP,F214) QW
              IF(OXYFLG.EQ.2.AND.RFGFLG.EQ.1)
     *          WRITE(IOUREP,F215) SHRMKT(2),SHRMKT(3),
     *          OXYCNT,NOYES(IGASHW),QW
              IF(OXYFLG.EQ.1)
     *          WRITE(IOUREP,F214) QW
              WRITE(IOUREP,F202) NAMREG(IREJN),QW
C
C  Only write Tier 1(I/M), LEV (I/M) or ZEV% when IV=1,2, and
C  ICY is greater than or equal to 1994.  And only write
C  Tier 1(I/M) when there is an I/M.  And only write LEV (I/M)
C  and ZEV% when there is an LEV (ICY is greater than or equal
C  to LEVYRS and LEVFLG=2).
C
              IF((I.EQ.1.OR.I.EQ.2).AND.IMFLAG.NE.1.AND.
     *        ICY.GE.IY1994) THEN
                WRITE(IOUREP,F202) NOYES(IMTIER+1),QW
              ELSE
                WRITE(IOUREP,F214) QW
              ENDIF
C
              IF((I.EQ.1.OR.I.EQ.2).AND.LEVFLG.NE.1.AND.
     *        ICY.GE.LEVYRS) THEN
                WRITE(IOUREP,F217) LEVYRS,IMFLG,QW
              ELSE
                WRITE(IOUREP,F214) QW
              ENDIF
C
              IF((I.EQ.1.OR.I.EQ.2).AND.LEVFLG.NE.1.AND.
     *        ICY.GE.LEVYRS) THEN
                WRITE(IOUREP,F204) ZEVTOT(I),QW
              ELSE
                WRITE(IOUREP,F214) QW
              ENDIF
              IF(ALHFLG.EQ.1) WRITE(IOUREP,F218) ALHFLG
              IF(ALHFLG.NE.1) WRITE(IOUREP,F218) 2
              WRITE(IOUREP,F218) DSFLAG
C
C  Write a 1 for input mileage accumulation rates if 1 or 3
C
              IF(MYMRFG.EQ.1.OR.MYMRFG.EQ.3) WRITE(IOUREP,F218) 1
              IF(MYMRFG.EQ.2.OR.MYMRFG.EQ.4) WRITE(IOUREP,F218) 2
C
C  Write a 1 for input registration distributions if 1 or 2
C
              IF(MYMRFG.EQ.1.OR.MYMRFG.EQ.2) WRITE(IOUREP,F218) 1
              IF(MYMRFG.EQ.3.OR.MYMRFG.EQ.4) WRITE(IOUREP,F218) 2
              IF(TAMFLG.EQ.2) WRITE(IOUREP,F218) 2
              IF(TAMFLG.NE.2) WRITE(IOUREP,F218) 1
              IF(VMFLAG.EQ.1) WRITE(IOUREP,F218) 1
              IF(VMFLAG.NE.1) WRITE(IOUREP,F218) 2
              IF(NEWFLG.EQ.2.OR.NEWFLG.EQ.4.OR.NEWFLG.EQ.6)
     *          WRITE(IOUREP,F219) 2,QW
              IF(NEWFLG.EQ.1.OR.NEWFLG.EQ.3.OR.NEWFLG.EQ.5)
     *          WRITE(IOUREP,F219) 1,QW
              IF(NEWFLG.EQ.3.OR.NEWFLG.EQ.4) WRITE(IOUREP,F220)
     *          IETPYR(1),ETPPIV(1),IETPYR(2),ETPPIV(2),
     *          IETPYR(3),ETPPIV(3),IETPYR(4),ETPPIV(4),QW
              IF(NEWFLG.NE.3.AND.NEWFLG.NE.4) WRITE(IOUREP,F214) QW
C
C  Field for Stage II VRS input record, RLFLAG=2,4
C
              IF(RLFLAG.EQ.2.OR.RLFLAG.EQ.4) WRITE(IOUREP,F221)
     *          IS2SY,NPHASE,S2EFF(1),S2EFF(4),QW
              IF(RLFLAG.NE.2.AND.RLFLAG.NE.4) WRITE(IOUREP,F214) QW
C
C  Field for Stage II VRS input record, RLFLAG=3,4
C
              IF(RLFLAG.EQ.3.OR.RLFLAG.EQ.4) WRITE(IOUREP,F222)
     *          IOBMY,IVOB,QW
              IF(RLFLAG.NE.3.AND.RLFLAG.NE.4) WRITE(IOUREP,F214) QW
              IF(SPDFLG.EQ.3.OR.SPDFLG.EQ.4) WRITE(IOUREP,F223)
     *          TLVMTU,QW
              IF(SPDFLG.NE.3.AND.SPDFLG.NE.4) WRITE(IOUREP,F214) QW
              WRITE(IOUREP,F216) SCNAME,Q
   15       CONTINUE
   20     CONTINUE
        ENDIF
   25 CONTINUE
C
      RETURN
      END
