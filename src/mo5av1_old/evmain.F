C Version @(#)$Id$ $Source$ $Date$ 

      SUBROUTINE EVMAIN(MY,KEYEVP,KEYEQ,IVGAS,ISTD,RIVAL,EVRATE,IY,ICY)
C
C  EVMAIN is the driver "main" for the block of functions that calculate
C  hot soak or diurnal rates, given the vehicle type, emission standard,
C  base/in-use/certification RVP/Index and the untampered/uncontrolled/cert
C  equation group specified in the parameter list.
C
C  Evap rates are adjusted for high altitude regions (IREJN=2) as follows:
C     Hot Soak High Alt Rate = Hot Soak Low Alt Rate * HIADJ
C     Diurnal  High Alt Rate = Diurnal  Low Alt Rate * HIADJ
C  EXCEPT:
C        LDGVs   1972-76my   Actual LDGV Test Data
C                1977my      Same as Low Alt
C                1984+my     Same as Low Alt
C        LDGT1s  1972-76my   Same as LDGVs
C                1977my      Same as Low Alt
C
C  Called by DIURNL and HOTSOK.
C
C  Calls EVPRE, EVSTD, EV81, EV85, EVUNC and EVHI.
C
C  Input on call:
C
C    parameter list: MY,KEYEVP,KEYEQ,IVGAS,ISTD,RIVAL
C    common blocks:
C    /EVAHIA/ HIADJ
C    /REGION/ IREJN
C
C  Output on return:
C
C    parameter list: EVRATE
C
C  Notes:
C
C  Certification fuels support (KEYEQ=3) from MOBILE4 was removed in MOBILE4.1.
C  Modified for MOBILE4.1 to add support for pressure/purge fail.
C
C
      COMMON /EVAHIA/ HIHS(3,2),HIADJ
      COMMON /REGION/ FEET(2),IREJN,ALT,INITPR
C
      JSTD=ISTD
C
C  Function call depends
C      first on the equation group key,
C            second on the vehicle type and
C                   third on the emission standard.
C
      GOTO (10,55,55), KEYEQ
C
C  KEYEQ = 1 = passing rates
C
   10 GOTO (15,15,35,50,25), IVGAS
C
C  LDGV/LDGT1
C
   15 GOTO (20,25,25,25,30,30,30), ISTD
C
   20 EVRATE=EVPRE(KEYEVP,KEYEQ,JSTD,IVGAS,RIVAL)
      GOTO 90
   25 EVRATE=EVSTD(KEYEVP,KEYEQ,JSTD,IVGAS,RIVAL,IY,ICY)
      IF(IVGAS.EQ.5) GOTO 100
      IF(IREJN.EQ.1.OR.MY.EQ.1977) GOTO 100
      IF(ISTD.EQ.2.OR.ISTD.EQ.4) GOTO 90
      IF(ISTD.EQ.3) EVRATE=EVHI(KEYEVP,KEYEQ,RIVAL)
      GOTO 100
   30 IFDS=ISTD-4
      EVRATE=EV81(KEYEVP,IFDS,IVGAS,RIVAL,IY,ICY)
      GOTO 90
C
C LDGT2
C
   35 GOTO (20,40,45,45,45), ISTD
   40 EVRATE=EVSTD(KEYEVP,KEYEQ,JSTD,IVGAS,RIVAL,IY,ICY)
      GOTO 90
C
   45 IFDS=ISTD-2
      EVRATE=EV81(KEYEVP,IFDS,IVGAS,RIVAL,IY,ICY)
      GOTO 90
C
C HDGV
C
   50 IF(ISTD.EQ.1) GOTO 20
      EVRATE=EV85(KEYEVP,KEYEQ,IVGAS,RIVAL,IY,ICY)
      GOTO 90
C
C  KEYEQ = 2 = failed purge uncontrolled emissions
C  KEYEQ = 3 = failed pressure uncontrolled emissions
C
   55 GOTO (60,60,70,50,25), IVGAS
C
C  Failed LDGV/LDGT1
C
   60 GOTO (20,25,25,25,80,80,80), ISTD
C
C  Failed LDGT2
C
   70 GOTO (20,40,80,80,80), ISTD
C
C  Failed 1981+ LDGV/T1/T2
C
   80 EVRATE=EVUNC(KEYEVP,KEYEQ,IVGAS,RIVAL)
   90 IF(IREJN.EQ.1) GOTO 100
      IF(IVGAS.LE.2.AND.MY.EQ.1977) GOTO 100
      IF(IVGAS.EQ.1.AND.MY.GE.1984) GOTO 100
      EVRATE=EVRATE*HIADJ
C
C  Limit Hot Soaks to 0-45 grams, inclusive, (low alt).
C  Limit Diurnals to 0-100 grams, inclusive, (low alt).
C  Account for hi-altitude accordingly.
C
  100 IF(KEYEVP.EQ.1) EVMAX=45.0
      IF(KEYEVP.EQ.2) EVMAX=100.0
C
      IF(IREJN.EQ.2) EVMAX=EVMAX*HIADJ
C
      EVRATE=AMIN1(AMAX1(EVRATE,0.0),EVMAX)
C
      RETURN
      END
