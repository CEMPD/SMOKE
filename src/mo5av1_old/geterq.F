C Version @(#)$Id$ $Source$ $Date$ 

      SUBROUTINE GETERQ(INERR,*)
C
C  GETERQ reads in, performs range checks on, converts to gpm if IV = HDV,
C  and stores replacement emission factor rates & their associated parameters.
C
C  Called by ONESEC.
C
C  Calls HDCGPM, QUITER, STOERQ and YRTEST.
C
C  Input on call:
C
C    parameter list: INERR
C    common blocks:
C    /BASEQ4/ MAXCT
C    /FLAGS1/ PROMPT
C    /IOUCOM/ IOUGEN,IOUASK
C    /YEARS4/ IY1941,IY2020
C
C  Output on return:
C
C    parameter list: INERR
C    common blocks:
C    /BASEQ3/ KINKS
C    /BASEQ4/ NEWPAR,BERNEW,NEWFIT,NEWCT
C
C  Local variable / array dictionary:
C
C   Name   Type              Description
C  ------  ----  -------------------------------------------------------
C  A50NEW   R    above 50K miles deterioration rate, nonzero for 1981+ LDGV
C                and 1991+ LDGT, HC/CO cases only (local variable: not used in
C                HDCGPM)
C  KINK50   I    kinked at 50K miles ef curve flag: 1 = no, 2 = yes
C  NEWPOL   I    new emission rate's pollutant code
C  NEWREG   I    new emission rate's region code
C
C  Notes:
C
C  Kinked EFs for LDGTs were added for MOBILE5.
C
C
      INTEGER PROMPT,TAMFLG,SPDFLG,VMFLAG,OXYFLG,DSFLAG
C
      COMMON /BASEQ3/ NEWVEH,NEWMYF,NEWMYL,ZMLNEW,DRNEW,KINKS
      COMMON /BASEQ4/NEWPAR(5,100),BERNEW(3,100),NEWFIT(100),NEWCT,MAXCT
      COMMON /FLAGS1/ PROMPT,TAMFLG,SPDFLG,VMFLAG,OXYFLG,DSFLAG
      COMMON /IOUCOM/ IOUIMD,IOUGEN,IOUREP,IOUERR,IOUASK
      COMMON /YEARS4/ IY1941,IY1960,IY1994,IY2020
C
C  Read in number of new ef rate records to be processed.
C
      IF(PROMPT.EQ.2) WRITE(IOUASK,200)
  200 FORMAT('&','Enter no. of new e.f. intercept/slope pairs (I3):')
      READ(IOUGEN,100,ERR=99,END=99) NEWCT
  100 FORMAT(I3)
      IF(NEWCT.LT.1.OR.NEWCT.GT.MAXCT) CALL QUITER(0.,NEWCT,23,INERR)
C
C  Read in and process one record at a time.
C
      KINKS=1
      DO 10 IFAC=1,NEWCT
      IF(PROMPT.EQ.2) WRITE(IOUASK,210)
  210 FORMAT(' ','Enter region,veh. type,pollutant,first my,'/
     *  ' ','last my,intercept,slope & above 50K miles slope:')
      READ(IOUGEN,110,ERR=99,END=99) NEWREG,NEWVEH,NEWPOL,NEWMYF,NEWMYL,
     *                        ZMLNEW,DRNEW,A50NEW
  110 FORMAT(I1,2(1X,I1),2(1X,I2),3(1X,F6.2))
      IF(NEWREG.LT.1.OR.NEWREG.GT.2) CALL QUITER(0.,NEWREG,62,INERR)
      IF(NEWVEH.LT.1.OR.NEWVEH.GT.8) CALL QUITER(0.,NEWVEH,25,INERR)
      IF(NEWPOL.LT.1.OR.NEWPOL.GT.3) CALL QUITER(0.,NEWPOL,26,INERR)
      CALL YRTEST(NEWMYF,27,IY1941,IY2020,INERR)
      CALL YRTEST(NEWMYL,27,IY1941,IY2020,INERR)
      IF(NEWMYL.LT.NEWMYF) CALL QUITER(0.,0,29,INERR)
      KINK50=1
      IF((NEWVEH.LE.3).AND.
     *   (NEWMYF.GE.1981.OR.NEWMYL.GE.1981)) KINK50=2
      IF(KINK50.EQ.2) KINKS=2
      IF(ZMLNEW.LE.0.0) CALL QUITER(ZMLNEW,0,30,INERR)
      IF(DRNEW.LT.0.0) CALL QUITER(DRNEW,0,31,INERR)
      IF(KINK50.EQ.2.AND.A50NEW.LT.0.0) CALL QUITER(A50NEW,0,31,INERR)
C
C  If the emission rate is for a heavy duty vehicle class, convert the rate's
C  units from grams per brakehorsepower to grams per mile.
C
      IF(NEWVEH.EQ.4.OR.NEWVEH.EQ.7) CALL HDCGPM
C
C  Save screened data in temporary storage areas.
C
      NEWPAR(1,IFAC)=NEWREG
      NEWPAR(2,IFAC)=NEWVEH
      NEWPAR(3,IFAC)=NEWPOL
      NEWPAR(4,IFAC)=NEWMYF
      NEWPAR(5,IFAC)=NEWMYL
      BERNEW(1,IFAC)=ZMLNEW
      BERNEW(2,IFAC)=DRNEW
      IF(KINK50.EQ.1) BERNEW(3,IFAC)=-1.0
      IF(KINK50.EQ.2) BERNEW(3,IFAC)=A50NEW
   10 CONTINUE
C
C  Store new ef equations for later access by BEF.
C
      CALL STOERQ(INERR)
C
      RETURN
C
C  Error or EOF on any attempted read => take alternate
C  return 1 => run aborts.
C
   99 CALL QUITER(0.,0,128,INERR)
C
      RETURN1
      END
