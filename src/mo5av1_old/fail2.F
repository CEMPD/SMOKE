C Version @(#)$Id$ $Source$ $Date$ 

      FUNCTION FAIL2(MYR,KDX,IV,IK)
C
C  FAIL2 handles the effect of Purge and Pressure failures as
C  well as evap system tampering, returning a fraction for
C  purge/pressure for evaluation years less than 1992.
C
C  Called by FAIL.
C
C  Calls ENFORC, ITAMPT, CHKFUN, and CHKATP.
C
C  Input on call:
C
C    parameter list: MYR,KDX,IV,IK
C    common blocks:
C    /CUMCOM/ CUMMIL
C    /EVPDAT/ PCTEQP
C    /IM240P/ DSIZE
C    /LOOKUP/ IVTAM
C    /MAXIMA/ MAXYRS
C    /MYCODE/ MY
C    /PRGCH2/ PRGPGM,PRGFQT,CRPRG
C    /PRSCH2/ PRSPGM,PRSFQT,CRPRS
C    /TAMOUT/ TDU
C
C  Output on return:
C
C    function: FAIL2
C
C  Local variable / array dictionary:
C
C   Name   Type              Description
C  ------  ----  -------------------------------------------------------
C  ADJ      R    Fraction of failures properly repaired.
C  DEFF     R    Discount of inspection effectiveness for decentralized
C                inspections.
C  DPRES    R    Pressure check discount.
C  DPURG    R    Purge check discount.
C  EADJ     R    Fraction of unrepaired failures remaining.
C  EPRES    R    EADJ value for the presure check.
C  EPURG    R    EADJ value for the purge check.
C  FACTOR   R    Fraction of tampered evap systems passing the purge
C                check.
C  FDR      R    Increase in the number of failures for each 10,000
C                miles.
C  FMILE    R    Accumulated mileage between zero and 50,000 miles.
C  FMILE5   R    Accumulated mileage greater than 50,000 miles.
C  FPRES    R    Combined EPRES and DPRES effect.
C  FPURG    R    Combined EPURG and DPURG effect.
C  FZML     R    Number of failures at zero mileage.
C  I        I    Loop counting variable, through 5 parts of Venn diag.
C  ICY      I    See top of program.
C  IDX      I    See top of program.
C  IHG      I    Indicates non-exhaust technology (for ITAMPT) where:
C                1=evaporative, 2=crankcase.
C  IVTAM    I    See top of program.
C  JCASE    I    Indicates the presence of functionality checks
C                for Model Year, Veh. and ICY where:
C                1=neither Pressure or Purge Checks
C                2=Pressure Checks not Purge Checks
C                3=Purge Checks not Pressure Checks
C                4=both Pressure and Purge Checks
C  JDX      I    See top of program.
C  KCASE    I    Indicates the presence of Anti-tampering Prog.
C                for Model Year, Veh. and ICY where:
C                1=no Anti-tampering Program
C                2=Anti-tampering Program exists without EVAP
C                3=Anti-tampering Program exists with EVAP
C  MY       I    See top of program.
C  OVLP     R    Initial size of pressure/purge overlap categories.
C  OVLR     R    Final size of pressure/purge overlap categories after
C                program effects.
C  PCTEQP   R    Fraction of the model year equipped with evaporative
C                control systems.
C  PRGR     R    Fraction of the model year fleet with purge system
C                failures.
C  PRSR     R    Fraction of the model year fleet with pressure system
C                failures.
C  PPEFF    R    Reduction in failures due to annual or biennial
C                inspection at each age.
C  PPRT     R    Fraction of the model year fleet with either purge or
C                pressure failures.
C  SUM      R    Maximum combined pressure and purge rate.
C  TEFF     R    Effect of ATP on overlaps 4 and 5
C  TRAT     R    Fraction of the model year with tampered evaporative
C                control systems.
C
C
C  Notes.
C
C  FAIL2 was added for MOBILE4.1 and was created from FAIL in
C  MOBILE5.
C  FAIL2 was modified in MOBILE5 to include ZML/DR.
C  FAIL2 was totally rewritten for MOBILE5 to call CHKFUN and
C  CHKATP.  CHKATP indicates that ATP either exists or exists
C  without EVAP: this critical information is only temporary.
C  FAIL2 does not incorporate changes that FAIL1 does.
C
C
      INTEGER PRGPGM,PRGFQT,PRSPGM,PRSFQT
      INTEGER PRGFLG,PRSFLG,ATPFLG
C
      COMMON /CUMCOM/ CUMMIL(25,8)
      COMMON /FLAGS3/ ATPFLG,TPDFLG,RLFLAG,LOCFLG,TEMFLG,OUTFMT
      COMMON /EVPDAT/ EVPAER(3,2),EVPTGS(3,2,4),CCEMI(7,4)
      COMMON /IM240P/ DSIZE(25,4)
      COMMON /LOOKUP/ IVTAM,IQG,IPG,JPGD,IHG,IGCSF
      COMMON /MAXIMA/ MAXVEH,MAXLTW,MAXPOL,MAXREG,MAXYRS
      COMMON /MYCODE/ MY,IDX,JDX,LDXSY,LMYRVT,IAY,IMDXSY,IMKINK
      COMMON /PRGCH2/ PRGPGM,PRGFQT,CRPRG,PRGFLG
      COMMON /PRSCH2/ PRSPGM,PRSFQT,CRPRS,PRSFLG
      COMMON /TAMOUT/ TAMBAG(3,3,25,4),THS(2,25,4),TDU(25,4),TCC(25,4)
C
      DIMENSION FZML(3),FDR(3,2),OVLP(5),OVLR(5),PPEFF(25,2)
C
C  Purge/Pressure Check Effectiveness, Annual, Biennial
C  For Mobile5, based on 95% detection,
C    .0649 ZML, 0.00096 & 0.03325 DRs
C  8/18/92...
C
      DATA PPEFF /
     * .187,.412,.947,.916,.740,.724,.787,.831,.862,.882,
     * .898,.911,.958,.937,.920,.925,.933,.939,.943,.948,
     * .952,.980,.964,.961,.960,
     * .187,.412,.947,.916,.727,.723,.575,.822,.711,.875,
     * .788,.906,.875,.933,.908,.888,.814,.899,.834,.914,
     * .858,.951,.926,.975,.953 /
C
C  Purge, Pressure and Purge/Pressure Failure ZMls
C
      DATA FZML/0.0265, 0.0474, 0.0649/
C
C  Purge, Pressure and Purge/Pressure Failure DRs
C  for 0-50,000 miles then 50,000+ miles.
C
      DATA FDR/0.00183, 0.00000, 0.00096,
     *         0.01813, 0.01913, 0.03325/
C
      DATA FACTOR,DEFF / 0.50, 0.50 /
C
C  Initialize variables
C
      DPRES=DEFF
      DPURG=DEFF
      FAIL2=0.0
      PRGR=0.0
      PRSR=0.0
      PPRT=0.0
      DO 10 I=1,5
        OVLP(I)=0.0
        OVLR(I)=0.0
   10 CONTINUE
C
C  For diesel & motorcycles, return 0.0 for FAIL2
C
      IF(IV.GT.4) RETURN
C
C  Determine percent of vehicles equipped with evap systems
C  (Use KDX and MYR to avoid problems with common block MYCODE)
C
      MY=MYR
      IDX=KDX
C
      IHG=1
      IVTAM=IV
      IG3=ITAMPT(3)
      PCTEQP=EVPTGS(IG3,1,IV)
C
C  Determine Mileage and Pointers
C
      JDX=MAXYRS+1-IDX
C
      FMILE=CUMMIL(JDX,IV)/10000.
      FMILE5=0.0
C
      IF(FMILE.GT.5.0) THEN
        FMILE5=FMILE-5.0
        FMILE=5.0
      ENDIF
C
C  Determine overall Pressure and Purge failure rates
C
      PRGR=(FZML(1)+FMILE*FDR(1,1)+FMILE5*FDR(1,2))*PCTEQP
      PRSR=(FZML(2)+FMILE*FDR(2,1)+FMILE5*FDR(2,2))*PCTEQP
      PPRT=(FZML(3)+FMILE*FDR(3,1)+FMILE5*FDR(3,2))*PCTEQP
      TRAT=DSIZE(IDX,IV)*PCTEQP
C
C  Set the Effects of ATP (TEFF) to 1.0 then compute it
C
      TEFF=1.0
      IF(ATPFLG.EQ.2.AND.TRAT.GT.0.0) TEFF=TDU(IDX,IV)/TRAT
C
      IF(TRAT.GT.PRSR) TRAT=PRSR
      SUM=PRGR+PRSR
      IF(PPRT.GT.SUM) PPRT=SUM
C
C  Determine proportion of overall in each overlap
C
      OVLP(4)=TRAT*FACTOR
      OVLP(5)=TRAT*(1.0-FACTOR)
      OVLP(1)=PPRT-PRGR-OVLP(4)
C
      IF(OVLP(1).LT.0.0) OVLP(1)=0.0
C
      OVLP(2)=PRSR-OVLP(1)-OVLP(4)-OVLP(5)
C
      IF(OVLP(2).LT.0.0) THEN
         OVLP(2)=0.0
         OVLP(5)=PRGR-(PPRT-PRSR)
         OVLP(4)=OVLP(4)+TRAT*(1.0-FACTOR)-OVLP(5)
      ENDIF
C
      OVLP(3)=PRGR-OVLP(2)-OVLP(5)
C
      IF(OVLP(3).LT.0.0) OVLP(3)=0.0
C
C  Initialize rates in overlaps
C
      DO 20 I=1,5
        OVLR(I)=OVLP(I)
   20 CONTINUE
C
C  ICY is needed for calls to CHKFUN and CHKATP.
C
      ICY=JDX-1+MY
C
C  One means none, two means pressure only, three is purge only
C  and four is both.
C
      JCASE=CHKFUN(MY,IVTAM,ICY)
      KCASE=CHKATP(MY,IVTAM,ICY)
      GOTO(30,40,50,60),JCASE
C
C  Neither Pressure or Purge is in effect
C  Anti-tampering program effects.
C
   30 IF(KCASE.GT.2) THEN
         OVLR(4)=OVLR(4)*TEFF
         OVLR(5)=OVLR(5)*TEFF
      ENDIF
      GOTO 70
C
C  Pressure Check Alone
C  Adjust rates to reflect programs
C
   40 EADJ=1.0-ENFORC(CRPRS,1)*PPEFF(JDX,PRSFQT)
      ADJ=1.0-EADJ
      IF(PRSPGM.EQ.1) DEFF=0.0
C
      OVLR(1)=DEFF*ADJ*OVLP(1)+EADJ*OVLP(1)
      OVLR(2)=DEFF*ADJ*OVLP(2)+EADJ*OVLP(2)
      OVLR(4)=DEFF*ADJ*OVLP(4)+EADJ*OVLP(4)
      OVLR(5)=DEFF*ADJ*OVLP(5)+EADJ*OVLP(5)
C
      OVLR(3)=OVLP(3)+OVLP(2)-OVLR(2)+OVLP(5)-OVLR(5)
      GOTO 70
C
C  Purge Check Alone
C  Adjust rates to reflect programs
C
   50 EADJ=1.0-ENFORC(CRPRG,1)*PPEFF(JDX,PRGFQT)
      ADJ=1.0-EADJ
      IF(PRGPGM.EQ.1) DEFF=0.0
C
      OVLR(2)=DEFF*ADJ*OVLP(2)+EADJ*OVLP(2)
      OVLR(3)=DEFF*ADJ*OVLP(3)+EADJ*OVLP(3)
      OVLR(5)=DEFF*ADJ*OVLP(5)+EADJ*OVLP(5)
C
      OVLR(1)=OVLP(1)+OVLP(2)-OVLR(2)
      OVLR(4)=OVLP(4)+OVLP(5)-OVLR(5)
C
C Anti-tampering program effects
C
      IF(KCASE.GT.2) OVLR(4)=OVLR(4)*TEFF
      GOTO 70
C
C  Combined Purge/Pressure Checks
C
   60 EPRES=1.0-ENFORC(CRPRS,1)*PPEFF(JDX,PRSFQT)
      IF(PRSPGM.EQ.1) DPRES=0.0
      FPRES=EPRES+DPRES-DPRES*EPRES
C
      EPURG=1.0-ENFORC(CRPRG,1)*PPEFF(JDX,PRGFQT)
      IF(PRGPGM.EQ.1) DPURG=0.0
      FPURG=EPURG+DPURG-DPURG*EPURG
C
      OVLR(1)=FPRES*OVLP(1)
      OVLR(3)=FPURG*OVLP(3)
      OVLR(4)=FPRES*OVLP(4)
C
      IF(FPRES.EQ.FPURG) THEN
         OVLR(2)=FPRES*OVLP(2)
         OVLR(5)=FPRES*OVLP(5)
      ENDIF
C
C     Assume compliance cars are the same.
C
      IF(FPRES.LT.FPURG) THEN
         OVLR(2)=FPURG*OVLP(2)
         OVLR(3)=OVLR(3)+OVLP(2)*(FPURG-FPRES)
     *                  +OVLP(5)*(FPURG-FPRES)
         OVLR(5)=FPURG*OVLP(5)
      ENDIF
C
      IF(FPRES.GT.FPURG) THEN
         OVLR(1)=OVLR(1)+OVLP(2)*(FPRES-FPURG)
         OVLR(2)=FPRES*OVLP(2)
         OVLR(4)=OVLR(4)+OVLP(5)*(FPRES-FPURG)
         OVLR(5)=FPRES*OVLP(5)
      ENDIF
C
C  Return adjusted overall Pressure/Purge rate
C
C  IK = 1 : Purge only rate
C       2 : Pressure rate
C
   70 IF(IK.EQ.1) FAIL2=OVLR(3)
      IF(IK.EQ.2) FAIL2=OVLR(1)+OVLR(2)+OVLR(4)+OVLR(5)
C
      RETURN
      END
