C Version "@(#)$Id$ $Source$ $Date$ 

      FUNCTION FAIL(MYR,KDX,IV,IK)
C
C  FAIL handles the effect of Purge and Pressure failures as
C  well as evap system tampering, returning a fraction for
C  purge/pressure.
C
C  Called by HCCALX.
C
C  Calls ENFORC, ITAMPT, CHKFUN, and CHKATP.
C
C  Input on call:
C
C    parameter list: MYR,KDX,IV,IK
C    common blocks:
C    /CUMCOM/ CUMMIL
C    /EVPDAT/ PCTEQP
C    /IMPAR5/ PPGUEF
C    /IM240P/ DSIZE
C    /LOOKUP/ IVTAM
C    /MAXIMA/ MAXYRS
C    /MYCODE/ MY
C    /PPEFFD/ PPEFF
C    /PRGCH1/ LPGSY,LPG1ST,LPGLST,LVPRG
C    /PRGCH2/ PRGPGM,PRGFQT,CRPRG,PRGFLG
C    /PRSCH1/ LPRSY,LPR1ST,LPRLST,LVPRS
C    /PRSCH2/ PRSPGM,PRSFQT,CRPRS,PRSFLG
C    /TAMOUT/ TDU
C
C  Output on return:
C
C    function: FAIL
C
C  Local variable / array dictionary:
C
C   Name   Type              Description
C  ------  ----  ------------------------------------------------
C  ADJ      R    Fraction of failures properly repaired.
C  DPRES    R    Pressure check discount.
C  DPURG    R    Purge check discount.
C  EADJ     R    Fraction of unrepaired failures remaining.
C  EPRES    R    EADJ value for the presure check.
C  EPURG    R    EADJ value for the purge check.
C  FACTOR   R    Fraction of tampered evap systems passing the
C                purge check.
C  FDR      R    Increase in the number of failures for each
C                10,000 miles.
C  FMILE    R    Accumulated mileage between zero and
C                50,000 miles.
C  FMILE5   R    Accumulated mileage greater than 50,000 miles.
C  FPRES    R    Combined EPRES and DPRES effect.
C  FPURG    R    Combined EPURG and DPURG effect.
C  FZML     R    Number of failures at zero mileage.
C  I        I    Loop counting variable, through 5 parts of Venn
C                diagram.
C  ICY      I    See top of program.
C  IDX      I    See top of program.
C  IHG      I    Indicates non-exhaust technology (for ITAMPT)
C                where: 1=evaporative, 2=crankcase.
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
C  OVLP     R    Initial size of pressure/purge overlap
C                categories.
C  OVLR     R    Final size of pressure/purge overlap categories after
C                program effects.
C  PCTEQP   R    Fraction of the model year equipped with
C                evaporative control systems.
C  PRGPIN   R    Correction factor for biennial phase in of
C                one year old purge check programs
C  PRGR     R    Fraction of the model year fleet with purge
C                system failures.
C  PRSPIN   R    Correction factor for biennial phase in of
C                one year old pressure check programs
C  PRSR     R    Fraction of the model year fleet with pressure
C                system failures.
C  PPEFF    R    Reduction in failures due to annual or biennial
C                inspection at each age.
C  PPRT     R    Fraction of the model year fleet with either
C                purge or pressure failures.
C  SUM      R    Maximum combined pressure and purge rate.
C  TEFF     R    Effect of ATP on overlaps 4 and 5
C  TRAT     R    Fraction of the model year with tampered
C                evaporative control systems.
C
C
C  Notes.
C
C  FAIL was added for MOBILE5 corrected and created from FAIL
C  and the old code is maintained in FAIL2 for evaluation years
C  less than 1992.
C  8-April-1993 @ ARC-bsg Subtask 212
C  to handle the biennial phase in of Pressure and Purge
C  15-June-1993 @ ARC-bwk Subtask 245
C  to fix a problem with the initialization of local variables
C  FAIL was modified in MOBILE5 to include ZML/DR.
C  FAIL was totally rewritten for MOBILE5 to call CHKFUN and
C  CHKATP. CHKATP indicates that ATP either exists or exists
C  without EVAP: this critical information is only temporary.
C  FAIL was modified to fix a bug in release MOBILE5, that was
C  causing evaporative emissions to go up when newer vehicles
C  were in use. This problem had been identified on
C  December 24, 1992.
C  FAIL was modified again for MOBILE5 corrected.  It was turned
C  into a driver subroutine and returns the value of FAIL1 or
C  FAIL2.
C  7-July-1994 @ CSC-tmm Subtask 438
C  Deleted FAIL2. Renamed FAIL1 back to FAIL. This eliminated
C  the old code for evaluation years of less than 1992.
C  Sep-27-1994 @ CSC-bsg Request 438 removed extra ICY definition
C
C  April-1995 @ CSC-zp Req. 2-495 "Add MOB5a-H to MOBILE5b"
C  Added: PPEFF is read in DRIVER as file PPEFF.D instead of
C  locally in FAIL.
C  23 Feb 96 @DynTel-MLA 2-611  Added include file IMPAR5 and removed the
C     local variable DEFF and the DATA statement that initialized it. The
C     code now uses PPGUEF, defined in COMMON area IMPAR5, which may 
C     contain a user-supplied value for the effectiveness of the purge/
C     pressure programs. The default is still 0.50 if the user does not 
C     supply a new value.
C
      INCLUDE 'CUMCOM.I'
      INCLUDE 'FLAGS3.I'
      INCLUDE 'EVPDAT.I'
      INCLUDE 'IMPAR5.I'
      INCLUDE 'IM240P.I'
      INCLUDE 'LOOKUP.I'
      INCLUDE 'MAXIMA.I'
      INCLUDE 'MYCODE.I'
      INCLUDE 'PPEFFD.I'
      INCLUDE 'PRGCH1.I'
      INCLUDE 'PRGCH2.I'
      INCLUDE 'PRSCH1.I'
      INCLUDE 'PRSCH2.I'
      INCLUDE 'TAMOUT.I'
C
      DIMENSION FZML(3),FDR(3,2),OVLP(5),OVLR(5)
C
C  Purge/Pressure Check Effectiveness, Annual, Biennial
C  For Mobile5, based on 95% detection,
C    .0649 ZML, 0.00096 & 0.03325 DRs
C
      DATA FZML/0.0265, 0.0474, 0.0649/
C
C  Purge, Pressure and Purge/Pressure Failure DRs
C  for 0-50,000 miles then 50,000+ miles.
C
      DATA FDR/0.00183, 0.00000, 0.00096,
     *         0.01813, 0.01913, 0.03325/
C
      DATA FACTOR / 0.50 /
C
C  Initialize variables. DPRES and DPURG are now set directly from the
C  COMMON variable PPGUEF (instead of from the local variable DEFF).
C
      DPRES=PPGUEF
      DPURG=PPGUEF
      FAILOC=0.0
      PRSPIN=1.0
      PRGPIN=1.0
      PRGR=0.0
      PRSR=0.0
      PPRT=0.0
      IF(IV.GT.4) GOTO 98
      DO 10 I=1,5
        OVLP(I)=0.0
        OVLR(I)=0.0
   10 CONTINUE
C
C  Sep-27-1994 @ CSC-bsg Request 438 removed extra ICY definition
C  (Use KDX and MYR to avoid problems with common block MYCODE)
C
      MY=MYR
      IDX=KDX
      ICY=MYR-IDX+MAXYRS
      IF(PRSFLG.EQ.2.AND.ICY.EQ.LPRSY+1) THEN
        IF(LVPRS(IV).EQ.2.AND.PRSFQT.EQ.2
     *  .AND.MYR.GE.LPR1ST.AND.MYR.LE.LPRLST) PRSPIN=0.5
      ENDIF
C
      IF(PRGFLG.EQ.2.AND.ICY.EQ.LPGSY+1) THEN
        IF(LVPRG(IV).EQ.2.AND.PRGFQT.EQ.2
     *  .AND.MYR.GE.LPG1ST.AND.MYR.LE.LPGLST) PRGPIN=0.5
      ENDIF
C
C  Determine percent of vehicles equipped with evap systems
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
      IF(TRAT.GT.PPRT) TRAT=PPRT
      IF(PRGR.GT.PPRT) PRGR=PPRT
      IF(PRSR.GT.PPRT) PRSR=PPRT
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
      IF(OVLP(2).LT.0.0) OVLP(2)=0.0
      OVLP(5)=PRSR+PRGR-PPRT-OVLP(2)
      OVLP(4)=TRAT-OVLP(5)
      OVLP(1)=PPRT-PRGR-OVLP(4)
C
      OVLP(3)=PRGR-OVLP(2)-OVLP(5)
C
      IF(OVLP(3).LT.0.0) THEN
         OVLP(3)=0.0
         OVLP(5)=PRGR
         OVLP(4)=TRAT-OVLP(5)
         OVLP(1)=PRSR-OVLP(4)-OVLP(5)-OVLP(2)
      ENDIF
C
C  Initialize rates in overlaps
C
      DO 20 I=1,5
        OVLR(I)=OVLP(I)
   20 CONTINUE
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
   40 EADJ=1.0-ENFORC(CRPRS,1)*PPEFF(JDX,PRSFQT)*PRSPIN
      ADJ=1.0-EADJ
      IF(PRSPGM.EQ.1) DPRES=0.0
C
      OVLR(1)=DPRES*ADJ*OVLP(1)+EADJ*OVLP(1)
      OVLR(2)=DPRES*ADJ*OVLP(2)+EADJ*OVLP(2)
      OVLR(4)=DPRES*ADJ*OVLP(4)+EADJ*OVLP(4)
      OVLR(5)=DPRES*ADJ*OVLP(5)+EADJ*OVLP(5)
C
      OVLR(3)=OVLP(3)+OVLP(2)-OVLR(2)+OVLP(5)-OVLR(5)
      GOTO 70
C
C  Purge Check Alone
C  Adjust rates to reflect programs
C
   50 EADJ=1.0-ENFORC(CRPRG,1)*PPEFF(JDX,PRGFQT)*PRGPIN
      ADJ=1.0-EADJ
      IF(PRGPGM.EQ.1) DPURG=0.0
C
      OVLR(2)=DPURG*ADJ*OVLP(2)+EADJ*OVLP(2)
      OVLR(3)=DPURG*ADJ*OVLP(3)+EADJ*OVLP(3)
      OVLR(5)=DPURG*ADJ*OVLP(5)+EADJ*OVLP(5)
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
   60 EPRES=1.0-ENFORC(CRPRS,1)*PPEFF(JDX,PRSFQT)*PRSPIN
      IF(PRSPGM.EQ.1) DPRES=0.0
      FPRES=EPRES+DPRES-DPRES*EPRES
C
      EPURG=1.0-ENFORC(CRPRG,1)*PPEFF(JDX,PRGFQT)*PRGPIN
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
   70 IF(IK.EQ.1) FAILOC=OVLR(3)
      IF(IK.EQ.2) FAILOC=OVLR(1)+OVLR(2)+OVLR(4)+OVLR(5)
C
   98 FAIL=FAILOC
      RETURN
C
      END
