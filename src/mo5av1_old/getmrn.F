C Version @(#)$Id$ $Source$ $Date$ 

      SUBROUTINE GETMRN(INERR,*)
C
C  GETMRN checks for and reads in and screens user supplied mileage accrual
C  rates and registration distributions.  Whether or not the distributions are
C  replaced, they are normalized.
C
C  Called by ONESEC.
C
C  Calls QUITER.
C
C  Input on call:
C
C    parameter list: INERR
C    common blocks:
C    /FLAGS1/ PROMPT
C    /FLAGS2/ MYMRFG
C    /IOUCOM/ IOUGEN,IOUASK
C    /MAXIMA/ MAXVEH,MAXYRS
C    /NAMES1/ NAMMMR
C
C  Output on return:
C
C    parameter list: INERR
C    common blocks:
C    /MYRSAV/ AMAR,JULMYR
C
C  Local variable / array dictionary:
C
C   Name   Type              Description
C  ------  ----  -------------------------------------------------------
C  CHKMYR   R    myr mix check: sum of the myr fractions = 1.0 +/- 1.E-6
C                (a normalized myr distribution has been entered for each IV)
C  MYF      I    my of first (JDX order) AMAR / JULMYR to be entered on READ
C  MYL      I    my of last (JDX order) AMAR / JULMYR to be entered on READ
C  SUMMYR   R    sum of myr for IV; used to normalize user entered JULMYR
C
C  Notes:
C
C  None.
C
C
      CHARACTER*8 NAMFLG,NAMMMR
      CHARACTER*8 NAMRVP
      CHARACTER*4 NAMTEM
      INTEGER PROMPT,TAMFLG,SPDFLG,VMFLAG,OXYFLG,DSFLAG
      INTEGER ALHFLG
      REAL JULMYR
C
      COMMON /FLAGS1/ PROMPT,TAMFLG,SPDFLG,VMFLAG,OXYFLG,DSFLAG
      COMMON /FLAGS2/ MYMRFG,NEWFLG,IMFLAG,ALHFLG,IMTIER
      COMMON /IOUCOM/ IOUIMD,IOUGEN,IOUREP,IOUERR,IOUASK
      COMMON /MAXIMA/ MAXVEH,MAXLTW,MAXPOL,MAXREG,MAXYRS
      COMMON /MYRSAV/ AMAR(25,8),JULMYR(25,8),NEWCUM
      COMMON /NAMES1/ NAMFLG(20),NAMMMR(2,3,8),NAMRVP(4),NAMTEM(2)
C
C  Check for user supplied my mileage accrual rates and / or registration
C  distributions.
C
      IF(MYMRFG.EQ.1) GOTO 65
      IF(MYMRFG.EQ.3) GOTO 20
C
C  Get mileage accrual rates.
C
      IF(PROMPT.EQ.2) WRITE(IOUASK,200)
  200 FORMAT(/)
C
      DO 10 IV=1,MAXVEH
      DO 10 IG=1,3
      IF(PROMPT.EQ.2) WRITE(IOUASK,210) (NAMMMR(ICH,IG,IV),ICH=1,2)
  210 FORMAT(' ','Enter MYM ages ',2A8,':')
      MYF=1
      IF(IG.EQ.2) MYF=11
      IF(IG.EQ.3) MYF=21
      MYL=10
      IF(IG.EQ.2) MYL=20
      IF(IG.EQ.3) MYL=25
      READ(IOUGEN,110,ERR=97,END=99) (AMAR(JDX,IV),JDX=MYF,MYL)
  110 FORMAT(10F7.5)
   10 CONTINUE
C
      DO 15 IV=1,MAXVEH
      DO 15 JDX=1,MAXYRS
      IF(AMAR(JDX,IV).LT.0.0) CALL QUITER(AMAR(JDX,IV),0,19,INERR)
   15 CONTINUE
C
      IF(MYMRFG.EQ.2) GOTO 45
C
C  Get vehicle registration distributions.
C
   20 IF(PROMPT.EQ.2) WRITE(IOUASK,200)
      DO 35 IV=1,MAXVEH
C
      DO 25 IG=1,3
      IF(PROMPT.EQ.2) WRITE(IOUASK,220) (NAMMMR(ICH,IG,IV),ICH=1,2)
  220 FORMAT(' ','Enter MYR ages ',2A8,':')
      MYF=1
      IF(IG.EQ.2) MYF=11
      IF(IG.EQ.3) MYF=21
      MYL=10
      IF(IG.EQ.2) MYL=20
      IF(IG.EQ.3) MYL=25
      READ(IOUGEN,120,ERR=98,END=99) (JULMYR(JDX,IV),JDX=MYF,MYL)
  120 FORMAT(10F5.3)
   25 CONTINUE
C
      CHKMYR=0.0
      DO 30 JDX=1,MAXYRS
      CHKMYR=JULMYR(JDX,IV)+CHKMYR
   30 CONTINUE
C
      IF(CHKMYR-.00001.GT.1..OR.CHKMYR+.00001.LT.1.)
     *  CALL QUITER(CHKMYR,0,49,INERR)
   35 CONTINUE
C
      DO 40 IV=1,MAXVEH
      DO 40 JDX=1,MAXYRS
      IF(JULMYR(JDX,IV).LT.0.0) CALL QUITER(JULMYR(JDX,IV),0,20,INERR)
   40 CONTINUE
C
   45 DO 55 IV=1,MAXVEH
C
C  Skip 0 with non-0 checks if user supplies MYR only and on HDDV pass:
C  the BD values are 0.0's for all 20 my's (HDDMYM has not yet been called).
C
      IF(MYMRFG.EQ.3.AND.IV.EQ.7) GOTO 55
      DO 50 JDX=1,MAXYRS
      IF(AMAR(JDX,IV).LE.0.0.AND.JULMYR(JDX,IV).GT.0.0)
     *    CALL QUITER(JULMYR(JDX,IV),0,21,INERR)
      IF(AMAR(JDX,IV).GT.0.0.AND.JULMYR(JDX,IV).EQ.0.0)
     *    CALL QUITER(AMAR(JDX,IV),0,22,INERR)
   50 CONTINUE
   55 CONTINUE
C
C  Check vehicle equalities (LDGV,LDDV), (LDGT1,LDDT) mod yr regs.
C
      DO 60 JDX=1,MAXYRS
      IF(JULMYR(JDX,1).NE.JULMYR(JDX,5)) CALL QUITER(0.,JDX,65,INERR)
      IF(JULMYR(JDX,2).NE.JULMYR(JDX,6)) CALL QUITER(0.,JDX,66,INERR)
   60 CONTINUE
C
C  Normalize JULMYR.
C
   65 DO 75 IV=1,MAXVEH
C
      SUMMYR=0.0
      DO 70 JDX=1,MAXYRS
      SUMMYR=JULMYR(JDX,IV)+SUMMYR
   70 CONTINUE
      IF(SUMMYR.EQ.0.0) SUMMYR=1.
C
      DO 75 JDX=1,MAXYRS
      JULMYR(JDX,IV)=JULMYR(JDX,IV)/SUMMYR
   75 CONTINUE
C
      RETURN
C
C  Error or EOF on any attempted read => take alternate
C  return 1 => run aborts.
C
   97 CALL QUITER(0.,0,122,INERR)
      GOTO 99
C
   98 CALL QUITER(0.,0,123,INERR)
   99 RETURN1
      END
