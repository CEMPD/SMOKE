C Version @(#)$Id$ $Source$ $Date$ 

      SUBROUTINE RDCITY(INERR,*)
C
C  RDCITY reads in the scenario's city specific parameters record.  The fields
C  are: city name, ASTM class, minimum/maximum fuel temperatures, two types of
C  fuel RVP levels, in-use RVP controls start year, oxy fuel parameters:
C  oxy flag, alcohol blend and ether blend market share, alcohol blend
C  and ether blend oxygen content, alcohol blend RVP waiver indicator.
C  diesel fuel flag, reformulated fuel flag and effect of industry average
C  fuel.  RDCITY assigns the RVP limits based on the max and min temperatures.
C
C  Called by ONESEC and LOCAL.
C
C  Calls QUITER.
C
C  Input on call:
C
C    common blocks:
C    /FLAGS1/ PROMPT
C    /FLAGS3/ LOCFLG
C    /IOUCOM/ IOUGEN
C    /RVPNAT/ RVPLIM
C
C  Output on return:
C
C    common block:
C    /CITPAR/ SCNAME
C    /CITRV1/ RVPBAS,RVPIUS,IUSESY
C    /FLAGS1/ OXYFLG,DSFLAG
C    /OXY1/   SHRMKT,OXYCNT,IGASHW
C    /REGISF/ USRGSF
C    /RFORM1/ RFGFLG,IAFCNT,WINFLG,ASTMCL
C    /RVPNAT/ RVPDIU
C    /TEMPS/  TEMMIN,TEMMAX
C
C  Local variable / array dictionary:
C
C  Name    Type              Description
C  ------  ----  ----------------------------------------------------
C  DSFCHR   C    input variable for DSFLAG
C  IAFLG    I    variable used to calculate IAFCNT
C  LOWUYR   I    lowest allowed in-use RVP controls start year
C  NUMBER   C    used to parse inputs
C  OXYCHR   C    input variable for OXYFLG
C  RFGCHR   C    input variable for RFGFLG
C  IAFCHR   C    input variable for IAFCNT
C  WINCHR   C    input variable for WINFLG
C
C  Notes:
C
C  RDCITY was changed in MOBILE4.1 to not read ASTM Class
C  Oxygenated fuel support was added for MOBILE4.1.
C  Support for user-entered gas/diesel sales fractions added for 4.1.
C  Reformulated fuel was added for MOBILE5.
C  Industry average fuel effect was added for MOBILE5.
C  ASTM class returned for MOBILE5.
C
      CHARACTER*4 SCNAME
      CHARACTER*1 DSFCHR,OXYCHR,RFGCHR,IAFCHR(3),WINCHR
      CHARACTER*1 ASTMCL
      INTEGER PROMPT,TAMFLG,SPDFLG,VMFLAG,OXYFLG,DSFLAG
      INTEGER RFGFLG,RFGSYR,WINFLG
      INTEGER ATPFLG,TPDFLG,RLFLAG,TEMFLG,OUTFMT
      INTEGER IAFLG(3)
      LOGICAL RFGON
C
      COMMON /CITPAR/ SCNAME(4)
      COMMON /CITRV1/ RVPBAS,RVPIUS,RVPICY,IUSESY,RVPUWX
      COMMON /FLAGS1/ PROMPT,TAMFLG,SPDFLG,VMFLAG,OXYFLG,DSFLAG
      COMMON /FLAGS3/ ATPFLG,TPDFLG,RLFLAG,LOCFLG,TEMFLG,OUTFMT
      COMMON /IOUCOM/ IOUIMD,IOUGEN,IOUREP,IOUERR,IOUASK
      COMMON /MAXIMA/ MAXVEH,MAXLTW,MAXPOL,MAXREG,MAXYRS
      COMMON /OXY1/ SHRMKT(3),OXYCNT(2),IGASHW,OXYFAC(3,2)
      COMMON /REGISF/ GSFVCT(8),USRGSF(25,2)
      COMMON /RFORM1/ IASTM,IPHASE,RFGFLG,IAFCNT,WINFLG,RFGON,ASTMCL
      COMMON /RFORM2/ RFGSYR(3),RFGRVP(5,3)
      COMMON /RVPNAT/ RVPLIM(2),RVPDIU(2),PFUL
      COMMON /TEMPS/ AMBT,TEMMIN,TEMMAX,TEMEXH(3),TEMEVP(6)
      COMMON /YEARS4/ IY1941,IY1960,IY1994,IY2020
C
      CHARACTER*1 NUMBER(10)
      DATA NUMBER/'0','1','2','3','4','5','6','7','8','9'/
      DATA LOWUYR/1970/
C
C  Acquire inputs.
C
      IF(PROMPT.EQ.2) WRITE(IOUASK,200)
  200 FORMAT(' ','Enter scenario name, ASTM class, min & max daily',
     *    ' temps,'/
     *    ' ','period 1 & 2 RVP, period 2 start year, oxy flag,'/
     *    ' ','diesel flag & reform flag:')
CC   *    ' ','diesel flag, reform flag & base fuel flag: ')
      READ(IOUGEN,100,ERR=96,END=96) (SCNAME(I),I=1,4),ASTMCL,TEMMIN,
     *    TEMMAX,RVPBAS,RVPIUS,IUSESY,OXYCHR,DSFCHR,RFGCHR,IAFCHR,WINCHR
  100 FORMAT(4A4,1X,A1,2F5.0,2F5.1,1X,I2,3(1X,A1),1X,3A1,1X,A1)
C
C  Parse OXYFLG, DSFLAG, RFGFLG, WINFLG and IAFCNT columns to prevent input
C  errors.
C
      OXYFLG=1
      DSFLAG=1
      RFGFLG=1
      WINFLG=0
      DO 60 I=1,3
        IAFLG(I)=9
   60 CONTINUE
      DO 70 I=1,10
      IF (OXYCHR.EQ.NUMBER(I)) OXYFLG=I-1
      IF (DSFCHR.EQ.NUMBER(I)) DSFLAG=I-1
      IF (RFGCHR.EQ.NUMBER(I)) RFGFLG=I-1
      IF (WINCHR.EQ.NUMBER(I)) WINFLG=I-1
      IF (IAFCHR(1).EQ.NUMBER(I)) IAFLG(1)=I-1
      IF (IAFCHR(2).EQ.NUMBER(I)) IAFLG(2)=I-1
      IF (IAFCHR(3).EQ.NUMBER(I)) IAFLG(3)=I-1
   70 CONTINUE
C
      IF(RFGFLG.EQ.0.OR.RFGFLG.GT.3)
     *  CALL QUITER(20.,RFGFLG,1,INERR)
C
C  Phase 2 reformulated gas starts in 2021
C
      IF(RFGFLG.EQ.3) THEN
        RFGSYR(3) = 2021
        RFGFLG = 2
      ENDIF
C
C  Determine base fuel characteristics.
C
      IAFCNT=339
      IF(IAFLG(1).EQ.2) IAFCNT=50
C
C  Round the minimum and maximum temperatures to the nearest 0.1F.
C
      TEMMIN=FLOAT(IFIX(TEMMIN*10.+.5))/10.
      TEMMAX=FLOAT(IFIX(TEMMAX*10.+.5))/10.
C
C  Check temperature ranges.  These bounds are for catching data entry errors
C  only.  What constitutes reasonable min & max temperature values depends on
C  local conditions, especially AMBT.
C
      IF(TEMMIN.LT.0.0.OR.TEMMIN.GT.100.0)
     *    CALL QUITER(TEMMIN,1,91,INERR)
C
      IF(TEMMAX.LT.0.0.OR.TEMMAX.GT.120.0)
     *    CALL QUITER(TEMMAX,2,91,INERR)
C
      IF(TEMMAX.LT.TEMMIN) CALL QUITER(TEMMAX,0,92,INERR)
C
C  Compute diurnal evap RVP limits.  Upper bound is dependent on max and min
C  temperatures.
C
      IF(TEMMIN.LE.95.) RVPDIU(2)=26.823-0.004781*TEMMIN-0.13407*TEMMAX
      IF(TEMMIN.GT.95.) RVPDIU(2)=26.000-0.004781*TEMMIN-0.13407*TEMMAX
      IF(RVPDIU(2).GT.RVPLIM(2)) RVPDIU(2)=RVPLIM(2)
C
C  Check RVPs and re-assign if necessary.
C
      IF(RVPBAS.LT.RVPLIM(1)) CALL QUITER(RVPLIM(1),1,90,INERR)
      IF(RVPBAS.LT.RVPLIM(1)) RVPBAS=RVPLIM(1)
C
      IF(RVPBAS.GT.RVPLIM(2)) CALL QUITER(RVPLIM(2),1,90,INERR)
      IF(RVPBAS.GT.RVPLIM(2)) RVPBAS=RVPLIM(2)
C
      IF(RVPIUS.LT.RVPLIM(1)) CALL QUITER(RVPLIM(1),2,90,INERR)
      IF(RVPIUS.LT.RVPLIM(1)) RVPIUS=RVPLIM(1)
C
      IF(RVPIUS.GT.RVPLIM(2)) CALL QUITER(RVPLIM(2),2,90,INERR)
      IF(RVPIUS.GT.RVPLIM(2)) RVPIUS=RVPLIM(2)
C
C  Convert in-use start year to 4 digit no. & then check if in range.
C
      CALL YRTEST(IUSESY,87,LOWUYR,IY2020,INERR)
C
C  Set defaults, then obtain and verify the oxy fuel inputs, if necessary.
C
      IF(OXYFLG.EQ.0.OR.OXYFLG.GT.2)
     *    CALL QUITER(18.,OXYFLG,1,INERR)
C
      SHRMKT(1)=1.0
      SHRMKT(2)=0.0
      SHRMKT(3)=0.0
C
      OXYCNT(1)=0.0
      OXYCNT(2)=0.0
C
      IGASHW=1
C
      IF(OXYFLG.NE.2) GOTO 10
C
      IF(PROMPT.EQ.2) WRITE(IOUASK,210)
  210 FORMAT(' ','Enter ether blend market share, alcohol blend market',
     *    ' share,'/' ','ether blend oxygen content, alcohol blend',
     *    ' oxygen content,'/' ','alcohol blend waiver flag:')
      READ(IOUGEN,110,ERR=97,END=97) SHRMKT(2),SHRMKT(3),OXYCNT(1),
     *                        OXYCNT(2),IGASHW
  110 FORMAT(4(F4.3,1X),I1)
C
      IF(SHRMKT(2).LT.0..OR.SHRMKT(2).GT.1.0)
     *    CALL QUITER(SHRMKT(2),0,2,INERR)
      IF(SHRMKT(3).LT.0..OR.SHRMKT(3).GT.1.0)
     *    CALL QUITER(SHRMKT(3),0,3,INERR)
C
C  Turn SHRMKT(IS) into NDP significant digits
C
      MKTSHR=0
      NDP=3
      DO 80 IS=2,3
         RVAL=SHRMKT(IS)*(10**NDP)
         NVAL=INT(RVAL)
         MKTSHR=MKTSHR+NVAL
         RVAL=1.0*RVAL/(10**NDP)
         SHRMKT(IS)=RVAL
   80 CONTINUE
      IF(MKTSHR.GT.1000) CALL QUITER(SHRMKT(2)+SHRMKT(3),0,4,INERR)
      SHRMKT(1)=1.0-SHRMKT(2)-SHRMKT(3)
      SHRMKT(1)=(NINT(SHRMKT(1)*1000))/1000.
C
      IF(OXYCNT(1).LT.0..OR.OXYCNT(1).GT..027001)
     *    CALL QUITER(OXYCNT(1),0,5,INERR)
      IF(OXYCNT(2).LT.0..OR.OXYCNT(2).GT..035001)
     *    CALL QUITER(OXYCNT(2),0,6,INERR)
C
      IF(IGASHW.LT.1.OR.IGASHW.GT.2) CALL QUITER(0.,IGASHW,7,INERR)
C
      IF(SHRMKT(2).GT.0..AND.OXYCNT(1).EQ.0. .OR.
     *   SHRMKT(2).EQ.0..AND.OXYCNT(1).GT.0.)
     *   CALL QUITER(0.,0,8,INERR)
      IF(SHRMKT(3).GT.0..AND.OXYCNT(2).EQ.0. .OR.
     *   SHRMKT(3).EQ.0..AND.OXYCNT(2).GT.0. .OR.
     *   SHRMKT(3).EQ.0..AND.IGASHW.EQ.2)
     *   CALL QUITER(0.,0,9,INERR)
C
C  Optional user-input of gas/diesel sales fractions (as diesel!)
C
   10 IF(DSFLAG.EQ.0.OR.DSFLAG.GT.2)
     *    CALL QUITER(19.,DSFLAG,1,INERR)
      IF(DSFLAG.NE.2) GOTO 50
      IF(LOCFLG.EQ.2) CALL QUITER(0.,0,118,INERR)
      IF(PROMPT.EQ.2) WRITE(IOUASK,220)
  220 FORMAT(' ','Enter diesel sales fractions (LDV,LDT pairs):')
      DO 20 IG=1,3
      MYF=1
      IF(IG.EQ.2) MYF=11
      IF(IG.EQ.3) MYF=21
      MYL=10
      IF(IG.EQ.2) MYL=20
      IF(IG.EQ.3) MYL=25
      READ(IOUGEN,120,ERR=98,END=98)
     *    ((USRGSF(IDX,IVTEMP),IVTEMP=1,2),IDX=MYF,MYL)
  120 FORMAT(20F4.3)
   20 CONTINUE
      DO 40 IDX=1,MAXYRS
      DO 40 IVTEMP=1,2
      IF(USRGSF(IDX,IVTEMP).GE.0.0.AND.USRGSF(IDX,IVTEMP).LE.1.0)
     *    GOTO 30
      CALL QUITER(USRGSF(IDX,IVTEMP),0,88,INERR)
      USRGSF(IDX,IVTEMP)=0.0
   30 USRGSF(IDX,IVTEMP)=1.0-USRGSF(IDX,IVTEMP)
   40 CONTINUE
C
   50 RETURN
C
C  Error or EOF on any attempted read => take alternate
C  return 1 => run aborts.
C
   96 CALL QUITER(0.,0,131,INERR)
      GOTO 99
   97 CALL QUITER(0.,0,137,INERR)
      GOTO 99
   98 CALL QUITER(0.,0,138,INERR)
      GOTO 99
C
   99 RETURN1
      END
