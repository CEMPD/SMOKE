C Version @(#)$Id$ $Source$ $Date$ 

      FUNCTION RSTLOS(MY,IV,TEMUSE)
C
C  RSTLOS computes the Resting Loss rate (g/hr) for each model year
C  and gas vehicle class at the given (hot soak) temperature
C
C  Called by HCCALX, HOTSOK, DIRUNL, and RULOSS.
C
C  Input on call:
C
C    parameter list: MY,IV,TEMUSE
C    common blocks:
C    /FLAGS3/ TEMFLG
C    /RESTLS/ PTGRST,RESTA,RESTB
C    /TEMPS/  TEMEVP
C
C  Output on return:
C
C    function: RSTLOS
C
C  Local variable / array dictionary:
C
C   Name   Type              Description
C  ------  ----  -------------------------------------------------------
C  MYRSTL   I    model year group number for PTGRST
C  RESTT    R    resting loss temperature (TEMEVP(1) or AMBT)
C                temperature and lower boundary RVP
C  RESTL    R    resting loss emission rate to be returned (g/hr)
C
C  Notes:
C
C  RSTLOS was added for MOBILE 4.1, version D03.
C
      INTEGER ATPFLG,TPDFLG,RLFLAG,TEMFLG,OUTFMT
C
      COMMON /FLAGS3/ ATPFLG,TPDFLG,RLFLAG,LOCFLG,TEMFLG,OUTFMT
      COMMON /RESTLS/ PTGRST(13,5,2),RESTA(2),RESTB
      COMMON /TEMPS/ AMBT,TEMMIN,TEMMAX,TEMEXH(3),TEMEVP(6)
C
C  Initialize
C
      RSTLOS=0.0
C
C  Determine vehicle index.  Currently, only LDGVs,LDGTs,and MCs have
C  resting losses, and LDGT1=LDGT2.  Also, MCs=open bottom only
C  Determine MYRSTL from MY.
C
      IF(IV.GE.5.AND.IV.LE.7) RETURN
      IVGAS=IV
      IF(IVGAS.EQ.8) IVGAS=5
      MYRSTL=MY-1979
      IF(MYRSTL.LT.1) MYRSTL=1
      IF(MYRSTL.GT.13) MYRSTL=13
C
C  Resting loss is zero if any evaporative temperature is 40F or less,
C  or if minimum temperature is 25F or less.
C
      IF(TEMMAX   .LE.40.0 .OR.
     *   AMBT     .LE.40.0.AND.TEMFLG.EQ.2 .OR.
     *   TEMEVP(1).LE.40.0 .OR.
     *   TEMEVP(2).LE.40.0 .OR.
     *   TEMEVP(6).LE.40.0 .OR.
     *   TEMMIN   .LE.25.0) RETURN
C
C  Choose Running Loss Temperature and reset if above limit (105 F)
C
      RESTT=TEMUSE
      IF(TEMFLG.EQ.2) RESTT=AMBT
      IF(RESTT.GT.105.0) RESTT=105.0
C
C  Loop through the three technologies and add up the resting loss
C
      RESTL=0
      DO 10 ITRSTL=1,2
      RESTL=RESTL+EXP(RESTA(ITRSTL)+RESTB*RESTT)
     *    *PTGRST(MYRSTL,IVGAS,ITRSTL)
   10 CONTINUE
C
      RSTLOS=RESTL
C
      RETURN
      END
