C Version @(#)$Id$ $Source$ $Date$ 

      SUBROUTINE BIGIN1
C
C  BIGIN1 sets up the BIGCFX coefficient calculations loops by initializing
C  parameters and variables / arrays used in those calculations.  For example,
C  BIGIN1 evaluates the 5th degree polynomial speed correction factor
C  equation for all 8 IV, 3 IP and 21 IG.  For HC, CO and 1 case of NOx, the
C  results are exponentiated.
C
C  Called by BIGCFX.
C
C  Input on call:
C
C    common blocks:
C    /FLAGS4/ COLDFG
C    /MAXIMA/ MAXVEH,MAXPOL,MAXYRS
C    /SCENE1/ SPD,PCCN,PCHC,PCCC
C    /SPEED1/ SK1
C    /SPEED4/ MAXSP1
C    /TEMPS/  TEMEXH
C    /TEMPC5/ NCC,CCSTD
C
C  Output on return:
C
C    common blocks:
C    /OFFSET/ OFFCO
C    /SCENE2/ FCC,FCN,FHC,FHN
C    /SPEED3/ SADJ,SCUNA1
C    /SPEED4/ SCUNA0
C    /TEMPC3/ MDLOHI,TDIFF
C    /TEMPC6/ CSUB,CMUL
C
C  Local array subscripts:
C
C  FTPSPD(3)  -  FTPSPD ( IB )
C
C  Local variable dictionary:
C
C   Name   Type              Description
C  ------  ----  -------------------------------------------------------
C  AVESPD   R    average speed for the vehicle type IV
C  FTPSPD   R    FTP standard speeds for the 3 operating modes
C  FTPT     R    FTP standard temperature
C  ICOEF5   I    pointer to term in 5th degree polynomial scf equation
C  ISADJ    I    integer location used in rounding the scenario's SADJ
C                for 1 used & then restored for a later use)
C  PASSRA   R    ratio of passing vehicle emission to the standard - assumed
C                independent of temperature for cold CO standard modeling
C  SCF      R    Low/Mid speed correction factor for the given IV, IP and IG
C  SCF1     R    16mph correction factor for the given IV, IP and IG
C  TMPSPD   R    tempory average speed for the vehicle type IV to keep SCF
C                from being corrected for speeds greater than 48 mph
C
C  Notes:
C
C  Although calculated for all 8 IV, SCUNA1 results are only used by BIGCFX for
C  subgroup IVA (LDGV, LDGT1/2, MC). SCUNA0 is calculated at 16.0.
C  BIGIN1 was modified for MOBILE4.1.
C
C
      COMMON /FLAGS4/ PRTFLG,IDLFLG,NMHFLG,HCFLAG,COLDFG
      INTEGER PRTFLG,HCFLAG,COLDFG
C
      COMMON /MAXIMA/ MAXVEH,MAXLTW,MAXPOL,MAXREG,MAXYRS
      COMMON /OFFSET/ OFFCO(25,3),OFFMTH(25,8,2)
      COMMON /SCENE1/ SPD(8),PCCN,PCHC,PCCC
      COMMON /SCENE2/ FCC,FCN,FHC,FHN
      COMMON /SPEED1/ SK1(6,18,3),SK2C(2,14,3,3),SK4C(2,14,3,3)
      COMMON /SPEED3/ SADJ,SCUNA1(18,3,8),SCADJ1(3,4),ISPGRP(3,4)
      COMMON /SPEED4/ MAXSP1,LB1STS(3),LBLAST,IVB,IVA,SCUNA0(18,3,8)
      COMMON /TEMPC3/ MDLOHI(3),TDIFF(3),CHFRAC(3),TCF(3)
      COMMON /TEMPC5/ NCC,MYCOL(2),CCSTD(3,2)
      COMMON /TEMPC6/ CSUB(3,2),CMUL(3,2),PIVPCT(3,3,2)
      COMMON /TEMPS/ AMBT,TEMMIN,TEMMAX,TEMEXH(3),TEMEVP(6)
C
      DIMENSION FTPSPD(3)
C
      DATA FTPSPD/26.0,16.0,26.0/,FTPT/75.0/
C
C  Compute operating mode fractions.
C
      FCN=PCCN*.01
      FHC=PCHC*.01
      FCC=PCCC*.01
      FHN=(FCC-FCN)+FHC
C
C  Setup for temperature correction factor loop: decide TCF case.
C
      DO 10 IP=1,MAXPOL
      TDIFF(IP)=TEMEXH(IP)-FTPT
      IF(TDIFF(IP).EQ.0) MDLOHI(IP)=0
      IF(TDIFF(IP).LT.0) MDLOHI(IP)=1
      IF(TDIFF(IP).GT.0) MDLOHI(IP)=2
   10 CONTINUE
C
C  If cold temperature and cold CO standard case, compute appropriate cfs.
C
      IF(MDLOHI(2).NE.1.OR.COLDFG.EQ.1) GOTO 30
C
      PASSRA=2.7/3.4
      DO 20 ICOLD=1,NCC
      DO 20 IVLDG=1,3
C
      IF(COLDFG.EQ.2) CSUB(IVLDG,ICOLD)=
     *    ((75.0-TEMEXH(2))/55.0)*(FCC/.206)
     *    *(13.2-(PASSRA*CCSTD(IVLDG,ICOLD)))
      IF(CSUB(IVLDG,ICOLD).LT.0.0) CSUB(IVLDG,ICOLD)=0.0
C
      IF(COLDFG.EQ.3) CMUL(IVLDG,ICOLD)=
     * ((PASSRA*CCSTD(IVLDG,ICOLD))-2.7)/10.5
      IF(CMUL(IVLDG,ICOLD).GT.1.0) CMUL(IVLDG,ICOLD)=1.0
C
   20 CONTINUE
C
C  Zero out OFFCO.
C
   30 DO 40 IV=1,3
      DO 40 IDX=1,MAXYRS
      OFFCO(IDX,IV)=0.0
   40 CONTINUE
C
C  Calculate average speed for the operating mode fractions.
C
      SADJ=PCCC/FTPSPD(1)
     -    +(100.-PCCC-PCHC)/FTPSPD(2)
     -    +PCHC/FTPSPD(3)
      SADJ=100./SADJ
      ISADJ=SADJ*10.+0.5
      SADJ=ISADJ/10.
C
C  Build the 5th degree polynomial ("old") equation unadjusted scf table.
C  Fill in SCUNA0 Array, it's the numerator in calculating SCIADJ;
C
      DO 60 IV=1,MAXVEH
      AVESPD=SPD(IV)
C
      DO 60 IP=1,MAXPOL
C
      DO 60 ISP=1,MAXSP1
      SCF=SK1(6,ISP,IP)
      SCF1=SCF
C
      TMPSPD=AVESPD
      IF((IV.LE.3.OR.IV.EQ.8).AND.AVESPD.GT.48.0) TMPSPD=48.0
      DO 55 ICOEF5=1,5
      SCF=SCF*TMPSPD+SK1(6-ICOEF5,ISP,IP)
      SCF1=SCF1*16.0+SK1(6-ICOEF5,ISP,IP)
   55 CONTINUE
C
      IF(IP.NE.3) THEN
       SCF=EXP(SCF)
       SCF1=EXP(SCF1)
      ENDIF
C
      SCUNA1(ISP,IP,IV)=SCF
      SCUNA0(ISP,IP,IV)=SCF1
C
  60  CONTINUE
C
C  Initialize ISPGRP.
C
      DO 70 IV=1,4
      DO 70 IP=1,3
      ISPGRP(IP,IV)=0
   70 CONTINUE
C
      RETURN
      END
