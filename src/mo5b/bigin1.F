C Version "@(#)$Id$ $Source$ $Date$ 

      SUBROUTINE BIGIN1
C
C  BIGIN1 sets up the BIGCFX coefficient calculations loops by initializing
C  parameters and variables / arrays used in those calculations.
C
C  Called by BIGCFX.
C
C  Calls SCFSUB.
C
C  Input on call:
C
C    common blocks:
C    /FLAGS4/ COLDFG
C    /MAXIMA/ MAXVEH,MAXPOL,MAXYRS
C    /SCENE1/ SPD,PCCN,PCHC,PCCC
C    /TEMPS/  TEMEXH
C    /TEMPC5/ NCC,CCSTD
C
C  Output on return:
C
C    common blocks:
C    /OFFSET/ OFFCO
C    /SCENE2/ FCC,FCN,FHC,FHN
C    /SPEED3/ SADJ
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
C  AVESPD   R    tempory average speed for the vehicle type IV to allow
C                SPD(IV) to be passed into SCFSUB
C  FTPSPD   R    FTP standard speeds for the 3 operating modes
C  FTPT     R    FTP standard temperature
C  ICOEF5   I    pointer to term in 5th degree polynomial scf equation
C  ISADJ    I    integer location used in rounding the scenario's SADJ
C                for 1 used & then restored for a later use)
C  PASSRA   R    ratio of passing vehicle emission to the standard - assumed
C                independent of temperature for cold CO standard modeling
C  SCF      R    Low/Mid speed correction factor for the given IV, IP and IG
C
C
C  Notes:
C
C  Although calculated for all 8 IV, SCUNA1 results are only used by BIGCFX for
C  subgroup IVA (LDGV, LDGT1/2, MC).
C  BIGIN1 was modified for MOBILE4.1.
C  SCF1, COMMON SPEED1 and SPEED4 was removed, SCUNA1 is calculated
C  in SCFSUB; July-14-1993 @ ARC-bsg Subtask 277
C
C
      INCLUDE 'FLAGS4.I'
      INCLUDE 'MAXIMA.I'
      INCLUDE 'OFFSET.I'
      INCLUDE 'SCENE1.I'
      INCLUDE 'SCENE2.I'
      INCLUDE 'SPEED3.I'
      INCLUDE 'TEMPC3.I'
      INCLUDE 'TEMPC5.I'
      INCLUDE 'TEMPC6.I'
      INCLUDE 'TEMPS.I'
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
C
      DO 60 IV=1,MAXVEH
        AVESPD=SPD(IV)
        CALL SCFSUB(AVESPD,IV)
  60  CONTINUE
C
      RETURN
      END
