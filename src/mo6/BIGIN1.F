      SUBROUTINE BIGIN1(IH) 
C 
C  BIGIN1 sets up the BIGCFX coefficient calculations loops by 
C  initializing parameters and variables / arrays used in those 
C  calculations. 
C 
C  Called by BIGCFX. 
C 
C  Calls 
C 
C  Changes : (Last change first) 
C 
C  14 Jun 00 AIR Task 03: Removed initialization of FCC, FCN, FHC and FHN 
C                         to BIGIN2. 
C  02 Jun 00 AIR Task 03: Modified MOBILE5B code to MOBILE6 useage. 
C  23 May 00 @EPA-elg   Removed TEMP_ST array. 
C  01 May 00 @EPA- BG    Removed use of DATA statement initializations       
C  14 Oct 99 @DynTel-ddj 2-714  Removed array OFFCO.  OFFCO is now a 
C            local variable in MAINTCF.  Removed array DIST_SUM. 
C            Variables in array were always equal to 1. 
C  26 May 99 @Dyntel-HXQ 2-698 delete code of including scene1.i 
C  12 MAY 99 @DynTel-BAG 2-650 Replaced call to MINUTES with use of 
C            precalculated soak times, SOK_TIME; SOAKDIST(IH,ST_INDX,IWEEK) 
C            now has hourly index IH increased to 24 (previously was 14) 
C  12 May 99 @DynTel-HQ 2-698 Delete the code to call SCFSUB 
C  15 Apr 99 @DynTel-HQ 2-698 operating mode data are from block data 
C  31 Mar 99 @DynTel-ddj #2-696  Initialized summation of hourly 
C              soak time distribution, array DIST_SUM 
C  14 Aug 98 @DynTel-zk 2-000  Explicit typing 
C 
C 
C  Input on call: 
C 
C    argument list: IH 
C    Constants: 
C    MAXVEH from include file IVTYPE.I 
C 
C    common blocks: 
C    /FLAGS4/ COLDFG 
C    /MAXIMA/ MAXVEH,MAXPOL,MAXYRS 
C    /TEMPS/  TEMHRLY 
C    /TEMPC5/ NCC,CCSTD 
C    /VVNAMS/ VVLDGAS   (IVTYPE.I) 
C    /WEEKLY/ IWEEK 
C 
C  Output on return: 
C 
C    common blocks: 
C    /OMTCOM/ TCFSTART,TCFRUN,TCF12 
C    /OFFSET/ OFFSET 
C    /SPEED3/ SADJ 
C    /TEMPC3/ MDLOHI,TDIFF 
C    /TEMPC6/ CSUB,CMUL 
C 
C  Local array subscripts: 
C 
C 
C  Local variable dictionary: 
C 
C   Name   Type              Description 
C  ------  ----  ------------------------------------------------------- 
C  AVESPD   R    temporary average speed to allow SPD to be passed into  
C                SCFSUB 
C  FTPT     R    FTP standard temperature 
C  ICOEF5   I    pointer to term in 5th degree polynomial scf equation 
C  ISADJ    I    integer location used in rounding the scenario's SADJ 
C                for 1 used & then restored for a later use) 
C 
C  SCF      R    Low/Mid speed correction factor for the given IV, IP and IG 
C  PASSRA   R    ratio of passing vehicle emission to the standard - assumed 
C                independent of temperature for cold CO standard modeling 
C 
C 
C  Notes: 
C 
C  Although calculated for all 8 IV, SCUNA1 results are only used by BIGCFX for 
C  subgroup IVA (LDGV, LDGT(1,2,3,4), MC). 
C  BIGIN1 was modified for MOBILE4.1. 
C  SCF1, COMMON SPEED1 and SPEED4 was removed, SCUNA1 is calculated 
C  in SCFSUB; July-14-1993 @ ARC-bsg Subtask 277 
C 
C  11 Jun 97 DynTel@MLA 2-622  Added IVTYPE.I for Mobile6 and modified the 
C            code to use the symbolic vehicle type names and vehicle type 
C            vectors that were introduced in Mobile6. 
C  12 November 1997 Dyntel@MMS 2-645 Removed index of TEMEXH(IP) 
C  18 November 1997 Dyntel@MMS 2-645 Removed IP index of TDIFF and MDKOHI 
C  16 December 1997 Dyntel@MMS 2-645 Removed multiple spd by modifying SPDFLG 
C                    and replacing it with RNLFLG, and SPD         
C  23 Jul 98 @DynTel-MLA 2-000  Changed LDGT names to conform to 
C            certification naming conventions. 
C 
      IMPLICIT NONE 
      INCLUDE 'FLAGS4.I' 
      INCLUDE 'IVTYPE.I' 
      INCLUDE 'MAXIMA.I' 
      INCLUDE 'OFFSET.I' 
      INCLUDE 'OMTCOM.I' 
      INCLUDE 'SCENE2.I' 
      INCLUDE 'SPEED3.I' 
      INCLUDE 'TEMPC3.I' 
      INCLUDE 'TEMPC5.I' 
      INCLUDE 'TEMPC6.I' 
      INCLUDE 'TEMPS.I' 
C 
C     Declare parameter list 
C 
      INTEGER, INTENT(IN) :: IH 
C 
C     Declare local variables 
C 
      INTEGER IVLDG 
      INTEGER ICOLD 
      REAL    TRISE 
      REAL    PASSRA 
C  Compute operating mode fractions. 
C 
C  In MOBILE5b, the exhaust temperature was calculated from the following equation: 
C  TEMEXH=2.2857+0.97674*TEMMIN+0.56881*TRISE+0.0024642*TEMMIN*TRISE. In MOBILE6 
C  set TEMEXH to the hourly ambient temperature. 
C 
      TEMEXH=TEMHRLY(IH) 
C 
C  From MOBILE5b: 
C  Calculate the hot soak temperature if TEMMAX is 
C  at or above 40F, otherwise use the daily average temperature. 
C 
      TEMEVP=(TEMMAX+TEMMIN)/2. 
C 
      TRISE=TEMMAX-TEMMIN 
C 
      IF(TEMMAX.GE.40..AND.TRISE.GT.0.) 
     *   TEMEVP=-1.7474+1.0290*TEMMIN+0.99202*TRISE 
     *          -0.0025173*TEMMIN*TRISE 
C 
      TEMEVP=FLOAT(IFIX(TEMEVP*10.+.5))/10. 
C 
C  Setup for temperature correction factor loop: decide TCF case. 
C 
      TDIFF=TEMEXH-75.0 
C 
      IF(TDIFF.EQ.0) MDLOHI=0 
      IF(TDIFF.LT.0) MDLOHI=1 
      IF(TDIFF.GT.0) MDLOHI=2 
C 
C  If cold temperature and cold CO standard case, compute appropriate cfs. 
C 
      IF(MDLOHI.NE.1.OR.COLDFG.EQ.1) GOTO 30 
C 
      PASSRA=2.7/3.4 
      DO ICOLD=1,NCC 
        DO IVLDG=1,5 
C 
          IF(COLDFG.EQ.2) CSUB(IVLDG,ICOLD)= 
     *        ((75.0-TEMEXH)/55.0)*(FCC/.206) 
     *        *(13.2-(PASSRA*CCSTD(IVLDG,ICOLD))) 
          IF(CSUB(IVLDG,ICOLD).LT.0.0) CSUB(IVLDG,ICOLD)=0.0 
C 
          IF(COLDFG.EQ.3) CMUL(IVLDG,ICOLD)= 
     *     ((PASSRA*CCSTD(IVLDG,ICOLD))-2.7)/10.5 
          IF(CMUL(IVLDG,ICOLD).GT.1.0) CMUL(IVLDG,ICOLD)=1.0 
C 
        END DO 
      END DO 
C 
C  Set running and start temperature correction factors and the offset 
C  to their neutral (no correction) values. 
C 
   30 TCFRUN   = 1.0 
C 
      TCFSTART = 1.0 
C 
      TCF12    = 1.0 
C 
      OFFSET   = 0.0 
C 
C 
      RETURN 
      END 
