      SUBROUTINE RVPEXH(ICY) 
C 
C  MOBILE4 calculates a correction factor array to adjust the base emission 
C  factor for fuel volatility (= Reid Vapor Pressure = RVP).  This 
C  multiplicative cf is neutral (= 1.0) only when (1) the actual RVP <= 9.0 
C  (= the RVP of the fuel used in the tests upon which the ber slopes and 
C  intercepts in Block Data are based) and / or (2) the ambient temperature 
C  (TEMEXH) is <= 45.0 (below 45 degrees there is no volatility effect) & / or 
C  the vehicle class is not equipped in the given MY with the technology known 
C  by testing to be subject to an RVP exhaust effect on the given pollutant. 
C  Thus the effect depends on IP and MY as well as RVP and TEMEXH. 
C 
C  There are 3 nonneutral cases.  RVPEXH computes the MYG 1 HC & CO case, 
C  as well as initializing the cf array to 1.0 to cover the neutral cases. 
C  BIGTCF computes the 2 operating mode ("bag") dependent cases: MYGs 2 & 3 
C  >= 75 degrees and MYGs 2 & 3 < 75 degrees. 
C 
C  Called by EFCALX. 
C 
C  Calls AMIN1 and IGRVPT. 
C 
C  Changes: (Last change first) 
C 
C  26 Feb 02 AIR Bug 17: set entire array of RVPCF to 1 to avoid problems with 
C            air toxic's use of FUELCF. 
C  12 Jun 00 @EPA-djb Bug99 Replaced RVPX with RVPICY. 
C  08 Jun 00 AIR Task 03: Removed nonblock DO constructs. 
C  08 Jun 00: Removed type declarations for intrinsic functions. 
C  10 Apr 00 AIR Task 02: Removed calls to OUTPOL by including DBSELPOL. 
C  22 Feb 99 @DynTel-DDJ 2-696 Added hourly indices to RVPCF initialization 
C            and hourly temperatures 
C  06 Jan 99 @DynTel-ZK 2-000 Explicit Typing 
C 
C  Input on call: 
C 
C    Constants:  
C    MAXVEH from include file IVTYPE.I 
C 
C    argument list: ICY 
C 
C    common blocks: 
C    /CITRV1/ RVPICY 
C    /CITRV2/ RVP090 
C    /MAXIMA/ MAXPOL,MAXYRS 
C    /RVPEX1/ R9G1,TLL 
C    /TEMPS/  TEMEXH 
C    /VVNAMS/ VVGAS, VVLDGAS  (IVTYPE.I) 
C 
C  Output on return: 
C 
C    common blocks: 
C    /RVPEX2/ RVPCF 
C 
C  Local variable dictionary: 
C 
C   Name   Type              Description 
C  ------  ----  ------------------------------------------------------- 
C  CFB      R    RVP cf at begin point of temperature interpolation interval 
C  CFL      R    RVP cf at end   point of temperature interpolation interval 
C  RADJCF   R    RVP cf assessed at FTP (75.0) temperature - accounts only for 
C                RVP deviation from norm, not temperature. 
C  REAT     R    RVP exhaust IP temperature 
C 
C  Notes: 
C 
C  Support for oxygenated and cert. fuels from MOBILE4 was removed in 4.1. 
C 
C  11 Jun 97 DynTel@MLA 2-622  Added IVTYPE.I for Mobile6 and modified the 
C            code to use the symbolic vehicle type names and vehicle type 
C            vectors that were introduced in Mobile6. 
C  12 November 1997 Dyntel@MMS 2-645 removed index of TEMEXH(IP) 
C  23 Jul 98 @DynTel-MLA 2-000  Changed LDGT names to conform to 
C            certification naming conventions. 
C 
      USE DATABASE, ONLY: DBSELPOL 
      IMPLICIT NONE 
C 
      INCLUDE 'CITRV1.I' 
      INCLUDE 'CITRV2.I' 
      INCLUDE 'IVTYPE.I' 
      INCLUDE 'MAXIMA.I' 
      INCLUDE 'RVPEX1.I' 
      INCLUDE 'RVPEX2.I' 
      INCLUDE 'TEMPS.I' 
C 
      INTEGER IV,IVTAM,IP,IDX,ICY,IG,IGRVPT,MY,IH 
      REAL RVP1,REAT,RIVAL,RADJCF,CFB,CFL 
C 
C  Initialize the RVP exhaust correction factor to the neutral 1.0 and depart, 
C  if no IVTAM x IP x IDX cases will be in the "not-by-bag" computed category. 
C 
      IVTAM = 0 
      RVPCF=1.0 
C 
C  The RVPCF algorithm applies only to HC & CO (not NOx) and MYG 1 (1970-79 
C  for LDGV, 1970-80 for LDGT12, 1978-80 for LDGT34, no MYG 1 for HDGV). 
C 
      IF((DBSELPOL(1).EQ.1 .AND. 
     *    DBSELPOL(2).EQ.1 .AND. 
     *    DBSELPOL(3).EQ.2).OR.ICY-MAXYRS.GT.1980) RETURN 
C 
C  Reduce supplied RVP to 11.7 for the RVP cf calculation only. 
C 
      RVP1=AMIN1(RVPICY,11.7) 
C 
C  Now can identify and compute the "not-by-bag" cases for LDGV/T. 
C 
      DO IH = 1, MAXIH 
        IVTAM = 0 
        DO 40 IV=1,MAXVEH 
C 
          IF (VVLDGAS(IV).EQ.0) CYCLE 
          IVTAM = IVTAM + 1 
C 
        DO 30 IP=1,MAXPOL 
         IF(DBSELPOL(IP).EQ.1) GO TO 30 
 
          IF(IP.EQ.3) GO TO 30 
C 
C  If IP's temperature is <= 45 degrees, then skip this IP iteration. 
C  If it is above 95 degrees, then reduce it to 95 for this computation. 
C 
 
          IF(TEMHRLY(IH).LE.45.0) GOTO 30 
          REAT=AMIN1(TEMHRLY(IH),95.0) 
C 
            DO 20 IDX=1,MAXYRS 
              MY=ICY+IDX-MAXYRS 
              IG=IGRVPT(MY,IVTAM) 
C 
C  If MYG 0, then no computation is needed for this MY.  Do next IDX iteration. 
C  If MYG 2 has been reached, further calculations for this IVTAM x IP case 
C  will be done in BIGTCF, so exit IDX loop. 
C 
              IF(IG.EQ.0) GOTO 20 
              IF(IG.GT.1) GOTO 30 
C 
C  Use base or in-use RVP.  If the selected RVP is <= 9.0, skip the 
C  computation.  Once RIVAL <= 9.0, exit IDX loop, not just 1 iteration. 
C 
              RIVAL=RVP1 
              IF(RIVAL.LE.RVP090) GOTO 30 
C 
C  Compute the cf.  There are 2 cases, regarding temperature. 
C 
C  (1) Temperature is neutral. 
C 
              RADJCF=(R9G1(1,IP)+R9G1(2,IP)*RIVAL)/R9G1(3,IP) 
              IF(REAT.EQ.75.0) RVPCF(IDX,IP,IVTAM,IH)=RADJCF 
              IF(REAT.EQ.75.0) GOTO 20 
C 
C  (2) Temperature is a variable in the RVP cf calculation.  Interpolate 
C      between 2 cf's computed at temperatures for which data exists. 
C 
C      (a) First assign the endpoints at 45 and 75 degrees. 
C 
              CFB=1.0 
              CFL=RADJCF 
C 
C      (b) Now interpolate between the 2 endpoints (CFB,CFL) or extrapolate 
C          beyond CFL (up to 95), using actual (REAT) and (TLL) temperatures. 
C 
              RVPCF(IDX,IP,IVTAM,IH)= 
     *          CFB+(CFL-CFB)*(REAT-TLL(1))/TLL(2) 
C 
   20       CONTINUE 
   30     CONTINUE 
   40   CONTINUE 
      END DO 
C 
      RETURN 
      END 
