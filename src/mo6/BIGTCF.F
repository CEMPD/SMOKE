      SUBROUTINE BIGTCF(IP,MYR,IV,IVTMP5) 
C 
C  BIGTCF computes the multiplicative temperature correction factor and, in 
C  some cases, the RVP cf and/or the additive temperature CO offset. 
C 
C  Called by BIGCFX. 
C 
C  Calls ICCOPT, IGRVPT, ITCPTR, IEVPTR, PIV and BIGRVP. 
C 
C  Changes : (Last change first) 
C 
C  27 Jul 01 AIR Task 15-21: Rmoved unused LIDX from parameter list. 
C  03 Jul 01 @EPA-djb Bug329, removed MYR and IV from argument list of BIGRVP 
C  02 Jun 00 AIR Task 03: Updated MOBILE5b code to MOBILE6 specifications. 
C  25 May 00 AIR Task 03: Updated code for Phase 1 CO. 
C  08 May 00 AIR Task 03: Replaced code with MOBILE5B TCF code. 
C  10 Apr 00 AIR Task 02: Substituted DBSELPOL for PRTFLG 
C  13 Jan 00 @DynTel-ddj 1-013  Removed call to BIGALH. 
C  10 Nov 99 @DynTel-ddj 1-003  Added initialization of OMTCF for 
C            diesel vehicles. 
C  28 Sep 99 @DynTel-ddj for EPA  Removed ISR from call to BIGMET. 
C  1 Mar  99 @DynTel-HQ 2-698 move the calculation of scf to efcalx, 
C                             hccalx and hrloop 
C  21 Apr 99 @DynTel-bg  0-632 Temporarily commented out call to BIGMET 
C   1 Apr 99 @DynTel-DDJ 2-696 Added hourly loop. 
C  30 Jul 98 @DynTel-zk 2-000 Explicit typing 
C 
C  30 Jul 98 @Dyntel -RJD 2-660 Disabled temperature corrections for 
C            vehicles with separated start and running emissions. 
C 
C  Input on call: 
C 
C    parameter list: IP,MYR,LIDX,IVTMP5 
C    common blocks: 
C    /CITRV2/ RVPX 
C    /FLAGS4/ COLDFG 
C    /FLAGS5/ TCF_FLAG,RVP_FLAG 
C    /INJECT/ PFI 
C    /SOAK  / SOK_TIME 
C    /TEMPC1/ TTFC,TTFD,TT4,TT8 
C    /TEMPC2/ MYGTFD,ISHIFT 
C    /TEMPC3/ MDLOHI,TDIFF 
C    /TEMPC4/ MYGCOO,LOWCO,COOTFC,COOTFD 
C    /TEMPC5/ NCC,MYCOL 
C    /TEMPC6/ CSUB,CMUL 
C    /TEMPS/  TEMEXH 
C 
C  Output on return: 
C 
C    common blocks: 
C    /OMTCOM/ TCFSTART,TCFRUN,TCF12 
C    /OFFSET/ OFFSET 
C    /TEMPC3/ TCF 
C 
C  Local array subscripts: 
C 
C  FRFDS(3)  -  FRFDS ( IFDS ) 
C  DTC1(3)   -  DTC1 (IFDS) 
C  DTC2(3)   -  DTC2 (IFDS) 
C 
C  Local variable / array dictionary: 
C 
C   Name   Type              Description 
C  ------  ----  ------------------------------------------------------- 
C  FRFDS    R    fraction of class IVTEMP given MY that is IFDS equipped 
C  HDGTCF   R    HDGV T CF - HDGV testing does not use the LDG operating modes, 
C                so multiply the 1 HDGTCF x each of the RVP CFs in TCF(IB). 
C  TC1      R    LDGV/T tcf before checked for cold CO std adjustment: 
C                no std applies = use as is, std applies = weight in std impact 
C  TC2      R    LDGV/T tcf checked & adjusted for cold CO std adjustment 
C  TFCB     R    tcf by bag for fuel delivery systems combined case 
C  TFDS     R    FDS weighted low tcf bag equation coefficient 
C  DTC1     R    TC1, but for differentiated FDS calculations 
C  DTC1     R    TC2, but for differentiated FDS calculations 
C 
C  Notes: 
C 
C  For low tcf, LDGV/T use the additive CO offset starting in MYGCOO(IVTEMP) & 
C  are differentiated by fuel delivery system (carbureted vs fuel injected) 
C  starting in MYGTFD(IVTEMP).  The CO offset start year precedes the fds start 
C  year in all 3 cases.  Hence, BIGTCF checks the LOWCO flag on both the LDG 
C  MYGTFD keyed branches. 
C 
C  BIGTCF was added for MOBILE4. 
C  BIGTCF was modified for MOBILE4, adding BIGRVP. 
C  Oxygenated and cert. fuels support from MOBILE4 was removed in MOBILE4.1. 
C  Sep-26-1994 @ CSC-bsg Request 426  Added IVTEMP parameter for calling 
C  P3WPOX, from BIGRVP. 
C  Oct-12-1994 @ CSC-bsg Request 426  Clarified documentation.  In specific, 
C  MY is output on return in MYCODE so that ITAMPT can use the value.  ITAMPT 
C  is called by P3WPOX which is called by BIGRVP.  IDX is set here to minimize 
C  code changes but called LIDX on the parameter list to avoid conflicts with 
C  Common. 
C 
C 
      IMPLICIT NONE 
      INCLUDE 'CITRV2.I' 
      INCLUDE 'FLAGS4.I' 
      INCLUDE 'FLAGS5.I' 
      INCLUDE 'INJECT.I' 
      INCLUDE 'IVPCOM.I' 
      INCLUDE 'IVTYPE.I' 
      INCLUDE 'OMTCOM.I' 
      INCLUDE 'OFFSET.I' 
      INCLUDE 'SOAK.I' 
      INCLUDE 'TEMPC1.I' 
      INCLUDE 'TEMPC2.I' 
      INCLUDE 'TEMPC3.I' 
      INCLUDE 'TEMPC4.I' 
      INCLUDE 'TEMPC5.I' 
      INCLUDE 'TEMPC6.I' 
      INCLUDE 'TEMPC7.I' 
      INCLUDE 'TEMPS.I' 
C 
      REAL TFCB(3),TFDS(3),FRFDS(3),DTC1(3),DTC2(3) 
      REAL TC1,TC2,HDGTCF,PIV,TCONST 
      INTEGER IB,IG,IVTEMP,IT,IGO,IG5,IFDS,IVRVP,IV,IVTMP5,IP, 
     *        MYR,IPHASE,ST_INDX 
      INTEGER IGRVPT,ITCPTR,IEVPTR 
C 
C  Determine the local pointer IVTEMP (see COMMENTS.I definition) 
C 
      IVTEMP=IVPTRT(IV) 
C 
C  Set the COLD CO Phase to 1 
C 
      IPHASE=1 
C 
C  Initialize TCF to the neutral (no correction) value. 
C 
      DO 10 IB=1,3 
      TCF(IB)=1.0 
   10 CONTINUE 
C 
C  First cases evaluated are the new combined high temperature and RVP cf 
C  algorithm and bag dependent low temperature RVP cf cases.  Check these first 
C  because their year ranges are mostly likely to be evaluated. 
C 
      IF(IVTMP5.GT.4) GOTO 30 
C 
      IVRVP=IVPTRR(IV) 
      IG=IGRVPT(MYR,IVRVP) 
      IF(IG.LE.1.OR.TEMEXH.LE.45.0) GOTO 30 
C 
C  Sep-26-1994 @ CSC-bsg Request 426  Added IVTEMP parameter for calling 
C  P3WPOX, from BIGRVP. 
C  Replace IVTEMP with IV so that P3WPOX is now correct. 
C  RVP Correction is applied according to the user option. 
C 
      IF(RVP_FLAG.EQ.1) CALL BIGRVP(IP,IG) 
C 
C  High temperature and RVP combined cf case is done. 
C 
      IF(MDLOHI.EQ.2) GOTO 999 
C 
C  All RVP cf calculations are done.  If T is neutral, or TCF is 
C  disabled, no T cf occurs. 
C  For the IG >= 2 case returned from BIGRVP, pass on the RVP cf in TCF. 
C 
   30 IF(MDLOHI.EQ.0.OR.TCF_FLAG.EQ.0) GOTO 999 
C 
C  At this point, TCF(IB) = 1.0, if IV = MC or IG <= 1 or (IG > 1 & T = low) 
C                           RVP cf, if IG >= 2 & T < 75.0 
C  Hence, the formulas below for the multiplicative tcf cases are all of the 
C  form "TCF(IB)=f(---)*TCF(IB)".  The result is either a tcf or the product 
C  by bag tcf * RVP cf. 
C 
      GOTO(40,40,40,70,80),IVTMP5 
C 
C  Low for all MY & high for RVP MYGs 0 & 1 TCF for LDGV, LDGT1, and LDGT2. 
C 
   40 IF(IP.EQ.2.AND.MDLOHI.EQ.1.AND.MYR.GE.MYGCOO(IVTEMP)) LOWCO=2 
      ADDCSW=1 
      MULCSW=1 
      IF(IP.EQ.2.AND.MDLOHI.EQ.1.AND.COLDFG.GT.1) CALL ICCOPT(MYR) 
      IF(MYR.GE.MYGTFD(IVTEMP)) GOTO 50 
C 
C  IVLDG multiplicative low & high tcf with FDS combined. 
C 
      IT=ITCPTR(MYR,IP,IVTEMP) 
      DO 45 IB=LOWCO,3 
      TFCB(IB)=TTFC(MDLOHI,IB,IT,IP) 
      TC1=EXP(TFCB(IB)*TDIFF) 
      IF(MULCSW.EQ.1) GOTO 43 
      TC2=((TC1-1.0)*CMUL(IVTEMP,ICOLD))+1.0 
      TC1=TC1*(1.0-PIV(MYR,IVTEMP,IPHASE))+TC2*PIV(MYR,IVTEMP,IPHASE) 
   43 TCF(IB)=TC1*TCF(IB) 
   45 CONTINUE 
      IF(LOWCO.EQ.1) GOTO 999 
C 
C  IVLDG additive low temperature CO offset with FDS combined. 
C  (For MOBILE6, removed CHFRAC from TC1 since soak correction is applied later) 
C 
      IGO=IT-ISHIFT 
      TC1=COOTFC(IGO)*TDIFF 
      IF(ADDCSW.EQ.1.AND.MULCSW.EQ.1) GOTO 47 
      IF(ADDCSW.EQ.2) TC2=TC1-CSUB(IVTEMP,ICOLD) 
      IF(MULCSW.EQ.2) TC2=TC1*CMUL(IVTEMP,ICOLD) 
      TC1=TC1*(1.0-PIV(MYR,IVTEMP,IPHASE))+TC2*PIV(MYR,IVTEMP,IPHASE) 
   47 OFFSET=TC1 
      GOTO 999 
C 
C  IVLDG multiplicative low tcf with FDS differentiated. 
C 
   50 IG5=IEVPTR(MYR,IVTEMP,2) 
      FRFDS(2)=TBI(IG5,IVTEMP) 
      FRFDS(3)=PFI(IG5,IVTEMP) 
      FRFDS(1)=1.0-FRFDS(2)-FRFDS(3) 
      DO 60 IB=LOWCO,3 
      TFDS(IB)=0.0 
      DO 55 IFDS=1,3 
      TFDS(IB)=TFDS(IB)+FRFDS(IFDS)*TTFD(IB,IP,IFDS) 
   55 CONTINUE 
      TC1=EXP(TFDS(IB)*TDIFF) 
      IF(MULCSW.EQ.1) GOTO 57 
      TC2=((TC1-1.0)*CMUL(IVTEMP,ICOLD))+1.0 
      TC1=TC1*(1.0-PIV(MYR,IVTEMP,IPHASE))+TC2*PIV(MYR,IVTEMP,IPHASE) 
   57 TCF(IB)=TC1*TCF(IB) 
   60 CONTINUE 
      IF(LOWCO.EQ.1) GOTO 999 
C 
C  IVLDG additive low temperature CO offset with FDS differentiated. 
C 
C  Calculate the uncontrolled offset (TC1) for each fuel delivery system. 
C  Also initialize the controlled offset (TC2) to neutral. Note that CHFRAC 
C  is no longer used here since soak correction is applied later. 
C 
      DO IFDS=1,3 
        DTC1(IFDS)=COOTFD(IFDS)*TDIFF 
        DTC2(IFDS)=0.0 
      END DO 
C 
C  When appropriate, compute the controlled offset (TC2) for carb and TBI 
C  fuel delivery systems only. (MOBILE5b method) 
C 
      IF(ADDCSW.EQ.1.AND.MULCSW.EQ.1) GOTO 65 
C 
      DO IFDS=1,2 
       IF(ADDCSW.EQ.2) DTC2(IFDS)=DTC1(IFDS)-CSUB(IVTEMP,ICOLD) 
       IF(MULCSW.EQ.2) DTC2(IFDS)=DTC1(IFDS)*CMUL(IVTEMP,ICOLD) 
      END DO 
C 
C  Compute the Phase 1 CO Offset for PFI LDGV and LDGT1-LDGT4. 
C  The values 3.82 and 3.24 are the Phase 1 CO Offsets at 20F for 
C  LDGV and LDGT1-LDGT4, respectively. The value 0.43 is the fraction 
C  of cold starts, 7.5 is the number of miles in the FTP, and 3.59 
C  is the FTP Bag 1 distance. (The Bag 1 distance is applied to put 
C  the DTC2 units in grams per mile, the same units as DTC1.). 
C 
      IF(IVTMP5.EQ.1) THEN 
        DTC2(3)=((3.82*7.5)/(0.43*3.59))*((75.0-TEMEXH)/(75.0-20.0)) 
      ELSE 
        DTC2(3)=((3.24*7.5)/(0.43*3.59))*((75.0-TEMEXH)/(75.0-20.0)) 
      ENDIF 
C 
C  The PFI Phase 1 CO Offset must be equal to or less than its 
C  uncontrolled counterpart. 
C 
      IF(DTC2(3).GT.DTC1(3)) DTC2(3)=DTC1(3) 
C 
C  Compute the fuel delivery system weighted average TC1 and TC2. 
C 
   65 TC1=0.0 
      TC2=0.0 
      DO IFDS=1,3 
       TC1=TC1+DTC1(IFDS)*FRFDS(IFDS) 
       TC2=TC2+DTC2(IFDS)*FRFDS(IFDS) 
      END DO 
C 
C  Finally, compute the CO offset by weighting the fraction of vehicles 
C  covered by the Cold CO Standards. 
  
      OFFSET=TC1*(1.0-PIV(MYR,IVTEMP,IPHASE)) 
     *      +TC2*PIV(MYR,IVTEMP,IPHASE) 
      GOTO 999 
C 
C 
C  Low for all MY & high for RVP MYG 0 (pre-1985) TCF for HDGV. 
C  Apply MOBILE5b BIGOMT code: HDV: a single RVP & T cf - assume 100% hot 
C  stabilized opmode = TCF's bag 2 
C 
   70 IT=ITCPTR(MYR,IP,IVTEMP) 
      HDGTCF=EXP(TT4(MDLOHI,IT,IP)*TDIFF) 
      HDGTCF=HDGTCF*TCF(2) 
      DO 75 IB=1,3 
      TCF(IB)=HDGTCF 
   75 CONTINUE 
      GOTO 999 
C 
C  Lo/hi TCF for MC. 
C 
   80 IT=ITCPTR(MYR,IP,IVTEMP) 
      DO 85 IB=1,3 
      TCF(IB)=EXP(TT8(MDLOHI,IB,IT,IP)*TDIFF) 
   85 CONTINUE 
C 
C  Compute the Multiplicative TCF for Exhaust Running Emissions by  
C  weighting the Bag 2 and Bag 3 factors. 
C 
 999  TCFRUN=TCF(2)*0.521+TCF(3)*(1.0-0.521) 
C 
C  Compute the Multiplication TCF for Basic Engine Start Emissions by 
C  weighting Bag 1 and Bag 3 by the soak time. Also store the 12-hour 
C  soak TCF for use when no soak correction is requested. 
C 
      TCONST=(TCF(1)-TCF(3))/(720.-10.) 
C 
      DO ST_INDX=1,69 
        TCFSTART(ST_INDX)=TCF(3)+TCONST*(SOK_TIME(ST_INDX)-10.) 
      END DO 
C 
      TCF12=TCF(1) 
C 
C  Convert the OFFSET from g/mi to g/start by multiplying 
C  by the FTP Bag 1 distance. 
C 
      OFFSET=OFFSET*3.59 
C 
      RETURN 
      END 
