C Version "@(#)$Id$ $Source$ $Date$ 

      SUBROUTINE HCCALX(ICY,VMLDGT,INERR,IY)
C
C  HCCALX calculates HC emission factors.
C
C  Called by EFCALX.
C
C  Calls BEF, CCEVRT, EFNFTP, FAIL, GSFIDX, IEVPTR, PCTLEV, QUITER,
C  RSTLOS, and RULOSS.
C
C  Input on call:
C
C    parameter list: ICY,VMLDGT,INERR
C    common blocks:
C    /CEVBMY/ BMYMPD
C    /CITCIN/ PFRATE,IPFCU
C    /CITRV2/ RVPHS,
C    /CUMCOM/ CUMMIL
C    /EVAPGR/ PARD
C    /FLAGS3/ RLFLAG
C    /FLAGS4/ HCFLAG
C    /GSFCOM/ GSFRAC
C    /MAXIMA/ MAXVEH,MAXYRS
C    /MYRCAL/ JANMYR,TF
C    /REGION/ IREJN
C    /RLCOM1/ ROADFE
C    /RLCOM3/ RLRATE
C    /RVPEX2/ RVPCF
C    /SPEED6/ SALHCF
C    /TEMPS/  TEMEVP
C    /VMXCOM/ REGMIX,VMTMIX
C
C  Output on return:
C
C    parameter list: INERR
C    common blocks:
C    /BYMYC2/ BYBEF4,BYTAM,BYFER
C    /BYMYC3/ BYEVAP,BYRUNL,BYREFL,BYRSTL
C    /BYMYC5/ BYIMRE,REDSUM
C    /CITCIN/ PFRATE
C    /RESUL1/ EFFTP,EFEXH,EFEVAP,EFREFL,EFRUNL
C    /RESUL2/ RLGGAL,RSTGPH,EFRSTL
C    /RESUL3/ VFTP,VEXH,VEVAP,VLOSS,VRUNLS
C    /EVAPGR/ GREVP,VGREVP,GRPD
C
C  Local array subscripts
C
C  EVPSUM ( 4 ) - EVPSUM ( IEVP )
C
C  Local variable / array dictionary:
C
C   Name   Type              Description
C  ------  ----  -------------------------------------------------------
C  COMPCC   R    composite evaporative & crankcase emission factor for IDX
C  COMPEF   R    composite emission factor for IDX: summed across IDX to get
C                EFFTP for the given IP & IV
C  EVPSUM   R    gram evap component summation variable
C  EXHHC    R    exhaust HC component of total HC, for given IV x IDX case
C  EXHWGT   R    multiplicative weights to be applied to the exh ef returned
C                by BEF: speed & "BIGALH" cf, travel fraction & RVP cf
C  GASCAP   R    refueling HC loss factor for given IDX and IV
C  MDX      I    model year index into road fuel economy arrays
C  RNGLOS   R    incremental running loss HC
C  RSLOSS   R    incremental resting loss HC
C  RVOPCF   R    exhaust HC & CO cf for RVP
C                given vehicle type and model year
C  VMTAGE   R    see parameter dictionary
C  Z        R    Discount evap emission by ZEV fraction.
C
C  Notes:
C
C  HCCALX was modified for MOBILE4 In-House Version 02 to save by model year
C  computations for the BYMY Tables.
C  HCCALX was modified for MOBILE4.1 to compute grams/mile hot soak,
C  diurnal, and partial diurnal values.  Also, oxygenated and cert.
C  fuels support from MOBILE4 was removed in MOBILE4.1.
C
C  (22-April-1994) @ CSC ljn Request 2-420 added LDGT2 (IV=3) to the LEV's
C  Aug-18-1994 @ CSC-bsg Request 462 Added calls to non-FTP emission factor 
C  calculation function.
C  Aug-25-1994 @ CSC-bsg Request 462 Made correction in non-FTP BYBEF4 
C  calculation to remove improperly making a SALHCF adjustment to by-model
C  Sep-19-1994 @ CSC-bsg Request 462 Made correction in HC BYMY FER
C  calculation, to properly account for tampering in the result printed.
C  Sep-21-1994 @ CSC-ked (for bsg) Request 438 Fix 1990 and earlier calendar year
C  removed IY from the call to CCEVRT
C
      INCLUDE 'BYMYC2.I'
      INCLUDE 'BYMYC3.I'
      INCLUDE 'BYMYC5.I'
      INCLUDE 'CEVBMY.I'
      INCLUDE 'CITCIN.I'
      INCLUDE 'CITRV1.I'
      INCLUDE 'CITRV2.I'
      INCLUDE 'CUMCOM.I'
      INCLUDE 'EGSCAL.I'
      INCLUDE 'EVAPGR.I'
      INCLUDE 'FLAGS3.I'
      INCLUDE 'FLAGS4.I'
      INCLUDE 'GSFCOM.I'
      INCLUDE 'LEVBLK.I'
      INCLUDE 'MAXIMA.I'
      INCLUDE 'MYRCAL.I'
      INCLUDE 'RLCOM1.I'
      INCLUDE 'RLCOM3.I'
      INCLUDE 'REGION.I'
      INCLUDE 'RESUL1.I'
      INCLUDE 'RESUL2.I'
      INCLUDE 'RESUL3.I'
      INCLUDE 'RVPEX2.I'
      INCLUDE 'SPEED6.I'
      INCLUDE 'TEMPS.I'
      INCLUDE 'VMXCOM.I'
C
      DIMENSION EVPSUM(4)
C
C  Initialize
C
      EFEXH(9)=0.0
      EFEVAP(9)=0.0
      EFREFL(9)=0.0
      RLGGAL(9)=0.0
      RSTGPH(9)=0.0
      EFRUNL(9)=0.0
      EFRSTL(9)=0.0
      VEXH=0.0
      VEVAP=0.0
      VLOSS=0.0
      VRUNLS=0.0
      VRSTLS=0.0
      DO 10 IEVP=1,4
      GREVP(IEVP,9)=0.0
      VGREVP(IEVP)=0.0
   10 CONTINUE
C
      IP=1
      EFFTP(IP,9)=0.0
      VFTP(IP)=0.0
C
      DO 12 IDU=2,4
      IPARTL=IDU-1
      GRPD(IPARTL,9)=0.0
   12 CONTINUE
C
C  For each vehicle class, calculate exhaust, evap and partial diurnals,
C  and running loss values
C
      DO 85 IV=1,MAXVEH
      EFFTP(IP,IV)=0.0
      IF(IV.LE.4) REDSUM(IP,IV)=0.0
      EFEXH(IV)=0.0
      EFEVAP(IV)=0.0
      EFREFL(IV)=0.0
      RLGGAL(IV)=0.0
      RSTGPH(IV)=0.0
      EFRUNL(IV)=0.0
      EFRSTL(IV)=0.0
      DO 15 IEVP=1,4
      GREVP(IEVP,IV)=0.0
   15 CONTINUE
C
      DO 16 IDU=2,4
      IPARTL=IDU-1
      GRPD(IPARTL,IV)=0.0
   16 CONTINUE
C
C  REGMIX is not affected by scenario input, only by BLOCK DATA / one-time
C  data values.  BYMY save areas can be initialized after this check,
C  since it would be technically in error to have a VMTMIX(i) > 0.0
C  when REGMIX(i) = 0.0, so the situation described in EFCALX for
C  CO and NOX can not happen in HCCALX for HC.
C
      IF(REGMIX(IV).EQ.0.0.AND.VMTMIX(IV).EQ.0.0) GOTO 85
      MY=0
C
      DO 70 IDX=1,MAXYRS
      JDX=(MAXYRS+1)-IDX
C
      BYBEF4(IP,JDX,IV)=0.0
      BYFER(IP,JDX,IV)=0.0
      BYEVAP(JDX,IV)=0.0
      BYRSTL(JDX,IV)=0.0
      IF(IV.GE.5) GOTO 18
      BYTAM(IP,JDX,IV)=0.0
      BYIMCR(IP,JDX,IV)=0.0
      BYIMRE(IP,JDX,IV)=0.0
      BYRUNL(JDX,IV)=0.0
      BYREFL(JDX,IV)=0.0
C
   18 IF(TF(IDX,IV).LE.0.0.AND.JANMYR(JDX,IV).EQ.0.0) GOTO 70
C
      COMPEF=0.0
      EXHHC=0.0
      COMPCC=0.0
      GASCAP=0.0
      RNGLOS=0.0
      RSLOSS=0.0
C
      DO 20 IEVP=1,4
      EVPSUM(IEVP)=0.0
   20 CONTINUE
C
      VMTAGE=CUMMIL(JDX,IV)/10000.
      MY=ICY+IDX-MAXYRS
C
C  Combine the speed and aluh correction factor and the travel fraction
C  into a single weight to correct the exhaust ef returned from BEF.
C  Aug-18-1994 @ CSC-bsg Request 462 Added calls to non-FTP emission factor
C  calculation function.  So had to move TF adjustment to end of HC exhaust 
C  section.
C
      EXHWGT=SALHCF(IDX,IP,IV)
C
C  Apply the exhaust weight to the BEF's computed ef.
C  BEF is used to get the operating mode corrections.
C
      EXHHC=BEF(MY,IDX,ICY,IP,IV,VMTAGE,IY)*EXHWGT
C
      IF(IV.LE.4) BYIMRE(IP,JDX,IV)=BYIMRE(IP,JDX,IV)*EXHWGT*TF(IDX,IV)
      Z = 1.0
C
C  (22-April-1994) @ CSC ljn Request 2-420 added LDGT2 (IV=3) to the LEV's
C
      IF(MY.GE.1994.AND.IV.LE.3.AND.(LEVFLG.EQ.2).AND.MY.GE.LEVYRS) THEN
        Z=1.0-PCTLEV(MY,IV,8)
      ENDIF
C
C  Exhaust HC ef is complete for diesels.
C
      GOTO(30,30,30,30,35,35,35,35),IV
C
C  For LDGV/T or HDGV, calculate and apply HC RVP exhaust correction
C  factor.  This cf is 1.0 in Release MOBILE4, except when IDX is
C  in RVP MYG 1 and RVP.NE.9.0
C
   30 RVOPCF=RVPCF(IDX,IP,IV)
      BYBEF4(IP,JDX,IV)=BYBEF4(IP,JDX,IV)*RVOPCF
      BYTAM(IP,JDX,IV)=BYTAM(IP,JDX,IV)*RVOPCF
      BYIMRE(IP,JDX,IV)=BYIMRE(IP,JDX,IV)*RVOPCF
      REDSUM(IP,IV)=REDSUM(IP,IV)+BYIMRE(IP,JDX,IV)
      EXHHC=EXHHC*RVOPCF
C
C  Aug-18-1994 @ CSC-bsg Request 462 Added calls to non-FTP emission factor
C  calculation function.  So had to move TF adjustment to end of HC exhaust
C  section.
C  Aug-25-1994 @ CSC-bsg Request 462 Made correction in non-FTP BYBEF4 
C  calculation to remove improperly making a SALHCF adjustment to by-model
C  Sep-19-1994 @ CSC-bsg Request 471 Made correction in HC BYMY FER 
C  calculation, to properly account for tampering in the result printed.
C
   35 BYBEF4(IP,JDX,IV)=BYBEF4(IP,JDX,IV)+EFNFTP(MY,IV,IP)
      EXHHC=(EXHHC+EFNFTP(MY,IV,IP))*TF(IDX,IV)
      EFEXH(IV)=EFEXH(IV)+EXHHC
      COMPEF=COMPEF+EXHHC
C
C  Diesels (IV=5-7) do not have evaporative, on-board or RVP
C  components, so branch accordingly.
C
      IF(IV.GE.5.AND.IV.LE.7) GOTO 65
C
C  Compute the rates of Failed Purge, Failed Pressure, and then Passed.
C  USE USER-SUPPLIED RATES WHEN NECESSARY
C
      IF(IPFCU.EQ.1) THEN
      DO 39 IK=1,2
      PFRATE(IK+1,IPFCU)=FAIL(MY,IDX,IV,IK)
   39 CONTINUE
      ENDIF
      PFRATE(1,IPFCU)=1.0-PFRATE(2,IPFCU)-PFRATE(3,IPFCU)
C
C  Now figure composite evaporative and crankcase factor for
C  current my, including any onboard effect.
C
      COMPCC=CCEVRT(MY,IDX,IV)
C
C  Call to CCEVRT also computes gram evap components.
C
      DO 40 IEVP=1,4
      EVPSUM(IEVP)=EVP(IEVP)
   40 CONTINUE
C
C  Calculate resting loss emissions for gas fueled vehicles
C
      RSLOSS=RSTLOS(MY,IV,TEMEVP(6))*Z
C
C  Evap. Test Procedure
C
      IF(IV.LE.4) RSLOSS=RSLOSS*(1.0-(ETPCAL(7,MY,1,0.,0.)+
     *            ETPCAL(7,MY,2,0.,0.)+ETPCAL(7,MY,3,0.,0.))/3.)
      IF(HCFLAG.EQ.3.OR.HCFLAG.EQ.5)
     *   RSTGPH(IV)=RSTGPH(IV)+RSLOSS*JANMYR(JDX,IV)
      IF(BMYMPD(JDX,IV).GT.0.0) RSLOSS=RSLOSS*24./BMYMPD(JDX,IV)
      IF(BMYMPD(JDX,IV).EQ.0.0) RSLOSS=0.0
      BYRSTL(JDX,IV)=RSLOSS
      RSLOSS=RSLOSS*TF(IDX,IV)
      EFRSTL(IV)=EFRSTL(IV)+RSLOSS
      COMPEF=COMPEF+RSLOSS
C
C  Weight by travel or registration fractions and summate.
C  Use January registration weighting for all grams output,
C  use travel fraction weighting for grams per mile output.
C  (Note: hot soak in g/trip mode has already be weighted by
C  JANMYR in the calculation of trips per day fractions.)
C
      BYEVAP(JDX,IV)=COMPCC
      COMPCC=COMPCC*TF(IDX,IV)
      EFEVAP(IV)=EFEVAP(IV)+COMPCC
      COMPEF=COMPEF+COMPCC
C
      IF(IV.EQ.8) GOTO 50
C
C  Calculate running loss emissions for gas fueled vehicles.
C  Motorcycles do not have running loss emissions.
C
      RNGLOS=RULOSS(MY,IV,RVPHS)*Z
C
C  Obtain corrected running losses for Oxy Fuel
C
      CALL FUEL2(3,RNGLOS)
C
      RNGLOS = RNGLOS
      BYRUNL(JDX,IV)=RNGLOS
      RNGLOS=RNGLOS*TF(IDX,IV)
      EFRUNL(IV)=EFRUNL(IV)+RNGLOS
      COMPEF=COMPEF+RNGLOS
C
   50 DO 60 IEVP=1,4
      IF(HCFLAG.EQ.4) THEN
         EVPSUM(IEVP)=EVPSUM(IEVP)*TF(IDX,IV)
      ELSE
         IF(IEVP.EQ.2.OR.IEVP.EQ.3) THEN
            EVPSUM(IEVP)=EVPSUM(IEVP)*JANMYR(JDX,IV)
         ELSEIF(IEVP.EQ.4) THEN
            EVPSUM(IEVP)=EVPSUM(IEVP)*TF(IDX,IV)
         ENDIF
      ENDIF
      GREVP(IEVP,IV)=GREVP(IEVP,IV)+EVPSUM(IEVP)
   60 CONTINUE
C
C  Set partial diurnals for LDGV, LDGT1/2. As before,
C  use January registration weighting for all grams output,
C  use travel fraction weighting for grams per mile output.
C
      IF(IV.GT.3) GOTO 64
      DO 62 IDU=2,4
      IPARTL=IDU-1
      IF(HCFLAG.NE.4)
     *     GRPD(IPARTL,IV)=GRPD(IPARTL,IV)+(PARD(IPARTL)*JANMYR(JDX,IV))
      IF(HCFLAG.EQ.4)
     *     GRPD(IPARTL,IV)=GRPD(IPARTL,IV)+(PARD(IPARTL)*TF(IDX,IV))
   62 CONTINUE
C
C  Refueling losses (GASCAP) only occur for gasoline vehicles (IV = 1-4)
C  under one of the nonzero refueling loss options (RLFLAG = 1-4).
C  For HCFLAG=3 case, do not include refueling losses in composite
C  emission factor or on refueling loss g/mi line. They will be
C  output in g/gal section (see below).
C
   64 IF(IV.GT.4.OR.RLFLAG.EQ.5) GOTO 65
      BYREFL( JDX,IV ) = RLRATE(IDX,IV) * Z
      GASCAP=(RLRATE(IDX,IV)*Z) * TF(IDX,IV)
C
C
C SET EFRREL unconditionally now, so the numbers will appear in the
C output, but don't adjust COMPEF when HCFLAG is 3. 7 Jan 93 - MLA
C
      EFREFL(IV)=EFREFL(IV)+GASCAP
      IF(HCFLAG.NE.3) THEN
        COMPEF=COMPEF+GASCAP
      ENDIF
C
C  Determine the refueling loss in g/gal by using the appropriate
C  fuel economy (only for HCFLAG=3 or 5).
C
      IF(HCFLAG.NE.3.AND.HCFLAG.NE.5) GOTO 65
C
      IF(MY.LE.1969)MDX=1
      IF(MY.GT.1969.AND.MY.LT.2000)MDX=MY-1969+1
      IF(MY.GE.2000)MDX=32
C
      RLGGAL(IV)=RLGGAL(IV)+GASCAP*ROADFE(MDX,IV)
C
C  Update FTP composite HC.
C
   65 BYFER(IP,JDX,IV)=COMPEF
      EFFTP(IP,IV)=EFFTP(IP,IV)+COMPEF
   70 CONTINUE
C
C  Sales fractions = 1.0 for all 4 HDV (HDG,HDD x L,H alt) groups and all MC's.
C
      IF(MY.EQ.0.OR.IV.EQ.4.OR.IV.GE.7) GOTO 75
      IF(EFFTP(IP,IV).GT.0.0.AND.GSFIDX(ICY,MY,IV,IREJN).EQ.0.0)
     *   CALL QUITER(0.,IV,67,INERR)
   75 IF(EFFTP(IP,IV).LE.0.0) CALL QUITER(0.,IV,68,INERR)
C
C  Weight by VMTMIX for mileage components and REGMIX for gram components.
C
      VFTP(IP)=VFTP(IP)+EFFTP(IP,IV)*VMTMIX(IV)
      VEXH=VEXH+EFEXH(IV)*VMTMIX(IV)
      VEVAP=VEVAP+EFEVAP(IV)*VMTMIX(IV)
      VLOSS=VLOSS+EFREFL(IV)*VMTMIX(IV)
      VRUNLS=VRUNLS+EFRUNL(IV)*VMTMIX(IV)
      VRSTLS=VRSTLS+EFRSTL(IV)*VMTMIX(IV)
C
      DO 80 IEVP=1,3
      VGREVP(IEVP)=VGREVP(IEVP)+GREVP(IEVP,IV)*REGMIX(IV)
   80 CONTINUE
      VGREVP(4)=VGREVP(4)+GREVP(4,IV)*VMTMIX(IV)
   85 CONTINUE
C
C  Assign composite light duty gas truck emission factors.
C
      IF(EFFTP(IP,2).GT.0.0.AND.EFFTP(IP,3).GT.0.0) EFFTP(IP,9)=
     *   (EFFTP(IP,2)*VMTMIX(2)+EFFTP(IP,3)*VMTMIX(3))/VMLDGT
      IF(EFEXH(2).GT.0.0.AND.EFEXH(3).GT.0.0)
     *   EFEXH(9)=(EFEXH(2)*VMTMIX(2)+EFEXH(3)*VMTMIX(3))/VMLDGT
      IF(EFEVAP(2).GT.0.0.AND.EFEVAP(3).GT.0.0)
     *   EFEVAP(9)=(EFEVAP(2)*VMTMIX(2)+EFEVAP(3)*VMTMIX(3))/VMLDGT
      IF(EFREFL(2).GT.0.0.AND.EFREFL(3).GT.0.0)
     *   EFREFL(9)=(EFREFL(2)*VMTMIX(2)+EFREFL(3)*VMTMIX(3))/VMLDGT
      IF(RLGGAL(2).GT.0.0.AND.RLGGAL(3).GT.0.0)
     *   RLGGAL(9)=(RLGGAL(2)*VMTMIX(2)+RLGGAL(3)*VMTMIX(3))/VMLDGT
      IF(RSTGPH(2).GT.0.0.AND.RSTGPH(3).GT.0.0)
     *   RSTGPH(9)=(RSTGPH(2)*VMTMIX(2)+RSTGPH(3)*VMTMIX(3))/VMLDGT
      IF(EFRUNL(2).GT.0.0.AND.EFRUNL(3).GT.0.0)
     *   EFRUNL(9)=(EFRUNL(2)*VMTMIX(2)+EFRUNL(3)*VMTMIX(3))/VMLDGT
      IF(EFRSTL(2).GT.0.0.AND.EFRSTL(3).GT.0.0)
     *   EFRSTL(9)=(EFRSTL(2)*VMTMIX(2)+EFRSTL(3)*VMTMIX(3))/VMLDGT
C
      DO 90 IEVP=1,3
      IF(GREVP(IEVP,2).GT.0.0.AND.GREVP(IEVP,3).GT.0.0)
     *   GREVP(IEVP,9)=(GREVP(IEVP,2)*REGMIX(2)
     *                 +GREVP(IEVP,3)*REGMIX(3))
     *                 /(REGMIX(2)+REGMIX(3))
   90 CONTINUE
      IF(GREVP(4,2).GT.0.0.AND.GREVP(4,3).GT.0.0)
     *   GREVP(4,9)=(GREVP(4,2)*VMTMIX(2)
     *              +GREVP(4,3)*VMTMIX(3))
     *              /VMLDGT
C
      DO 100 IDU=2,4
      IPARTL=IDU-1
      IF(GRPD(IPARTL,2).LE.0.0.AND.GRPD(IPARTL,3).LE.0.0) GOTO 100
      IF(HCFLAG.NE.4)
     *   GRPD(IPARTL,9)=(GRPD(IPARTL,2)*REGMIX(2)
     *                 + GRPD(IPARTL,3)*REGMIX(3))
     *                 /(REGMIX(2)+REGMIX(3))
      IF(HCFLAG.EQ.4)
     *   GRPD(IPARTL,9)=(GRPD(IPARTL,2)*VMTMIX(2)
     *                 + GRPD(IPARTL,3)*VMTMIX(3))
     *                 /(VMTMIX(2)+VMTMIX(3))
  100 CONTINUE
C
      RETURN
      END
