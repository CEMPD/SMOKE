C Version @(#)$Id$ $Source$ $Date$ 

      SUBROUTINE IDLCAL(ICY,VMLDGT,INERR)
C
C  IDLCAL calculates the idle emission factors.
C
C  Called by EFCALX.
C
C  Calls BEFIDL, BEFID1, BEFID2, BEFID3, BEFID4, GSFIDX, IDLPTR,
C  QUITER and VOCFID.
C
C  Input on call:
C
C    parameter list: ICY,VMLDGT,INERR
C    common blocks:
C    /BASEQ6/ NUMERU
C    /BASEQ7/ IPMOD
C    /CUMCOM/ CUMMIL
C    /FLAGS2/ NEWFLG
C    /FLAGS4/ NMHFLG,PRTFLG
C    /GSFCOM/ GSFRAC
C    /IDLEQ1/ ZMLIDL,DRIDL
C    /IDLEQ3/ LDVIDL
C    /MAXIMA/ MAXVEH,MAXPOL,MAXYRS
C    /MYRCAL/ TF
C    /REGION/ IREJN
C    /VMXCOM/ VMTMIX
C
C  Output on return:
C
C    parameter list: INERR
C    common blocks:
C    /RESUL2/ EFIDLE
C    /RESUL3/ VIDLE
C
C  Local array subscripts:
C
C  IDLREV(4)  -  IDLREV ( IVREV )
C
C  Local variable / array dictionary:
C
C   Name   Type              Description
C  ------  ----  -----------------------------------------------------
C  IDLFAC   R    idle emission factor for IDX: summed across IDX to
C                get EFIDLE for the given IP & IV
C  IDLREV   I    lower bound MYs of REVised IDLe EF calculation
C                algorithm MYGs
C  IGIDL    I    pointer variable to Heavy Duty Diesel Vehicle array.
C
C
C  Notes:
C
C  IDLCAL was modified in MOBILE4.1 to include TOC, TOG and oxy fuels.
C  IDLCAL was mofified in MOBILE5 to move methane calculations to
C  BEFIDL and to break BEFIDL into several smaller functions.
C
C
      COMMON /BASEQ6/ MYGERU(12,2,3,8,2),MAXERU,NUMERU(3,8,2),KEYER,IGER
      COMMON /BASEQ7/ IPMOD(3)
      COMMON /CUMCOM/ CUMMIL(25,8)
      COMMON /FLAGS2/ MYMRFG,NEWFLG,IMFLAG,ALHFLG,IMTIER
      COMMON /FLAGS4/ PRTFLG,IDLFLG,NMHFLG,HCFLAG,COLDFG
      COMMON /GSFCOM/ MAXGSF,GSFRAC(22,8,2),MYGSF(22,8,2)
      COMMON /IDLEQ1/ ZMLIDL(15,3,8,2),DRIDL(15,3,8,2)
      COMMON /IDLEQ3/ UTIDOF(14,3,2),LDVIDL(25)
      COMMON /MAXIMA/ MAXVEH,MAXLTW,MAXPOL,MAXREG,MAXYRS
      COMMON /MYRCAL/ XMYM(25,8),JANMYR(25,8),TF(25,8),TFMYM(25,8)
      COMMON /REGION/ FEET(2),IREJN,ALT,INITPR
      COMMON /RESUL2/ RLGGAL(9),RSTGPH(9),EFRSTL(9),EFIDLE(3,9)
      COMMON /RESUL3/ VFTP(3),VEXH,VEVAP,VLOSS,VRUNLS,VRSTLS,VIDLE(3)
      COMMON /VMXCOM/ REGMIX(8),TFNORM(8),VMTMIX(8),VCOUNT(39,8)
      COMMON /IOUCOM/ IOUIMD,IOUGEN,IOUREP,IOUERR,IOUASK
C
      REAL IDLFAC
      INTEGER ALHFLG
      INTEGER PRTFLG,HCFLAG,COLDFG
      REAL JANMYR,LDVIDL
C
      DIMENSION IDLREV(4)
      DATA IDLREV/1977,1979,1979,1984/
C
      DO 10 IDX=1,MAXYRS
        LDVIDL(IDX)=0.0
   10 CONTINUE
C
      DO 70 IP=1,MAXPOL
      IF(PRTFLG.NE.IP.AND.PRTFLG.NE.4) GOTO 70
      VIDLE(IP)=0.0
      EFIDLE(IP,9)=0.0
C
      DO 50 IV=1,MAXVEH
      EFIDLE(IP,IV)=0.0
      IF(VMTMIX(IV).EQ.0.0) GOTO 50
C
      DO 30 IDX=1,MAXYRS
C
C  See EFCALX explanation for this TF check.
C
      IF(TF(IDX,IV).LE.0.0) GOTO 30
C
C  Initialize emissions to zero and figure parameters.
C
      IDLFAC=0.0
C
C  If pollutant IP -- vehicle class IV has one or more sets of
C  replacement BEF equation coefficients, the idle emission factor
C  for that IP/IV combination is set as zero.  This is because
C  IF((ICY-IDLREV(IVLDG)+1).GE.25) THEN all MYs of ICY's 25 MY
C  window use the revised idle algorithm.  In the future it can
C  be reprogrammed if anyone can confidently produce an algorithm
C  that can be shown to handle NEWFLG for ICYs passing the IF.
C  Note: at the present time IDLREV ( 1 ) = 1977; ( 2/3 ) = 1979
C  and (4) is 1984.
C
      IF(NUMERU(IP,IV,IREJN).GE.1) GOTO 30
C
C  Initialize variables JDX counter, MY and VMT age.
C
      JDX=MAXYRS-IDX+1
      VMTAGE=CUMMIL(JDX,IV)/10000.
      MY=ICY+IDX-MAXYRS
C
C  Heavy Duty Diesel Exit for Idle emissions calculations.
C
      IF(IV.EQ.7) THEN
         IGIDL=IDLPTR(MY,IP,IV)
         IDLFAC=ZMLIDL(IGIDL,IP,IV,IREJN)
     *         +DRIDL(IGIDL,IP,IV,IREJN)*VMTAGE
         IDLFAC=IDLFAC*TF(IDX,IV)
        GOTO 20
      ENDIF
C
C  Compute idle emissions.
C
      IF(IV.LT.4) THEN
C
C  NOx gets a special consideration
C
         IF(IP.EQ.3.AND.MY.GE.IDLREV(IV)) THEN
            IDLFAC=BEFID1(MY,IDX,IP,IV)*TF(IDX,IV)
            GOTO 20
         ENDIF
C
C  HC, CO LDGV and LDGT are the original BEFIDL
C
         IF(MY.GE.IDLREV(IV)) THEN
            IDLFAC=BEFIDL(MY,IDX,IP,IV,VMTAGE)*TF(IDX,IV)
            GOTO 20
         ENDIF
C
C  Cars and Trucks before IDLREV(IV) are in the same group.
C
         IDLFAC=BEFID2(MY,IDX,IP,IV,VMTAGE)*TF(IDX,IV)
         GOTO 20
      ENDIF
C
C  HDGV get a special consideration
C
      IF(IV.EQ.4.AND.IP.EQ.1) THEN
         IDLFAC=BEFID3(MY,IDX,VMTAGE)*TF(IDX,IV)
         GOTO 20
      ENDIF
      IF(IV.EQ.4.AND.IP.EQ.2) THEN
         IDLFAC=BEFID5(MY,IDX,VMTAGE)*TF(IDX,IV)
         GOTO 20
      ENDIF
      IF(IV.EQ.4.AND.IP.EQ.3) THEN
         IDLFAC=BEFID6(MY,IDX,VMTAGE)*TF(IDX,IV)
         GOTO 20
      ENDIF
C
C  Special consideration, if IV=6, 7 OR 8.
C
      IDLFAC=BEFID4(MY,IDX,IP,IV,VMTAGE)*TF(IDX,IV)
C
C  Summate emissions.
C
   20 EFIDLE(IP,IV)=EFIDLE(IP,IV)+IDLFAC
   30 CONTINUE
C
C  NEWFLG versus calculated emissions check.  Improve by checking all 25 IDX
C  for use/non-use of NEWFLG & err exit only on non-use of alternate equations.
C  Skip for IV for which MOBILE4 does not store gas-diesel combined
C  registrations (splitting them with GSFRAC during execution).
C
      IF(IV.EQ.4.OR.IV.GE.7) GOTO 40
      IF(NEWFLG.EQ.1.AND.EFIDLE(IP,IV).GT.0.
     *  .AND.GSFIDX(ICY,MY,IV,IREJN).EQ.0.0) CALL QUITER(0.,IV,69,INERR)
   40 IF(NEWFLG.EQ.1.AND.EFIDLE(IP,IV).LE.0.0)
     *   CALL QUITER(0.,IV,70,INERR)
C
   50 CONTINUE
C
C  Assign composite light duty gas truck emission factor.
C
      IF(EFIDLE(IP,2).GT.0.0.AND.EFIDLE(IP,3).GT.0.0) EFIDLE(IP,9)=
     *   (EFIDLE(IP,2)*VMTMIX(2)+EFIDLE(IP,3)*VMTMIX(3))/VMLDGT
C
C  If one or more of the IV for a given IP have had 1 or more replacement basic
C  FTP ef parameters records accepted via NEWFLG and therefore had all the
C  corresponding idle ef parameters zeroed out, then it does not make sense to
C  calculate a weighted by IV idle ef (VIDLE).
C
      IF(IPMOD(IP).EQ.1) GOTO 70
      DO 60 IV=1,MAXVEH
      VIDLE(IP)=VIDLE(IP)+EFIDLE(IP,IV)*VMTMIX(IV)
   60 CONTINUE
   70 CONTINUE
C
      RETURN
      END
