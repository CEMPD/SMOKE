C Version "@(#)$Id$ $Source$ $Date$ 

      SUBROUTINE IDLCAL(ICY)
C
C  IDLCAL calculates the idle emission factors.
C
C  Called by MOBILE.
C
C  Calls BIGSAL, ISPPTR and SCFSUB.
C
C  Input on call:
C
C    parameter list: ICY
C    common blocks:
C    /BYMYC2/ BYBEF4,BYTAM
C    /IVPCOM/ IVPTRL,IVPTRA,IVPTRB
C    /MAXIMA/ MAXVEH,MAXPOL,MAXYRS
C    /MYRCAL/ TF
C    /SPEED4/ LB1STS
C    /SPEED6/ SCFIDL
C    /VMXCOM/ VMTMIX
C
C  Output on return:
C
C    common blocks:
C    /RESUL2/ EFIDLE
C    /RESUL3/ VIDLE
C    /SPEED4/ IVB,IVA
C    /SPEED6/ IDLYES
C
C  Local array subscripts:
C
C
C  Local variable / array dictionary:
C
C   Name   Type              Description
C  ------  ----  -----------------------------------------------------
C  IDLFAC   R    idle emission factor for IDX: summed across IDX to
C                get EFIDLE for the given IP & IV
C  VMLDGT   R    mileage for composite light duty gas truck
C
C  Notes:
C
C  IDLCAL was modified in MOBILE4.1 to include TOC, TOG and oxy fuels.
C  IDLCAL was modified for MOBILE5A in house to use the corrected BEF
C  multiplied by the speed correction factor @2.5 mph multiplied by 2.5
C  to convert to gph.  13 July 1993 @ ARC-rwk Task 2-277
C
CC    IMPLICIT NONE
C
      INTEGER ICY,ISPPTR
C
      INCLUDE 'BYMYC2.I'
      INCLUDE 'FLAGS4.I'
      INCLUDE 'IVPCOM.I'
      INCLUDE 'MAXIMA.I'
      INCLUDE 'MYRCAL.I'
      INCLUDE 'RESUL2.I'
      INCLUDE 'RESUL3.I'
      INCLUDE 'SPEED4.I'
      INCLUDE 'SPEED6.I'
      INCLUDE 'VMXCOM.I'
C
      INTEGER IP,IV,IVLM,IDX,JDX,ISP,MY
      REAL IDLFAC,VMLDGT
C
C  Build the 5th degree polynomial ("old") equation unadjusted scf table.
C
      DO 10 IV=1,MAXVEH
        CALL SCFSUB(2.5,IV)
   10 CONTINUE
C
      IDLYES=.TRUE.
C
      DO 50 IP=1,MAXPOL
        IF(PRTFLG.NE.IP.AND.PRTFLG.NE.4) GOTO 50
        VIDLE(IP)=0.0
        EFIDLE(IP,9)=0.0
C
        DO 30 IV=1,MAXVEH
          EFIDLE(IP,IV)=0.0
          IF(VMTMIX(IV).EQ.0.0) GOTO 30
          IVA=IVPTRA(IV)
          IVB=IVPTRB(IV)
          IVLM=IVPTRL(IV)
C
          DO 20 IDX=1,MAXYRS
C
C  See EFCALX explanation for this TF check.
C
            IF(TF(IDX,IV).LE.0.0) GOTO 20
C
C  Initialize emissions to zero and MY, JDX parameters.
C
            IDLFAC=0.0
            MY=ICY+IDX-MAXYRS
            JDX=MAXYRS+1-IDX
C
C  Set up speed pointer for BIGSAL which sets up SCFIDL.
C
            IF(IVLM.LE.3.AND.IVLM.NE.0.) THEN
              IF(MY.LT.LB1STS(IVLM)) ISP=ISPPTR(MY,0,IVLM,1)
            ENDIF
C
            CALL BIGSAL(IP,IV,MY,IDX,ISP)
C
C  Compute idle emissions.  Use the corrected BEF multiplied by the speed
C  correction factor @2.5 mph multiplied by 2.5 to convert to gph.
C
            IF(IV.LE.4) THEN
              IDLFAC=
     *        (BYBEF4(IP,JDX,IV)+BYTAM(IP,JDX,IV))*SCFIDL*TF(IDX,IV)*2.5
            ELSE
              IDLFAC=BYBEF4(IP,JDX,IV)*SCFIDL*TF(IDX,IV)*2.5
            ENDIF
C
C  Summate emissions.
C
            EFIDLE(IP,IV)=EFIDLE(IP,IV)+IDLFAC
   20     CONTINUE
   30   CONTINUE
C
C  Assign composite light duty gas truck emission factor.
C
        VMLDGT=VMTMIX(2)+VMTMIX(3)
        IF(VMLDGT.GT.0.0) EFIDLE(IP,9)=
     *    (EFIDLE(IP,2)*VMTMIX(2)+EFIDLE(IP,3)*VMTMIX(3))/VMLDGT
C
C  Calculate a weighted by IV idle ef (VIDLE).
C
        DO 40 IV=1,MAXVEH
          VIDLE(IP)=VIDLE(IP)+EFIDLE(IP,IV)*VMTMIX(IV)
   40   CONTINUE
C
   50 CONTINUE
C
      IDLYES=.FALSE.
C
      RETURN
      END
