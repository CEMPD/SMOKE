C Version "@(#)$Id$ $Source$ $Date$ 

      SUBROUTINE BIGIN2(IV,MY,JDX,IDX,IVLM,IVTEMP,ISP)
C
C  BIGIN2 initializes the vehicle class specific parameters and flags used
C  subsequently by BIGCFX in the calculation of emission factor correction
C  factors and offsets.
C
C  Called by BIGCFX.
C
C  Calls ISPPTR.
C
C  Input on call:
C
C    constant: NEYR
C    parameter list: IV,MY,JDX,IDX
C    common blocks:
C    /CUMCOM/ CUMMIL
C    /IVPCOM/ IVPTRL,IVPTRT,IVPTRA,IVPTRB
C    /SCENE2/ FCC,FCN,FHC,FHN
C    /SPEED4/ LB1STS
C    /TAMEQ2/ F50K
C
C  Output on return:
C
C    parameter list: IVLM,IVTEMP,ISP
C    common blocks:
C    /OPMOD2/ MYSTEP
C    /OPMOD3/ BFMILE,KINKBF,BFGT50
C    /SPEED4/ IVB,IVA
C    /SPEED6/ SALHCF,SCFIDL
C    /TEMPC3/ CHFRAC,TCF
C
C  Local array subscripts:
C
C  NCATYR(8)       -  NCATYR ( IV )
C
C  Local variable / array dictionary:
C
C   Name   Type              Description
C  ------  ----  -------------------------------------------------------
C  NCATYR   I    last non-catalyst MY for the given IV
C
C  Notes:
C
C  SALHCF is reinitialized to 0.0, starting in Version 02.
C  BIGIN2 was modified for MOBILE4.1.
C
C
      INCLUDE 'CUMCOM.I'
      INCLUDE 'IVPCOM.I'
      INCLUDE 'MAXIMA.I'
      INCLUDE 'OPMOD2.I'
      INCLUDE 'OPMOD3.I'
      INCLUDE 'PARAME.I'
      INCLUDE 'SCENE2.I'
      INCLUDE 'SPEED4.I'
      INCLUDE 'SPEED6.I'
      INCLUDE 'TAMEQ2.I'
      INCLUDE 'TEMPC3.I'
C
      DIMENSION NCATYR(8)
      DATA NCATYR/1974,1974,1978,1983,4*NEYR/
C
      IVLM=IVPTRL(IV)
      IVTEMP=IVPTRT(IV)
      IVA=IVPTRA(IV)
      IVB=IVPTRB(IV)
      IF(IVLM.EQ.0) RETURN
C
C  Reinitialize TCF to neutral (no effect) value for diesel IV (Wipe HDGV NOx
C  TCF, else it would be used for IP=1-3 / IV=5-7, for which TCF is not calc.).
C
      IF(IV.NE.5) GOTO 20
      DO 10 IB=1,3
      TCF(IB)=1.0
   10 CONTINUE
C
   20 GOTO(30,30,30,40,40,40),IVLM
   30 IF(MY.LT.LB1STS(IVLM)) ISP=ISPPTR(MY,0,IVLM,1)
   40 BFMILE=CUMMIL(JDX,IV)*.0001
C
C  Select bag coefficients based on vehicle class and model year
C
      MYSTEP=1
      IF(MY.LT.1981.OR.IVLM.GT.3) GOTO 50
      IF(MY.LT.1984.AND.IVLM.GT.1) GOTO 50
C
C  Use yearly stepped bag coefficients for 1981+ LDGV and 1984+ LDGT1/2.
C  1981+ LDGV class has bag fraction coefficients curve kink at 50K miles.
C
C  1981+ LDGV and 1984+ LDGT1/2 have kinks at 50+K miles.
C
      MYSTEP=2
      KINKBF=1
      IF(IVLM.LE.3.AND.BFMILE.GT.F50K) KINKBF=2
      IF(KINKBF.EQ.2) BFGT50=BFMILE-F50K
C
C  Assign cold/hot fractions by model year group.
C
   50 CHFRAC(1)=FCN
      CHFRAC(3)=FHN
      IF(MY.LE.NCATYR(IV)) GOTO 60
      CHFRAC(1)=FCC
      CHFRAC(3)=FHC
   60 CHFRAC(2)=1.-CHFRAC(1)-CHFRAC(3)
C
C  Reinitialize SALHCF for the given IDX x IV case to 0.0.  This is done so
C  that the prior scenario's SALHCF does not show up in the current
C  scenario's (a possibility with the TF=0.0 bypass of the BIGSAL call).
C  Also, initialize Idle Speed Correction Factor Adjustments to 1.0.
C
      IF(IDLYES) THEN
        SCFIDL=0.0
      ELSE
        DO 70 IP=1,3
          SALHCF(IDX,IP,IV)=0.0
   70   CONTINUE
      ENDIF
C
      RETURN
      END
