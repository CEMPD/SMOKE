C Version "@(#)$Id$ $Source$ $Date$ 

      SUBROUTINE BIGALH(ICY)
C
C  BIGALH computes emission correction factors for air conditioning, extra
C  loading, trailer towing and humidity for vehicle types  LDGV, LDGT1, LDGT2
C  and LDGT (composite of LDGT1 & LDGT2).  Additionally, the humidity
C  correction is applied to MC vehicles.
C
C  Called by BIGCFX.
C
C  Calls IGEPTR, IGLPTR, IGTPTR and IGUPTR.
C
C  Input on call:
C
C    parameter list: ICY
C    common blocks:
C    /ALUHAC/ ACEQIP,ACCF
C    /ALUHIN/ XLOAD,TRAILR,ABSHUM,DB,WB
C    /ALUHLT/ XLCF,TWCF
C    /MAXIMA/ MAXLTW,MAXPOL,MAXYRS
C    /REGION/ IREJN
C    /TEMPS/  TEMEXH
C
C  Output on return:
C
C    common blocks:
C    /ALUHIN/ AC
C    /ALUOUT/ ALHRET
C
C  Local array subscripts:
C
C  TDIFF(3)  -  TDIFF ( IP )
C
C  Local variable dictionary:
C
C   Name   Type              Description
C  ------  ----  -------------------------------------------------------
C  ACU      R    fraction of vehicles that are equipped with air conditioning
C                and are using it
C  DI       R    discomfort index computed from the user supplied (default)
C                values of dry bulb (db) and wet bulb (wb) temperatures
C  DIHI     R    high reference temperature for DI: above this temperature
C                AC is on (100% use)
C  DILO     R    low reference temperature for DI: below this temperature
C                AC is off (0% use)
C  HH       R    HC or CO => 1.; NOx => humidity correction factor
C  HHH      R    HC or CO => 1.; NOx => humidity correction factor
C  TDIFF    R    deviation from FTP standard temperature (75 degrees F)
C  TWU      R    fraction of vehicles towing a trailer
C  XLU      R    fraction of vehicle s with an extra load
C
C  Notes:
C
C  None.
C
C
      INCLUDE 'ALUHAC.I'
      INCLUDE 'ALUHIN.I'
      INCLUDE 'ALUHLT.I'
      INCLUDE 'ALUOUT.I'
      INCLUDE 'MAXIMA.I'
      INCLUDE 'REGION.I'
      INCLUDE 'SCENE1.I'
      INCLUDE 'TEMPS.I'
C
      DIMENSION TDIFF(3)
C
      DATA DILO/70./
      DATA DIHI/80./
C
      IF(AC.EQ.0.0) GOTO 10
C
C  Compute AC based on discomfort index DI.
C
      DI=(DB+WB)*.4+15.
      AC=(DI-DILO)/(DIHI-DILO)
      IF(AC.LT.0.0) AC=0.0
      IF(AC.GT.1.) AC=1.
C
C  Setup humidity correction factors.
C
   10 DO 20 IP=1,MAXPOL
      TDIFF(IP)=TEMEXH(IP)-75.
   20 CONTINUE
C
      HHH=1.-0.00380*(ABSHUM-75.)
C
C  Assign MC only for IP=3.
C
      DO 30 IDX=1,MAXYRS
      ALHRET(IDX,4,1)=1.
      ALHRET(IDX,4,2)=1.
      ALHRET(IDX,4,3)=HHH
   30 CONTINUE
C
C  Assign other vehicle correction factors.
C
      DO 50 IVA=1,MAXLTW
      XLU=XLOAD(IVA)
      TWU=TRAILR(IVA)
C
      DO 40 IDX=1,MAXYRS
      MY=ICY-(MAXYRS-IDX)
      IGU=IGUPTR(MY,IVA)
      IGE=IGEPTR(MY,IVA)
      IGL=IGLPTR(MY,IVA)
      IGT=IGTPTR(MY,IVA)
      ACU=ACEQIP(IGE,IVA,IREJN)*AC
C
      DO 40 IP=1,MAXPOL
      HH=1.
      IF(IP.EQ.3) HH=HHH
      ALHRET(IDX,IVA,IP)=
     *   ((ACCF(2,IP,IGU,IVA,IREJN)*TDIFF(IP)
     *    +ACCF(1,IP,IGU,IVA,IREJN)-1.)*ACU+1.)
     *  *((XLCF(IP,IGL,IVA,IREJN)-1.)*XLU+1.)
     *  *((TWCF(IP,IGT,IVA,IREJN)-1.)*TWU+1.)
     *  *HH
   40 CONTINUE
C
   50 CONTINUE
C
      RETURN
      END
