C Version "@(#)$Id$ $Source$ $Date$ 

      SUBROUTINE CALCUM(ICY)
C
C  CALCUM adjusts model year mileage arrays.
C
C  Called by ONESEC.
C
C  Calls HDDMYM.
C
C  Input on call:
C
C    parameter list: ICY,INERR
C    common blocks:
C    /FLAGS2/ MYMRFG
C    /MAXIMA/ MAXVEH,MAXYRS
C    /MYRSAV/ AMAR
C
C  Output on return:
C
C    parameter list: INERR
C    common blocks:
C    /MYRCAL/ XMYM
C
C  Local array subscripts:
C
C  CUM(25)  -  CUM ( JDX )
C
C  Local variable / array dictionary:
C
C   Name   Type              Description
C  ------  ----  -------------------------------------------------------
C  AMC      R    mileage accrual rate for vehicle of type IV in its JDX'th
C                year of operation
C  AMN      R    mileage accrual rate for vehicle of type IV in its JDX-1th
C                year of operation (a vehicle in JDX-1 is 1 year "newer" than
C                a JDX year old vehicle)
C  CUM      R    cumulative mileage distribution as of January 1
C  FRC      R    fraction of January 1 registrations in current year JDX
C                which are in their JDX'th year of operation
C  FRN      R    fraction of January 1 registrations in current year JDX
C                which are in their JDX-1'th year of operation (1 year newer)
C
C  Notes:
C
C  CALCUM was created for MOBILE4.1v6.
C
C
      INCLUDE 'FLAGS2.I'
      INCLUDE 'MAXIMA.I'
      INCLUDE 'MYRCAL.I'
      INCLUDE 'MYRSAV.I'
C
      DIMENSION CUM(25)
      REAL AMC,AMN
C
C  Generate the mileage accrual rates.
C
C  If MOBILE4 defaults used, determine the HDDV calendar year mileage
C  by registration weighting. Otherwise, use user input mileage.
C
      IF(MYMRFG.EQ.1.OR.MYMRFG.EQ.3) CALL HDDMYM(ICY)
C
C  Compute the annual mileage / 100,000  distributions, XMYM.
C
      DO 30 IV=1,MAXVEH
      FRC=.25
      IF(IV.EQ.4.OR.IV.EQ.7.OR.IV.EQ.8) FRC=0.0
      FRN=1.-FRC
      AMN=AMAR(1,IV)
      AMC=AMAR(2,IV)
      CUM(1)=.5*FRC*AMN
      CUM(2)=FRC*(.5*FRC*AMC+AMN)+FRN*((1.+FRC)*.5*AMN)
      AMC=AMC+AMN
C
      DO 10 JDX=3,MAXYRS
      CUM(JDX)=FRC*(.5*FRC*AMAR(JDX,IV)+AMC)+FRN*(.5*(1.+FRC)
     *   *AMAR(JDX-1,IV)+AMN)
      AMN=AMC
      AMC=AMC+AMAR(JDX,IV)
   10 CONTINUE
C
C  Subroutine GETCUM produces the cumulative mileage partial sums array,
C  CUMMIL, from XMYM.
C
      XMYM(1,IV)=CUM(1)
      DO 20 JDX=2,MAXYRS
      XMYM(JDX,IV)=CUM(JDX)-CUM(JDX-1)
   20 CONTINUE
   30 CONTINUE
C
      RETURN
      END
