C Version "@(#)$Id$ $Source$ $Date$ 

      SUBROUTINE REGMOD(ICY,INERR)
C
C  REGMOD adjusts model year mileage/registration arrays.
C
C  Called by MOBILE.
C
C  Calls CEVMPD, CVTPDF, HDDMYM, GSFIDX and TFCALX.
C
C  Input on call:
C
C    parameter list: ICY,INERR
C    common blocks:
C    /FLAGS2/ MYMRFG
C    /MAXIMA/ MAXVEH,MAXYRS
C    /MYRSAV/ AMAR,JULMYR
C    /REGION/ IREJN
C
C  Output on return:
C
C    parameter list: INERR
C    common blocks:
C    /MYRCAL/ XMYM,JANMYR,TFMYM
C    /REGISF/ GSFVCT
C
C  Local array subscripts:
C
C  CUM(25)  -  CUM ( JDX )
C  DAF(8)   -  DAF ( IV )
C
C  Local variable / array dictionary:
C
C   Name   Type              Description
C  ------  ----  -------------------------------------------------------
C  CUM      R    cumulative mileage distribution as of January 1
C  DAF      R    unnormalized sum by IV of registrations across my's
C                (acronym is from ef derivation document)
C  FRC      R    fraction of January 1 registrations in current year JDX
C                which are in their JDX'th year of operation
C  FRN      R    fraction of January 1 registrations in current year JDX
C                which are in their JDX-1'th year of operation (1 year newer)
C  AMC      R    mileage accrual rate for vehicle of type IV in its JDX'th
C                year of operation
C  AMN      R    mileage accrual rate for vehicle of type IV in its JDX-1th
C                year of operation (a vehicle in JDX-1 is 1 year "newer" than
C                a JDX year old vehicle)
C  MY       I    see parameter dictionary
C  SUMLDV   R    LDV DAF's combined (both gas & diesel registrations)
C  SUMLDT   R    LDT DAF's combined (both gas & diesel registrations)
C
C  Notes:
C
C  None.
C
C
      INCLUDE 'FLAGS2.I'
      INCLUDE 'FLAGS3.I'
      INCLUDE 'GSFCOM.I'
      INCLUDE 'MAXIMA.I'
      INCLUDE 'MYRCAL.I'
      INCLUDE 'MYRSAV.I'
      INCLUDE 'REGION.I'
      INCLUDE 'REGISF.I'
C
      DIMENSION CUM(25),DAF(8)
      REAL AMC,AMN
C
C  Compute January 1 registration distributions, JANMYR.
C
      DO 20 IV=1,MAXVEH
C
C  Step 1: for each model year of each vehicle type, select and apply the
C          corresponding GSF to the (gas & diesel combined) July 1
C          registrations.
C
      DO 10 JDX=1,MAXYRS
      MY=ICY-JDX+1
      JANMYR(JDX,IV)=GSFIDX(ICY,MY,IV,IREJN)*JULMYR(JDX,IV)
   10 CONTINUE
C
C  Step 2: backup ICY's registrations from 7/1 to 1/1 -
C          LDV/T => my starts 10/1 prior yr => 1/1 sales are 1/3 of 7/1 sales.
C          HDV & MC => my starts 1/1 => 1/1 sales are 0.
C
      JANMYR(1,IV)=JANMYR(1,IV)/3.
      IF(IV.EQ.4.OR.IV.GE.7) JANMYR(1,IV)=0.0
   20 CONTINUE
C
C  Step 3: normalize JANMYR.
C
      DO 40 IV=1,MAXVEH
      DAF(IV)=0.0
      DO 30 JDX=1,MAXYRS
      DAF(IV)=DAF(IV)+JANMYR(JDX,IV)
   30 CONTINUE
      DO 40 JDX=1,MAXYRS
      IF(DAF(IV).GT.0.0) JANMYR(JDX,IV)=JANMYR(JDX,IV)/DAF(IV)
   40 CONTINUE
C
C  Given DAF (JAN 1 unnormalized total (MAXYRS my's combined) registrations),
C  calculate the <gas | diesel> / (gas + diesel) fractions for LDV & LDT.
C  HDV & MC are assigned unity, since there is no combined gas + diesel
C  number in VCOUNT in TFCALX to split by fuel type fraction for these
C  vehicle types.  An attempt was made to model the HDV gas vs diesel sales,
C  but the results proved unsatisfactory.
C
      SUMLDV=DAF(1)+DAF(5)
      SUMLDT=DAF(2)+DAF(6)
C
      GSFVCT(1)=DAF(1)/SUMLDV
      GSFVCT(2)=DAF(2)/SUMLDT
      GSFVCT(3)=1.
      GSFVCT(4)=1.
      GSFVCT(5)=DAF(5)/SUMLDV
      GSFVCT(6)=DAF(6)/SUMLDT
      GSFVCT(7)=1.
      GSFVCT(8)=1.
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
      DO 60 IV=1,MAXVEH
      FRC=.25
      IF(IV.EQ.4.OR.IV.EQ.7.OR.IV.EQ.8) FRC=0.0
      FRN=1.-FRC
      AMN=AMAR(1,IV)
      AMC=AMAR(2,IV)
      CUM(1)=.5*FRC*AMN
      CUM(2)=FRC*(.5*FRC*AMC+AMN)+FRN*((1.+FRC)*.5*AMN)
      AMC=AMC+AMN
C
      DO 50 JDX=3,MAXYRS
      CUM(JDX)=FRC*(.5*FRC*AMAR(JDX,IV)+AMC)+FRN*(.5*(1.+FRC)
     *   *AMAR(JDX-1,IV)+AMN)
      AMN=AMC
      AMC=AMC+AMAR(JDX,IV)
   50 CONTINUE
C
C  Subroutine GETCUM produces the cumulative mileage partial sums array,
C  CUMMIL, from XMYM.
C
      XMYM(1,IV)=CUM(1)
      DO 60 JDX=2,MAXYRS
      XMYM(JDX,IV)=CUM(JDX)-CUM(JDX-1)
   60 CONTINUE
C
C  Now figure the January 1 model year mileage accumulation rates, TFMYM.
C
      DO 80 IV=1,MAXVEH
      FRC=.25
      TFMYM(1,IV)=AMAR(1,IV)
      IF(IV.NE.4.AND.IV.NE.7.AND.IV.NE.8) GOTO 70
      FRC=0.0
      TFMYM(1,IV)=0.0
   70 FRN=1.-FRC
C
      DO 80 JDX=2,MAXYRS
      TFMYM(JDX,IV)=AMAR(JDX-1,IV)*FRN+AMAR(JDX,IV)*FRC
   80 CONTINUE
C
C  Use TFMYM to generate miles per day by model year group.
C
      CALL CEVMPD
C
C  Use JANMYR to generate trips per day fractions by model year gorup.
C
      CALL CVTPDF
C
C  Use JANMYR, TFMYM, and GSFVCT to generate the vmt mix and travel fractions.
C
      CALL TFCALX(ICY,INERR)
C
      RETURN
      END
