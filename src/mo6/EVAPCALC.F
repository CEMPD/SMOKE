      SUBROUTINE EVAPCALC(ICY,IV) 
C 
C  This routine drives the calculation for the evap weighting factor arrays for  
C  without and with i/m and without and with obd. 
C 
C  Called by HCCALC. 
C 
C  Calls EVNOIMCALC, EVFAILCALC, FINDIM, EVIMCALC, SAWTOOTH. 
C 
C  Changes: (Last change first) 
C 
C  13 Dec 00 @epa-bag t-x6 added ACTUAL_MY and clean air act switch parameters 
C  03 Oct 00 AIR Task 08: Added constant IEXHAUST=1 and IEVAP=2 for call to FINDIM 
C  26 Sep 00 AIR Task 08: Corrected block data. 
C  15 Sep 00 AIR Task 08: New routine. 
C 
C  
C  Input on call: 
C 
C    argument list : ICY, IV 
C    common blocks: None 
C 
C  Output on return: 
C 
C    argument list : None 
C    common blocks: 
C    /EVIM/ EVIM 
C 
C  Local variable / array dictionary:] 
C 
C   Name     Type                      Description 
C  ------    ----  ---------------------------------------------------------- 
C  IV         I    Vehicle class 
C  JDX        I    MOBILE6 JDX 
C  IEVAP      I    Switch to obtain evap I/M programs from FINDIM 
C  IEXHAUST   I    Switch to obtain exhaust I/M programs from FINDIM 
C  IMPGM      I    Applicable Exhaust I/M program number 
C  EVIMPGM    I    Applicable Evap I/M program number 
C  IOBD       I    1=no OBD, 2=with OBD 
C  IPF        I    Index: 1=pressure/pass purge, 2=pass pressure/fail purge, 
C                  3=fail pressure/pass purge, 4=fail pressure/fair purge 
C  ICY        I    Calendar year 
C  MY         I    Model year 
C 
      IMPLICIT NONE 
C 
      INCLUDE 'EVIM.I' 
      INCLUDE 'FLAGS3.I' 
C 
      INTEGER EVIMPGM 
      INTEGER FINDIM 
      INTEGER ICY 
      INTEGER IMPGM 
      INTEGER IOBD 
      INTEGER IPF 
      INTEGER IV 
      INTEGER JDX 
      INTEGER MY 
      INTEGER ACTUAL_MY 
      INTEGER, SAVE :: IEXHAUST = 1 
      INTEGER, SAVE :: IEVAP    = 2 
C 
C  Compute the no I/M rates, with and without OBD, for all 26 ages. 
C 
      CALL EVNOIMCALC 
C 
C  Compute the failure rates, for all 26 ages. 
C 
      CALL EVFAILCALC(IV) 
C 
C  Loop for each age 
C 
      DO JDX=1,26 
C 
C  Determine the model year 
C 
        ACTUAL_MY=ICY-JDX+1 
        MY=ACTUAL_MY 
C 
        IF(CAAFLG.EQ.1 .AND. ICY.GE.CAAMY) THEN 
          MY = MIN(CAAMY,ACTUAL_MY) 
        ENDIF 
C 
C  Determine if there is an EVAP I/M program available for this vehicle 
C  class, calendar year and model year. 
C 
        EVIMPGM=FINDIM(ICY,MY,IV,IEVAP) 
C 
C  If no I/M program, then assign no I/M rates, without and with OBD, to 
C  the with-I/M arrays. Otherwise, determine if an exhaust I/M program 
C  is active and compute the with I/M evap rates. 
C 
        IF(EVIMPGM.EQ.0) THEN 
C 
          DO IOBD=1,2 
            DO IPF=1,4 
              EVIM(IPF,JDX,IOBD)=EVNOIM(IPF,JDX,IOBD) 
            END DO 
          END DO 
C 
        ELSE 
C 
          IMPGM=FINDIM(ICY,MY,IV,IEXHAUST) 
C 
          CALL EVIMCALC(MY,IMPGM,EVIMPGM,JDX) 
C 
        ENDIF 
C 
      END DO 
C 
C  Apply the sawtooth weighting method for without and with OBD. 
C 
      CALL SAWTOOTH(IV,ICY) 
C 
      RETURN 
      END 
       
