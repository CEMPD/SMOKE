      SUBROUTINE SAWTOOTH(IV,ICY) 
C 
C  This subroutine computes the sawthooth weighted I/M evap pass/fail rates, 
C  without and with OBD. 
C 
C  Called by EVAPCALC. 
C 
C  Calls FINDIM. 
C 
C  Changes: (Last change first) 
C 
C  13 Dec 00 @epa-bag t-x6 added ACTUAL_MY and clean air act switch parameters 
C  26 Sep 00 AIR Task 08: Reinstalled SAWFLAG processing. 
C  15 Sep 00 AIR Task 08: New routine. 
C 
C  
C  Input on call: 
C 
C    argument list : IV,ICY 
C    common blocks: 
C     /LDGOBD/ MIL, OBDF 
C     /IMPAR1/ IMGRACE 
C     /EVIM/   EVNOIM, EVIM 
C     /IMPAR6/ IFREQ 
C     /EVAL/   MEVAL 
C 
C  Output on return: 
C 
C    argument list : None 
C    common blocks: 
C     /EVIM/ EVSAW 
C 
C  Local variable / array dictionary:] 
C 
C   Name   Type                      Description 
C  ------  ----  ---------------------------------------------------------- 
C  ADJIM    R    I/M Adjustment (see M6IMFIXC.XLS spreadsheet for element origins.) 
C  AGE      R    Vehicle Age with fractional part 
C  AIM1     I    Average I/M age, earlier 
C  AIM2     I    Average I/M age  later 
C  DELTAT1  R    Fractional age, earlier 
C  DELTAT2  R    Fractional age, later 
C  EVIMPGM  I    Applicable evap I/M program number 
C  FFRC     R    Fraction of model year after MEVAL 
C  FFRN     R    Fraction of model year before MEVEL 
C  FRACT1   R    Fraction only part of DELTAT1   
C  FRACT2   R    Fraction only part of DELTAT1   
C  IAGE     I    Vehicle Age, integer portion only 
C  IOBD     I    Index, 1=no OBD, 2=with OBD 
C  IPF      I    Index for 1=fail pressure, 2=fail either, 3=fail both, 4=leaker 
C  ICY      I    Calendar year 
C  JDX      I    MOBILE6 JDX 
C  MY       I    Model Year 
C  PNEX     R    Next sawtooth portion 
C  PREV     R    Previous sawtooth portion 
C  X        R    Temporary variable to hold without or with sawtooth weighting 
C 
      IMPLICIT NONE 
C 
      INCLUDE 'EVIM.I' 
      INCLUDE 'FLAGS3.I' 
      INCLUDE 'IMPAR1.I' 
      INCLUDE 'EVAL.I' 
      INCLUDE 'IMPAR6.I' 
      INCLUDE 'EVAPACT.I' 
C 
      INTEGER AIM1 
      INTEGER AIM2 
      INTEGER EVIMPGM 
      INTEGER FINDIM 
      INTEGER ICY 
      INTEGER IOBD 
      INTEGER IPF 
      INTEGER IV 
      INTEGER JDX 
      INTEGER MY 
      INTEGER ACTUAL_MY 
      REAL    DELTAT1 
      REAL    DELTAT2 
      REAL    FFRC 
      REAL    FFRN 
      REAL    FRACT1 
      REAL    FRACT2 
      REAL    PNEX 
      REAL    PREV 
      REAL    X 
C 
C  Determine FRC and FRN from the evaulation month 
C  (Do NOT use the FRC and FRN values in EVAL.I) 
C 
      FFRC=(FLOAT(MEVAL)-10.)/12. 
      IF(FFRC.LT.0.) FFRC=1.+FFRC 
C       
      FFRN=1.-FFRC 
C 
C  Loop for all ages 
C 
      DO JDX=1,25 
C 
C  Determine if there is an active evap I/M program for this model year  
C 
        ACTUAL_MY=ICY-JDX+1 
        MY=ACTUAL_MY 
C 
        IF(CAAFLG.EQ.1 .AND. ICY.GE.CAAMY) THEN 
          MY = MIN(CAAMY,ACTUAL_MY) 
        ENDIF 
C 
        EVIMPGM=FINDIM(ICY,MY,IV,2) 
C 
C  If an evap I/M program is active, set up the with-I/M values 
C 
        IF(EVIMPGM.GT.0) THEN 
C 
C Compute the average I/M ages 
C 
          AIM1=0 
          IF(JDX.GE.2.AND.JDX.GT.IMGRACE(EVIMPGM)+1) 
     *      AIM1=JDX-2-MOD(JDX-2-IMGRACE(EVIMPGM),IFREQ(EVIMPGM)) 
C 
          AIM2=0 
          IF(JDX.GE.IMGRACE(EVIMPGM)+1) 
     *      AIM2=JDX-1-MOD(JDX-1-IMGRACE(EVIMPGM),IFREQ(EVIMPGM)) 
C 
C  Compute the DELTAT's and their fractions 
C 
          IF(JDX.LT.2) THEN 
            DELTAT1=0. 
          ELSE 
            DELTAT1=1.0-FFRN/2.0+JDX-2.0-AIM1 
          ENDIF 
          DELTAT2=      FFRC/2.0+JDX-1.0-AIM2 
C 
          FRACT1=DELTAT1-INT(DELTAT1) 
          FRACT2=DELTAT2-INT(DELTAT2) 
C 
C  Map AIM1 and AIM2 into JDX space 
C 
          AIM1=AIM1+1 
          AIM2=AIM2+1 
C 
        ENDIF 
C 
C  Loop for without and with OBD 
C 
        DO IOBD=1,2 
C 
C  Loop for each failure type 
C 
          DO IPF=1,4 
C 
C  Compute no-I/M or I/M values depending upon I/M activity 
C 
            IF(EVIMPGM.EQ.0) THEN 
C 
C  **** COMPUTE NO IM VALUES***** 
C 
C   If user has disabled sawtooth averaging, use JDX value, 
C   otherwise, perform sawtooth calculations. 
C 
              IF(SAWFLAG.EQ.0) THEN 
C 
                X=EVNOIM(IPF,JDX,IOBD) 
C 
              ELSE 
C 
C   Remember to treat the first row differently (age zero) 
C 
                IF(JDX.EQ.1) THEN 
C 
                  X = 0.5*FFRC *EVNOIM(IPF,JDX+1,IOBD) 
     *                + (1.-0.5*FFRC)*EVNOIM(IPF,JDX,IOBD) 
C 
C       Treat later JDX-ages using the full  
C       no IM interpolation formula 
C         
                ELSE 
C 
C          Interpolate for the FRC segment 
C 
                  PNEX = 0.5*FFRC *EVNOIM(IPF,JDX+1,IOBD) 
     *                   + (1.-0.5*FFRC)*EVNOIM(IPF,JDX,IOBD) 
C 
C          and for the FRN segment 
C 
                  PREV = 0.5*FFRN *EVNOIM(IPF,JDX-1,IOBD) 
     *                   + (1.-0.5*FFRN)*EVNOIM(IPF,JDX,IOBD) 
C 
C          and combine for the sawtooth average 
C  
                  X=FFRC*PNEX + FFRN*PREV 
C  
                END IF   ! JDX 
              END IF   ! SAWFLAG 
C 
              EVSAW(IPF,JDX,IOBD)=X 
C 
C  ***** COMPUTE WITH IM VAULUES ***** 
C 
            ELSE 
C 
C  Handle first year 
C 
              IF (JDX.EQ.1) THEN 
C 
                X = 0.5*FFRC *EVNOIM(IPF,JDX+1,IOBD) 
     *              + (1.-0.5*FFRC)*EVNOIM(IPF,JDX,IOBD) 
C 
              ELSE 
C 
C          Interpolate for the FRC segment 
C 
                PNEX = EVIM(IPF,AIM2,IOBD) 
     *                 - EVNOIM(IPF,AIM2,IOBD) 
     *                 + (1.-FRACT2) * EVNOIM(IPF,JDX,IOBD) 
     *                 + FRACT2 * EVNOIM(IPF,JDX+1,IOBD) 
C 
C          and for the FRN segment 
C 
                PREV = EVIM(IPF,AIM1,IOBD) 
     *                 - EVNOIM(IPF,AIM1,IOBD) 
     *                 + (1.-FRACT1) * EVNOIM(IPF,JDX-1,IOBD) 
     *                 + FRACT1 * EVNOIM(IPF,JDX,IOBD) 
C 
C          and combine for the sawtooth average 
C  
                X=FFRC*PNEX + FFRN*PREV 
C  
              END IF !JDX 
C 
              EVSAW(IPF,JDX,IOBD)=X 
C 
            ENDIF  !EVIMPGM 
C 
          END DO  ! IPF 
        END DO  ! IOBD 
      END DO  ! JDX 
C 
      RETURN 
      END 
