      SUBROUTINE EVNOIMCALC
C
C  This subroutine computes the no I/M evap pass/fail rates,
C  without and with OBD.
C
C  Called by EVAPCALC.
C
C  Calls None.
C
C  Changes: (Last change first)
C
C  02 May 02 ARI Task 32: Modified for OBD EVAP RESPONSE
C  15 Sep 00 AIR Task 08: New routine.
C
C 
C  Input on call:
C
C    argument list : None
C    common blocks:
C     /LDGOBD/  WMIL, OBDVA, OBDVF, OVDVLIN
C
C  Output on return:
C
C    argument list : None
C    common blocks:
C     /EVIM/ EVNOIM
C
C  Local array subscripts:
C
C  FAIL(4,26,2) - FAIL (IPF, JDX, 1=without OBD, 2=with OBD)
C
C  Local variable / array dictionary:]
C
C   Name   Type                      Description
C  ------  ----  ----------------------------------------------------------
C                for 1=durability only, 2=durabilitY without I/M
C  A        R    1st term coefficients for failure rates
C  AGE      R    Vehicle Age with fractional part
C  B        R    2nd term coefficients for failure rates
C  C        R    3rd term coefficients for failure rates
C  D        R    4th term coefficients for failure rates
C  DAGE     R    Durability age (=age/2.)
C  FAIL     R    Failure rates for 1=fail pressure, 2=fail either, 3=fail both, 4=leaker,
C  IAGE     I    Vehicle Age, integer portion only
C  INTR     R    Intercept for OBD interpolation
C  IPF      I    Index for 1=fail pressure, 2=fail either, 3=fail both, 4=leaker
C  JDX      I    MOBILE6 JDX
C  KDX      I    Index for evap fail age
C  OBD      R    Interpolated OBD
C  SLOPE    R    Slope for OBD interpolation 
C
      IMPLICIT NONE
C
      INCLUDE 'EVIM.I'
      INCLUDE 'LDGOBD.I'
C
C  Local variables
C
      INTEGER IAGE
      INTEGER IPF
      INTEGER JDX
      INTEGER KDX
      REAL    AGE
      REAL    DAGE
      REAL    FAIL(4,26,2)
      REAL    INTR
      REAL    OBD
      REAL    SLOPE
C
C  Define the fail coefficients
C                             FAIL     FAIL      FAIL      
C                             PRES   EITHER      BOTH        LEAKER
      REAL, 
     &DIMENSION(4),
     &SAVE         :: A=(/ 0.60450, 0.72000,  0.3536 ,  0.089018051/)
      REAL,
     &DIMENSION(4),
     &SAVE         :: B=(/17.73300,13.40000,414.96   ,414.61285    /)
      REAL,
     &DIMENSION(4),
     &SAVE         :: C=(/-0.01362,-0.01450, -0.32955, -0.36836226 /)
      REAL,
     &DIMENSION(4),
     &SAVE         :: D=(/ 2.0    , 2.0    ,  1.0    ,  1.0        /)
C
C  Compute the no I/M rates without OBD
C
C  Loop for all JDX
C
      DO JDX=1,26
C
C  Assign age
C
        AGE=JDX-1
C
C  Determine the basic rates for fail pressure, fail either, fail both, leaker
C
        DO IPF=1,4
          FAIL(IPF,JDX,1)=A(IPF)/(1.0+B(IPF)*EXP(C(IPF)*(AGE**D(IPF))))
        END DO
C
C  Compute the pass/pass, pass pressure, pass purge, fail/fail rates
C  without durability or OBD
C
        EVNOIM(1,JDX,1)=1.0-FAIL(2,JDX,1)
        EVNOIM(2,JDX,1)=FAIL(2,JDX,1)-FAIL(1,JDX,1)
        EVNOIM(3,JDX,1)=FAIL(1,JDX,1)-FAIL(3,JDX,1)
        EVNOIM(4,JDX,1)=FAIL(3,JDX,1)
C
      END DO
C
C  Compute the no I/M rates with OBD
C
C  Determine fail pressure, fail either, fail both, leaker rates
C  Account for durability only, then durability without I/M
C
      DO JDX=1,26
C
        IAGE=JDX-1
        DAGE=FLOAT(IAGE)/2.
C
        DO IPF=1,4
C
          FAIL(IPF,JDX,1)=A(IPF)/(1.0+B(IPF)*EXP(C(IPF)*(DAGE**D(IPF))))
C
          IF(IAGE.LE.OBDVA(1)) THEN
            FAIL(IPF,JDX,2)=FAIL(IPF,JDX,1)*(1.0-WMIL*OBDVF(1))
C
          ELSE IF(IAGE.GT.OBDVA(1).AND.IAGE.LE.OBDVA(2)) THEN
            KDX=OBDVA(1)+1
            OBD=OBDVF(2)
            IF(OBDVLIN) THEN  ! Apply linear interpolation (Y=AX+B)
              SLOPE=(OBDVF(2)-OBDVF(1))/FLOAT(OBDVA(2)-OBDVA(1))
              INTR=(FLOAT(OBDVA(2))*OBDVF(1)-OBDVF(2)*FLOAT(OBDVA(1)))
     &             /FLOAT(OBDVA(2)-OBDVA(1))
              OBD=SLOPE*IAGE+INTR
             ENDIF          
            FAIL(IPF,JDX,2)=FAIL(IPF,JDX,1)*(1.0-WMIL*OBD)
     *                     -FAIL(IPF,KDX,1)*(1.0-WMIL*OBD)
     *                     +FAIL(IPF,KDX,2)
C
          ELSE
            KDX=OBDVA(2)+1
            FAIL(IPF,JDX,2)=FAIL(IPF,JDX,1)*(1.0-WMIL*OBDVF(3))
     *                     -FAIL(IPF,KDX,1)*(1.0-WMIL*OBDVF(3))
     *                     +FAIL(IPF,KDX,2)
          END IF
        END DO
C
C  Compute the pass/pass, pass pressure, pass purge, fail/fail rates
C  with OBD, durability, but no I/M
C
        EVNOIM(1,JDX,2)=1.0-FAIL(2,JDX,2)
        EVNOIM(2,JDX,2)=FAIL(2,JDX,2)-FAIL(1,JDX,2)
        EVNOIM(3,JDX,2)=FAIL(1,JDX,2)-FAIL(3,JDX,2)
        EVNOIM(4,JDX,2)=FAIL(3,JDX,2)
C
      END DO
C
      RETURN
      END
