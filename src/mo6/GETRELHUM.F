      SUBROUTINE GETRELHUM(INPREC,INERR,RC)
C
C  GETRELHUM reads in the hourly relative humidity values, in percent.
C
C  Called by DSPTCH.
C
C  Calls QUITER.
C
C  Changes : 
C
C  02 Jul 02 AIR Task 36, new subroutine.
C       
C  Input on call:
C
C    argument list: INPREC, INERR
C
C
C  Output on return:
C
C    argument list: INERR, RC
C
C    common block:
C    /ALUHIN/   RHMAX, RHMIN, RH_HRLY, RH_DIST
C
C  Local variable / array dictionary:
C
C  Name                 Type              Description
C  ------               ----  ----------------------------------------------------
C  IH                    I   Hour index
C
C
      IMPLICIT NONE
      INCLUDE 'ALUHIN.I'
      INCLUDE 'IOUCOM.I'
      INCLUDE 'MAXIMA.I'
C
      INTEGER, INTENT(INOUT) :: INERR
      INTEGER, INTENT(OUT) :: RC
      CHARACTER INPREC*(*)
C
      INTEGER IH
C
      RC = 0
C
      RH_DIST = .FALSE.
C
C  Read and verify the hourly relative humidities.
C
        READ(INPREC(21:),*,ERR=95,END=95) (RH_HRLY(IH), IH = 1,12)
        READ(IOUGEN,*,ERR=95,END=95) (RH_HRLY(IH), IH = 13,24)
C
C Find the maximum and minimum relative humidities. Also, check relative humidity
C ranges. These bounds are for catching data entry errors only.
C
        RHMIN = RH_HRLY(1)
        RHMAX = RH_HRLY(1)
C
        DO IH = 1,MAXIH
C
         IF(RH_HRLY(IH).LT.0.0.OR.RH_HRLY(IH).GT.100.0) THEN
           CALL QUITER(RH_HRLY(IH),IH,89,INERR)
           RC=1
           GOTO 99
         ENDIF
C
          IF( RH_HRLY(IH) .GT. RHMAX ) RHMAX = RH_HRLY(IH)
          IF( RH_HRLY(IH) .LT. RHMIN ) RHMIN = RH_HRLY(IH)
        END DO
C
      RH_DIST = .TRUE.
C
      GOTO 99
C
C  Error or EOF on any attempted read.
C
   95 RC = 1
      CALL QUITER(0.,0,102,INERR)
C
   99 RETURN
      END
