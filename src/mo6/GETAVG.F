      SUBROUTINE GETAVG(INPREC,INERR,RC) 
C 
C  GETAVG obtains the user specified average speed for a specific 
C  roadway type and the roadway type distribution. 
C 
C 
C  Called by DSPTCH. 
C 
C 
C  Calls NXTTOK, QUITER, UCCOMP. 
C 
C 
C  Changes: (Last change first) 
C 
C  22 Jan 02 @EPA-djb  Bug379   Added check for Case #2 (freeway) and 
C                               changed case error to device IOUREP. 
C  18 Jan 02 @EPA-djb  Bug381   Initialize USRVMT after error 580, 
C                               changed writes to IOUREP to QUITER calls. 
C  14 Jan 02 @EPA-bag  Bug375   Added calls to QUITER to avoid user input of 
C                            VMT fractions which are inconsistent with facility 
C                            options, FREEWAY and AREAWIDE 
C  19 Jul 01 @EPA-djb Wish235, Add Average Speed command. 
C 
C  Input on call: 
C 
C    argument list: INPREC,NSCEN,INERR,RC 
C    common blocks: 
C      /SPEED0/ DEFAVG 
C 
C  Output on return: 
C 
C    argument list: INPREC, NSCEN, INERR, RC 
C    common blocks: 
C      /SPEED0/ USRSPD,USRVMT,USRFAC 
C 
C  Local variable / array dictionary: 
C 
C   Name   Type              Description 
C  ------  ----  ------------------------------------------------------- 
C  INPSPD    R   Variable to temporarily store the user input average speed. 
C  MAXSPD    R   Maximum speed that can be calculated from the combination 
C                of average speed and roadway type VMT fractions. 
C  MINSPD    R   Minimum speed that can be calculated from the combination 
C                of average speed and roadway type VMT fractions. 
C  VMTSUM    R   Sum of user supplied roadway type VMT fractions. 
C 
C 
C  Notes: 
C 
C 
      IMPLICIT NONE 
C 
      INCLUDE 'IOUCOM.I'      ! IOUOUT,IOUREP 
      INCLUDE 'SPEED0.I'      ! USRSPD,USRVMT,USRFAC,DEFAVG,AVGFLG 
C 
C     Declare variables. 
C 
      CHARACTER(80)   :: DATAFIELD 
      CHARACTER(*)    :: INPREC 
      CHARACTER*10    :: TESTFIELD 
C 
      INTEGER, INTENT(INOUT)  :: INERR 
      INTEGER, INTENT(INOUT)  :: RC 
C 
      INTEGER         :: IROAD    ! Roadway type index. 
      INTEGER         :: PTR1 
      INTEGER         :: PTR2 
C 
      REAL            :: INPSPD   ! Storage of user input average speed. 
      REAL            :: MAXSPD   ! Given user inputs, maximum allowed area wide average speed. 
      REAL            :: MINSPD   ! Given user inputs, maximum allowed area wide average speed. 
      REAL            :: VMTSUM   ! Sum of user supplied roadway type VMT fractions. 
      REAL            :: VMT1     ! Freeway portion of freeway and ramp VMT. 
      REAL            :: VMT2     ! Ramp portion of freeway and ramp VMT. 
C 
      PTR2 = 21 
C 
C     Read the user supplied average speed value. 
C 
      CALL NXTTOK(INPREC,DATAFIELD,PTR1,PTR2) 
      READ (DATAFIELD,*,IOSTAT=RC) USRSPD 
      IF(USRSPD.LT.2.5) THEN 
        CALL QUITER(USRSPD,0,52,INERR) 
        USRSPD=2.5 
      ELSEIF(USRSPD.GT.65.0) THEN 
        CALL QUITER(USRSPD,0,96,INERR) 
        USRSPD=65.0 
      ENDIF 
C 
C     Read the user specified average speed option. 
C 
      CALL NXTTOK(INPREC,DATAFIELD,PTR1,PTR2) 
      CALL UCCOMP(DATAFIELD,TESTFIELD) 
      IF (TESTFIELD.EQ.'NON-RAMP') THEN 
        USRFAC=1 
      ELSEIF (TESTFIELD.EQ.'FREEWAY') THEN 
        USRFAC=2 
      ELSEIF (TESTFIELD.EQ.'ARTERIAL') THEN 
        USRFAC=3 
      ELSEIF (TESTFIELD.EQ.'AREAWIDE') THEN 
        USRFAC=4 
      ELSEIF (TESTFIELD.EQ.'LOCAL') THEN 
        USRFAC=5 
      ELSEIF (TESTFIELD.EQ.'RAMP') THEN 
        USRFAC=6 
      ELSE 
        CALL QUITER(0.0,0,512,INERR) 
        RC = 1 
        GOTO 99 
      ENDIF 
C 
C     Read the user supplied distribution of roadway VMT. 
C 
      DO IROAD=1,4 
        CALL NXTTOK(INPREC,DATAFIELD,PTR1,PTR2) 
        READ (DATAFIELD,*,IOSTAT=RC) USRVMT(IROAD) 
      END DO 
      VMTSUM=USRVMT(1)+USRVMT(2)+USRVMT(3)+USRVMT(4) 
      IF(VMTSUM.LT.0.0001 .OR. RC.NE.0) THEN   ! Set values to defaults. 
        USRVMT(1) = 0.342 
        USRVMT(2) = 0.498 
        USRVMT(3) = 0.130 
        USRVMT(4) = 0.030 
        RC=0 
      ELSEIF(     (VMTSUM.GT.0.999 .AND. VMTSUM.LT.1.001) 
     *       .OR. (VMTSUM.GT. 99.9 .AND. VMTSUM.LT.100.1) ) THEN 
        IF(VMTSUM.NE.1.0) THEN  ! Normalize user input. 
          USRVMT(1) = USRVMT(1)/VMTSUM 
          USRVMT(2) = USRVMT(2)/VMTSUM 
          USRVMT(3) = USRVMT(3)/VMTSUM 
          USRVMT(4) = USRVMT(4)/VMTSUM 
        ENDIF 
      ELSE 
        CALL QUITER(VMTSUM,0,16,INERR) 
        USRVMT(1) = USRVMT(1)/VMTSUM 
        USRVMT(2) = USRVMT(2)/VMTSUM 
        USRVMT(3) = USRVMT(3)/VMTSUM 
        USRVMT(4) = USRVMT(4)/VMTSUM 
      ENDIF 
C 
C     Check for inconsistent freeway and arterial vehicle fractions 
      IF(USRFAC==4.AND.USRVMT(1)<=0.00001.AND.USRVMT(2)<=0.00001) THEN 
        CALL QUITER(0.0,0,577,INERR) 
        RC = 1 
        GOTO 99 
      ENDIF 
C     Check for inconsistent freeway and ramp vehicle fractions 
      IF(USRFAC==2) THEN 
        IF(USRVMT(1)<=0.00001.AND.USRVMT(4)<=0.00001) THEN 
          CALL QUITER(0.0,0,578,INERR) 
          RC = 1 
          GOTO 99 
        ENDIF 
        IF(USRVMT(1)<=0.00001.AND.USRVMT(2)<=0.00001) THEN 
          CALL QUITER(0.0,0,580,INERR) 
          RC = 1 
          USRVMT(1) = 0.342  !Reset VMT to defaults. 
          USRVMT(2) = 0.498 
          USRVMT(3) = 0.130 
          USRVMT(4) = 0.030 
          GOTO 99 
        ENDIF 
      ENDIF 
C 
C  If the freeway option is chosen, check the combination of 
C  user supplied values to determine if this option can be used. 
C 
      IF(USRFAC.EQ.2) THEN 
C     Minimum and maximum combined freeway and freeway ramp speeds 
C       based on the user supplied (or default) VMT fractions for both. 
        VMT1=0.0 
        VMT2=0.0 
        IF(USRVMT(1).GT.0.0 .OR. USRVMT(4).GT.0.0) THEN 
          VMT1 = USRVMT(1) / (USRVMT(1) + USRVMT(4)) 
          VMT2 = USRVMT(4) / (USRVMT(1) + USRVMT(4)) 
        ENDIF 
        MINSPD = 1.0/((1.0/2.5)*VMT1+(1.0/DEFAVG(4))*VMT2) 
        MAXSPD = 1.0/((1.0/65.0)*VMT1+(1.0/DEFAVG(4))*VMT2) 
        IF(USRSPD.LT.MINSPD) THEN 
          CALL QUITER(MINSPD,0,514,INERR) 
          USRSPD=MINSPD 
        ELSEIF(USRSPD.GT.MAXSPD) THEN 
          CALL QUITER(MAXSPD,0,515,INERR) 
          USRSPD=MAXSPD 
        ENDIF 
        INPSPD=USRSPD 
      ENDIF 
C 
C  If the areawide option is chosen, check the combination of 
C  user supplied values to determine if this option can be used. 
C 
      IF(USRFAC.EQ.4) THEN 
C     Minimum and maximum freeway and arterial speeds assume 
C       than the same adjustment will be applied to both. 
        MINSPD = 1.0/((1.0/2.92)*USRVMT(1)+(1.0/2.5)*USRVMT(2) 
     *    +(1.0/DEFAVG(3))*USRVMT(3)+(1.0/DEFAVG(4))*USRVMT(4)) 
        MAXSPD = 1.0/((1.0/65.0)*USRVMT(1)+(1.0/55.56)*USRVMT(2) 
     *    +(1.0/DEFAVG(3))*USRVMT(3)+(1.0/DEFAVG(4))*USRVMT(4)) 
        IF(USRSPD.LT.MINSPD) THEN 
          CALL QUITER(MINSPD,0,14,INERR) 
          USRSPD=MINSPD 
        ELSEIF(USRSPD.GT.MAXSPD) THEN 
          CALL QUITER(MAXSPD,0,15,INERR) 
          USRSPD=MAXSPD 
        ENDIF 
C     Determine and apply adjustment to default arterial/collector average speed. 
        INPSPD=USRSPD 
        USRSPD=DEFAVG(2) * 
     *    ((1/DEFAVG(1))*USRVMT(1)+(1/DEFAVG(2))*USRVMT(2)) / 
     *    ((1/USRSPD)-((1/DEFAVG(3))*USRVMT(3)+(1/DEFAVG(4))*USRVMT(4))) 
      ENDIF 
C 
C  Alert user in descriptive output that an average speed will be used 
C  and that the facility VMT will be changed. 
C 
      SELECT CASE (USRFAC) 
          CASE (1) 
            CALL QUITER(USRSPD,0,581,INERR) 
          CASE (2) 
            CALL QUITER(INPSPD,0,582,INERR) 
          CASE (3) 
            CALL QUITER(USRSPD,0,583,INERR) 
          CASE (5) 
            CALL QUITER(DEFAVG(3),0,585,INERR) 
          CASE (6) 
            CALL QUITER(DEFAVG(4),0,586,INERR) 
          CASE DEFAULT     ! Case 4 is the default. 
            CALL QUITER(INPSPD,0,584,INERR) 
      END SELECT 
C 
C  Error or EOF on any attempted read. 
C 
      IF(RC.NE.0)THEN 
         INERR = INERR + 1 
      END IF 
C 
  99  RETURN 
      END 
