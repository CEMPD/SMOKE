      SUBROUTINE AVGSPD(INERR) 
C 
C  AVGSPD uses the average speed user input parameters to create the 
C  appropriate average speed distributios for each hour of the day that 
C  will match the average speed intended by the user.  This also 
C  will require adjusting the VMT by facility in each hour to match the 
C  intended mix of roadway types. 
C 
C  Called by INISCN 
C 
C  Calls *NONE* 
C 
C  Changes: (Last change first) 
C 
C  22 Jan 02 @EPA-djb  Bug379   Set limit on maximum freeway speed for 
C                               Case #2 and added comments. 
C  14 Jan 02 @EPA-bag  Bug375   Added call to QUITER to avoid VMT fractions 
C                      which are inconsistent with freeway facility option 
C  19 Jul 01 @EPA-djb Wish235, Add Average Speed command. 
C 
C  Input on call: 
C 
C    common blocks: 
C    /SPEED0/ USRSPD,USRVMT,USRFAC,DEFAVG 
C 
C  Output on return: 
C 
C    common blocks: 
C    /SPEED9/ SVMT, FVMT 
C 
C  Notes: 
C 
C  The fourteen speed bins are: 
C 
C    1 : (2.5 mph) 0-2.5 mph       8 : (35 mph) >32.5-37.5 mph 
C    2 : (5 mph)  >2.5-7.5 mph     9 : (40 mph) >37.5-42.5 mph 
C    3 : (10 mph) >7.5-12.5 mph   10 : (45 mph) >42.5-47.5 mph 
C    4 : (15 mph) >12.5-17.5 mph  11 : (50 mph) >47.5-52.5 mph 
C    5 : (20 mph) >17.5-22.5 mph  12 : (55 mph) >52.5-57.5 mph 
C    6 : (25 mph) >22.5-27.5 mph  13 : (60 mph) >57.5-62.5 mph 
C    7 : (30 mph) >27.5-32.5 mph  14 : (65 mph) >62.5 mph 
C 
C  The six AVERAGE SPEED roadway scenarios (USRFAC) are: 
C 
C    1 : NON-RAMP 
C    2 : FREEWAY 
C    3 : ARTERIAL 
C    4 : AREAWIDE 
C    5 : LOCAL 
C    6 : RAMP 
C 
      IMPLICIT NONE  
C 
      INCLUDE 'SPEED0.I'      ! USRSPD,USRVMT,USRFAC,DEFAVG,AVGFLG 
      INCLUDE 'SPEED9.I'      ! SVMT, FVMT 
C 
C     Declare parrameter list 
C 
      INTEGER, INTENT(INOUT)  :: INERR 
C 
C     Declare local variables. 
C 
      INTEGER                :: I         ! Speed bin index. 
      INTEGER                :: J         ! Roadway type (1:freeway, 2:arterial) 
      INTEGER                :: IH        ! Hour index. 
      INTEGER                :: IV        ! Vehicle class index. 
      REAL                   :: BIN1      ! Lower bin fraction. 
      REAL                   :: BIN2      ! Upper bin fraction. 
      REAL                   :: FWYSPD    ! Non-ramp freeway average speed. 
      REAL, DIMENSION(14)    :: NSPD=(/   ! Nominal speed of bin. 
C      Bin   1     2     3     4     5     6     7 
     *      2.5,  5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 
C            8     9    10    11    12    13    14 
     *     35.0, 40.0, 45.0, 50.0, 55.0, 60.0, 65.0/) 
C 
C  Initialize the speed VMT array to zero before calculations. 
C 
      FWYSPD = 0.0 
      IF(USRFAC.LE.4) THEN                      ! If local or ramp, don't zero distributions. 
        DO J=1,2                                ! Loop, 1:Freeway, 2:Arterial 
          IF(USRFAC.LE.2 .AND. J.EQ.2) CYCLE    ! If freeay, don't zero arterial distributions. 
          IF(USRFAC.EQ.3 .AND. J.EQ.1) CYCLE    ! If arterial, don't zero freeway distributions. 
          DO I=1,14                             ! Loop over fourteen average speed bins. 
            DO IH=1,24                          ! Loop over 24 hours. 
              SVMT(J,I,IH) = 0.0                ! Initialize speed VMT to zero. 
            END DO 
          END DO 
        END DO 
      ENDIF 
C 
C  Select the average speed case specified by the user: 
C 
      SELECT CASE (USRFAC) 
C 
C         Case 1 : All VMT occurs on freeways, not including ramps. 
C 
          CASE (1) 
C           Calculate speed bin fractions which give the user 
C           supplied average speed for non-ramp freeways. 
            DO I=1,13 
C           All VMT is at 2.5 mph. 
              IF(USRSPD.EQ.2.5) THEN 
                DO IH=1,24 
                  SVMT(1,1,IH)=1.0 
                END DO 
                CYCLE 
              ENDIF 
C           General case for speeds >2.5 mph. 
              IF(USRSPD.GT.NSPD(I) .AND. 
     *           USRSPD.LE.NSPD(I+1)) THEN 
                BIN1=((1.0/USRSPD)-(1.0/NSPD(I+1))) / 
     *               ((1.0/NSPD(I))-(1.0/NSPD(I+1))) 
                BIN2=1.0-BIN1 
C           Set speed bin fractions for freeways. 
                DO IH=1,24 
                  SVMT(1,I+1,IH) = BIN2 
                  SVMT(1,I,IH) = BIN1 
                END DO 
                CYCLE 
              ENDIF 
            END DO 
C           Set roadway type to 100% freeway (no ramps). 
            DO IH=1,24 
              DO IV=1,28 
                FVMT(1,IH,IV) = 1.0 
                FVMT(2,IH,IV) = 0.0 
                FVMT(3,IH,IV) = 0.0 
                FVMT(4,IH,IV) = 0.0 
              END DO 
            END DO 
C 
C         Case 2 : All VMT occurs on freeways, including ramps. 
C 
          CASE (2) 
C           Check for inconsistent freeway and ramp vehicle fractions 
            IF( (USRVMT(1)<=0.00001.AND.USRVMT(4)<=0.00001) .OR. 
     *          (USRVMT(1)<=0.00001.AND.USRVMT(2)<=0.00001) .OR. 
     *          (1/USRSPD)==(USRVMT(4)/DEFAVG(4))           ) THEN 
              CALL QUITER(0.0,0,579,INERR) 
              USRVMT(1) = 0.342 
              USRVMT(2) = 0.498 
              USRVMT(3) = 0.130 
              USRVMT(4) = 0.030 
            ENDIF 
 
C           Determine the non-ramp freeway average speed. 
C           User supplied VMT fractions used in calculations. 
            FWYSPD=(1-USRVMT(4)/(USRVMT(1)+USRVMT(4))) / 
     *             ((1/USRSPD)-(USRVMT(4)/DEFAVG(4))) 
            IF(FWYSPD.GT.65.0) FWYSPD=65.0 
C           Calculate speed bin fractions which give the user 
C           supplied average speed for non-ramp freeways. 
            DO I=1,13 
C           All VMT is at 2.5 mph. 
              IF(FWYSPD.EQ.2.5) THEN 
                DO IH=1,24 
                  SVMT(1,1,IH)=1.0 
                END DO 
                CYCLE 
              ENDIF 
C           General case for speeds >2.5 mph. 
              IF(FWYSPD.GT.NSPD(I) .AND. 
     *           FWYSPD.LE.NSPD(I+1)) THEN 
                BIN1=((1.0/FWYSPD)-(1.0/NSPD(I+1))) / 
     *               ((1.0/NSPD(I))-(1.0/NSPD(I+1))) 
                BIN2=1.0-BIN1 
C           Set speed bin fractions for freeways. 
                DO IH=1,24 
                  SVMT(1,I+1,IH) = BIN2 
                  SVMT(1,I,IH) = BIN1 
                END DO 
                CYCLE 
              ENDIF 
            END DO 
C           Set roadway type to 100% only freeways and ramps. 
            DO IH=1,24 
              DO IV=1,28 
                FVMT(1,IH,IV) = USRVMT(1)/(USRVMT(1)+USRVMT(4)) 
                FVMT(2,IH,IV) = 0.0 
                FVMT(3,IH,IV) = 0.0 
                FVMT(4,IH,IV) = USRVMT(4)/(USRVMT(1)+USRVMT(4)) 
              END DO 
            END DO 
C 
C         Case 3 : All VMT occurs on arterial/collector roadways. 
C 
          CASE (3) 
C           Calculate speed bin fractions which give the user 
C           supplied average speed for arterial/collector roadways. 
            DO I=1,13 
C           All VMT is at 2.5 mph. 
              IF(USRSPD.EQ.2.5) THEN 
                DO IH=1,24 
                  SVMT(2,1,IH)=1.0 
                END DO 
                CYCLE 
              ENDIF 
C           General case for speeds >2.5 mph. 
              IF(USRSPD.GT.NSPD(I) .AND. 
     *           USRSPD.LE.NSPD(I+1)) THEN 
                BIN1=((1.0/USRSPD)-(1.0/NSPD(I+1))) / 
     *               ((1.0/NSPD(I))-(1.0/NSPD(I+1))) 
                BIN2=1.0-BIN1 
C           Set speed bin fractions for arterial/collector roadways. 
                DO IH=1,24 
                  SVMT(2,I+1,IH) = BIN2 
                  SVMT(2,I,IH) = BIN1 
                END DO 
                CYCLE 
              ENDIF 
            END DO 
C           Set roadway type to 100% arterial/collector roadways. 
            DO IH=1,24 
              DO IV=1,28 
                FVMT(1,IH,IV) = 0.0 
                FVMT(2,IH,IV) = 1.0 
                FVMT(3,IH,IV) = 0.0 
                FVMT(4,IH,IV) = 0.0 
              END DO 
            END DO 
C 
C         Case 5 : All VMT occurs on local roadways. 
C 
          CASE (5) 
C           Set roadway type to 100% local roadways. 
            DO IH=1,24 
              DO IV=1,28 
                FVMT(1,IH,IV) = 0.0 
                FVMT(2,IH,IV) = 0.0 
                FVMT(3,IH,IV) = 1.0 
                FVMT(4,IH,IV) = 0.0 
              END DO 
            END DO 
C 
C         Case 6 : All VMT occurs on freeway ramps. 
C 
          CASE (6) 
C           Set roadway type to 100% freeway ramps. 
            DO IH=1,24 
              DO IV=1,28 
                FVMT(1,IH,IV) = 0.0 
                FVMT(2,IH,IV) = 0.0 
                FVMT(3,IH,IV) = 0.0 
                FVMT(4,IH,IV) = 1.0 
              END DO 
            END DO 
C 
C         Case 4 : Average speed includes a mix of roadway types 
C                  specified by the user in USRVMT. 
C 
          CASE DEFAULT    ! Case 4 is the default. 
C           Determine the non-ramp freeway average speed adjustment from the 
C           adjusted arterial/collector speed stored in USRSPD. 
            FWYSPD=DEFAVG(1)*(USRSPD/DEFAVG(2)) 
C           Calculate speed bin fractions which give the user 
C           supplied average speed for non-ramp freeways. 
            DO I=1,13 
C           All VMT is at 2.5 mph. 
              IF(FWYSPD.EQ.2.5) THEN 
                DO IH=1,24 
                  SVMT(1,1,IH)=1.0 
                END DO 
                CYCLE 
              ENDIF 
C           General case for speeds >2.5 mph. 
              IF(FWYSPD.GT.NSPD(I) .AND. 
     *           FWYSPD.LE.NSPD(I+1)) THEN 
                BIN1=((1.0/FWYSPD)-(1.0/NSPD(I+1))) / 
     *               ((1.0/NSPD(I))-(1.0/NSPD(I+1))) 
                BIN2=1.0-BIN1 
C           Set speed bin fractions for freeways. 
                DO IH=1,24 
                  SVMT(1,I+1,IH) = BIN2 
                  SVMT(1,I,IH) = BIN1 
                END DO 
                CYCLE 
              ENDIF 
            END DO 
C           Calculate speed bin fractions which give the user 
C           supplied average speed for arterial/collector roadways. 
            DO I=1,13 
C           All VMT is at 2.5 mph. 
              IF(USRSPD.EQ.2.5) THEN 
                DO IH=1,24 
                  SVMT(2,1,IH)=1.0 
                END DO 
                CYCLE 
              ENDIF 
C           General case for speeds >2.5 mph. 
              IF(USRSPD.GT.NSPD(I) .AND. 
     *           USRSPD.LE.NSPD(I+1)) THEN 
                BIN1=((1.0/USRSPD)-(1.0/NSPD(I+1))) / 
     *               ((1.0/NSPD(I))-(1.0/NSPD(I+1))) 
                BIN2=1.0-BIN1 
C           Set speed bin fractions for arterial/collector roadways. 
                DO IH=1,24 
                  SVMT(2,I+1,IH) = BIN2 
                  SVMT(2,I,IH) = BIN1 
                END DO 
                CYCLE 
              ENDIF 
            END DO 
C           Set roadway type to 100% freeway (no ramps). 
            DO IH=1,24 
              DO IV=1,28 
                FVMT(1,IH,IV) = USRVMT(1) 
                FVMT(2,IH,IV) = USRVMT(2) 
                FVMT(3,IH,IV) = USRVMT(3) 
                FVMT(4,IH,IV) = USRVMT(4) 
              END DO 
            END DO 
      END SELECT 
C 
      RETURN 
      END 
