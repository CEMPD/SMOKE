      SUBROUTINE GETFUEL(INPREC,NSCEN,INERR,RC) 
C 
C  GETFUEL sets the Fuel Program index (FUEL_P) and for RFG fuel, sets 
C  the RFG region (RFG_REGION). 
C 
C 
C  Called by DSPTCH. 
C 
C 
C  Calls NXTTOK, QUITER. 
C 
C 
C  Changes: (Last change first) 
C 
C  24 Aug 01 AIR Task 15: Added trap for FUEL PROGRAM Option 2 and  
C            Air Toxic run. 
C  11 Sep 00 @EPA-elg Added min/max check on sulfur level input. 
C  08 Aug 00 @EPA-djb bug216 Removed the GLS_PPM storage. 
C  07 aug 00 @EPA-BAG Bug217 Call to QUITER with message #614(user supplied 
C            diesel sales fractions) changed to message #616(user supplied 
C            sulfur content) 
C  10 Jul 00 @EPA-ddj Bug184 Added assignment to GLS_PPM. 
C  24 May 00 @EPA-ddj E9 Added error messages for invalid user input 
C            of CURRENT and MAXIMUM sulfur contents (FUEL PROGRAM = 4). 
C  20 Apr 00 @EPA-ddj E7 Renamed from RFGFLG and added Fuel Program 
C            index. 
C  21 Oct 1998 @DynTel-ZK 2-000 Explicit Typing 
C  Input on call: 
C 
C 
C   argument list: INPREC 
C 
C  Output on return: 
C 
C    argument list: INPREC, NSCEN, INERR, RC 
C    common blocks: 
C    /RFORM1/ GLRFG,RFGFLG 
C    /RFORM2/ GLRFGSYR,RFGSYR 
C    /YEARS4/ IYEND 
C 
C  Local variable / array dictionary: 
C 
C   Name   Type              Description 
C  ------  ----  ------------------------------------------------------- 
C  DATA     C    Holds the data value from the input record. 
C  S_FUEL   R    Holds the current/maximum gasoline sulfur content 
C                by calendar/model year. 
C 
C 
C  Notes: 
C 
C  27 Mar 98 @DynTel-MLA 2-659  New routine for Mobile6 (the transitional 
C            labeled inputs version). This subroutine sets the RFGFLG. 
C  24 Jun 98 @DynTel-bag 2-000 Changed QUITER argument RVALUE (28 --> 5)  
C            to represent the changes made in the NAMFLG array. 
C 
      IMPLICIT NONE 
C 
      INCLUDE 'ATOX1.I' 
      INCLUDE 'IOUCOM.I' 
      INCLUDE 'RFORM1.I' 
      INCLUDE 'RFORM2.I' 
      INCLUDE 'SULFUR.I' 
      INCLUDE 'YEARS4.I' 
C 
      CHARACTER(80)   :: DATAFIELD 
      CHARACTER(*)    :: INPREC 
      CHARACTER       :: RFG_LOCAL 
C 
      INTEGER, INTENT(INOUT)  :: INERR 
      INTEGER, INTENT(IN)     :: NSCEN 
      INTEGER, INTENT(INOUT)  :: RC 
C 
      INTEGER         :: I 
      INTEGER         :: PTR1 
      INTEGER         :: PTR2 
C 
      PTR2 = 21 
      CALL NXTTOK(INPREC,DATAFIELD,PTR1,PTR2) 
      READ(DATAFIELD,*, IOSTAT=RC) FUEL_P 
      IF(FUEL_P.LT.1 .OR. FUEL_P.GT.4) THEN 
         INERR = INERR + 1 
         CALL QUITER(4.,FUEL_P,1,INERR) 
      ELSE 
         CALL QUITER(0.0,0,616,INERR) 
      ENDIF 
C 
C  Prevent using FUEL PROGRAM Option 2 (RFG) with any Air Toxic Runs. 
C 
      IF(FUEL_P.EQ.2 .AND. TOXFLG) THEN 
        CALL QUITER(0.0,0,115,INERR) 
        RC = 2 
        RETURN 
      END IF 
C 
C  Process the input 
C 
      IF(RC.EQ.0 .AND. FUEL_P.EQ.2) THEN 
         RFG_REGION = 0 
         CALL NXTTOK(INPREC,DATAFIELD,PTR1,PTR2) 
         READ (DATAFIELD,*,IOSTAT=RC) RFG_LOCAL 
         IF(RFG_LOCAL.EQ.'s' .OR. RFG_LOCAL.EQ.'S') THEN 
            RFG_REGION = 1 
         ELSEIF(RFG_LOCAL.EQ.'n' .OR. RFG_LOCAL.EQ.'N') THEN 
            RFG_REGION = 2 
         END IF 
         IF(RFG_REGION.NE.1 .AND. RFG_REGION.NE.2) THEN 
            CALL QUITER(5.,IACHAR(RFG_LOCAL),1,INERR) 
         ELSE 
            RFGFLG = 2                ! RFG implemented 
            IF(NSCEN .EQ. 0) THEN 
              GLRFG = 2 
              GLFUEL_P = FUEL_P 
            END IF 
         END IF 
      ELSEIF(RC.EQ.0 .AND. FUEL_P.EQ.4) THEN !User supplied sulfur levels. 
            READ(IOUGEN, *, IOSTAT=RC) 
     &         S_FUEL(1, 4, 1:16, 2)   ! Read average S ppm levels 
            IF(RC.EQ.0) THEN 
               READ(IOUGEN, *, IOSTAT=RC) 
     &         S_FUEL(2, 4, 1:16, 2)    ! Read maximum S ppm levels 
               DO I = 1, 16             ! Set No Tier2 ppm = Tier2 ppm 
                 IF(S_FUEL(1, 4, I, 2).LT.30.0 .OR. 
     &              S_FUEL(1, 4, I, 2).GT.1000.0) 
     &                CALL QUITER(S_FUEL(1,4,I,2),0,554,INERR) 
                 IF(S_FUEL(2, 4, I, 2).LT.30.0 .OR. 
     &              S_FUEL(2, 4, I, 2).GT.1000.0) 
     &                CALL QUITER(S_FUEL(2,4,I,2),0,554,INERR) 
                 IF(S_FUEL(1, 4, I, 2).GT.S_FUEL(2, 4, I, 2)) 
     &                CALL QUITER(S_FUEL(2,4,I,2),0,555,INERR) 
                 S_FUEL(1, 4, I, 1) = S_FUEL(1, 4, I, 2) 
                 S_FUEL(2, 4, I, 1) = S_FUEL(2, 4, I, 2) 
               END DO 
C 
               IF(NSCEN .EQ. 0) THEN 
                 GLFUEL_P = FUEL_P 
               END IF 
            ELSE 
               WRITE(IOUERR,200) 
               INERR = INERR + 1 
               RC = 0 
            END IF 
            IF(RC.NE.0) THEN 
               WRITE(IOUERR,210) 
               INERR = INERR + 1 
            END IF 
      ELSEIF(RC.EQ.0 .AND. (FUEL_P.EQ.1 .OR. FUEL_P.EQ.3)) THEN 
         IF(NSCEN .EQ. 0) THEN 
           GLFUEL_P = FUEL_P 
         END IF 
      END IF 
  200 FORMAT(9x,' Error reading Current sulfur levels - FUEL PROGRAM 4') 
  210 FORMAT(9X,' Error reading Maximum sulfur levels - FUEL PROGRAM 4') 
C 
C  Error or EOF on any attempted read. 
C 
      IF(RC.NE.0)THEN 
         INERR = INERR + 1 
      END IF 
C 
      RETURN 
      END 
