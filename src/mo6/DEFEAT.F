      SUBROUTINE DEFEAT(MY,ICY,IV,IROAD,DDNOX) 
C 
C  Called by : HRLOOP 
C 
C  Changes: (Last change first) 
C 
C 29 Oct 01 @EPA-djb    Changed OR to AND in MY conditional for rebuild effects. 
C 25 May 01 @EPA-elg    QUITER message 191 changed to 567 
C 13 Feb 01 @EPA-elg    Added bus category to IVV determination. 
C 08 Feb 01 @EPA-elg    Changed logic for IMY mapping for the NO PULL AHEAD 
C                       case and removed RETURN. Added logic to include 
C                       IRBLD, IRBSTRT and IRBEND to rebuild calculation. 
C 01 May 00 @EPA- BG    Removed use of DATA statement initializations 
C 22 Dec 99 @Dyntel-nh  1-004 removed diagnostic parameters and references. 
C 07 Dec 99 @DynTel-NH  1-004  added subroutine DEFEAT_STAT 
C 25 Oct 99 @DynTel-NH  1-004  New for the HDDV defeat task 
C 
C..Rebuild modeling options using variable IBUILD. 
C    IBUILD = 0  -  Default option - Rebuild effects are included at a 
C                                    90% penetration for MY 1994-98 
C    IBUILD = 1  -  No Rebuild effects are included. 
C    IBUILD = 2  -  Rebuild effects are included for MY 1994-98 at a 
C                   user specified level (RBLEV). 
C 
C    IPULL = 0   -  Default option - Pull ahead effects are included for 
C                   MY 2002 and 2003. Negative emission increases. 
C    IPULL = 1   -  Pull ahead effects are not included for MY 2002 and 
C                   2003.  Defeat device effects are equal to MY 2001 effects. 
C 
C 
      USE DDDATA 
C 
      IMPLICIT NONE 
C 
      INTEGER, INTENT(IN) :: MY 
      INTEGER, INTENT(IN) :: ICY 
      INTEGER, INTENT(IN) :: IV 
      INTEGER, INTENT(IN) :: IROAD 
 
      INTEGER :: IMY, IVV                      ! Local variables 
C 
      REAL, INTENT(OUT) :: DDNOX 
C 
      REAL :: DDNOXX                           ! Local variables 
C 
C 
      DDNOX=0.0 
      IF(MY.LT.1988.OR.MY.GT.2003.OR.ICY.LT.1988.OR. 
     *   IROAD.EQ.4) THEN 
        DDNOX = 0.0                             ! Out of Range: No defeat device effects 
        RETURN 
      ENDIF 
C 
      IF(MY.GE.2004) THEN 
         IMY = 1 
      ELSEIF(MY.LE.1987) THEN 
         IMY = 18 
      ELSE 
         IMY = 2004 - MY + 1 
      ENDIF 
C 
      IF((MY.EQ.2002.OR.MY.EQ.2003).AND.IPULL.EQ.1) THEN 
         IMY = 2004 - 2001 + 1       ! Always takes 2001 value for no pull ahead. 
      ENDIF 
C 
      IF(IV.GE.16.AND.IV.LE.19) THEN !was OR 
         IVV = 1 
      ELSEIF(IV.EQ.20.OR.IV.EQ.21) THEN 
         IVV = 2 
      ELSEIF(IV.EQ.22) THEN 
         IVV = 3 
      ELSEIF(IV.EQ.23) THEN 
         IVV = 4 
      ELSEIF(IV.EQ.26.OR.IV.EQ.27) THEN 
         IVV = 5 
      ENDIF 
C 
C     The base NOx excess does not contain rebuild effects. 
C 
      DDNOXX=NEFFCT(IMY,IVV) 
C 
C     Effect of rebuilds: 
C     No rebuilds can occur before calendar year 2000 because 
C     the defeat device litigation had not yet been concluded. 
C 
      IF(ICY.GT.2000) THEN 
        IF(MY.GE.IRBSTRT.AND.MY.LE.IRBEND) THEN 
          IF(ICY.GT.MY+IRBLD(IVV)) THEN 
            IF(IBUILD.EQ.2) THEN 
              DDNOXX = RBLEV * 
     *                ((RBLD(IMY,IVV)-0.10*NEFFCT(IMY,IVV))/0.90) 
     *                + (1-RBLEV)*NEFFCT(IMY,IVV) 
            ELSEIF(IBUILD.EQ.1) THEN 
              DDNOXX = NEFFCT(IMY,IVV) 
            ELSE 
              DDNOXX = RBLD(IMY,IVV) 
            ENDIF 
          ENDIF 
        ENDIF 
      ENDIF 
C 
C     Effects are adjusted for defeat device activation by roadway and percentage 
C     of defeat device engines in the fleet. 
C 
      IF(IROAD.EQ.1) THEN 
         DDNOX = DDNOXX * FLEET(IMY,IVV) * URBAN(IMY,IVV) 
      ELSEIF(IROAD.EQ.2) THEN 
         DDNOX = DDNOXX * FLEET(IMY,IVV) * ARTERL(IMY,IVV) 
      ELSEIF(IROAD.EQ.3) THEN 
         DDNOX = DDNOXX * FLEET(IMY,IVV) * EXPRES(IMY,IVV) 
      ENDIF 
C 
C 
      RETURN 
      END 
 
      SUBROUTINE RBEFFCTS(INPREC,INERR,RC) 
C 
C Changes last first 
C 
C 18 Oct 1999 @Duntel-nh. New subroutine added to set DDDATA module variable IBUILD 
C                         when "REBUILD EFFECTS" label or "NO REBUILD" label is 
C                         found in input file. 
C 
C Called by DSPTCH 
C 
C Calls NXTTOK 
C 
C Argument List 
C 
C CHARACTER*(*) INPREC   record from input file 
C INTEGER       INERR 
C INTEGER       RC 
C MODULE DDDATA variables 
C 
C INTEGER       IBUILD 
C REAL          RBLEV 
C 
C Local variables 
C 
C CHARACTER*10  STR1  Return string from NXTTOK 
C INTEGER       PTR1  Returned index from NXTTOK 
C INTEGER       PTR2  INDEX for search in NXTTOK 
C REAL          DELTA Used for REAL xero comparisons 
C 
      USE DDDATA 
      IMPLICIT NONE 
      CHARACTER*(*) INPREC 
      INTEGER INERR 
      INTEGER RC 
      CHARACTER*10 STR1 
      INTEGER PTR1,PTR2 
      REAL, PARAMETER :: DELTA=0.000001 
C 
      IF(INDEX(INPREC,'REBUILD EFFECTS')>0) THEN 
        IBUILD=2 
        PTR2=21 
        CALL NXTTOK(INPREC,STR1,PTR1,PTR2) 
        READ(STR1,'(F4.0)')RBLEV 
        IF ((ABS(RBLEV).LT.DELTA).OR.(RBLEV.GT.0.90)) THEN 
           CALL QUITER(RBLEV,0,567,INERR) 
        END IF 
      ELSEIF(INDEX(INPREC,'NO REBUILD')>0)THEN 
        IBUILD=1 
      ELSEIF(INDEX(INPREC,'NO DEFEAT')>0)THEN 
        DEFEAT_FLAG=.FALSE. 
      ENDIF 
      RC=0 
      RETURN 
      END SUBROUTINE 
 
      SUBROUTINE NOXPUL(RC) 
C 
C CHANGES LAST FIRST 
C 
C 07 Jul 00 AIR Task 05: Removed unused INPREC and INERR from parameter list. 
C 18 OCT 1999 @DYNTEL-NH. NEW SUBROUTINE ADDED TO SET DDDATA MODULE VARIABLE IPULL 
C                         WHEN "NO NOX PULL AHEAD" LABEL FOUND IN INPUT FILE 
C CALLED BY DSPTCH 
C 
C ARGUMENT LIST 
C 
C INTEGER       RC 
C 
C MODULE DDDATA VARIABLES 
C 
C INTEGER       IPULL 
C REAL          RBLEV 
C 
      USE DDDATA 
      IMPLICIT NONE 
      INTEGER RC 
      IPULL=1 
      RC=0 
      RETURN 
      END SUBROUTINE 
 
      SUBROUTINE INIDEFEAT() 
C 
C CALLED BY REINIT TO RESET DEFAULT VALUES FOR EACH RUN. 
C MODULE DDDATA VARIABLES 
C 
C INTEGER       IPULL 
C INTEGER       IBUILD 
C REAL          RBLEV 
      USE DDDATA 
      IMPLICIT NONE 
      IPULL=0 
      IBUILD=0 
      RBLEV=0.90 
      DEFEAT_FLAG=.TRUE. 
      FIRST_DEFEAT=.TRUE. 
      RETURN 
      END SUBROUTINE 
 
      LOGICAL FUNCTION DEFEAT_DEVICE() 
      USE DDDATA 
      IMPLICIT NONE 
      DEFEAT_DEVICE=DEFEAT_FLAG 
      END 
 
      SUBROUTINE DEFEAT_STAT(IOUREP) 
C 
C 
C CHANGES LAST FIRST 
C 
C 22 DEC 99 @DYNTEL-NH. INCREASED VERBOSITY OF MESSAGES TO REPORT FILE. 
C 08-DEC-1999 @DYNTEL-NH. NEW SUBROUTINE ADDED TO REPORT DEFEAT OPTIONS. 
C 
C CALLED BY HRLOOP TO REPORT DEFEAT OPTIONS TO REPORT FILE 
C 
C MODULE DDDATA VARIABLES 
C 
C INTEGER       IPULL 
C INTEGER       IBUILD 
C REAL          RBLEV 
C 
C ARGUMENT LIST 
C 
C INTEGER IOUREP LOGICAL DEVICE NUMBER FOR REPORT FILE 
C 
C OUTPUTS 
C 
C STATUS TO REPORT FILE 
C 
      USE DDDATA 
      IMPLICIT NONE 
      INTEGER IOUREP 
      IF (IPULL.EQ.1) THEN 
         WRITE(IOUREP, 
     C    '(" NO HDDV DEFEAT DEVICE PULL AHEAD EFFECTS ARE PRESENT.")') 
      END IF 
      SELECT CASE (IBUILD) 
          CASE (1) 
            WRITE(IOUREP, 
     C       '( " NO HDDV DEFEAT DEVICE REBUILD PROGRAM IN EFFECT." )') 
          CASE (2) 
            WRITE(IOUREP,'( " HDDV DEFEAT DEVICE EFFECTS ARE PRESENT. ", 
     C             "THE REBUILD FRACTION IS ",F4.2,"." )')RBLEV 
          CASE DEFAULT 
      END SELECT 
      IF (.NOT. DEFEAT_FLAG) THEN 
         WRITE (IOUREP, 
     C    '(" NO HDDV DEFEAT DEVICE EFFECTS ARE PRESENT.")') 
      END IF 
      END SUBROUTINE 
