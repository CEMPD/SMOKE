      SUBROUTINE CKOBDIM(INERR) 
C 
C  CKOBDIM checks for existence of an Exhaust I/M program operating on 
C  1996 and later model year vehicles, and the existence of an OBD I/M 
C  and/or an EVAP OBD operating on pre-1996 model year vehicles. The 
C  routine also checks that no OBD tests are made on pre-2007 HD vehicles. 
C 
C  Called by TESTIM 
C 
C  Changes: (Last change first) 
C 
C 09 Jan 02 @EPA-djb   : Bug 356, No OBD I/M for pre-2007 HDGV. 
C 09 May 01 @EPA-elg   : POST96IM Flag has been effectively disabled. 
C 16 Apr 00 @EPA-elg   : Allows POST 1996 I/M Input 
C 28 Nov 00 @EPA-elg   : Allows NO CAA I/M Input 
C 15 Sep 00 AIR Task 08: Updated code to handle evap I/M. 
C 25 Aug 00 AIR Task 08: Added code to trap invalid Evap I/M inputs. 
C 07 Jul 00 AIR Task 05: Removed unused RC from parameter list. 
C 18 May 00 @EPA-ELG  This is a new routine for Mobile6. 
C 
C 
C  Input on call: 
C 
C    common blocks: 
C    /CHKLST/ IMCHECK,MAXIMPGM 
C    /FLAGS2/ IMFLAG,EVIMFLAG 
C    /FLAGS3/ CAAFLG 
C    /IMPAR1/ ICYIM,MODYR 
C    /IMPAR2/ ILDT 
C    /OIUCOM/ IOUERR 
C    /IVTYPE/ MAXVEH 
C    /VVNAMS/ VVGAS     (IVTYPE.I) 
C    /VVTYPS/ VTGASBUS  (IVTYPE.I) 
C 
C  Output on return: 
C    argument list: INERR 
C 
C    Function: None. (Error messages are written to the error file.) 
C 
C 
C  Notes: 
C 
      IMPLICIT NONE 
C 
      INCLUDE 'CHKLST.I' 
      INCLUDE 'FLAGS2.I' 
      INCLUDE 'FLAGS3.I' 
      INCLUDE 'IMPAR1.I' 
      INCLUDE 'IMPAR2.I' 
      INCLUDE 'IOUCOM.I' 
      INCLUDE 'IVTYPE.I' 
C 
      INTEGER, INTENT(INOUT) :: INERR 
C 
      INTEGER IVIM   !Local gasoline vehicle classes (no MC). 
      INTEGER PGM1   !Local program number. 
      INTEGER JHDGV  !Are HDGV included in OBD Evap for pre-2007 MYRs? (0:No,1:Yes) 
C 
C 
C 
      IF (IMFLAG.LE.1 .AND. EVIMFLAG.LE.1) GOTO 99 
C 
      DO PGM1 = 1, MAXIMPGM - 1 
C 
C       If no data was entered for this program, skip it. 
C 
        IF (IMCHECK(1,PGM1).EQ.1) CYCLE 
C 
C       Exhaust I/M programs are NOT supported in MOBILE6 for 1996 and later 
C       model year vehicles.  Exception is if the No Clean Air Act Flag is set 
C       to 1. In these cases only the Idle test and the 2500/Idle tests are 
C       supported.  Evap I/M tests are not supported. 
C 
        IF(CAAFLG.NE.1) THEN 
C 
C       The POST96IM flag has been effectively disabled by always setting it to 1. 
C       This allows non-OBD tests for 1996 and newer vehicles. 
C 
          IF(ITEST(PGM1).LE.11 .AND. MODYR(1,PGM1).GT.1995) THEN 
            IF(POST96IM.NE.1) CALL QUITER(0.0,0,546,INERR)   !Exhaust I/M not supported 
          ELSEIF(ITEST(PGM1).LE.11 .AND. MODYR(2,PGM1).GT.1995) THEN 
            IF(POST96IM.NE.1) CALL QUITER(0.0,0,546,INERR)   !Exhaust I/M not supported 
C 
C 
C       OBD I/M programs are NOT supported in MOBILE6 for pre-1996 
C       model year vehicles. 
C 
          ELSEIF(ITEST(PGM1).EQ.12 .AND. MODYR(1,PGM1).LT.1996) THEN 
           CALL QUITER(0.0,0,545,INERR)    !OBD not supported for pre-1996 MYs 
          ELSEIF(ITEST(PGM1).EQ.12 .AND. MODYR(2,PGM1).LT.1996) THEN 
           CALL QUITER(0.0,0,545,INERR)    !OBD not supported for pre-1996 MYs 
          ENDIF 
C 
C       EVAP OBD and EVAP OBD & GC programs are NOT supported in MOBILE6 for 
C       pre-1996 model year vehicles. 
C 
          IF(ITEST(PGM1).GE.20 .AND. ITEST(PGM1).LE.21 .AND.  
     &           MODYR(1,PGM1).LT.1996) THEN 
           CALL QUITER(0.0,0,556,INERR)    !Evap OBD not supported for pre-1996 MYs 
          ELSEIF(ITEST(PGM1).GE.20 .AND. ITEST(PGM1).LE.21 .AND. 
     &           MODYR(2,PGM1).LT.1996) THEN 
           CALL QUITER(0.0,0,556,INERR)    !Evap OBD not supported for pre-1996 MYs 
          ENDIF 
C 
C       OBD systems for Evap do not appear on HDGV until the 2007 model year. 
C       Do not allow Evap OBD tests for HDGV before the 2007 model year. 
C 
          IF(ITEST(PGM1).EQ.20 .OR. ITEST(PGM1).EQ.21) THEN 
            JHDGV=0 
            IF(MODYR(1,PGM1).LT.2007) THEN 
              DO IVIM=6,14 
                IF(ILDT(IVIM,PGM1).NE.1) JHDGV=1          !Evap OBD I/M tests are NOT supported 
              END DO                                      !for pre-2007 model year HDGV. 
            ENDIF 
            IF(JHDGV.EQ.1) CALL QUITER(0.0,PGM1,531,INERR) 
          ENDIF 
C 
        ELSEIF(CAAFLG.EQ.1) THEN 
C 
C       If the CAA is disabled, then non-OBD tests are allowed for 
C       1996 and newer vehicles, however OBD tests are not allowed. 
C 
          IF(ITEST(PGM1) .EQ. 12 .OR. ITEST(PGM1).EQ. 20  .OR. 
     &       ITEST(PGM1) .EQ. 21) THEN 
            IF(POST96IM.NE.1) CALL QUITER(0.0,0,565,INERR)   !OBD I/M tests are NOT supported if CAAFLG=1 (Inerr=1) 
          ELSEIF(ITEST(PGM1) .GE.4 .AND. ITEST(PGM1).LE. 11) THEN 
            IF(POST96IM.NE.1) CALL QUITER(0.0,0,566,INERR)   !Program produces a warning if these tests are selected. 
          ENDIF 
C 
        ENDIF 
C 
      END DO       !End of loop over PGM1 
C 
C 
   99 CONTINUE 
      RETURN 
      END 
