      FUNCTION DGADJF(IP,ICY,IFDS) 
C 
C  DGADJF computes the detergent gasoline adjustment to the basic emission factor 
C 
C  Called by FUEL. 
C 
C  Calls QUITER 
C 
C  Changes: (Last change first) 
C   
C  27 Nov 00 @EPA-elg    Turn off Detergent gas benefit if CAAFLG = 1. 
C  13 Jun 00 @DynTel-ddj Bug140, Bypassed calculation if ICY less than 
C            initial year (DGPISY) for implementing detergent gasoline 
C            rule. 
C  27 Jan 99 @DynTel-ZK  2-000 Explicit Typing 
C  16 Sep 98 @DynTel-ZK  2-000 Explicit Typing 
C  03 Sep 98 @DynTel-RJD 0-632  
C  24 Aug 98 @DynTel-RJD 0-632 Created DJADJF to simplfy fuel. 
C 
C  Input on call: 
C 
C    argument list: IP,ICY,IFDS 
C    common blocks: 
C     /DTGGAS/ DGPISY,DGPIEY,DGPHSF,FIDFRC,FIDEFF,IVDFRC,IVDEFF 
C     /FLAGS3/  CAAFLG, CAAMY 
C 
C 
C  Output on return: 
C 
C    function: DGADJF 
C    common blocks: 
C 
C  Local variable dictionary: 
C 
C   Name   Type              Description 
C  ------  ----  ------------------------------------------------------- 
C  DPHASE   R    local detergent gasoline phase-in factor 
C 
      IMPLICIT NONE 
C 
      INCLUDE 'DTGGAS.I' 
      INCLUDE 'FLAGS3.I' 
C 
C     declare parameter list 
C 
      INTEGER ICY, IFDS, IP 
C 
C     declare local variables 
C 
      REAL    DPHASE 
      INTEGER INERR 
C 
C     declare functions 
C 
      REAL DGADJF 
C 
C  Calculate the phase-in factor 
C 
C  13-Sep-1996 @DynTel-yc Request #2-630 The detergent gas phase-in year was 
C              divided into phase-in start year(DGPISY) and end year(DGPIEY). 
C 
C  Since call to QUITER in this routine stops processing INERR is a 
C  dummy variable for QUITER and can be set to an arbitrary value. 
C 
      INERR = 1 
C 
      IF (ICY.LT.DGPISY .OR. (CAAFLG.EQ.1 .AND. ICY.GE.CAAMY)) THEN 
          DGADJF = 1.0 
      ELSE 
        IF (ICY.GE.DGPISY .AND. ICY.LT.DGPIEY) THEN 
          DPHASE = DGPHSF 
        ELSE 
          DPHASE = 1.0 
        END IF 
C 
C  Compute the detergent gasoline adjustment to the basic emission factor 
C 
        IF (IFDS.EQ.1) THEN 
C         *** Carbureted 
          DGADJF = 1.0 - IVDFRC*IVDEFF(IP)*DPHASE 
        ELSE 
C         *** TBI or PFI 
          DGADJF = 1.0 - (IVDFRC*IVDEFF(IP) + FIDFRC*FIDEFF(IP))*DPHASE 
        END IF 
C 
C  If DGADJF is not reasonable, call QUITER and abort the run. 
C 
        IF (DGADJF .GT. 1.0 .OR. DGADJF .LT. 0.001) THEN 
          CALL QUITER(DGADJF,IP,85,INERR) 
        ENDIF 
C 
      ENDIF 
C 
      RETURN 
      END 
 
