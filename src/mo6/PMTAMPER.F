      SUBROUTINE PMTAMPER(ICY, MY, IV) 
C 
C     This subroutine finds the proper misfueling rate for vehicle 
C     class IV, age IAGE, and whether or not there is an IM program,  
C     determined by IMFLAG.    
C 
C     This subroutine also finds the catalyst removal rate  
C 
C Changes: 
C 
C 20 Dec 01 AIR Task 15-28: Corrected a long-standing PART5 logic error which 
C                        allowed misfueling to continue even though leaded 
C                        fuel availability ended after CY1994.         
C 21 Feb 01 AIR Task 11: Added OBD effect to tampering. 
C 10 Feb 01 AIR Task 11: New MOBILE6 routine. 
C 10 Feb 01 AIR Task 11: Updated and merged into MOBILE6.1. Explicit 
C                        use of PART5 vehicle class IVPART. Uses 
C                        MOBILE6 CUMMIL array. 
C 10 Nov 00 AIR Task 06: Updated code to MOBILE6 standards. 
C 
C     Called by: PMLEADED, PMGSO4SO2, PMGASSL 
C 
C     Calls: none 
C 
C     Input on call: 
C     /TAMBLK/ TAMRAT(2,3,3,2,2) 
C     /CUMCOM/ CUMMIL(12,25) 
C 
C     Output on return: TAMFRC, RMIS 
C 
C  Local variable dictionary: 
C 
C   Name     Type             Description 
C  ------    ----  ----------------------------------------------------- 
C  IAGE       I    age 
C  ICY        I    calendar year 
C  IG         I    tampering group pointer 
C  IMFLAG     I    tampering I/M (hardcoded to 1 (no I/M) for MOBILE6 
C  IV         I    vehicle class 
C  IVA        I    vehicle class switch, LDGV=1, all others=2 
C  IVPART     I    MOBILE6 vehicle class mapped into PART6 space 
C  MY         I    model year 
C  XMILE      R    temporary mileage 
C  XXMILE     R    temporary mileage, above 50,000 miles 
C 
      IMPLICIT NONE 
C 
      INCLUDE 'CUMCOM.I'   ! CUMMIL 
      INCLUDE 'IVTYPE.I'   ! VM6M5B, vtldgt4 
      INCLUDE 'LEAD.I'     ! TAMRAT 
      INCLUDE 'PART2.I'    ! TAMFRC, RMIS 
C 
C  Declare parameter list. 
C 
      INTEGER ICY 
      INTEGER IV 
      INTEGER MY 
C 
C  Declare local variables/arrays. 
C 
      INTEGER IAGE 
      INTEGER IG 
      INTEGER IMFLAG 
      INTEGER IVA 
      INTEGER IVPART 
      REAL    XMILE 
      REAL    XXMILE 
C 
C 
      TAMFRC = 0.0 
      RMIS = 0.0 
C 
C  Map MOBILE6 vehicle classes into PART5 space. GasBusses are 
C  treated like HDGV. 
C 
      IVPART=VM6M5B(IV) 
C 
C  Special cases return 0.0 
C  EPA added 1996+MY LDGV-LDGT4 (MOBILE6) to account for OBD. 
C 
      IF (IVPART .GT. 4) RETURN 
      IF (MY .LT. 1975) RETURN 
      IF (IVPART .EQ. 4 .AND. MY .LT. 1982) RETURN 
C 
      IF (IV.LE.VTLDGT4.AND.MY.GE.1996) RETURN 
C 
C  Set PART5 IMFLAG to 1 (no I/m). I/M has no effect on PM in MOBILE6.1 
C 
      IMFLAG=1 
C 
C  Need to calculate age to index into mileage array 
C 
      IAGE = ICY - MY + 1 
C 
C  Mileage is measure in terms of 10,000 miles 
C  (Use MOBILE6 CUMMIL array) 
C 
      XMILE = CUMMIL(IAGE,IV)/10000. 
C 
C  XXMILE is the above 50,000 miles amount 
C 
      XXMILE = XMILE - 5.0 
C 
C  Set the second mileage to zero when the mileage is less than 50,000 miles 
C 
      IF(XXMILE.LT.0.0) THEN  
        XXMILE=0.0 
C 
C  Capping the mileage at 130,000 miles 
C 
      ELSEIF(XXMILE.GT.8.0) THEN 
        XXMILE=8.0 
      ENDIF 
      IF(XXMILE.GT.0.0) XMILE=5.0 
C 
C  Discover the model year group and set the IG pointer 
C 
      IF(MY.LE.1980) THEN 
        IG=1 
      ELSEIF(MY.GE.1984) THEN 
        IG=3 
      ELSE 
        IG=2 
      ENDIF 
C 
C  Choose the vehicle class (PART5 definition) cars or trucks 
C 
      IVA=IVPART 
      IF(IVA.NE.1) IVA=2 
C 
C  Calculate the tampering rates 
C 
      TAMFRC = TAMRAT(IMFLAG,1,IG,IVA,1) +  
     2         TAMRAT(IMFLAG,2,IG,IVA,1) * XMILE + 
C 
C  Multiple the second rate times the first rate 
C 
     3   TAMRAT(IMFLAG,2,IG,IVA,1) * TAMRAT(IMFLAG,3,IG,IVA,1) * XXMILE 
C 
C  Calculate the misfueling rates. Leaded fuel is not available for 
C  CY1995+. Therefore, only allow misfueling to occur for model years 
C  which could have been misfueled. 
C 
      IF(ICY-IAGE.LT.1995)  
     1 RMIS =   TAMRAT(IMFLAG,1,IG,IVA,2) +  
     2          TAMRAT(IMFLAG,2,IG,IVA,2) * XMILE + 
C 
C  Multiple the second rate times the first rate 
C 
     3   TAMRAT(IMFLAG,2,IG,IVA,2) * TAMRAT(IMFLAG,3,IG,IVA,2) * XXMILE 
C 
C  If rates are less than zero reset the rates to 0.0 
C 
      IF (RMIS .LT. 0.0) RMIS = 0.0 
      IF (TAMFRC .LT. 0.0) TAMFRC = 0.0 
C 
C  If rates are greater than one reset the rates to 1.0 
C 
      IF (RMIS .GT. 1.0) RMIS = 1.0 
      IF (TAMFRC .GT. 1.0) TAMFRC = 1.0 
C 
C  Only need to make this check once here in PMTAMPER 
C 
      IF(RMIS-TAMFRC .LT. 0.0) RMIS=TAMFRC 
C 
      RETURN  
      END  
