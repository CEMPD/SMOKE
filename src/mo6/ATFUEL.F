      SUBROUTINE ATFUEL(MYR,IV,ICY,IFDS,INH,IDX,ISTD,IVT,EF,FUEL_CF) 
C 
C  ATFUEL computes the detergent gas, sulfur, RVP and oxy fuel, effect of 
C  industry average fuel and reformulated gas adjustments.  Due to the 
C  introduction of oxidation and 3-way catalyst on LDTs, special 
C  handling of those technologies was needed. (Air Toxic version of 
C  FUEL.) 
C 
C  Called by ATHRLOOP. 
C 
C  Calls NMOGFID 
C 
C  Changes: (Last change first) 
C 
C  09 Jul 01 AIR Task 15: New Subroutine (Based on FUEL) 
C 
C  Input on call: 
C 
C    argument list: MYR,IV,ICY,IFDS,INH,IDX,ISTD,EF,FUEL_CF 
C    common blocks: 
C     /FLAGS1/ DTGFLG, OXYFLG 
C     /FLAGS2/ NEWFLG 
C     /OXY1/   SHRMKT,OXYCNT 
C     /RFORM1/ JPHASE,SEAFLG,RFGON 
C     /RFORM2/ RFG_CF 
C     /VVNAMS/ VVDSL, VVHDGV, VVLDGT12, VVLDGT34  (IVTYPE.I) 
C     /VVTYPS/ VTLDGV, VTMC  (IVTYPE.I) 
C 
C  Output on return: 
C 
C    common blocks: 
C     None 
C 
C  Local array subscripts: 
C 
C     COHDGV(10)  - COHDGV(EMY) 
C     COMC(10)    - COMC(EMY) 
C     EMY depends on the emitter type (Normals and High) and 
C          model years  
C 
C 
C  Local variable dictionary: 
C 
C   Name   Type              Description 
C  ------  ----  ------------------------------------------------------- 
C  COHDGV   R    Hardcoded oxygen fuel correction factors for HDGV 
C  COMC     R    Hardcoded oxygen fuel correction factors for MC 
C  OXYCT    R    local variable for the oxygen content 
C 
      IMPLICIT NONE 
C 
      INCLUDE 'BASEQ9.I' 
      INCLUDE 'FLAGS1.I' 
      INCLUDE 'FLAGS2.I' 
      INCLUDE 'IVTYPE.I' 
      INCLUDE 'MAXIMA.I' 
      INCLUDE 'OXY1.I' 
      INCLUDE 'RFORM1.I' 
      INCLUDE 'RFORM2.I' 
      INCLUDE 'RVPEX2.I' 
C 
      INTEGER, INTENT(IN)      :: ICY 
      INTEGER, INTENT(IN)      :: IDX 
      INTEGER, INTENT(IN)      :: IFDS 
      INTEGER, INTENT(IN)      :: INH 
      INTEGER, INTENT(IN)      :: ISTD 
      INTEGER, INTENT(IN)      :: IV 
      INTEGER, INTENT(IN)      :: IVT 
      INTEGER, INTENT(IN)      :: MYR 
C 
      REAL, DIMENSION(2), INTENT(IN)                :: EF 
      REAL, DIMENSION(2,2,10,3,24,2), INTENT(INOUT) :: FUEL_CF 
C 
      REAL C_FACTOR 
      REAL, DIMENSION(2) :: F_CF 
C 
      INTEGER IH 
      INTEGER ISR 
      INTEGER IVTAM 
      INTEGER, PARAMETER :: IP=1   !HC 
C 
C    External Functions 
C 
      REAL DGADJF 
      REAL HCOXY 
      REAL SUL_COR 
C 
C  If vehicle is diesel fueled, no adjustment is necessary. 
C 
      IF(VVDSL(IV) .NE. 1) THEN 
C 
C  Initialize 
C 
          C_FACTOR = 1.0 
C 
          IVTAM = IV 
          IF (IV.EQ.VTGASBUS) IVTAM = 14 
C 
C  If reformulated or oxygenated fuel, compute oxy effect. 
C 
        IF(RFGON .OR. OXYFLG.NE.1) C_FACTOR = HCOXY(MYR,IV,EF) 
C 
C  Begin detergent gasoline effects calculations 
C    22-April-1994 @ CSC-tmm Request #2-421 
C  Modified to only calculate when the new DTGFLG is equal to 1 or 2 
C    7-May-1994 @ CSC-tmm Request #2-441 
C 
C  8-Feb-1996 @ DynTel Request 2-615 Disables detergent gasoline effects when 
C               NEWFLG is 5 or 6. 
C 
        IF (DTGFLG.GT.0 .AND. NEWFLG.LT.5)  
     &    C_FACTOR = C_FACTOR*DGADJF(IP,ICY,IFDS) 
C 
C  Apply RFG corrections to light duty gasoline vehicles 
C 
        IF(VVLDGAS(IV).EQ.1 .AND. ICY.GE.RFGSYR(1) 
     &                      .AND. RFGFLG.EQ.2)  
     &    C_FACTOR = C_FACTOR*RFG_CF(IP,RFG_REGION,JPHASE,IV) 
C 
C  Compute Sulfur effect 
C 
        IF(VVGAS(IV) .EQ. 1) THEN !No sulfur correction for motorcycles 
           F_CF(1) = C_FACTOR *      ! Sulfur CF running 
     &        SUL_COR(IP,ISTD,ICY,MYR,IV,INH,1) 
           IF(VVLDGAS(IV) .EQ. 1) THEN 
             F_CF(2) = C_FACTOR *    ! Sulfur CF start 
     &          SUL_COR(IP,ISTD,ICY,MYR,IV,INH,2) 
           ELSE 
             F_CF(2) = C_FACTOR 
           END IF 
        END IF 
C 
C     For LDGV/T & HDGV, the RVPCF applied may be 1.0, either because there 
C     is no correction or because  the correction is by bag (& has already 
C     been applied via OMTCF).  Only RVP MYG 1 uses RVPCF. Diesel classes 
C     and motorcycles have no RVPCF corrections. 
C 
        DO ISR = 1, 2 
          DO IH = 1, MAXIH 
            IF(IV .NE. VTMC) THEN 
              FUEL_CF(ISR,INH,ISTD,IFDS,IH,IVT) =  
     &          F_CF(ISR) * RVPCF(IDX,IP,IVTAM,IH) 
            ELSE 
              FUEL_CF(ISR,INH,ISTD,IFDS,IH,IVT) = C_FACTOR 
            END IF 
          END DO    ! IH loop 
        END DO      ! ISR loop 
C 
      ELSE          ! Diesel Vehicles have no fuel adjustment. 
        DO ISR = 1, 2 
          DO IH = 1, MAXIH 
            FUEL_CF(ISR,INH,ISTD,IFDS,IH,IVT) = 1.0 
          END DO 
        END DO 
      END IF 
C 
      RETURN 
      END 
 
