      SUBROUTINE
     *   FUEL(MYR,IP,IV,ICY,IFDS,INH,IDX,ISTD,HREF,FUEL_CF,INERR)
C
C  FUEL computes the detergent gas, sulfur, RVP and oxy fuel, effect of
C  industry average fuel and reformulated gas adjustments.  Due to the
C  introduction of oxidation and 3-way catalyst on LDTs, special
C  handling of those technologies was needed.
C
C  Called by HRLOOP.
C
C  Calls NMOGFID
C
C  Changes: (Last change first)
C
C  11 Jul 02 AIR Task 36: Hourly humidity makes some EF values hour dependent,
C            modified code to accept hourly EF values.
C  15 Apr 02  AIR Task Bug Fix 388: modification of constants via QUITER.
C  26 Jun 00 @EPA-elg  Changed IPHASE to JPHASE.
C  20 Jun 00 @EPA-ddj Bug127 Added conditional to bypass sulfur
C            correction if not LDGV.
C  15 May 00 @EPA-ddj BUG27, Added all fuel corrections.
C  02 May 00 @EPA-ddj E-7, Removed DATA statements and initialized
C            contants in type statements.
C  04 Apr 00 AIR Task 02: Removed IASTM reference in RFORM1 common.
C  11 Apr 00 @DynTel-BG 2-668 removed ISR and ISTD from argument list;
C            fixed CO 1981 model problem (replaced a .le. with a .lt.)
C  22 Mar 00 @EPA-eg           hrloop revisions
C  25 Aug 99 @DynTel-ddj 2-668 Simplified calculation of BENFAC to
C            reduce calculation time.  Added effect of non P3W
C            technology for post 1981 vehicles and effect of P3W for
C            post 1994 Tier0 vehicles.
C  12 Aug 99 @DynTel-HXQ 2-667 Moved the effect for HC to HCOXY.FOR
C   9 Aug 99 @dynTel-ddj 2-694 Added ISTD parameter to Fuel function
C            list and to call to SUL_COR.
C  13 Jul 99 @DynTel-ddj 2-699 Changed sulfur correction factor to a
C            call to new function SUL_COR.
C  09 Feb 99 @DynTel-RWK&HXQ 2-667 Made some requested cosemetic changes
C  11 Dec 98 @DynTel-HXQ 2-667  Fuel Oxygen Effects on Exchaust Co Emission
C  02 Oct 98 @DynTel-ZK 2-000 Explicit Typing
C  03 Sep 98 @DynTel-RJD 2-663
C  24 Aug 98 @DynTel-RJD 0-632 Added IFDS parameter to argument list.
C            Moved detergent gasoline adjustment to fuel and eliminated
C            some dead code.
C
C  Input on call:
C
C    argument list: MYR,IP,IV,ICY,IFDS,INH,IDX,ISTD,EF,FUEL_CF
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
C              F_CF
C
C  Local array subscripts:
C
C     COLDGV(10)  - COLDGV(EMY)
C     COLDGT(10)  - COLDGT(EMY)
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
C  BENFAC   R    benefit factor for determining oxy adjustment
C  CAT      R    fraction of vehicles that are catalyst equipped
C  FUEL1    R    1st oxy fuel adjusment factor: pre-81 LDV, HDGV, MC and
C                pre-91 LDTs
C  COLDGV   R    Hardcoded oxygen fuel correction factors for LDGV
C  COLDGT   R    Hardcoded oxygen fuel correction factors for LDGT
C  COHDGV   R    Hardcoded oxygen fuel correction factors for HDGV
C  COMC     R    Hardcoded oxygen fuel correction factors for MC
C  NOCAT    R    fraction of vehicles that are not catalyst equipped
C  OXYCT    R    local variable for the oxygen content
C  P3W      R    percent of vehicles with 3-way catalyst
C  POX      R    percent of vehicles with ox catalyst
C
C
C  Notes:
C
C  FUEL was added for MOBILE4.1 as OXYADJ and renamed for MOBILE5.
C  FUEL was MODIFIED for MOBILE5v03 to include the effect of industry average
C  fuel and reformulated gasoline fuel adjustments.
C
C  22-April-1994 @ CSC-tmm Request #2-421 Added code for detergent gasoline adjustment
C     Added ICY to argument list
C  7-May-1994 @ CSC-tmm Request #2-441 Modified code to permit the detergent gasoline
C     parameters to be read in and stored in COMMON.
C
C  29-July-1994 @ CSC-tmm Request #2-457 Replaced WINFLG with SEAFLG
C  18-August-1994 @ CSC-rwk for tmm Request #0-441, Added a modification for
C     July runs for detergent gas.  Also added a call to QUITER if the User
C     inputs calculated an unreasonable detergent gas factor.
C  24-August-1994 @CSC-pme request 2-446. The third dimension of OFFMTH array and
C                 ICAP argument were removed.
C  22-Sep-1994 @ CSC-bsg Request 426  Modularized the calculation of the
C  3-Way oxidation values, by calling a function to calculate them.  Updated
C  PHSFAC.
C  12-Dec-1994 @ CSC-ked Request 0-441 Detergent Gas Fix for July and smooth
C     error control.
C  08-Feb-1996 @ DynTel Request 2-615 Detergent gasoline effects are disabled
C                when NEWFLG is set to 5 or 6.
C  22-Jul-1996 @ DynTel-wcs Request 2-627 Alter PHSFAC values in order
C                meet VOC HC Exhaust target values for LDGV & LDGT
C  13-Sep-1996 @DynTel-yc Request #2-630 The detergent gas phase-in effect
C               was updated by refining phase-in period
C  11 Jun 97 DynTel@MLA 2-622  Added IVTYPE.I for Mobile6 and modified the
C            code to use the symbolic vehicle type names and vehicle type
C            vectors that were introduced in Mobile6.
C  26 Feb 98 @Dyntel-rjd (after gjr) Request 2-204 Removed IY from 
C            parameter list. Removed code that was used for July 
C            emissions that is no longer needed.  That is because July 
C            emissions will be calculated from extending the VMT mix by 
C            6 months rather than interpolating it from two consecutive 
C            Januaries.
C  21 Apr 98 @DynTel-RJD 2-654 Added the ISR variable to the input argument 
C            list. This is used in the OFFMTH array.    
C  23 Jul 98 @DynTel-MLA 2-000  Changed LDGT names to conform to
C            certification naming conventions.
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
      INTEGER, INTENT(IN)      :: IP
      INTEGER, INTENT(IN)      :: IV
      INTEGER, INTENT(IN)      :: INH
      INTEGER, INTENT(IN)      :: ICY
      INTEGER, INTENT(IN)      :: MYR
      INTEGER, INTENT(IN)      :: IFDS
      INTEGER, INTENT(IN)      :: IDX
      INTEGER, INTENT(IN)      :: ISTD
      INTEGER, INTENT(INOUT)   :: INERR
C
      REAL, DIMENSION(24,2), INTENT(IN)           :: HREF
      REAL, DIMENSION(2,2,10,3,24), INTENT(INOUT) :: FUEL_CF
C
      REAL DGADJF
      REAL C_FACTOR
      REAL, DIMENSION(2,2,10,3)      :: F_CF
      REAL, DIMENSION(2)             :: EF
C
      INTEGER IH
      INTEGER ISR
      INTEGER IVTAM
C
C    External Functions
C
      REAL COOXY
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
          IVTAM = IV
          IF (IV.EQ.VTGASBUS) THEN
            IVTAM = 14
          ENDIF
C
C  If reformulated or oxygenated fuel effects.
C
        IF(RFGON .OR. OXYFLG.NE.1) THEN
C
C  Compute Oxy effect
C  (Note: Due to hourly humidity CF, only NOx basic emission factors are
C  actually hour dependent. Therefore, for HC, put first hour's HREF values
C  into the local EF array so that code modifications are kept to a minimum.)
C
          IF(IP .EQ. 1) THEN
             EF(1)=HREF(1,1)  ! HC basic emission factors are hour-independent
             EF(2)=HREF(1,2)  ! HC basic emission factors are hour-independent
             C_FACTOR = HCOXY(MYR,IV,EF)
          ELSEIF(IP .EQ. 2)THEN
             C_FACTOR = COOXY(MYR,IV,INH)
          END IF
        END IF
C
C  Begin detergent gasoline effects calculations
C    22-April-1994 @ CSC-tmm Request #2-421
C  Modified to only calculate when the new DTGFLG is equal to 1 or 2
C    7-May-1994 @ CSC-tmm Request #2-441
C
C  8-Feb-1996 @ DynTel Request 2-615 Disables detergent gasoline effects when
C               NEWFLG is 5 or 6.
C
        IF (DTGFLG.GT.0 .AND. NEWFLG.LT.5) THEN
          C_FACTOR = C_FACTOR*DGADJF(IP,ICY,IFDS,INERR)
        END IF
C
C  Apply RFG corrections to light duty gasoline vehicles
C
        IF(VVLDGAS(IV).EQ.1 .AND. ICY.GE.RFGSYR(1)
     &                      .AND. RFGFLG.EQ.2) THEN
           C_FACTOR = C_FACTOR*RFG_CF(IP,RFG_REGION,JPHASE,IV)
        END IF
C
C  Compute Sulfur effect
C
        IF(VVGAS(IV) .EQ. 1) THEN !No sulfur correction for motorcycles
           F_CF(1,INH,ISTD,IFDS) = C_FACTOR *      ! Sulfur CF running
     &        SUL_COR(IP,ISTD,ICY,MYR,IV,INH,1)
           IF(VVLDGAS(IV) .EQ. 1) THEN
             F_CF(2,INH,ISTD,IFDS) = C_FACTOR *    ! Sulfur CF start
     &          SUL_COR(IP,ISTD,ICY,MYR,IV,INH,2)
           ELSE
             F_CF(2,INH,ISTD,IFDS) = C_FACTOR
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
              FUEL_CF(ISR,INH,ISTD,IFDS,IH) = F_CF(ISR,INH,ISTD,IFDS)*
     &                                        RVPCF(IDX,IP,IVTAM,IH)
            ELSE
              FUEL_CF(ISR,INH,ISTD,IFDS,IH) = C_FACTOR
            END IF
          END DO    ! IH loop
        END DO      ! ISR loop
C
      ELSE          ! Diesel Vehicles have no fuel adjustment.
        DO ISR = 1, 2
          DO IH = 1, MAXIH
            FUEL_CF(ISR,INH,ISTD,IFDS,IH) = 1.0
          END DO
        END DO
      END IF
C
      RETURN
      END

