      SUBROUTINE PMBERPAR(ICY,INERR) 
C 
C  This subroutine determines the associated PM BER sulfur levels and the 
C  the maximum BER allowed if various regulations are enforced. The maximum 
C  BERs are at PM30. 
C 
C  Changes : 
C 
C  08 Mar 02 AIR Task 25: Fixed TIER2 LDDV/LDDT Bug 
C  20 Dec 01 AIR Task 25: The associated diesel sulfur level is now by-model-year. 
C  04 Dec 01 AIR Task 27-28: Removed 15 ppm diesel sulfur limit from  
C            HDV 2007 Rule. Did not remove warning message, however. 
C  20 Oct 01 AIR Task 25: New M6.2T subroutine 
C 
C     Called by:  PMCALX 
C 
C     Calls: NONE 
C 
C  Input on call: 
C 
C    Parameter list: ICY 
C 
C    Common blocks:  See Includes below 
C 
C 
      IMPLICIT NONE 
C 
      INCLUDE 'FLAGS5.I'   ! SOAK_FLG,TIER2_FLAG,HDR2007 
      INCLUDE 'IVTYPE.I'   ! MAXVEH, VVDSL 
      INCLUDE 'MAXIMA.I'   ! MAXYRS 
      INCLUDE 'MYRCAL.I'   ! TF 
      INCLUDE 'PART1.I'    ! PRTCHK, D_PPM, BERSULF, BERMAX 
      INCLUDE 'PART4.I'    ! BSULFLVL 
      INCLUDE 'SULFUR.I'   ! S_PPM 
      INCLUDE 'TIER2PM.I'  ! T2PM_50K, T2PM_FRAC 
      INCLUDE 'VMXCOM.I'   ! VMTMIX 
C 
C  Declare external functions. 
C 
      REAL   , EXTERNAL :: HDCPIC 
C 
C  Declare parameter list 
C 
      INTEGER ICY     !  Calendar Year 
      INTEGER INERR   !  Errors 
C 
C  Declare local variables/arrays 
C 
      INTEGER IBIN    ! BIN number 
      INTEGER IDX     ! IDX=(MAXYRS+1)-JDX 
      INTEGER IMY     ! Temporary model year 
      INTEGER IV      ! Vehicle Class 
      INTEGER IVTEMP  ! temporary IV 
      INTEGER JDX     ! Age 
      INTEGER MY      ! Model year 
C 
C 
C  Initialize the BER associated suflur levels. 
C 
      BERSULF = 0. 
C 
C  Initialize the BER limit. The negative value means no limit will be applied. 
C 
      BERMAX = -1. 
C 
C  If HDD2007 rule is in effect and the calendar year diesel 
C  sulfur is above 15 ppm, warn user. Do Not reset the sulfur level. 
C  This warning is issued because the diesel ECARBON and OCARBON 
C  values may become zero. 
C 
      IF(HDR2007.EQ.1.AND.ICY.GE.2007.AND.D_PPM.GT.15.) 
     *   CALL QUITER(D_PPM,0,111,INERR) 
C 
C  Loop for each vehicle class 
C 
      DO IV=1,MAXVEH 
C 
C  Skip if VMTMIX is 0. 
C 
      IF(VMTMIX(IV).EQ.0.) CYCLE 
C 
C  Loop for all ages 
C 
      DO JDX=1,MAXYRS 
C 
C  Set IDX 
C 
      IDX=(MAXYRS+1)-JDX 
C 
C  If the travel fraction for this age/vehicle combination is zero, skip. 
C 
      IF(TF(IDX,IV) .EQ. 0.) CYCLE 
C 
C  Set model year 
C 
      MY=ICY-JDX+1 
C 
C  Determine the associated fuel sulfur content of the BER. (For gas, the  
C  BERs are based on the GAS SULFUR command for CY1999- or MOBILE6 CY2000+. 
C  For diesel, the BERs are based on the sulfur level input via PMDZML.CSV.) 
C  Store the sulfur values for use in determining the amount of associated 
C  SO4 to remove later. 
C 
      IF(VVDSL(IV).EQ.1) THEN 
        IMY=MY-1950+1 
        IF(IMY.LT.1) IMY=1 
        IF(IMY.GT.71) IMY=71 
        BERSULF(JDX,IV)=BSULFLVL(IMY) 
      ELSE 
        BERSULF(JDX,IV)=S_PPM 
      ENDIF 
C 
C  Gas BINs are applicable to Tier2 LDGV/T MY2004+ when the TIER2 flag is on. 
C  Diesel BINs are applicable to LDDV/T MY2004+ when the TIER2 flag is on. 
C  Note: The associated fuel sulfur content for TIER2 BERs is 30 ppm for 
C  gas and 8 ppm diesel. 
C 
      IF(TIER2_FLAG .EQ. 1         .AND.       ! TIER 2 LDGT/LDGT and LDDV/LDDT applicability 
     *           MY .GE. 2004      .AND. 
     *          (IV .LE. VTLDGT4   .OR. 
     *  VVLDDSL(IV) .EQ. 1)) THEN 
C 
        IVTEMP=IV                               ! Map IV into BINs 
        IF(IV.EQ.VTLDDV .OR. 
     *     IV.EQ.VTLDDT12) IVTEMP=IV-8 
        IF(IV.EQ.VTLDDT34) IVTEMP=8 
C 
        IMY=MY-2004+1                           ! Map MY into BINs 
        IF(IMY.GT.12) IMY=12 
C 
        BERMAX(JDX,IV)=0. 
C 
        DO IBIN = 1, 12                         ! Compute BIN Limit 
          BERMAX(JDX,IV) = BERMAX(JDX,IV)  
     *                   + T2PM_50K(IVTEMP,IBIN) 
     *                   * T2PM_FRAC(IMY,IBIN,IVTEMP) 
        END DO 
C 
          IF(VVDSL(IV).EQ.1) THEN               ! Store associated sulfur level 
            BERSULF(JDX,IV)=8. 
          ELSE 
            BERSULF(JDX,IV)=30. 
          ENDIF 
C 
      END IF 
C 
C  Determine if 2007 HDV Rule is in effect. If so, set the limit to 
C  0.0092 g/bhp-hr (start + running).  Apply same limit to HDGV.  
C  Note: The associated fuel sulfur content for 2007 HDV BERs is 30 ppm for 
C  gas and 8 ppm diesel. 
C  (Programming Note: currently the deterioration rate for the 0.0092 
C   g/bhp-hr standard is set to 0.0.) 
C 
      IF(HDR2007.EQ.1.AND.MY.GE.2007) THEN             ! HDV 2007 applicability 
C 
      IF(VVHDDV(IV).EQ.1 .OR. VVHDGV(IV).EQ.1) THEN 
C 
        BERMAX(JDX,IV) = 0.0092 * HDCPIC(MY,IV)        ! Compute 2007 HDV BER limit in g/mi, at PM30. 
C 
          IF(VVDSL(IV).EQ.1) THEN                      ! Store associated sulfur level 
            BERSULF(JDX,IV)=8. 
          ELSE 
            BERSULF(JDX,IV)=30. 
          ENDIF 
C 
      END IF 
      END IF 
C 
      END DO ! JDX 
C 
      END DO ! IV 
C 
      RETURN 
      END 
       
