      REAL FUNCTION AFTRPREF(MY,ACTUAL_MY,IDX,ICY,JDX,IP,IV,ISR,IFDS, 
     *                       ISTD,IH,IROAD,CAT,VMTAGE,IMPGM) 
C 
C    This module determines the emission rates for I&M saw tooth 
C 
C 
C  Called by IMPROG. 
C 
C  Calls EF_CASE, ARBEF, IMCF 
C 
C  17 Dec 01 @EPA-BAG        Remove I&M sawtooth 
C  12 Dec 01 @EPA-BAG bug    Use minimum value of ARBEF and EF_NRML0 
C                            rather than the maximum of the two values 
C  04 Dec 01 @EPA-BAG bug    Use maximmum value of ARBEF and EF_NRML0 
C                            rather than the minimum of the two values 
C  17 Oct 01 @EPA-BAG bug331 Biennial/ annual I&M program discrepancy 
C 
C    argument list: MY,ACTUAL_MY,IDX,ICY,IP,IV,ISR,IFDS,ISTD,IH,IROAD, 
C                   CAT,VMTAGE,VMTAGE2M,VMTAGE1M,VMTAGE1P 
C 
C    common blocks: 
C    /IMPAR6/ IFREQ 
C    /LDGBER/ CRT050K 
C    /TIER2/ T2CRT50K 
C 
C  Output on return: 
C 
C    modules: 
C 
C    function : 
C 
C    common blocks: 
C    /LDGOBD/ AR94B 
C 
C 
C  Local variable / array dictionary: 
C 
C   Name        Type              Description 
C  ------       ----  ----------------------------------------------------- 
C   MEF          R     Temporary storage for variable AR94B 
C   NORMAL       I     Constant(=1) indicating normal emission calculations 
C   HIGH         I     Constant(=2) indicating high emission calculations 
C 
C 
      IMPLICIT NONE 
      INCLUDE 'IMPAR6.I' 
      INCLUDE 'LDGBER.I' 
      INCLUDE 'LDGOBD.I' 
      INCLUDE 'TIER2.I' 
C 
      INTEGER, INTENT(IN) :: MY 
      INTEGER, INTENT(IN) :: ACTUAL_MY 
      INTEGER, INTENT(IN) :: IDX 
      INTEGER, INTENT(IN) :: ICY 
      INTEGER, INTENT(IN) :: JDX 
      INTEGER, INTENT(IN) :: IH 
      INTEGER, INTENT(IN) :: IROAD 
      INTEGER, INTENT(IN) :: IP 
      INTEGER, INTENT(IN) :: IV 
      INTEGER, INTENT(IN) :: IFDS 
      INTEGER, INTENT(IN) :: ISR 
      INTEGER, INTENT(IN) :: ISTD 
      INTEGER, INTENT(IN) :: IMPGM 
      REAL,    INTENT(IN) :: CAT 
      REAL,    INTENT(IN) :: VMTAGE 
C 
C     Declare Local Variables 
C 
      REAL EF_AR 
      REAL EF_NRML 
      REAL MEF 
      INTEGER, SAVE :: NORMAL=1 
C 
C  Declare external functions 
C 
      REAL EF_CASE 
      REAL ARBEF 
C 
      AFTRPREF =0.0 
      AR94B   =0.0 
      EF_AR   =0.0 
      MEF     =0.0 
C 
C 
C     REPAIRED EMISSIONS RATES 
C 
C     Model Years Less Than 1994 
      SELECT CASE (MY) 
      CASE(:1993) 
        EF_AR = ARBEF(MY,JDX,IP,IV,VMTAGE,IFDS,ISR,IMPGM,1) 
        CALL IMCF(ICY,IV,IP,ISR,IDX,MY,ACTUAL_MY,ISTD,IFDS,NORMAL,IH, 
     *            IROAD,CAT,EF_AR) 
C 
C     Determine the Post '94 Repaired Emission Rates 
C     Model Years Greater Than 1993 
      CASE(1994:) 
C 
        IF((MY==1994.OR.MY==1995).AND.ISTD==1) THEN 
          EF_AR = ARBEF(MY,JDX,IP,IV,VMTAGE,IFDS,ISR,IMPGM,1) 
          CALL IMCF(ICY,IV,IP,ISR,IDX,MY,ACTUAL_MY,ISTD,IFDS,NORMAL,IH, 
     *              IROAD,CAT,EF_AR) 
C 
        ELSE 
          IF(ISTD.NE.4) THEN 
            AR94B   = OBDREP*CRT050K(ISTD,IP,IV)*R_S_CONV(ISTD,IP,ISR) 
          ELSEIF(ISTD.EQ.4) THEN      !If Tier2 use T2CRT50K 
            AR94B   = OBDREP*T2CRT50K(IV,IP,JDX)*R_S_CONV(ISTD,IP,ISR) 
          ENDIF 
C 
C         After repair emissions are either Normals or capped at 1.5 times 
C         Cert Standards. 
C 
C         Mileage at Current year : 
          MEF = AR94B 
          CALL IMCF(ICY,IV,IP,ISR,IDX,MY,ACTUAL_MY,ISTD,IFDS,NORMAL,IH, 
     *              IROAD,CAT,MEF) 
          AR94B = MEF 
C         NORMAL Emissions Rates 
          EF_NRML = EF_CASE(MY,JDX,IP,IV,ISR,VMTAGE,ISTD,NORMAL,IFDS) 
          CALL IMCF(ICY,IV,IP,ISR,IDX,MY,ACTUAL_MY,ISTD,IFDS,NORMAL,IH, 
     *              IROAD,CAT,EF_NRML) 
C 
          AR94B=AMIN1(AR94B,EF_NRML) 
          EF_AR = AR94B 
        ENDIF 
C 
      END SELECT 
C 
      AFTRPREF=EF_AR 
C 
      RETURN 
      END 
