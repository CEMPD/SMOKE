      SUBROUTINE CLCIMPAR(MY,JDX,IP,IV,ISR,IFDS,ISTD,IMPGM,VMTAGE, 
     *         PRCNT_RPRD,PRCNT_NTIDD,PRCNT_NCOMP,PRCNT_WVRS,IDR,INERR) 
C 
C  This module determines the vehicle fractions used for I&M saw tooth 
C 
C  Called by IMPROG. 
C 
C  Calls ARBEF, IDRATE 
C 
C  17 Dec 01 @EPA-BAG        Remove I&M sawtooth 
C  04 Dec 01 @EPA-BAG        Removed unused variable, AFTR_RPR_HI1P, and associated 
C                            1P variables 
C  30 Nov 01 @EPA-BAG bug331 added checks for ratios larger than 1. and 
C            less than 0.; fixed mileages, VMTAGE1M and VMTAGE2M for JDX=2 
C  17 Oct 01 @EPA-BAG bug331 Biennial/ annual I&M program discrepancy 
C 
C    argument list: MY,JDX,IP,IV,ISR,IFDS,ISTD,IMPGM,VMTAGE,PRCNT_RPRD, 
C                   PRCNT_NTIDD,PRCNT_NCOMP,PRCNT_WVRS,IDR,INERR 
C 
C    modules: 
C 
C    common blocks: 
C    /IMPAR1/ CRIM,WAV81PLUS 
C    /IMPAR5/ DISCNT 
C    /IMPAR6/ INTYP 
C    /LDGIMC/ WAVRDC 
C    /LDGOBD/ OBDHI_IM, OBDHINIM 
C 
C  Output on return: 
C 
C    function: IMBENEFIT 
C    common blocks: 
C 
C  Local variable / array dictionary: 
C 
C   Name         Type              Description 
C  ------        ----  ----------------------------------------------------- 
C  IDR            R    I&M identification fraction of vehicles for model 
C                      year under consideration 
C  NONCOM         R    fraction of non-compliant vehicles 
C  WAIVER         R    fraction of waivered vehicles determined from 
C                      input/default data 
C  BETA           R    variable used for model year dependence OBD waiver and 
C                      non-compliant fractions (current model year) 
C  WAVR_NCM       R    fraction of waivered and non-compliant vehicles 
C 
      IMPLICIT NONE 
      INCLUDE 'IMPAR1.I' 
      INCLUDE 'IMPAR5.I' 
      INCLUDE 'IMPAR6.I' 
      INCLUDE 'LDGIMC.I' 
      INCLUDE 'LDGOBD.I' 
C 
      INTEGER, INTENT(IN)    :: MY 
      INTEGER, INTENT(IN)    :: JDX 
      INTEGER, INTENT(IN)    :: IP 
      INTEGER, INTENT(IN)    :: IV 
      INTEGER, INTENT(IN)    :: IFDS 
      INTEGER, INTENT(IN)    :: ISR 
      INTEGER, INTENT(IN)    :: ISTD 
      INTEGER, INTENT(IN)    :: IMPGM 
      INTEGER, INTENT(INOUT) :: INERR 
      REAL,    INTENT(IN)    :: VMTAGE 
      REAL,    INTENT(OUT)   :: PRCNT_RPRD 
      REAL,    INTENT(OUT)   :: PRCNT_NTIDD 
      REAL,    INTENT(OUT)   :: PRCNT_NCOMP 
      REAL,    INTENT(OUT)   :: PRCNT_WVRS 
      REAL,    INTENT(OUT)   :: IDR 
C 
C  Declare Local Variables 
C 
      REAL BETA 
      REAL WAVR_NCM 
      REAL NONCOM 
      REAL WAIVER 
C 
C  Declare external functions 
C 
      REAL IDRATE 
C 
C     Initializations 
      WAVR_NCM   = 0.0 
      WAIVER     = 0.0 
      NONCOM     = 0.0 
      IDR        = 0.0 
      BETA       = 1.0 
      PRCNT_RPRD = 0.0 
      PRCNT_NTIDD= 0.0 
      PRCNT_NCOMP= 0.0 
      PRCNT_WVRS = 0.0 
C 
C     Assign waiver rates 
C     Control WAIVER (waver rates) for INTYP (I/M program inspection 
C     type), where 2 = I/M test and repair (Computerized), 
C     & 3 = test and repair (Manual), if INTYP is not a 1 = test 
C     only I/M and the WAIVER is LE 50%, WAIVER is set to zero. 
C 
      WAIVER=WAV81PLUS(IMPGM) 
      IF (WAIVER<=0.50.AND.INTYP(IMPGM)>1) THEN 
        IF (DISCNT(IP, INTYP(IMPGM)).LT.0.50) THEN 
          WAIVER = 0.0 
        ELSE 
          WAIVER=WAIVER*(2*DISCNT(IP,INTYP(IMPGM))-1) 
        ENDIF 
      ENDIF 
C 
C   determine the fraction of vehicles which are I&M non-compliant 
C 
      NONCOM = 1.0 - (CRIM(IMPGM)/100.0) 
C 
C   check for OBD 
C 
      WAVR_NCM=(1.0-WAIVER)*(1.0-NONCOM) 
      IF(WAVR_NCM<0.) WAVR_NCM=0. 
      IF(WAVR_NCM>1.) WAVR_NCM=1. 
C 
      IF(MY>=1996.AND.VMTAGE<=3.6) THEN 
        WAVR_NCM=1.0 
        BETA    =0.0 
      ENDIF 
C 
C  Determine I/M  program identification fractions. 
C 
      IF(MY<=1995 ) THEN 
        IDR =IDRATE(MY,JDX,IP,IV,ISR,IFDS,ISTD,IMPGM,VMTAGE,INERR) 
      ELSEIF(MY >= 1996) THEN 
        IF(OBDHINIM(JDX,IP,IV) < 0.000001 ) THEN 
          IDR  = 0.0 
        ELSE 
          IDR = 1.0 - (OBDHI_IM(JDX,IP,IV)/OBDHINIM(JDX,IP,IV))  ! 
          IF(IDR<0.) IDR=0. 
          IF(IDR>1.) IDR=1. 
        ENDIF 
      ENDIF 
C 
      PRCNT_RPRD =IDR*WAVR_NCM 
      PRCNT_NTIDD=(1.-NONCOM*BETA)*(1.-IDR) 
      PRCNT_NCOMP=NONCOM*BETA 
      PRCNT_WVRS =IDR*WAIVER*(1.-NONCOM)*BETA 
C 
      RETURN 
      END SUBROUTINE 
