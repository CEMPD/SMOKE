      REAL FUNCTION EFIVT(IV,IEF,IP) 
C 
C  EFIVT returns emission factors for individual vehicle  
C  types.        
C 
C  Called by EFCVT 
C 
C  Changes: (Last change first) 
C 
C  10 May 01 AIR Task 11/5: IDLE PM 
C  19 Feb 01 AIR Task 11: Expanded to include particules, IEFPM 
C  03 Nov 99 @DynTel-bag 1-006 Removed partial diurnals; renamed refueling 
C            and resting loss emissions; added G_MI_EVAP array for the 
C            evaporative emissions with units of grams/mile 
C  19 Mar 99 @DynTel-MLA 2-701 Added cases for start (ISTARTEF), running 
C            (IRUNEF), and total evaporative (ITEVAP) emissions. 
C  22 Sep 98 @DynTel-ZK  2-000 Explicit Typing 
C 
C  Input on call or from calls: 
C 
C    parameter list: IV,IEF,IP 
C    common blocks: 
C    /EVAPGR/ GREVP,IGREVP 
C    /RESUL1/ EFFTP,EFEXH,EFEVAP,REFL_GM_MI,EFRUNL, 
C             IEFFTP,IEFEXH,IEFEVAP,IEFREFL,IEFRUNL 
C    /RESUL2/ EFIDLE,G_MI_EVP,IRLGGAL,IRSTGPH,IEFIDLE,IGMIEVP 
C    /RESUL4/ IRUNEF, ISTARTEF, ITEVAP, RUNEF, STARTEF, TOTEVAP 
C 
C  Notes. 
C   
C  10 Mar 98 @DynTel-RJD 2-622 Created this function for Mobile6 
C            to assist EFCVT in the calculation of emission factors  
C            for composite vehicle types.  
C 
      IMPLICIT NONE  
      INCLUDE 'EVAPGR.I' 
      INCLUDE 'PART3.I' 
      INCLUDE 'RESUL1.I' 
      INCLUDE 'RESUL2.I' 
      INCLUDE 'RESUL4.I' 
C 
      INTEGER, INTENT(IN) :: IV 
      INTEGER, INTENT(IN) :: IEF 
      INTEGER, INTENT(IN) :: IP 
C 
      IF(IEF.EQ.IEFFTP) THEN 
           EFIVT=EFFTP(IP,IV) 
      ELSE IF (IEF.EQ.IEFEXH) THEN 
           EFIVT=EFEXH(IV) 
      ELSE IF (IEF.EQ.IEFEVAP) THEN  
           EFIVT=EFEVAP(IV) 
      ELSE IF (IEF.EQ.IEFREFL) THEN 
           EFIVT=REFL_GM_MI(IV) 
      ELSE IF (IEF.EQ.IEFRUNL) THEN 
           EFIVT=EFRUNL(IV) 
      ELSE IF (IEF.EQ.IEFIDLE) THEN  
           EFIVT=EFIDLE(IP,IV) 
      ELSE IF (IEF.EQ.IRLGGAL) THEN 
           EFIVT=RLGGAL(IV) 
      ELSE IF (IEF.EQ.IRSTGPH) THEN  
           EFIVT=RSTGPH(IV) 
      ELSE IF (IEF.EQ.IGREVP) THEN  
           EFIVT=GREVP(IP,IV) 
      ELSE IF (IEF.EQ.IRUNEF) THEN 
           EFIVT=RUNEF(IP,IV) 
      ELSE IF (IEF.EQ.ISTARTEF) THEN 
           EFIVT=STARTEF(IP,IV) 
      ELSE IF (IEF.EQ.ITEVAP) THEN 
           EFIVT=TOTEVAP(IV) 
      ELSE IF (IEF.EQ.IGREVP) THEN 
           EFIVT=GREVP(IP,IV)     ! IP here is the evaporative emissions type/index 
      ELSE IF (IEF.EQ.IGMIEVP) THEN 
           EFIVT=G_MI_EVP(IP,IV)  ! IP here is the evaporative emissions type/index 
                                  ! IHOTSOAK=1 
                                  ! IDIURNAL=2 
                                  ! IRUNLOSS=3 
                                  ! ICRNCASE=4 
                                  ! IRESTLOSS=5 
                                  ! IREFUEL=6 
      ELSE IF (IEF.EQ.IEFPM) THEN 
           EFIVT=EFPM(IP,IV)      ! IP here is the PM emissions type/index 
                                  ! SO4=1 
                                  ! OCARBON=2 
                                  ! ECARBON=3 
                                  ! GASPM=4 
                                  ! LEAD=5 
                                  ! SO2=6 
                                  ! NH3=7 
                                  ! BRAKE=8 
                                  ! TIRE=9 
                                  ! EMFACP_PM=10 
                                  ! IDLE PM=11 
      ENDIF 
C 
      RETURN 
      END 
