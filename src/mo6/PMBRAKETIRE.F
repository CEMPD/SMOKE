      SUBROUTINE PMBRAKETIRE 
C 
C  Computes the PART5-based BRAKE and TIRE WEAR emission factors. 
C 
C  Changes : 
C 
C  11 Feb 01 AIR Task 11: New Routine 
C 
C  Called by: PMCALX 
C 
C  Calls: None 
C 
C  Input on call: 
C 
C    Common blocks: 
C    /MAXIMA/     MAXYRS,MAXIH 
C     IVTYPE.I    MAXVEH 
C    /PSIZE/      PSTIRE,PSBRK 
C    /MYRCAL/     TF 
C    /PART1/      IVEHWL 
C    /SPEED9/     FVMT,HVMT  
C    /VMXCOM/     VMTMIX     
C 
C  Output on return: 
C 
C    Common blocks: 
C    /PART3.I/    EFPM,VPM 
C 
C  Local variable / array dictionary: 
C 
C   Name        Type              Description 
C  ------       ----  ----------------------------------------------------- 
C  CBRAKE        R    Temporary brake wear calculation 
C  CTIRE         R    Temporary tire wear calculation 
C  IDX           I    IDX=MAXYRS+1-JDX 
C  IH            I    Hour 
C  IROAD         I    Road tyep 
C  IV            I    Vehicle class 
C  JDX           I    Age 
C  SUMBIH        R    Temporary summation for brake 
C  SUMBJDX       R    Temporary summation for brake 
C  SUMTIH        R    Temporary summation for tire 
C  SUMTJDX       R    Temporary summation for tire 
C 
C 
      IMPLICIT NONE 
C 
      INCLUDE 'IVTYPE.I'  !MAXVEH 
      INCLUDE 'MAXIMA.I' 
      INCLUDE 'MYRCAL.I' 
      INCLUDE 'PART1.I' 
      INCLUDE 'PART3.I' 
      INCLUDE 'PSIZE.I' 
      INCLUDE 'SPEED9.I' 
      INCLUDE 'VMXCOM.I' 
C 
      INTEGER IDX 
      INTEGER IH 
      INTEGER IROAD 
      INTEGER IV 
      INTEGER JDX 
      REAL    CBRAKE 
      REAL    CTIRE 
      REAL    SUMBIH 
      REAL    SUMBJDX 
      REAL    SUMTIH 
      REAL    SUMTJDX 
C 
C 
C  Calculate the BRAKE emission factor, in g/mi. BRAKE emissions are 
C  dependent only on particle size. 
C 
      CBRAKE = 0.0128 * PSBRK 
C 
C  Loop for each vehcile class 
C 
      DO IV=1,MAXVEH 
C 
C  Skip if VMTMIX is zero 
C 
      IF(VMTMIX(IV).EQ.0.) CYCLE 
C 
C  Calculate the TIRE WEAR emission factor, in g/mi. TIRE WEAR emmisions 
C  only depend on particle size and number of wheels on the vehicle. 
C 
       CTIRE = 0.002 * PSTIRE * FLOAT(IVEHWL(IV)) 
C 
C  Loop for each age 
C 
       DO JDX=1,MAXYRS 
C 
C  Set IDX 
C 
        IDX=MAXYRS+1-JDX 
C 
C  If travel fraction is 0, skip 
C 
        IF(TF(IDX,IV).EQ.0) CYCLE 
C 
C  Initialize 
C 
        SUMBJDX=0.0 
        SUMTJDX=0.0 
C 
C  Load arrays for Hour and Road. Currently, the EFs are independent 
C  of Road Type. Also compute by-model-year and composite emission values. 
C 
        DO IH=1,MAXIH 
C 
         SUMBIH=0. 
         SUMTIH=0. 
C 
         DO IROAD=1,4 
C 
          DBPM(8,IROAD,IH,JDX,IV) = CBRAKE 
          DBPM(9,IROAD,IH,JDX,IV) = CTIRE    
C 
          SUMBIH = SUMBIH + DBPM(8,IROAD,IH,JDX,IV) * FVMT(IROAD,IH,IV)  !BRAKE 
          SUMTIH = SUMTIH + DBPM(9,IROAD,IH,JDX,IV) * FVMT(IROAD,IH,IV)  !TIRE 
C 
         END DO !IROAD 
C 
         SUMBJDX = SUMBJDX + SUMBIH * HVMT(IH)  !BRAKE 
         SUMTJDX = SUMTJDX + SUMTIH * HVMT(IH)  !TIRE 
C 
        END DO !IH 
C 
        EFPM(8,IV) = EFPM(8,IV) + SUMBJDX * TF(IDX,IV)  !BRAKE 
        EFPM(9,IV) = EFPM(9,IV) + SUMTJDX * TF(IDX,IV)  !TIRE 
C 
       END DO !JDX 
C 
       VPM(8) = VPM(8) + EFPM(8,IV) * VMTMIX(IV)  !BRAKE 
       VPM(9) = VPM(9) + EFPM(9,IV) * VMTMIX(IV)  !TIRE 
C 
      END DO !IV 
C 
      RETURN 
      END 
