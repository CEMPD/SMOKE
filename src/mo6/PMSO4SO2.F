      SUBROUTINE PMSO4SO2(ICY) 
C 
C  Computes the final, in-use fuel adjusted SO4 and SO2 
C  values for gas, diesel and NGV.  SO2 is gaseous, 
C  so no particle size adjustment is necessary. SO4 
C  is assumed to be PM2.5 or less, so no particle 
C  size adjustment can be applied. 
C  The original SO4 and SO2 values which come into this 
C  routine are based on GSO4SO2 and DSO4SO2. They are 
C  based on BER associated sulfur levels, hence the need 
C  to adjust them to in-use fuel levels. 
C 
C 
C  Changes : 
C 
C  03 May 01 AIR Task 12: Apply BER associated sulfur level. 
C  17 Apr 01 AIR Task 11: Changed code to handle by model year 
C            sulfur correction for SO4. 
C  02 Mar 01 AIR Task 11: New Routine 
C 
C  Input on call: 
C 
C    Constants:                None 
C    argument list:            ICY 
C    common blocks:            (See include files below) 
C 
C  Output on return:           DBPM, EFPM, VPM 
C 
C  common blocks:              (See include files below) 
C 
C  Local variable / array subscript/dictionary: 
C 
C   SLFCF(2)  --  SLFCF(IGD) 
C 
C   Name     Type             Description 
C  ------    ----  ----------------------------------------------------- 
C  CADJ        R   Combined NGV and ZEV sulfur correction factor 
C  IDX         I   Maxyrs - Age 
C  IH          I   Hour index 
C  IROAD       I   Road type index 
C  IV          I   Vehicle Loop 
C  JDX         I   Age    
C  MY          I   Model year 
C  NGVFRAC     R   Fleet NGV fraction (after ZEVs are removed) 
C  NGVPT       I   NGV model year index 
C  SULFCF      R   Sulfur correction factors 
C  SLFNGV      R   NGV sulfur correction factor 
C  SUMSO2IH    R   Summation of SO2 over IH,IROAD 
C  SUMSO2JDX   R   Summation of SO2 over IH,IROAD,JDX 
C  SUMSO4IH    R   Summation of SO4 over IH,IROAD 
C  SUMSO4JDX   R   Summation of SO4 over IH,IROAD,JDX 
C  ZEVFRAC     R   Fleet ZEV fraction 
C        
C 
      IMPLICIT NONE 
C 
      INCLUDE 'IVTYPE.I' ! /VVTYPS / VVDSL 
      INCLUDE 'MAXIMA.I' ! /MAXIMA / MAXYRS,MAXIH 
      INCLUDE 'MYRCAL.I' ! /MYRCAL / TF 
      INCLUDE 'NGV.I'    ! /NGVDAT1/ NGVFLG1,NGV_START,NGVPHIN,NGVFRAC 
      INCLUDE 'PART1.I'  ! /PART1  / D_PPM, BERSULF 
      INCLUDE 'PART3.I'  ! /PART3  / DBPM,EFPM,VPM 
      INCLUDE 'SPEED9.I' ! /SPEED9 / FVMT,HVMT 
      INCLUDE 'SULFUR.I' ! /S_DATA / S_PPM 
      INCLUDE 'VMXCOM.I' ! /VMXCOM / VMTMIX 
C 
C  Declare parameter list 
C 
      INTEGER ICY 
C 
C  Declare external functions. 
C 
      REAL, EXTERNAL :: PCTLEV 
C 
C  Declare local variables/arrays 
C 
      INTEGER IDX 
      INTEGER IH 
      INTEGER IROAD 
      INTEGER IV 
      INTEGER JDX 
      INTEGER MY 
      INTEGER NGVPT 
      REAL    CADJ 
      REAL    SULFCF 
      REAL    SLFNGV 
      REAL    SUMSO2IH 
      REAL    SUMSO2JDX 
      REAL    SUMSO4IH 
      REAL    SUMSO4JDX 
      REAL    ZEVFRAC 
C 
C 
C  Reinitialize to hold final adjusted values 
C 
       VPM(1) = 0.  !SO4 
       VPM(6) = 0.  !SO2 
C 
C  Compute the sulfur correction for NGV. Assume the NG sulfur 
C  level to be 0.01 ppm and that the corresponding gas sulfur 
C  level is 50 ppm. 
C 
      SLFNGV = 0.01 / 50.  !NGV 
C 
C  Loop for all vehicles 
C 
      DO IV=1,MAXVEH 
C 
C  Reinitialize to hold final adjusted values 
C 
       EFPM(1,IV) = 0.  !SO4 
       EFPM(6,IV) = 0.  !SO2 
C 
C  Skip if VMTMIX is zero 
C 
      IF(VMTMIX(IV).EQ.0.) CYCLE 
C 
C  Loop for all ages 
C 
       DO JDX=1,MAXYRS 
C 
C  Set IDX 
C 
        IDX = (MAXYRS + 1) - JDX 
C 
C  If the travel fraction for this age/vehicle combination is zero, skip. 
C 
        IF(TF(IDX,IV) .EQ. 0.) CYCLE 
C 
C  Set model year 
C 
        MY=ICY-JDX+1 
C 
C  Determine the BER associated sulfur level to in-use sulfur level 
C  correction factor. 
C 
        IF(VVGASMC(IV).EQ.1) THEN 
          SULFCF = S_PPM / BERSULF(JDX,IV) 
        ELSE 
          SULFCF = D_PPM / BERSULF(JDX,IV) 
        END IF 
C 
C  Compute the NGV fraction 
C 
        IF(NGVFLG1.EQ.1) THEN 
           IF(MY.LT.NGV_START) THEN 
              NGVFRAC = 0.0 
           ELSE 
              NGVPT = MY - NGV_START + 1 
              NGVFRAC = NGVPHIN(NGVPT,IV) 
           END IF 
        ELSE 
           NGVFRAC = 0.0 
        END IF 
C 
C  Determine the ZEV fraction for LDGV, LDGT1-LDGT4. 
C 
        ZEVFRAC = PCTLEV(MY,IV,11) 
C 
C  Combine the NGV and ZEV adjustments. This form assumes that 
C  the ZEVs are removed first, and the NGVs are from the remainder. 
C 
        CADJ = SULFCF * (1.0 - ZEVFRAC) * (1.0 - NGVFRAC)  !Gas portion 
     *       + SLFNGV * NGVFRAC                            !NG  portion 
C 
C  Initialize 
C 
        SUMSO4JDX = 0. 
        SUMSO2JDX = 0. 
C 
C  Loop for hour 
C 
        DO IH=1,MAXIH 
C 
         SUMSO4IH = 0. 
         SUMSO2IH = 0. 
C 
C  Loop for road type 
C 
         DO IROAD=1,4 
C 
C  Adjust for sulfur, ngv, zev 
C 
          DBPM(1,IROAD,IH,JDX,IV) = DBPM(1,IROAD,IH,JDX,IV) * CADJ  !SO4 
          DBPM(6,IROAD,IH,JDX,IV) = DBPM(6,IROAD,IH,JDX,IV) * CADJ  !SO2 
C 
          SUMSO4IH = SUMSO4IH+DBPM(1,IROAD,IH,JDX,IV)*FVMT(IROAD,IH,IV)  !SO4 
          SUMSO2IH = SUMSO2IH+DBPM(6,IROAD,IH,JDX,IV)*FVMT(IROAD,IH,IV)  !SO2 
C   
         END DO  !IROAD 
C 
         SUMSO4JDX = SUMSO4JDX + SUMSO4IH * HVMT(IH)  !SO4 
         SUMSO2JDX = SUMSO2JDX + SUMSO2IH * HVMT(IH)  !SO2 
C 
        END DO  !IH 
C 
        EFPM(1,IV) = EFPM(1,IV) + SUMSO4JDX * TF(IDX,IV)   !SO4 
        EFPM(6,IV) = EFPM(6,IV) + SUMSO2JDX * TF(IDX,IV)   !SO2 
C 
       END DO  !JDX 
C 
       VPM(1) = VPM(1) + EFPM(1,IV) * VMTMIX(IV)   !SO4 
       VPM(6) = VPM(6) + EFPM(6,IV) * VMTMIX(IV)   !SO2 
C 
      END DO  !IV 
C 
      RETURN 
      END 
