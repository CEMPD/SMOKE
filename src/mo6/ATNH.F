      SUBROUTINE ATNH(IV,MY,MKTSHR,CAT) 
C 
C  This subroutine weights the 9 normal emission conversions 
C  into one Normal conversion, and places the market-share 
C  weighted normal and high values into the ATCONV array. 
C  Complex Model based conversion rates are based on VOC, 
C  so IVT is always 1. 
C 
C  Changes (Last listed first) 
C  20 Sep 01: AIR Task 15-22: Comments 
C  09 Sep 01: AIR Task 15-22: Disabled Case 2 per EPA direction 
C  10 Jul 01: AIR Task 15-22: New MOBILE6 Subroutine 
C 
C  Called by: ATCALX 
C 
C  Calls: None 
C   
C 
      IMPLICIT NONE 
C 
      INCLUDE 'ATOX3.I'  ! GROUP,VOCCOMPLEX,BENZRATION,MTBERATION,BUTARATIO,FORMRATIO,ACETRATIO 
      INCLUDE 'ATOX4.I'  ! ATNHWT,ATCONV 
      INCLUDE 'IVTYPE.I' ! VTLDGV 
C 
C  Declare parameter list 
C  
      INTEGER IV      ! Vehicle class 
      INTEGER MY      ! Model year 
      REAL    CAT(3)  ! Catalyst fractions (1=No Cat, 2=Ox, 3=3-way) 
      REAL    MKTSHR  ! Market share 
C 
C  Declare parameter list 
C  
      INTEGER INH     ! Normal/High toggle 
      INTEGER ITG     ! Tech group index 
      INTEGER ITGUSE  ! Temporary tech group index 
      REAL    VOC     ! VOC from Complex Model, in g/mi) 
      REAL    WT      ! Weighting factor 
C 
C  Loop for the 10 Tech Groups 
C 
      DO ITG=1,10 
C 
C  Set normal and high emitter index 
C 
        INH=1                !Normals 1-9 
        IF(ITG.EQ.10) INH=2  !High 10 
C 
C  Set weighting factor. Weight factor must include fraction of 
C  3-ways so that VOC-based and TOG-based values will properly 
C  sum to 100% 
C 
        WT=GROUP(ITG) * MKTSHR * CAT(3) 
C 
C  For 1982+ LDGV only: 
C  Case 1) If the Normal Emitter VOC computed by the Complex Model  
C          exceeds 0.6 g/mi, then use the High Emitter Ratio. Only  
C          do this if the MOBILE6 High Emitter Fraction is greater  
C          than 0. Do NOT reset the MOBILE6-computed Emitter Fractions. 
C  Case 2: If the High Emitter VOC computed by the Complex Model is 
C          equal to or less than 0.6 g/mi, use the only the Normal Emitter 
C          Ratios. Only do this if the MOBILE6 Normal Emitter Fraction 
C          is greater than 0. Reset the MOBILE6-Emitter Fractions. 
C  NOTE: CASE 2 IS DISABLED BECAUSE THE COMPLEX MODEL IS LOOPED OVER 
C        ALL 5 POSSIBLE OXYGENATE COMBINATIONS. 
C 
        ITGUSE=ITG 
C 
        IF(IV.EQ.VTLDGV .AND. MY.GE.1982) THEN 
          VOC = VOCCOMPLEX(ITG)/1000. ! Convert from mg/mi to g/mi 
          IF(ITG.LE.9 .AND. VOC.GT.0.6 .AND. ATNHWT(2).GT.0.) ITGUSE=10 
C        IF(ITG.EQ.10.AND. VOC.LE.0.6 .AND. ATNHWT(1).GT.0.) THEN 
C        WT=0. 
C        ATNHWT(1)=1. 
C        ATNHWT(2)=0. 
C        END IF 
        END IF 
C 
C  Compute the Exhaust AT Conversion Factors, weighted by Tech Group 
C  and Market Share. 
C 
        ATCONV(INH,1,1) = ATCONV(INH,1,1) + BENZRATIO(ITGUSE) * WT  !Benzene 
        ATCONV(INH,1,2) = ATCONV(INH,1,2) + MTBERATIO(ITGUSE) * WT  !MTBE 
        ATCONV(INH,1,3) = ATCONV(INH,1,3) + BUTARATIO(ITGUSE) * WT  !Butadiene 
        ATCONV(INH,1,4) = ATCONV(INH,1,4) + FORMRATIO(ITGUSE) * WT  !Formaldehyde 
        ATCONV(INH,1,5) = ATCONV(INH,1,5) + ACETRATIO(ITGUSE) * WT  !Acetaldehyde 
C 
      END DO 
C 
      RETURN 
      END 
