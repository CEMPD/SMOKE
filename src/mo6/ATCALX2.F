      SUBROUTINE ATCALX2(CAT,IV,IOXY) 
C 
C  This subroutine compute the fraction of TOG which is emitted as 
C  benzene, formaldehyde, acetaldehyde, 1,3 butadiene and MTBE for 
C  all vehicle/catalyst combinations not covered by the Complex Model.  
C  Changes:  (Last change first) 
C 
C  13 Aug 01 AIR Task 20: Updated equations to table sent 09Aug01. 
C  26 Jul 01 AIR Task 20: New Subroutine. 
C 
C  Input on call: 
C 
      IMPLICIT NONE 
C 
      INCLUDE 'ATOX1.I' 
      INCLUDE 'ATOX4.I' 
      INCLUDE 'IVTYPE.I' 
C 
      INTEGER ICAT 
      INTEGER INH 
      INTEGER IOXY 
      INTEGER IV 
      REAL    ACET1 
      REAL    BENZ1 
      REAL    BUTA1 
      REAL    CAT(3) 
      REAL    FORM1 
      REAL    MS 
      REAL    MTBE1 
      REAL    WT 
C 
C Programming Note: Due to lack of data, TAME is assumed to behave like 
C MTBE, and ETBE is assumed to behave like ETOH. 
C 
C Programming Note: All oxygenate levels are in wt% (via ATOXY(1,IOXY)) 
C 
C  Determine the market share. Reset Diesels to 100% 
C 
      MS=0.0 
C 
      IF(VVDSL(IV).EQ.1) THEN 
        IF(IOXY.EQ.5) MS=1.0 
      ELSE 
        MS=ATOXY(2,IOXY) 
      END IF 
C 
C  Skip if fuel market share is zero. 
C 
       IF(MS .EQ. 0.) RETURN 
C 
C  Loop for catalyst 
C 
        DO ICAT=1,3  !1:NoCat/Diesel, 2:Ox, 3:3-way 
C 
C  Skip if catalyst fraction is zero. 
C 
        IF(CAT(ICAT) .EQ. 0.) CYCLE 
C 
C******************** BENZENE ******************** 
C 
C Programming Note: The original Benzene equations for gasoline vehicles 
C were found to produce percents, not fractions, like all other equations. 
C Therefore, the percent to fraction conversions were needed. 
C 
      BENZ1 = 0. 
C 
C  LDGV,LDGT1-4/oxcat, MC/nocat, LDGV,LDGT1-4/nocat, HDGV/nocat 
C  All baseline gasoline and oxygenated fuels. 
C 
      IF(VVLDGAS(IV) .EQ. 1    .AND. ICAT .LE. 2 .OR.   !LDGV,LDGT1-4/NoCat,Ox 
     *           IV  .EQ. VTMC .AND. ICAT .EQ. 1 .OR.   !MC/NoCat 
     *    VVHDGV(IV) .EQ. 1    .AND. ICAT .EQ. 1) THEN  !HDGV/NoCat 
        
        BENZ1 = 0.8551  * BENZ 
     *        + 0.12198 * AROM 
     *        - 1.1626 
C 
        BENZ1 = BENZ1 / 100.   ! Convert from percent to fraction.  
C 
      END IF   !LDGV,HDGV,MC 
C 
C  HDGV/cat 
C  All baseline gasoline and oxygenated fuels. 
C 
      IF(VVHDGV(IV) .EQ. 1 .AND. ICAT .GE. 2) THEN  !HDGV/cat 
C 
        BENZ1 = 0.7732 * BENZ 
     *        + 0.0987 * (AROM - BENZ) 
     *        + 1.077 
C 
        BENZ1 = BENZ1 / 100.   ! Convert from percent to fraction.  
C 
      END IF   !HDGV/cat 
C 
C  LDDV/LDDT 
C 
      IF(VVLDDSL(IV).EQ.1) BENZ1 = 0.0200 
C 
C  HDDV 
C 
      IF(VVHDDV(IV).EQ.1) BENZ1 = 0.0105 
C 
C******************** FORMALDEHYDE ******************** 
C 
      FORM1 = 0. 
C 
C  LDGV,LDGT1-4,MC 
C 
      IF(VVLDGAS(IV) .EQ. 1     .OR. 
     *           IV  .EQ. VTMC) THEN   
C 
        IF(ICAT .EQ. 2) THEN     !OxCat      
C 
        SELECT CASE (IOXY) 
          CASE (1,4)                                 !MTBE,TAME                   
            FORM1 = 0.0151 
     *            + 0.0151 * 1.2082 * ATOXY(1,IOXY) / 2.7 
          CASE (2,3)                                 !ETOH,ETBE 
            FORM1 = 0.0151 
     *            + 0.0151 * 0.3350 * ATOXY(1,IOXY) / 3.5 
          CASE (5)                                   !NOOXY 
            FORM1 = 0.0151 
          END SELECT 
C 
        ELSE IF (ICAT .EQ. 1) THEN  !NoCat 
C 
        SELECT CASE (IOXY)  
          CASE (1,4)                                 !MTBE,TAME 
            FORM1 = 0.0224 
     *            + 0.0224 * 0.4336 * ATOXY(1,IOXY) / 2.7 
          CASE (2,3)                                 !ETOH,ETBE 
            FORM1 = 0.0224 
     *            + 0.0224 * 0.1034 * ATOXY(1,IOXY) / 3.5 
           CASE (5)                                  !NOOXY 
            FORM1 = 0.0224 
          END SELECT 
C 
        END IF  !Nocat,OxCat 
C 
      END IF  !LDGV,LDGT1-4,MC 
C 
C  HDGV 
C 
      IF(VVHDGV(IV) .EQ. 1) THEN 
C   
        IF(ICAT .EQ. 1) THEN     !NoCat  
C     
        SELECT CASE (IOXY) 
          CASE (1,4)                                 !MTBE,TAME 
            FORM1 = 0.0347 
     *            + 0.0347 * 0.1259 * ATOXY(1,IOXY) / 2.7 
          CASE (2,3)                                 !ETOH,ETBE 
            FORM1 = 0.0347 
     *            + 0.0347 * 0.1034 * ATOXY(1,IOXY) / 3.5 
          CASE (5)                                   !NOOXY 
            FORM1 = 0.0347 
          END SELECT 
C 
        ELSE IF (ICAT .GE. 2) THEN  !Cat 
C 
        SELECT CASE (IOXY) 
          CASE (1,4)                                 !MTBE,TAME 
            FORM1 = 0.0151 
     *            + 0.0151 * 1.2082 * ATOXY(1,IOXY) / 2.7 
          CASE (2,3)                                 !ETOH,ETBE 
            FORM1 = 0.0151 
     *            + 0.0151 * 0.3350 * ATOXY(1,IOXY) / 3.5 
          CASE (5)                                   !NOOXY 
            FORM1 = 0.0151 
          END SELECT 
C 
        END IF  !Nocat,Cat 
C 
      END IF  !HDGV 
C 
C  LDDV/LDDT 
C 
      IF(VVLDDSL(IV).EQ.1) FORM1 = 0.0386 
C 
C  HDDV 
C 
      IF(VVHDDV(IV).EQ.1) FORM1 = 0.0782 
C 
C******************** ACETALDEHYDE ******************** 
C 
      ACET1 = 0. 
C 
C  LDGV,LDGT1-4,MC 
C 
      IF(VVLDGAS(IV) .EQ. 1     .OR. 
     *           IV  .EQ. VTMC) THEN   
C 
        IF(ICAT .EQ. 2) THEN     !OxCat      
C 
        SELECT CASE (IOXY) 
          CASE (1,4)                                 !MTBE,TAME                   
            ACET1 = 0.0047 
     *            + 0.0047 * 0.2556 * ATOXY(1,IOXY) / 2.7 
          CASE (2,3)                                 !ETOH,ETBE 
            ACET1 = 0.0047 
     *            + 0.0047 * 2.1074 * ATOXY(1,IOXY) / 3.5 
          CASE (5)                                   !NOOXY 
            ACET1 = 0.0047 
          END SELECT 
C 
        ELSE IF (ICAT .EQ. 1) THEN  !NoCat 
C 
        SELECT CASE (IOXY)  
          CASE (1,4)                                 !MTBE,TAME 
            ACET1 = 0.0060 
     *            + 0.0060 * 0.2303 * ATOXY(1,IOXY) / 2.7 
          CASE (2,3)                                 !ETOH,ETBE 
            ACET1 = 0.0060 
     *            + 0.0060 * 1.1445 * ATOXY(1,IOXY) / 3.5 
          CASE (5)                                   !NOOXY 
            ACET1 = 0.0060 
          END SELECT 
C 
        END IF  !Nocat,OxCat 
C 
      END IF  !LDGV,LDGT1-4,MC 
C 
C  HDGV 
C 
      IF(VVHDGV(IV) .EQ. 1) THEN 
C   
        IF(ICAT .EQ. 1) THEN     !NoCat  
C     
        SELECT CASE (IOXY) 
          CASE (1,4,5)                               !MTBE,TAME,NOOXY 
            ACET1 = 0.0067 
          CASE (2,3)                                 !ETOH,ETBE 
            ACET1 = 0.0067 
     *            + 0.0067 * 1.1445 * ATOXY(1,IOXY) / 3.5 
          END SELECT 
C 
        ELSE IF (ICAT .GE. 2) THEN  !Cat 
C 
        SELECT CASE (IOXY) 
          CASE (1,4)                                 !MTBE,TAME                   
            ACET1 = 0.0047 
     *            + 0.0047 * 0.2556 * ATOXY(1,IOXY) / 2.7 
          CASE (2,3)                                 !ETOH,ETBE 
            ACET1 = 0.0047 
     *            + 0.0047 * 2.1074 * ATOXY(1,IOXY) / 3.5 
          CASE (5)                                   !NOOXY 
            ACET1 = 0.0047 
          END SELECT 
C 
        END IF  !Nocat,Cat 
C 
      END IF  !HDGV 
C 
C  LDDV/LDDT 
C 
      IF(VVLDDSL(IV).EQ.1) ACET1 = 0.0123 
C 
C  HDDV 
C 
      IF(VVHDDV(IV).EQ.1) ACET1 = 0.0288 
C 
C******************** BUTADIENE ******************** 
C 
      BUTA1 = 0. 
C 
C  LDGV,LDGT1-4,MC 
C 
      IF(VVLDGAS(IV) .EQ. 1     .OR. 
     *           IV  .EQ. VTMC) THEN   
C 
        IF(ICAT .EQ. 2) THEN     !OxCat      
C 
        SELECT CASE (IOXY) 
          CASE (1,4)                                 !MTBE,TAME                   
            BUTA1 = 0.0044 
     *            - 0.0044 * 0.2227 * ATOXY(1,IOXY) / 2.7 
          CASE (2,3)                                 !ETOH,ETBE 
            BUTA1 = 0.0044 
     *            - 0.0044 * 0.2804 * ATOXY(1,IOXY) / 3.5 
          CASE (5)                                   !NOOXY 
            BUTA1 = 0.0044 
          END SELECT 
C 
        ELSE IF (ICAT .EQ. 1) THEN  !NoCat 
C 
        SELECT CASE (IOXY)  
          CASE (1,4)                                 !MTBE,TAME 
            BUTA1 = 0.0092 
     *            + 0.0092 * 0.1517 * ATOXY(1,IOXY) / 2.7 
          CASE (2,3)                                 !ETOH,ETBE 
            BUTA1 = 0.0092 
     *            + 0.0092 * 0.1233 * ATOXY(1,IOXY) / 3.5 
           CASE (5)                                  !NOOXY 
            BUTA1 = 0.0092 
          END SELECT 
C 
        END IF  !Nocat,OxCat 
C 
      END IF  !LDGV,LDGT1-4,MC 
C 
C  HDGV 
C 
      IF(VVHDGV(IV) .EQ. 1) THEN 
C   
        IF(ICAT .EQ. 1) THEN     !NoCat  
C     
        SELECT CASE (IOXY) 
          CASE (1,4)                                 !MTBE,TAME 
            BUTA1 = 0.0074 
     *            - 0.0074 * 0.2172 * ATOXY(1,IOXY) / 2.7 
          CASE (2,3)                                 !ETOH,ETBE 
            BUTA1 = 0.0074 
     *            + 0.0074 * 0.1233 * ATOXY(1,IOXY) / 3.5 
          CASE (5)                                   !NOOXY 
            BUTA1 = 0.0074 
          END SELECT 
C 
        ELSE IF (ICAT .GE. 2) THEN   !Cat 
C 
        SELECT CASE (IOXY) 
          CASE (1,4)                                 !MTBE,TAME 
            BUTA1 = 0.0029 
     *            - 0.0029 * 0.3233 * ATOXY(1,IOXY) / 2.7 
          CASE (2,3)                                 !ETOH,ETBE 
            BUTA1 = 0.0029 
     *            - 0.0029 * 0.1188 * ATOXY(1,IOXY) / 3.5 
          CASE (5)                                   !NOOXY 
            BUTA1 = 0.0029 
          END SELECT 
C 
        END IF  !Nocat,Cat 
C 
      END IF  !HDGV 
C 
C  LDDV/LDDT 
C 
      IF(VVLDDSL(IV).EQ.1) BUTA1 = 0.0090 
C 
C  HDDV 
C 
      IF(VVHDDV(IV).EQ.1) BUTA1 = 0.0061 
C 
C************************** MTBE ************************** 
C 
      MTBE1 = 0.0 
C 
C  Only MTBE fuels emit this toxic. 
C 
      IF(IOXY .EQ. 1) THEN 
C 
C  LDGV,LDGT1-4,MC 
C 
      IF(VVLDGAS(IV) .EQ. 1     .OR. 
     *           IV  .EQ. VTMC) THEN   
C 
        IF(ICAT .EQ. 2) MTBE1 = 0.0464 * ATOXY(1,IOXY) / 2.7   !OxCat 
        IF(ICAT .EQ. 1) MTBE1 = 0.0333 * ATOXY(1,IOXY) / 2.7   !NoCat 
C 
      END IF  !LDGV,LDGT1-4,MC 
C 
C  HDGV 
C 
      IF(VVHDGV(IV) .EQ. 1) THEN 
C   
        IF(ICAT .EQ. 1) MTBE1 = 0.0209 * ATOXY(1,IOXY) / 2.7   !NoCat  
        IF(ICAT .GE. 2) MTBE1 = 0.0155 * ATOXY(1,IOXY) / 2.7   !Cat 
C 
        END IF  !HDGV 
C 
      END IF  !IOXY 
C 
C********************* WEIGHTING ************************** 
C 
C  Compute cat and market share weighted percentages.  
C  Low and High Emitters are the same. These values 
C  are TOG based, so IVT is always 2. 
C 
      WT = MS * CAT(ICAT) 
C 
      DO INH = 1,2 
        ATCONV(INH,2,1) = ATCONV(INH,2,1) + BENZ1 * WT 
        ATCONV(INH,2,2) = ATCONV(INH,2,2) + MTBE1 * WT 
        ATCONV(INH,2,3) = ATCONV(INH,2,3) + BUTA1 * WT 
        ATCONV(INH,2,4) = ATCONV(INH,2,4) + FORM1 * WT 
        ATCONV(INH,2,5) = ATCONV(INH,2,5) + ACET1 * WT 
      END DO ! INH 
C 
      END DO ! ICAT 
C 
      RETURN 
      END 
