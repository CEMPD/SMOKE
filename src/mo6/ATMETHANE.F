      SUBROUTINE ATMETHANE(IV,MY,ISTD,IFDS,IDX,ISR,LEFVOC,LEFTOG) 
C 
C   This subroutine calculates the VOC and TOG HC speciation for 
C   Air Toxics calculations. (Based on METHANE subroutine.) 
C 
C   Called By:  ATHRLOOP 
C 
C   CALLs:      VOCFID, NMOGFID 
C 
C   Changes: 
C 
C   09 Jul 01 AIR Task 15:  New Subroutine Created 
C 
C  Input on call: 
C    argument list: IV,MY,ISTD,IFDS,IDX,ISR,INH,LEF 
C 
C    common blocks: 
C    /FLAGS1/   OXYFLG 
C    /VVNAMS/   VVLDDSL, VVLDGAS, VVHDGV, VVHDDV  (IVTYPE.I) 
C    /VVTYPS/   VTMC  (IVTYPE.I) 
C    /OFFSET/   OFFMTH(IDX,IV,ISR,ISTD,IFDS), MOFFMTH(IDX,IV,ISR,ISTD,IFDS) 
C    /RFORM1/   RFGON 
C 
C  Output on return: 
C    common blocks: 
C 
C  Local variable dictionary: 
C 
C   Name        Type              Description 
C  ------       ----  --------------------------------------------------- 
C   HDGMET       R     Methane percentage offset for all HDG classes. 
C   HDDMET       R     Methane percentage offset for all HDD classes. 
C   LDDMET       R     Methane percentage offset for LDDV/T. 
C   LEF          R     Basic emission factor sent to HRLOOP. 
C   OFFMTH       R     Additive Methane emission adjustment factor. 
C   MCMET        R     Methane start/running ratio for motorcycles. 
C   MOFFMTH      R     Multiplicative Methane adjustment factor for IV=1 thru 5. 
C 
C   Notes: 
C 
C 
      IMPLICIT NONE 
C 
      INCLUDE 'FLAGS1.I' 
      INCLUDE 'IVTYPE.I' 
      INCLUDE 'OFFSET.I' 
      INCLUDE 'RFORM1.I' 
C 
C     Declare called functions 
C 
      REAL NMOGFID 
      REAL VOCFID 
C 
C     Declare parameter list 
C 
      REAL,  DIMENSION(2) :: LEFVOC 
      REAL,  DIMENSION(2) :: LEFTOG 
      INTEGER, INTENT(IN) :: IV 
      INTEGER, INTENT(IN) :: MY 
      INTEGER, INTENT(IN) :: ISTD 
      INTEGER, INTENT(IN) :: IFDS 
      INTEGER, INTENT(IN) :: IDX 
      INTEGER, INTENT(IN) :: ISR 
C 
C     Declare local variables 
C 
      REAL, PARAMETER    :: LDDMET = 0.01 
      REAL, PARAMETER    :: HDGMET = 0.116 
      REAL, PARAMETER    :: HDDMET = 0.046 
      REAL, DIMENSION(2), SAVE :: MCMET  = (/ 0.792, 2.1225 /) 
C 
C     Determine the methane values for 1981+ LDV and LDT 
C     MOFFMTH is calculated for the THC case and must be adjusted if 
C     the basic emission rate is not in THC units. 
C     (Programming Note: At this point ATVOC=ATTOG, so OK to use 
C      either array in the following OFFMTH calculations.) 
C 
        IF(VVLDGAS(IV).EQ.1 .AND. MY.GE.1981) THEN 
          IF(ISTD.EQ.2.OR.ISTD.EQ.3) THEN           !Base NMHC units 
            OFFMTH(IDX,IV,ISR,ISTD,IFDS) = LEFVOC(ISR) / 
     *        ((1.0 /MOFFMTH(IDX,IV,ISR,ISTD,IFDS))-1.0) 
          ELSEIF(ISTD.GE.4) THEN                    !Base NMOG units 
            IF(RFGON .AND. OXYFLG.EQ.1) THEN 
              OFFMTH(IDX,IV,ISR,ISTD,IFDS) = LEFVOC(ISR) / 
     *         (NMOGFID(MY,IV,2) * 
     *         ((1.0 /MOFFMTH(IDX,IV,ISR,ISTD,IFDS))-1.0)) 
            ELSE 
              OFFMTH(IDX,IV,ISR,ISTD,IFDS) = LEFVOC(ISR) / 
     *         (NMOGFID(MY,IV,OXYFLG) * 
     *         ((1.0 /MOFFMTH(IDX,IV,ISR,ISTD,IFDS))-1.0)) 
            ENDIF 
          ELSE                                      !Base THC units 
            OFFMTH(IDX,IV,ISR,ISTD,IFDS) = LEFVOC(ISR) * 
     *        MOFFMTH(IDX,IV,ISR,ISTD,IFDS) 
          ENDIF  !ISTD Conditional 
        ENDIF    !VVLDGAS/MY Conditional 
C 
C     Methane calculation for IV = 14, 15, and 28 (LDDV/T) uses a 
C     constant 1.0% methane offset (LDDMET).  The values in array ZDMTH3 
C     for the methane offset are not used for 1987 and newer LDDV/T. 
C 
C     1997 and newer HDGV and HDDV methane is calculated here and 
C     do not use ZDMTH values.  These calculations must take into 
C     account that the HDGMET and HDDMET value are based on basic 
C     emission rates in THC units and 2005+ HDG are in NMOG units and 
C     2004+ HDD are in NMHC units. 
C 
C     MCMET for start (2.1225 = 0.283 * 7.5) both estimates the start 
C     fraction and converts the results to a methane start offset 
C     in grams/start. 
C     (Programming Note: At this point ATVOC=ATTOG, so OK to use 
C      either array in the following OFFMTH calculations.) 
C 
           IF(VVLDDSL(IV).EQ.1 .AND. MY.GE.1988) THEN 
              OFFMTH(IDX,IV,ISR,1,IFDS)  = LDDMET * LEFVOC(ISR) 
           ELSEIF(VVHDGV(IV).EQ.1 .AND. MY.GE.1988.AND.MY.LE.2004) THEN 
              OFFMTH(IDX,IV,ISR,1,IFDS)  = HDGMET * LEFVOC(ISR) 
           ELSEIF(VVHDGV(IV).EQ.1 .AND. MY.GE.2005) THEN   !NMOG units 
             IF(RFGON .AND. OXYFLG.EQ.1) THEN 
               OFFMTH(IDX,IV,ISR,1,IFDS)  = 
     *           LEFVOC(ISR)/(NMOGFID(MY,IV,2)*((1.0 /HDGMET)-1.0)) 
             ELSE 
               OFFMTH(IDX,IV,ISR,1,IFDS)  = 
     *           LEFVOC(ISR)/(NMOGFID(MY,IV,OXYFLG)*((1.0 /HDGMET)-1.0)) 
             ENDIF 
           ELSEIF(VVHDDV(IV).EQ.1 .AND. MY.GE.1988.AND.MY.LE.2003) THEN 
              OFFMTH(IDX,IV,ISR,1,IFDS)  = HDDMET * LEFVOC(ISR) 
           ELSEIF(VVHDDV(IV).EQ.1 .AND. MY.GE.2004) THEN   !NMHC units 
               OFFMTH(IDX,IV,ISR,1,IFDS)  =  
     *            LEFVOC(ISR)/((1.0 /HDDMET)-1.0) 
           ELSEIF(IV .EQ. VTMC) THEN 
              OFFMTH(IDX,IV,ISR,1,IFDS)  = OFFMTH(IDX,IV,ISR,1,IFDS) * 
     *                                     MCMET(ISR) 
           ENDIF 
C 
        IF( (ISTD.EQ.1  .AND.  VVLDGAS(IV).EQ.1 )   .OR. 
     *      (MY.LE.2004 .AND.  VVHDGV(IV) .EQ.1 )   .OR. 
     *      (MY.LE.2003 .AND.  VVHDDV(IV) .EQ.1 )   .OR. 
     *      (                  VVLDDSL(IV).EQ.1 )   .OR. 
     &      (                  IV .EQ. VTMC     ) ) THEN   !THC Base 
C 
C     Adjust the HC species for all vehicles with base emission THC units. 
C 
C     Base HC units for LDGV & LDGT are always THC prior to 1994. 
C     Base HC units for LDDV & LDDT are always THC. 
C     Base HC units for HDGV are THC for all model years prior to 2005. 
C     Base HC units for HDDV are THC for all model years prior to 2004. 
C     Base HC units for MC are always THC. 
C 
C     VOC Units 
          LEFVOC(ISR) = (LEFVOC(ISR) - OFFMTH(IDX,IV,ISR,1,IFDS) ) * 
     *                       VOCFID(MY,IV) 
C 
C     TOG Units, Effect of RFG/OXY on TOG is controlled in NMOGFID. 
          IF(RFGON .AND. OXYFLG.EQ.1) THEN 
            LEFTOG(ISR) = (LEFTOG(ISR) - OFFMTH(IDX,IV,ISR,1,IFDS) ) * 
     *                NMOGFID(MY,IV,2) + OFFMTH(IDX,IV,ISR,1,IFDS) 
          ELSE 
            LEFTOG(ISR) = (LEFTOG(ISR) - OFFMTH(IDX,IV,ISR,1,IFDS) ) * 
     *           NMOGFID(MY,IV,OXYFLG) + OFFMTH(IDX,IV,ISR,1,IFDS) 
          ENDIF 
C 
C     Adjust the HC species for all vehicles with base emission in non-THC units. 
C 
C     Base HC units for LDGV & LDGT are always THC for Tier 0. 
C     Base HC units for Tier 1 ISTD's are Non-Methane Hydrocarbons. 
C     Base HC units for LEV, ULEV and Tier2 ISTD's are NMOG. 
C     Base HC units for HDGV are THC for all model years prior to 2005. 
C     Bace HC units for HDGV are NMOG for 2005+ model years. 
C     Base HC units for LDDV/T are always THC. 
C     Base HC units for MC are always THC. 
C     Base HC units for HDDV are THC for all model years prior to 2004. 
C     Base HC units for HDDV and Urban buses for 2004+ MY is NMHC. 
C 
        ELSEIF( (ISTD .LT. 4) .AND. 
     &        ( (MY.GE.1994 .AND. VVLDGAS(IV).EQ.1 )    .OR. 
     &          (MY.GE.2004 .AND. VVHDDV(IV) .EQ.1 ) )) THEN  !NMHC Base 
C 
C     VOC Units 
          LEFVOC(ISR) = LEFVOC(ISR) * VOCFID(MY,IV) 
C 
C     TOG Units, Effect of RFG/OXY on TOG is controlled in NMOGFID. 
          IF(  VVHDDV(IV).EQ.1 .OR. (RFGON .AND. OXYFLG.EQ.1) ) THEN 
            LEFTOG(ISR) = (LEFTOG(ISR) * NMOGFID(MY,IV,2)) + 
     *                    OFFMTH(IDX,IV,ISR,ISTD,IFDS) 
          ELSE 
            LEFTOG(ISR) = (LEFTOG(ISR) * NMOGFID(MY,IV,OXYFLG)) + 
     *                    OFFMTH(IDX,IV,ISR,ISTD,IFDS) 
          ENDIF 
C 
        ELSE 
C 
C     NMOG Base 
C 
C     VOC Units 
        IF(RFGON .AND. OXYFLG.EQ.1) THEN 
          LEFVOC(ISR) = (LEFVOC(ISR) / NMOGFID(MY,IV,2)) * 
     *                  VOCFID(MY,IV) 
        ELSE 
          LEFVOC(ISR) = (LEFVOC(ISR) / NMOGFID(MY,IV,OXYFLG)) * 
     *                     VOCFID(MY,IV) 
        ENDIF 
C 
C     TOG Units, Effect of RFG/OXY on TOG is controlled in NMOGFID. 
            LEFTOG(ISR) = LEFTOG(ISR) + OFFMTH(IDX,IV,ISR,ISTD,IFDS) 
C 
        ENDIF                          !HC Base NMOG ENDIF 
C 
      RETURN 
      END 
