      SUBROUTINE METHANE(IV,IP,MY,ISTD,IFDS,IDX,ISR,LEF) 
C 
C   This subroutine calculates the correct HC speciation based on 
C   user inputs. 
C 
C   Called By:  HRLOOP 
C 
C   CALLs:      VOCFID, NMOGFID 
C 
C   Changes: 
C 
C  04 Oct 00 @EPA-djb Bug276 replaced HDGMET with HDDMET in line 143, lddv&t 
C                     calculations 
C  04 Oct 00 @EPA-djb Bug276 Speciated HDGMET, HDDMET & MOFFMTH application. 
C  27 Jun 00 @EPA-ddj Bug174 Completed implementation of bug116 fix 
C            for HDGV. 
C  15 Jun 00 @EPA-BAG Removed INH from parameter list; INH is not 
C                     used in the METHANE subroutine 
C  15 Jun 00 @EPA-ddj Bug149 Added SAVE to PARAMETER type declaration. 
C   30 May 00 @EPA-ddj Bug116, HC species for HDGVs now in NMOG for 
C             2005 model years. 
C   28 Apr 00 @EPA-elg    Methane for light-duty gas now ratioed to 
C             composite High / Normal emitter.  NGV methane code moved 
C             to subroutine NGV_EF.  New methane ratios for 1988+ HDGV, 
C             HDDT, LDDV/T and MC vehicles. 
C   22 Feb 00 @EPA-elg    New Subroutine Created 
C 
C  Input on call: 
C    argument list: IV,IP,MY,ISTD,IFDS,IDX,ISR,INH,LEF 
C 
C    common blocks: 
C    /FLAGS1/   OXYFLG 
C    /FLAGS4/   NMHFLG 
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
      INCLUDE 'FLAGS4.I' 
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
      REAL,  DIMENSION(2) :: LEF 
      INTEGER, INTENT(IN) :: IV 
      INTEGER, INTENT(IN) :: IP 
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
C     Methane applies the proper HC Speciation Correction based on user input. 
C 
      IF (IP.EQ.1) THEN 
C 
C     Determine the methane values for 1981+ LDV and LDT 
C     MOFFMTH is calculated for the THC case and must be adjusted if 
C     the basic emission rate is not in THC units. 
C 
        IF(VVLDGAS(IV).EQ.1 .AND. MY.GE.1981) THEN 
          IF(ISTD.EQ.2.OR.ISTD.EQ.3) THEN           !Base NMHC units 
            OFFMTH(IDX,IV,ISR,ISTD,IFDS) = LEF(ISR) / 
     *        ((1.0 /MOFFMTH(IDX,IV,ISR,ISTD,IFDS))-1.0) 
          ELSEIF(ISTD.GE.4) THEN                    !Base NMOG units 
            IF(RFGON .AND. OXYFLG.EQ.1) THEN 
              OFFMTH(IDX,IV,ISR,ISTD,IFDS) = LEF(ISR) / 
     *         (NMOGFID(MY,IV,2) * 
     *         ((1.0 /MOFFMTH(IDX,IV,ISR,ISTD,IFDS))-1.0)) 
            ELSE 
              OFFMTH(IDX,IV,ISR,ISTD,IFDS) = LEF(ISR) / 
     *         (NMOGFID(MY,IV,OXYFLG) * 
     *         ((1.0 /MOFFMTH(IDX,IV,ISR,ISTD,IFDS))-1.0)) 
            ENDIF 
          ELSE                                      !Base THC units 
            OFFMTH(IDX,IV,ISR,ISTD,IFDS) = LEF(ISR) * 
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
C 
           IF(VVLDDSL(IV).EQ.1 .AND. MY.GE.1988) THEN 
              OFFMTH(IDX,IV,ISR,1,IFDS)  = LDDMET * LEF(ISR) 
           ELSEIF(VVHDGV(IV).EQ.1 .AND. MY.GE.1988.AND.MY.LE.2004) THEN 
              OFFMTH(IDX,IV,ISR,1,IFDS)  = HDGMET * LEF(ISR) 
           ELSEIF(VVHDGV(IV).EQ.1 .AND. MY.GE.2005) THEN   !NMOG units 
             IF(RFGON .AND. OXYFLG.EQ.1) THEN 
               OFFMTH(IDX,IV,ISR,1,IFDS)  = 
     *           LEF(ISR)/(NMOGFID(MY,IV,2)*((1.0 /HDGMET)-1.0)) 
             ELSE 
               OFFMTH(IDX,IV,ISR,1,IFDS)  = 
     *           LEF(ISR)/(NMOGFID(MY,IV,OXYFLG)*((1.0 /HDGMET)-1.0)) 
             ENDIF 
           ELSEIF(VVHDDV(IV).EQ.1 .AND. MY.GE.1988.AND.MY.LE.2003) THEN 
              OFFMTH(IDX,IV,ISR,1,IFDS)  = HDDMET * LEF(ISR) 
           ELSEIF(VVHDDV(IV).EQ.1 .AND. MY.GE.2004) THEN   !NMHC units 
              OFFMTH(IDX,IV,ISR,1,IFDS)  = LEF(ISR)/((1.0 /HDDMET)-1.0) 
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
C          IF (NMHFLG .EQ. 1) THEN                   !THC UNITS 
C              No effect 
C 
           IF (NMHFLG .EQ. 2) THEN                   !NMHC UNITS 
             LEF(ISR) = LEF(ISR) - OFFMTH(IDX,IV,ISR,1,IFDS) 
C 
           ELSEIF (NMHFLG .EQ. 3) THEN               !VOC UNITS 
             LEF(ISR) = (LEF(ISR) - OFFMTH(IDX,IV,ISR,1,IFDS) ) * 
     *                     VOCFID(MY,IV) 
C 
C     Effect of RFG/OXY on TOG is controlled in NMOGFID. 
C 
           ELSEIF (NMHFLG .EQ. 4) THEN               !TOG UNITS 
             IF(RFGON .AND. OXYFLG.EQ.1) THEN 
               LEF(ISR) = (LEF(ISR) - OFFMTH(IDX,IV,ISR,1,IFDS) ) * 
     *                    NMOGFID(MY,IV,2) + OFFMTH(IDX,IV,ISR,1,IFDS) 
             ELSE 
               LEF(ISR) = (LEF(ISR) - OFFMTH(IDX,IV,ISR,1,IFDS) ) * 
     *               NMOGFID(MY,IV,OXYFLG) + OFFMTH(IDX,IV,ISR,1,IFDS) 
             ENDIF 
C 
           ELSEIF (NMHFLG .EQ. 5) THEN               !NMOG UNITS 
             IF(RFGON .AND. OXYFLG.EQ.1) THEN 
               LEF(ISR) = (LEF(ISR) - OFFMTH(IDX,IV,ISR,1,IFDS) ) * 
     *                       NMOGFID(MY,IV,2) 
             ELSE 
               LEF(ISR) = (LEF(ISR) - OFFMTH(IDX,IV,ISR,1,IFDS) ) * 
     *                       NMOGFID(MY,IV,OXYFLG) 
             ENDIF 
           ENDIF                      !NMHFLG ENDIF 
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
           IF (NMHFLG .EQ. 1) THEN                   !TOTAL HC UNITS 
               LEF(ISR) = LEF(ISR) + OFFMTH(IDX,IV,ISR,ISTD,IFDS) 
C 
C          ELSEIF (NMHFLG .EQ. 2) THEN               !NMHC UNITS 
C              No effect 
C 
           ELSEIF (NMHFLG .EQ. 3) THEN               !VOC UNITS 
               LEF(ISR) = LEF(ISR) * VOCFID(MY,IV) 
C 
C     Effect of RFG/OXY on TOG is controlled in NMOGFID. 
C 
           ELSEIF (NMHFLG .EQ. 5) THEN               !NMOG UNITS 
              IF( VVHDDV(IV).EQ.1 .OR. (RFGON .AND. OXYFLG.EQ.1) ) THEN 
                 LEF(ISR) = (LEF(ISR) * NMOGFID(MY,IV,2)) 
              ELSE 
                 LEF(ISR) = (LEF(ISR) * NMOGFID(MY,IV,OXYFLG)) 
              ENDIF 
           ELSEIF (NMHFLG .EQ. 4) THEN               !TOG UNITS 
              IF(  VVHDDV(IV).EQ.1 .OR. (RFGON .AND. OXYFLG.EQ.1) ) THEN 
                 LEF(ISR) = (LEF(ISR) * NMOGFID(MY,IV,2)) + 
     *                       OFFMTH(IDX,IV,ISR,ISTD,IFDS) 
              ELSE 
                 LEF(ISR) = (LEF(ISR) * NMOGFID(MY,IV,OXYFLG)) + 
     *                       OFFMTH(IDX,IV,ISR,ISTD,IFDS) 
              ENDIF 
           ENDIF                      !NMHFLG ENDIF 
C 
        ELSE 
C        IF( (MY.GE.2005 .AND. VVHDGV(IV).EQ.1 )   .OR. 
C            (ISTD.GE.4 .AND. VVLDGAS(IV).EQ.1 ) ) THEN   !NMOG Base 
C 
C     Effect of RFG/OXY on TOG is controlled in NMOGFID. 
C 
          IF (NMHFLG .EQ. 1) THEN                   !TOTAL HC UNITS 
            IF(RFGON .AND. OXYFLG.EQ.1) THEN 
               LEF(ISR) = (LEF(ISR) / NMOGFID(MY,IV,2)) + 
     *                     OFFMTH(IDX,IV,ISR,ISTD,IFDS) 
            ELSE 
               LEF(ISR) = (LEF(ISR) / NMOGFID(MY,IV,OXYFLG)) + 
     *                     OFFMTH(IDX,IV,ISR,ISTD,IFDS) 
            ENDIF 
C 
          ELSEIF (NMHFLG .EQ. 2) THEN               !NMHC UNITS 
            IF(RFGON .AND. OXYFLG.EQ.1) THEN 
               LEF(ISR) = (LEF(ISR) / NMOGFID(MY,IV,2)) 
            ELSE 
               LEF(ISR) = (LEF(ISR) / NMOGFID(MY,IV,OXYFLG)) 
            ENDIF 
C 
          ELSEIF (NMHFLG .EQ. 3) THEN               !VOC UNITS 
            IF(RFGON .AND. OXYFLG.EQ.1) THEN 
               LEF(ISR) = (LEF(ISR) / NMOGFID(MY,IV,2)) * 
     *                     VOCFID(MY,IV) 
            ELSE 
               LEF(ISR) = (LEF(ISR) / NMOGFID(MY,IV,OXYFLG)) * 
     *                     VOCFID(MY,IV) 
            ENDIF 
C 
          ELSEIF (NMHFLG .EQ. 4) THEN               !TOG UNITS 
            LEF(ISR) = LEF(ISR) + OFFMTH(IDX,IV,ISR,ISTD,IFDS) 
C 
C         ELSEIF (NMHFLG .EQ. 5) THEN               !NMOG UNITS 
C           No effect 
          ENDIF                        !NMHFLG ENDIF 
        ENDIF                          !HC Base NMOG ENDIF 
      ENDIF                            !IP = HC ENDIF 
C 
      RETURN 
      END 
