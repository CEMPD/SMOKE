      REAL FUNCTION SPEED_CF(IH,IV,INH,IP,MY,IDX,ISTD,IROAD, 
     *                       BEFIN,AC_DUM,TCF) 
C 
C  SPEED_CF calculates the speed correction for different road types 
C  (facilities), applies the temperature correction factor and 
C  determines the air conditioning effect. 
C 
C 12 Oct 00 AIR Task 01: Optimizations 
C 25 Sep 00 AIR Task 07: Converted all code to upper case to remove 
C           problems between platforms (Mac, PC, Unix). 
C 30 Aug 00 @EPA-djb Bug 243, Tier 2 uses LEV for calculations of cfs. 
C           Added BASEQ5 include file to get EMFAC. 
C 16 Aug 00 @EPA-djb Bug 223, Added TCF to the argument and apply the 
C           TCF to the speed corrected emissions before air conditioning. 
C 20 Jul 00 @EPA-ddj Bug194  Removed IDX from call to ACCF. 
C 01 May 00 @EPA-elg   Changed the name of this subroutine from 
C           BEFRDSPAD to SPEED_CF. 
C 18 Jan 00 @DynTel-ddj  1-020  Added INH and ISTD to parameter 
C           list.  Removed old Offcycle code and added call to 
C           OFFCYCLE.  Added call to SFTP_SPD. 
C 13 Jan 00 @DynTel-ddj  1-013  Placed AVSPR and AVSPL into common 
C           block (SPEED9.I). 
C 17 Aug 99 @DynTel-ddj  Streamlined calculation for OEAD and PLEVELn. 
C           Also changed value for PLEVEL1 from 1.0 to 0.0 when FEFN 
C           greater than or equal to LEVEL(3,IP).  HC calculation 
C           changed from THC to NMHC. 
C 12 May  1999 @DynTel-HQ 2-698 initialzed AVSPR=34.9 rather than 34.6 
C  1 March 1999 @DynTel-HQ 2-698 calculate SCF table 
C                          using MOBILE5 method. 
C 
C    Called by:  HRLOOP 
C 
C 
C    Calls:  OFFCYCLE, SFTP_SPD 
C 
C 
C Input on call: 
C 
C   argument list : IH, IV, INH, IP, MY, IDX, ISTD, IROAD, BEF, AC_DUM 
C 
C 
C Output on return: 
C 
C   function : SPEED_CF 
C 
C   Variables: AC_DUM 
C 
C 
C    Constants: 
C    MAXVEH from include file IVTYPE.I 
C 
C     common blocks: 
C     /BASEQ5/  EMFAC 
C     /BDSPEED1/ TVMT, HVMT, FVMT, SVMT, LOCAL, RAMP, 
C                OCECO, OCY, OCX, LEVEL, EMSFC 
C 
      IMPLICIT NONE 
C 
      INCLUDE 'BASEQ5.I' 
      INCLUDE 'FLAGS5.I' 
      INCLUDE 'SPEED9.I' 
      INCLUDE 'IVTYPE.I' 
C 
      INTEGER, INTENT(IN)    :: IH 
      INTEGER, INTENT(IN)    :: IV 
      INTEGER, INTENT(IN)    :: INH 
      INTEGER, INTENT(IN)    :: IP 
      INTEGER, INTENT(IN)    :: IDX 
      INTEGER, INTENT(IN)    :: IROAD 
      INTEGER, INTENT(IN)    :: ISTD 
      INTEGER, INTENT(IN)    :: MY 
C 
      REAL, INTENT(IN)       :: BEFIN 
      REAL, INTENT(OUT)      :: AC_DUM 
      REAL, INTENT(IN)       :: TCF 
C 
      INTEGER I 
      INTEGER IPI 
C 
      REAL ACCF 
      REAL AVSP 
      REAL BEF 
      REAL BEFF 
      REAL BEFAC 
      REAL BEFL 
      REAL BEFR 
      REAL BEFRC 
      REAL BEF_SPD 
      REAL BEFT2 
      REAL BEFT2N 
      REAL OFFCYCLE 
      REAL PLEVEL1 
      REAL PLEVEL2 
      REAL PLEVEL3 
      REAL BEFN 
      REAL BEF_AVSP 
      REAL SFTP_SPD 
C 
C  initialize 
      AC_DUM = 0.0 
      AVSP  = 19.6 
      BEFR  =  0.0 
      BEFL  =  0.0 
      BEFF  =  0.0 
      BEFAC =  0.0 
      BEFIDL(IP,IDX,IV,IH) = 0.0 
      SPEED_CF = BEFIN 
C 
C  Tier 2 vehicles use LEV BEF for calculating speed correction 
C  factors and air conditioning effects. Store the Tier 2 emission 
C  rate for later use in BEFT2. 
C 
      BEF=BEFIN 
      BEFT2=BEFIN 
      IF(ISTD.EQ.4 .AND. BEFIN.GT.0.0000001) BEF=EMFAC(1,8,3,INH) 
C 
      IF(VVLDGAS(IV).NE.1)THEN 
C 
C  Adjustments for light duty diesel vehicles, all heavy duty vehicles 
C  and motorcycles using SCF calculated from the MOBILE5 method. 
C 
        SELECT CASE (IROAD) 
         CASE(1)       ! Freeways 
            DO I = 1, 14 
               BEFF = BEFF+BEFT2*M5SCF(I,IP)*SVMT(1,I,IH) 
            END DO 
            SPEED_CF = BEFF 
C 
         CASE(2)       ! Arterial and Collector Roadways 
            DO I = 1, 14 
               BEFAC = BEFAC+BEFT2*(M5SCF(I,IP))*SVMT(2,I,IH) 
            END DO 
            SPEED_CF = BEFAC 
C 
         CASE(3)       ! Local Roadways 
            BEFL = BEFT2*SCFL 
            SPEED_CF = BEFL 
C 
         CASE(4)        !Ramp Roadways 
            BEFR = BEFT2*SCFR 
            SPEED_CF = BEFR 
         END SELECT 
C 
C     Apply the temperature correction factor 
C 
         SPEED_CF = SPEED_CF * TCF 
C 
CDFK 
CDFK      WRITE(99,99) IV,MY,SPEED_CF,TCF 
CDFK   99 FORMAT(' SPEED CF: ',I3,I5,2F10.5) 
C 
C calculate idl bef  
C 
         BEFIDL(IP,IDX,IV,IH)=BEFT2*M5SCF(1,IP) 
C 
C     Adjustments for light duty gasoline vehicles 
C 
      ELSE    ! VVLDGAS(IV) = 1 
C 
C     Add the off-cycle emissions. 
C 
         IF(IROAD.LE.2 .AND. SFTP_SP_FLAG.EQ.1) THEN 
C     Offcycle for Freeways & Collectors 
            BEFN = BEF + OFFCYCLE(INH,IP,BEF)* 
     &                   SFTP_SPD(INH,MY,IV,IP,ISTD) 
            BEFT2N = BEFT2 + OFFCYCLE(INH,IP,BEF)* 
     &                   SFTP_SPD(INH,MY,IV,IP,ISTD) 
         ELSE 
            BEFN = MAX(BEF, 0.000001) 
            BEFT2N = MAX(BEFT2, 0.000001) 
         END IF 
C 
C Here we determine the weight factor for three levels 
C 
         IF (IP.EQ.1) THEN 
            IPI = 4 
         ELSE 
            IPI = IP 
         END IF 
C 
         IF(BEFN.LE.LEVEL(1,IPI)) THEN 
           PLEVEL1=1.0 
           PLEVEL2=0.0 
           PLEVEL3=0.0 
         ELSE 
           IF(BEFN.LE.LEVEL(2,IPI)) THEN 
             PLEVEL1=(LEVEL(2,IPI)-BEFN)/(LEVEL(2,IPI)-LEVEL(1,IPI)) 
             PLEVEL2=1.0-PLEVEL1 
             PLEVEL3=0.0 
           ELSE 
             IF(BEFN.LT.LEVEL(3,IPI)) THEN 
               PLEVEL1=0.0 
               PLEVEL2=(LEVEL(3,IPI)-BEFN)/(LEVEL(3,IPI)-LEVEL(2,IPI)) 
               PLEVEL3=1.0-PLEVEL2 
             ELSE 
               PLEVEL1=0.0 
               PLEVEL2=0.0 
               PLEVEL3=1.0 
             END IF 
           END IF 
         END IF 
C 
         SELECT CASE (IROAD) 
          CASE(1)       ! Freeways 
C 
C    Calculate the BEF for Freeways Roadway.  First add speeds from 
C    2.5 to 65. 
C 
            DO I=1,14 
C 
C    Calculate the speed/temperature corrected BEF for air conditioning 
C    calculations.  For Tier 2, this will be based on the LEV BEF and LEV-based 
C    speed correction factors. 
C 
              IF(AC_FLAG .EQ. 1) THEN    ! Calculate A/C Correction 
              BEF_SPD = BEFN*(PLEVEL1*SCFF(1,I,IPI) + 
     &          PLEVEL2*SCFF(2,I,IPI)+PLEVEL3*SCFF(3,I,IPI)) 
                BEF_SPD = BEF_SPD * TCF  ! Apply temperature correction 
                IF(BEF_SPD.GT.0..OR.MY.GE.2001)  
     &           AC_DUM = AC_DUM + 
     &             ACCF(IV,IP,INH,MY,IROAD,ISTD,BEF_SPD, 
     &                  AVESPD(I,IROAD)) * SVMT(1,I,IH) 
              END IF 
C    Calculate the speed/temperature corrected BEF for return to HRLOOP. 
C    For Tier 2, this will be based on the Tier 2 BEF and the LEV-based 
C    speed correction factors. 
C 
              BEF_SPD = BEFT2N*(PLEVEL1*SCFF(1,I,IPI) + 
     &        PLEVEL2*SCFF(2,I,IPI)+PLEVEL3*SCFF(3,I,IPI)) 
              BEF_SPD = BEF_SPD * TCF    ! Apply temperature correction 
              BEFF = BEFF + BEF_SPD*SVMT(1,I,IH) 
            END DO 
            SPEED_CF = BEFF 
C 
          CASE(2)       ! Arterial and Collector Roadways 
C 
C    Calculate the BEF for Arterials & Collectors Roadway. 
C    First add speeds from 2.5 to 65. 
C 
            DO I=1,14 
C 
C    Calculate the speed/temperature corrected BEF for air conditioning 
C    calculations.  For Tier 2, this will be based on the LEV BEF and LEV-based 
C    speed correction factors. 
C 
              IF(AC_FLAG .EQ. 1) THEN    ! Calculate A/C Correction 
              BEF_SPD = BEFN*(PLEVEL1*SCFAC(1,I,IPI) + 
     &          PLEVEL2*SCFAC(2,I,IPI)+PLEVEL3*SCFAC(3,I,IPI)) 
                BEF_SPD = BEF_SPD * TCF  ! Apply temperature correction 
                IF(BEF_SPD.GT.0..OR.MY.GE.2001)  
     &          AC_DUM = AC_DUM + 
     &             ACCF(IV,IP,INH,MY,IROAD,ISTD,BEF_SPD, 
     &                  AVESPD(I,IROAD)) * SVMT(2,I,IH) 
              END IF 
C    Calculate the speed/temperature corrected BEF for return to HRLOOP. 
C    For Tier 2, this will be based on the Tier 2 BEF and the LEV-based 
C    speed correction factors. 
C 
              BEF_SPD = BEFT2N*(PLEVEL1*SCFAC(1,I,IPI) + 
     &        PLEVEL2*SCFAC(2,I,IPI)+PLEVEL3*SCFAC(3,I,IPI)) 
              BEF_SPD = BEF_SPD * TCF    ! Apply temperature correction 
              BEFAC = BEFAC + BEF_SPD*SVMT(2,I,IH) 
            END DO 
            SPEED_CF = BEFAC 
C 
          CASE(3)       ! Local Roadways 
C 
C  Here We calculate the BEF for Locals Roadway 
C 
C 
C    Calculate the speed/temperature corrected BEF for air conditioning 
C    calculations.  For Tier 2, this will be based on the LEV BEF and LEV-based 
C    speed correction factors. 
C 
            IF(AC_FLAG .EQ. 1) THEN    ! Calculate A/C Correction 
              BEF_AVSP = BEF*AVSP 
              BEFL = (LOCAL(1,IPI)+(BEF_AVSP)*(LOCAL(2,IPI) + 
     *               LOCAL(3,IPI)*(BEF_AVSP)))/AVSPL 
              BEFL = BEFL * TCF        ! Apply temperature correction 
                IF(BEFL.GT.0..OR.MY.GE.2001)  
     &         AC_DUM = AC_DUM + 
     &                  ACCF(IV,IP,INH,MY,IROAD,ISTD,BEFL,AVSPL) 
            END IF 
C    Calculate the speed/temperature corrected BEF for return to HRLOOP. 
C    For Tier 2, this will be based on the Tier 2 BEF and the LEV-based 
C    speed correction factors. 
C 
            BEF_AVSP = BEFT2*AVSP 
            BEFL = (LOCAL(1,IPI)+(BEF_AVSP)*(LOCAL(2,IPI) + 
     *             LOCAL(3,IPI)*(BEF_AVSP)))/AVSPL 
            BEFL = BEFL * TCF        ! Apply temperature correction 
            SPEED_CF = BEFL 
C 
          CASE(4)        !Ramp Roadways 
C 
C  Here We calculate the BEF for Ramps Roadway 
C 
C    Calculate the BEF emissions in grams/hour. 
C    NOx is capped at 6 g/mi. 
C 
            IF(IP .EQ. 3 .AND. BEF .GT. 6.0) THEN 
               BEF_AVSP = 6.0*AVSP 
            ELSE 
               BEF_AVSP = BEF*AVSP 
            END IF 
C 
C      SFTP reduction applied to difference between ramp emissions 
C      and LA4 emissions.  First compute offcycle increase. 
C 
            BEFRC = ((RAMP(1,IPI)+(BEF_AVSP)*(RAMP(2,IPI) + 
     &                  RAMP(3,IPI)*(BEF_AVSP)))/AVSPR) - BEF 
C 
C    Calculate the speed/temperature corrected BEF for air conditioning 
C    calculations.  For Tier 2, this will be based on the LEV BEF and LEV-based 
C    speed correction factors (BEFRC). 
C 
            IF(AC_FLAG .EQ. 1) THEN    ! Calculate A/C Correction 
C 
C    Apply SFTP reduction to offcycle increase and add to BEF 
C    emissions. 
C 
            BEFR = BEF + BEFRC*SFTP_SPD(INH,MY,IV,IP,ISTD) 
            BEFR = BEFR * TCF          ! Apply temperature correction 
                IF(BEFR.GT.0..OR.MY.GE.2001)  
     &      AC_DUM = AC_DUM + 
     &                 ACCF(IV,IP,INH,MY,IROAD,ISTD,BEFR,AVSPR) 
            END IF 
C 
C    Calculate the speed/temperature corrected BEF for return to HRLOOP. 
C    For Tier 2, this will be based on the Tier 2 BEF and the LEV-based 
C    speed correction factors. 
C 
C    Apply SFTP reduction to offcycle increase and add to BEF 
C    emissions.  For non-Tier 2 vehicles, BEF and BEFT2 will be the same. 
C 
            BEFR = BEFT2 + BEFRC*SFTP_SPD(INH,MY,IV,IP,ISTD) 
            BEFR = BEFR * TCF          ! Apply temperature correction 
            SPEED_CF = BEFR 
C 
         END SELECT 
C 
         SPEED_CF = MAX(SPEED_CF, 0.0) 
C 
C   Calculate idl bef 
C 
         BEFIDL(IP,IDX,IV,IH)=BEFT2N*(PLEVEL1*SCFF(1,1,IPI) + 
     &        PLEVEL2*SCFF(2,1,IPI) + PLEVEL3*SCFF(3,1,IPI)) 
C 
      END IF 
C 
      BEFIDL(IP,IDX,IV,IH)=BEFIDL(IP,IDX,IV,IH)*HVMT(IH) 
C 
      RETURN 
      END 
