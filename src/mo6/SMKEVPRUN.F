      SUBROUTINE SMKEVPRUN(IEVAP,IV,JDX,GRAMS_MI,IROAD)
      
C  SMKEVPRUN aggregates evaporative running emissions by vehicle age and
C  vehicle category (8 types).  It then writes the emmisions to the database
C  file.

      USE SMKSAVE,  ONLY : EVPRUN, ZMILES_EVPRUN, ZGMHOUR_EVPRUN, 
     &                     ZGMMILE_EVPRUN
      USE DATABASE, ONLY : DBSELPOL, DBSELVEH, DBSELEFT, DBSCTITLE,
     &                     TABCHAR, DBUNIT, DBSELFAC                
#ifdef M6LIB
      USE MODMBSET, ONLY : EMISSIONS
#endif

      IMPLICIT NONE

      INCLUDE 'CEVBMY.I' ! BMYMPD      
      INCLUDE 'MYRCAL.I' ! MEVMYR
      INCLUDE 'SPEED9.I' ! HVMT
      INCLUDE 'VMXCOM.I' ! VMTMIX,VMTGT12,VMTGT34,VMTHDG,VMTLDDT,VMTHDD
      
      INTEGER STR2INT
      INTEGER EFPOSITION
      
      EXTERNAL STR2INT, EFPOSITION
      
      INTEGER,INTENT(IN) ::  IEVAP         ! evap emission process
      INTEGER,INTENT(IN) ::  IV            ! vehicle type (1 - 28)
      INTEGER,INTENT(IN) ::  JDX           ! age
      INTEGER,INTENT(IN) ::  IROAD         ! road type
      REAL, INTENT(IN)   ::  GRAMS_MI(24)  ! gm/mile emission factors
      
      INTEGER :: IP = 1    ! pollutant will always be HC
      INTEGER :: IH        ! hour
      INTEGER :: IFAC      ! facility type
      INTEGER :: IVEH      ! vehicle type (1 - 8)
      INTEGER :: EFTYPE    ! emission factor type
      INTEGER :: EMISPOS   ! position in master emission factor array
      
      REAL :: MILES        ! hourly miles
      REAL :: SUMVMT       ! vmt fraction for one of 8 vehicle types
      
      LOGICAL, SAVE :: INITIAL = .TRUE.  ! true: first time through subroutine
      
C  First time, initialize emission factor array
      IF (INITIAL) THEN
         EVPRUN = 0.

         ZMILES_EVPRUN  = 0.
         ZGMHOUR_EVPRUN = 0.
         ZGMMILE_EVPRUN = 0.
         
         INITIAL = .FALSE.
      END IF

C  Select appropriate vehicle class (8 types) for vehicle type (28)
ccs  (check if this is available as a function somewhere)
      SELECT CASE (IV)
      CASE (1)
         IVEH = 1
         SUMVMT = VMTMIX(1)
      CASE (2,3)
         IVEH = 2
         SUMVMT = VMTGT12
      CASE (4,5)
         IVEH = 3
         SUMVMT = VMTGT34
      CASE (6:13,25)
         IVEH = 4
         SUMVMT = VMTHDG

C  The remaining types are not used for evaporative running emissions.
      CASE (14)
         RETURN      
      CASE (15,28)
         RETURN 
      CASE (16:23,26,27)
         RETURN  
      CASE (24)
         RETURN 
      CASE DEFAULT
ccs  Print some error message
      END SELECT

      DO IH = 1,24
         MILES = BMYMPD(JDX,IV)*HVMT(IH)*MEVMYR(JDX,IV)
         ZMILES_EVPRUN(IH,IROAD) = ZMILES_EVPRUN(IH,IROAD) + MILES

         ZGMHOUR_EVPRUN(IH,IROAD) = ZGMHOUR_EVPRUN(IH,IROAD) + 
     &                           GRAMS_MI(IH)*MILES
           
C.....   If age = 0, we've finished summing and can get the aggregated
C        gm/mile over ages, then aggregate by vehicle type
         IF (JDX-1 == 0) THEN
            IF (ZMILES_EVPRUN(IH,IROAD) < 0.000001) THEN
               ZGMMILE_EVPRUN(IH,IROAD) = 0.0
            ELSE
               ZGMMILE_EVPRUN(IH,IROAD) = 
     &                ZGMHOUR_EVPRUN(IH,IROAD)/ZMILES_EVPRUN(IH,IROAD)
            END IF

            EVPRUN(IH,IROAD,IVEH) = EVPRUN(IH,IROAD,IVEH) + 
     &            (ZGMMILE_EVPRUN(IH,IROAD)*VMTMIX(IV))/SUMVMT    

C.....      Reinitialize arrays after last hour and facility
            IF (IH == 24 .AND. IROAD == 4) THEN
               ZMILES_EVPRUN  = 0.
               ZGMHOUR_EVPRUN = 0.
               ZGMMILE_EVPRUN = 0.
            END IF
         END IF
      END DO

C  If vehicle = 25, age = 0, and road type = 4, we're done, so output array
      IF (IV == 25 .AND. JDX-1 == 0 .AND. IROAD == 4) THEN
         DO IVEH = 1,4
            IF (DBSELPOL(IP) == 2) THEN
               IF (DBSELEFT(IEVAP) == 2) THEN
                  EFTYPE = IEVAP
                  DO IFAC = 1,4
                     IF (DBSELFAC(IFAC) == 2) THEN
#ifdef M6LIB
                        DO IH = 1,24
                           EMISPOS = EFPOSITION(STR2INT(DBSCTITLE),
     &                                          IP, IVEH, EFTYPE, IFAC, 
     &                                          IH, .TRUE. )
                           EMISSIONS(EMISPOS) = EVPRUN(IH,IFAC,IVEH)                           
                        END DO
#else
                        WRITE(DBUNIT(1),110) (DBSCTITLE,TABCHAR,IP,
     &                  TABCHAR,IVEH,TABCHAR,EFTYPE,TABCHAR,IFAC,
     &                  TABCHAR,IH,TABCHAR,EVPRUN(IH,IFAC,IVEH),IH=1,24)
  110 FORMAT(A10,A1,I1,A1,I2,A1,I1,A1,I1,A1,I2,A1,F9.4)
#endif
                     END IF
                  END DO
               END IF
            END IF
         END DO
         
C.....   Reinitialize emission factors array         
         EVPRUN = 0.
      END IF
      
      END SUBROUTINE SMKEVPRUN
      
      
