      SUBROUTINE SMKEVPOTHER(IEVAP,IV,JDX,GRAMS_MI,IROAD)
      
C  SMKEVPOTHER aggregates evaporative emissions (except running loss) by 
C  vehicle age and vehicle category (8 types).  It then writes the emmisions 
C  to the database file.

      USE SMKSAVE,  ONLY : EVPOTHER, ZMILES_EVP, ZGMHOUR_EVP, 
     &                     ZGMMILE_EVP
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
      INTEGER,INTENT(IN) ::  IV            ! vehicle type
      INTEGER,INTENT(IN) ::  JDX           ! age
      INTEGER,INTENT(IN) ::  IROAD         ! road type
      REAL, INTENT(IN)   ::  GRAMS_MI(24)  ! gm/mile emission factors
            
      INTEGER :: IP = 1        ! pollutant will always be HC
      INTEGER :: IFAC = 5      ! road type will always be none
      
      INTEGER :: IH            ! hour
      INTEGER :: VTYPE, IVX    ! vehicle type and index
      INTEGER :: EFTYPE, EFX   ! emission factor type and index
      INTEGER :: EMISPOS       ! position in master emission factor array
      
      REAL :: MILES            ! hourly miles
      REAL :: SUMVMT           ! vmt fraction for one of 8 vehicle types
      
      LOGICAL, SAVE :: INITIAL = .TRUE.   ! true: first time through subroutine
      
C  First time, initialize emission factor array
      IF (INITIAL) THEN
         EVPOTHER    = 0.

         ZMILES_EVP  = 0.
         ZGMHOUR_EVP = 0.
         ZGMMILE_EVP = 0.
         
         INITIAL = .FALSE.
      END IF

C  Convert database emission type number to indexing numbers
      SELECT CASE (IEVAP)
      CASE (3)        ! hot soak
         EFX = 1
      CASE (4)        ! diurnal
         EFX = 2
      CASE (5)        ! resting loss
         EFX = 3
      CASE (7)        ! crankcase
         EFX = 4
      CASE (8)        ! refueling
         EFX = 5
      CASE DEFAULT
ccs  Print some error message
      END SELECT

C  Select appropriate vehicle class (8 types) for vehicle type (28)
ccs  (check if this is available as a function somewhere)
      SELECT CASE (IV)
      CASE (1)
         IVX = 1
         SUMVMT = VMTMIX(1)
      CASE (2,3)
         IVX = 2
         SUMVMT = VMTGT12
      CASE (4,5)
         IVX = 3
         SUMVMT = VMTGT34
      CASE (6:13,25)
         IVX = 4
         SUMVMT = VMTHDG

C  The following types are not used for evaporative emissions.
      CASE (14)
         RETURN
      CASE (15,28)
         RETURN
      CASE (16:23,26,27)
         RETURN

C  Motorcycles are a special case, but will need to be stored in the arrays
      CASE (24)                
         IVX = 5              ! set to 5 for arrays, but use 8 for output 
         SUMVMT = VMTMIX(24)
      CASE DEFAULT
ccs  Print some error message
      END SELECT

      DO IH = 1,24
         MILES = BMYMPD(JDX,IV)*HVMT(IH)*MEVMYR(JDX,IV)
         ZMILES_EVP(IH,EFX) = ZMILES_EVP(IH,EFX) + MILES

         ZGMHOUR_EVP(IH,EFX) = ZGMHOUR_EVP(IH,EFX) + GRAMS_MI(IH)*MILES
           
C.....   If age = 0, we've finished summing and can get the aggregated
C        gm/mile over ages, then aggregate by vehicle type
         IF (JDX-1 == 0) THEN
            IF (ZMILES_EVP(IH,EFX) < 0.000001) THEN
               ZGMMILE_EVP(IH,EFX) = 0.0
            ELSE
               ZGMMILE_EVP(IH,EFX) = 
     &                ZGMHOUR_EVP(IH,EFX)/ZMILES_EVP(IH,EFX)
            END IF

            EVPOTHER(IH,IVX,EFX) = EVPOTHER(IH,IVX,EFX) + 
     &            (ZGMMILE_EVP(IH,EFX)*VMTMIX(IV))/SUMVMT

C.....      Reinitialize arrays after last hour and emission type
            IF (IH == 24 .AND. EFX == 5) THEN
               ZMILES_EVP  = 0.
               ZGMHOUR_EVP = 0.
               ZGMMILE_EVP = 0.
            END IF
         END IF
      END DO

C  If vehicle = 25, age = 0, and emission type = 5, we're done, so output array
      IF (IV == 25 .AND. JDX-1 == 0 .AND. EFX == 5) THEN
         DO EFX = 1,5
         
C.....      Convert emission index back to factor type
            SELECT CASE (EFX)
               CASE (1)
                  EFTYPE = 3
               CASE (2)
                  EFTYPE = 4 
               CASE (3)
                  EFTYPE = 5
               CASE (4)
                  EFTYPE = 7
               CASE (5)
                  EFTYPE = 8
            END SELECT
            
            DO IVX = 1,5
            
C....          Convert vehicle index back to vehicle class
               IF (IVX == 5) THEN
                  VTYPE = 8
               ELSE
               	  VTYPE = IVX
               END IF
                           
               IF (DBSELPOL(IP) == 2) THEN
                  IF (DBSELEFT(EFTYPE) == 2) THEN
                     IF (DBSELFAC(IFAC) == 2) THEN
#ifdef M6LIB  
                        DO IH = 1,24
                           EMISPOS = EFPOSITION(STR2INT(DBSCTITLE), 
     &                                          IP, VTYPE, EFTYPE, IFAC, 
     &                                          IH, .TRUE. )
                           EMISSIONS(EMISPOS) = EVPOTHER(IH,IVX,EFX)                                                       
                        END DO
#else
                        WRITE(DBUNIT(1),110) (DBSCTITLE,TABCHAR,IP,
     &                     TABCHAR,VTYPE,TABCHAR,EFTYPE,TABCHAR,IFAC,
     &                     TABCHAR,IH,TABCHAR,EVPOTHER(IH,IVX,EFX),
     &                     IH=1,24)
  110 FORMAT(A10,A1,I1,A1,I2,A1,I1,A1,I1,A1,I2,A1,F9.4)
#endif  
                     END IF
                  END IF
               END IF
            END DO
         END DO
         
C.....   Reinitialize emission factors array         
         EVPOTHER = 0.
      END IF
      
      END SUBROUTINE SMKEVPOTHER
