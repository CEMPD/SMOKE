      SUBROUTINE SMKEVPOTHER(IEVAP,IV,JDX,GRAMS_MI,IROAD)
      
C  SMKEVPOTHER aggregates evaporative emissions (except running loss) by 
C  vehicle age and vehicle category.  It then writes the emissions 
C  to the database file.

      USE DATABASE, ONLY : DBSELPOL, DBSELVEH, DBSELEFT, DBSCTITLE,
     &                     TABCHAR, DBUNIT, DBSELFAC                
#ifdef M6LIB
      USE MODEMFAC, ONLY : EMISSIONS, NVTYPE, SMKVEH2EF
#endif
      
      IMPLICIT NONE

      INCLUDE 'CEVBMY.I' ! BMYMPD    
      INCLUDE 'MAXIMA.I' ! MAXIH    
      INCLUDE 'MYRCAL.I' ! MEVMYR
      INCLUDE 'SPEED9.I' ! HVMT
      INCLUDE 'VMXCOM.I' ! VMTMIX,VMTGT12,VMTGT34,VMTHDG,VMTLDDT,VMTHDD
      
      INCLUDE 'M6CNST3.EXT'   !  Mobile6 constants
      
      INTEGER STR2INT
      INTEGER EFPOSITION
      
      EXTERNAL STR2INT, EFPOSITION
      
      INTEGER,INTENT(IN) ::  IEVAP         ! evap emission process
      INTEGER,INTENT(IN) ::  IV            ! vehicle type
      INTEGER,INTENT(IN) ::  JDX           ! age
      INTEGER,INTENT(IN) ::  IROAD         ! road type
      REAL, INTENT(IN)   ::  GRAMS_MI(24)  ! gm/mile emission factors
            
      REAL, ALLOCATABLE, SAVE :: EVPOTHER (:,:,:) ! aggregated evap efs
                                                  ! dimensions: (IH,IVEH,IEVAP)
      REAL, ALLOCATABLE, SAVE :: ZMILES   (:,:)   ! summed miles
                                                  ! dimensions: (IH,IEVAP)
      REAL, ALLOCATABLE, SAVE :: ZGMHOUR  (:,:)   ! summed gm/hour
                                                  ! dimensions: (IH,IEVAP)
                                                  
      INTEGER :: IP = 1        ! pollutant will always be HC
      INTEGER :: IFAC = 5      ! road type will always be none
      
      INTEGER :: IH            ! hour
      INTEGER :: IOS           ! I/O status
      INTEGER :: IVEH          ! vehicle type
      INTEGER :: EFTYPE, EFX   ! emission factor type and index
      INTEGER :: EMISPOS       ! position in master emission factor array
      
      REAL :: MILES            ! hourly miles
      REAL :: SUMVMT           ! vmt fraction for one of 8 vehicle types
      REAL :: ZGMMILE          ! summed gm/mile
      
      LOGICAL, SAVE :: INITIAL = .TRUE.   ! true: first time through subroutine
      
      CHARACTER(300) :: MESG                      ! message buffer
      CHARACTER(16)  :: PROGNAME = 'SMKEVPOTHER'  ! program name
      
C  First time, initialize emission factor array
      IF (INITIAL) THEN
         ALLOCATE( EVPOTHER( MAXIH,NVTYPE,5 ), STAT=IOS )
         CALL CHECKMEM( IOS, 'EVPOTHER', PROGNAME )
         ALLOCATE( ZMILES  ( MAXIH,5   ), STAT=IOS )
         CALL CHECKMEM( IOS, 'ZMILES', PROGNAME )
         ALLOCATE( ZGMHOUR ( MAXIH,5   ), STAT=IOS )
         CALL CHECKMEM( IOS, 'ZGMHOUR', PROGNAME )
 
         EVPOTHER = 0.
         ZMILES   = 0.
         ZGMHOUR  = 0.
         
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
         MESG = 'INTERNAL ERROR: Unexpected emission process'
         CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
      END SELECT

C  Skip vehicle types not used for evaporative emissions
      IF (IV > 13 .AND. IV /= 24 .AND. IV /= 25) RETURN
      
C  Map MOBILE6 vehicle type to SMOKE vehicle type
      CALL SMKCALCVMT( IV, IVEH, SUMVMT )

      DO IH = 1,24
         MILES = BMYMPD(JDX,IV)*HVMT(IH)*MEVMYR(JDX,IV)
         ZMILES(IH,EFX) = ZMILES(IH,EFX) + MILES

         ZGMHOUR(IH,EFX) = ZGMHOUR(IH,EFX) + GRAMS_MI(IH)*MILES
           
C.....   If age = 0, we've finished summing and can get the aggregated
C        gm/mile over ages, then aggregate by vehicle type
         IF (JDX-1 == 0) THEN
            IF (ZMILES(IH,EFX) < 0.000001) THEN
               ZGMMILE = 0.0
            ELSE
               ZGMMILE = ZGMHOUR(IH,EFX)/ZMILES(IH,EFX)
            END IF

            EVPOTHER(IH,IVEH,EFX) = EVPOTHER(IH,IVEH,EFX) + 
     &            (ZGMMILE*VMTMIX(IV))/SUMVMT

C.....      Reinitialize arrays after last hour and emission type
            IF (IH == 24 .AND. EFX == 5) THEN
               ZMILES  = 0.
               ZGMHOUR = 0.
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

C....       Check if this is a valid pollutant/process combo
            IF (M6POL2EF( EFTYPE, IP ) == -1) CYCLE

C....       Check if this is a valid facility/process combo
            IF (M6FAC2EF( EFTYPE, IFAC ) == -1) CYCLE
            
            DO IVEH = 1,NVTYPE

C....          Check if this is a valid vehicle/process combo
               IF (SMKVEH2EF( EFTYPE, IVEH ) == -1) CYCLE

#ifdef M6LIB  
                  DO IH = 1,24
                     EMISSIONS( EFTYPE )%PTR( 
     &                   STR2INT( DBSCTITLE ),
     &                   M6POL2EF( EFTYPE, IP ),
     &                   SMKVEH2EF( EFTYPE, IVEH ),
     &                   M6FAC2EF( EFTYPE, IFAC ), IH ) = 
     &               EVPOTHER( IH, IVEH, EFX )
                  END DO
#else
                  WRITE(DBUNIT(1),110) (DBSCTITLE,TABCHAR,IP,
     &               TABCHAR,IVEH,TABCHAR,EFTYPE,TABCHAR,IFAC,
     &               TABCHAR,IH,TABCHAR,EVPOTHER(IH,IVEH,EFX),
     &               IH=1,24)
  110 FORMAT(A10,A1,I1,A1,I2,A1,I1,A1,I1,A1,I2,A1,E12.5)
#endif  

            END DO
         END DO
         
C.....   Reinitialize emission factors array         
         EVPOTHER = 0.
      END IF
      
      END SUBROUTINE SMKEVPOTHER
