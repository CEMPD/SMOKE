      SUBROUTINE SMKEVPOTHER(IEVAP,IV,JDX,GRAMS_MI,IROAD)
      
C  SMKEVPOTHER aggregates evaporative emissions (except running loss) by 
C  vehicle age and vehicle category (8 types).  It then writes the emmisions 
C  to the database file.

      USE DATABASE, ONLY : DBSELPOL, DBSELVEH, DBSELEFT, DBSCTITLE,
     &                     TABCHAR, DBUNIT, DBSELFAC                
#ifdef M6LIB
      USE MODEMFAC, ONLY : EMISSIONS
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
                                                  ! dimensions: (IH,IV,IEVAP)
      REAL, ALLOCATABLE, SAVE :: ZMILES   (:,:)   ! summed miles
                                                  ! dimensions: (IH,IEVAP)
      REAL, ALLOCATABLE, SAVE :: ZGMHOUR  (:,:)   ! summed gm/hour
                                                  ! dimensions: (IH,IEVAP)
                                                  
      INTEGER :: IP = 1        ! pollutant will always be HC
      INTEGER :: IFAC = 5      ! road type will always be none
      
      INTEGER :: IH            ! hour
      INTEGER :: IOS           ! I/O status
      INTEGER :: VTYPE, IVX    ! vehicle type and index
      INTEGER :: EFTYPE, EFX   ! emission factor type and index
      INTEGER :: EMISPOS       ! position in master emission factor array
      
      REAL :: MILES            ! hourly miles
      REAL :: SUMVMT           ! vmt fraction for one of 8 vehicle types
      REAL :: ZGMMILE          ! summed gm/mile
      
      LOGICAL, SAVE :: INITIAL = .TRUE.   ! true: first time through subroutine
      
      CHARACTER(LEN=300) :: MESG                      ! message buffer
      CHARACTER(LEN=16)  :: PROGNAME = 'SMKEVPOTHER'  ! program name
      
C  First time, initialize emission factor array
      IF (INITIAL) THEN
      	 ALLOCATE( EVPOTHER( MAXIH,5,5 ), STAT=IOS )
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
         MESG = 'INTERNAL ERROR: Unexpected vehicle type'
         CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
      END SELECT

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

            EVPOTHER(IH,IVX,EFX) = EVPOTHER(IH,IVX,EFX) + 
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
            
            DO IVX = 1,5
            
C....          Convert vehicle index back to vehicle class
               IF (IVX == 5) THEN
                  VTYPE = 8
               ELSE
               	  VTYPE = IVX
               END IF

C....          Check if this is a valid vehicle/process combo
               IF (M6VEH2EF( EFTYPE, VTYPE ) == -1) CYCLE

#ifdef M6LIB  
                  DO IH = 1,24
                     EMISSIONS( EFTYPE )%PTR( 
     &                   STR2INT( DBSCTITLE ),
     &                   M6POL2EF( EFTYPE, IP ),
     &                   M6VEH2EF( EFTYPE, VTYPE ),
     &                   M6FAC2EF( EFTYPE, IFAC ), IH ) = 
     &               EVPOTHER( IH, IVX, EFX )
                  END DO
#else
                  WRITE(DBUNIT(1),110) (DBSCTITLE,TABCHAR,IP,
     &               TABCHAR,VTYPE,TABCHAR,EFTYPE,TABCHAR,IFAC,
     &               TABCHAR,IH,TABCHAR,EVPOTHER(IH,IVX,EFX),
     &               IH=1,24)
  110 FORMAT(A10,A1,I1,A1,I2,A1,I1,A1,I1,A1,I2,A1,E12.5)
#endif  

            END DO
         END DO
         
C.....   Reinitialize emission factors array         
         EVPOTHER = 0.
      END IF
      
      END SUBROUTINE SMKEVPOTHER
