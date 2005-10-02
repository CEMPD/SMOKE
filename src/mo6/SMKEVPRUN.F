      SUBROUTINE SMKEVPRUN(IEVAP,IV,JDX,GRAMS_MI,IROAD)
      
C  SMKEVPRUN aggregates evaporative running emissions by vehicle age and
C  vehicle category.  It then writes the emissions to the database
C  file.

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
      INTEGER,INTENT(IN) ::  IV            ! vehicle type (1 - 28)
      INTEGER,INTENT(IN) ::  JDX           ! age
      INTEGER,INTENT(IN) ::  IROAD         ! road type
      REAL, INTENT(IN)   ::  GRAMS_MI(24)  ! gm/mile emission factors

      REAL, ALLOCATABLE, SAVE :: EVPRUN  (:,:,:)  ! aggregated evap running efs
                                                  ! dimensions: (IH,IFAC,IVEH)
      REAL, ALLOCATABLE, SAVE :: ZMILES  (:,:)    ! summed miles
                                                  ! dimensions: (IH,IFAC)
      REAL, ALLOCATABLE, SAVE :: ZGMHOUR (:,:)    ! summed gm/hour
                                                  ! dimensions: (IH,IFAC)
      
      INTEGER :: IP = 1    ! pollutant will always be HC
      INTEGER :: IH        ! hour
      INTEGER :: IFAC      ! facility type
      INTEGER :: IVEH      ! vehicle type
      INTEGER :: IOS       ! I/O status
      INTEGER :: EFTYPE    ! emission factor type
      INTEGER :: EMISPOS   ! position in master emission factor array
      
      REAL :: MILES        ! hourly miles
      REAL :: SUMVMT       ! vmt fraction for one of 8 vehicle types
      REAL :: ZGMMILE      ! summed gm/mile
      
      LOGICAL, SAVE :: INITIAL = .TRUE.  ! true: first time through subroutine
      
      CHARACTER(300) :: MESG                    ! message buffer
      CHARACTER(16)  :: PROGNAME = 'SMKEVPRUN'  ! program name
      
C  First time, initialize emission factor array
      IF (INITIAL) THEN
         ALLOCATE( EVPRUN ( MAXIH,4,NVTYPE ), STAT=IOS )
         CALL CHECKMEM( IOS, 'EVPRUN', PROGNAME )
         ALLOCATE( ZMILES ( MAXIH,4   ), STAT=IOS )
         CALL CHECKMEM( IOS, 'ZMILES', PROGNAME )
         ALLOCATE( ZGMHOUR( MAXIH,4   ), STAT=IOS )
         CALL CHECKMEM( IOS, 'ZGMHOUR', PROGNAME )
 
         EVPRUN  = 0.
         ZMILES  = 0.
         ZGMHOUR = 0.
         
         INITIAL = .FALSE.
      END IF

C  Skip vehicle types not used for evaporative running emissions
      IF (IV > 13 .AND. IV /= 24 .AND. IV /= 25 ) RETURN
      
C  Map MOBILE6 vehicle type to SMOKE vehicle type
      CALL SMKCALCVMT( IV, IVEH, SUMVMT )

      DO IH = 1,24
         MILES = BMYMPD(JDX,IV)*HVMT(IH)*MEVMYR(JDX,IV)
         ZMILES(IH,IROAD) = ZMILES(IH,IROAD) + MILES

         ZGMHOUR(IH,IROAD) = ZGMHOUR(IH,IROAD) + 
     &                           GRAMS_MI(IH)*MILES
           
C.....   If age = 0, we've finished summing and can get the aggregated
C        gm/mile over ages, then aggregate by vehicle type
         IF (JDX-1 == 0) THEN
            IF (ZMILES(IH,IROAD) < 0.000001) THEN
               ZGMMILE = 0.0
            ELSE
               ZGMMILE = ZGMHOUR(IH,IROAD)/ZMILES(IH,IROAD)
            END IF

            EVPRUN(IH,IROAD,IVEH) = EVPRUN(IH,IROAD,IVEH) + 
     &            (ZGMMILE*VMTMIX(IV))/SUMVMT    

C.....      Reinitialize arrays after last hour and facility
            IF (IH == 24 .AND. IROAD == 4) THEN
               ZMILES  = 0.
               ZGMHOUR = 0.
            END IF
         END IF
      END DO

C  If vehicle = 25, age = 0, and road type = 4, we're done, so output array
      IF (IV == 25 .AND. JDX-1 == 0 .AND. IROAD == 4) THEN
         EFTYPE = IEVAP

C....    Check if this is a valid pollutant/process combo
         IF (M6POL2EF( EFTYPE, IP ) /= -1) THEN
               
            DO IVEH = 1,NVTYPE
            
C....          Check if this is a valid vehicle/process combo
               IF (SMKVEH2EF( EFTYPE, IVEH ) == -1) CYCLE
            
               DO IFAC = 1,4
               
C....             Check if this is a valid facility/process combo
                  IF (M6FAC2EF( EFTYPE, IFAC ) == -1) CYCLE
#ifdef M6LIB
                  DO IH = 1,24
                     EMISSIONS( EFTYPE )%PTR( STR2INT( DBSCTITLE ),
     &                                        M6POL2EF( EFTYPE, IP ),
     &                                        SMKVEH2EF( EFTYPE, IVEH ),
     &                                        M6FAC2EF( EFTYPE, IFAC ), 
     &                                        IH ) = 
     &               EVPRUN( IH, IFAC, IVEH )
                  END DO  ! hour loop
#else         
                  WRITE(DBUNIT(1),110) (DBSCTITLE,TABCHAR,IP,
     &            TABCHAR,IVEH,TABCHAR,EFTYPE,TABCHAR,IFAC,
     &            TABCHAR,IH,TABCHAR,EVPRUN(IH,IFAC,IVEH),IH=1,24)
  110 FORMAT(A10,A1,I1,A1,I2,A1,I1,A1,I1,A1,I2,A1,E12.5)
#endif   
              
               END DO  ! facility loop
            END DO  ! vehicle loop
         END IF
         
C.....   Reinitialize emission factors array         
         EVPRUN = 0.
      END IF
      
      END SUBROUTINE SMKEVPRUN
      
      
