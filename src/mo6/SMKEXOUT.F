      SUBROUTINE SMKEXOUT(HEF,IP,IV,JDX)
      
C  SMKEXOUT aggregates exhaust running and start emissions by vehicle age and
C  vehicle category (8 types).  It then writes the emmisions to the database
C  file.

      USE DATABASE, ONLY : DBSELPOL, DBSELVEH, DBSELEFT, 
     &                     DBSCTITLE, TABCHAR, DBUNIT,
     &                     DBSELFAC
#ifdef M6LIB
      USE MODMBSET, ONLY : EMISSIONS
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
      
      REAL,DIMENSION(24,5),INTENT(IN)  ::  HEF  ! gm/mile emission factors
      INTEGER,INTENT(IN)               ::  IP   ! pollutant
      INTEGER,INTENT(IN)               ::  IV   ! vehicle type (1 - 28)
      INTEGER,INTENT(IN)               ::  JDX  ! age
      
      REAL, ALLOCATABLE, SAVE :: EXHAUST (:,:,:)  ! aggregated exhaust emission factors
                                                  ! dimensions: (IH,IFAC,IV)
      REAL, ALLOCATABLE, SAVE :: ZMILES (:)       ! summed miles
                                                  ! dimensions: (IH)
      REAL, ALLOCATABLE, SAVE :: ZGMHOUR (:,:)    ! summed gm/hour
                                                  ! dimensions: (IH,IFAC)
      
      INTEGER :: IH        ! hour
      INTEGER :: IFAC      ! facility type
      INTEGER :: IVEH      ! vehicle type (1 - 8)
      INTEGER :: IOS       ! I/O status
      INTEGER :: EFTYPE    ! emission factor type
      INTEGER :: EMISPOS   ! position in master emission factor array
      
      REAL :: MILES        ! hourly miles
      REAL :: SUMVMT       ! vmt fraction for one of 8 vehicle types
      REAL :: ZGMMILE      ! summed gm/mile
      
      LOGICAL, SAVE :: INITIAL = .TRUE.  ! true: first time through subroutine
      
      CHARACTER(LEN=300) :: MESG                   ! message buffer
      CHARACTER(LEN=16)  :: PROGNAME = 'SMKEXOUT'  ! program name
      
C  First time, initialize emission factor array
      IF (INITIAL) THEN
      	 ALLOCATE( EXHAUST( MAXIH,5,8 ), STAT=IOS )
      	 CALL CHECKMEM( IOS, 'EXHAUST', PROGNAME )
      	 ALLOCATE( ZMILES ( MAXIH     ), STAT=IOS )
      	 CALL CHECKMEM( IOS, 'ZMILES', PROGNAME )
      	 ALLOCATE( ZGMHOUR( MAXIH,5   ), STAT=IOS )
      	 CALL CHECKMEM( IOS, 'ZGMHOUR', PROGNAME )
      	
         EXHAUST = 0.
         ZMILES  = 0.
         ZGMHOUR = 0.

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
      CASE (14)
         IVEH = 5
         SUMVMT = VMTMIX(14)
      CASE (15,28)
         IVEH = 6
         SUMVMT = VMTLDDT
      CASE (16:23,26,27)
         IVEH = 7
         SUMVMT = VMTHDD
      CASE (24)
         IVEH = 8
         SUMVMT = VMTMIX(24)
      CASE DEFAULT
         MESG = 'INTERNAL ERROR: Unexpected vehicle type'
         CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
      END SELECT

      DO IH = 1,24
         MILES = BMYMPD(JDX,IV)*HVMT(IH)*MEVMYR(JDX,IV)
         ZMILES(IH) = ZMILES(IH) + MILES
      
         DO IFAC = 1,5
            ZGMHOUR(IH,IFAC) = ZGMHOUR(IH,IFAC)+HEF(IH,IFAC)*MILES
  
C.....      If age = 0, we've finished summing and can get the aggregated
C           gm/mile over ages, then aggregate by vehicle type
            IF (JDX-1 == 0) THEN
               IF (ZMILES(IH) < 0.000001) THEN
                  ZGMMILE = 0.0
               ELSE
                  ZGMMILE = ZGMHOUR(IH,IFAC)/ZMILES(IH)
               END IF
            
               EXHAUST(IH,IFAC,IVEH) = EXHAUST(IH,IFAC,IVEH) + 
     &            (ZGMMILE*VMTMIX(IV))/SUMVMT    
     
C.....         Reinitialize arrays after last hour and facility
               IF (IH == 24 .AND. IFAC == 5) THEN
                  ZMILES  = 0.
                  ZGMHOUR = 0.
               END IF
            END IF
         END DO
      END DO

C  If vehicle type = 28 and age = 0, we're done with this pollutant, so output array
      IF (IV == 28 .AND. JDX-1 == 0) THEN
      	 DO EFTYPE = 1,2
      	 
C.....      Check if this is a valid pollutant/process combo
            IF (M6POL2EF( EFTYPE, IP ) == -1) CYCLE
         
            DO IVEH = 1,8
            
C.....         Check if this is a valid vehicle/process combo
               IF (M6VEH2EF( EFTYPE, IVEH ) == -1) CYCLE
                           
               DO IFAC = 1,5

C.....             Check if this is a valid facility/process combo
                   IF (M6FAC2EF( EFTYPE, IFAC ) == -1) CYCLE
#ifdef M6LIB     
                   DO IH = 1,24
                      EMISSIONS( EFTYPE )%PTR( STR2INT( DBSCTITLE ),
     &                                         M6POL2EF( EFTYPE, IP ),
     &                                         M6VEH2EF( EFTYPE, IVEH ),
     &                                         M6FAC2EF( EFTYPE, IFAC ),
     &                                         IH ) = 
     &                EXHAUST( IH,IFAC,IVEH )
                   END DO
#else
                   WRITE(DBUNIT(1),110) (DBSCTITLE,TABCHAR,IP,TABCHAR,
     &                  IVEH,TABCHAR,EFTYPE,TABCHAR,IFAC,TABCHAR,IH,
     &                  TABCHAR,EXHAUST(IH,IFAC,IVEH),IH=1,24)
  110 FORMAT(A10,A1,I1,A1,I2,A1,I1,A1,I1,A1,I2,A1,E12.5)
#endif
               END DO  ! facility loop
            END DO  ! vehicle loop
         END DO  ! ef type loop
         
C.....   Reinitialize emission factors array         
         EXHAUST = 0.
      END IF
      
      END SUBROUTINE SMKEXOUT
      
      
