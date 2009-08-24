      SUBROUTINE SMKATOUT
       
C  SMKATOUT aggregates air toxics emissions by vehicle age and 
C  vehicle category.  It then writes the emissions to the
C  database file.

      USE DATABASE, ONLY : DBSELTOX,DBSELEFT,DBSELFAC, 
     &                     DBUNIT,TABCHAR,DBSCTITLE 
      USE ATHEAP, ONLY : DBAT
      
#ifdef M6LIB      
      USE MODEMFAC, ONLY : EMISSIONS, NVTYPE, SMKVEH2EF
#endif
     
      IMPLICIT NONE

      INCLUDE 'ATOX1.I'  ! MAXIAT     
      INCLUDE 'CEVBMY.I' ! BMYMPD 
      INCLUDE 'IVTYPE.I' ! VVDSL,VVHEAVY 
      INCLUDE 'MAXIMA.I' ! MAXYRS 
      INCLUDE 'MYRCAL.I' ! MEVMYR
      INCLUDE 'PART1.I'  ! MAXIPPM
      INCLUDE 'SPEED9.I' ! HVMT
      INCLUDE 'VMXCOM.I' ! VMTMIX,VMTGT12,VMTGT34,VMTHDG,VMTLDDT,VMTHDD

      INCLUDE 'M6CNST3.EXT'   !  Mobile6 constants
        
      INTEGER STR2INT
      
      EXTERNAL STR2INT

      INTEGER  ::  I        ! loop variable
      INTEGER  ::  ETYPE    ! emission factor type 
      INTEGER  ::  IXV      ! eXhaust/eVap index 
      INTEGER  ::  IFAC     ! facility type
      INTEGER  ::  STFAC    ! starting facility type
      INTEGER  ::  ENDFAC   ! ending facility type
      INTEGER  ::  IH       ! hour of day
      INTEGER  ::  IAT      ! air toxic pollutant number (1 - MAXIAT)
      INTEGER  ::  IP       ! actual pollutant number (16 - 21)
      INTEGER  ::  JDX      ! age
      INTEGER  ::  IDX      ! year
      INTEGER  ::  IV       ! vehicle type (1 - 28)
      INTEGER  ::  IVEH     ! SMOKE vehicle class
      INTEGER  ::  SMKPOL   ! SMOKE pollutant number based on M6 index
      INTEGER  ::  IOS      ! I/O status
      
      REAL              ::  SUMVMT   ! vmt fraction for one of SMOKE vehicle types 
      REAL              ::  MILES    ! hourly miles                            
      REAL              ::  ZMILES   ! used to sum hourly miles                
      REAL              ::  ZGMHOUR  ! used to sum gm/hour                     
      REAL              ::  ZGMMILE  ! used to calculate hourly gm/mile
      
      REAL, ALLOCATABLE, SAVE ::  ATAGGR(:) ! fully aggregated AT output
                                      ! dimensions: (IVEH)

      LOGICAL, SAVE :: INITIAL = .TRUE.   ! true: first time through subroutine

      CHARACTER(300) :: MESG                   ! message buffer
      CHARACTER(16)  :: PROGNAME = 'SMKATOUT'  ! program name

C  Begin body of subroutine SMKATOUT

C  Allocate array first time
      IF (INITIAL) THEN
          ALLOCATE( ATAGGR( NVTYPE ), STAT=IOS )
          CALL CHECKMEM( IOS, 'ATAGGR', PROGNAME )
          
          INITIAL = .FALSE.
      END IF

C  Loop over air toxic pollutants
      DO IAT = 1,MAXIAT

C  Check that we need this pollutant
         IF (DBSELTOX(IAT) == 1) CYCLE
         IP = IAT + MAXIPPM + 4  ! AT pollutants start after classic and PM pollutants 

C  Find SMOKE pollutant number based on M6 number
         SMKPOL = 0
         DO I = 1,MXM6POLS
             IF( IP == SMK2M6POL( I ) ) THEN
                 SMKPOL = I
                 EXIT
             END IF
         END DO

         IF( SMKPOL == 0 ) THEN
             MESG = 'INTERNAL ERROR: Unexpected pollutant'
             CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
         END IF

C  Loop over emission process types
         DO ETYPE = 1,8

C  Check that we need this process type
            IF (DBSELEFT(ETYPE) == 1) CYCLE

C  Only output evap emissions for Benzene and MTBE
            IF (IAT >= 3 .AND. ETYPE >= 3) EXIT 

C  Check if this is a valid pollutant/process combo
            IF (M6POL2EF( ETYPE, SMKPOL ) == -1) CYCLE
            
C  Map emission factor type to air toxic EF type
            IF (ETYPE <= 2) IXV=1  ! Exhaust Running and Starts 
            IF (ETYPE == 3) IXV=2  ! Evap Hot Soak 
            IF (ETYPE == 4) IXV=3  ! Evap Diurnal 
            IF (ETYPE == 5) IXV=5  ! Evap Resting 
            IF (ETYPE == 6) IXV=4  ! Evap Running 
            IF (ETYPE == 7) CYCLE  ! Evap Crankcase (Not Available) 
            IF (ETYPE == 8) IXV=6  ! Evap Refueling 

C  Loop over facility types
            DO IFAC = 1,MXM6FACS

C  Check that this is a valid facility/process combo
               IF (M6FAC2EF( ETYPE, IFAC ) == -1) CYCLE

C  Check that we need this facility type
               IF (DBSELFAC(IFAC) /= 2) CYCLE
                     
               DO IH = 1,24
C  Reinitialize output array
                  ATAGGR = 0.
            
                  DO IV = 1,MAXVEH 

C  Skip unused vehicle, pollutant, and emission factor combinations

C.....               MTBE is output for gas vehicles only
                     IF (VVDSL(IV) == 1 .AND. IAT == 2) CYCLE   
                     
C.....               Heavy duty vehicles and buses do not have engine start emissions
                     IF (VVHEAVY(IV) == 1 .AND. ETYPE == 2) CYCLE
                     
C.....               Diesels do not have evaporative emissions
                     IF (VVDSL(IV) == 1 .AND. ETYPE >= 3) CYCLE
                                 
C  Map MOBILE6 vehicle type to SMOKE vehicle type  
                     CALL SMKCALCVMT( IV, IVEH, SUMVMT )
      
                     ZMILES  = 0.
                     ZGMHOUR = 0.
                     ZGMMILE = 0.
                  
                     DO JDX = 1,MAXYRS
                        IDX = (MAXYRS+1) - JDX
                        MILES = BMYMPD(JDX,IV)*HVMT(IH)*MEVMYR(JDX,IV)
                        ZMILES = ZMILES + MILES
                        ZGMHOUR = ZGMHOUR + 
     &                            DBAT(IDX,IH,IFAC,IXV,IAT,IV)*MILES
                     END DO  ! loop over years
               
                     IF (ZMILES < 0.000001 ) THEN
                        ZGMMILE = 0.0
                     ELSE
                        ZGMMILE = ZGMHOUR/ZMILES
                     END IF
 
                     ATAGGR(IVEH) = ATAGGR(IVEH) + 
     &                              (ZGMMILE*VMTMIX(IV))/SUMVMT
 
                  END DO  ! loop over 28 vehicle types

#ifdef M6LIB
                  DO IVEH = 1,NVTYPE
                  
C  Check if this is a valid process/vehicle combo
                      IF( SMKVEH2EF( ETYPE, IVEH ) == -1 ) CYCLE
                                        
                      EMISSIONS( ETYPE )%PTR( STR2INT( DBSCTITLE), 
     &                                        M6POL2EF( ETYPE, SMKPOL ), 
     &                                        SMKVEH2EF( ETYPE, IVEH ), 
     &                                        M6FAC2EF( ETYPE, IFAC ), 
     &                                        IH ) = ATAGGR( IVEH )
                  END DO  ! loop over SMOKE vehicle types
#else
                  WRITE(DBUNIT(1),110) (DBSCTITLE,TABCHAR,IP, 
     &               TABCHAR,IVEH,TABCHAR,ETYPE,TABCHAR,IFAC,
     &               TABCHAR,IH,TABCHAR,ATAGGR(IVEH),IVEH=1,NVTYPE)
  110 FORMAT(A10,A1,I2,A1,I2,A1,I2,A1,I1,A1,I2,A1,E12.5)
#endif

               END DO  ! loop over hours
            END DO  ! loop over facilities
         END DO  ! loop over emission factor types
      END DO  ! loop over pollutants

      END SUBROUTINE SMKATOUT 
      
