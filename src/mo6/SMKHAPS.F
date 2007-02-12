      SUBROUTINE SMKHAPS(IV,JDX,IP,EFIDX,IH,IFAC,XEF,OUTPUT)

C  SMKHAPS aggregates user-defined toxic emission factors by vehicle
C  age and type.

      USE DATABASE, ONLY : DBUNIT, TABCHAR, DBSCTITLE
      USE MODEMFAC, ONLY : NUMSCEN, NTOTHAPS, PMHAPS, HAPEFS, NVTYPE
      
      IMPLICIT NONE
      
      INCLUDE 'ATOX5.I'  ! NMAP
      INCLUDE 'CEVBMY.I' ! BMYMPD      
      INCLUDE 'MAXIMA.I' ! MAXIH
      INCLUDE 'MYRCAL.I' ! MEVMYR
      INCLUDE 'PART1.I'  ! PSZCUT
      INCLUDE 'SPEED9.I' ! HVMT
      INCLUDE 'VMXCOM.I' ! VMTMIX,VMTGT12,VMTGT34,VMTHDG,VMTLDDT,VMTHDD
      
      INTEGER, EXTERNAL :: STR2INT
       
      INTEGER, INTENT(IN) :: IV      ! vehicle type (1 - 28)
      INTEGER, INTENT(IN) :: JDX     ! age
      INTEGER, INTENT(IN) :: IP      ! pollutant
      INTEGER, INTENT(IN) :: EFIDX   ! emission factor type (1 - 6)
      INTEGER, INTENT(IN) :: IH      ! hour
      INTEGER, INTENT(IN) :: IFAC    ! facility type (1 - 5)
      REAL,    INTENT(IN) :: XEF     ! mg/mile emission factor
      LOGICAL, INTENT(IN) :: OUTPUT  ! true: output data to file or memory
      
      REAL, ALLOCATABLE, SAVE :: HAPAGGR (:,:,:,:,:)  ! aggregated toxics emission factors
                                                      ! dimensions: (IVEH,IP,EFIDX,IH,IFAC)
      REAL, ALLOCATABLE, SAVE :: ZMILES (:,:,:,:)     ! summed miles
                                                      ! dimensions: (IP,EFIDX,IH,IFAC)
      REAL, ALLOCATABLE, SAVE :: ZGMHOUR(:,:,:,:)     ! summed gm/hour
                                                      ! dimensions: (IP,EFIDX,IH,IFAC)
      
      INTEGER :: IVEH        ! SMOKE vehicle type
      INTEGER :: EFTYPE      ! emission factor type (1 - 8)
      INTEGER :: IPNEW       ! pollutant number adjusted for PM species
      INTEGER :: IPOL        ! pollutant number from ADD HAPS file (50 - 99)
      INTEGER :: IOS         ! I/O status
      INTEGER :: IPTMP       ! pollutant loop variable
      INTEGER :: IHTMP       ! hour loop variable
      INTEGER :: IFTMP       ! facility loop variable
      INTEGER :: IEF         ! emission factor type for output
      INTEGER :: SFAC, EFAC  ! start and end facility types
      INTEGER, SAVE :: PREVVEH     ! previous vehicle type
      INTEGER :: NPMHAPS     ! number of PM HAPS

      REAL :: MILES          ! hourly miles
      REAL :: SUMVMT         ! vmt fraction for one of 8 vehicle types
      REAL :: ZGMMILE        ! summed gm/mile
      
      LOGICAL, SAVE :: INITIAL = .TRUE.  ! true: first time through subroutine

      CHARACTER(16) :: PROGNAME = 'SMKHAPS'   ! program name

C  Begin body of subroutine SMKHAPS

C  Allocate and initialize arrays first time
      IF (INITIAL) THEN
          ALLOCATE( HAPAGGR( NVTYPE,NTOTHAPS,6,MAXIH,5 ), STAT=IOS )
          CALL CHECKMEM( IOS, 'HAPAGGR', PROGNAME )

          ALLOCATE( HAPEFS ( NUMSCEN,NVTYPE,NTOTHAPS,8,MAXIH,5 ), 
     &                       STAT=IOS )
          CALL CHECKMEM( IOS, 'HAPEFS', PROGNAME )

          ALLOCATE(  ZMILES( NTOTHAPS,6,MAXIH,5 ), STAT=IOS )
          CALL CHECKMEM( IOS, 'ZMILES', PROGNAME )

          ALLOCATE( ZGMHOUR( NTOTHAPS,6,MAXIH,5 ), STAT=IOS )
          CALL CHECKMEM( IOS, 'ZGMHOUR', PROGNAME )
 
          HAPAGGR = 0.0
          HAPEFS  = 0.0
          ZMILES  = 0.0
          ZGMHOUR = 0.0

          PREVVEH = IV
          INITIAL = .FALSE.
      END IF

      IPNEW = IP

C  If this is the second time through, only bother with PM ratio pollutants
      IF( PSZCUT == 10 ) THEN
          IF( .NOT. OUTPUT ) THEN
              IF( PMHAPS( IP ) == 0 ) RETURN
              IPNEW = PMHAPS( IP )
          END IF
      END IF

C  If the vehicle type has changed, need to aggregate previous values
      IF( PREVVEH /= IV ) THEN

C  Map MOBILE6 vehicle type to SMOKE vehicle type
          CALL SMKCALCVMT( PREVVEH, IVEH, SUMVMT ) 
 
          DO IPTMP = 1,NTOTHAPS
              DO IEF = 1,6
                  
C  Set start and end facility based on emission factor type
                  SFAC = 5
                  EFAC = 5
                      
                  IF( IEF == 1 .OR. IEF == 4 ) THEN
                      SFAC = 1
                  END IF
                      
                  IF( IEF == 4 ) THEN
                      EFAC = 4
                  END IF
                  
                  DO IHTMP = 1,MAXIH
                  
                      DO IFTMP = SFAC,EFAC
                          IF (ZMILES(IPTMP,IEF,IHTMP,IFTMP) < 
     &                                            0.000001) THEN
                              ZGMMILE = 0.0
                          ELSE
                              ZGMMILE =
     &                            ZGMHOUR(IPTMP,IEF,IHTMP,IFTMP)/
     &                            ZMILES( IPTMP,IEF,IHTMP,IFTMP)
                          END IF
                        
                          HAPAGGR(IVEH,IPTMP,IEF,IHTMP,IFTMP) =
     &                        HAPAGGR(IVEH,IPTMP,IEF,IHTMP,IFTMP) +
     &                        (ZGMMILE*VMTMIX(PREVVEH))/SUMVMT
     
                          ZMILES (IPTMP,IEF,IHTMP,IFTMP) = 0.0
                          ZGMHOUR(IPTMP,IEF,IHTMP,IFTMP) = 0.0
                      END DO
                  END DO
              END DO
          END DO
      END IF
          
      PREVVEH = IV

      IF (.NOT. OUTPUT) THEN
                          
C  Calculate and store miles and gm/hour
          MILES = BMYMPD(JDX,IV)*HVMT(IH)*MEVMYR(JDX,IV)
          ZMILES (IPNEW,EFIDX,IH,IFAC) = 
     &                          ZMILES (IPNEW,EFIDX,IH,IFAC) + MILES
          ZGMHOUR(IPNEW,EFIDX,IH,IFAC) = 
     &                          ZGMHOUR(IPNEW,EFIDX,IH,IFAC) + XEF*MILES

      ELSE
   
          DO IPTMP = 1,NMAP
              IPNEW = IPTMP
              
              DO IVEH = 1,NVTYPE
                  DO EFTYPE = 1,8
                  
C  Skip evap crankcase since there are toxics are not available
                      IF( EFTYPE == 7 ) CYCLE
                      
C  Map emission factor type to toxics EF types
C    EFTYPE = 1 -> EFIDX = 1  ! ex. running
C    EFTYPE = 2 -> EFIDX = 1  ! ex. start
C    EFTYPE = 3 -> EFIDX = 2  ! evp. hot soak
C    EFTYPE = 4 -> EFIDX = 3  ! evp. diurnal
C    EFTYPE = 5 -> EFIDX = 5  ! evp. resting
C    EFTYPE = 6 -> EFIDX = 4  ! evp. running
C    EFTYPE = 7 -> N/A        ! evp. crankcase
C    EFTYPE = 8 -> EFIDX = 6  ! evp. refueling

                      SELECT CASE( EFTYPE )
                      CASE( 1 )
                           IEF = 1
                      CASE( 2 )
                           IEF = 1
                      CASE( 3 )
                           IEF = 2
                      CASE( 4 )
                           IEF = 3
                      CASE( 5 )
                           IEF = 5
                      CASE( 6 )
                           IEF = 4
                      CASE( 8 )
                           IEF = 6
                      END SELECT

C  If this is the second time through, only bother with PM ratio pollutants
                      IF( PSZCUT == 10 ) THEN
                          IF( IEF /= 1 ) CYCLE
                          IF( PMHAPS( IPTMP ) == 0 ) CYCLE
                          IPNEW = PMHAPS( IPTMP )
                      END IF
          
C  Set facility types based on EF type
                      IF( EFTYPE == 1 .OR. EFTYPE == 6 ) THEN
                          SFAC = 1
                          EFAC = 4
                      ELSE
                          SFAC = 5
                          EFAC = 5
                      END IF
                      
                      DO IFTMP = SFAC,EFAC

#ifdef M6LIB
                         DO IHTMP = 1,24
                           IF( HAPEFS(STR2INT(DBSCTITLE),IVEH,IPNEW,
     &                                EFTYPE,IHTMP,IFTMP) == 0. ) THEN

C  Divide value by 1000 to account for milligram to gram conversion
                               HAPEFS(STR2INT(DBSCTITLE),IVEH,IPNEW,
     &                                EFTYPE,IHTMP,IFTMP) =
     &                             HAPAGGR(IVEH,IPNEW,IEF,IHTMP,IFTMP)
     &                             / 1000.
                           END IF

C  Subtract PM2.5 values from PM10 emission factors                        
                           IF( IPNEW > NMAP ) THEN
                               HAPEFS(STR2INT(DBSCTITLE),IVEH,IPNEW,
     &                                EFTYPE,IHTMP,IFTMP) =
     &                         HAPEFS(STR2INT(DBSCTITLE),IVEH,IPNEW,
     &                                EFTYPE,IHTMP,IFTMP) -
     &                         HAPEFS(STR2INT(DBSCTITLE),IVEH,IPTMP,
     &                                EFTYPE,IHTMP,IFTMP)
                           END IF
                        
                         END DO
#else                    
                         IPOL = IPNEW + 50  
                         WRITE(DBUNIT(1),110) (DBSCTITLE,TABCHAR,IPOL,
     &                   TABCHAR,IVEH,TABCHAR,EFTYPE,TABCHAR,IFTMP,
     &                   TABCHAR,IHTMP,TABCHAR,
     &                   HAPAGGR(IVEH,IPNEW,IEF,IHTMP,IFTMP)/1000,
     &                   IHTMP=1,24)
  110 FORMAT(A10,A1,I2,A1,I2,A1,I2,A1,I1,A1,I2,A1,E12.5)
#endif
                      END DO  ! facility loop
                  END DO  ! ef type loop
              END DO  ! vehicle loop 
          END DO  ! pollutant loop
          
          PREVVEH = 1
          
          HAPAGGR = 0.0
          ZMILES  = 0.0
          ZGMHOUR = 0.0
      END IF  ! check for output

      END SUBROUTINE SMKHAPS
