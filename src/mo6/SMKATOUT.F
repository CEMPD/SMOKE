      SUBROUTINE SMKATOUT
       
C  SMKATOUT aggregates air toxics emissions by vehicle age and 
C  vehicle category (8 types).  It then writes the emissions to the
C  database file.

      USE DATABASE, ONLY : DBSELTOX,DBSELEFT,DBSELFAC, 
     &                     DBUNIT,TABCHAR,DBSCTITLE 
      USE ATHEAP, ONLY : DBAT
      USE MODMBSET, ONLY : EMISSIONS
     
      IMPLICIT NONE

      INCLUDE 'ATOX1.I'  ! MAXIAT     
      INCLUDE 'CEVBMY.I' ! BMYMPD 
      INCLUDE 'IVTYPE.I' ! VVDSL,VVHEAVY 
      INCLUDE 'MAXIMA.I' ! MAXYRS 
      INCLUDE 'MYRCAL.I' ! MEVMYR
      INCLUDE 'PART1.I'  ! MAXIPPM
      INCLUDE 'SPEED9.I' ! HVMT
      INCLUDE 'VMXCOM.I' ! VMTMIX,VMTGT12,VMTGT34,VMTHDG,VMTLDDT,VMTHDD

      INTEGER STR2INT
      INTEGER EFPOSITION
      
      EXTERNAL STR2INT, EFPOSITION

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
      INTEGER  ::  IVEH     ! vehicle class (1 - 8)
      INTEGER  ::  EMISPOS  ! position in master emission factor array
      
      REAL              ::  SUMVMT   ! vmt fraction for one of 8 vehicle types 
      REAL              ::  MILES    ! hourly miles                            
      REAL              ::  ZMILES   ! used to sum hourly miles                
      REAL              ::  ZGMHOUR  ! used to sum gm/hour                     
      REAL              ::  ZGMMILE  ! used to calculate hourly gm/mile        
      REAL,DIMENSION(8) ::  ATAGGR   ! fully aggregated AT output              

C  Begin body of subroutine SMKATOUT

      DO IAT = 1,MAXIAT

C  Check that we need this pollutant
         IF (DBSELTOX(IAT) == 1) CYCLE
         IP = IAT + MAXIPPM + 4  ! AT pollutants start after classic and PM pollutants 

         DO ETYPE = 1,8

C  Check that we need this pollutant
            IF (DBSELEFT(ETYPE) == 1) CYCLE

C  Only output evap emissions for Benzene and MTBE
            IF (IAT >= 3 .AND. ETYPE >= 3) EXIT 

C  Map emission factor type to air toxic EF type
            IF (ETYPE <= 2) IXV=1  ! Exhaust Running and Starts 
            IF (ETYPE == 3) IXV=2  ! Evap Hot Soak 
            IF (ETYPE == 4) IXV=3  ! Evap Diurnal 
            IF (ETYPE == 5) IXV=5  ! Evap Resting 
            IF (ETYPE == 6) IXV=4  ! Evap Running 
            IF (ETYPE == 7) CYCLE  ! Evap Crankcase (Not Available) 
            IF (ETYPE == 8) IXV=6  ! Evap Refueling 

C  Set starting and ending facility types based on EF type -
C  Exhaust Running and Evap Running Loss are road dependent (Roads 1-4),
C  Starts and the remaining Evaps are road independent (Road 5)
            IF (ETYPE == 1 .OR. ETYPE == 6) THEN
               STFAC = 1
               ENDFAC = 4
            ELSE
               STFAC = 5
               ENDFAC = 5
            END IF

            DO IFAC = STFAC,ENDFAC

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
                                 
C  Select appropriate vehicle class (8 types) for vehicle type (28)               
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
ccs  Print some error message
                     END SELECT
      
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

C                  WRITE(DBUNIT(1),110) (DBSCTITLE,TABCHAR,IP, 
C     &            TABCHAR,IVEH,TABCHAR,ETYPE,TABCHAR,IFAC,
C     &            TABCHAR,IH,TABCHAR,ATAGGR(IVEH),IVEH=1,8)
C  110 FORMAT(A10,A1,I2,A1,I2,A1,I2,A1,I1,A1,I2,A1,F9.4)

                  DO IVEH = 1,8
                     EMISPOS = EFPOSITION(STR2INT(DBSCTITLE),
     &                                    IP, IVEH, ETYPE, IFAC, 
     &                                    IH, .TRUE. )
                     EMISSIONS(EMISPOS) = ATAGGR(IVEH)
                  END DO  ! loop over 8 vehicle types
               
               END DO  ! loop over hours
            END DO  ! loop over facilities
         END DO  ! loop over emission factor types
      END DO  ! loop over pollutants

      END SUBROUTINE SMKATOUT 
      