      SUBROUTINE SMKPMOUT
       
C  SMKPMOUT aggregates particulate emissions by vehicle age and 
C  vehicle category (8 types).  It then writes the emissions to the
C  database file.

      USE DATABASE, ONLY : DBSELPART,DBSELEFT,DBSELFAC, 
     &                     DBUNIT,TABCHAR,DBSCTITLE 
      USE MODMBSET, ONLY : EMISSIONS
           
      IMPLICIT NONE
        
      INCLUDE 'CEVBMY.I' ! BMYMPD 
      INCLUDE 'IVTYPE.I' ! VVDSL,VVGASMC
      INCLUDE 'MAXIMA.I' ! MAXYRS 
      INCLUDE 'MYRCAL.I' ! MEVMYR
      INCLUDE 'PART1.I'  ! MAXIPPM 
      INCLUDE 'PART3.I'  ! DBPM 
      INCLUDE 'SPEED9.I' ! HVMT
      INCLUDE 'VMXCOM.I' ! VMTMIX,VMTGT12,VMTGT34,VMTHDG,VMTLDDT,VMTHDD

      INTEGER STR2INT
      INTEGER EFPOSITION
      
      EXTERNAL STR2INT, EFPOSITION

      INTEGER  ::  EFTYPE   ! emission factor type 
      INTEGER  ::  IFAC     ! facility type
      INTEGER  ::  IROAD    ! road type (only used for output)
      INTEGER  ::  IH       ! hour of day
      INTEGER  ::  IPPM     ! particulate pollutant number (1 - MAXIPPM)
      INTEGER  ::  IP       ! actual pollutant number (7 - 15)
      INTEGER  ::  JDX      ! age
      INTEGER  ::  IV       ! vehicle type (1 - 28)
      INTEGER  ::  IVEH     ! vehicle class (1 - 8)
      INTEGER  ::  EMISPOS  ! position in master emission factor array
      
      REAL              ::  SUMVMT  ! vmt fraction for one of 8 vehicle types
      REAL              ::  MILES   ! hourly miles
      REAL              ::  ZMILES  ! used to sum hourly miles
      REAL              ::  ZGMHOUR ! used to sum gm/hour
      REAL              ::  ZGMMILE ! used to calculate hourly gm/mile
      REAL,DIMENSION(8) ::  PMAGGR  ! fully aggregated PM output

C  Begin body of subroutine SMKPMOUT

      DO IPPM = 1,MAXIPPM-2

C  Check that we need this pollutant
         IF (DBSELPART(IPPM) == 1) CYCLE
         IP = IPPM+6  ! PM pollutants start after classic pollutants

C  Set emission factor type based on pollutant 
C  (always exhaust running except for brake and tire pollutants)
         IF (IPPM <= 7) THEN
            EFTYPE = 1
         ELSE
            EFTYPE = IPPM + 1
         END IF

         DO IFAC = 1,4

C  Only need to do brake and tire pollutants once
            IF (IFAC > 1 .AND. IPPM > 7) CYCLE
                     
            DO IH = 1,24
C  Reinitialize output array
               PMAGGR = 0.
            
               DO IV = 1,MAXVEH

C  Skip unused vehicle and emission factor combinations

C.....            Gas vehicles do not have OCARBON or ECARBON output
                  IF ((IPPM == 2 .OR. IPPM == 3) .AND. 
     &                              VVGASMC(IV) == 1) CYCLE

C.....            No GASPM output for diesel vehicles
                  IF (IPPM == 4 .AND. VVDSL(IV) == 1) CYCLE
               
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
                     MILES = BMYMPD(JDX,IV)*HVMT(IH)*MEVMYR(JDX,IV)
                     ZMILES = ZMILES + MILES
                     ZGMHOUR = ZGMHOUR + DBPM(IPPM,IFAC,IH,JDX,IV)*MILES
                  END DO  ! loop over ages
               
                  IF (ZMILES < 0.000001 ) THEN
                     ZGMMILE = 0.0
                  ELSE
                     ZGMMILE = ZGMHOUR/ZMILES
                  END IF
 
                  PMAGGR(IVEH) = PMAGGR(IVEH) + 
     &                           (ZGMMILE*VMTMIX(IV))/SUMVMT
 
               END DO  ! loop over 28 vehicle types

C.....  Check that this EF type, and road type should be output
               IF (DBSELEFT(EFTYPE) == 2 ) THEN
                     	
C.....  Set road type for output
                  IF (IPPM <= 7) THEN
                     IROAD = IFAC
                  ELSE
                     IROAD = 5
                  END IF
                     
                  IF (DBSELFAC(IROAD) == 2 ) THEN

C                     WRITE(DBUNIT(1),110) (DBSCTITLE,TABCHAR,IP, 
C     &               TABCHAR,IVEH,TABCHAR,EFTYPE,TABCHAR,IROAD,
C     &               TABCHAR,IH,TABCHAR,PMAGGR(IVEH),IVEH=1,8)
C  110 FORMAT(A10,A1,I2,A1,I2,A1,I2,A1,I1,A1,I2,A1,F9.4)

                     DO IVEH = 1,8
                        EMISPOS = EFPOSITION(STR2INT(DBSCTITLE),
     &                               IP, IVEH, EFTYPE, IROAD, IH)
                        EMISSIONS(EMISPOS) = PMAGGR(IVEH)
                     END DO  ! loop over 8 vehicle types

                  END IF  ! road check
               END IF  ! emission factor check
               
            END DO  ! loop over hours
         END DO  ! loop over facilities
      END DO  ! loop over pollutants

      END SUBROUTINE SMKPMOUT 
      