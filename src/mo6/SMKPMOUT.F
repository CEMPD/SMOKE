      SUBROUTINE SMKPMOUT
       
C  SMKPMOUT aggregates particulate emissions by vehicle age and 
C  vehicle category (8 types).  It then writes the emissions to the
C  database file.

      USE DATABASE, ONLY : DBSELPART,DBSELEFT,DBSELFAC, 
     &                     DBUNIT,TABCHAR,DBSCTITLE 
#ifdef M6LIB
      USE MODMBSET, ONLY : EMISSIONS
#endif

      IMPLICIT NONE
        
      INCLUDE 'CEVBMY.I' ! BMYMPD 
      INCLUDE 'IVTYPE.I' ! VVDSL,VVGASMC
      INCLUDE 'MAXIMA.I' ! MAXYRS 
      INCLUDE 'MYRCAL.I' ! MEVMYR
      INCLUDE 'PART1.I'  ! MAXIPPM, PSZCUT 
      INCLUDE 'PART3.I'  ! DBPM 
      INCLUDE 'SPEED9.I' ! HVMT
      INCLUDE 'VMXCOM.I' ! VMTMIX,VMTGT12,VMTGT34,VMTHDG,VMTLDDT,VMTHDD

      INCLUDE 'M6CNST3.EXT'   !  Mobile6 constants
        
      INTEGER STR2INT
      INTEGER EFPOSITION
      
      EXTERNAL STR2INT, EFPOSITION

      INTEGER  ::  I        ! loop variable
      INTEGER  ::  EFTYPE   ! emission factor type 
      INTEGER  ::  IFAC     ! facility type
      INTEGER  ::  IH       ! hour of day
      INTEGER  ::  IPPM     ! particulate pollutant number (1 - MAXIPPM)
      INTEGER  ::  IP       ! actual pollutant number (7 - 15)
      INTEGER  ::  JDX      ! age
      INTEGER  ::  IV       ! vehicle type (1 - 28)
      INTEGER  ::  IVEH     ! vehicle class (1 - 8)
      INTEGER  ::  SMKPOL   ! SMOKE pollutant number based on M6 index
      INTEGER  ::  POL25    ! SMOKE pollutant number for 2.5 pollutant
      
      REAL              ::  SUMVMT  ! vmt fraction for one of 8 vehicle types
      REAL              ::  MILES   ! hourly miles
      REAL              ::  ZMILES  ! used to sum hourly miles
      REAL              ::  ZGMHOUR ! used to sum gm/hour
      REAL              ::  ZGMMILE ! used to calculate hourly gm/mile
      REAL,DIMENSION(8) ::  PMAGGR  ! fully aggregated PM output
      REAL              ::  SUBVAL  ! 2.5 value to subtract from 10 value

      CHARACTER(LEN=300) :: MESG                   ! message buffer
      CHARACTER(LEN=16)  :: PROGNAME = 'SMKPMOUT'  ! program name

C  Begin body of subroutine SMKPMOUT

C  Loop over particulate pollutants
      DO IPPM = 1,MAXIPPM

C  Skip PM idle pollutants
         IF (IPPM == 10 .OR. IPPM == 11) CYCLE

C  Set M6 pollutant index; PM pollutants start after classic pollutants
         IP = IPPM + 6

C  Skip lead since it is not used by SMOKE
         IF (IP == 11) CYCLE

C  If particle size is 10, this is the second time through this routine -
C  skip pollutants that don't depend on size (SO4, SO2, NH3)
         IF (PSZCUT == 10) THEN
             IF (IP == 7 .OR. IP == 12 .OR. IP == 13) CYCLE
             
C  Reset pollutant number to not interfere with 2.5 pollutants
             IP = IP + 20
         END IF

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

C  Set emission factor type based on pollutant 
C  (always exhaust running except for brake and tire pollutants)
         IF (IPPM <= 7) THEN
            EFTYPE = 1
         ELSE
            EFTYPE = IPPM + 1
         END IF

C  Check if this pollutant/process combo is valid
         IF (M6POL2EF( EFTYPE, SMKPOL ) == -1) CYCLE

         DO IFAC = 1,4

C  Check if this facility/process combo is valid
            IF (M6FAC2EF( EFTYPE, IFAC ) == -1) CYCLE
                     
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
                     SUMVMT = VMTMIX(IV)
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
                     SUMVMT = VMTMIX(IV)
                  CASE (15,28)
                     IVEH = 6
                     SUMVMT = VMTLDDT
                  CASE (16:23,26,27)
                     IVEH = 7
                     SUMVMT = VMTHDD
                  CASE (24)
                     IVEH = 8
                     SUMVMT = VMTMIX(IV)
                  CASE DEFAULT
                     MESG = 'INTERNAL ERROR: Unexpected vehicle type'
                     CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
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

#ifdef M6LIB
               DO IVEH = 1,8
C  Check that this vehicle/process combo is valid
                  IF (M6VEH2EF( EFTYPE, IVEH ) == -1) CYCLE
     
C  Subtract 2.5 values from 10 values if needed
                  IF (PSZCUT == 10) THEN

C  Find SMOKE pollutant number for 2.5 pollutant
                     POL25 = 0
                     DO I = 1,MXM6POLS
                         IF( IP - 20 == SMK2M6POL( I ) ) THEN
                             POL25 = I
                             EXIT
                         END IF
                     END DO
                     
                     IF( POL25 == 0 ) THEN
                         MESG = 'INTERNAL ERROR: Unexpected pollutant'
                         CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                     END IF

C  Check that new pollutant is valid for current process
                     IF (M6POL2EF( EFTYPE, POL25 ) == -1) THEN
                         SUBVAL = 0.
                     ELSE
                         SUBVAL = EMISSIONS( EFTYPE )%PTR(
     &                       STR2INT( DBSCTITLE ),
     &                       M6POL2EF( EFTYPE, POL25 ),
     &                       M6VEH2EF( EFTYPE, IVEH ),
     &                       M6FAC2EF( EFTYPE, IFAC ), IH )
                     END IF
                  ELSE
                     SUBVAL = 0.
                  END IF
                   
                  EMISSIONS( EFTYPE )%PTR( 
     &                STR2INT( DBSCTITLE ),
     &                M6POL2EF( EFTYPE, SMKPOL ),
     &                M6VEH2EF( EFTYPE, IVEH ),
     &                M6FAC2EF( EFTYPE, IFAC ), IH ) = 
     &            PMAGGR( IVEH ) - SUBVAL
                        
               END DO  ! loop over 8 vehicle types
#else
               WRITE(DBUNIT(1),110) (DBSCTITLE,TABCHAR,IP, 
     &            TABCHAR,IVEH,TABCHAR,EFTYPE,TABCHAR,IFAC,
     &            TABCHAR,IH,TABCHAR,PMAGGR(IVEH),IVEH=1,8)
  110 FORMAT(A10,A1,I2,A1,I2,A1,I2,A1,I1,A1,I2,A1,E12.5)
#endif
               
            END DO  ! loop over hours
         END DO  ! loop over facilities
      END DO  ! loop over pollutants

      END SUBROUTINE SMKPMOUT 
