      SUBROUTINE SMKEXOUT(HEF,IP,IV,JDX)
      
C  SMKEXOUT aggregates exhaust running and start emissions by vehicle age and
C  vehicle category (8 types).  It then writes the emmisions to the database
C  file.

      USE SMKSAVE,  ONLY : EXHAUST, ZMILES_EX, ZGMMILE_EX,
     &                     ZGMHOUR_EX
      USE DATABASE, ONLY : DBSELPOL, DBSELVEH, DBSELEFT, 
     &                     DBSCTITLE, TABCHAR, DBUNIT,
     &                     DBSELFAC
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
      
      REAL,DIMENSION(24,5),INTENT(IN)  ::  HEF  ! gm/mile emission factors
      INTEGER,INTENT(IN)               ::  IP   ! pollutant
      INTEGER,INTENT(IN)               ::  IV   ! vehicle type (1 - 28)
      INTEGER,INTENT(IN)               ::  JDX  ! age
      
      INTEGER :: IH        ! hour
      INTEGER :: IFAC      ! facility type
      INTEGER :: IVEH      ! vehicle type (1 - 8)
      INTEGER :: EFTYPE    ! emission factor type
      INTEGER :: EMISPOS   ! position in master emission factor array
      
      REAL :: MILES        ! hourly miles
      REAL :: SUMVMT       ! vmt fraction for one of 8 vehicle types
      
      LOGICAL, SAVE :: INITIAL = .TRUE.  ! true: first time through subroutine
      
C  First time, initialize emission factor array
      IF (INITIAL) THEN
         EXHAUST    = 0.
         ZMILES_EX  = 0.
         ZGMHOUR_EX = 0.
         ZGMMILE_EX = 0.

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
ccs  Print some error message
      END SELECT

      DO IH = 1,24
         MILES = BMYMPD(JDX,IV)*HVMT(IH)*MEVMYR(JDX,IV)
         ZMILES_EX(IH) = ZMILES_EX(IH) + MILES
      
         DO IFAC = 1,5
            ZGMHOUR_EX(IH,IFAC) = ZGMHOUR_EX(IH,IFAC)+HEF(IH,IFAC)*MILES
  
C.....      If age = 0, we've finished summing and can get the aggregated
C           gm/mile over ages, then aggregate by vehicle type
            IF (JDX-1 == 0) THEN
               IF (ZMILES_EX(IH) < 0.000001) THEN
                  ZGMMILE_EX(IH,IFAC) = 0.0
               ELSE
                  ZGMMILE_EX(IH,IFAC) = 
     &                            ZGMHOUR_EX(IH,IFAC)/ZMILES_EX(IH)
               END IF
            
               EXHAUST(IH,IFAC,IVEH) = EXHAUST(IH,IFAC,IVEH) + 
     &            (ZGMMILE_EX(IH,IFAC)*VMTMIX(IV))/SUMVMT    
     
C.....         Reinitialize arrays after last hour and facility
               IF (IH == 24 .AND. IFAC == 5) THEN
                  ZMILES_EX  = 0.
                  ZGMHOUR_EX = 0.
                  ZGMMILE_EX = 0.
               END IF
            END IF
         END DO
      END DO

C  If vehicle type = 28 and age = 0, we're done with this pollutant, so output array
      IF (IV == 28 .AND. JDX-1 == 0) THEN
         DO IVEH = 1,8
         IF (DBSELPOL(IP) == 2) THEN
            IF (DBSELEFT(1) == 2) THEN
               EFTYPE = 1
               DO IFAC = 1,4
                  IF (DBSELFAC(IFAC) == 2) THEN
#ifdef M6LIB     
                     DO IH = 1,24
                        EMISPOS = EFPOSITION( STR2INT( DBSCTITLE ),
     &                                       IP, IVEH, EFTYPE, IFAC, 
     &                                       IH, .TRUE. )
                        EMISSIONS( EMISPOS ) = EXHAUST( IH,IFAC,IVEH )
                     END DO
#else
                     WRITE(DBUNIT(1),110) (DBSCTITLE,TABCHAR,IP,TABCHAR,
     &                  IVEH,TABCHAR,EFTYPE,TABCHAR,IFAC,TABCHAR,IH,
     &                  TABCHAR,EXHAUST(IH,IFAC,IVEH),IH=1,24)
#endif
                  END IF
               END DO
            END IF
                             
            IF (DBSELEFT(2) == 2) THEN
               EFTYPE = 2               
               IFAC = 5
               IF (DBSELFAC(IFAC) == 2) THEN
#ifdef M6LIB  
                  DO IH = 1,24
                     EMISPOS = EFPOSITION( STR2INT( DBSCTITLE ), 
     &                                    IP, IVEH, EFTYPE, IFAC, 
     &                                    IH, .TRUE. )
                     EMISSIONS( EMISPOS ) = EXHAUST( IH,IFAC,IVEH )
                  END DO
#else
                  WRITE(DBUNIT(1),110) (DBSCTITLE,TABCHAR,IP,TABCHAR,
     &               IVEH,TABCHAR,EFTYPE,TABCHAR,IFAC,TABCHAR,IH,
     &               TABCHAR,EXHAUST(IH,IFAC,IVEH),IH=1,24)
  110 FORMAT(A10,A1,I1,A1,I2,A1,I1,A1,I1,A1,I2,A1,F9.4)
#endif
               END IF
            END IF
         END IF
         END DO
         
C.....   Reinitialize emission factors array         
         EXHAUST = 0.
      END IF
      
      END SUBROUTINE SMKEXOUT
      
      
