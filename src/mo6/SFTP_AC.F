      REAL FUNCTION SFTP_AC(INH,MY,IV,IP,STD,SPEED,FRACT) 
C 
C    Computes the off-cycle emission rate for light duty gasoline 
C    vehicles and trucks.  The off-cycle emission rate is added to 
C    the base emissions prior to applying the speed correction. 
C 
C    Called by: ACCF 
C 
C    Calls:  MIN, MAX, STANDARD 
C 
C 
C    Changes: (Last change first) 
C 
C  20 Jul 00 @EPA-ddj Bug194 Removed NOF from CO SFTP calculation. 
C  16 May 00 @EPA- BG BUG#88 changed NOF from INTEGER to REAL 
C   02 May 00 @EPA-ddj E-7, Changed DATA statement for P to a 
C             type statement. 
C    7 Jan 00 @DynTel-ddj 1-020, New Function 
C 
C 
C  Input on call: 
C    argument list: INH, MY, IV, IP, STD, SPEED, FRACT 
C 
C    common blocks: 
C    /SFTP/   SFTP_FRAC, SFTP_BENE 
C    /IVTYPS/ VVLDGT34 
C    /VTNAMS/ VTLDGV, VTLDGT1 
C 
C 
C  Output on return: 
C    function:  SFTP_AC 
C 
C    variable:  FRACT 
C 
C 
C  Local variable dictionary: 
C 
C   Name   Type              Description 
C  ------  ----  ------------------------------------------------------- 
C   FRACT   R    Fraction of fleet by model year and vehicles type 
C                subject to the SFTP regulations. 
C   BENEFIT R    AC benefit of the supplemental FTP requirements by 
C                pollutant, emitter type, vehicles type and standard. 
C   P(3,2)  R    Polynomial constants for CO Multiplicative Effect. 
C 
C 
      IMPLICIT NONE 
C 
      INCLUDE 'AC.I' 
      INCLUDE 'IVTYPE.I' 
      INCLUDE 'SFTP.I' 
C 
      INTEGER, INTENT(IN) :: INH 
      INTEGER, INTENT(IN) :: IP 
      INTEGER, INTENT(IN) :: STD 
      INTEGER, INTENT(IN) :: IV 
      INTEGER, INTENT(IN) :: MY 
C 
      REAL, INTENT(IN)    :: SPEED 
      REAL, INTENT(OUT)   :: FRACT 
C 
      INTEGER CO 
      INTEGER IMY 
      INTEGER IPI 
      INTEGER ISTD 
      INTEGER STANDARD 
C 
      REAL BENEFIT 
C 
      PARAMETER (CO = 2) 
C 
      REAL, DIMENSION(3,2)  :: P = RESHAPE 
     & ( (/1.34,  -0.006134,  0.000053, 1.27,  -0.004939,  0.000048/), 
     &   (/3,2/) ) 
C 
C     initialize function value 
C 
      FRACT = 0.0 
      SFTP_AC = 1.0 
C 
      IMY = MY - 2000 
      ISTD = STANDARD(STD) - 1 
C 
      IF(MY .GT. 2000 ) THEN 
         IF(IMY .GT. 4) THEN 
            FRACT = 1.0 
         ELSE IF(VVLDGT34(IV) .EQ. 0) THEN 
               FRACT = SFTP_FRAC(1,IMY) 
            ELSE                     ! LDGT3 or LDGT4 
               FRACT = SFTP_FRAC(2,IMY) 
         END IF 
         IF(IP .NE. CO) THEN 
            IF(ISTD .LE. 0 .OR. ISTD .EQ. 3) THEN 
               BENEFIT = 0.0         ! Tier0 or ULEV 
            ELSE 
               IF(ISTD .GE. 3) THEN  ! ISTD index is 1 for Tier1, 
                  ISTD = 2           ! 2 for LEV and Tier2. 
               END IF 
               IPI = MIN(2,IP)       ! IPI index is 1 for HC, 2 for NOx 
               BENEFIT = SFTP_BENE_AC(INH,IPI,IV,ISTD) 
            END IF 
            SFTP_AC = 1.0 - FRACT*BENEFIT 
C 
         ELSE                        ! For pollutant CO. 
C 
C     Light Duty Gasoline vehicle/truck and CO (multiplicative effect) 
C 
            IF(IV .EQ. VTLDGV .OR. IV .EQ. VTLDGT1) THEN 
               SFTP_AC = (P(1,1) + SPEED*(P(2,1) + SPEED*P(3,1)) ) 
            ELSE                     !  LDGT2, LDGT3 or LDGT4 
               SFTP_AC = (P(1,2) + SPEED*(P(2,2) + SPEED*P(3,2)) ) 
            END IF 
         END IF 
C 
         SFTP_AC = MAX(SFTP_AC, 0.0) 
      END IF 
C 
      RETURN 
      END FUNCTION                   ! SFTP_AC 
