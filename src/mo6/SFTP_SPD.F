      REAL FUNCTION SFTP_SPD(INH,MY,IV,IP,STD) 
C 
C    Computes the SFTP_SPEED(emissions using a table lookup to obtain 
C    the fraction of vehicles affected and the benefit of the SFTP. 
C    requirement. 
C 
C    Called by: SPEED_CF 
C 
C    Calls: 
C 
C 
C    Changes: (Last change first) 
C 
C    7 Jan 00 ADynTel-ddj 1-020, New Function 
C 
C 
C  Input on call: 
C    argument list: INH, MY, IV, IP, ISTD 
C 
C    common blocks: 
C    /SFTP/   SFTP_FRAC, SFTP_BENE 
C    /VVTYPS/ VVLDGT34 
C 
C 
C  Output on return: 
C   function: 
C             SFTP_SPD 
C 
C 
C  Local variable dictionary: 
C 
C   Name   Type              Description 
C  ------  ----  ------------------------------------------------------- 
C   FRACT   R    Fraction of fleet by model year and vehicles type 
C                subject to the SFTP regulations. 
C   BENEFIT R    Benefit of the supplemental FTP requirements by 
C                pollutant, emitter type, vehicles type and standard. 
C 
      IMPLICIT NONE 
C 
      INCLUDE 'IVTYPE.I' 
      INCLUDE 'SFTP.I' 
C 
      INTEGER, INTENT(IN) :: INH 
      INTEGER, INTENT(IN) :: IP 
      INTEGER, INTENT(IN) :: STD 
      INTEGER, INTENT(IN) :: IV 
      INTEGER, INTENT(IN) :: MY 
C 
      INTEGER IMY 
      INTEGER ISTD 
      INTEGER IVI 
C 
      REAL BENEFIT 
      REAL FRACT 
C 
      IMY = MY - 2000 
C 
      IF(STD .GE. 3 .AND. STD .LE. 8 .AND. IMY .GT. 0) THEN 
          IF(IMY .GT. 4) THEN      ! If Model Year > 2004 
             FRACT = 1.0           ! Full Phase-In. 
          ELSE 
             IF(VVLDGT34(IV) .EQ. 0) THEN 
                   IVI = 1 
                ELSE               ! LDGT3 or LDGT4 
                   IVI = 2 
             END IF 
             FRACT = SFTP_FRAC(IVI,IMY) 
          END IF 
          IF(STD .EQ. 3) THEN      ! Tier1 Standard 
             ISTD = 1 
          ELSE                     ! LEV or Tier2 Standard 
             ISTD = 2 
          END IF 
          BENEFIT = SFTP_BENE(INH,IP,IV,ISTD) 
      ELSE 
          FRACT = 0.0 
          BENEFIT = 0.0 
      END IF 
C 
      SFTP_SPD = 1.0 - FRACT*BENEFIT 
C 
      RETURN 
      END FUNCTION 
