      BLOCK DATA BD_SFTP 
C 
C   BD_SFTP contains coefficients used in determining the speed 
C   and A/C correction factors for light Duty Gasoline Vehicles 
C   and Trucks. 
C  
C   Changes. 
C    7 Jan 00 @DynTel-ddj 1-020  New block data. 
C 
C   Common block array subscripts: 
C 
C       SFTP_FRAC     (2, 4)           (IV, IMY) 
C       SFTP_BENE_AC  (2, 2, 5, 2)     (INH, IP, IV, ISTD) 
C       SFTP_BENE     (2, 3, 5, 2)     (INH, IP, IV, ISTD) 
C 
C 
C   Common block dictionary: 
C 
C     Name       Type              Description 
C    ------      ----  --------------------------------------------- 
C /SFTP/ 
C  SFTP_FRAC      R    Fraction of vehicles subject to the SFTP 
C                      Regulations. 
C  SFTP_BENE      R    Fractional reduction in off-cycle emissionf due to 
C                      SFTP regulations by pollutant, vehicle type, 
C                      emission standard and emitter category. 
C  SFTP_BENE_AC   R    Fractional reduction in air conditioning emissions 
C                      due to SFTP regulations by pollutant (HC and NOx), 
C                      vehicle type, emission standard and emitter 
C                      category. 
C 
C 
C  Local array subscripts: 
C 
C 
C 
C 
C  Notes: 
C 
C 
C 
      IMPLICIT NONE 
C 
      INCLUDE 'SFTP.I' 
C 
      DATA SFTP_FRAC/ 
c 
C 200X   (phase-in year, 1=2001, 2=2002, 3=2003, 4=2004 
C 
C      LDGV, LDGT1&2       LDGT3&4 
     1     0.25,            0.00, 
     2     0.50,            0.40, 
     3     0.85,            0.60, 
     4     1.00,            1.00/ 
C 
      DATA SFTP_BENE/ 
C 
C   Continuation code:  V=LDV, 1=LDT1, 2=LDT2, 3=LDT3, 4=LDT4 
C 
C Tier1      HC            CO           NOx 
C       normal high   normal high   normal high 
     V  0.88, 0.44,   0.78, 0.39,   0.78,  0.39, 
     1  0.88, 0.44,   0.78, 0.39,   0.78,  0.39, 
     2  0.88, 0.44,   0.78, 0.39,   0.78,  0.39, 
     3  0.88, 0.44,   0.78, 0.39,   0.78,  0.39, 
     4  0.88, 0.44,   0.78, 0.39,   0.78,  0.39, 
C LEV 
     V  0.93, 0.47,   0.78, 0.39,   0.92,  0.46, 
     1  0.93, 0.47,   0.78, 0.39,   0.92,  0.46, 
     2  0.94, 0.47,   0.78, 0.39,   0.93,  0.47, 
     3  0.90, 0.45,   0.78, 0.39,   0.83,  0.42, 
     4  0.90, 0.45,   0.78, 0.39,   0.81,  0.41/ 
C 
      DATA SFTP_BENE_AC/ 
C Tier1      HC            NOx 
C       normal high   normal high 
     V  1.00, 1.00,   0.50,  0.25, 
     1  1.00, 1.00,   0.50,  0.25, 
     2  1.00, 1.00,   0.50,  0.25, 
     3  1.00, 1.00,   0.50,  0.25, 
     4  1.00, 1.00,   0.50,  0.25, 
C LEV 
     V  1.00, 1.00,   0.70,  0.35, 
     1  1.00, 1.00,   0.70,  0.35, 
     2  1.00, 1.00,   0.85,  0.43, 
     3  1.00, 1.00,   0.79,  0.39, 
     4  1.00, 1.00,   0.78,  0.39/ 
C 
      END 
