      MODULE SMKSAVE
      
      IMPLICIT NONE
      
C.... Aggregated exhaust running and start emission factors
      REAL, DIMENSION(24,5,8) :: EXHAUST     ! (hours, facility types, vehicle types)

C.... Aggregated evaporative emission factors
      REAL, DIMENSION(24,4,4) :: EVPRUN      ! (hours, facility types, vehicle types)      
      REAL, DIMENSION(24,5,5) :: EVPOTHER    ! (hours, vehicle types, emission types)
      
C.... Temporary arrays for summing values

C     Hourly miles for each emission factor type
      REAL, DIMENSION(24)     :: ZMILES_EX      ! (hours)
      REAL, DIMENSION(24,5)   :: ZMILES_EVPRUN  ! (hours, facility types)
      REAL, DIMENSION(24,5)   :: ZMILES_EVP     ! (hours, emission types)
      
C     Summed gm/hour emission factors for each emission factor type
      REAL, DIMENSION(24,5)   :: ZGMHOUR_EX     ! (hours, facility types)
      REAL, DIMENSION(24,5)   :: ZGMHOUR_EVPRUN ! (hours, facility types)
      REAL, DIMENSION(24,5)   :: ZGMHOUR_EVP    ! (hours, emission types)

C     Summed gm/mile emission factors for each emission factor type
      REAL, DIMENSION(24,5)   :: ZGMMILE_EX     ! (hours, facility types)
      REAL, DIMENSION(24,5)   :: ZGMMILE_EVPRUN ! (hours, facility types)
      REAL, DIMENSION(24,5)   :: ZGMMILE_EVP    ! (hours, emission types)
      
      END MODULE SMKSAVE
      
