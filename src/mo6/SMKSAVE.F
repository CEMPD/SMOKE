      MODULE SMKSAVE
      
      IMPLICIT NONE
      
      REAL, DIMENSION(24,5,8) :: EXHAUST     ! (hours, facility types, vehicle types)
      
      REAL, DIMENSION(24,4,4) :: EVPRUN      ! (hours, facility types, vehicle types)      
      REAL, DIMENSION(24,5,5) :: EVPOTHER    ! (hours, vehicle types, emission types)
      
C.... Temporary arrays for summing values      
      REAL, DIMENSION(24)     :: ZMILES_EX      ! (hours)
      REAL, DIMENSION(24,5)   :: ZMILES_EVPRUN  ! (hours, facility types)
      REAL, DIMENSION(24,5)   :: ZMILES_EVP     ! (hours, emission types)
      
      REAL, DIMENSION(24,5)   :: ZGMHOUR_EX     ! (hours, facility types)
      REAL, DIMENSION(24,5)   :: ZGMHOUR_EVPRUN ! (hours, facility types)
      REAL, DIMENSION(24,5)   :: ZGMHOUR_EVP    ! (hours, emission types)
      
      REAL, DIMENSION(24,5)   :: ZGMMILE_EX     ! (hours, facility types)
      REAL, DIMENSION(24,5)   :: ZGMMILE_EVPRUN ! (hours, facility types)
      REAL, DIMENSION(24,5)   :: ZGMMILE_EVP    ! (hours, emission types)
      
      END MODULE SMKSAVE
      
