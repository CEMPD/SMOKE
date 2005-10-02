      SUBROUTINE SMKCALCVMT( M6TYPE, SMKTYPE, SUMVMT )

C  SMKCALCVMT determines the matching SMOKE vehicle type for a given MOBILE6
C  vehicle type. It also calculates and returns the VMT fraction to use for
C  each SMOKE vehicle type.
      
      USE MODEMFAC, ONLY: M6VEHMAP
      
      IMPLICIT NONE
      
      INCLUDE 'VMXCOM.I'    ! VMTMIX
      
      INCLUDE 'M6CNST3.EXT' ! MXM6VTYP
      
      INTEGER, INTENT(IN)  :: M6TYPE   ! MOBILE6 vehicle type
      INTEGER, INTENT(OUT) :: SMKTYPE  ! matching SMOKE vehicle type
      REAL,    INTENT(OUT) :: SUMVMT   ! VMT fraction for SMOKE vehicle type
      
      LOGICAL, SAVE :: CALCVMT( MXM6VTYP ) ! true: calculate VMT fraction
      REAL,    SAVE :: VMTFRAC( MXM6VTYP ) ! VMT fraction
      
      INTEGER TYPE    ! vehicle type loop variable
      
      LOGICAL, SAVE :: INITIAL = .TRUE.   ! true: first time through subroutine
      
      CHARACTER(16) :: PROGNAME = 'SMKCALCVMT'  ! program name

C  First time, initialize storage arrays      
      IF( INITIAL ) THEN
          CALCVMT = .TRUE.
          VMTFRAC = 0.
        
          INITIAL = .FALSE.
      END IF
      
      SMKTYPE = M6VEHMAP( M6TYPE )

C  Check if VMT fraction needs to be calculated for this type      
      IF( CALCVMT( M6TYPE ) ) THEN
          SUMVMT = 0.
          DO TYPE = 1, MXM6VTYP
              IF( M6VEHMAP( TYPE ) == SMKTYPE ) THEN
                  SUMVMT = SUMVMT + VMTMIX( TYPE )
              END IF
          END DO
        
          VMTFRAC( M6TYPE ) = SUMVMT
          CALCVMT( M6TYPE ) = .FALSE.
      END IF
      
      SUMVMT = VMTFRAC( M6TYPE )
      
      END SUBROUTINE SMKCALCVMT
