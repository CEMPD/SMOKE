
        REAL FUNCTION CALCRELHUM( TEMPK, PRESSURE, MIXRATIO )

C***********************************************************************
C  function body starts at line 69
C
C  DESCRIPTION:
C       Calculates relative humidity as a percentage based on air temperature,
C       pressure, and water vapor mixing ratio. Uses Lowe's approximation
C       to calculate saturation vapor pressure, then calculates saturation
C       water vapor mixing ratio.
C
C  PRECONDITIONS REQUIRED: 
C       Temperature in Kelvin
C       Pressure in pascals
C       Water vapor mixing ratio in kg/kg
C
C  SUBROUTINES AND FUNCTIONS CALLED:  none
C
C  REVISION  HISTORY:
C     12/03: Created by C. Seppanen
C
C***********************************************************************
C
C Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
C                System
C File: @(#)$Id$
C
C COPYRIGHT (C) 2004, Environmental Modeling for Policy Development
C All Rights Reserved
C 
C Carolina Environmental Program
C University of North Carolina at Chapel Hill
C 137 E. Franklin St., CB# 6116
C Chapel Hill, NC 27599-6116
C 
C smoke@unc.edu
C
C Pathname: $Source$
C Last updated: $Date$ 
C
C***********************************************************************
        
        IMPLICIT NONE
        
C.........  Function arguments
        REAL, INTENT(IN) :: TEMPK     ! temperature in Kelvin
        REAL, INTENT(IN) :: PRESSURE  ! pressure in pascals
        REAL, INTENT(IN) :: MIXRATIO  ! mixing ratio in kg/kg

C.........  Local parameters
        REAL, PARAMETER :: A0 = 6984.505294         ! constants for Lowe's
        REAL, PARAMETER :: A1 = -188.9039310        ! approximation of
        REAL, PARAMETER :: A2 = 2.133357675         ! saturation vapor
        REAL, PARAMETER :: A3 = -1.288580973E-2     ! pressure
        REAL, PARAMETER :: A4 = 4.393587233E-5
        REAL, PARAMETER :: A5 = -8.023923082E-8
        REAL, PARAMETER :: A6 = 6.136820929E-11

        REAL, PARAMETER :: WVRATIO = 0.622      ! ratio of MW of water vapor to dry air

C.........  Local variables
        REAL    SVP         ! saturation vapor pressure
        REAL    SMR         ! saturation mixing ratio
        REAL    RELHUM      ! relative humidity
        
        CHARACTER(16) :: PROGNAME = 'CALCRELHUM' ! program name
        
C***************************************************************
C   begin body of function CALCRELHUM

C.........  Calculate saturation vapor pressure; uses Lowe's (1977)
C           polynomial approximation to the Clausius Clayperon equation
        SVP = A0 + 
     &          TEMPK * ( A1 +
     &              TEMPK * ( A2 +
     &                  TEMPK * ( A3 +
     &                      TEMPK * ( A4 +
     &                          TEMPK * ( A5 +
     &                              TEMPK * A6 )))))

C.........  Convert saturation vapor pressure from millibars to pascals
        SVP = SVP * 100
     
C.........  Calculate saturation mixing ratio
        SMR = WVRATIO * ( SVP / ( PRESSURE - SVP ))
        
C.........  Calculate relative humidity %
        RELHUM = ( MIXRATIO / SMR ) * 100
        
C.........  Make sure relative humidity is not more than 100%
C           or less than 0%
        RELHUM = MIN( RELHUM, 100.0 )
        RELHUM = MAX( RELHUM, 0.0 )
        
        CALCRELHUM = RELHUM
        
        RETURN
        
        END FUNCTION CALCRELHUM
                
