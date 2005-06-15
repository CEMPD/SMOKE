
        SUBROUTINE GETSPDFMT( INTFMT, REALFMT )

C***********************************************************************
C  subroutine body starts at line 56
C
C  DESCRIPTION:
C       Creates the format strings needed for reading and writing the 
C       SPDSUM file. Uses the parameter SPDLEN3 to determine how wide
C       the speed / speed profile field should be.
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:  none
C
C  REVISION  HISTORY:
C     06/05: Created by C. Seppanen
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
        
C.........  INCLUDES
        INCLUDE 'EMCNST3.EXT'   ! emissions constant parameters
        
C.........  SUBROUTINE ARGUMENTS
        CHARACTER(*), INTENT(OUT) :: INTFMT     ! format string with integer
        CHARACTER(*), INTENT(OUT) :: REALFMT    ! format string with real
        
C.........  Other local variables
        CHARACTER(50) BUFFER        ! temporary buffer
        
        CHARACTER(16) :: PROGNAME = 'GETSPDFMT'
        
C***********************************************************************
C   begin body of subroutine GETSPDFMT

        INTFMT = '(I6, 1X, I1, 1X,'
        REALFMT = INTFMT
        
        IF( SPDLEN3 >= 10 ) THEN
            WRITE(BUFFER, '(I2)') SPDLEN3
        ELSE
            WRITE(BUFFER, '(I1)') SPDLEN3
        END IF
        
        INTFMT = TRIM( INTFMT ) // ' I' // TRIM( BUFFER )
        REALFMT = TRIM( REALFMT ) // ' F' // TRIM( BUFFER ) // '.2'
        
        BUFFER = ', 7( 1X, I6 ), 1X, A1 )'
        
        INTFMT = TRIM( INTFMT ) // TRIM( BUFFER )
        REALFMT = TRIM( REALFMT ) // TRIM( BUFFER )
        
        END SUBROUTINE GETSPDFMT
