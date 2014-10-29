
        LOGICAL FUNCTION CHKEXPSIC( CSIC )

C***********************************************************************
C  subroutine body starts at line 
C
C  DESCRIPTION:
C      Returns true if SIC is expanded.
C
C  REVISION  HISTORY:
C     Created 10/13 by C. Seppanen
C
C****************************************************************************/
C
C Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
C                System
C File: @(#)$Id$
C
C COPYRIGHT (C) 2005, Environmental Modeling for Policy Development
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
C***************************************************************************

        IMPLICIT NONE

C...........   INCLUDES:
        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters

C...........   EXTERNAL FUNCTIONS and their descriptions:
        LOGICAL, EXTERNAL :: CHKINT

C........  Function arguments
        CHARACTER(SICLEN3), INTENT(IN) :: CSIC   ! SIC code

C***********************************************************************
C   begin body of function CHKEXPSIC

        CHKEXPSIC = .FALSE.

C.........  Check if any of the expanded characters are not zero        
        IF( CSIC( 1:SICEXPLEN3 ) /= REPEAT( '0', SICEXPLEN3 ) ) THEN
            CHKEXPSIC = .TRUE.
            RETURN
        END IF

C.........  Check if last 4 characters are not digits
        IF( .NOT. CHKINT( CSIC( SICLEN3-3:SICLEN3 ) ) ) THEN
            CHKEXPSIC = .TRUE.
            RETURN
        END IF
        
        RETURN
        
        END FUNCTION CHKEXPSIC
