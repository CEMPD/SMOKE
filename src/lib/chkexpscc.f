
        LOGICAL FUNCTION CHKEXPSCC( TSCC )

C***********************************************************************
C  subroutine body starts at line 
C
C  DESCRIPTION:
C      Returns true if SCC is expanded.
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

C........  Function arguments
        CHARACTER(SCCLEN3), INTENT(IN) :: TSCC   ! SCC code

C***********************************************************************
C   begin body of function CHKEXPSCC

        CHKEXPSCC = .FALSE.
        
        IF( TSCC( 1:SCCEXPLEN3 ) /= '0000000000' ) THEN
            CHKEXPSCC = .TRUE.
        END IF
        
        RETURN
        
        END FUNCTION CHKEXPSCC
