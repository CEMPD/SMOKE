
        LOGICAL FUNCTION SETSCCTYPE( TSCC )

C***********************************************************************
C  function body starts at line 68
C
C  DESCRIPTION:
C       Checks SCC code and resets parameters based on type
C
C  PRECONDITIONS REQUIRED:
C       CATEGORY type must be set in MODINFO
C       SCC must be 10-digits long and right-justified
C       8-digit SCCs must start with '00'
C
C  SUBROUTINES AND FUNCTIONS CALLED: none
C
C  REVISION  HISTORY:
C     7/03: Created by C. Seppanen
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

C.........  MODULES for public variables
C.........  This module contains the information about the source category
        USE MODINFO, ONLY: LSCCEND, RSCCBEG, SCCLEV1, SCCLEV2,
     &                     SCCLEV3, CATEGORY
            
        IMPLICIT NONE

C...........   INCLUDES:
        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
       
C........  Function arguments
        CHARACTER(*), INTENT (IN) :: TSCC   ! SCC code

C........  Local variables and their descriptions:
        
        CHARACTER(16) :: PROGNAME = 'SETSCCTYPE' ! program name
        
C***********************************************************************
C   begin body of function SETSCCTYPE

        SETSCCTYPE = .FALSE.

C.........  Don't change any parameters if category is mobile
        IF( CATEGORY == 'MOBILE' ) RETURN

C.........  Check if first two digits of SCC are zero
        IF( TSCC( SCCEXPLEN3+1:SCCEXPLEN3+2 ) == '00' ) THEN

C.............  Only set new values if needed and set flag
            IF( LSCCEND /= SCCEXPLEN3 + 5 ) THEN
                SETSCCTYPE = .TRUE.  ! flag indicates that values have been changed
                LSCCEND = SCCEXPLEN3 + 5
                RSCCBEG = SCCEXPLEN3 + 6
                SCCLEV1 = SCCEXPLEN3 + 3
                SCCLEV2 = SCCEXPLEN3 + 5
                SCCLEV3 = SCCEXPLEN3 + 8
            END IF
        ELSE
            IF( LSCCEND /= SCCEXPLEN3 + 7 ) THEN
                SETSCCTYPE = .TRUE.
                LSCCEND = SCCEXPLEN3 + 7
                RSCCBEG = SCCEXPLEN3 + 8
                SCCLEV1 = SCCEXPLEN3 + 2
                SCCLEV2 = SCCEXPLEN3 + 4
                SCCLEV3 = SCCEXPLEN3 + 7
            END IF
        END IF

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Formatted file I/O formats............ 93xxx

93000   FORMAT( A )  
      
C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )
       
        END FUNCTION SETSCCTYPE
