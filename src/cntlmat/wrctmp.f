
        SUBROUTINE WRCTMP( IDEV, POLID, IDX, VIDX, LOUTANY )

C***********************************************************************
C  subroutine body starts at line
C
C  DESCRIPTION:
C      This subroutine writes control packet data table indices
C      to a temporary file.
C
C      
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C
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
C***************************************************************************

C.........  MODULES for public variables
C.........  This module contains the information about the source category
        USE MODINFO, ONLY: NSRC, NIPPA

        IMPLICIT NONE

C...........   SUBROUTINE ARGUMENTS:

        INTEGER     , INTENT (IN) :: IDEV           ! logical file name
        INTEGER     , INTENT (IN) :: POLID          ! pollutant number
        INTEGER     , INTENT (IN) :: IDX ( NSRC )   ! index to data tables
        INTEGER     , INTENT (IN) :: VIDX( NIPPA )  ! pollutant/act flags
        LOGICAL     , INTENT(OUT) :: LOUTANY        ! true: at least one pollutant output

C...........   Other local variables

        INTEGER   S   ! indices

	LOGICAL, SAVE :: FIRSTTIME = .TRUE.

C***********************************************************************
C   Begin body of subroutine WRCTMP

        IF ( FIRSTTIME ) THEN
            LOUTANY =  .FALSE.
            FIRSTTIME = .FALSE.
        ENDIF

C.............. Write indices to control factor packets to a temporary file
C               for only those pollutants that have controls
        IF ( VIDX( POLID ) .EQ. 1 ) THEN
            DO S = 1, NSRC

                WRITE( IDEV, '(I8)' ) IDX( S )

            END DO   ! end source loop
            LOUTANY = .TRUE.

        END IF
    
        RETURN

C******************  FORMAT  STATEMENTS   ******************************

        END SUBROUTINE WRCTMP
