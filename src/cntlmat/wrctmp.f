
        SUBROUTINE WRCTMP( IDEV, POLID, IDX, VIDX )

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
C COPYRIGHT (C) 2002, MCNC Environmental Modeling Center
C All Rights Reserved
C
C See file COPYRIGHT for conditions of use.
C
C Environmental Modeling Center
C MCNC
C P.O. Box 12889
C Research Triangle Park, NC  27709-2889
C
C smoke@emc.mcnc.org
C
C Pathname: $Source$
C Last updated: $Date$ 
C
C***************************************************************************

C.........  MODULES for public variables
C.........  This module contains the information about the source category
        USE MODINFO

        IMPLICIT NONE

C...........   SUBROUTINE ARGUMENTS:

        INTEGER     , INTENT (IN) :: IDEV           ! logical file name
        INTEGER     , INTENT (IN) :: POLID          ! pollutant number
        INTEGER     , INTENT (IN) :: IDX ( NSRC )   ! index to data tables
        INTEGER     , INTENT(OUT) :: VIDX( NIPPA )  ! pollutant/act flags

C...........   Other local variables

        INTEGER   S   ! indices

C***********************************************************************
C   Begin body of subroutine WRCTMP

C.............. Write indices to control factor packets to a temporary file
C               for only those pollutants that have controls
        IF ( VIDX( POLID ) .EQ. 1 ) THEN
            DO S = 1, NSRC

                WRITE( IDEV, '(I8)' ) IDX( S )

            END DO   ! end source loop
        END IF

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

        END SUBROUTINE WRCTMP
