
        SUBROUTINE WRCTMP( IDEV, IGRP, NGSZ, IDX, VIDX )

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
C****************************************************************************/
C
C Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
C                System
C File: @(#)$Id$
C
C COPYRIGHT (C) 1999, MCNC--North Carolina Supercomputing Center
C All Rights Reserved
C
C See file COPYRIGHT for conditions of use.
C
C Environmental Programs Group
C MCNC--North Carolina Supercomputing Center
C P.O. Box 12889
C Research Triangle Park, NC  27709-2889
C
C env_progs@mcnc.org
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

        INTEGER     , INTENT (IN) :: IDEV             ! logical file name
        INTEGER     , INTENT (IN) :: IGRP             ! pollutant group
        INTEGER     , INTENT (IN) :: NGSZ             ! # pollutants/group
        INTEGER     , INTENT (IN) :: IDX(NSRC,NGSZ)   ! index to data tables
        INTEGER     , INTENT(OUT) :: VIDX(NIPOL)      ! pollutant flags


C...........   Other local variables

        INTEGER   I,J,K,S   ! indices

C***********************************************************************
C   Begin body of subroutine WRCTMP

        DO I = 1,NGSZ

        K = NGSZ*IGRP - NGSZ + I  ! compute index to master pollutant
                                  ! list for current pollutant

C............. If pollutant has controls applied, write indices to control
C              data tables to a temporary file
        IF ( VIDX(K) .EQ. 1 ) THEN

           DO S = 1,NSRC

              WRITE( IDEV, * ) IDX( S, I )

           END DO   ! end source loop

        END IF

        END DO   ! end pollutant group loop

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

        END SUBROUTINE WRCTMP
