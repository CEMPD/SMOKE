
        SUBROUTINE SRCGRPCNT( NSRC, NMAT1, NX, IX )

C***********************************************************************
C  subroutine body starts at line 
C
C  DESCRIPTION:
C     This subroutine matches source groups to grid cells to count the
C     total number of data points. This subroutine only exists because
C     of the way the gridding matrix is stored in memory.
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C     Created 7/2013 by C. Seppanen
C
C***************************************************************************
C
C Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
C                System
C File: @(#)$Id$
C
C COPYRIGHT (C) 2013, Environmental Modeling for Policy Development
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
C**************************************************************************

C.........  MODULES for public variables
        USE MODMERGE, ONLY: ISRCGRP, NSGOUTPUT, GRPCNT, PFLAG

C.........  This module contains the global variables for the 3-d grid
        USE MODGRID, ONLY: NGRID

C.........  This module contains arrays for plume-in-grid and major sources
        USE MODELEV, ONLY: LMAJOR, LPING
        
        IMPLICIT NONE

C...........   SUBROUTINE ARGUMENTS
        INTEGER, INTENT(IN) :: NSRC         ! number of sources
        INTEGER, INTENT(IN) :: NMAT1        ! dim 1 for gridding matrix
        INTEGER, INTENT(IN) :: NX( NGRID )  ! no. of sources per cell
        INTEGER, INTENT(IN) :: IX( NMAT1 )  ! list of sources per cell

C...........   Other local variables
        INTEGER C, J, K, SRC, GIDX, CNT   ! counters and indices

***********************************************************************
C   begin body of subroutine SRCGRPCNT

        K = 0
        DO C = 1, NGRID
            DO J = 1, NX( C )
                K = K + 1
                SRC = IX( K )

C.................  Skip elevated sources
                IF( PFLAG ) THEN
                    IF( LMAJOR( SRC ) .OR. LPING( SRC ) ) CYCLE
                END IF

                GIDX = ISRCGRP( SRC )
                CNT = GRPCNT( C, GIDX )
                IF( CNT == 0 ) THEN
                    NSGOUTPUT = NSGOUTPUT + 1
                END IF
                GRPCNT( C, GIDX ) = CNT + 1
            END DO
        END DO
        
        RETURN

        END SUBROUTINE SRCGRPCNT
