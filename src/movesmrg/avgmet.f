
        SUBROUTINE AVGMET( NSRC, TSTEP ) 

C***********************************************************************
C  subroutine body starts at line 78
C
C  DESCRIPTION:
C       Averages hourly meteorology data based on number of sources
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:  none
C
C  REVISION  HISTORY:
C     10/01: Created by C. Seppanen
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

C...........   MODULES for public variables

C...........   This module is the derived meteorology data for emission factors
        USE MODMET, ONLY: TKHOUR, RHHOUR, NDAYSRC
        
        IMPLICIT NONE

C...........   INCLUDES

        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters

C...........   EXTERNAL FUNCTIONS 
        CHARACTER(2) CRLF
        INTEGER      FIND1FIRST

        EXTERNAL     CRLF, FIND1FIRST

C...........   SUBROUTINE ARGUMENTS
        INTEGER, INTENT   (IN) :: NSRC                   ! no. sources
        INTEGER, INTENT   (IN) :: TSTEP                  ! current time step

C...........   Other local variables
        INTEGER I, J, K, L, S                ! counters and indices                     
        
        INTEGER IOS                       ! I/O status

        CHARACTER(300) MESG           ! message buffer

        CHARACTER(16) :: PROGNAME = 'AVGMET' ! program name

C***********************************************************************
C   begin body of subroutine AVGMET
        
C.........  Loop through all counties
        DO S = 1, NSRC
        
C.............  Skip sources with no days; this can happen when the
C               gridding surrogates do not contain data for all counties
            IF( NDAYSRC( S,TSTEP ) == 0 ) CYCLE

            TKHOUR( S,TSTEP ) = TKHOUR( S,TSTEP ) / NDAYSRC( S,TSTEP )
            RHHOUR( S,TSTEP ) = RHHOUR( S,TSTEP ) / NDAYSRC( S,TSTEP )

        END DO
        
        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I6, :, 1X ) )

94020   FORMAT( A, 4( 1X, F8.2, 1X, A ) )
 
        END SUBROUTINE AVGMET
