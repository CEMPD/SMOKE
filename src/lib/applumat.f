
        SUBROUTINE APPLUMAT( NSRC, NMATX, VAL, NU, IU, CU, VALBYSRC )

C***********************************************************************
C  subroutine APPLUMAT body starts at line < >
C
C  DESCRIPTION:
C      Applies the "ungridding" matrix to gridded data to compute a per-source
C      value of the data.  If the ungridding matrix has no factors for a source,
C      a missing value is returned for that source (i.e., BADVAL3)
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION HISTORY:
C
C***************************************************************************
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
C****************************************************************************

        IMPLICIT NONE

C...........   INCLUDES:
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters

C...........   SUBROUTINE ARGUMENTS
        INTEGER     , INTENT (IN) :: NSRC          ! no. sources
        INTEGER     , INTENT (IN) :: NMATX         ! ugridding factors dimension
        REAL        , INTENT (IN) :: VAL( * )      ! gridded data
        INTEGER     , INTENT (IN) :: NU( NSRC  )   ! no. cells per source
        INTEGER     , INTENT (IN) :: IU( NMATX )   ! grid cell IDs 
        REAL        , INTENT (IN) :: CU( NMATX )   ! ungridding coefficients
        REAL        , INTENT(OUT) :: VALBYSRC( NSRC ) ! ungridded data

C...........   Other local variables
        INTEGER     J, K, S     ! counters and indices

        REAL        RDUM        ! tmp value for summing over cells

        CHARACTER*16 :: PROGNAME = 'APPLUMAT' ! program name

C***********************************************************************
C   begin body of subroutine APPLUMAT

C.........  Apply ungridding matrix 
        K = 0
        DO S = 1, NSRC

            IF( NU( S ) .GT. 0 ) THEN
                RDUM = 0.0
            ELSE
                RDUM = BADVAL3
            END IF

            DO J = 1, NU( S )
                K = K + 1
                RDUM = RDUM + VAL( IU( K ) ) * CU( K )
            END DO

            VALBYSRC( S ) = RDUM

        END DO

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I9, :, 1X ) )
 
        END SUBROUTINE APPLUMAT
