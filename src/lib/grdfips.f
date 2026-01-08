
        SUBROUTINE GRDFIPS( NSRC, CNTY, VAL, VALBYSRC, TFLAG )

C***********************************************************************
C  subroutine APPLUMAT body starts at line 72
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
C       Updated with USE M3UTILIO by Huy Tran UNC-IE on 2026-01
C****************************************************************************

C...........   Modules for public variables
C...........   This module contains the gridding surrogates tables
        USE M3UTILIO

        USE MODSURG, ONLY: NSRGFIPS, SRGFIPS, NCELLS, FIPCELL
     
C.........  This module contains the global variables for the 3-d grid
        USE MODGRID, ONLY: NCOLS

C.........  This module is the derived meteorology data for emission factors
        USE MODMET, ONLY: MINTSRC, MAXTSRC, MAXTDAY, MINTDAY
 
        IMPLICIT NONE

C...........   INCLUDES:
C        INCLUDE 'PARMS3.EXT'    !  I/O API parameters

C...........   EXTERNAL FUNCTIONS and their descriptions:        
C       INTEGER      FINDC
C        EXTERNAL     FINDC

C...........   SUBROUTINE ARGUMENTS
        INTEGER     , INTENT (IN) :: NSRC             ! no. sources
        CHARACTER(*), INTENT (IN) :: CNTY( NSRC  )    ! no. counties
        REAL        , INTENT (IN) :: VAL( * )         ! gridded data
        REAL        , INTENT(OUT) :: VALBYSRC( NSRC ) ! ungridded data
        LOGICAL     , INTENT (IN) :: TFLAG            ! processing temp. variable

C...........   Other local variables
        INTEGER     C, I, J, K, L, S  ! counters and indices
        INTEGER     COL            ! subgrid column number
        INTEGER     ROW            ! subgrid row number

        REAL        CNTOT          ! tmp value for summing over cells
        REAL        TEMPVAL        ! temperature value in Farenheight

        CHARACTER(16) :: PROGNAME = 'GRDFIPS' ! program name

C***********************************************************************
C   begin body of subroutine APPLUMAT

C.........  Apply ungridding matrix from a (possible) subgrid to data on base 
C           grid.  If no subgrid, then XOFF and YOFF will be 1 and no problem.
        DO S = 1, NSRC

            L = FINDC( CNTY( S ), NSRGFIPS, SRGFIPS )

            IF( L < 1 ) CYCLE
            
            IF( NCELLS( L ) .GT. 0 ) THEN
                CNTOT = 0.0
            ELSE
                CNTOT = BADVAL3
            END IF
            
            K = 0
            DO J = 1, NCELLS( L )
                K = K + 1

C.................  Get column and row from subgrid
                C = FIPCELL( J,L )

                ROW = C / NCOLS          ! note: integer math
                IF( MOD( C, NCOLS ) .GT. 0. ) ROW = ROW + 1
                COL = C - ( ROW-1 ) * NCOLS

C bbh           C = (ROW-1)*NCOLS + COL   ! org cell equations using row,col

C.................  Convert K to F
                IF( TFLAG ) THEN
                    TEMPVAL = 1.8 * VAL( C ) - 459.67

                    MAXTSRC( S ) = MAX( TEMPVAL, MAXTSRC( S ) )
                    MINTSRC( S ) = MIN( TEMPVAL, MINTSRC( S ) )

                    MAXTDAY( S ) = MAX( TEMPVAL, MAXTDAY( S ) ) 
                    MINTDAY( S ) = MIN( TEMPVAL, MINTDAY( S ) ) 
                END IF

                CNTOT = CNTOT + VAL( C )
                
            END DO

            VALBYSRC( S ) = CNTOT / K    ! averaged

        END DO

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I9, :, 1X ) )
 
        END SUBROUTINE GRDFIPS
