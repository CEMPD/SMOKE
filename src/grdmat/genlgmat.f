
        SUBROUTINE GENLGMAT( GNAME, FDEV, NCOEF, CMAX, CMIN )

C***********************************************************************
C  FUNCTION:
C      Compute gridding matrix going from Lat-Lon grid to a
C      specified modeling grid, converting Kg/(M^2 S) to Kg/(cell S)
C
C      Also supports non-Standard WMO "lat-lon" that claims
C      longitude-range is [0,360], contrary to ISO Standard 6709.
C      Does NOT apply units-conversion (e.g., Kg/(M^2 S) ~> Ton/(cell year) )
C
C  PRECONDITIONS:
C      I/O API 3.2 or later (for MODGCTP/GRID2XY())
C      Initialized MODULE MODGRID
C      setenv GRIDMASK     <path-name> for Lat-Lon emissions-input file
C      setenv GRID_CRO_2D  <path-name> for model-grid GRIDCRO2 file
C
C  REVISION  HISTORY:
C      prototype 2/2016 by Carlie J. Coats, Jr., UNC IE, adapted from
C      "m3tools/mtxcalc.f90"
C
C**************************************************************************
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
C  Version "$Id: genlgmat.f 2016-02-05 15:53:52Z coats $"
C  Copyright (C)  2016 UNC Institute for the Environment.
C
C***********************************************************************

C...........   MODULES for public variables
C...........   This module is the source inventory arrays
        USE MODSOURC, ONLY: CIFIP, CELLID, CSOURC

C.........  This module contains the information about the source category
        USE MODINFO, ONLY: NCHARS, NSRC

C.........  This module contains the global variables for the 3-d grid
        USE MODGRID, ONLY : GDTYP, P_ALP, P_BET, P_GAM, XCENT, YCENT,
     &                      XORIG, YORIG, XCELL, YCELL, VGTYP, VGTOP, 
     &                      NCOLS, NROWS, NGRID

C.........  MODULES for I/O API INTERFACEs, geo-transform codes:
        USE M3UTILIO
        USE MODGCTP

        IMPLICIT NONE

C.........  SUBROUTINE ARGUMENTS and their descriptions:
        CHARACTER(*),         INTENT (IN) :: GNAME  ! gridding mtx logical name
        INTEGER,              INTENT (IN) :: FDEV   ! surg codes report file
        INTEGER,              INTENT(OUT) :: NCOEF  !  no. of gridding coeffs
        INTEGER,              INTENT(OUT) :: CMAX   !  max no. of sources per cell
        INTEGER,              INTENT(OUT) :: CMIN   !  min no. of sources per cell

C.........  PARAMETERS
        REAL*8, PARAMETER    :: PI     = 3.14159265358979324D0 !  pi
        REAL*8, PARAMETER    :: PI180  = PI / 180.0D0          !  degrees-to-radians
        REAL*8, PARAMETER    :: REARTH = 6367333.0D0           !  mean radius of the earth (meters)
        REAL*8, PARAMETER    :: DG2M   = REARTH * PI / 180.0D0 !  degrees-to-meters at equator

C.........  LOCAL VARIABLES and their descriptions:
        INTEGER     I, J, K, L, M, N, C, R, S, V
        INTEGER     GDTYPE, NCOLSE, NROWSE, NGRIDE
        INTEGER     NCOLS2, NROWS2, NGRID2, RATIO
        INTEGER     COL, ROW
        INTEGER     COLMAX, ROWMAX
        INTEGER     COLMIN, ROWMIN
        INTEGER     NSCC
        INTEGER     IOS

        REAL        DDXY, AREA
        REAL*8      XCELL2, YCELL2, DDXE, DDYE, X, Y
        REAL*8      XORIGE, YORIGE, XCELLE, YCELLE    ! for input grid (assumed Lat-Lon)
        REAL*8      P_ALPE, P_BETE, P_GAME, XCENTE, YCENTE
        REAL*8      COSFAC, DXRAT, DYRAT

        INTEGER, ALLOCATABLE :: NS(:)  !  no. srcs per Edgar-cell

        INTEGER, ALLOCATABLE :: NX(:)  !  no. srcs per model-cell
        INTEGER, ALLOCATABLE :: IX(:)  !  source IDs
        REAL   , ALLOCATABLE :: CX(:)  !  gridding coefficients
        INTEGER, ALLOCATABLE :: ICEL( : )
        INTEGER, ALLOCATABLE :: OCEL( : )
        INTEGER, ALLOCATABLE :: ICNT( :,: )
        REAL,    ALLOCATABLE :: FRAC( : )
        REAL,    ALLOCATABLE :: MSFX2( : )  ! square of map-scale factor (single-indexed)
        REAL*8,  ALLOCATABLE :: YLAT( :,: ) ! latitude (degrees)
        REAL*8,  ALLOCATABLE :: XLON( :,: ) ! longitude (degrees)
        REAL*8,  ALLOCATABLE :: XLOC( :,: )
        REAL*8,  ALLOCATABLE :: YLOC( :,: )

        CHARACTER(256)          MESG                  ! tmp line buffer
        
        CHARACTER(NAMLEN3) ::   PROGNAME = 'GENLGMAT'   ! program name

C***********************************************************************
C   begin body of subroutine GENLGMAT

C.........   Get "raw-netCDF" input-grid description from "GRIDMASK" header:

        IF ( .NOT. OPEN3( 'GRIDMASK', FSREAD3, PROGNAME ) ) THEN
            MESG   = 'Could not open "GRIDMASK"'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        ELSE IF ( .NOT. DESC3( 'GRIDMASK' ) ) THEN
            MESG   = 'Could not DESC3( "GRIDMASK" )'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        ELSE IF ( XCELLE * YCELLE .LT. 0.0D0 ) THEN
            MESG   = 'Negative input-grid XCELL or YCELL:  unsupported case'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        ELSE
            NCOLSE = NCOLS3D
            NROWSE = NROWS3D
            NGRIDE = NCOLSE * NROWSE
            XORIGE = XORIG3D
            YORIGE = YORIG3D
            XCELLE = XCELL3D
            YCELLE = YCELL3D
            GDTYPE = GDTYP3D
            P_ALPE = P_ALP3D
            P_BETE = P_BET3D
            P_GAME = P_GAM3D
            XCENTE = XCENT3D
            YCENTE = YCENT3D
        END IF
        
        NSCC = NSRC / NGRIDE

        IF ( MOD( NSRC, NGRIDE ) .NE. 0 ) THEN
            MESG   = 'NSRC is not a multiple of NGRIDE'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

C.........   Open and check dims/proj/grid for GRID_CRO_2D file:

        IF( .NOT.OPEN3( 'GRID_CRO_2D', FSREAD3, PROGNAME ) ) THEN
            MESG   = 'Could not open "GRID_CRO_2D"'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

        ELSE IF( .NOT. DESC3( 'GRID_CRO_2D' ) ) THEN
            MESG   = 'Could not DESC3( "GRID_CRO_2D" )'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

        ELSE 
C.............   Allocate local variables for output modeling domain
C.............   Save modeling-domain description in MODGRID variables:
            NCOLS = NCOLS3D
            NROWS = NROWS3D
            NGRID = NCOLS * NROWS
            XORIG = XORIG3D
            YORIG = YORIG3D
            XCELL = XCELL3D
            YCELL = YCELL3D
            GDTYP = GDTYP3D
            P_ALP = P_ALP3D
            P_BET = P_BET3D
            P_GAM = P_GAM3D
            XCENT = XCENT3D
            YCENT = YCENT3D

        END IF        

        ALLOCATE( MSFX2( NGRID ), STAT=IOS )    !!  ratio  dx*dy/true-Earth-area( cell )
        CALL CHECKMEM( IOS, 'MSFX2', PROGNAME )

        IF ( GDTYP .EQ. LATGRD3 ) THEN

!$OMP       PARALLEL DO DEFAULT( NONE ),
!$OMP&                   SHARED( NCOLS, NROWS, MSFX2 ),
!$OMP&                  PRIVATE( C, R, K )
                DO R = 1, NROWS
                DO C = 1, NCOLS
                    Y = YORIG + ( DBLE( R ) - 0.5d0 )*YCELL
                    K = C + NCOLS*( R - 1 )
                    MSFX2( K ) = 1.0 / ( DG2M**2 * COS( PI180*Y ) )
                END DO
                END DO

            IF ( GDTYPE .EQ. LATGRD3 ) THEN

                DXRAT  = ABS( XCELL / XCELLE )
                DYRAT  = ABS( YCELL / YCELLE )

            ELSE        !!  gdtype not latgrd3:
            
                X = YORIG + 0.5d0 * YCELL
                Y = YORIG + ( DBLE( NROWS ) - 0.5d0 )*YCELL
                IF ( X * Y .LT. 0.0d0 ) THEN
                    COSFAC = 1.0d0      !!  grid crosses equator, where cosfac==1
                ELSE
                    COSFAC = MAX( COS( PI180*X ), COS( PI180*Y ) )
                END IF

                DXRAT  = ABS( DG2M * XCELL * COSFAC / XCELLE )
                DYRAT  = ABS( DG2M * YCELL          / YCELLE )

            END IF      !!  if gdtype is latgrd3, or not
        
        ELSE

            IF( .NOT. READ3('GRID_CRO_2D', 'MSFX2', 1,0,0, MSFX2 ) ) THEN
                MESG  = 'Could not read "MSFX2" from "GRID_CRO_2D"'
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

            END IF

            IF ( GDTYPE .EQ. LATGRD3 ) THEN

C.........   Compute output-grid Lat-Lon:

                ALLOCATE( XLON( NCOLS,NROWS ),
     &                    YLAT( NCOLS,NROWS ),  STAT=IOS )
                CALL CHECKMEM( IOS, 'XLON.YLAT', PROGNAME )

                CALL GRID2XY( GDTYPE,P_ALPE, P_BETE, P_GAME, XCENTE, YCENTE,
     &                        GDTYP, P_ALP,  P_BET,  P_GAM,  XCENT,  YCENT,
     &                        NCOLS, NROWS,  XORIG,  YORIG,  XCELL,  YCELL,
     &                        XLON, YLAT )

C  compute grid refinement ratio (at least 2x finer than input-grid)
C  Need min of cos(lat) on the modeling grid in order to find the
C  smallest meters "dx" for that grid
C  This min will occur either on the top row or on the bottom row,
C  depending upon (N or S) hemisphere coverage of modeling grid...
C  Needs XCELLE,YCELLE converted to meters...

                COSFAC = COS( PI180*YLAT(1,1) )

!$OMP       PARALLEL DO DEFAULT( NONE ),
!$OMP&                   SHARED( NCOLS, NROWS, YLAT ),
!$OMP&                  PRIVATE( C ),
!$OMP&                REDUCTION( MIN: COSFAC )

                DO C = 1, NCOLS
                    COSFAC = MIN( COSFAC, COS( PI180 * YLAT( C,1 ) ),
     &                                    COS( PI180 * YLAT( C,NROWS ) ) )
                END DO

                DXRAT  = ABS( XCELL / ( DG2M * XCELLE * COSFAC ) )
                DYRAT  = ABS( YCELL / ( DG2M * YCELLE ) )

                DEALLOCATE( XLON, YLAT )

            ELSE        !!  gdtype not latgrd3:

                DXRAT  = ABS( XCELL / XCELLE )
                DYRAT  = ABS( YCELL / YCELLE )

            END IF      !!  if gdtype is latgrd3, or not

        END IF      !!  if gdtyp  is latgrd3, or not

        RATIO  = CEILING( 3.0D0 * MAX( DXRAT, DYRAT ) ) !  samples at least 3x per input-cell
        NCOLS2 = RATIO * NCOLS
        NROWS2 = RATIO * NROWS
        NGRID2 = NCOLS2 * NROWS2
        XCELL2 = XCELL / DBLE( RATIO )
        YCELL2 = YCELL / DBLE( RATIO )

C.........   Allocate  and compute Edgar-grid-to-source matrix:
        ALLOCATE( NS( NGRIDE ), STAT = IOS )
        NS = 0

        IF ( IOS .NE. 0 ) THEN
            MESG   = 'Matrix-buffer allocation failure'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF
        
        DO S = 1, NGRIDE    ! assign source ID to first set of Cell IDs
            C = CELLID( S )
            NS( C ) = S
        END DO

C.........   Allocate scratch buffers:
        ALLOCATE( ICNT( NCOLSE,NROWSE ), 
     &            ICEL( NGRID2 ), 
     &            OCEL( NGRID2 ), 
     &            FRAC( NGRID2 ),
     &            XLOC( NCOLS2, NROWS2 ), 
     &            YLOC( NCOLS2, NROWS2 ), STAT = IOS )

        IF ( IOS .NE. 0 ) THEN
            MESG   = 'Work-buffer allocation failure'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

C.........   Compute fractions:
C.........   First, find input-grid cell-centers for refined modeling grid:
        CALL GRID2XY( GDTYPE,  P_ALPE, P_BETE, P_GAME, XCENTE, YCENTE,
     &                GDTYP,   P_ALP,  P_BET,  P_GAM,  XCENT,  YCENT,
     &                NCOLS2, NROWS2, XORIG, YORIG, XCELL2, YCELL2,
     &                XLOC, YLOC )

C  Fix-up for WMO-style [0,360] input-grid longitudes

        IF ( GDTYPE .EQ. LATGRD3 ) THEN
            IF ( XORIGE + DBLE(NCOLSE)*XCELLE > 180.0D0 ) THEN

                MESG = 'This file seems to use WMO-style "longitude" ' //
     &                 'that violates ISO Standard 6709'
                CALL M3MESG( MESG )

!$OMP           PARALLEL DO
!$OMP&             DEFAULT( NONE ),
!$OMP&              SHARED( NCOLS2, NROWS2, XLOC ),
!$OMP&             PRIVATE( C, R )

                DO R = 1, NROWS2
                DO C = 1, NCOLS2
                    XLOC( C,R ) = MOD( XLOC( C,R ) + 360.0D0 , 360.0D0 )    !  now in range [0,360]
                END DO
                END DO

            END IF      !  if WMO-style "lon" or not
        END IF      !  if LATLON input grid

!$OMP       PARALLEL DO
!$OMP&         DEFAULT( NONE ),
!$OMP&          SHARED( NCOLSE, NROWSE, ICNT ),
!$OMP&         PRIVATE( C, R )

        DO  R = 1, NROWSE
        DO  C = 1, NCOLSE
            ICNT( C,R ) = 0
        END DO
        END DO

        DDXE = 1.0D0 / XCELLE
        DDYE = 1.0D0 / YCELLE
        DDXY = 1.0   / FLOAT( RATIO**2 )
        N    = 0

        DO  R = 1, NROWS         ! traverse output grid:
        DO  C = 1, NCOLS         ! serial loop net:  dependency on ICNT, N

            COLMAX = 0
            ROWMAX = 0
            COLMIN = NCOLSE + 1
            ROWMIN = NROWSE + 1

            DO  J = 1, RATIO       ! loop nest sub-sampling this model-grid cell
            DO  I = 1, RATIO

                K = RATIO*( C - 1 ) + I
                L = RATIO*( R - 1 ) + J
                X = XLOC( K,L ) - XORIGE
                Y = YLOC( K,L ) - YORIGE
                IF ( X .GE. 0.0D0  .AND. Y .GE. 0.0D0 ) THEN
                    COL = 1 + INT( DDXE * X )
                    ROW = 1 + INT( DDYE * Y )
                    IF ( COL .LE. NCOLSE .AND. ROW .LE. NROWSE ) THEN
                        ICNT( COL,ROW ) = ICNT( COL,ROW ) + 1
                        IF ( COL .GT. COLMAX ) COLMAX = COL
                        IF ( COL .LT. COLMIN ) COLMIN = COL
                        IF ( ROW .GT. ROWMAX ) ROWMAX = ROW
                        IF ( ROW .LT. ROWMIN ) ROWMIN = ROW
                    END IF
                END IF

            END DO      !  end loop on sub-sampled cols
            END DO      !  end loop on sub-sampled rows

            DO ROW = ROWMIN, ROWMAX     !  loop nest on impacted input-grid cells
            DO COL = COLMIN, COLMAX
                IF ( ICNT( COL,ROW ) .GT. 0 ) THEN
                    N = N + 1
                    FRAC( N ) = DDXY * FLOAT( ICNT( COL,ROW ) )
                    ICEL( N ) = COL  + ( ROW - 1 ) * NCOLSE  ! input  cell #
                    OCEL( N ) = C    + ( R   - 1 ) * NCOLS   ! output cell #
                    ICNT( COL,ROW ) = 0                      ! reset to zero for next pass
                END IF
            END DO
            END DO

        END DO
        END DO          !  end traversal of input grid

        NCOEF = N * NSCC

        IF ( N .GT. NGRID2 ) THEN
            WRITE( MESG, 94010 )
     &          'ERROR:  Overflow allocation. Estimated', NGRID2, ':: Actual:', N
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

        ELSE IF ( N .EQ. 0 ) THEN
            MESG = 'ERROR:  NO INTERSECTION FOUND:  Number of coeffs=0'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )     !  if input-grid is global, should never happen

        END IF

        WRITE( MESG, 94010 ) 'Number of coeffs', N
        CALL M3MESG( MESG )

C.........   Allocate and organize the output sparse matrix:
        ALLOCATE( NX( NGRID ),
     &            IX( NCOEF ),
     &            CX( NCOEF ), STAT = IOS )

        IF ( IOS .NE. 0 ) THEN
            WRITE( MESG, 94010 ) 'ERROR: Buffer allocation failure. STATUS=', IOS
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

        K    = 0
        CMAX = 0
        CMIN = NCOEF + 1

        DO I = 1, NGRID         !  SERIAL LOOP on model-grid:  dependency on K, L

            L    = 0
            AREA = XCELL * YCELL / MSFX2( I )    !  true cell-area (M^2)

            DO J = 1, NCOEF

                IF ( OCEL( J ) .EQ. I ) THEN
                    C = ICEL( J )                   !!  EDGAR cell-ID 1 ... NGRIDE
                    DO N = 1, NSCC 
                        K = K + 1
                        L = L + 1
                        IX( K ) = NS( C ) + NGRIDE * ( N-1 )
                        CX( K ) = FRAC( J ) * AREA
                    END DO
                END IF

            END DO          !  end loop on coeffs J

            NX( I ) = L
            CMAX = MAX( CMAX, NX( I ) )
            CMIN = MIN( CMIN, NX( I ) )

        END DO

C.........  Output output gridding matrix

        CALL WRITELLG( GNAME, NGRID, NCOEF, NSRC, NX, IX, CX )

C.........  Write output surrogates codes
        DO S = 1, NGRIDE
            WRITE( FDEV,93360 ) S, 0, 0
        END DO

        DEALLOCATE( NS, ICNT, ICEL, OCEL, FRAC, XLOC, YLOC, MSFX2, NX, IX, CX )

        RETURN

C******************  FORMAT  STATEMENTS   ******************************
C...........   Formatted file I/O formats............ 93xxx

93360   FORMAT( I8, 1X, I8, 1X, I8 )

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I10, :, 1X ) )

      END SUBROUTINE GENLGMAT

C---------------------------------------------------------------------
C     Assemble contiguous matrix in memory and write it to GNAME:
C     (needs to be a separate routine because of type-checking)
C---------------------------------------------------------------------
      SUBROUTINE  WRITELLG( GNAME, NGRID, NCOEF, NSRC, NX, IX, CX )

        USE M3UTILIO
        USE MODGRID, ONLY : NCOLS, NROWS

C.........   ARGUMENTS and their descriptions:
        CHARACTER(*),   INTENT(IN) :: GNAME  ! gridding mtx logical name
        INTEGER,        INTENT(IN) :: NGRID, NCOEF, NSRC
        INTEGER,        INTENT(IN) :: NX(NGRID)  !!  no. srcs per cell
        INTEGER,        INTENT(IN) :: IX(NCOEF)  !!  src IDs
        INTEGER,        INTENT(IN) :: CX(NCOEF)  !!  gridding coefficients

        INTEGER     MTXBUF( NGRID + 2*NCOEF )
        INTEGER     I, J
        CHARACTER(256) MESG
C***********************************************************************

        I = 0
        DO J = 1, NGRID
            I = I + 1
            MTXBUF(I) = NX(J)
        END DO

        DO J = 1, NCOEF
            I = I + 1
            MTXBUF(I) = IX(J)
        END DO

        DO J = 1, NCOEF
            I = I + 1
            MTXBUF(I) = CX(J)
        END DO

C.........  Initialize header by using the DESC3(GRID_CRO_2D) results.
        NCOLS3D = NCOEF
        NROWS3D = NGRID
        NTHIK3D = NSRC
        NVARS3D = 1
        FTYPE3D = SMATRX3
        VNAME3D( 1 ) = 'AGRDMAT'
        UNITS3D( 1 ) = 'n/a'
        VDESC3D( 1 ) = 'Area source re-gridding coefficients'
        VTYPE3D( 1 ) = M3REAL
        FDESC3D      = ' '
        FDESC3D( 1 ) = 'Area source re-gridding coefficients'
        FDESC3D( 2 ) = '/FROM/ GENLGMAT/WRITELLG'
        FDESC3D( 3 ) = '/VERSION/ '         !  // VERCHAR( CVSW )
        FDESC3D( 4 ) = '/GDESC/ '           !  // GDESC
        WRITE( FDESC3D(5), 94010 ) '/NCOLS3D/ ', NCOLS
        WRITE( FDESC3D(6), 94010 ) '/NROWS3D/ ', NROWS

        FDESC3D( 11 ) = '/INVEN FROM/ '     !  // INVPROG
        FDESC3D( 12 ) = '/INVEN VERSION/ '  !  // INVVERS
        
        IF ( .NOT.OPEN3( GNAME, FSUNKN3, 'GENLGMAT/WRITELL' ) ) THEN
            MESG = 'Error opening GRIDDING MATRIX file.'
            CALL M3EXIT( 'GENLGMAT/WRITELL', 0, 0, MESG, 2 )
        ELSE IF( .NOT. WRITE3( GNAME, 'ALL', 0, 0, MTXBUF ) ) THEN
            MESG = 'Error writing GRIDDING MATRIX file.'
            CALL M3EXIT( 'GENLGMAT/WRITELL', 0, 0, MESG, 2 )
        END IF
        
        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

      END SUBROUTINE WRITELLG 
