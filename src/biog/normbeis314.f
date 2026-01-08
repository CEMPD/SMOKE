
        SUBROUTINE NORMBEIS312(  )

C***********************************************************************
C
C  DESCRIPTION:  Produces normalized biogenic emissions for use with
C                SMOKE-BEIS3 v3.12.
C
C  SUBROUTINES AND FUNCTIONS CALLED: RDB3FAC
C
C  REVISION  HISTORY: 4/00 Prototype, Jeff Vukovich
C                     1/03 changes to NO, George Pouliot
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
C       Updated with USE M3UTILIO by Huy Tran UNC-IE on 2026-01
C***********************************************************************

C.........  MODULES for public variables
C.........  This module contains biogenic variables
        USE M3UTILIO

        USE MODBEIS3, ONLY: NVEG, VEGID, AVGEMIS, AVGLAI, NOEMIS, 
     &                      EMFAC, LAI, SLW, WFAC, LFBIO

        USE MODGRDLIB
 
        IMPLICIT NONE

C.........  INCLUDES
C        INCLUDE 'PARMS3.EXT'      ! I/O API constants
C        INCLUDE 'FDESC3.EXT'      ! I/O API file description data structure
C        INCLUDE 'IODECL3.EXT'     ! I/O API function declarations
        INCLUDE 'B3V14DIMS3.EXT'     ! BEIS3-related declarations
        INCLUDE 'IOSTRG3.EXT'     ! I/O API string declarations
  
C.........  EXTERNAL FUNCTIONS and their descriptions
        INTEGER         GETFLINE
C       INTEGER         PROMPTFFILE
C       CHARACTER(16)   PROMPTMFILE
        CHARACTER(16)   VERCHAR
C       LOGICAL         ENVYN

C        EXTERNAL        GETFLINE, PROMPTFFILE, PROMPTMFILE, VERCHAR, ENVYN
        EXTERNAL     GETFLINE, VERCHAR

C.........  ARGUMENTS and their descriptions
C       CHARACTER(50), INTENT(IN) :: CVSW    ! CVS release tag

C.........  LOCAL VARIABLES and their descriptions
        INTEGER         B, C, R, I, J, K, L, M, N ! loop counters and subscripts
        INTEGER         IOS     !  I/O status result

        INTEGER         FDEV    !  unit number for emissions factor file
        INTEGER         LDEV    !  unit number for log file

        CHARACTER(16)   ENAME   !  logical name for normalized emissions output
        CHARACTER(16)   GNAME1  !  unit number for gridded land use file
        CHARACTER(16)   GNAME2  !  unit number for gridded land use file
        CHARACTER(16)   DOTNAME !  logical name for gridded 2D dot file
        CHARACTER(16)   GRDNM   !  grid name
        CHARACTER(16), ALLOCATABLE  :: LUNAMA( : )  ! land use type names
        CHARACTER(16), ALLOCATABLE  :: LUNAMB( : )  ! land use type names
        CHARACTER(16)   LUNAM   ! temporary tag for land use names    
        CHARACTER(16)   GNAMET  !  unit number for gridded land use totals file

        CHARACTER(256)  MESG    !  message buffer for M3EXIT()
        CHARACTER(5)    BTMP    ! temporary tag for naming output variables

        INTEGER, ALLOCATABLE      :: LUINDX( : )   ! index for land use types
        INTEGER         NCOLS   ! no. of grid columns
        INTEGER         NROWS   ! no. of grid rows
        INTEGER         NCDOT   ! no. of columns in dot file
        INTEGER         NRDOT   ! no. of rows in dot file
        INTEGER         NVARSA  ! no. of land use types in land use file A
        INTEGER         NVARSB  ! no. of land use types in land use file B
        INTEGER         NVARST  ! total land use types = A + B
        INTEGER         IFOUND  ! used for checking land use vs. emis facs
        INTEGER         IUSDEC  ! USGS decid forest
        INTEGER         IUSBRD  ! USGS evbrdleaf
        INTEGER         IUSCON  ! USGS coniferfor
        INTEGER         IUSMIX  ! USGS mixed forest                         
        INTEGER         IUSSHR  ! USGS shrubland
        INTEGER         IUSCGS  ! USGS Cropgrass
        INTEGER         IUSCWD  ! USGS Cropwoodland
        INTEGER         IUSCDY  ! USGS Drycrop
        INTEGER         IUSCIR  ! USGS Irrcrop
        INTEGER         IALFAL  ! Alfalfa
        INTEGER         IBARLE  ! Barley
        INTEGER         ICORN   ! Corn
        INTEGER         ICOTTO  ! Cotton
        INTEGER         IGRASS  ! Grass
        INTEGER         IHAY    ! Hay
        INTEGER         IMISCC  ! Misc_crop
        INTEGER         IOATS   ! Oats
        INTEGER         IPASTU  ! Pasture
        INTEGER         IPEANU  ! Peanuts
        INTEGER         IPOTAT  ! Potatoes
        INTEGER         IRICE   ! Rice
        INTEGER         IRYE    ! Rye
        INTEGER         ISORGH  ! Sorghum
        INTEGER         ISOYBE  ! Soybeans
        INTEGER         ITOBAC  ! Tobacco
        INTEGER         IWHEAT  ! Wheat

        REAL, ALLOCATABLE    :: XVALS ( :, : )     ! x values for grid cell boundaries
        REAL, ALLOCATABLE    :: YVALS ( :, : )     ! y values for grid cell boundaries

        REAL, ALLOCATABLE    :: LUSE ( :, :, :  )  ! BELD3 land use data
        REAL, ALLOCATABLE    :: FIA ( :, : )       ! Forest inventory data

        REAL, ALLOCATABLE    ::   SUMEM( : )       ! Summer emissions 
        REAL, ALLOCATABLE    ::   SUMEMW( : )      ! Winter emissions
        REAL, ALLOCATABLE    ::   NOEM( : )        ! NO emissions 
        REAL, ALLOCATABLE    ::   SUMLAI( : )      ! Summer LAIs
        REAL, ALLOCATABLE    ::   SUMLAIW( : )     ! Winter LAIs
        REAL  VEGAREA                              ! Veg. area for 
        REAL  TOTFOR                               ! USGS forest
        REAL  EFTMP                                ! Emission factors
        
        REAL  XCELL                                ! x cell size
        REAL  YCELL                                ! y cell size

        DOUBLE PRECISION, ALLOCATABLE :: PRCNT2KM2( :, : ) ! Prcnt to km**2

        LOGICAL         USE_SHRUB
        LOGICAL         VFLAG                      ! variable grid flag
        LOGICAL       :: EFLAG = .FALSE.           ! Error flag

        CHARACTER(16) :: PROGNAME = 'NORMBEIS312'  ! Program name

C***********************************************************************
C   begin body of subroutine NORMBEIS312

C.........  Check for variable grid data
        VFLAG = ENVYN( 'USE_VARIABLE_GRID',
     &                 'Use variable grid definition',
     &                 .FALSE., IOS )
   
C.........  Get file name; open emission factors file
        FDEV = PROMPTFFILE( 
     &           'Enter logical name for EMISSION FACTORS file',
     &           .TRUE., .TRUE., 'B3FAC', PROGNAME )

C.........  Open gridded landuse files 
        GNAMET = PROMPTMFILE(
     &           'Enter logical name for GRIDDED LANDUSE totals file',
     &           FSREAD3, 'BELD3_TOT', PROGNAME )

        IF ( .NOT. DESC3( GNAMET ) ) THEN
            MESG = 'Could not get description of file "' //
     &             TRIM( GNAMET ) // '"'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

C.........  Initialize grid definition
        CALL CHKGRID( GNAMET, 'GRID' , 0, EFLAG )

C.........  Grid cell resolution assumed to be in meters
C           Input land use will be in percentages
C           Compute conversion factor to go from percentages
C           to km**2
        NCOLS = NCOLS3D
        NROWS = NROWS3D
        GRDNM = GDNAM3D

        IF ( VFLAG ) THEN
C.............  Open GRIDDOT2D file
            DOTNAME = PROMPTMFILE(
     & 'Enter logical name for DOT-POINT SURFACE GRID file',
     &                FSREAD3, 'GRID_DOT_2D', PROGNAME )

            If ( .NOT. DESC3( DOTNAME ) ) THEN
                MESG = 'Could not get description of file "' //
     &                 TRIM( DOTNAME ) // '"'
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            ENDIF

C.............  Check grid definition
            CALL CHKGRID( DOTNAME, 'DOT', 0, EFLAG )
        
            IF ( EFLAG ) THEN
                MESG = 'Grid in file "' // TRIM( DOTNAME ) //
     &                 '" does not match previously set grid.'
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF
        
            NCDOT = NCOLS + 1
            NRDOT = NROWS + 1

C.............  Allocate memory for grid cell coordinates
            ALLOCATE( XVALS( NCDOT, NRDOT ), STAT=IOS )
            CALL CHECKMEM( IOS, 'XVALS', PROGNAME )
            ALLOCATE( YVALS( NCDOT, NRDOT ), STAT=IOS )
            CALL CHECKMEM( IOS, 'YVALS', PROGNAME )

C.............  Read grid cell coordinates
            IF( .NOT. READ3( DOTNAME, 'LON', 1, 0, 0, XVALS ) ) THEN
                MESG = 'Could not read LON from file "' //
     &                 TRIM( DOTNAME ) // '"'
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF

            IF( .NOT. READ3( DOTNAME, 'LAT', 1, 0, 0, YVALS ) ) THEN
                MESG = 'Could not read LAT from file "' //
     &                 TRIM( DOTNAME ) // '"'
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF

C.............  Convert coordinates to map projection units
            CALL CONVRTXY( NCDOT, NRDOT, GDTYP3D, GRDNM,
     &                     P_ALP3D, P_BET3D, P_GAM3D,
     &                     XCENT3D, YCENT3D, XVALS, YVALS )

C.............  Calculate cell size for each cell and conversion factor
            ALLOCATE( PRCNT2KM2( NCOLS, NROWS ), STAT=IOS )
            CALL CHECKMEM( IOS, 'PRCNT2KM2', PROGNAME )
            DO I = 1, NCOLS
                DO J = 1, NROWS
                    XCELL = ABS( XVALS( I + 1, J ) - XVALS( I, J ) )
                    YCELL = ABS( YVALS( I, J + 1 ) - YVALS( I, J ) )
                    PRCNT2KM2( I,J ) = XCELL * YCELL * 1E-08
                END DO
            END DO
        ELSE
            ALLOCATE( PRCNT2KM2( 1, 1 ), STAT=IOS )
            CALL CHECKMEM( IOS, 'PRCNT2KM2', PROGNAME )
            PRCNT2KM2( 1, 1 ) = XCELL3D * YCELL3D * 1E-08
        END IF

        GNAME1 = PROMPTMFILE(
     &           'Enter logical name for GRIDDED LANDUSE A file',
     &           FSREAD3, 'BELD3_A', PROGNAME )

        IF ( .NOT. DESC3( GNAME1 ) ) THEN
            MESG = 'Could not get description of file "' //
     &             TRIM( GNAME1 ) // '"'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        ENDIF

C.........  Check grid definition
        CALL CHKGRID( GNAME1, 'GRID' , 0, EFLAG )

C.........  If grid definition does not match first landuse file then stop
        IF ( EFLAG ) THEN
            MESG = 'Grid in file "' // TRIM( GNAME1 ) //
     &             '" does not match previously set grid.'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

        NVARSA = NVARS3D

C.........  Store landuse variable names from first file
        ALLOCATE( LUNAMA( NVARSA ), STAT=IOS )
        CALL CHECKMEM( IOS, 'LUNAMA', PROGNAME )

        DO I = 1, NVARSA
            LUNAMA ( I ) = VNAME3D ( I ) 
        END DO
 
        GNAME2 = PROMPTMFILE(
     &           'Enter logical name for GRIDDED LANDUSE B file',
     &           FSREAD3, 'BELD3_B', PROGNAME )

        IF ( .NOT. DESC3( GNAME2 ) ) THEN
            MESG = 'Could not get description of file "' //
     &             TRIM( GNAME2 ) // '"'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

C.........  Check grid definition
        CALL CHKGRID( GNAME2, 'GRID' , 0, EFLAG )

C.........  If grid definition does not match first landuse file then stop
        IF ( EFLAG ) THEN
            MESG = 'Grid in file "' // TRIM( GNAME2 ) //
     &             '" does not match previously set grid.'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        ENDIF

C.........  Store landuse variable names from second file
        NVARSB = NVARS3D

        ALLOCATE( LUNAMB( NVARSB ), STAT=IOS )
        CALL CHECKMEM( IOS, 'LUNAMB', PROGNAME )

        DO I = 1, NVARSB
            LUNAMB ( I ) = VNAME3D ( I ) 
        END DO

C.........  Set up header variables for output file B3GRD
        NROWS3D = NROWS
        NCOLS3D = NCOLS
        GDNAM3D = GRDNM
        FTYPE3D = GRDDED3

        SDATE3D = 0       !  n/a
        STIME3D = 0       !  n/a
        TSTEP3D = 0       !  time independent
        NVARS3D = ( (NSEF-1) + NLAI ) * NSEASONS + NNO    ! special treatment of NO
        NLAYS3D = 1
        NTHIK3D = 1
        VGTYP3D = IMISS3
        VGTOP3D = AMISS3

        FDESC3D = ' '   ! array

        FDESC3D( 1 ) = 'BEIS3 normalized emissions values.'
        FDESC3D( 2 ) = '/FROM/ '    // PROGNAME
        FDESC3D( 3 ) = '/VERSION/ ' // VERCHAR( CVSW )
        FDESC3D( 4 ) = '/LANDUSE/ SMOKE TOOL '

        I = 0

C.........  Set up variable names and output units
        DO M = 1, NSEASONS 
            DO B = 1, NSEF
            
                BTMP = BIOTYPES( B ) 

C.................  Handle types except NO
                IF( TRIM( BTMP ) /= 'NO' ) THEN
                    I = I + 1
                    
                    VNAME3D( I ) = 'AVG_' // TRIM( BTMP ) // SEASON( M )
                    VDESC3D( I ) = 'normalized emissions'
                    UNITS3D( I ) = 'gramsC/hour'
                    VTYPE3D( I ) = M3REAL
                END IF
            END DO
            
            DO N = 1, NLAI
            
                BTMP = LAITYPES( N )
                I = I + 1
                
                VNAME3D( I ) = 'LAI_' // TRIM( BTMP ) // SEASON( M )
                VDESC3D( I ) = 'normalized emissions'
                UNITS3D( I ) = 'index'
            END DO
        END DO

C.............  Handle NO types (not dependent on season)            
        DO B = 1, NSEF
        
            BTMP = BIOTYPES( B )
            
            IF( TRIM( BTMP ) == 'NO' ) THEN
            
                I = I + 1
                VNAME3D( I ) = 'AVG_' // TRIM( BTMP ) // 'AG_GROW'
                VDESC3D( I ) = 'normalized emissions for NO AG_GROW'
                UNITS3D( I ) = 'gramsN/hour'
                VTYPE3D( I ) = M3REAL
                
                I = I + 1
                VNAME3D( I ) = 'AVG_' // TRIM( BTMP ) // 'AG_NONGROW'
                VDESC3D( I ) = 'normalized emissions for NO AG_NONGROW'
                UNITS3D( I ) = 'gramsN/hour'
                VTYPE3D( I ) = M3REAL

                I = I + 1
                VNAME3D( I ) = 'AVG_' // TRIM( BTMP ) // 'NONAG'
                VDESC3D( I ) = 'normalized emissions for NO NONAG'
                UNITS3D( I ) = 'gramsN/hour'
                VTYPE3D( I ) = M3REAL

            END IF
        END DO

C.........  Open output file
        ENAME = PROMPTMFILE(  
     &        'Enter logical name for NORMALIZED emissions output file',
     &        FSUNKN3, 'B3GRD', PROGNAME )

C.........  Get length of BFAC file
        NVEG = GETFLINE( FDEV, 'Emissions factor file' )

C.........  Allocate memory for emission factor variables   
        ALLOCATE( VEGID ( NVEG ), STAT=IOS )
        CALL CHECKMEM( IOS, 'VEGID', PROGNAME )

        ALLOCATE( EMFAC ( NVEG, NSEF ), STAT=IOS )
        CALL CHECKMEM( IOS, 'EMFAC', PROGNAME )

        ALLOCATE( LAI ( NVEG ), STAT=IOS )
        CALL CHECKMEM( IOS, 'LAI', PROGNAME )

        ALLOCATE( SLW ( NVEG ), STAT=IOS )
        CALL CHECKMEM( IOS, 'SLW', PROGNAME )

        ALLOCATE( WFAC ( NVEG ), STAT=IOS )
        CALL CHECKMEM( IOS, 'WFAC', PROGNAME )

        ALLOCATE( LFBIO ( NVEG ), STAT=IOS )
        CALL CHECKMEM( IOS, 'LFBIO', PROGNAME )

C.........  Read emissions factor file
        MESG = 'Reading emissions factor file...'
        CALL M3MSG2( MESG )
        
        WRITE( MESG,94010 ) 'Number of landuse types in factor file: ',
     &                      NVEG
        CALL M3MSG2( MESG )

C.........  This routine reads in emission factors 
        CALL RDB3FAC( .FALSE., NSEF, FDEV, NVEG, VEGID, LAI, LFBIO,
     &                WFAC, SLW, EMFAC ) 

C.........  Total number of land use types (should be 230)
        NVARST = NVARSA + NVARSB

        ALLOCATE( LUINDX ( NVARST ), STAT=IOS )
        CALL CHECKMEM( IOS, 'LUINDX', PROGNAME )
        ALLOCATE( LUSE ( NCOLS, NROWS, NVEG ), STAT=IOS )
        CALL CHECKMEM( IOS, 'LUSE', PROGNAME )

        LUINDX = -9   ! array

C.........  Check to see if there are emissions factors for all landuse types
        DO I = 1, NVARST

            IFOUND = 0

            IF ( I .LE. NVARSA ) THEN
                LUNAM = LUNAMA( I ) 
            ELSE
                K = I - NVARSA
                LUNAM = LUNAMB( K ) 
            END IF 
 
            DO J = 1, NVEG

                IF ( VEGID ( J ) .EQ. LUNAM ) THEN
                    IFOUND = 1  
                    LUINDX( I ) = J
                END IF

C.................  Store vegids for certain USGS categories for use later
                IF ( VEGID( J ) .EQ. 'USGS_decidforest' ) IUSDEC = J
                IF ( VEGID( J ) .EQ. 'USGS_evbrdleaf  ' ) IUSBRD = J
                IF ( VEGID( J ) .EQ. 'USGS_coniferfor ' ) IUSCON = J
                IF ( VEGID( J ) .EQ. 'USGS_mxforest   ' ) IUSMIX = J
                IF ( VEGID( J ) .EQ. 'USGS_shrubland  ' ) IUSSHR = J

                IF ( VEGID( J ) .EQ. 'USGS_cropgrass  ' ) IUSCGS = J
                IF ( VEGID( J ) .EQ. 'USGS_cropwdlnd  ' ) IUSCWD = J

                IF ( VEGID( J ) .EQ. 'USGS_drycrop    ' ) IUSCDY = J
                IF ( VEGID( J ) .EQ. 'USGS_irrcrop    ' ) IUSCIR = J

                IF ( VEGID( J ) .EQ. 'Alfalfa         ' ) IALFAL = J
                IF ( VEGID( J ) .EQ. 'Barley          ' ) IBARLE = J
                IF ( VEGID( J ) .EQ. 'Corn            ' ) ICORN  = J
                IF ( VEGID( J ) .EQ. 'Cotton          ' ) ICOTTO = J
                IF ( VEGID( J ) .EQ. 'Grass           ' ) IGRASS = J
                IF ( VEGID( J ) .EQ. 'Hay             ' ) IHAY   = J
                IF ( VEGID( J ) .EQ. 'Misc_crop       ' ) IMISCC = J
                IF ( VEGID( J ) .EQ. 'Oats            ' ) IOATS  = J
                IF ( VEGID( J ) .EQ. 'Pasture         ' ) IPASTU = J
                IF ( VEGID( J ) .EQ. 'Peanuts         ' ) IPEANU = J
                IF ( VEGID( J ) .EQ. 'Potatotes       ' ) IPOTAT = J
                IF ( VEGID( J ) .EQ. 'Rice            ' ) IRICE  = J
                IF ( VEGID( J ) .EQ. 'Rye             ' ) IRYE   = J
                IF ( VEGID( J ) .EQ. 'Sorghum         ' ) ISORGH = J
                IF ( VEGID( J ) .EQ. 'Soybeans        ' ) ISOYBE = J
                IF ( VEGID( J ) .EQ. 'Tobacco         ' ) ITOBAC = J
                IF ( VEGID( J ) .EQ. 'Wheat           ' ) IWHEAT = J

            END DO

            IF( IFOUND .EQ. 0 ) THEN
               MESG =   LUNAM // 
     &                ' does NOT have emissions factors in BFAC file'
               CALL M3WARN( PROGNAME, 0, 0, MESG )
            END IF 
         
        END DO

C.........  Read the gridded landuse from the landuse files
        DO M = 1, NVARST 
     
            N = LUINDX( M ) 

            IF( N > 0 ) THEN

                IF( M <= NVARSA ) THEN
                    MESG = 'Reading ' // LUNAMA( M )
                    CALL M3MESG( MESG )

                    IF( .NOT. READ3( GNAME1, LUNAMA( M ), 1, 0, 0, 
     &                               LUSE( 1, 1, N ) ) ) THEN
                        MESG = 'Could not find "' // LUNAMA( M ) // 
     &                         '" in file "' // TRIM( GNAME1 ) // '"'
                        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                    END IF
                ELSE
                    K = M - NVARSA
                    MESG = 'Reading ' // LUNAMB( K )
                    CALL M3MESG( MESG )

                    IF( .NOT. READ3( GNAME2, LUNAMB( K ), 1, 0, 0, 
     &                               LUSE( 1, 1, N ) ) ) THEN
                        MESG = 'Could not find "' // LUNAMB( K ) // 
     &                         '" in file "' // TRIM( GNAME2 ) // '"'
                        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                    END IF
                END IF
            END IF
        END DO

C.........  Allocate memory and read forest inventory data
        ALLOCATE( FIA( NCOLS, NROWS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'FIA', PROGNAME )

        IF ( .NOT. READ3( GNAMET, 'FOREST', 1, 0, 0, FIA ) ) THEN
            MESG = 'Could not read FOREST from file "' //
     &              TRIM( GNAMET ) // '"'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

C.........  Allocate memory for output normalized fluxes
        ALLOCATE( AVGEMIS( NCOLS, NROWS, NSEF, NSEASONS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'AVGEMIS', PROGNAME )

        ALLOCATE( NOEMIS( NCOLS, NROWS, NNO ), STAT=IOS )
        CALL CHECKMEM( IOS, ' NOEMIS', PROGNAME )

        ALLOCATE( AVGLAI( NCOLS, NROWS, NLAI, NSEASONS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'AVGLAI', PROGNAME )

        ALLOCATE( SUMEM( NSEF ), STAT=IOS )
        CALL CHECKMEM( IOS, 'SUMEM', PROGNAME )
       
        ALLOCATE( NOEM( NNO ), STAT=IOS )
        CALL CHECKMEM( IOS, 'NOEM', PROGNAME )

        ALLOCATE( SUMEMW( NSEF ), STAT=IOS )
        CALL CHECKMEM( IOS, 'SUMEMW', PROGNAME )

        ALLOCATE( SUMLAI( NLAI ), STAT=IOS )
        CALL CHECKMEM( IOS, 'SUMLAI', PROGNAME )

        ALLOCATE( SUMLAIW( NLAI ), STAT=IOS )
        CALL CHECKMEM( IOS, 'SUMLAIW', PROGNAME )

        AVGEMIS = 0.0  !  array
        AVGLAI  = 0.0  !  array
        NOEMIS  = 0.0  !  array

C.........  Calculate normalized fluxes 
        DO I = 1, NCOLS
            DO J = 1, NROWS

C.................  Initialize variables
                SUMEM = 0.0   ! array
                SUMEMW = 0.0  ! array
                NOEM  = 0.0   ! array
                SUMLAI = 0.0  ! array
                SUMLAIW = 0.0 ! array

C.................  Sum USGS forest
                TOTFOR = LUSE( I, J, IUSDEC ) + LUSE( I, J, IUSBRD ) +
     &                   LUSE( I, J, IUSCON ) + LUSE( I, J, IUSMIX )
         
                DO M = 1, NVEG

C.....................  Assuming that land use is in percentages
C                       Converting to units of area (km**2)
                    IF ( VFLAG ) THEN
                        VEGAREA = LUSE( I, J, M ) * PRCNT2KM2( I, J )
                    ELSE
                        VEGAREA = LUSE( I, J, M ) * PRCNT2KM2( 1, 1 )
                    END IF

C.....................  If FIA and USGS data greater than 0, then use shrubland factors
                    USE_SHRUB = .FALSE.

                    IF ( TOTFOR > 0. .AND. FIA( I, J ) > 0 ) THEN

                        IF ( M .EQ. IUSDEC .OR. M .EQ. IUSBRD .OR.
     &                       M .EQ. IUSCON .OR. M .EQ. IUSMIX ) THEN

                            USE_SHRUB = .TRUE.
                  
                        END IF

                    END IF

                    DO N = 1, NSEF
                    
                        BTMP = BIOTYPES( N ) 

C.........................  Special handling for NO emissions
                        IF ( TRIM( BTMP ) == 'NO' ) THEN
                        
                            IF ( VEGAREA > 0. ) THEN
                                
                                IF( IS_AG( M,IUSCGS,IUSCWD,IUSCDY,
     &                                       IUSCIR,IALFAL,IBARLE,
     &                                       ICORN ,ICOTTO,IGRASS,
     &                                       IHAY  ,IMISCC,IOATS ,
     &                                       IPASTU,IPEANU,IPOTAT,
     &                                       IRICE ,IRYE  ,ISORGH,
     &                                       ISOYBE,ITOBAC,IWHEAT)) THEN

C.....................................  Compute NO emissions for agriculture regions 
C                                      during growing season
                                    IF( IS_HAG (M,IUSCGS,IUSCWD) ) THEN
                                        NOEM( 1 ) = NOEM( 1 ) 
     &                                      + VEGAREA * EMFAC(M,N)*0.5
 
                                        NOEM( 3 ) = NOEM( 3 ) 
     &                                      + VEGAREA * EMFAC(M,N)*0.5
                                    ELSE
                                        NOEM( 1 ) = NOEM( 1 )
     &                                      + VEGAREA * EMFAC(M,N)
                                    END IF

C.....................................  Compute NO emissions for agriculture regions 
C                                       outside growing season
                                    NOEM( 2 ) = NOEM( 2 )
     &                                  + VEGAREA * EMFAC(IGRASS,N) 

                                ELSE
                                
C.....................................  Compute NO emissions for Non-Agriculture regions 
                                    NOEM( 3 ) = NOEM( 3 ) 
     &                                  + VEGAREA * EMFAC(M,N)

                                END IF

                            END IF

                        ELSE

                            IF ( .NOT. USE_SHRUB ) THEN
                                EFTMP = EMFAC( M, N )
                            ELSE
                                EFTMP = EMFAC( IUSSHR, N )
                            ENDIF

                            IF ( VEGAREA > 0. ) THEN

C.................................  Compute summer emissions
                                SUMEM( N ) = SUMEM( N ) 
     &                               + VEGAREA * EFTMP

C.................................  Compute winter emissions
                                SUMEMW( N ) = SUMEMW( N )
     &                               + VEGAREA * EFTMP * WFAC( M )
                            END IF

C.............................  Compute LAI on ISOP and MBO and METH
C                               Note: assumption that these are the first three species
C                               in the BIOTYPES array
                            IF ( N .LE. NLAI ) THEN
                                SUMLAI( N ) = SUMLAI( N ) + VEGAREA
     &                                * LAI( M ) * EFTMP
                                SUMLAIW( N ) = SUMLAIW( N ) + VEGAREA
     &                                * LAI( M ) * EFTMP * WFAC( M )
                            END IF
                            
                        END IF  ! check if NO emissions
                    END DO  ! end of emis fac loop
                END DO  ! end of veg land use loop

                DO K = 1, NLAI

                    IF ( SUMLAI( K ) <= 1E-06 ) THEN
                        AVGLAI( I, J, K, 1 ) = 0.
                    ELSE IF ( SUMEM( K ) == 0. ) THEN
                        AVGLAI(  I, J, K, 1 ) = 0.
                    ELSE
                        AVGLAI( I, J, K, 1 ) = SUMLAI( K )/SUMEM( K ) 
                    ENDIF
                    
                    IF ( SUMLAIW( K ) <= 1E-06 ) THEN
                        AVGLAI( I, J, K, 2 ) = 0.
                    ELSE IF ( SUMEMW( K ) == 0. ) THEN
                        AVGLAI(  I, J, K, 2 ) = 0.
                    ELSE
                        AVGLAI( I, J, K, 2) = SUMLAIW( K ) /SUMEMW( K )
                    ENDIF

                END DO

                DO N = 1, NSEF

                    BTMP = BIOTYPES( N ) 

C.....................  Check for NO emissions
                    IF ( TRIM( BTMP ) == 'NO' ) THEN    
                        NOEMIS( I, J, 1 ) = NOEM( 1 ) 
                        NOEMIS( I, J, 2 ) = NOEM( 2 )
                        NOEMIS( I, J, 3 ) = NOEM( 3 )
                    ELSE
                        AVGEMIS( I, J, N, 1 ) = SUMEM( N ) 
                        AVGEMIS( I, J, N, 2 ) = SUMEMW( N )
                    END IF

                END DO  ! end loop over emission factors
            END DO  ! end loop over rows
        END DO  ! end loop over columns

C.........  Write output file
        I = 0
        DO M = 1, NSEASONS
            DO B = 1, NSEF
                BTMP = BIOTYPES( B ) 

C.................  Handle types other than NO
                IF ( TRIM( BTMP ) /= 'NO' ) THEN  

                    I = I + 1 
                    IF ( .NOT. WRITE3( ENAME, VNAME3D( I ), 0, 0,
     &                                 AVGEMIS(1, 1, B, M ) ) ) THEN
                        MESG = 'Could not write "' //
     &                          TRIM( VNAME3D( I ) ) //
     &                          '" to "' // TRIM( ENAME ) // '"'
                        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                    END IF
                END IF

            END DO

            DO N = 1, NLAI

                I = I + 1
                IF ( .NOT. WRITE3( ENAME, VNAME3D( I ), 0, 0,
     &                             AVGLAI(1, 1, N, M ) ) ) THEN
                    MESG = 'Could not write "' //
     &                      TRIM( VNAME3D( I ) ) //
     &                      '" to "' // TRIM( ENAME ) // '"'
                    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                END IF

            END DO

        END DO

        I = I + 1 
        IF ( .NOT. WRITE3( ENAME, VNAME3D( I ), 0, 0,
     &                     NOEMIS(1, 1, 1 ) ) ) THEN
            MESG = 'Could not write "' //
     &              TRIM( VNAME3D( I ) ) //
     &              '" to "' // TRIM( ENAME ) // '"'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

        I = I + 1 
        IF ( .NOT. WRITE3( ENAME, VNAME3D( I ), 0, 0,
     &                     NOEMIS(1, 1, 2 ) ) ) THEN
            MESG = 'Could not write "' //
     &              TRIM( VNAME3D( I ) ) //
     &              '" to "' // TRIM( ENAME ) // '"'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

        I = I + 1 
        IF ( .NOT. WRITE3( ENAME, VNAME3D( I ), 0, 0,
     &                     NOEMIS(1, 1, 3 ) ) ) THEN
            MESG = 'Could not write "' //
     &              TRIM( VNAME3D( I ) ) //
     &              '" to "' // TRIM( ENAME ) // '"'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

C.........  End of subroutine
        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10 ( A, :, I5, :, 2X ) )

C***************** CONTAINS ********************************************

        CONTAINS

C.............  This internal function checks for "half" agricultural areas
            LOGICAL FUNCTION IS_HAG( M,IUSCGS, IUSCWD )

            IMPLICIT NONE
            
C.............  Function arguments            
            INTEGER, INTENT(IN) :: M
            INTEGER, INTENT(IN) :: IUSCGS
            INTEGER, INTENT(IN) :: IUSCWD

C-----------------------------------------------------------------------------

            IS_HAG = .FALSE.

            IF( M == IUSCGS ) IS_HAG = .TRUE.
            IF( M == IUSCWD ) IS_HAG = .TRUE.

            RETURN
            
            END FUNCTION IS_HAG

C-----------------------------------------------------------------------------

C.............  This internal function checks for agricultural areas
            LOGICAL FUNCTION IS_AG( M, IUSCGS, IUSCWD, IUSCDY,
     &                                 IUSCIR, IALFAL, IBARLE,
     &                                 ICORN , ICOTTO, IGRASS,
     &                                 IHAY  , IMISCC, IOATS ,
     &                                 IPASTU, IPEANU, IPOTAT,
     &                                 IRICE , IRYE  , ISORGH,
     &                                 ISOYBE, ITOBAC, IWHEAT )
     
            IMPLICIT NONE
            
C.............  Function arguments  
            INTEGER, INTENT(IN) :: M
            INTEGER, INTENT(IN) :: IUSCGS
            INTEGER, INTENT(IN) :: IUSCWD
            INTEGER, INTENT(IN) :: IUSCDY
            INTEGER, INTENT(IN) :: IUSCIR
            INTEGER, INTENT(IN) :: IALFAL
            INTEGER, INTENT(IN) :: IBARLE
            INTEGER, INTENT(IN) :: ICORN
            INTEGER, INTENT(IN) :: ICOTTO
            INTEGER, INTENT(IN) :: IGRASS
            INTEGER, INTENT(IN) :: IHAY
            INTEGER, INTENT(IN) :: IMISCC
            INTEGER, INTENT(IN) :: IOATS
            INTEGER, INTENT(IN) :: IPASTU
            INTEGER, INTENT(IN) :: IPEANU
            INTEGER, INTENT(IN) :: IPOTAT
            INTEGER, INTENT(IN) :: IRICE
            INTEGER, INTENT(IN) :: IRYE
            INTEGER, INTENT(IN) :: ISORGH
            INTEGER, INTENT(IN) :: ISOYBE
            INTEGER, INTENT(IN) :: ITOBAC
            INTEGER, INTENT(IN) :: IWHEAT

C-----------------------------------------------------------------------------

            IS_AG = .FALSE.

            IF( M == IUSCGS ) IS_AG = .TRUE.
            IF( M == IUSCWD ) IS_AG = .TRUE.
            IF( M == IUSCDY ) IS_AG = .TRUE.
            IF( M == IUSCIR ) IS_AG = .TRUE.
            IF( M == IALFAL ) IS_AG = .TRUE.
            IF( M == IBARLE ) IS_AG = .TRUE.
            IF( M == ICORN  ) IS_AG = .TRUE.
            IF( M == ICOTTO ) IS_AG = .TRUE.
            IF( M == IGRASS ) IS_AG = .TRUE.
            IF( M == IHAY   ) IS_AG = .TRUE.
            IF( M == IMISCC ) IS_AG = .TRUE.
            IF( M == IOATS  ) IS_AG = .TRUE.
            IF( M == IPASTU ) IS_AG = .TRUE.
            IF( M == IPEANU ) IS_AG = .TRUE.
            IF( M == IPOTAT ) IS_AG = .TRUE.
            IF( M == IRICE  ) IS_AG = .TRUE.
            IF( M == IRYE   ) IS_AG = .TRUE.
            IF( M == ISORGH ) IS_AG = .TRUE.
            IF( M == ISOYBE ) IS_AG = .TRUE.
            IF( M == ITOBAC ) IS_AG = .TRUE.
            IF( M == IWHEAT ) IS_AG = .TRUE.

            RETURN
            
            END FUNCTION IS_AG

C-----------------------------------------------------------------------------
        
        END SUBROUTINE NORMBEIS312
