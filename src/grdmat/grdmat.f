
        PROGRAM GRDMAT

C***********************************************************************
C  program body starts at line 148
C
C  DESCRIPTION:
C     Creates the gridding matrix for any source category and creates the
C     "ungridding" matrix for mobile sources.
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
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
C***************************************************************************

C...........   MODULES for public variables
C...........   This module is the source inventory arrays
        USE MODSOURC, ONLY: XLOCA, YLOCA, XLOC1, YLOC1, XLOC2, YLOC2,
     &                      CELLID   

C...........   This module contains the cross-reference tables
        USE MODXREF, ONLY: SRGIDPOS, SGFIPPOS

C.........  This module contains the global variables for the 3-d grid
        USE MODGRID, ONLY: GRDNM, NCOLS, NROWS, COORD, GDTYP, 
     &                     P_ALP, P_BET, P_GAM, XCENT, YCENT, NGRID,
     &                     OFFLAG

C.........  This module contains the information about the source category
        USE MODINFO, ONLY: CATEGORY, CRL, NSRC, NIPPA, EANAM

C.........  This module is required by the FileSetAPI
        USE MODFILESET

        IMPLICIT NONE

C...........   INCLUDES:
        
        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
        INCLUDE 'SETDECL.EXT'   !  FileSetAPI variables
c      INCLUDE 'PARMS3.EXT'    !  I/O API parameters (in modfileset)
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations 
c      INCLUDE 'FDESC3.EXT'    !  I/O API file description data structures (in modfileset)

C...........   EXTERNAL FUNCTIONS and their descriptions:
        
        CHARACTER(2)       CRLF
        LOGICAL            DSCM3GRD
        LOGICAL            ENVYN
        CHARACTER(IODLEN3) GETCFDSC
        INTEGER            INDEX1
        LOGICAL            INGRID
        INTEGER            PROMPTFFILE
        CHARACTER(16)      VERCHAR
   
        EXTERNAL  CRLF, ENVYN, DSCM3GRD, GETCFDSC, INDEX1, INGRID, 
     &            PROMPTFFILE, PROMPTMFILE, VERCHAR

C...........   LOCAL PARAMETERS
        CHARACTER(50), PARAMETER :: CVSW = '$Name$' ! CVS release tag

C...........   LOCAL VARIABLES and their descriptions:

C...........   Gridding Matrix

        INTEGER, ALLOCATABLE :: GMAT( : ) ! Contiguous gridding matrix

C...........   Ungridding Matrix

        INTEGER, ALLOCATABLE :: UMAT( : ) ! Contiguous ungridding matrix

C.........  Array that contains the names of the inventory variables needed 
C           for this program
        CHARACTER(IOVLEN3) IVARNAMS( MXINVARR )

C...........   File units and logical/physical names
c        INTEGER         ADEV    !  for adjustments file
        INTEGER         IDEV    ! tmp unit number if ENAME is map file
        INTEGER         LDEV    !  log-device
        INTEGER         KDEV    !  for link defs file
        INTEGER         GDEV    !  for surrogate coeff file
        INTEGER      :: MDEV = 0!  mobile sources codes file
        INTEGER      :: RDEV = 0!  report file for surrogates usage
        INTEGER      :: SDEV = 0!  ASCII part of inventory unit no.
        INTEGER         XDEV    !  for surrogate xref  file

        CHARACTER(16)   ANAME   !  logical name for ASCII inventory input file
        CHARACTER(16)   ENAME   !  logical name for i/o api inventory input file
        CHARACTER(16)   INAME   !  tmp name for inven file of unknown fmt
        CHARACTER(16)   GNAME   !  logical name for grid matrix output file
        CHARACTER(16)   UNAME   !  logical name for ungrid matrix output file

C...........   Other local variables
        
        INTEGER         I, L1, L2, K, S !  indices and counters.

        INTEGER         COL     ! tmp column
        INTEGER         CMAX    ! max number srcs per cell
        INTEGER         CMIN    ! min number srcs per cell
        INTEGER         ENLEN   ! length of the emissions inven name
        INTEGER         IOS     ! i/o status
        INTEGER         NK      ! Number of gridding coefficients 
        INTEGER         NKU     ! Number of ungridding coefficients
        INTEGER         NINVARR ! no. of inventory characteristics
        INTEGER         NMATX   ! no cell-source intersections
        INTEGER         NMATXU  ! no county-source intrsctns for all sources
        INTEGER         MXCCL   ! max no cells per county or link
        INTEGER         MXCSRC  ! max no cells per source
        INTEGER         MXSCEL  ! max no sources per cell
        INTEGER         ROW     ! tmp row
        INTEGER         SRGNROWS! no. rows in surrogates file
        INTEGER         SRGNCOLS! no. cols in surrogates file

        REAL            CAVG   ! average number sources per cell

        LOGICAL      :: A2PFLAG = .FALSE.  ! true: inv has ar-to-pt locations
        LOGICAL      :: AFLAG   = .FALSE.  ! true: use grid adjustments file
        LOGICAL      :: DFLAG   = .FALSE.  ! true: use link defs file
        LOGICAL      :: EFLAG   = .FALSE.  ! true: error found
        LOGICAL      :: SRGFLAG = .FALSE.  ! true: surrogates are needed
        LOGICAL      :: UFLAG   = .FALSE.  ! true: create ungridding matrix
        LOGICAL      :: VFLAG   = .FALSE.  ! true: use variable grid

        CHARACTER(16)       COORUNIT !  coordinate system projection units
        CHARACTER(16)    :: INVGRDNM  = ' '  !  inventory grid name
        CHARACTER(16)    :: SRGGRDNM  = ' '  !  surrogates file grid name
        CHARACTER(16)       SRGFMT   !  surrogates format
        CHARACTER(80)       GDESC    !  grid description
        CHARACTER(300)      MESG     !  message buffer

        CHARACTER(IODLEN3)  IFDESC2, IFDESC3 !  fields 2 & 3 from PNTS FDESC

        CHARACTER(16) :: PROGNAME = 'GRDMAT'   !  program name

C***********************************************************************
C   begin body of program GRDMAT
        
        LDEV = INIT3()

C.........  Write out copywrite, version, web address, header info, and prompt
C           to continue running the program.
        CALL INITEM( LDEV, CVSW, PROGNAME )

C.........  Get environment variables that control this program
        AFLAG = ENVYN( 'GRDMAT_ADJUST',
     &                 'Use grid adjustments file or not', 
     &                 .FALSE., IOS )

        DFLAG = ENVYN( 'GRDMAT_LINKDEFS',
     &                 'Use link definitions file or not', 
     &                 .FALSE., IOS )

        VFLAG = ENVYN( 'USE_VARIABLE_GRID',
     &                 'Use variable grid definition',
     &                 .FALSE., IOS )

C.........  Temporary section for disallowing optional files
        IF( AFLAG ) THEN
            MESG = 'NOTE: Grid adjustments file is not supported yet!'
            CALL M3MSG2( MESG )
            AFLAG = .FALSE.
        END IF

        IF( DFLAG ) THEN
            MESG = 'NOTE: Link definitions file is not supported yet!'
            CALL M3MSG2( MESG )
            DFLAG = .FALSE.
        END IF

C.........  Set source category based on environment variable setting
        CALL GETCTGRY

C.........  Get inventory file names given source category
        CALL GETINAME( CATEGORY, ENAME, ANAME )

C.........  Prompt for and open inventory file 
        INAME = ENAME
        MESG = 'Enter logical name for the MAP INVENTORY file'
        IDEV = PROMPTFFILE( MESG, .TRUE., .TRUE., INAME, PROGNAME )

C.........  Open and read map file
        CALL RDINVMAP( INAME, IDEV, ENAME, ANAME, SDEV )

c        IF( AFLAG ) 
c     &  ADEV = PROMPTFFILE( 
c     &           'Enter logical name for ADJUSTMENT FACTORS file',
c     &           .TRUE., .TRUE., CRLF // 'ADJUST', PROGNAME )

C.............  Store source specific information based on header
        CALL GETSINFO( ENAME )

        IFDESC2  = GETCFDSC( FDESC3D, '/FROM/', .TRUE. )
        IFDESC3  = GETCFDSC( FDESC3D, '/VERSION/', .TRUE. )
        INVGRDNM = GETCFDSC( FDESC3D, '/GRIDNAME/', .FALSE. )

C.........  Set inventory variables to read for all source categories
        IVARNAMS( 1 ) = 'IFIP'

        SELECT CASE ( CATEGORY )

        CASE ( 'AREA' )
            NINVARR = 4
            IVARNAMS( 2 ) = 'CSCC'
            IVARNAMS( 3 ) = 'CELLID'
            IVARNAMS( 4 ) = 'CSOURC'

C............  Check to see if point locations are in the AREA
C              file (i.e. we have area-to-point sources)
            K = INDEX1( 'XLOCA', NVARSET, VNAMESET )
            IF ( K .GT. 0 ) THEN
            
C.................  Make sure we're not using a variable grid
                IF( VFLAG ) THEN
                    MESG = 'Cannot use area-to-point sources ' //
     &                     'with a variable grid.'
                    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                END IF
                
                A2PFLAG = .TRUE.
                NINVARR = 6
                IVARNAMS( 5 ) = 'XLOCA'
                IVARNAMS( 6 ) = 'YLOCA'
            END IF

        CASE ( 'MOBILE' )
            NINVARR = 10
            IVARNAMS( 2 ) = 'IRCLAS'
            IVARNAMS( 3 ) = 'IVTYPE'
            IVARNAMS( 4 ) = 'CSOURC'
            IVARNAMS( 5 ) = 'CSCC'
            IVARNAMS( 6 ) = 'CLINK'
            IVARNAMS( 7 ) = 'XLOC1'
            IVARNAMS( 8 ) = 'YLOC1'
            IVARNAMS( 9 ) = 'XLOC2'
            IVARNAMS( 10 ) = 'YLOC2'

        CASE ( 'POINT' )
            NINVARR = 3
            IVARNAMS( 2 ) = 'XLOCA'
            IVARNAMS( 3 ) = 'YLOCA'

        END SELECT

C.........  Allocate memory for and read in required inventory characteristics
        CALL RDINVCHR( CATEGORY, ENAME, SDEV, NSRC, NINVARR, IVARNAMS )

C.........  Define source-category-specific settings
        SELECT CASE( CATEGORY )

        CASE( 'AREA' )

C.............  Determine if surrogates are needed by checking whether all cell
C               values are defined or not
            DO S = 1, NSRC
                IF( CELLID( S ) .LE. 0 ) SRGFLAG = .TRUE.
                IF( SRGFLAG ) EXIT
            END DO

        CASE( 'MOBILE' )
        
C.............  Check for link-based sources when using a variable grid        
            IF( MAXVAL( XLOC1 ) .GT. AMISS3 .AND. VFLAG ) THEN
                MESG = 'Cannot use link-based data ' //
     &                 'with a variable grid.'
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF

C.............  If all sources are link-based, don't need surrogates (needs to be tested)
c            IF( MINVAL( XLOC1 ) .GT. AMISS3 ) THEN
c                SRGFLAG = .FALSE.
c            ELSE
               SRGFLAG = .TRUE.
c            END IF

        CASE( 'POINT' )

        END SELECT

C.........  Get gridding surrogates and x-ref file, if needed
        IF( SRGFLAG ) THEN

            XDEV = PROMPTFFILE( 
     &           'Enter logical name for GRIDDING SURROGATE XREF file',
     &           .TRUE., .TRUE., CRL // 'GREF', PROGNAME )

            GDEV = PROMPTFFILE( 
     &           'Enter logical name for SURROGATE COEFFICIENTS file',
     &           .TRUE., .TRUE., CRL // 'GPRO', PROGNAME )
        END IF

C.........  Get mobile-specific files
        IF( CATEGORY .EQ. 'MOBILE' ) THEN
            MESG = 'Enter logical name for MOBILE CODES file'
            MDEV = PROMPTFFILE( MESG, .TRUE., .TRUE., 'MCODES',
     &                          PROGNAME )

            IF( DFLAG ) KDEV = PROMPTFFILE( 
     &               'Enter logical name for LINK DEFINITIONS file',
     &               .TRUE., .TRUE., CRL // 'GLNK', PROGNAME )

        END IF  ! End of mobile file opening

C.........  If surrogates are in use, read the gridding cross-reference,
C           and the header of the surrogates file to get the grid on which
C           the surrogates are available.
        IF( SRGFLAG ) THEN

C.............  Build unique lists of SCCs and country/state/county codes
C               from the inventory arrays
            CALL GENUSLST

C.............  For mobile sources, read the mobile codes
            IF( MDEV .GT. 0 ) CALL RDMVINFO( MDEV )

            CALL M3MSG2( 'Reading gridding cross-reference file...' )

C.............  Read the gridding cross-reference
            CALL RDGREF( XDEV )

            CALL M3MSG2( 'Reading gridding surrogates header...' )

C.............  Read the surrogates header and initialize the grid description
C.............  Also, obtain the format of the file.
C.............  Save the name of the input grid
            CALL RDSRGHDR( VFLAG, GDEV, SRGFMT )
            SRGGRDNM = GRDNM
            SRGNCOLS = NCOLS
            SRGNROWS = NROWS

        END IF   ! If surrogates are being used or not

C.........  Get grid name from the environment and read grid parameters
        IF ( .NOT. DSCM3GRD( GDNAM3D, GDESC, COORD, GDTYP3D, COORUNIT,
     &                     P_ALP3D, P_BET3D, P_GAM3D, XCENT3D, 
     &                     YCENT3D, XORIG3D, YORIG3D, XCELL3D,
     &                     YCELL3D, NCOLS3D, NROWS3D, NTHIK3D)) THEN

            MESG = 'Could not get Models-3 grid description.'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

        END IF

C.........  Check or initialize the output grid grid settings (depends on
C           if a surrogates file is being used). For variable grids,
C           do not allow subgrids.
        IF( VFLAG ) THEN
            CALL CHKGRID( GDNAM3D, 'GRIDDESC', 0, EFLAG )
            
            IF( EFLAG ) THEN
                MESG = 'Problem with variable grid input data.'
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 ) 
            END IF
        ELSE
            CALL CHKGRID( GDNAM3D, 'GRIDDESC', 1, EFLAG )
            
            IF ( EFLAG ) THEN
                MESG = 'Problem with gridded input data.'
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF
        END IF

C.........  If surrogates are needed, read the gridding surrogates,  
C           allocate memory for the surrogate assignments, and assign
C           surrogates to each source.
        IF( SRGFLAG ) THEN

            CALL M3MSG2( 'Reading gridding surrogates...' )

C.............  Allocate memory for and read the gridding surrogates file,
C               extracting data for a subgrid, if necessary
            CALL RDSRG( VFLAG, GDEV, SRGFMT, SRGNROWS, SRGNCOLS )
            
C..............  Read the link definition file
c            CALL RDLNKDEF( )

C.............  Allocate memory for indices to surrogates tables for each source
            ALLOCATE( SRGIDPOS( NSRC ), STAT=IOS )
            CALL CHECKMEM( IOS, 'SRGIDPOS', PROGNAME )
            ALLOCATE( SGFIPPOS( NSRC ), STAT=IOS )
            CALL CHECKMEM( IOS, 'SGFIPPOS', PROGNAME )

C.............  Assigns the index of the surrogate to each source (stored
C               in SRGIDPOS passed through MODXREF)
            CALL ASGNSURG

        END IF   ! If surrogates are being used or not

C.........  Ensure that the output grid is consistent with the input grid
        IF ( INVGRDNM .NE. ' ' .AND.
     &       INVGRDNM .NE. GRDNM     ) THEN

            L1 = LEN_TRIM( GRDNM )
            L2 = LEN_TRIM( INVGRDNM )
            MESG = 'Output grid "' // GRDNM( 1:L1 ) // 
     &             '" is inconsistent with inventory input grid "'//
     &             INVGRDNM( 1:L2 ) // '".'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

        END IF

C.........  Write message stating grid name and description
        L1 = LEN_TRIM( GRDNM )
        MESG = 'NOTE: Output grid "' // GRDNM( 1:L1 ) // 
     &         '" set; described as' // CRLF() // BLANK10 // GDESC
        CALL M3MSG2( MESG )

C.........  If output grid is different from surrogates, write message
        IF ( OFFLAG ) THEN
            L1 = LEN_TRIM( SRGGRDNM )
            MESG = 'NOTE: gridding surrogates extracted for output '//
     &             'grid from grid "' // SRGGRDNM( 1:L1 ) // '"'
            CALL M3MSG2( MESG )
        END IF

C.........  If area or mobile inventory has point source locations,
C           convert point source coordinates from lat-lon to output grid
        IF( A2PFLAG ) THEN
            CALL CONVRTXY( NSRC, GDTYP, GRDNM, P_ALP, P_BET, P_GAM, 
     &                     XCENT, YCENT, XLOCA, YLOCA )
        END IF

C.........  Determine if ungridding matrix is needed
        I = INDEX1( 'VMT', NIPPA, EANAM )
        IF( I .GT. 0 ) THEN
        
C.............  Make sure we're not using a variable grid
            IF( VFLAG ) THEN
                MESG = 'Cannot create ungridding matrix ' //
     &                 'when using a variable grid.'
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF
            
            UFLAG = .TRUE.
            MESG = 'NOTE: VMT detected in inventory, ungridding matrix '
     &             //'will be created'
        ELSE
            MESG = 'NOTE: VMT not detected in inventory, ungridding ' //
     &             'matrix will NOT be created'
        END IF
        CALL M3MSG2( MESG )

C.........  Depending on source category, convert coordinates, determine size
C           of gridding matrix, and allocate gridding matrix.

        SELECT CASE( CATEGORY )

        CASE( 'AREA' )

C.............  Determine sizes for allocating area gridding matrix 
            CALL SIZGMAT( CATEGORY, NSRC, MXSCEL, MXCSRC, 
     &                    MXCCL, NMATX, NMATXU )

C.............  Allocate memory for mobile source gridding matrix
            ALLOCATE( GMAT( NGRID + 2*NMATX ), STAT=IOS )
            CALL CHECKMEM( IOS, 'GMAT', PROGNAME )

        CASE( 'MOBILE' )

C.............  Convert mobile source coordinates from lat-lon to output grid
            CALL CONVRTXY( NSRC, GDTYP, GRDNM, P_ALP, P_BET, P_GAM,
     &                     XCENT, YCENT, XLOC1, YLOC1 )
            CALL CONVRTXY( NSRC, GDTYP, GRDNM, P_ALP, P_BET, P_GAM, 
     &                     XCENT, YCENT, XLOC2, YLOC2 )

C.............  Determine sizes for allocating mobile gridding matrix 
            CALL SIZGMAT( CATEGORY, NSRC, MXSCEL, MXCSRC,  
     &                    MXCCL, NMATX, NMATXU )
 
C.............  Allocate memory for mobile source gridding matrix
            ALLOCATE( GMAT( NGRID + 2*NMATX ), STAT=IOS )
            CALL CHECKMEM( IOS, 'GMAT', PROGNAME )

C.............  Allocate memory for mobile source ungridding matrix
            IF( .NOT. UFLAG ) NMATXU = 1
            ALLOCATE( UMAT( NSRC + 2*NMATXU ), STAT=IOS )
            CALL CHECKMEM( IOS, 'UMAT', PROGNAME )

        CASE( 'POINT' )

C.............  Convert point source coordinates from lat-lon to output grid
            CALL CONVRTXY( NSRC, GDTYP, GRDNM, P_ALP, P_BET, P_GAM, 
     &                     XCENT, YCENT, XLOCA, YLOCA )

C.............  Set the number of source-cell intersections
            DO S = 1, NSRC
                IF( INGRID( XLOCA( S ), YLOCA( S ), 
     &                      NCOLS, NROWS, COL, ROW  ) ) THEN
                    NMATX = NMATX + 1
                END IF
            END DO

C.............  Allocate memory for point source gridding matrix
            ALLOCATE( GMAT( NGRID + NMATX ), STAT=IOS )
            CALL CHECKMEM( IOS, 'GMAT', PROGNAME )

        END SELECT

C.........  Abort of there are no source-cell intersections
        IF( NMATX .EQ. 0 ) THEN
            MESG = 'No source-cell intersections found.'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

C.........  Get file names; open output gridding matrix (and ungridding matrix
C           for mobile) using grid characteristics from DSCM3GRD() above        
C.........  Also open report file 
        CALL OPENGMAT( NMATX, NMATXU, UFLAG, IFDESC2, IFDESC3, VFLAG,
     &                 GNAME, UNAME, RDEV )

        CALL M3MSG2( 'Generating gridding matrix...' )

C.........  Generate gridding matrix for given source category, and write it
C           out.  It is necessary to write it out while in the subroutine,
C           because of the type transformation from real to integer that
C           is done so the sparse i/o api format can be used.

        SELECT CASE( CATEGORY )

        CASE( 'AREA' )

            CALL GENAGMAT( GNAME, RDEV, MXSCEL, NSRC, NMATX, 
     &                     GMAT( 1 ), GMAT( NGRID+1 ), 
     &                     GMAT( NGRID+NMATX+1 ), NK, CMAX, CMIN )

        CASE( 'MOBILE' )

            CALL GENMGMAT( ENAME, GNAME, UNAME, RDEV, MXSCEL, MXCSRC, 
     &                     MXCCL, NSRC, NMATX, NMATXU, UFLAG, GMAT( 1 ),
     &                     GMAT( NGRID+1 ), GMAT( NGRID+NMATX+1 ), 
     &                     UMAT( 1 ), UMAT( NSRC+1 ), 
     &                     UMAT( NSRC+NMATXU+1 ),NK, CMAX, CMIN, NKU )

        CASE( 'POINT' )
   
            CALL GENPGMAT( GNAME, NSRC, NGRID, XLOCA, YLOCA, VFLAG,
     &                     GMAT( 1 ), GMAT( NGRID+1 ), NK, CMAX, CMIN )

        END SELECT

C.........  Report statistics for gridding matrix

        CAVG = FLOAT( NK ) / FLOAT( NGRID )
        CALL M3MSG2( 'GRIDDING-MATRIX statistics:' )

        WRITE( MESG,94010 ) 
     &         'Total number of coefficients    :', NK   ,
     &         CRLF() // BLANK5 //
     &         'Max  number of sources per cell :', CMAX,
     &         CRLF() // BLANK5 //
     &         'Min  number sources per cell > 0:', CMIN

        L1 = LEN_TRIM( MESG )
        WRITE( MESG,94020 ) MESG( 1:L1 ) // CRLF() // BLANK5 //
     &         'Mean number of sources per cell :', CAVG

        CALL M3MSG2( MESG )

C.........  Report statistics for ungridding matrix
        IF( UFLAG ) THEN

            CAVG = FLOAT( NKU ) / FLOAT( NSRC )
            CALL M3MSG2( 'UNGRIDDING-MATRIX statistics:' )

            WRITE( MESG,94010 ) 
     &        'Total number of coefficients   :', NKU 

            WRITE( MESG, 94020 ) MESG( 1:LEN_TRIM( MESG ) ) //
     &        CRLF() // BLANK5 //
     &        'Mean number of cells per source:', CAVG

            CALL M3MSG2( MESG )

        END IF 

C.........  End of program
      
        CALL M3EXIT( PROGNAME, 0, 0, ' ', 0 )

C******************  FORMAT  STATEMENTS   ******************************

C...........   Informational (LOG) message formats... 92xxx

92000   FORMAT( 5X, A )

92010   FORMAT( 5X, A, :, I12 )

92020   FORMAT( 5X, A, :, F17.4 )


C...........   Formatted file I/O formats............ 93xxx

93000   FORMAT( A )

93010   FORMAT( A16 )


C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I10, :, 1X ) )

94020   FORMAT( A, :, F10.2 )

        END PROGRAM GRDMAT


