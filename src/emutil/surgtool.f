
        PROGRAM SURGTOOL

C***********************************************************************
C  program body starts at line 
C
C  DESCRIPTION: 
C      Inputs SMOKE-formatted gridding surrogates
C      for a "fine" input grid and a grid definition for a "coarse" output
C      grid.  It produces approximate "coarse" grid surrogates file, the 
C      accuracy of which depends on how fine the resolution of the input grid 
C      is relative to that of the output grid.
C
C  PRECONDITIONS REQUIRED:
C      - Input & output surrogates on Lat-Lon, Lambert Conformal, or UTM map
C        projections (can perform Lambert-to-Lambert and UTM zone-to-zone
C        transformations).
C      - Input grid resolution is much finer than output grid resolution
C      - Supports SMOKE-formatted surrogate coefficient files only.
C      - Input logical file names are defined.
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C             10/00 : Prototype by JMV
C              1/02 : Updated for Models-3 surrogate format and included
C                     in SMOKE emutil module
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
C*************************************************************************

C...........   MODULES for public variables
C...........   This module contains the gridding surrogates tables
        USE MODSURG, ONLY: NSRGREC, IDXSRGA, SCELLA, SFIPSA, SSRGIDA,
     &                     SFRACA

C.........  This module contains the global variables for the 3-d grid
        USE MODGRID, ONLY: GDTYP, P_ALP, P_BET, P_GAM, YCENT,
     &                     NCOLS, NROWS, XORIG, YORIG, XCELL, YCELL,
     &                     GRDNM

       IMPLICIT NONE

C...........   INCLUDES:
        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
        INCLUDE 'FDESC3.EXT'    !  I/O API file description data structures.

C...........   EXTERNAL FUNCTIONS and their descriptions:

        CHARACTER(2)    CRLF
        LOGICAL         DSCM3GRD
        INTEGER         PROMPTFFILE
        INTEGER         STR2INT
        
        EXTERNAL   CRLF, DSCM3GRD, PROMPTFFILE, STR2INT

C...........   LOCAL PARAMETERS
        CHARACTER(50), PARAMETER :: 
     &  CVSW = '$Name SMOKEv4.8_Oct2020$' ! CVS release tag

C............   Allocatable arrays for computing grid cell intersections
        INTEGER, ALLOCATABLE::   C2( :, : )     !  output grid cell
        REAL,    ALLOCATABLE::   XX( :, : )     !  output grid X
        REAL,    ALLOCATABLE::   YY( :, : )     !  output grid Y

C...........   Allocatable arrays for unsorted output surrogates
        INTEGER                 OUTNSRG       !  number of output entries
        INTEGER, ALLOCATABLE :: INDXA( : )    !  sorting index
        INTEGER, ALLOCATABLE :: CELLA( : )    !  cell number
        INTEGER, ALLOCATABLE :: FIPSA( : )    !  region code
        INTEGER, ALLOCATABLE :: SSCSA( : )    !  surrogate ID
        REAL   , ALLOCATABLE :: FRACA( : )    !  surrogate fraction

C...........   Logical file names and unit numbers

        INTEGER         LDEV    !  log-device
        INTEGER         SDEV    !  for  input surrogate coeff file
        INTEGER         XDEV    !  for output surrogate coeff  file

C...........   Arguments for GTPZ0:
        REAL(8)           DSCR            !  scratch variables
        INTEGER           DEG, MNT        !  scratch variables
        REAL(8)           CRDIN( 2 )      !  input coordinates x,y
        INTEGER(4)        INSYS           !  input projection code
        INTEGER(4)        INZONE          !  input utm zone, etc.
        REAL(8)           TPAIN( 15 )     !  input projection parameters
        INTEGER(4)        INUNIT          !  input units code
        INTEGER(4)        INSPH           !  spheroid code
        INTEGER(4)        IPR             !  error print flag
        INTEGER(4)        JPR             !  projection parameter print flag
        INTEGER(4)        LEMSG           !  error message unit number
        INTEGER(4)        LPARM           !  projection parameter unit number
        REAL(8)           CRDIO( 2 )      !  output coordinates x,y
        INTEGER(4)        IOSYS           !  output projection code
        INTEGER(4)        IOZONE          !  output utm zone, etc.
        REAL(8)           TPOUT( 15 )     !  output projection parameters
        INTEGER(4)        IOUNIT          !  output units code
        INTEGER(4)        IOSPH           !  spheroid code
        INTEGER(4)        LN27            !  NAD1927 file unit number
        INTEGER(4)        LN83            !  NAD1983 file unit number
        CHARACTER(128)    FN27            !  NAD1927 file name
        CHARACTER(128)    FN83            !  NAD1983 file name
        INTEGER(4)        LENGTH          !  NAD* record-length
        INTEGER(4)        IFLG            !  error flag

C...........    Output grid variables
        INTEGER      NCOLS2      ! number of grid columns
        INTEGER      NROWS2      ! number of grid rows
        INTEGER      NTHIK2      ! BOUNDARY:  perim thickness (cells)
        INTEGER      GDTYP2      ! grid type:  1=LAT-LON, 2=UTM, ...
        REAL(8)      P_ALP2      ! first, second, third map
        REAL(8)      P_BET2      ! projection descriptive
        REAL(8)      P_GAM2      ! parameters.
        REAL(8)      XCENT2      ! lon for coord-system X=0
        REAL(8)      YCENT2      ! lat for coord-system Y=0
        REAL(8)      XORIG2      ! X-coordinate origin of grid (map units)
        REAL(8)      YORIG2      ! Y-coordinate origin of grid
        REAL(8)      XCELL2      ! X-coordinate cell dimension
        REAL(8)      YCELL2      ! Y-coordinate cell dimension
        CHARACTER(IOVLEN3) :: GDNAM2 = ' '    ! output grid name
        CHARACTER(IOVLEN3) :: COORUNIT2 = ' ' !  coord sys projection units
        CHARACTER(IOVLEN3) :: COORD2 = ' '    ! coord system name
        CHARACTER(IODLEN3) :: GDESC2 = ' '    ! output grid desc (if any)

C...........   Other local variables
        INTEGER          C, S, F, J, K, L, M, N, R    !  counters, subscripts
        INTEGER          L1, L2                 !  tmp lengths
        INTEGER          COL, ROW               !  tmp row and column
        INTEGER          IOS                    !  I/O status
        INTEGER          IREC                   !  input line (record) number
        INTEGER          LCEL                   !  previous cell number
        INTEGER          LFIP                   !  previous region code
        INTEGER          LSSC                   !  previous surrogate code

        REAL             DDX, DDY               !  inverse of cell widths
        REAL             FRAC                   !  tmp output surrogate fraction
        REAL             X0, Y0                 !  surrogate file header
        REAL             XZ, YZ

        LOGICAL       :: EFLAG = .FALSE.        !  true: error found

        CHARACTER(16)    OUTTYPE
        CHARACTER(16)    OUTUNIT
        CHARACTER(16)    SRGFMT                 !  surrogates format
        CHARACTER(256)   MESG                   !  message buffer

        CHARACTER(16) :: PROGNAME = 'SURGTOOL'  !  program name

C***********************************************************************
C   begin body of program SURGTOOL

        LDEV = INIT3()

C.........  Write out copyright, version, web address, header info, and prompt
C           to continue running the program.
        CALL INITEM( LDEV, CVSW, PROGNAME )

C.........  Open input surrogates file
        SDEV = PROMPTFFILE( 
     &         'Enter logical name for input GRIDDING SURROGATE file',
     &         .TRUE., .TRUE., 'INFILE', PROGNAME )

C.........  Read the surrogates header and initialize the grid description
C.........  Also get the format of the surrogates file.
        CALL RDSRGHDR( .FALSE., SDEV, SRGFMT )

        L = LEN_TRIM( SRGFMT )
        MESG = 'NOTE: Input surrogates are ' // SRGFMT( 1:L ) // 
     &         ' format.'
        CALL M3MSG2( MESG )

C.........  Initialize gridding module
        CALL CHKGRID( GDNAM3D, 'SURROGATES', 1, EFLAG )

C.........  Precompute intermediate grid overlay variables
        SELECT CASE( GDTYP )
        CASE( LAMGRD3 )

           TPAIN( 1 ) = 0.0D0
           TPAIN( 2 ) = 0.0D0

           DSCR = P_ALP
           DEG  = INT( DSCR )                              !  int degrees
           DSCR = 60.0D0 * ( DSCR - DBLE( DEG ) )          !  minutes
           MNT  = INT( DSCR )                              !  int minutes
           DSCR = 60.0D0 * ( DSCR - DBLE( MNT ) )          !  seconds
           TPAIN( 3 ) = DSCR + 1000.0D0 * ( MNT + 1000*DEG ) !  dddmmmsss.sssD0

           DSCR = P_BET
           DEG  = INT( DSCR )                              !  int degrees
           DSCR = 60.0D0 * ( DSCR - DBLE( DEG ) )          !  minutes
           MNT  = INT( DSCR )                              !  int minutes
           DSCR = 60.0D0 * ( DSCR - DBLE( MNT ) )          !  seconds
           TPAIN( 4 ) = DSCR + 1000.0D0*( MNT + 1000*DEG ) !  dddmmmsss.sssD0

           DSCR = P_GAM
           DEG  = INT( DSCR )                              !  int degrees
           DSCR = 60.0D0 * ( DSCR - DBLE( DEG ) )          !  minutes
           MNT  = INT( DSCR )                              !  int minutes
           DSCR = 60.0D0 * ( DSCR - DBLE( MNT ) )          !  seconds
           TPAIN( 5 ) = DSCR + 1000.0D0*( MNT + 1000*DEG ) !  dddmmmsss.sssD0

           DSCR = YCENT
           DEG  = INT( DSCR )                              !  int degrees
           DSCR = 60.0D0 * ( DSCR - DBLE( DEG ) )          !  minutes
           MNT  = INT( DSCR )                              !  int minutes
           DSCR = 60.0D0 * ( DSCR - DBLE( MNT ) )          !  seconds
           TPAIN( 6 ) = DSCR + 1000.0D0*( MNT + 1000*DEG ) !  dddmmmsss.sssD0
           TPAIN( 7 ) = 0.0D0
           TPAIN( 8 ) = 0.0D0
           INSYS  = 4       !  Lambert conformal conic
           INZONE = 67
           INUNIT = 2       !  input units:  meters
           INSPH  = 19      !  normal sphere
           IPR    = 0       !  print error messages, if any
           JPR    = 1       !  do NOT print projection parameters
           LEMSG  = INIT3() !  unit number for log file
           LPARM  = LEMSG   !  projection parameters file

        CASE( UTMGRD3 ) 

           TPAIN = 0.0D0
           INSYS  = 1       !  Universal Transverse Mercator
           INZONE = NINT( P_ALP )
           INUNIT = 2       !  input units:  meters
           INSPH  = 19      !  normal sphere 
           IPR    = 0       !  print error messages, if any
           JPR    = 1       !  do NOT print projection parameters
           LEMSG  = INIT3() !  unit number for log file
           LPARM  = LEMSG   !  projection parameters file

        CASE( LATGRD3 )

           TPAIN = 0.0D0
           INSYS  = 0       !  geographic (=Lat-Lon)
           INZONE = 0
           INUNIT = 4       !  input units:  meters
           INSPH  = 19      !  normal sphere
           IPR    = 0       !  print error messages, if any
           JPR    = 1       !  do NOT print projection parameters
           LEMSG  = INIT3() !  unit number for log file
           LPARM  = LEMSG   !  projection parameters file

        CASE DEFAULT

            WRITE( MESG, '( A, I4, 2X, A )' ) 
     &      'Grid type', GDTYP, 
     &      'not supported (does Lat-Lon, Lambert, and UTM only)'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

        END SELECT  ! if grid not selected

C.........  Get grid parameters from grid description file for output grid.        
        IF ( .NOT. DSCM3GRD( GDNAM2, GDESC2, COORD2, GDTYP2, COORUNIT2,
     &                       P_ALP2, P_BET2, P_GAM2, XCENT2, 
     &                       YCENT2, XORIG2, YORIG2, XCELL2,
     &                       YCELL2, NCOLS2, NROWS2, NTHIK2 ) ) THEN

            MESG = 'Could not get Models-3 grid description.'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

        END IF

C.........  Compute intermediate grid-overlay variables
        IF ( GDTYP2 .EQ. LAMGRD3 ) THEN

            TPOUT( 1 ) = 0.0D0
            TPOUT( 2 ) = 0.0D0

            DSCR = P_ALP2
            DEG  = INT( DSCR )                              !  int degrees
            DSCR = 60.0D0 * ( DSCR - DBLE( DEG ) )          !  minutes
            MNT  = INT( DSCR )                              !  int minutes
            DSCR = 60.0D0 * ( DSCR - DBLE( MNT ) )          !  seconds
            TPOUT( 3 ) = DSCR + 1000.0D0*( MNT + 1000*DEG ) !  dddmmmsss.sssD0

            DSCR = P_BET2
            DEG  = INT( DSCR )                              !  int degrees
            DSCR = 60.0D0 * ( DSCR - DBLE( DEG ) )          !  minutes
            MNT  = INT( DSCR )                              !  int minutes
            DSCR = 60.0D0 * ( DSCR - DBLE( MNT ) )          !  seconds
            TPOUT( 4 ) = DSCR + 1000.0D0*( MNT + 1000*DEG ) !  dddmmmsss.sssD0

            DSCR = P_GAM2
            DEG  = INT( DSCR )                              !  int degrees
            DSCR = 60.0D0 * ( DSCR - DBLE( DEG ) )          !  minutes
            MNT  = INT( DSCR )                              !  int minutes
            DSCR = 60.0D0 * ( DSCR - DBLE( MNT ) )          !  seconds
            TPOUT( 5 ) = DSCR + 1000.0D0*( MNT + 1000*DEG ) !  dddmmmsss.sssD0

            DSCR = YCENT2
            DEG  = INT( DSCR )                              !  int degrees
            DSCR = 60.0D0 * ( DSCR - DBLE( DEG ) )          !  minutes
            MNT  = INT( DSCR )                              !  int minutes
            DSCR = 60.0D0 * ( DSCR - DBLE( MNT ) )          !  seconds
            TPOUT( 6 ) = DSCR + 1000.0D0*( MNT + 1000*DEG ) !  dddmmmsss.sssD0

            TPOUT( 7 ) = 0.0D0
            TPOUT( 8 ) = 0.0D0
            IOSYS  = 4       !  Lambert conformal conic
            IOZONE = 69
            IOUNIT = 2       !  output units:  meters
            IOSPH  = 19      !  normal sphere

            OUTTYPE = 'LAMBERT'
            OUTUNIT = 'METERS'

        ELSE IF ( GDTYP2 .EQ. UTMGRD3 ) THEN

            TPOUT  = 0.0D0
            IOSYS  = 1       !  Universal Transverse Mercator
            IOZONE = NINT( P_ALP2 )
            IOUNIT = 2       !  input units:  meters
            IOSPH  = 19      !  normal sphere

            OUTTYPE = 'LAT-LON'
            OUTUNIT = 'DEGREES'

        ELSE IF ( GDTYP2 .EQ. LATGRD3 ) THEN

            TPOUT  = 0.0D0
            IOSYS  = 0       !  geographic coords (=Lat-Lon)
            IOZONE = 0
            IOUNIT = 4       !  output units:  degrees
            IOSPH  = 19      !  normal sphere

            OUTTYPE = 'UTM'
            OUTUNIT = 'METERS'

        ELSE

            WRITE( MESG, '( A, I4, 2X, A )' ) 
     &      'Grid type', GDTYP2,
     &      'not supported (does Lat-Lon, Lambert, and UTM only)'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

        END IF                  ! if dscgrid() failed, or if non-Lambert grid

C........  Compute X-Y coordinates of fine grid relative to coarse grid
C........  and then mapping C2 from fine cells to coarse cells:

        ALLOCATE( C2( NCOLS, NROWS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'C2', PROGNAME )
        ALLOCATE( XX( NCOLS, NROWS ), STAT = IOS )
        CALL CHECKMEM( IOS, 'XX', PROGNAME )
        ALLOCATE( YY( NCOLS, NROWS ), STAT = IOS )
        CALL CHECKMEM( IOS, 'YY', PROGNAME )

        X0 = XORIG - 0.5 * XCELL
        Y0 = YORIG - 0.5 * YCELL
        DDX = 1.0 / XCELL2
        DDY = 1.0 / YCELL2

        DO  R = 1, NROWS       !  loop through cells of input grid
        DO  C = 1, NCOLS

            CRDIN( 1 ) = X0 + FLOAT( C ) * XCELL
            CRDIN( 2 ) = Y0 + FLOAT( R ) * YCELL
            CALL GTPZ0( CRDIN, INSYS, INZONE, TPAIN, INUNIT, INSPH, 
     &                  IPR, JPR, LEMSG, LPARM, CRDIO, IOSYS, IOZONE, 
     &                  TPOUT, IOUNIT, LN27, LN83, FN27, FN83, LENGTH,
     &                  IFLG )

            IF ( IFLG .NE. 0 ) THEN
                IFLG = MAX( MIN( 9, IFLG ), 1 )     !  trap between 1 and 9
                MESG = 'Failure in LAM2LL  for ' // GRDNM
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF
            
            XX( C,R ) = SNGL( CRDIO( 1 ) )
            YY( C,R ) = SNGL( CRDIO( 2 ) )

            XZ = XX( C,R ) - XORIG2
            YZ = YY( C,R ) - YORIG2
            COL = 1 + INT( DDX * XZ )
            ROW = 1 + INT( DDY * YZ )

C.............  If center of input cell (col,row) is outside output grid,
C               then store missing intersection.
            IF ( XZ  .LT. 0.0     .OR.  YZ  .LT. 0.0    .OR.  
     &           COL .GT. NCOLS2  .OR.  ROW .GT. NROWS2 ) THEN
                C2( C,R ) = IMISS3

C.............  Store which cell in output grid for current input cell (col,row)
            ELSE
                C2( C,R ) = COL + NCOLS2 * ( ROW - 1 )

            END IF

        END DO
        END DO          !  end looping through input grid

C.........  Open output file
        XDEV = PROMPTFFILE( 
     &         'Enter logical name for output GRIDDING SURROGATE file',
     &           .FALSE., .TRUE., 'OUTFILE', PROGNAME )

        CALL M3MSG2( 'Reading gridding surrogates...' )

C.........  Allocate memory for and read the gridding surrogates file,
C           extracting data for a subgrid, if necessary
        CALL RDSRG( .FALSE., SDEV, SRGFMT, NROWS, NCOLS )

C.........  Allocate memory for output surrogates assuming that number of
C           records from input file will > output records.
        ALLOCATE( INDXA( NSRGREC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'INDXA', PROGNAME )
        ALLOCATE( CELLA( NSRGREC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'CELLA', PROGNAME )
        ALLOCATE( FIPSA( NSRGREC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'FIPSA', PROGNAME )
        ALLOCATE( SSCSA( NSRGREC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'SSCSA', PROGNAME )
        ALLOCATE( FRACA( NSRGREC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'FRACA', PROGNAME )

        CALL M3MSG2( 'Regridding surrogates...' )

C.........  Regrid the surrogates based on previously computed factors
        OUTNSRG = 0
        DO J = 1, NSRGREC

            K = IDXSRGA( J )
            C = SCELLA ( K )

            ROW = 1 + INT( (C-1) / NCOLS )
            COL = C - NCOLS * ( ROW-1 )

C.............  Do not store output cells if input cell does not intersect
C               any output cells

            IF ( C2( COL,ROW ) .EQ. IMISS3 ) THEN
                CYCLE

            ELSE

                OUTNSRG = OUTNSRG + 1

                IF ( OUTNSRG .GT. NSRGREC ) CYCLE

                INDXA( OUTNSRG ) = OUTNSRG
                CELLA( OUTNSRG ) = C2( COL,ROW )
                FIPSA( OUTNSRG ) = STR2INT( SFIPSA ( K ) )
                SSCSA( OUTNSRG ) = SSRGIDA( K )
                FRACA( OUTNSRG ) = SFRACA ( K )

            END IF
            
        END DO

C...........   Give error if couldn't store all of the output intersections
        IF( OUTNSRG .GT. NSRGREC ) THEN

            MESG = 'INTERNAL ERROR: Assumption about memory ' //
     &             'allocation for output surrogates is false.'
            CALL M3MSG2( MESG )
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

        END IF

C...........   Sort and process the surrogates:
        CALL M3MSG2( 'Sorting output gridding surrogates...' )

        CALL SORTI3( OUTNSRG, INDXA, SSCSA, CELLA, FIPSA )

        CALL M3MSG2 ( 'Creating output gridding surrogate file...' )

C.........  Write header line
        MESG = '#GRID ' // GDNAM2
        L = LEN_TRIM( MESG )

        L1 = LEN_TRIM( OUTTYPE )
        L2 = LEN_TRIM( OUTUNIT )
        WRITE( XDEV, 93980 ) MESG( 1:L ), XORIG2, YORIG2, XCELL2,
     &         YCELL2, NCOLS2, NROWS2, NTHIK2, OUTTYPE( 1:L1 ), 
     &         OUTUNIT( 1:L2 ), P_ALP2, P_BET2, P_GAM2, XCENT2, YCENT2

C.........  Loop through output entries and write them
        LCEL = IMISS3
        LFIP = IMISS3
        LSSC = IMISS3
        M    = 0
        N    = 0
        FRAC = 0.0

        DO  J = 1, OUTNSRG

            K = INDXA ( J )
            C = CELLA ( K )
            F = FIPSA ( K )
            S = SSCSA ( K )

C.............  If same FIPS code, cell, and surrogate, sum surrogate fraction
            IF ( F .EQ. LFIP  .AND.   
     &           C .EQ. LCEL  .AND.   
     &           S .EQ. LSSC ) THEN

                FRAC = FRAC + FRACA( K )
                N    = N + 1

C.............  Otherwise...
            ELSE
            
                IF ( N .GT. 0 ) THEN

                    ROW = 1 + INT( (LCEL-1) / NCOLS2 )
                    COL = LCEL - NCOLS2 * ( ROW - 1 )

                    WRITE( XDEV, 93020, IOSTAT=IOS )    
     &                  LSSC, LFIP, COL, ROW, FRAC

                    IF ( IOS .GT. 0 ) THEN
                        EFLAG = .TRUE.
                        WRITE( MESG, 94010 ) 'ERROR: I/O error', IOS, 
     &                         'writing output file at line', M
                    END IF

                END IF
                
                M = M + 1
                N = 1
                LCEL = C
                LFIP = F
                
                LSSC = S
                FRAC = FRACA( K )
            
            END IF            !  if new record encountered

        END DO

C.........  Write the last entry in the output file
        IF ( N .GT. 0 ) THEN
            ROW = 1 + INT( (LCEL-1) / NCOLS2 )
            COL = LCEL - NCOLS2 * ( ROW-1 )

            FRAC = FRAC / FLOAT( N )

            WRITE( XDEV, 93020, IOSTAT=IOS )     
     &              LSSC, F, COL, ROW, FRAC

C.............  Write error message, if needed
            IF ( IOS .GT. 0 ) THEN
                EFLAG = .TRUE.
                WRITE( MESG, 94010 ) 'ERROR: I/O error', IOS,     
     &                 'writing output file at line', M
            END IF

        END IF

C.........  Abort if error
        IF( EFLAG ) THEN

            MESG = 'Problem writing output file.'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

        END IF
                
C..........   Normal completion
        MESG = ' ' 
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 0 )

C******************  FORMAT  STATEMENTS   ******************************

C...........   Formatted file I/O formats............ 93xxx

93000   FORMAT( A )

93020   FORMAT( I3, 1X, I6.5, I8, I8, 1X, F11.8 )

93980   FORMAT( A, 1X, 4( E15.8, 1X ), 3( I8, 1X ), 2( A, 1X ), 
     &          5( E15.8, 1X ) )

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I10, :, 2X ) )

        END PROGRAM SURGTOOL
