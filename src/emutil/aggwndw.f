        PROGRAM AGGWNDW

C***********************************************************************
C  program body starts at line 135
C
C  DESCRIPTION: 
C      Aggregates or extracts a window of emissions data from an input
C      (fine) to an output (coarse) grid.
C      Inputs SMOKE-formatted emissions data for a "fine" input grid
C      and a grid definition for a "coarse" output grid.
C      It produces an aggregate (sum) "coarse" grid emissions file.
C
C  PRECONDITIONS REQUIRED:
C      - Input & output map projections (alpha, beta, gamma,
C        x center, y center) must be the same.
C      - Grid cell size (x cell, y cell) of the coarse output grid
C        must be an integer multiple of the fine input grid:
C        4 km to 12 km is okay; 4 km to 10 km is not.
C      - Grid cell boundaries must be the concentric:
C        difference in x origin and y origin on the coarse grid must be
C        an integer multiple of the original fine grid's cell size.
C      - The coarse output grid must be entirely contained within the
C        fine input grid.
C      - Input logical file names are defined.
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C       CHECKMEM
C       DBLERR
C       DESC3
C       DSCM3GRD
C       INIT3
C       INITEM
C       M3EXIT
C       M3MSG2
C       NEXTIME
C       PROMPTMFILE
C       READ3
C       WRITE3
C
C  REVISION  HISTORY:
C       Created 4/2005 by C. Mattocks
C
C***********************************************************************
C
C Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
C                System
C File: @(#)$Id$
C
C COPYRIGHT (C) 2005, Environmental Modeling for Policy Development
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

C.........  MODULES for public variables
C.........  This module contains the global variables for the 3-d grid
        USE MODGRID, ONLY: GRDNM, COORD, GDTYP, P_ALP, P_BET, P_GAM,
     &                     XCENT, YCENT, XORIG, YORIG, XCELL, YCELL,
     &                     NCOLS, NROWS

        IMPLICIT NONE

C.........  INCLUDES
        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters
        INCLUDE 'FDESC3.EXT'    !  I/O API file description data structures.
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations

C.........  EXTERNAL FUNCTIONS
        CHARACTER(2)    CRLF
        LOGICAL         DBLERR
        LOGICAL         DSCM3GRD
        INTEGER         PROMPTFFILE
        CHARACTER(16)   PROMPTMFILE
        LOGICAL         SETENVVAR
        LOGICAL         ENVYN

        EXTERNAL        CRLF, DBLERR, DSCM3GRD, PROMPTFFILE, 
     &                  PROMPTMFILE, SETENVVAR, ENVYN

C.........  LOCAL PARAMETERS
        CHARACTER(50), PARAMETER :: 
     &  CVSW = '$Name SMOKEv4.8.1_Jan2021$'  ! CVS release tag

C.........  LOCAL VARIABLES

C.........  Allocatable arrays
        INTEGER, ALLOCATABLE :: OUTROW( :,: )       ! output row for each input cell
        INTEGER, ALLOCATABLE :: OUTCOL( :,: )       ! output column for each input cell
        
        REAL,    ALLOCATABLE :: VDATA_IN ( :,: )    ! input data
        REAL,    ALLOCATABLE :: VDATA_OUT( :,:,: )  ! output data

C.........  Logical file names and unit numbers
        INTEGER         LDEV              !  unit number for log file
        INTEGER         RDEV              !  unit number of report file

        CHARACTER(16)   INAME             !  logical name for input  emissions data file
        CHARACTER(16)   ONAME             !  logical name for output emissions data file

C.........  Other local variables
        INTEGER     IOS                         ! I/O status
        INTEGER     I, J, K, L, T               ! counters, subscripts
        INTEGER     COL, ROW                    ! tmp row and column
        INTEGER     JDATE                       ! loop date
        INTEGER     JTIME                       ! loop time
        INTEGER     NTHIK                       ! thickness of grid boundary
        INTEGER     NROWS_IN                    ! number of rows in input grid
        INTEGER     NCOLS_IN                    ! number of columns in input grid
        INTEGER     NSTEPS                      ! number of time steps in input
        INTEGER     POLLWIDTH                   ! width of pollutant names
        INTEGER     UNITWIDTH                   ! width of emissions units
        INTEGER     IDIVSR                      ! No of aggregating cells (divisor)

        REAL(8)     DDX, DDY                    ! inverse of grid cell size
        REAL(8)     X0, Y0                      ! shifted origin
        REAL(8)     XX, YY                      ! grid cell center coordinates
        REAL(8)     XE_IN, XE_OUT               ! right boundary of input and output grids
        REAL(8)     YN_IN, YN_OUT               ! upper boundary of input and output grids
        REAL(8)     CHK_X                       ! check subgrid even with grid in x
        REAL(8)     CHK_Y                       ! check subgrid even with grid in y
        
        DOUBLE PRECISION  INTOTAL               ! input emissions total for reporting
        DOUBLE PRECISION  OUTTOTAL              ! output emissions total for reporting
        DOUBLE PRECISION  PCTDIFF               ! percent difference in totals

        LOGICAL     EMISFLAG                    ! true: aggregating mass emissions across grid cells
        LOGICAL ::  EFLAG = .FALSE.             ! true: error found

        CHARACTER(80)    GDESC                  ! grid description
        CHARACTER(80)    COORUNIT               ! coordinate units
        CHARACTER(100)   HDRFMT                 ! format for report header
        CHARACTER(100)   RPTFMT                 ! format for report output
        CHARACTER(100)   TMPFMT                 ! temporary format
        CHARACTER(256)   MESG                   ! message buffer

        CHARACTER(16) :: PROGNAME = 'AGGWNDW'   ! program name

C***********************************************************************
C   begin body of program AGGWNDW

        LDEV = INIT3()

C.........  Write out copyright, version, web address, header info,
C           and prompt to continue running the program
        CALL INITEM( LDEV, CVSW, PROGNAME )

C.........  Open input emissions file
        MESG = 'Enter logical name for INPUT DATA file'
        INAME = PROMPTMFILE( MESG, FSREAD3, 'INFILE', PROGNAME )

C.........  Aggregating emission mass or not
        MESG = 'Aggregating gridded emissions or not?'
        EMISFLAG = ENVYN( 'AGGREGATE_EMIS_YN', MESG, .TRUE., IOS )

C.........  Get name of output grid and read grid parameters
        IF( .NOT. DSCM3GRD( GRDNM, GDESC, COORD, GDTYP, COORUNIT,
     &                      P_ALP, P_BET, P_GAM, XCENT, YCENT,
     &                      XORIG, YORIG, XCELL, YCELL, NCOLS,
     &                      NROWS, NTHIK ) ) THEN
            MESG = 'Could not get grid description'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

C.........  Get description in input file
        IF( .NOT. DESC3( INAME ) ) THEN
            MESG = 'Could not get description of file ' // 
     &             TRIM( INAME )
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

C.........  Compare input and output grids
        CALL M3MSG2( 'Checking input and output grids...' )

C.........  Check projection parameters
        IF ( DBLERR( P_ALP3D, P_ALP ) .OR.
     &       DBLERR( P_BET3D, P_BET ) .OR.
     &       DBLERR( P_GAM3D, P_GAM ) .OR.
     &       DBLERR( XCENT3D, XCENT ) .OR.
     &       DBLERR( YCENT3D, YCENT )      ) THEN

            EFLAG = .TRUE.
            MESG = 'ERROR: Input and output grid projections do not match'
            CALL M3MSG2( MESG )
        END IF

C.........  Check that output grid cell size is a multiple of input cell size
        IF( .NOT. EFLAG ) THEN
            CHK_X = XCELL / XCELL3D
            CHK_Y = YCELL / YCELL3D

            IDIVSR = INT( CHK_X * CHK_X )
   
            IF ( DBLERR( CHK_X, DNINT( CHK_X ) ) .OR.
     &           DBLERR( CHK_Y, DNINT( CHK_Y ) )      ) THEN
                EFLAG = .TRUE.
                MESG = 'ERROR: Output grid cell size is not an ' //
     &                 'integer multiple of input grid cell size'
                CALL M3MSG2( MESG )
            END IF
        END IF

C.........  Check that the grid boundaries are aligned
        IF( .NOT. EFLAG ) THEN
            CHK_X = ( XORIG - XORIG3D ) / XCELL3D
            CHK_Y = ( YORIG - YORIG3D ) / YCELL3D
            
            IF( CHK_X /= 0. ) THEN
                IF( ABS( (DNINT(CHK_X) - CHK_X) / CHK_X ) > 0.001 ) THEN
                    EFLAG = .TRUE.
                    MESG = 'ERROR: Input and output grid cell ' //
     &                'boundaries are not aligned in the x-direction'
                    CALL M3MSG2( MESG )
                END IF
            END IF
            
            IF( CHK_Y /= 0. ) THEN
                IF( ABS( (DNINT(CHK_Y) - CHK_Y) / CHK_Y ) > 0.001 ) THEN
                    EFLAG = .TRUE.
                    MESG = 'ERROR: Input and output grid cell ' //
     &                'boundaries are not aligned in the y-direction'
                    CALL M3MSG2( MESG )
                END IF
            END IF
        END IF

C.........  Check that the output grid is entirely within the input grid
        IF( .NOT. EFLAG ) THEN
            XE_IN  = XORIG3D + NCOLS3D * XCELL3D   ! East edges
            XE_OUT = XORIG   + NCOLS   * XCELL
    
            YN_IN  = YORIG3D + NROWS3D * YCELL3D   ! North edges
            YN_OUT = YORIG   + NROWS   * YCELL

            IF( XORIG3D > XORIG  .OR.
     &          YORIG3D > YORIG  .OR.
     &          XE_IN   < XE_OUT .OR.
     &          YN_IN   < YN_OUT      ) THEN

                EFLAG = .TRUE.
                MESG = 'ERROR: Output grid is not entirely ' //
     &                 'contained within input grid.'
                CALL M3MSG2( MESG )
            END IF
        END IF

        IF( EFLAG ) THEN
            MESG = 'Input and output grids are not compatible'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

C.........  For each input grid cell, determine the corresponding output
C           grid cell
        CALL M3MSG2( 'Calculating output grid information...' )

        ALLOCATE( OUTROW( NCOLS3D, NROWS3D ), STAT=IOS )
        CALL CHECKMEM( IOS, 'OUTROW', PROGNAME )
        ALLOCATE( OUTCOL( NCOLS3D, NROWS3D ), STAT=IOS )
        CALL CHECKMEM( IOS, 'OUTCOL', PROGNAME )

        DDX = 1. / XCELL
        DDY = 1. / YCELL

C.........  Calculate origin offset        
        X0 = XORIG3D - XCELL3D * 0.5
        Y0 = YORIG3D - YCELL3D * 0.5
        
        DO J = 1, NROWS3D
        
            DO I = 1, NCOLS3D

C.................  Calculate center of current grid cell            
                XX = X0 + XCELL3D * I
                YY = Y0 + YCELL3D * J

C.................  Calculate column in output grid
                CHK_X = DDX * ( XX - XORIG )
                IF( CHK_X < 0. ) THEN
                    OUTROW( I,J ) = 0
                    OUTCOL( I,J ) = 0
                    CYCLE
                END IF
                
                COL = 1 + INT( CHK_X )
                IF( COL > NCOLS .OR. COL <= 0 ) THEN
                    OUTROW( I,J ) = 0
                    OUTCOL( I,J ) = 0
                    CYCLE
                ELSE
                    OUTCOL( I,J ) = COL
                END IF

C.................  Calculate row in output grid
                CHK_Y = DDY * ( YY - YORIG )
                IF( CHK_Y < 0. ) THEN
                    OUTROW( I,J ) = 0
                    OUTCOL( I,J ) = 0
                    CYCLE
                END IF

                ROW = 1 + INT( CHK_Y )
                IF( ROW > NROWS .OR. ROW <= 0 ) THEN
                    OUTROW( I,J ) = 0
                    OUTCOL( I,J ) = 0
                    CYCLE
                ELSE
                    OUTROW( I,J ) = ROW
                END IF
            
            END DO
        
        END DO

C.........  Allocate memory for storing input and output data
        ALLOCATE( VDATA_IN( NCOLS3D, NROWS3D ), STAT=IOS )
        CALL CHECKMEM( IOS, 'VDATA_IN', PROGNAME )

        ALLOCATE( VDATA_OUT( NCOLS, NROWS, NLAYS3D ), STAT=IOS )
        CALL CHECKMEM( IOS, 'VDATA_OUT', PROGNAME )

C.........  Save input grid information
        NROWS_IN = NROWS3D
        NCOLS_IN = NCOLS3D
        NSTEPS   = MXREC3D
        
C.........  Build header for output file; use info from input emissions file
C           and just change grid information
        XORIG3D = XORIG
        YORIG3D = YORIG
        XCELL3D = XCELL
        YCELL3D = YCELL
        NCOLS3D = NCOLS
        NROWS3D = NROWS
        GDNAM3D = GRDNM

C.........  Open output emissions file
        MESG = 'Enter logical name for OUTPUT DATA file'
        ONAME = PROMPTMFILE( MESG, FSUNKN3, 'OUTFILE', PROGNAME )

C.........  Open report file
        MESG = 'Enter logical name for OUTPUT REPORT file'
        RDEV = PROMPTFFILE( MESG, .FALSE., .TRUE., 'REPFILE', PROGNAME )

C.........  Build report header and format string for report output

C.........  Determine maximum width of pollutant names and units
        POLLWIDTH = LEN_TRIM( 'Pollutant' )
        UNITWIDTH = 0
        DO I = 1, NVARS3D
            POLLWIDTH = MAX( POLLWIDTH, LEN_TRIM( VNAME3D( I ) ) )
            UNITWIDTH = MAX( UNITWIDTH, LEN_TRIM( UNITS3D( I ) ) )
        END DO

        HDRFMT = "(A7, '; ', A6, '; ', A5, '; ', A"
        RPTFMT = "(I7, '; ', I6, '; ', I5, '; ', A"

        IF( POLLWIDTH < 10 ) THEN
            TMPFMT = '(A,I1,A)'
        ELSE
            TMPFMT = '(A,I2,A)'
        END IF
        
        WRITE( HDRFMT, TMPFMT ) TRIM( HDRFMT ), POLLWIDTH, ", '; ', A"
        WRITE( RPTFMT, TMPFMT ) TRIM( RPTFMT ), POLLWIDTH,
     &                          ", '; ', F15.2, ' ', A"

        IF( UNITWIDTH < 10 ) THEN
            TMPFMT = '(A,I1,A,I1,A)'
        ELSE
            TMPFMT = '(A,I2,A,I2,A)'
        END IF

        WRITE( HDRFMT, '(A,I2,A,I2,A)' ) TRIM( HDRFMT ), UNITWIDTH + 16,
     &                                 ", '; ', A", UNITWIDTH + 16, ")"
        
        WRITE( RPTFMT, TMPFMT ) TRIM( RPTFMT ), UNITWIDTH,
     &                          ", '; ', F15.2, ' ', A", UNITWIDTH, ")"

        WRITE( RDEV, HDRFMT ) 'Date', 'Hour', 'Layer', 'Pollutant',
     &                        'Input emissions', 'Output emissions'

C.....  For each variable in the original input emissions file,
C       at each time step:
C       - Read the emissions data from the original file using READ3()
C       - Loop through the cells of the input grid and
C         calculate emissions for the new grid
C       - When multiple input grid cells map to the same output grid
C         cell, sum the emissions
C       - Write the emissions data to the output file using WRITE3()

C.........  Loop through time steps of emissions data
        JDATE = SDATE3D          ! Initialize date
        JTIME = STIME3D          ! Initialize time
        DO T = 1, NSTEPS         ! Loop over time steps

            DO L = 1, NVARS3D    ! Loop over variables

C.................  Initialize output data array
                VDATA_OUT = 0.00000

                DO K = 1, NLAYS3D

C.....................  Initialize input data array
                    VDATA_IN = 0.

C.....................  Read variable from input emissions file
                    IF( .NOT. READ3( INAME, VNAME3D( L ), K,
     &                               JDATE, JTIME, VDATA_IN ) ) THEN
                        MESG = 'Could not read variable ' // 
     &                         TRIM( VNAME3D( I ) ) // ' from file ' // 
     &                         TRIM( INAME )
                        CALL M3EXIT( PROGNAME, JDATE, JTIME, MESG, 2 )
                    END IF

C.....................  Loop over grid cells
                    DO J = 1, NROWS_IN
                        DO I = 1, NCOLS_IN
                            IF( OUTROW( I,J ) == 0 .OR.
     &                          OUTCOL( I,J ) == 0      ) CYCLE
                            
                            ROW = OUTROW( I,J )
                            COL = OUTCOL( I,J )
                            
                            VDATA_OUT( COL,ROW,K ) = 
     &                          VDATA_OUT( COL,ROW,K ) + VDATA_IN( I,J )
     
                        END DO
                    END DO

C.....................  Calculate report information
                    INTOTAL = 0.00000
                    DO J = 1, NROWS_IN
                        DO I = 1, NCOLS_IN
                            IF( OUTROW( I,J ) == 0 .OR.
     &                          OUTCOL( I,J ) == 0      ) CYCLE

                            INTOTAL = INTOTAL + VDATA_IN( I,J )
                        END DO
                    END DO
                    
                    OUTTOTAL = 0.00000
                    DO J = 1, NROWS
                        DO I = 1, NCOLS
                            OUTTOTAL = OUTTOTAL + VDATA_OUT( I,J,K )
                        END DO
                    END DO
                    
C.....................  Output report information
                    WRITE( RDEV, RPTFMT ) JDATE, JTIME, K, VNAME3D( L ),
     &                  INTOTAL, UNITS3D( L ), OUTTOTAL, UNITS3D( L )
                    
C.....................  Check percent difference between totals and flag
C                       if difference is too big; "too big" is pretty much
C                       a guess here (0.01%)
                    IF( INTOTAL /= 0. ) THEN
                        PCTDIFF = 100. * ( OUTTOTAL-INTOTAL ) / INTOTAL
                    
                        IF( ABS( PCTDIFF ) > 0.01 ) THEN
                            EFLAG = .TRUE.
                            WRITE( MESG, 94020 )
     &                        'ERROR: Input and output emissions ' //
     &                        'differ by ', PCTDIFF, ' percent for' //
     &                        CRLF() // BLANK10 // 'pollutant ' //
     &                        TRIM( VNAME3D( L ) ) // ' at time ',
     &                        JDATE, ':', JTIME
                            CALL M3MESG( MESG )
                        END IF
                    END IF

C.....................   Aggregate non-mass values like percentage or ratio
                    IF( .NOT. EMISFLAG ) THEN
                        DO J = 1, NROWS
                            DO I = 1, NCOLS
                               VDATA_OUT( I,J,K) = VDATA_OUT( I,J,K )/
     &                                             FLOAT( IDIVSR )
                            END DO
                        END DO
                    END IF

                END DO  ! End loop over layers

C.................  Output variable to new emissions file
                IF( .NOT. WRITE3( ONAME, VNAME3D( L ), 
     &                            JDATE, JTIME, VDATA_OUT ) ) THEN
                    MESG = 'Could not write variable ' // 
     &                     TRIM( VNAME3D( I ) ) // ' to file ' //
     &                     TRIM( ONAME )
                    CALL M3EXIT( PROGNAME, JDATE, JTIME, MESG, 2 )
                END IF

            END DO   ! End loop over variables
            
C.............  Increment date & time by time step
            CALL NEXTIME( JDATE, JTIME, TSTEP3D )
        
        END DO   ! End loop over time steps

C.........  Check if any in and out totals did not match
        IF( EFLAG ) THEN
            MESG = 'Input and output emissions differ by more ' //
     &             'than 0.01%'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF
 
        CALL M3EXIT( PROGNAME, 0, 0, ' ', 0 )

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats............ 94xxx

94020   FORMAT( A, E12.5, A, I7, A, I6.6 )

        END PROGRAM AGGWNDW
