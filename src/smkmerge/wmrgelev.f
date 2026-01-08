
        SUBROUTINE WMRGELEV( VNAME, NSRC, NMAJOR, JDATE, JTIME )

C***********************************************************************
C  subroutine WMRGELEV body starts at line
C
C  DESCRIPTION:
C      The purpose of this subroutine is to write out an ASCII elevated point
C      source emissions file. The default format is the UAM/CAMx elevated
C      sources file.
C
C  PRECONDITIONS REQUIRED:  
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C       Created 4/2000 by M. Houyoux
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
C****************************************************************************

C.........  MODULES for public variables
C.........  This module contains the major data structure and control flags
        USE M3UTILIO

        USE MODMERGE, ONLY: EXPLFLAG, EMLAYS, PVSDATE, PVSTIME, PVNAME,
     &                      PENAME, PSDEV, NPSRC, NMSPC, EMNAM,
     &                      SDATE, STIME, EDATE, ETIME, TSTEP,
     &                      JSTACK, EVDEV, LFRAC

C...........   This module is the source inventory arrays
        USE MODSOURC, ONLY: CSOURC

C.........  This module contains arrays for plume-in-grid and major sources
        USE MODELEV, ONLY: NHRSRC, GRPCOL, GRPROW, GRPXX, GRPYY, 
     &                     GRPHT, GRPDM, GRPTK, GRPVE,
     &                     NGROUP, GRPGID, ELEVSRC, LPING,
     &                      ELEVSIDX, GROUPID, INDXH, ELEVEMIS

C.........  This module contains the global variables for the 3-d grid
        USE MODGRID, ONLY: NCOLS, NROWS, XCELL, YCELL, XORIG, YORIG,
     &                     GRDNM, P_ALP, GDTYP, VGTYP, VGLVS

        IMPLICIT NONE

C.........  INCLUDES:
        
        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
C        INCLUDE 'PARMS3.EXT'    !  I/O API parameters
C        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
C        INCLUDE 'FDESC3.EXT'    !  I/O API file desc. data structures

C.........  EXTERNAL FUNCTIONS and their descriptions:
        
C       CHARACTER(2)    CRLF
C       INTEGER         ENVINT
C       REAL            ENVREAL
C       INTEGER         FIND1
C       INTEGER         INDEX1
C       INTEGER         STR2INT

C        EXTERNAL        CRLF, ENVINT, ENVREAL, FIND1, INDEX1, STR2INT

C.........  SUBROUTINE ARGUMENTS
        CHARACTER(*), INTENT (IN) :: VNAME   ! variable name to output
        INTEGER     , INTENT (IN) :: NSRC    ! no. sources
        INTEGER     , INTENT (IN) :: NMAJOR  ! no. elevated sources
        INTEGER     , INTENT (IN) :: JDATE   ! Julian date to output (YYYYDDD)
        INTEGER     , INTENT (IN) :: JTIME   ! time to output (HHMMSS)

C.........  Local parameters
        INTEGER, PARAMETER :: STKLEN3 = FIPLEN3 + PLTLEN3 + CHRLEN3

C.........  Variables allocated by module settings...
        INTEGER, ALLOCATABLE, SAVE :: ESRTIDX( : ) ! major srcs sorting index
        INTEGER, ALLOCATABLE, SAVE :: EIDXSRC( : ) ! Source ID per major source
        INTEGER, ALLOCATABLE, SAVE :: EIDX2  ( : ) ! another index
        INTEGER, ALLOCATABLE, SAVE :: ELAYER ( : ) ! addt'l srcs for explicit

        REAL, ALLOCATABLE :: EOUTHT ( : ) ! output stack heights [m]
        REAL, ALLOCATABLE :: EOUTDM ( : ) ! output stack diameters [m]
        REAL, ALLOCATABLE :: EOUTTK ( : ) ! output exit tmprs [K]
        REAL, ALLOCATABLE :: EOUTVE ( : ) ! output exit velocities [m/hr]
        REAL, ALLOCATABLE :: EOUTXL ( : ) ! output x-loc [units of grid]
        REAL, ALLOCATABLE :: EOUTYL ( : ) ! output y-loc [units of grid]

        CHARACTER(STKLEN3), ALLOCATABLE, SAVE :: ECSRCA  (:)! FIPS//plt//stk
        CHARACTER(STKLEN3), ALLOCATABLE       :: EOUTCSRC(:)! FIPS//plt//stk

C.........  Allocatable array for fake stack heights for explicit plums
        REAL, ALLOCATABLE :: LAYRMID( : )

C.........  Fixed size arrays
        REAL            ESUM   ( MXLAYS3 )  ! emissions data sum

C.........  UAM-format specific variables
        INTEGER         DDEVOUT             ! diffbreak file unit number
        INTEGER         ESWITCH             ! 0=no, 1=yes - print vert mthds tbl
        INTEGER         GSWITCH             ! 0=no, 1=yes - print output grid
        INTEGER         LSWITCH             ! 0=no, 1=yes - print src locs tbl
        INTEGER      :: MDEVOUT = 0         ! metscalars file unit number
        INTEGER         MSWITCH             ! 0=no, 1=yes - print methods tbl
        INTEGER         NH                  ! tmp no. explicit sources
        INTEGER         NPARAM              ! no. paramaters for control pkt
        INTEGER         NULAYS              ! no. UAM model layers
        INTEGER         NZLOWR              ! no. layers below Diffbreak height
        INTEGER         NZUPPR              ! no. layers above Diffbreak height
        INTEGER      :: PDEVOUT = 0         ! PTSRCE output unit no.
        INTEGER      :: RDEVOUT = 0         ! regiontop file unit number
        INTEGER      :: TDEVOUT = 0         ! temperatur file unit number
        INTEGER         USWITCH             ! 0=no, 1=yes - print units table
        INTEGER         VSWITCH             ! 0=no, 1=yes - print values tbl
        INTEGER      :: WDEVOUT = 0         ! wind file unit number

        REAL            HTSUR               ! height of surface layer [m]
        REAL            HTLOWR              ! min cell ht b/w sfc and diffbr [m]
        REAL            HTUPPR              ! min cell ht b/w diffbr and top [m]

        CHARACTER(10)       :: SPCNAM       ! UAM-format species name
        CHARACTER(10), SAVE :: VTYPE        ! User-spec vertical method type
        CHARACTER(44)          NOTEDEF      ! Default note
        CHARACTER(44)          UNOTE        ! UAM file note from env variable
        CHARACTER(60)       :: FNOTE        ! UAM file header note

C.........  Other local variables
        INTEGER          I, J, K, KK, L, LL, LM, LN, M, N, S 

        INTEGER          COL                 ! tmp column
        INTEGER          EML                 ! tmp emissions layers
        INTEGER          FIP                 ! tmp co/st/cy code
        INTEGER          GID                 ! tmp stack group ID
        INTEGER          IBD                 ! beg 5-digit Julian date of step
        INTEGER          IBT                 ! beginning time of time step
        INTEGER          IED                 ! end 5-digit Julian date of step
        INTEGER          IET                 ! ending time of time step
        INTEGER          IOS                 ! i/o status
        INTEGER          JDATEP1             ! julian date, plus 1 time step
        INTEGER          JTIMEP1             ! time, plus 1 time step
        INTEGER          NEXPLOOP            ! no. for explicit plume loop
        INTEGER, SAVE :: NOUT                ! number output stacks
        INTEGER, SAVE :: PTIME = -9          ! previous call's time (HHMMSS)
        INTEGER          ROW                 ! tmp row

        REAL             DM                  ! tmp stack diameter [m]
        REAL             HT                  ! tmp stack height [m]
        REAL             XLOC                ! tmp x location
        REAL             XMAX                ! rightmost x location
        REAL             YLOC                ! tmp y location
        REAL             YMAX                ! topmost y location

        LOGICAL       :: EFLAG    = .FALSE.  ! true: error occurred
        LOGICAL       :: FIRSTIME = .TRUE.   ! true: first time routine called

        CHARACTER(16),SAVE :: OUTFMT       ! output format for elevated ASCII
        CHARACTER(100)     EFMT         ! output emissions foamat
        CHARACTER(200)     BUFFER       ! source chars buffer
        CHARACTER(300)     MESG         ! message buffer

        CHARACTER(FIPLEN3) CFIP         ! tmp country/state/county code
        CHARACTER(PLTLEN3) FCID         ! tmp facility ID
        CHARACTER(CHRLEN3) SKID         ! tmp stack ID
        CHARACTER(STKLEN3) ECS          ! stack elevated source chars
        CHARACTER(STKLEN3) PECS         ! tmp previous ECS
        CHARACTER(IOULEN3) GRDENV       ! gridded output units from envrmt

        CHARACTER(16) :: PROGNAME = 'WMRGELEV' ! program name

C***********************************************************************
C   begin body of subroutine WMRGELEV

C.........  For the first time the subroutine is called...
        IF( FIRSTIME ) THEN

            MESG = 'Setting up to output ASCII elevated file...'
            CALL M3MSG2( MESG )

C.............  Get environment variable settings
            
C.............  Get type of plume rise calculation from environment
            MESG = 'Elevated ASCII output file format'
            CALL ENVSTR( 'SMK_ASCIIELEV_FMT', MESG, 'UAM', OUTFMT, IOS )
            IF( VTYPE .EQ. 'UAM4' ) VTYPE = 'UAM'

C.............  Set size for allocating output elevated arrays
C.............  Adjustment for EXPLFLAG so that one record for each layer of
C               each explicit source can be inserted.
            I = NMAJOR
            IF ( EXPLFLAG ) I = I - NHRSRC + EMLAYS * NHRSRC

C.............  Allocate memory for local elevated sources arrays
            ALLOCATE( ESRTIDX( I ), STAT=IOS )
            CALL CHECKMEM( IOS, 'ESRTIDX', PROGNAME )
            ALLOCATE( EIDXSRC( I ), STAT=IOS )
            CALL CHECKMEM( IOS, 'EIDXSRC', PROGNAME )
            ALLOCATE( ECSRCA( I ), STAT=IOS )
            CALL CHECKMEM( IOS, 'ECSRCA', PROGNAME )
            ALLOCATE( EOUTCSRC( I ), STAT=IOS )
            CALL CHECKMEM( IOS, 'EOUTCSRC', PROGNAME )
            ALLOCATE( EOUTHT( I ), STAT=IOS )
            CALL CHECKMEM( IOS, 'EOUTHT', PROGNAME )
            ALLOCATE( EOUTDM( I ), STAT=IOS )
            CALL CHECKMEM( IOS, 'EOUTDM', PROGNAME )
            ALLOCATE( EOUTTK( I ), STAT=IOS )
            CALL CHECKMEM( IOS, 'EOUTTK', PROGNAME )
            ALLOCATE( EOUTVE( I ), STAT=IOS )
            CALL CHECKMEM( IOS, 'EOUTVE', PROGNAME )
            ALLOCATE( EOUTXL( I ), STAT=IOS )
            CALL CHECKMEM( IOS, 'EOUTXL', PROGNAME )
            ALLOCATE( EOUTYL( I ), STAT=IOS )
            CALL CHECKMEM( IOS, 'EOUTYL', PROGNAME )
            ALLOCATE( EIDX2( I ), STAT=IOS )
            CALL CHECKMEM( IOS, 'EIDX2', PROGNAME )
            ALLOCATE( ELAYER( I ), STAT=IOS )
            CALL CHECKMEM( IOS, 'ELAYER', PROGNAME )

            EIDX2 = 0   ! array
            ELAYER = 0  ! array

C.............  If needed for explicit plume rise, allocate and set local array
C               with the mid-point of layers to use as fake stack heights.
            IF ( EXPLFLAG ) THEN

                ALLOCATE( LAYRMID( EMLAYS ), STAT=IOS )
                CALL CHECKMEM( IOS, 'LAYRMID', PROGNAME )

C.................  Method of midpoint depends on vertical structure
                SELECT CASE ( VGTYP ) 
                CASE ( VGHVAL3 ) 

                    DO L = 1, EMLAYS
                        LAYRMID( L ) = VGLVS( L-1 ) + 0.5 *
     &                               ( VGLVS( L ) - VGLVS( L-1 ) )
                    END DO

                CASE DEFAULT
                    WRITE( MESG,94010 ) 'Do not know ' //
     &                     'how to set layer midpoints for layer '//
     &                     'structure type', VGTYP
                    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

                END SELECT

            END IF

C.............  Read stack groups file variables. Read in group ID using x-loc
C               array because SAFE_READ3 is expecting a real array to be
C               passed to it.
            I = PVSDATE
            J = PVSTIME
            CALL SAFE_READ3_I( PVNAME, 'COL' , 1, I, J, GRPCOL )
            CALL SAFE_READ3_I( PVNAME, 'ROW' , 1, I, J, GRPROW )
            CALL SAFE_READ3( PVNAME, 'XLOCA' , 1, I, J, GRPXX )
            CALL SAFE_READ3( PVNAME, 'YLOCA' , 1, I, J, GRPYY )
            CALL SAFE_READ3( PVNAME, 'STKHT' , 1, I, J, GRPHT )
            CALL SAFE_READ3( PVNAME, 'STKDM' , 1, I, J, GRPDM )
            CALL SAFE_READ3( PVNAME, 'STKTK' , 1, I, J, GRPTK )
            CALL SAFE_READ3( PVNAME, 'STKVE' , 1, I, J, GRPVE )

C.............  Read source characteristics from the point source inventory
C               file
            CALL RDINVCHR( 'POINT', PENAME, PSDEV, NPSRC, 1, 'CSOURC' )

C.............  Copy source characteristics from sources to elevated list with
C               only country/state/county, plant, and stack.  For IDA 
C               inventories, the stack ID comes after the point ID, so need to
C               use position of stack ID in source definition to properly
C               build the elevated sources arrays
            DO S = 1, NSRC 

                I    = ELEVSIDX( S )
                IF ( I .EQ. 0 ) CYCLE                
                CFIP = CSOURC( S )( PTBEGL3( 1 ):PTENDL3( 1 ) )
                FCID = CSOURC( S )( PTBEGL3( 2 ):PTENDL3( 2 ) )
                SKID = CSOURC( S )( PTBEGL3(JSTACK):PTENDL3(JSTACK) )
                ESRTIDX( I ) = I
                EIDXSRC( I ) = S
                ECSRCA ( I ) = CFIP // FCID // SKID                

            END DO

C.............  Sort elevated sources
            CALL SORTIC( NMAJOR, ESRTIDX, ECSRCA )

            XMAX = XORIG + XCELL * NCOLS
            YMAX = YORIG + YCELL * NROWS

C.............  Create indices from major-sources list to output list
C.............  Eliminate sources that are not in the domain
C.............  Add sources that have explicit plume rise and require fake
C               stacks to ensure that their emissions get put in the correct
C               model layers.
            K  = 0
            KK = 0
            M  = 0
            PECS = ' '
            DO I = 1, NMAJOR 

                J   = ESRTIDX ( I )
                S   = EIDXSRC ( J )
                ECS = ECSRCA  ( J )

C.................  Find group for this record in stack groups list
                GID = GROUPID( S )
                N = FIND1( GID, NGROUP, GRPGID )

C.................  If group not found, error...
                IF( N .LE. 0 ) THEN
                    EFLAG = .TRUE.
c                    CALL FMTCSRC( ECS, 3, BUFFER, L )
c note: FMTCSRC doesn't work here because MODINFO is not populated
                    WRITE( MESG,94010 )
     &                     'ERROR: Group ID', GID, 'in elevated ' //
     &                     'source ASCII input not in stack ' //
     &                     'groups file for:' // CRLF() // BLANK10 // 
     &                     ECS
                    CALL M3MSG2( MESG )
                    CYCLE
                END IF

                COL  = GRPCOL( N )
                ROW  = GRPROW( N )
                XLOC = GRPXX ( N )
                YLOC = GRPYY ( N )

C.................  Skip records that are outside the grid
                IF ( COL .EQ. 0 .OR. ROW .EQ. 0 ) THEN

c                    CALL FMTCSRC( ECS, 3, BUFFER, L )
                    WRITE( MESG,94010 ) 'NOTE: Source is outside the '//
     &                     'domain boundaries - excluding from UAM ' //
     &                     CRLF() // BLANK10 // 'elevated file for:' //
     &                     CRLF()// BLANK10// ECS
                    CALL M3MESG( MESG )
                    CYCLE

C.................  As per UAM's PTSRCE specifications, eliminate sources
C                   where the lat/lon coordinates are on the domain boundaries
                ELSE IF ( OUTFMT .EQ. 'UAM' .AND.
     &                  ( XLOC   .EQ. XMAX  .OR.
     &                    XLOC   .EQ. XORIG .OR.
     &                    YLOC   .EQ. YMAX  .OR.
     &                    YLOC   .EQ. YORIG       ) ) THEN

c                    CALL FMTCSRC( ECS, 3, BUFFER, L )
                    WRITE( MESG,94010 ) 'NOTE: Source is on the ' //
     &                     'domain boundaries - excluding from UAM ' //
     &                     CRLF() // BLANK10 // 'elevated file for:' //
     &                     CRLF()// BLANK10// ECS
                    CALL M3MESG( MESG )
                    CYCLE

                END IF

C...............  Check to see if source has explicit plume rise and 
C                 needs fake stacks
                IF ( EXPLFLAG ) M = FIND1( S, NHRSRC, ELEVSRC )
                IF ( M .GT. 0 ) THEN
                    NEXPLOOP = EMLAYS
                ELSE
                    NEXPLOOP = 1
                END IF

C...............  Loop through layers for explicit plume srcs (1 for
C                 standard sources)
                DO L = 1, NEXPLOOP
                    
C.........................  Set stack height depending on whether the source
C                           has explicit plume rise or not
                    IF( M .GT. 0 ) THEN
                        HT = LAYRMID( L )
                    ELSE
                        HT = GRPHT ( N )
                    END IF

C.........................  Set stack diameter depending on whether the source
C                           is a PinG source (for UAM-V & CAMx), or whether
C                           it's an explicit plume rise source (REMSAD).
                    DM = GRPDM ( N )
                    IF( ALLOCATED( LPING ) ) THEN
                        IF( LPING( S ) ) DM = -DM
                    END IF
                    IF( M .GT. 0 ) DM = -DM  ! Explicit source

C..................  If FIPS/plant/stack is not the same as the
C                    previous FIP/plant/stack, then store this unique
C                    list as needed for output to UAM format header
                    IF( ECS .NE. PECS ) THEN

                        KK = KK + 1
                        EOUTCSRC( KK ) = ECS

C.........................  Set output stack parameters and coordinates
                        EOUTHT ( KK ) = HT
                        EOUTDM ( KK ) = DM
                        EOUTTK ( KK ) = GRPTK ( N )
                        EOUTVE ( KK ) = GRPVE ( N ) * 3600.    ! m/s -> m/hr
                        EOUTXL ( KK ) = XLOC
                        EOUTYL ( KK ) = YLOC

                    END IF  ! end new FIP/plant/stack

C..................  For PinG sources, may need to overwrite stored diameter;
C                    first source in group may not be PinG, but later sources
C                    can be
                    IF( DM < 0 ) EOUTDM ( KK ) = DM

C..................  Store arrays for writing emissions
                    ELAYER ( I ) = L
                    EIDX2  ( I ) = KK

                END DO      ! end loop over layers or 1

                PECS = ECS

            END DO          ! end of major sources

            NOUT = KK

            IF( EFLAG ) THEN
                MESG = 'Problem processing elevated sources.'
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF

            IF( OUTFMT .EQ. 'UAM' ) THEN

C.................  Get settings from the environment...

C.................  Get type of plume rise calculation from environment
                MESG = 'Vertical method: "PLUMERISE" or "STACKHGT"'
                CALL ENVSTR( 'UAM_VERTICAL_METHOD', MESG, 'STACKHGT', 
     &                       VTYPE, IOS )

C.................  Get UAM file note
                NOTEDEF = 'UAM elev pt emis from ' // 
     &                    PROGNAME
                MESG = 'UAM file note for output file'
                CALL ENVSTR( 'UAM_NOTE', MESG, NOTEDEF , UNOTE, IOS )

C.................  Insert grid name
                L = LEN_TRIM( GRDNM )
                FNOTE  = GRDNM( 1:L ) // ' ' // ADJUSTL( UNOTE )

C.................  Vertical structure information
                EML = EMLAYS
                MESG = 'Number of emission layers'
                IF( EML .EQ. 1 ) 
     &              EML = ENVINT( 'SMK_EMLAYS', MESG, 8, IOS )

                MESG = 'Number of UAM vertical layers'
                NULAYS = ENVINT( 'UAM_LAYERS', MESG,
     &                           EML , IOS )

                NZLOWR = ENVINT( 'UAM4_LAYBELOW', 
     &                           'Layers below diffbreak', 3, IOS )
                NZUPPR = NULAYS - NZLOWR
                NZUPPR = ENVINT( 'UAM4_LAYABOVE', 
     &                           'Layers above diffbreak', NZUPPR, IOS )

                MESG = 'Height of surface layer [meters]'
                HTSUR = ENVREAL( 'UAM4_HTSFC', MESG, 0., IOS )

                MESG= 'Min. ht. of cells b/w sfc and diffbreak [meters]'
                HTLOWR = ENVREAL( 'UAM4_HTLOWR', MESG, 20., IOS )

                MESG= 'Min. ht. of cells b/w diffbreak and top [meters]'
                HTUPPR = ENVREAL( 'UAM4_HTUPPR', MESG, 100., IOS )

C.................  Control packet settings
                MESG = 'Number of parameters for UAM control packet'
                NPARAM = ENVINT( 'UAM_NPARAM', MESG, 20, IOS )

                MESG = 'PTSRCE output file unit number'
                PDEVOUT = ENVINT( 'UAM_PTSRCE_OUTUNIT', MESG, 20, IOS )

                MESG = 'Print output grid, 0=no, 1=yes'
                GSWITCH = ENVINT( 'UAM_PRINT_OUTGRD', MESG, 0, IOS )

                MESG = 'Print units table, 0=no, 1=yes'
                USWITCH = ENVINT( 'UAM_PRINT_UNITS', MESG, 0, IOS )

                MESG = 'Print source locations table, 0=no, 1=yes'
                LSWITCH = ENVINT( 'UAM_PRINT_LOCATIONS', MESG, 0, IOS )

                MESG = 'Print methods table, 0=no, 1=yes'
                MSWITCH = ENVINT( 'UAM_PRINT_METHODS', MESG, 0, IOS )

                MESG = 'Print source values table, 0=no, 1=yes'
                VSWITCH = ENVINT( 'UAM_PRINT_VALUES', MESG, 0, IOS )

                MESG = 'Print vertical methods table, 0=no, 1=yes'
                ESWITCH = ENVINT( 'UAM_PRINT_VERTMETH', MESG, 0, IOS )

                MESG = 'DIFFBREAK file unit number'
                DDEVOUT = ENVINT( 'UAM4_DIFFBREAK_UNIT', MESG, 0, IOS )

                MESG = 'REGIONTOP file unit number'
                RDEVOUT = ENVINT( 'UAM4_REGIONTOP_UNIT', MESG, 0, IOS )

C.................  If PLUMERISE method has been selected, get file unit
C                   assignments
                IF( VTYPE .EQ. 'PLUMERISE' ) THEN

                    MESG = 'TEMPERATUR file unit number'
                    TDEVOUT= ENVINT( 'UAM4_TEMPERATUR_UNIT', 
     &                               MESG, 14, IOS ) 

                    MESG = 'METSCALARS file unit number'
                    MDEVOUT= ENVINT( 'UAM4_METSCALARS_UNIT', 
     &                               MESG, 15, IOS )

                    MESG = 'WIND file unit number'
                    WDEVOUT = ENVINT( 'UAM4_WIND_UNIT', MESG, 16, IOS )

                END IF

C.................  Get output units from the environment
                BUFFER = 'Units for output gridded emissions'
                CALL ENVSTR( 'MRG_GRDOUT_UNIT', BUFFER,
     &                         ' ', GRDENV, IOS)

C.................  Write header information to the file
                WRITE( EVDEV,93010 ) 'CONTROL', GRDENV
                WRITE( EVDEV,93000 ) 'PTSOURCE'
                WRITE( EVDEV,93000 ) FNOTE
    
C.................  Insert crrect number of sources into 1st line of ESTRING
                WRITE( EVDEV,93015 ) NMSPC, 0, NOUT, 1, NPARAM
                WRITE( EVDEV,93015 ) PDEVOUT, 0, GSWITCH
                WRITE( EVDEV,93015 ) USWITCH, LSWITCH, 0, MSWITCH,
     &                               VSWITCH
                WRITE( EVDEV,93015 ) 1, 0, ESWITCH
                WRITE( EVDEV,93015 ) DDEVOUT, RDEVOUT, 0, TDEVOUT,
     &                               MDEVOUT, WDEVOUT

C.................  Write species names
                DO I = 1, NMSPC
                    SPCNAM = EMNAM( I )( 1:LEN( SPCNAM ) )

                    IF( SPCNAM .NE. EMNAM( I ) ) THEN
                        L = LEN_TRIM( EMNAM( I ) ) 
                        MESG = 'WARNING: Species name "' // 
     &                         EMNAM( I )( 1:L ) // 
     &                         '" is being truncated to "' // SPCNAM //
     &                         CRLF() // BLANK10 // '" on output to ' //
     &                         'UAM elevated point sources file'
                        CALL M3MSG2( MESG )

                    END IF
 
                    WRITE( EVDEV, 93020, IOSTAT=IOS ) SPCNAM

                END DO

C.................  Output the beginning and ending dates and times
                IBD = REMOVE_4DIGIT_YEAR( SDATE )
                IED = REMOVE_4DIGIT_YEAR( EDATE )
                IBT = STIME / 100
                IET = ETIME / 100
                WRITE( EVDEV, 93030 ) IBD, IBT, IED, IET

                WRITE( EVDEV, 93000 ) 'END'

C.................  Region packet:
                WRITE( EVDEV, 93000 ) 'REGION'
                WRITE( EVDEV, 93040 ) 0., 0., INT( P_ALP )
                WRITE( EVDEV, 93040 ) XORIG, YORIG

                IF( GDTYP .EQ. LATGRD3 ) THEN           ! Print out for Lat-Lon
                    WRITE( EVDEV, 93045 ) XCELL, YCELL
                ELSE                                    ! other coordinate sys
                    WRITE( EVDEV, 93040 ) XCELL, YCELL
                END IF

                WRITE( EVDEV, 93050 ) NCOLS, NROWS, NULAYS
                WRITE( EVDEV, 93060 ) NZLOWR, NZUPPR, HTSUR, 
     &                                HTLOWR, HTUPPR 
                WRITE( EVDEV, 93000 ) 'END' 

C.................  Point Sources packet:

                WRITE( EVDEV, 93000 ) 'POINT SOURCES'

                DO I = 1, NOUT

                    ECS  = EOUTCSRC( I )
                    CFIP = ECS( 1:FIPLEN3 )
                    FCID = ADJUSTL( ECS( PTBEGL3( 2 ):PTENDL3( 2 ) ) )
                    SKID = ADJUSTL( ECS( PTBEGL3( 3 ):PTENDL3( 3 ) ) )

                    L = LEN_TRIM( FCID )
                    IF( L .GT. 10 ) THEN
                        MESG = 'WARNING: Plant name truncated from "' //
     &                         FCID( 1:L ) // '" to "' // FCID( 1:10 )//
     &                         '" on output to ASCII elevated file'
                        CALL M3MESG( MESG )
                    END IF

                    L = LEN_TRIM( SKID )
                    IF( L .GT. 10 ) THEN
                        MESG = 'WARNING: Stack name truncated from "' //
     &                         SKID( 1:L ) // '" to "' // SKID( 1:10 )//
     &                         '" on output to ASCII elevated file'
                        CALL M3MESG( MESG )
                    END IF

                    FIP  = STR2INT( CFIP )
                    FCID = ADJUSTR( FCID( 1:10 ) )
                    SKID = ADJUSTR( SKID( 1:10 ) )

C.....................  Output source informat in a different format if 
C                       the grid is lat-lon or not
                    IF( GDTYP .EQ. LATGRD3 ) THEN
                        WRITE( EVDEV, 93065 )
     &                      I, 'STD       ', EOUTXL( I ), EOUTYL( I ),
     &                      FCID( 1:10 ), SKID( 1:10 ), FIP

                    ELSE 
                        WRITE( EVDEV, 93070, IOSTAT=IOS )
     &                      I, 'STD       ', EOUTXL( I ), EOUTYL( I ),
     &                      FCID( 1:10 ), SKID( 1:10 ), FIP

                    END IF

                    WRITE( EVDEV, 93080, IOSTAT=IOS )
     &                     EOUTHT ( I ), EOUTDM ( I ), 
     &                     EOUTTK ( I ), EOUTVE ( I )

                END DO

                WRITE( EVDEV, 93000 ) 'END' 

            END IF    ! End of output format for UAM file

            FIRSTIME = .FALSE.

        END IF        ! end if firstime routine is called

C.........  Write out date/time-specific data...
C.........  Data previously merged in MRGELEV
        IF( JTIME .NE. PTIME .AND. OUTFMT .EQ. 'UAM' ) THEN

C.............  Time Interval packet:
C.............  Set required time parameters
            JDATEP1 = JDATE
            JTIMEP1 = JTIME
            CALL NEXTIME( JDATEP1, JTIMEP1, TSTEP )

            IBT = JTIME   / 100
            IET = JTIMEP1 / 100            
            IBD = REMOVE_4DIGIT_YEAR( JDATE )
            IED = REMOVE_4DIGIT_YEAR( JDATEP1 )

            WRITE( EVDEV, 93000 ) 'TIME INTERVAL'
            WRITE( EVDEV, 93050 ) IBD, IBT, IED, IET

C.............  Method packet: 
C.............  Provide same output as PSTPNT

            WRITE( EVDEV, 93000 ) 'METHOD'
            WRITE( EVDEV, 93000 ) 
     &            'STD       ALL       EMVALUES  0.        50000.'
            WRITE( EVDEV, 93000 ) 'END' 

C.............  Vertical Method packet:
C.............  Provide same output as PSTPNT with variable plume-rise method

            WRITE( EVDEV, 93000 ) 'VERTICAL METHOD'
            WRITE( EVDEV, 93000 ) 
     &             'STD       ALL       ' // VTYPE // ' 0.       10000.'
            WRITE( EVDEV, 93000 ) 'END' 

        END IF    ! End UAM/CAMx elevated point sources format

C.........  But first determine how many explicit sources
C           there are for the current hour.
        DO N = 1, NHRSRC
            IF ( INDXH( N ) .EQ. 0 ) EXIT
        END DO
        NH = N - 1

C.........  Write out emissions data for UAM/CAMx elevated point sources

        IF( OUTFMT .EQ. 'UAM' ) THEN

C.............  Emissions Values packet.  Only write this if the time is 
C               new, otherwise, it will write each time routine is called
C               for the next species.
            IF( JTIME .NE. PTIME ) THEN
                WRITE( EVDEV, 93000 ) 'EMISSIONS VALUES'
                WRITE( EVDEV, 93000 ) 'ALL       ALL            0.000'
            END IF

            LN = 0
            LM = 0
            LL = 0
            M  = 0
            ESUM = 0.  ! array
            PECS = ' '
            DO I = 1, NMAJOR

                J   = ESRTIDX ( I )  ! sort for UAM output index
                L   = ELAYER  ( I )  ! layer, 1 or 0 if outside grid
                N   = EIDX2   ( I )  ! UAM elev src no. or 0 if outside grd

                IF ( N .LE. 0 ) CYCLE   ! skip because not in grid

                S   = EIDXSRC( J )                
                ECS = ECSRCA ( J )

                IF ( EXPLFLAG ) M = FIND1( S, NHRSRC, ELEVSRC )

C NOTE: This algorithm will not work if there are explicit elevated sources
C    N: that have duplicates.

C.................  If FIP/plant/stack is different from last time, 
C                   write emissions for previous elevated output src 
C                   and initialize
                IF ( PECS .NE. ECS .OR. L .NE. LL ) THEN

C.....................  Write emissions for standard output srcs
                    IF( LM .LE. 0 ) THEN

C......................  Only write emissions if not zero
                        IF( ESUM( 1 ) .NE. 0. ) THEN 
                            CALL GET_ESUM_FORMAT( VNAME, ESUM(1), EFMT )
                            WRITE( EVDEV, EFMT ) LN, VNAME, ESUM(1)
                        END IF
                        ESUM( 1 ) = 0.

C.....................  Write emissions for explicit source records
                    ELSE

                        IF( ESUM( LL ) .NE. 0. ) THEN 
                            CALL GET_ESUM_FORMAT( VNAME,ESUM(LL),EFMT )
                            WRITE( EVDEV,EFMT ) LN, VNAME, ESUM( LL )
                        END IF
                        ESUM( L ) = 0.

                    END IF

C.....................  Initialize for standard FIP/plant/stack
                    IF( M .LE. 0 ) THEN
                        ESUM( 1 ) = ELEVEMIS( J )

C.....................  Initialize for explicit sources
                    ELSE
                        ESUM( L ) = ELEVEMIS( J ) * LFRAC( M,L )

                    END IF

C.................  If FIP/plant/stack is the same as the previous 
C                   iteration, sum emissions
                ELSE

C.....................  For standard sources
                    IF( M .LE. 0 ) THEN

                        ESUM( 1 ) = ESUM( 1 ) + ELEVEMIS( J )

C.....................  For explicit sources
                    ELSE
                        ESUM(L) = ESUM(L) + ELEVEMIS(J) * LFRAC(M,L)

                    END IF

                END IF

                LN = N       ! index to grouped UAM source number
                LL = L       ! layer number or 1
                LM = M
                PECS = ECS

            END DO

C.............  Output for final record(s)...
C.............  For standard sources
            IF( M .LE. 0 ) THEN

                IF( ESUM( 1 ) .NE. 0. ) THEN 
                    CALL GET_ESUM_FORMAT( VNAME, ESUM(1), EFMT )
                    WRITE( EVDEV, EFMT ) NOUT, VNAME, ESUM(1)
                END IF

C.............  For explicit sources
            ELSE
                IF( ESUM( L ) .NE. 0. ) THEN
                    CALL GET_ESUM_FORMAT( VNAME, ESUM( L ), EFMT )
                    WRITE( EVDEV,EFMT ) NOUT, VNAME, ESUM( L )
                END IF

            END IF

C.............  After last species, write packet end fields
            IF( VNAME .EQ. EMNAM( NMSPC ) ) THEN
                WRITE( EVDEV, 93000 ) 'END' 
                WRITE( EVDEV, 93000 ) 'ENDTIME'
            END IF

        END IF   ! End UAM/CAMx elevated point sources format

        PTIME = JTIME

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Formatted file I/O formats............ 93xxx

93000   FORMAT( A )

93010   FORMAT( A7, 3X, A )

93015   FORMAT( 6I10 )

93020   FORMAT( A10 )

93030   FORMAT( 4I10 )

93040   FORMAT( 2F10.0, I10 )

93045   FORMAT( 2F10.7 )

93050   FORMAT( 4I10 )

93060   FORMAT( 2I10, 3F10.0 )

93065   FORMAT( I10, A10, F10.5, F10.5, 2A10, I10.5 )

93070   FORMAT( I10, A10, F10.0, F10.0, 2A10, I10.5 )

93080   FORMAT( F10.1, F10.2, F10.1, F10.0 )

93090   FORMAT( I10, A10, F10.3 )

C...........   Internal buffering formats.............94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

94300   FORMAT( A, I10, A )

C*****************  INTERNAL SUBPROGRAMS  ******************************

        CONTAINS

C.............  This internal subprogram tries to read a variable from an
C               I/O API file, and aborts if not successful.
            SUBROUTINE SAFE_READ3( FILNAM, VARNAM, LAYER,
     &                             JDATE, JTIME, XBUF     )

C.............  Subprogram arguments
            CHARACTER(*) FILNAM    ! logical file name
            CHARACTER(*) VARNAM    ! variable name
            INTEGER      LAYER     ! layer number (or ALLAYS3)
            INTEGER      JDATE     ! Julian date
            INTEGER      JTIME     ! time
            REAL         XBUF( * ) ! read buffer

C.............  Local variables
            INTEGER L1, L2

C----------------------------------------------------------------------

            IF ( .NOT. READ3( FILNAM, VARNAM, LAYER,
     &                        JDATE, JTIME, XBUF ) ) THEN

                L1 = LEN_TRIM( VARNAM )
                L2 = LEN_TRIM( FILNAM )
                MESG = 'Could not read "' // VARNAM( 1:L1 ) //
     &                 '" from file "' // FILNAM( 1:L2 ) // '."'
                CALL M3EXIT( PROGNAME, JDATE, JTIME, MESG, 2 )

            END IF

            END SUBROUTINE SAFE_READ3

C----------------------------------------------------------------------
C----------------------------------------------------------------------

C.............  This internal subprogram tries to read a variable from an
C               I/O API file, and aborts if not successful.
            SUBROUTINE SAFE_READ3_I( FILNAM, VARNAM, LAYER,
     &                              JDATE, JTIME, IBUF     )

C.............  Subprogram arguments
            CHARACTER(*) FILNAM    ! logical file name
            CHARACTER(*) VARNAM    ! variable name
            INTEGER      LAYER     ! layer number (or ALLAYS3)
            INTEGER      JDATE     ! Julian date
            INTEGER      JTIME     ! time
            INTEGER      IBUF( * ) ! read buffer

C.............  Local variables
            INTEGER L1, L2

C----------------------------------------------------------------------

            IF ( .NOT. READ3( FILNAM, VARNAM, LAYER,
     &                        JDATE, JTIME, IBUF ) ) THEN

                L1 = LEN_TRIM( VARNAM )
                L2 = LEN_TRIM( FILNAM )
                MESG = 'Could not read "' // VARNAM( 1:L1 ) //
     &                 '" from file "' // FILNAM( 1:L2 ) // '."'
                CALL M3EXIT( PROGNAME, JDATE, JTIME, MESG, 2 )

            END IF

            END SUBROUTINE SAFE_READ3_I

C----------------------------------------------------------------------
C----------------------------------------------------------------------

C.............  This internal subprogram removes the first two digits of the
C               year from a 7-digit Julian date
            INTEGER FUNCTION REMOVE_4DIGIT_YEAR( JDATE )

C.............  Subroutine arguments
            INTEGER, INTENT (IN) :: JDATE

C.............  Local variables
            INTEGER    YRREMOVE 

C----------------------------------------------------------------------

            YRREMOVE = ( JDATE / 100000 ) * 100000  ! integer math
            REMOVE_4DIGIT_YEAR = JDATE - YRREMOVE

            RETURN

            END FUNCTION REMOVE_4DIGIT_YEAR

C----------------------------------------------------------------------
C----------------------------------------------------------------------

C.............  This internal subprogram determines the format to output
C               emission values
            SUBROUTINE GET_ESUM_FORMAT( VBUF, VAL, FMT )

C.............  Subroutine arguments
            CHARACTER(*), INTENT (IN) :: VBUF
            REAL        , INTENT (IN) :: VAL
            CHARACTER(*), INTENT(OUT) :: FMT

C----------------------------------------------------------------------

            IOS = 0

C.............  Value is too large for 
            IF( VAL .GT. 999999999. ) THEN
                FMT = '( I10, A10, E10.3 )'

                L = LEN_TRIM( VBUF )
                WRITE( MESG,94020 ) 'WARNING: "' // VBUF( 1:L ) // 
     &            '"Emissions value of', ESUM, CRLF()// BLANK10// 
     &            '" is too large for file format, so writing ' //
     &            'in scientific notation for source:',
     &            CRLF() // BLANK10 // ECSRCA( J )

                CALL M3MESG( MESG )

            ELSE IF( VAL .GT. 99999999. ) THEN
                FMT = '( I10, A10, F10.0 )'

            ELSE IF( VAL .GT. 9999999. ) THEN
                FMT = '( I10, A10, F10.1 )'

            ELSE IF( VAL .GT. 999999. ) THEN
                FMT = '( I10, A10, F10.2 )'

            ELSE IF( VAL .GT. 0. .AND. VAL .LT. 1. ) THEN
                FMT = '( I10, A10, E10.4 )'

            ELSE
                FMT = '( I10, A10, F10.3 )'

            END IF

            RETURN

94020       FORMAT( 10( A, :, E12.5, :, 1X ) )

            END SUBROUTINE GET_ESUM_FORMAT

        END SUBROUTINE WMRGELEV
