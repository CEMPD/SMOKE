
        SUBROUTINE WMRGELEV( VNAME, NMAJOR, JDATE, JTIME )

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
C COPYRIGHT (C) 2000, MCNC--North Carolina Supercomputing Center
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

C.........  MODULES for public variables
C.........  This module contains the major data structure and control flags
        USE MODMERGE

C...........   This module is the source inventory arrays
        USE MODSOURC

C.........  This module contains arrays for plume-in-grid and major sources
        USE MODELEV

C.........  This module contains the global variables for the 3-d grid
        USE MODGRID

        IMPLICIT NONE

C.........  INCLUDES:
        
        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
        INCLUDE 'FDESC3.EXT'    !  I/O API file desc. data structures

C.........  EXTERNAL FUNCTIONS and their descriptions:
        
        CHARACTER*2     CRLF
        INTEGER         ENVINT
        REAL            ENVREAL
        INTEGER         FIND1
        INTEGER         INDEX1
        INTEGER         STR2INT

        EXTERNAL        CRLF, ENVINT, ENVREAL, FIND1, INDEX1, STR2INT

C.........  SUBROUTINE ARGUMENTS
        CHARACTER(*), INTENT (IN) :: VNAME   ! variable name to output
        INTEGER     , INTENT (IN) :: NMAJOR  ! no. elevated sources
        INTEGER     , INTENT (IN) :: JDATE   ! Julian date to output (YYYYDDD)
        INTEGER     , INTENT (IN) :: JTIME   ! time to output (HHMMSS)

C.........  Local parameters
        INTEGER, PARAMETER :: STKLEN3 = FIPLEN3 + PLTLEN3 + CHRLEN3

C.........  Variables allocated by module settings...
        INTEGER, ALLOCATABLE, SAVE :: ESRTIDX( : ) ! major srcs sorting index
        INTEGER, ALLOCATABLE, SAVE :: EGRPIDX( : ) ! group index

        REAL   , ALLOCATABLE :: EOUTHT ( : ) ! output stack heights [m]
        REAL   , ALLOCATABLE :: EOUTDM ( : ) ! output stack diameters [m]
        REAL   , ALLOCATABLE :: EOUTTK ( : ) ! output exit tmprs [K]
        REAL   , ALLOCATABLE :: EOUTVE ( : ) ! output exit velocities [m/hr]
        REAL   , ALLOCATABLE :: EOUTXL ( : ) ! output x-loc [units of grid]
        REAL   , ALLOCATABLE :: EOUTYL ( : ) ! output y-loc [units of grid]

        CHARACTER(LEN=STKLEN3), ALLOCATABLE :: ECSRCA   ( : ) ! FIPS//plt//stk
        CHARACTER(LEN=STKLEN3), ALLOCATABLE :: EOUTCSRC ( : ) ! FIPS//plt//stk

C.........  UAM-format specific variables
        INTEGER         DDEVOUT             ! diffbreak file unit number
        INTEGER         ESWITCH             ! 0=no, 1=yes - print vert mthds tbl
        INTEGER         GSWITCH             ! 0=no, 1=yes - print output grid
        INTEGER         LSWITCH             ! 0=no, 1=yes - print src locs tbl
        INTEGER         MDEVOUT             ! metscalars file unit number
        INTEGER         MSWITCH             ! 0=no, 1=yes - print methods tbl
        INTEGER         NPARAM              ! no. paramaters for control pkt
        INTEGER         NULAYS              ! no. UAM model layers
        INTEGER         NZLOWR              ! no. layers below Diffbreak height
        INTEGER         NZUPPR              ! no. layers above Diffbreak height
        INTEGER         PDEVOUT             ! PTSRCE output unit no.
        INTEGER         RDEVOUT             ! regiontop file unit number
        INTEGER         TDEVOUT             ! temperatur file unit number
        INTEGER         USWITCH             ! 0=no, 1=yes - print units table
        INTEGER         VSWITCH             ! 0=no, 1=yes - print values tbl
        INTEGER         WDEVOUT             ! wind file unit number

        REAL            HTSUR               ! height of surface layer [m]
        REAL            HTLOWR              ! min cell ht b/w sfc and diffbr [m]
        REAL            HTUPPR              ! min cell ht b/w diffbr and top [m]

        CHARACTER*10        :: SPCNAM       ! UAM-format species name
        CHARACTER*10 , SAVE :: VTYPE        ! User-spec vertical method type
        CHARACTER*44           NOTEDEF      ! Default note
        CHARACTER*44           UNOTE        ! UAM file note from env variable
        CHARACTER*60        :: FNOTE        ! UAM file header note

C.........  Other local variables
        INTEGER          I, J, K, L, LK, N, S 

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
        INTEGER          NOUT                ! number output stacks
        INTEGER, SAVE :: PTIME = -9          ! previous call's time (HHMMSS)
        INTEGER          ROW                 ! tmp row

        REAL             ESUM                ! emissions data sum
        REAL             XLOC                ! tmp x location
        REAL             XMAX                ! rightmost x location
        REAL             YLOC                ! tmp y location
        REAL             YMAX                ! topmost y location

        LOGICAL       :: EFLAG    = .FALSE.  ! true: error occurred
        LOGICAL       :: FIRSTIME = .TRUE.   ! true: first time routine called

        CHARACTER*16           OUTFMT       ! output format for elevated ASCII
        CHARACTER*200          BUFFER       ! source chars buffer
        CHARACTER*300          MESG         ! message buffer

        CHARACTER(LEN=FIPLEN3) CFIP         ! tmp country/state/county code
        CHARACTER(LEN=PLTLEN3) FCID         ! tmp facility ID
        CHARACTER(LEN=CHRLEN3) SKID         ! tmp stack ID
        CHARACTER(LEN=STKLEN3) ECS          ! stack elevated source chars
        CHARACTER(LEN=STKLEN3) PECS         ! tmp previous ECS

        CHARACTER*16 :: PROGNAME = 'WMRGELEV' ! program name

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

C.............  Allocate memory for local elevated sources arrays
            ALLOCATE( ESRTIDX( NMAJOR ), STAT=IOS )
            CALL CHECKMEM( IOS, 'ESRTIDX', PROGNAME )
            ALLOCATE( ECSRCA( NMAJOR ), STAT=IOS )
            CALL CHECKMEM( IOS, 'ECSRCA', PROGNAME )
            ALLOCATE( EOUTCSRC( NMAJOR ), STAT=IOS )
            CALL CHECKMEM( IOS, 'EOUTCSRC', PROGNAME )
            ALLOCATE( EGRPIDX( NMAJOR ), STAT=IOS )
            CALL CHECKMEM( IOS, 'EGRPIDX', PROGNAME )
            ALLOCATE( EOUTHT( NMAJOR ), STAT=IOS )
            CALL CHECKMEM( IOS, 'EOUTHT', PROGNAME )
            ALLOCATE( EOUTDM( NMAJOR ), STAT=IOS )
            CALL CHECKMEM( IOS, 'EOUTDM', PROGNAME )
            ALLOCATE( EOUTTK( NMAJOR ), STAT=IOS )
            CALL CHECKMEM( IOS, 'EOUTTK', PROGNAME )
            ALLOCATE( EOUTVE( NMAJOR ), STAT=IOS )
            CALL CHECKMEM( IOS, 'EOUTVE', PROGNAME )
            ALLOCATE( EOUTXL( NMAJOR ), STAT=IOS )
            CALL CHECKMEM( IOS, 'EOUTXL', PROGNAME )
            ALLOCATE( EOUTYL( NMAJOR ), STAT=IOS )
            CALL CHECKMEM( IOS, 'EOUTYL', PROGNAME )

C.............  Initialize local elevated sources arrays
            EGRPIDX = 0   ! array

C.............  Read stack groups file variables. Read in group ID using x-loc
C               array because SAFE_READ3 is expecting a real array to be
C               passed to it.
            I = PVSDATE
            J = PVSTIME
            CALL SAFE_READ3_I( PVNAME, 'COL' , 1, I, J, GRPCOL )
            CALL SAFE_READ3_I( PVNAME, 'ROW' , 1, I, J, GRPROW )
            CALL SAFE_READ3( PVNAME, 'XLOCA' , 1, I, J, GRPXL )
            CALL SAFE_READ3( PVNAME, 'YLOCA' , 1, I, J, GRPYL )
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
            DO I = 1, NMAJOR 

                S    = ELEVSIDX( I )                
                CFIP = CSOURC( S )( PTBEGL3( 1 ):PTENDL3( 1 ) )
                FCID = CSOURC( S )( PTBEGL3( 2 ):PTENDL3( 2 ) )
                SKID = CSOURC( S )( PTBEGL3(JSTACK):PTENDL3(JSTACK) )
                ESRTIDX( I ) = I
                ECSRCA ( I ) = CFIP // FCID // SKID                

            END DO

C.............  Sort elevated sources
            CALL SORTIC( NMAJOR, ESRTIDX, ECSRCA )

            XMAX = XORIG + XCELL * NCOLS
            YMAX = YORIG + YCELL * NROWS

C.............  Create indices from major-sources list to output list
C.............  Eliminate sources that are not in the domain
            PECS = ' '
            K = 0
            DO I = 1, NMAJOR

                J   = ESRTIDX ( I )
                S   = ELEVSIDX( J )                
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
                XLOC = GRPXL ( N )
                YLOC = GRPYL ( N )

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

                IF( ECS .NE. PECS ) THEN
                    K = K + 1

                    EOUTCSRC( K ) = ECS

C.....................  Set stack parameters and coordinates
                    EOUTHT ( K ) = GRPHT ( N )
                    EOUTDM ( K ) = GRPDM ( N )
                    EOUTTK ( K ) = GRPTK ( N )
                    EOUTVE ( K ) = GRPVE ( N ) * 3600.    ! m/s -> m/hr
                    EOUTXL ( K ) = XLOC
                    EOUTYL ( K ) = YLOC

                    PECS = ECS


                END IF

                EGRPIDX( J ) = K

            END DO

            NOUT = K

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
                NOTEDEF = 'UAM elevated point emissions from ' // 
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

C.................  Write header information to the file
        	WRITE( EVDEV,93000 ) 'CONTROL'
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

        END IF

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

C.........  Write out emissions data for UAM/CAMx elevated point sources

        IF( OUTFMT .EQ. 'UAM' ) THEN

C.............  Emissions Values packet.  Only write this if the time is 
C               new, otherwise, it will write each time routine is called
C               for the next species.
            IF( JTIME .NE. PTIME ) THEN
        	WRITE( EVDEV, 93000 ) 'EMISSIONS VALUES'
        	WRITE( EVDEV, 93000 ) 'ALL       ALL            0.000'
            END IF

            LK = -9
            ESUM = 0.
            DO I = 1, NMAJOR

                J   = ESRTIDX ( I )
                K   = EGRPIDX ( J )

C.................  If K is zero, then source is being skipped
                IF( K .EQ. 0 ) THEN
                    CYCLE

C.................  If current record is new, print out previous value and
C                   initialize emissions
                ELSE IF( K .NE. LK ) THEN
                    
                    IF( ESUM .NE. 0. ) 
     &                  WRITE( EVDEV, 93090 ) LK, VNAME, ESUM
                    ESUM = ELEVEMIS( I )

C.................  If current record is for the same plant and stack, sum
C                   emissions
                ELSE IF( K .EQ. LK ) THEN
                    ESUM = ESUM + ELEVEMIS( I )

                END IF

                LK = K

            END DO

C.............  Output for final record
            IF( ESUM .NE. 0. ) WRITE( EVDEV, 93090 ) K, VNAME, ESUM

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

93015   FORMAT( 6I10 )

93020   FORMAT( A10 )

93030   FORMAT( 4I10 )

93040   FORMAT(	2F10.0, I10 )

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

        END SUBROUTINE WMRGELEV
