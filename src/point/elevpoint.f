
        PROGRAM ELEVPOINT

C***********************************************************************
C  program body starts at line 
C
C  DESCRIPTION:
C       Identifies sources as elevated (major or plume-in-grid) or minor.
C       Major sources will get plume rise and minor sources will not, however,
C       the major/minor distinction is not required because SMOKE will compute
C       layer fractions for all sources efficiently when needed.  If desired,
C       the program can use an analytical computation (the PLUMRIS routine)
C       along with a cutoff height to determine the major sources.
C
C       NOTE - The initial version uses input files to identify the sources
C       or can use the PLUMERIS and cutoff height.  Future versions will
C       permit more flexible run-time identification of major and plume-in-grid
C       sources based on source characteristics and hourly emissions.
C       
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C       Copied from elevpoint.F 4.2 by M Houyoux
C
C************************************************************************
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
C Last updated: %G 
C  
C***********************************************************************

C...........   MODULES for public variables
C...........   This module is the source inventory arrays
        USE MODSOURC

C.........  This module contains arrays for plume-in-grid and major sources
        USE MODELEV

C.........  This module contains the information about the source category
        USE MODINFO

        IMPLICIT NONE

C...........   INCLUDES:
        
        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
        INCLUDE 'CONST3.EXT'    !  physical and mathematical constants
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
        INCLUDE 'FDESC3.EXT'    !  I/O API file description data structures.

C...........   EXTERNAL FUNCTIONS and their descriptions:

        CHARACTER*2     CRLF
        LOGICAL         DSCM3GRD
        REAL            ENVREAL
        LOGICAL         ENVYN
        INTEGER         FIND1
        INTEGER         FINDC
        INTEGER         GETFLINE
        REAL            PLUMRIS
        INTEGER         PROMPTFFILE
        CHARACTER*16    PROMPTMFILE

        EXTERNAL        CRLF, DSCM3GRD, ENVREAL, ENVYN, FIND1, FINDC, 
     &                  GETFLINE, PLUMRIS, PROMPTFFILE, PROMPTMFILE

C...........  LOCAL PARAMETERS and their descriptions:
        CHARACTER*50, PARAMETER :: SCCSW = '@(#)$Id$'

C...........   Indicator for which public inventory arrays need to be read
        INTEGER               , PARAMETER :: NINVARR = 9
        CHARACTER(LEN=IOVLEN3), PARAMETER :: IVARNAMS( NINVARR ) = 
     &                                 ( / 'IFIP           '
     &                                   , 'XLOCA          '
     &                                   , 'YLOCA          '
     &                                   , 'STKHT          '
     &                                   , 'STKDM          '
     &                                   , 'STKTK          '
     &                                   , 'STKVE          '
     &                                   , 'CSCC           '
     &                                   , 'CSOURC         ' / )

C...........   Allocateable arrays for using GENPTCEL routine to get grid-cell
C              numbers based on current projection
        INTEGER, ALLOCATABLE :: INDX ( : )  ! sorting index (unused)
        INTEGER, ALLOCATABLE :: GN   ( : )  ! cell numbers
        INTEGER, ALLOCATABLE :: SN   ( : )  ! stack group pos in list (unused)
        INTEGER, ALLOCATABLE :: NX   ( : )  ! no. stack groups per cell (unused)

C...........   Allocatable arrays for reading in stack splits def'n file
        INTEGER, ALLOCATABLE :: SPTINDX ( : )  ! sorting index
        INTEGER, ALLOCATABLE :: SPTGIDA( : )  ! unsorted stack group ID

        REAL   , ALLOCATABLE :: SPTLAT ( : )  ! unsorted splits file latitude
        REAL   , ALLOCATABLE :: SPTLON ( : )  ! unsorted splits file longitude
        REAL   , ALLOCATABLE :: SPTDM  ( : )  ! unsorted splits file diameter 
        REAL   , ALLOCATABLE :: SPTHT  ( : )  ! unsorted splits file height   
        REAL   , ALLOCATABLE :: SPTTK  ( : )  ! unsorted splits file tmpr     
        REAL   , ALLOCATABLE :: SPTVE  ( : )  ! unsorted splits file velocity
        REAL   , ALLOCATABLE :: SPTFL  ( : )  ! unsorted splits file flow

        LOGICAL, ALLOCATABLE :: SPTMMSA( : )  ! true: Major stack (unsorted)
        LOGICAL, ALLOCATABLE :: SPTMPSA( : )  ! true: PinG stack (unsorted)
        LOGICAL, ALLOCATABLE :: FOUND  ( : )  ! true: entry found in inven

        CHARACTER(LEN=ALLLEN3), ALLOCATABLE :: SPTCSRCA( : ) ! src info (unsrt)
        CHARACTER(LEN=ALLLEN3), ALLOCATABLE :: SPTCSRC ( : ) ! src info (sorted)

C...........   File units and logical/physical names
        INTEGER         GDEV    !  stack groups file
        INTEGER         LDEV    !  log-device
        INTEGER         PDEV    !  for output major/mepse src ID file
        INTEGER         SDEV    !  ASCII part of inventory unit no.
        INTEGER         TDEV    !  stack splits file

        CHARACTER*16    ANAME   !  logical name for ASCII inventory input file
        CHARACTER*16    ENAME   !  logical name for i/o api inventory input file
        CHARACTER*16    MNAME   !  plume-in-grid srcs stack groups output file

C...........   Other local variables
        INTEGER         I, J, K, S, L, L2, N      ! indices and counters

        INTEGER         COID          ! tmp country ID
        INTEGER         COL           ! tmp column number
        INTEGER         CYID          ! tmp county ID
        INTEGER         ENLEN         ! inventory file name length
        INTEGER         FIP           ! tmp FIPS code
        INTEGER         GID           ! tmp group ID
        INTEGER         IOS           ! i/o status
        INTEGER         IOSCUT        ! i/o status for cutoff E.V.
        INTEGER         IREC          ! record counter
        INTEGER         NCOLS         ! no. grid column
        INTEGER      :: NEXCLD = 0    ! no. stack groups exlcuded from the grid
        INTEGER         NGRID         ! no. grid cells
        INTEGER         NROWS         ! no. grid rows
        INTEGER         NGLINES       ! no. lines in stack group file
        INTEGER      :: NMAJOR = 0    ! no. major sources
        INTEGER         NPG           ! tmp number per group
        INTEGER      :: NPING  = 0    ! no. plume-in-grid sources
        INTEGER         NSLINES       ! no. lines in stack splits file
        INTEGER         NSTEPS        ! no. time steps
        INTEGER         MS            ! tmp src ID for major sources
        INTEGER         PS            ! tmp src ID for plume in grid sources
        INTEGER         ROW           ! tmp row number
        INTEGER         STID          ! tmp state ID
        INTEGER         SDATE         ! Julian start date
        INTEGER         STIME         ! start time
        INTEGER         TZONE         ! output time zone

        REAL            CUTOFF        ! plume rise cutoff for elev pts
        REAL            DM            ! tmp inside stack diameter [m]
        REAL            FL            ! tmp stack exit flow rate [m^3/s]
        REAL            HT            ! tmp inside stack diameter [m]
        REAL            LAT           ! tmp latitude [degrees]
        REAL            LON           ! tmp longitude [degrees]
        REAL            MINDM         ! min stack group diam
        REAL            MINFL         ! min stack group flow
        REAL            MINHT         ! min stack group height
        REAL            MINTK         ! min stack group temperature
        REAL            MINVE         ! min stack group velocity
        REAL            RISE          ! calculated plume rise
        REAL            TK            ! tmp stack exit temperature [K]
        REAL            VE            ! tmp stack exit velocity diameter [m/s]

        LOGICAL :: CFLAG    = .FALSE. ! true: convert from English to metric
        LOGICAL :: EFLAG    = .FALSE. ! true: error detected
        LOGICAL :: MAJRFLAG = .FALSE. ! true: use major/minor specifier
        LOGICAL :: PINGFLAG = .FALSE. ! true: output for plume-in-grid
        LOGICAL :: WFLAG    = .FALSE. ! true: convert lon to western

        CHARACTER*1     CSWITCH1  ! major/minor (TDEV) or PinG switch (GDEV)
        CHARACTER*1     CSWITCH2  ! PinG switch       
        CHARACTER*8     FMTFIP    ! format for writing co/st/cy code
        CHARACTER*80    GDESC               !  grid description
        CHARACTER*300   BUFFER
        CHARACTER*300   MESG

        CHARACTER(LEN=FIPLEN3) CFIP     !  char FIPS code
        CHARACTER(LEN=IOVLEN3) COORD3D  !  coordinate system name
        CHARACTER(LEN=IOVLEN3) COORUN3D !  coordinate system units 
        CHARACTER(LEN=ALLLEN3) CSRC     !  buffer for source char, incl pol/act
        CHARACTER(LEN=CHRLEN3) CHAR1    !  tmp plant characteristic 1
        CHARACTER(LEN=CHRLEN3) CHAR2    !  tmp plant charactersitic 2
        CHARACTER(LEN=PLTLEN3) PLT      !  tmp plant ID

        CHARACTER*16 :: PROGNAME = 'ELEVPOINT'   !  program name

C***********************************************************************
C   begin body of program ELEVPOINT

        LDEV = INIT3()

C.........  Write out copywrite, version, web address, header info, and prompt
C           to continue running the program.
        CALL INITEM( LDEV, SCCSW, PROGNAME )

C.........  Get environment variables that control this program
        MESG = 'Plume height elevated source cutoff [m]'
        CUTOFF = ENVREAL( 'SMK_CUTOFF_HT', MESG, 75., IOSCUT )

        MESG = 'Indicator for create plume-in-grid outputs'
        PINGFLAG = ENVYN( 'SMK_PING_YN', MESG, .FALSE., IOS )

        MESG = 'Indicator for defining major/minor sources'
        MAJRFLAG = ENVYN( 'SMK_SPECELEV_YN', MESG, .FALSE., IOS )

        MESG = 'Indicator for converting all longitudes to Western'
        WFLAG = ENVYN( 'WEST_HSPHERE', MESG, .TRUE., IOS )

        MESG = 'Indicator for English to metric units conversion'
        CFLAG = ENVYN( 'SMK_ENG2METRIC_YN', MESG, .FALSE., IOS )

C.........  Set source category based on environment variable setting
        CALL GETCTGRY

C.........  Get inventory file names given source category
        CALL GETINAME( CATEGORY, ENAME, ANAME )

C.........  Make sure only run for point sources
        IF( CATEGORY .NE. 'POINT' ) THEN
            MESG = 'ERROR: ' // PROGNAME( 1:LEN_TRIM( PROGNAME ) ) //
     &             ' is not valid for ' // CATEGORY( 1:CATLEN ) // 
     &             ' sources'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

C.........  Write note if using specific elevated sources from the stack splits
C           file, or if using the cutoff method
        IF( MAJRFLAG ) THEN
            MESG = 'NOTE: Major sources will be identified based on ' //
     &             'the stack splits file.'

        ELSE IF ( CUTOFF .GT. 0. ) THEN
            MESG = 'NOTE: Major sources will be defined using'
            L = LEN_TRIM( MESG )
            IF( IOSCUT .LT. 0 ) THEN
                WRITE( MESG,94020 ) MESG( 1:L ) // 
     &                 ' the default cutoff ' // CRLF() // BLANK10//
     &                 'height of', CUTOFF, '[m]'
            ELSE
                WRITE( MESG,94020 ) MESG( 1:L ) // 
     &                 ' a user-defined cutoff ' // CRLF() // BLANK10//
     &                 'height of', CUTOFF, '[m]'
            END IF

        ELSE
            WRITE( MESG,94020 ) 'WARNING: Major sources will be ' //
     &             'defined using a cutoff height of 0.'

        END IF

        CALL M3MSG2( MESG )

C.........  Create format for country/state/county code
        WRITE( FMTFIP, 94300 ) '(I', FIPLEN3, '.', FIPLEN3, ')'

C.......   Get file name; open input point source and output
C.......   elevated points files; get plume-rise cutoff for
C.......   elevated points file

C.........   Get file names and open inventory files
        ENAME = PROMPTMFILE( 
     &          'Enter logical name for the I/O API INVENTORY file',
     &          FSREAD3, ENAME, PROGNAME )
        ENLEN = LEN_TRIM( ENAME )

        SDEV = PROMPTFFILE( 
     &         'Enter logical name for the ASCII INVENTORY file',
     &         .TRUE., .TRUE., ANAME, PROGNAME )

C.............  For plume-in-grid inputs using any method of major source
C               identification
        IF( PINGFLAG ) THEN

C.............  Get stack split groups file
            GDEV = PROMPTFFILE( 
     &         'Enter logical name for the STACK SPLIT GROUPS file',
     &         .TRUE., .TRUE., CRL // 'GROUP', PROGNAME )

C.............  Get stack split definitions file
            TDEV = PROMPTFFILE( 
     &         'Enter logical name for the STACK SPLIT DEFINITION file',
     &         .TRUE., .TRUE., CRL // 'SPLIT', PROGNAME )

        END IF

C.........  Get header description of inventory file, error if problem
        IF( .NOT. DESC3( ENAME ) ) THEN
            MESG = 'Could not get description of file "' //
     &             ENAME( 1:ENLEN ) // '"'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

C.........  Otherwise, store source-category-specific header information, 
C           including the inventory pollutants in the file (if any).  Note that 
C           the I/O API head info is passed by include file and the
C           results are stored in module MODINFO.
        ELSE

            CALL GETSINFO

        END IF

C.........  Get episode information for setting date and time of STACK_PING file
        MESG = 'NOTE: Getting date/time information for use in ' //
     &         'STACK_PING file'
        CALL M3MSG2( MESG )

        CALL GETM3EPI( -9, SDATE, STIME, -9 )

C.........  Allocate memory for and read in required inventory characteristics
        CALL RDINVCHR( CATEGORY, ENAME, SDEV, NSRC, NINVARR, IVARNAMS )

C.........  Allocate memory for source status arrays and group numbers
        ALLOCATE( LMAJOR( NSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'LMAJOR', PROGNAME )
        ALLOCATE( LPING( NSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'LPING', PROGNAME )
        ALLOCATE( GROUPID( NSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'GROUPID', PROGNAME )

C.........  Initialize source status and group number arrays
        LMAJOR  = .FALSE.   ! array
        LPING   = .FALSE.   ! array
        GROUPID = 0         ! array

C.........  Get grid description for converting the stack group coordinates
C           if using PinG or cutoff method
        IF( .NOT. MAJRFLAG .OR. PINGFLAG ) THEN
            IF( .NOT. DSCM3GRD( 
     &                GDNAM3D, GDESC, COORD3D, GDTYP3D, COORUN3D,
     &                P_ALP3D, P_BET3D, P_GAM3D, XCENT3D, YCENT3D,
     &                XORIG3D, YORIG3D, XCELL3D, YCELL3D,
     &                NCOLS3D, NROWS3D, NTHIK3D ) ) THEN

        	MESG = 'Could not get Models-3 grid description.'
        	CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

            ELSE
                NCOLS = NCOLS3D
                NROWS = NROWS3D
                NGRID = NCOLS * NROWS
            END IF            

C.............  Allocate memory so that we can use the GENPTCEL
            ALLOCATE( NX( NGRID ), STAT=IOS )
            CALL CHECKMEM( IOS, 'NX', PROGNAME )

        END IF

C.............  If processing for models-3...
C.........  If using plume-and-grid inputs, allocate memory and read in 
C           stack groups file
        IF( PINGFLAG ) THEN

C.............  Get the number of lines in the stack groups file
            NGLINES = GETFLINE( GDEV, 'Stack groups' )

C.............  Scan the groups file to determine the number of PinG groups
            DO I = 1, NGLINES

                READ( GDEV, 93500, IOSTAT=IOS, END=999 ) 
     &                GID, CSWITCH1
                IREC = IREC + 1

C.................  Check read error status
        	IF( IOS .GT. 0 ) THEN
                    EFLAG = .TRUE.
                    WRITE( MESG,94010 ) 'I/O error', IOS, 
     &                 'reading stack groups file at line', IREC
                    CALL M3MESG( MESG )
                    CYCLE
        	END IF

                IF( CSWITCH1 .NE. ' ' ) NGROUP = NGROUP + 1

            END DO

            REWIND( GDEV )

C.............  Allocate memory for stack groups based on the number of lines
            NGROUP = MAX( NGROUP, 1 )
            CALL ALLOCATE_GROUP_VARIABLES

            MESG = 'Reading stack split groups file...'
            CALL M3MSG2( MESG )

            IF( CFLAG ) THEN
                MESG = 'NOTE: Converting stack parameters from ' //
     &                 'English to metric units'
        	CALL M3MSG2( MESG )
            END IF

C.............  Read stack groups file
            IREC = 0
            J    = 0
            DO I = 1, NGLINES

                READ( GDEV, 93500, IOSTAT=IOS, END=999 ) 
     &                GID, CSWITCH1, NPG, LON, LAT, DM, HT, TK, VE, FL
                IREC = IREC + 1

C.................  Check read error status
        	IF( IOS .GT. 0 ) THEN
                    EFLAG = .TRUE.
                    WRITE( MESG,94010 ) 'I/O error', IOS, 
     &                 'reading stack groups file at line', IREC
                    CALL M3MESG( MESG )
                    CYCLE
        	END IF

C.................  Skip entry if not a plume-in-grid source
                IF( CSWITCH1 .EQ. ' ' ) CYCLE

C.................  Convert longitude to western hemisphere if needed
                IF( WFLAG .AND. LON .GT. 0 ) LON = -LON

C.................  Convert from English to metric, if needed...
                IF( CFLAG ) THEN
                    DM = DM * FT2M
                    HT = HT * FT2M
                    TK = ( TK - 32. ) * FTOC + CTOK
                    VE = VE * FT2M
                    FL = FL * FLWE2M
                END IF

C.................  When flow is not defined, set it with the vel & diam
                IF( FL .LE. 0. .AND. VE .GT. 0. ) THEN
                   FL = VE * PI * ( 0.25 * DM * DM )
                END IF

C.................  Store data
                J = J + 1
                IF( J .LE. NGROUP ) THEN
                    GRPIDX ( J ) = J
                    GRPGIDA( J ) = GID
                    GRPCNT ( J ) = NPG 
                    IF( LON .NE. 0. ) GRPLON ( J ) = LON
                    IF( LON .NE. 0. ) GRPXL  ( J ) = LON
                    IF( LAT .NE. 0. ) GRPLAT ( J ) = LAT
                    IF( LAT .NE. 0. ) GRPYL  ( J ) = LAT
                    IF( DM  .GT. 0. ) GRPDM  ( J ) = DM
                    IF( HT  .GT. 0. ) GRPHT  ( J ) = HT
                    IF( TK  .GT. 0. ) GRPTK  ( J ) = TK
                    IF( VE  .GT. 0. ) GRPVE  ( J ) = VE
                    IF( FL  .GT. 0. ) GRPFL  ( J ) = FL
                END IF

            END DO    ! End loop on input file lines

C.............  Abort if overflow
            IF( J .GT. NGROUP ) THEN
                EFLAG = .TRUE.
                WRITE( MESG,94010 ) 
     &                  'INTERNAL ERROR: Number of stack groups ' //
     &                  'J=', J, 
     &                  'exceeds dimension NGROUP=', NGROUP
                CALL M3MSG2( MESG ) 

            END IF

C.............  Sort stack group information
            CALL SORTI1( NGROUP, GRPIDX, GRPGIDA )

C.............  Store sorted stack groups for lookups in reading stack splits 
C               file
            DO I = 1, NGROUP
                GRPGID( I ) = GRPGIDA( GRPIDX( I ) )
            END DO

        END IF   ! End if plume-in-grid

C.........  If major/minor definitions are to be used...
C.........  If doing cutoff method, but still want to identify PinG sources,
C           then the splits file is still needed (PINGFLAG = T)
        IF( MAJRFLAG .OR. PINGFLAG ) THEN

C.............  Allocate memory for reading stack splits file
            NSLINES = GETFLINE( TDEV, 'Stack splits' )

            ALLOCATE( SPTINDX( NSLINES ), STAT=IOS )
            CALL CHECKMEM( IOS, 'SPTINDX', PROGNAME )
            ALLOCATE( SPTGIDA( NSLINES ), STAT=IOS )
            CALL CHECKMEM( IOS, 'SPTGIDA', PROGNAME )
            ALLOCATE( SPTMMSA( NSLINES ), STAT=IOS )
            CALL CHECKMEM( IOS, 'SPTMMSA', PROGNAME )
            ALLOCATE( SPTMPSA( NSLINES ), STAT=IOS )
            CALL CHECKMEM( IOS, 'SPTMPSA', PROGNAME )
            ALLOCATE( SPTLON( NSLINES ), STAT=IOS )
            CALL CHECKMEM( IOS, 'SPTLON', PROGNAME )
            ALLOCATE( SPTLAT( NSLINES ), STAT=IOS )
            CALL CHECKMEM( IOS, 'SPTLAT', PROGNAME )
            ALLOCATE( SPTDM( NSLINES ), STAT=IOS )
            CALL CHECKMEM( IOS, 'SPTDM', PROGNAME )
            ALLOCATE( SPTHT( NSLINES ), STAT=IOS )
            CALL CHECKMEM( IOS, 'SPTHT', PROGNAME )
            ALLOCATE( SPTTK( NSLINES ), STAT=IOS )
            CALL CHECKMEM( IOS, 'SPTTK', PROGNAME )
            ALLOCATE( SPTVE( NSLINES ), STAT=IOS )
            CALL CHECKMEM( IOS, 'SPTVE', PROGNAME )
            ALLOCATE( SPTFL( NSLINES ), STAT=IOS )
            CALL CHECKMEM( IOS, 'SPTFL', PROGNAME )
            ALLOCATE( SPTCSRCA( NSLINES ), STAT=IOS )
            CALL CHECKMEM( IOS, 'SPTCSRCA', PROGNAME )
            ALLOCATE( SPTCSRC( NSLINES ), STAT=IOS )
            CALL CHECKMEM( IOS, 'SPTCSRC', PROGNAME )
            ALLOCATE( FOUND  ( NSLINES ), STAT=IOS )
            CALL CHECKMEM( IOS, 'FOUND', PROGNAME )

C.............  Initialize
            SPTINDX = 1

C.............  Initialize status of PSPLIT entries found in inventory
            FOUND = .FALSE.    ! array 
        
            MESG = 'Reading stack splits file...'
            CALL M3MSG2( MESG )

C.............  Read stack splits file
            IREC = 0
            DO I = 1, NSLINES

                READ( TDEV, 93550, IOSTAT=IOS, END=999 ) 
     &                GID, CSWITCH1, CSWITCH2, COID, STID, CYID, PLT,
     &                CHAR1, CHAR2, LON, LAT, DM, HT, TK, VE, FL
                IREC = IREC + 1

C.................  Check read error status
        	IF( IOS .GT. 0 ) THEN
                    EFLAG = .TRUE.
                    WRITE( MESG,94010 ) 'I/O error', IOS, 
     &                 'reading stack splits file at line', IREC
                    CALL M3MESG( MESG )
                    CYCLE
        	END IF

C.................  Store contents of splits file entry
                SPTINDX( I ) = I
                SPTGIDA( I ) = GID
                SPTMMSA( I ) = ( CSWITCH1 .NE. ' ' )
                SPTMPSA( I ) = ( CSWITCH2 .NE. ' ' )
                SPTLON ( I ) = LON
                SPTLAT ( I ) = LAT
                SPTDM  ( I ) = DM
                SPTHT  ( I ) = HT
                SPTTK  ( I ) = TK
                SPTVE  ( I ) = VE
                SPTFL  ( I ) = FL

                FIP = COID * 100000 + STID * 1000 + CYID
                WRITE( CFIP, FMTFIP ) FIP

                CSRC = ' '
                CALL BLDCSRC( CFIP, PLT, CHAR1, CHAR2, CHRBLNK3,
     &                        CHRBLNK3, CHRBLNK3, POLBLNK3, CSRC )

                SPTCSRCA( I ) = CSRC

            END DO    ! End loop on input file lines

            MESG = 'Processing splits data with inventory...'
            CALL M3MSG2( MESG )

C.............  Sort splits file source characteristics
            CALL SORTIC( NSLINES, SPTINDX, SPTCSRCA ) 

C.............  Store sorted splits file source characteristics for searching
            DO I = 1, NSLINES
                J = SPTINDX( I )
                SPTCSRC( I ) = SPTCSRCA( J )
            END DO

C.............  Loop through sources and match with records in splits file.
C.............  Flag and sources as major or plume-in-grid, and store group 
C               number for plume-in-grid sources
            DO S = 1, NSRC

                CSRC = CSOURC( S )( 1:PTENDL3( 4 ) )
                I = FINDC( CSRC, NSLINES, SPTCSRC )

                IF( I .LE. 0 ) THEN
                    CSRC = CSOURC( S )( 1:PTENDL3( 3 ) )
                    I = FINDC( CSRC, NSLINES, SPTCSRC )
                END IF

                IF( I .GT. 0 ) THEN
                    J = SPTINDX( I )
                    FOUND  ( I ) = .TRUE.
                    GROUPID( S ) = SPTGIDA( J )

C.....................  Find group ID in stack groups file to make sure that
C                       output is desired for this group
                    GID = SPTGIDA( J )
                    K = FIND1( GID, NGROUP, GRPGID )

C.....................  Skip source if stack group is not in stack group file
                    IF( K .LE. 0 ) CYCLE

C.....................  Store per-source major and PinG source info
                    LMAJOR ( S ) = ( MAJRFLAG .AND. SPTMMSA( J ) )
                    LPING  ( S ) = ( PINGFLAG .AND. SPTMPSA( J ) )

                    IF( LPING( S ) ) THEN
                        NPING  = NPING  + 1
                    ELSE IF( LMAJOR( S ) ) THEN
                        NMAJOR = NMAJOR + 1
                    END IF

C.....................  Check the stack groups file for stack parameters and
C                       coordinates. If the group is missing any information,
C                       update it with the data from the stack splits file.
C                       If the stack splits file is missing data, update with
C                       data from the inventory...

C.....................  Reassign index from position in sorted list to unsorted
                    K = GRPIDX( K )

                    CALL VALID_GRP_INFO( 'longitude', CSRC, GID, 
     &                         SPTLON( J ), XLOCA( S ), GRPLON( K ) )
                    CALL VALID_GRP_INFO( 'x-location', CSRC, GID, 
     &                         SPTLON( J ), XLOCA( S ), GRPXL ( K ) )
                    CALL VALID_GRP_INFO( 'latitude', CSRC, GID, 
     &                         SPTLAT( J ), YLOCA( S ), GRPLAT( K ) )
                    CALL VALID_GRP_INFO( 'y-location', CSRC, GID, 
     &                         SPTLAT( J ), YLOCA( S ), GRPYL( K ) )
                    CALL VALID_GRP_INFO( 'diameter', CSRC, GID, 
     &                         SPTDM( J ), STKDM( S ), GRPDM( K ) )
                    CALL VALID_GRP_INFO( 'height', CSRC, GID,
     &                         SPTHT( J ), STKHT( S ), GRPHT( K ) )
                    CALL VALID_GRP_INFO( 'exit temperature', CSRC, GID, 
     &                         SPTTK( J ), STKTK( S ), GRPTK( K ) )
                    CALL VALID_GRP_INFO( 'exit velocity', CSRC, GID, 
     &                         SPTVE( J ), STKVE( S ), GRPVE( K ) )
                    CALL VALID_GRP_INFO( 'exit flow', CSRC, GID, 
     &                         SPTFL( J ), -9., GRPFL( K ) )

                    IF( GRPFL( K ) .LT. 0. .AND. 
     &                  GRPVE( K ) .GT. 0.       ) THEN
                       DM = GRPDM( K )
                       GRPFL( K ) = GRPVE( K ) * PI * ( 0.25*DM*DM )

                    ELSE IF( GRPVE( K ) .LE. 0. ) THEN
                        EFLAG = .TRUE.
                        WRITE( MESG,94010 ) 'INTERNAL ERROR: ' //
     &                         'Bad velocity set for stack group',
     &                         GID, 'and could not compute flow.'
                        CALL M3MESG( MESG )

                    END IF

                             

                END IF

            END DO

C.............  Give warnings if any entries in the SPLITS file are not in the
C               inventory
            DO I = 1, NSLINES

                IF( .NOT. FOUND( I ) ) THEN

                    J    = SPTINDX ( I )
                    CSRC = SPTCSRCA( J )
                    CALL FMTCSRC( CSRC, NCHARS, BUFFER, L2 )

                    MESG = 'WARNING: Entry from PSPLIT not found ' //
     &                     'in inventory:' // CRLF() // BLANK10 // 
     &                     BUFFER( 1:L2 )              
                    CALL M3MESG( MESG )

                END IF

            END DO

        END IF  ! End of section for major/minor split file input

C.........  Abort if error found while reading stack group or stack splits file
        IF( EFLAG ) THEN
            MESG = 'Problem reading stack groups and/or stack ' //
     &             'splits file.'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

C.........  If using cutoff method...
        IF( .NOT. MAJRFLAG ) THEN

C.............  Write note 
            MESG = 'Computing plume rise and comparing to cutoff...'
            CALL M3MSG2( MESG )

C.............  Process the stacks to determine elevated sources
            DO S = 1, NSRC

C.................  Get stack parameters from stack groups file, if they
C                   are available. This will only be true for PinG source
C                   identification.
                GID = GROUPID( S )
                K = FIND1( GID, NGROUP, GRPGID )
                IF( K .GT. 0 ) THEN

                    K = GRPIDX( K )
                    HT = GRPHT( K )
                    TK = GRPTK( K )
                    DM = GRPDM( K )
                    VE = GRPVE( K )

                ELSE
                    HT = STKHT( S )
                    TK = STKTK( S )
                    DM = STKDM( S )
                    VE = STKVE( S )

                END IF

C.................  Check stack parameters so PLUMRIS doesn't blow up
C.................  If parameters are bad, skip plume rise calculation
                IF( HT .LT. 0. .OR. 
     &              TK .LE. 0. .OR.
     &              VE .LE. 0. .OR.
     &              DM .LE. 0.      ) THEN

                    EFLAG = .TRUE.
                    CALL FMTCSRC( CSRC, NCHARS, BUFFER, L2 )

                    WRITE( MESG,94030 ) HT, DM, TK, VE

                    L = LEN_TRIM( MESG )
                    MESG = 'ERROR: Invalid stack parameters for:' //
     &                     CRLF() // BLANK10 // 
     &                     BUFFER( 1:L2 )// ' with'// CRLF()// BLANK10//
     &                     MESG( 1:L )                
                    CALL M3MESG( MESG )

C.................  When stack parameters are okay...
                ELSE

C.....................  Calculate estimated plume rise
                    RISE = PLUMRIS( HT, TK, VE, DM )

C.....................  Identify sources as major when the plume rise is
C                       greater than the cutoff
                    IF( RISE .GT. CUTOFF ) THEN

                	NMAJOR = NMAJOR + 1
                	LMAJOR ( S ) = .TRUE.

                    END IF    ! if rise > cutoff

                END IF        ! end bad stack parms or not
            END DO            ! end loop on sources S

C.............  Abort if no sources are elevated
            IF( NMAJOR .EQ. 0 ) THEN

                MESG = 'WARNING: No sources found that meet the ' //
     &                 'elevated source cutoff height'
                CALL M3MSG2( MESG )
                CALL M3EXIT( PROGNAME, 0, 0, ' ', 0 )

            END IF

        END IF                ! end groups/splits file or not



C.........  If using cutoff method but not PinG, create the data for the stack
C           groups file.
        IF( .NOT. MAJRFLAG .AND. .NOT. PINGFLAG ) THEN

            NGROUP = NMAJOR
            CALL ALLOCATE_GROUP_VARIABLES

            J = 0            
            DO S = 1, NSRC

                IF( .NOT. LMAJOR( S ) ) CYCLE

                J = J + 1

                IF( J .LE. NGROUP ) THEN

                    GROUPID( S ) = J
                    GRPGID ( J ) = J
                    GRPGIDA( J ) = J
                    GRPIDX ( J ) = J
                    GRPXL  ( J ) = XLOCA( S )
                    GRPLON ( J ) = XLOCA( S )
                    GRPYL  ( J ) = YLOCA( S )
                    GRPLAT ( J ) = YLOCA( S )
                    GRPDM  ( J ) = STKDM( S ) 
                    GRPHT  ( J ) = STKHT( S )
                    GRPTK  ( J ) = STKTK( S )
                    GRPVE  ( J ) = STKVE( S )
                    GRPCNT ( J ) = 1

                    DM = GRPDM( J )
                    GRPFL( J ) = GRPVE( J ) * PI * ( 0.25*DM*DM )

                END IF

            END DO

C.............  Write internal error if dimensions exceeded
            IF( J .GT. NGROUP ) THEN

                WRITE( MESG,94010 ) 'INTERNAL ERROR: Stack group ' //
     &                 'arrays dimensioned to', NGROUP, 'but needed', J
                CALL M3MSG2( MESG )
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

            END IF

        END IF

C.........  Abort if an error occurred
        IF( EFLAG ) THEN
            MESG = 'Problem selecting major/plume-in-grid sources'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

C.........  Write status of processing
        IF( NGROUP .GT. 0 ) THEN
            WRITE( MESG,94010 ) 'Number of plume-in-grid source ' //
     &             '(MEPSE) groups:', NGROUP
            CALL M3MSG2( MESG )
        END IF

        IF( NMAJOR .GT. 0 ) THEN
            WRITE( MESG,94010 ) 'Number of major sources:', NMAJOR
            CALL M3MSG2( MESG )
        END IF

        IF( NPING .GT. 0 ) THEN
            WRITE( MESG,94010 ) 'Number of plume-in-grid sources:',NPING
            CALL M3MSG2( MESG )
        END IF

C.........  Open output files
        CALL OPENEOUT( NGROUP, SDATE, STIME, ENAME, PDEV, MNAME )

C.........  Write ASCII file
        MESG = 'Writing ELEVATED POINT SOURCE output file'
        CALL M3MSG2( MESG )

        DO S = 1, NSRC

            IF( LMAJOR( S ) .OR. LPING( S ) ) THEN
 
                MS  = 0
                PS  = 0
                GID = 0
                IF( LMAJOR( S ) ) MS = S
                IF( LPING ( S ) ) THEN
                    MS = 0
                    PS = S
                END IF
                GID = GROUPID( S )

                WRITE( PDEV, 93620 ) MS, PS, GID

            END IF

        END DO  

C.........  Sort and write plume-in-grid output file for Models-3 processing
C           or for elevated source identidied by cutoff method
        IF( .NOT. MAJRFLAG .OR. PINGFLAG ) THEN

C.............  Process the stack group coordinates for the current grid

C.............  Convert x,y location to coordinates of the projected grid
            CALL CONVRTXY( NGROUP, GDTYP3D, P_ALP3D, P_BET3D, 
     &                     P_GAM3D, XCENT3D, YCENT3D, GRPXL, GRPYL )

C.............  Determine grid cells for these coordinate locations
            CALL GENPTCEL( NGROUP, NGRID, GRPXL, GRPYL, NEXCLD, NX, 
     &                     INDX, GN, SN )

C.............  Convert grid cells to row and columns numbers
            DO I = 1, NGROUP

               ROW = 0
               COL = 0
               N   = GN( I )

               IF( N .GT. 0 ) THEN
                   ROW = N / NCOLS          ! note: integer math
                   IF( MOD( N, NCOLS ) .GT. 0. ) ROW = ROW + 1
                   COL = N - ( ROW-1 ) * NCOLS
               END IF

               GRPROW( I ) = ROW
               GRPCOL( I ) = COL

            END DO

C.............  Make sure that stack parameters are set for all groups
C.............  This is a simplistic way of doing this for now, later
C               add call to FIXSTK routine
            MINDM = MINVAL( GRPDM )
            MINHT = MINVAL( GRPHT )
            MINTK = MINVAL( GRPTK )
            MINVE = MINVAL( GRPVE )
            MINFL = MINVAL( GRPFL )

            IF( MIN( MINDM, MINHT, MINTK, MINVE, MINFL ) .LE. 0. ) THEN
                MESG = 'Bad stack group or stack split file. ' //
     &                 'Unable to assign stack ' // CRLF()//BLANK10//
     &                 'parameters to all stack groups. Could be '//
     &                 'a source matching problem.'
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF

C.............  Give warning if any plume-in-grid stack groups are outside the
C               grid
            IF( NEXCLD .GT. 0 ) THEN
                WRITE( MESG,94010 ) 'WARNING: ', NEXCLD, 'stack ' //
     &                 'groups are outside of grid "' // 
     &                 GDNAM3D( 1:LEN_TRIM( GDNAM3D ) ) // '"'
                CALL M3MSG2( MESG )
            END IF

            MESG='Writing ELEVATED/PING STACK PARAMETERS output file...'
            CALL M3MSG2( MESG )

            CALL WPINGSTK( MNAME, SDATE, STIME )

        END IF

C.........  Normal completion of program
        CALL M3EXIT( PROGNAME, 0, 0, ' ', 0 )

999     MESG = 'End of file reached unexpectedly'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

C******************  FORMAT  STATEMENTS   ******************************

C...........   Formatted file I/O formats............ 93xxx

93500   FORMAT( I6, A1, 21X, I5, F9.0, F9.0, 3X, F8.0, F7.0, F7.0, 
     &          F7.0, F10.0 )

93550   FORMAT( 6X, I6, A1, A1, I1, I2, I3, A15, A15, A11, 7X,
     &          F9.0, F9.0, F8.0, F7.0, F7.0, F7.0, F10.0 )

93620   FORMAT( 3(I8,1X) )

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10 ( A, :, I8, :, 2X  ) )

94020   FORMAT( A, 1X, F8.2, 1X, A )

94030   FORMAT( 'H[m]:', 1X, F6.2, 1X, 'D[m]:'  , 1X, F4.2, 1X,
     &          'T[K]:', 1X, F7.1, 1X, 'V[m/s]:', 1X, F10.1 )

94300   FORMAT( A, I2.2, A, I2.2, A )

C******************  INTERNAL SUBPROGRAMS  *****************************

        CONTAINS

C.............  This internal subprogram is for verifying that the group stack
C               parameters are valid and resetting them with the inventory
C               stack parameters if they are not.
            SUBROUTINE VALID_GRP_INFO( DESC, CSRC, GID, SPLVAL, 
     &                                 SRCVAL, GRPVAL )

C.............  Subroutine arguments 
            CHARACTER(*), INTENT (IN) :: DESC    ! description of data
            CHARACTER(*), INTENT (IN) :: CSRC    ! description of source chars
            INTEGER     , INTENT (IN) :: GID     ! group ID
            REAL        , INTENT (IN) :: SPLVAL  ! stack split data value
            REAL        , INTENT (IN) :: SRCVAL  ! source data value
            REAL        , INTENT(OUT) :: GRPVAL  ! stack group data value

C.............  Local variables
            INTEGER        L, L2
            CHARACTER*300  MESG

C----------------------------------------------------------------------

            L = LEN_TRIM( DESC )
            IF( GRPVAL .LT. AMISS3 ) THEN

                CALL FMTCSRC( CSRC, 4, BUFFER, L2 )

                IF( SPLVAL .GT. 0 ) THEN

                   GRPVAL = SPLVAL
                   WRITE( MESG,94010 ) 'WARNING: Stack group', GID, 
     &                    DESC( 1:L ) // ' set using stack splits'//
     &                    CRLF() // BLANK10 // 'file data for:'//
     &                    CRLF() // BLANK10 // BUFFER( 1:L2 )
                   CALL M3MESG( MESG )

                ELSE IF( SRCVAL .GE. 0 ) THEN

                   GRPVAL = SRCVAL
                   WRITE( MESG,94010 ) 'WARNING: Stack group', GID, 
     &                    DESC( 1:L ) // ' set using SMOKE'//
     &                    CRLF() // BLANK10 // 'inventory data for:'//
     &                    CRLF() // BLANK10 // BUFFER( 1:L2 )
                   CALL M3MESG( MESG )

                ELSE IF( DESC .EQ. 'exit flow' ) THEN
                    WRITE( MESG,94010 ) 'WARNING: Stack group', GID,
     &                    DESC( 1:L ) // ' set by computing'// CRLF()//
     &                    BLANK10// 'from velocity and diameter for:'//
     &                    CRLF()// BLANK10// BUFFER( 1:L2 )
                    CALL M3MESG( MESG )

                ELSE

                   EFLAG = .TRUE.
                   WRITE( MESG,94010 ) 'INTERNAL ERROR: ' //
     &                    DESC( 1:L ) // ' data in inventory less '// 
     &                    'than zero for:'//
     &                    CRLF() // BLANK10 // BUFFER( 1:L2 )
                   CALL M3MSG2( MESG )

                END IF

            END IF

            RETURN

C------------------- SUBPROGRAM FORMAT STATEMENTS ----------------------

C...........   Internal buffering formats............ 94xxx

94010       FORMAT( 10 ( A, :, I8, :, 2X  ) )

            END SUBROUTINE VALID_GRP_INFO 

C----------------------------------------------------------------------
C----------------------------------------------------------------------

C.............  This internal subprogram allocates memory for and initializes
C               the stack group arrays
            SUBROUTINE ALLOCATE_GROUP_VARIABLES

C----------------------------------------------------------------------

            ALLOCATE( GRPIDX( NGROUP ), STAT=IOS )
            CALL CHECKMEM( IOS, 'GRPIDX', PROGNAME )
            ALLOCATE( GRPGIDA( NGROUP ), STAT=IOS )
            CALL CHECKMEM( IOS, 'GRPGIDA', PROGNAME )
            ALLOCATE( GRPGID( NGROUP ), STAT=IOS )
            CALL CHECKMEM( IOS, 'GRPGID', PROGNAME )
            ALLOCATE( GRPLAT( NGROUP ), STAT=IOS )
            CALL CHECKMEM( IOS, 'GRPLAT', PROGNAME )
            ALLOCATE( GRPLON( NGROUP ), STAT=IOS )
            CALL CHECKMEM( IOS, 'GRPLON', PROGNAME )
            ALLOCATE( GRPDM( NGROUP ), STAT=IOS )
            CALL CHECKMEM( IOS, 'GRPDM', PROGNAME )
            ALLOCATE( GRPHT( NGROUP ), STAT=IOS )
            CALL CHECKMEM( IOS, 'GRPHT', PROGNAME )
            ALLOCATE( GRPTK ( NGROUP ), STAT=IOS )
            CALL CHECKMEM( IOS, 'GRPTK', PROGNAME )
            ALLOCATE( GRPVE( NGROUP ), STAT=IOS )
            CALL CHECKMEM( IOS, 'GRPVE', PROGNAME )
            ALLOCATE( GRPFL( NGROUP ), STAT=IOS )
            CALL CHECKMEM( IOS, 'GRPFL', PROGNAME )
            ALLOCATE( GRPCNT( NGROUP ), STAT=IOS )
            CALL CHECKMEM( IOS, 'GRPCNT', PROGNAME )
            ALLOCATE( GRPCOL( NGROUP ), STAT=IOS )
            CALL CHECKMEM( IOS, 'GRPCOL', PROGNAME )
            ALLOCATE( GRPROW( NGROUP ), STAT=IOS )
            CALL CHECKMEM( IOS, 'GRPROW', PROGNAME )
            ALLOCATE( GRPXL( NGROUP ), STAT=IOS )
            CALL CHECKMEM( IOS, 'GRPXL', PROGNAME )
            ALLOCATE( GRPYL( NGROUP ), STAT=IOS )
            CALL CHECKMEM( IOS, 'GRPYL', PROGNAME )
            ALLOCATE( INDX( NGROUP ), STAT=IOS )
            CALL CHECKMEM( IOS, 'INDX', PROGNAME )
            ALLOCATE( GN( NGROUP ), STAT=IOS )
            CALL CHECKMEM( IOS, 'GN', PROGNAME )
            ALLOCATE( SN( NGROUP ), STAT=IOS )
            CALL CHECKMEM( IOS, 'SN', PROGNAME )

C.............  Initialize stack group variables
            GRPIDX  = 0
            GRPGIDA = 0
            GRPGID  = 0
            GRPLAT  = BADVAL3
            GRPLON  = BADVAL3
            GRPDM   = BADVAL3
            GRPHT   = BADVAL3
            GRPTK   = BADVAL3
            GRPVE   = BADVAL3
            GRPFL   = BADVAL3
            GRPCNT  = 0
            GRPROW  = 0
            GRPCOL  = 0
            GRPXL   = BADVAL3
            GRPYL   = BADVAL3
            INDX    = 0
            GN      = 0
            SN      = 0

            END SUBROUTINE ALLOCATE_GROUP_VARIABLES

        END PROGRAM ELEVPOINT

