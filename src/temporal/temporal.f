
        PROGRAM TMPPOINT

C***********************************************************************
C  program body starts at line 
C
C  DESCRIPTION:
C
C  PRECONDITIONS REQUIRED:  
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C    copied by: M. Houyoux 01/99
C    origin: tmppoint.F 4.3
C
C***********************************************************************
C
C Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
C                System
C File: @(#)$Id$
C
C COPYRIGHT (C) 1998, MCNC--North Carolina Supercomputing Center
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
C*************************************************************************

C.........  MODULES for public variables
C.........  This module contains the inventory arrays
        USE MODSOURC

C.........  This module contains the temporal cross-reference tables
        USE MODXREF

C.........  This module contains the temporal profile tables
        USE MODTPRO

C.........  This module contains the information about the source category
        USE MODINFO

        IMPLICIT NONE
 
C.........  INCLUDES:
        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
        INCLUDE 'FDESC3.EXT'    !  I/O API file description data structures

C..........  EXTERNAL FUNCTIONS and their descriptions:

        LOGICAL         CHKEMEPI    ! checks the emissions episode parameters
        CHARACTER*2     CRLF
        INTEGER         ENVINT
        LOGICAL         ENVYN
        INTEGER         GETDATE
        INTEGER         GETFLINE
        INTEGER         GETNUM
        INTEGER         INDEX1
        CHARACTER*14    MMDDYY
        INTEGER         PROMPTFFILE
        CHARACTER*16    PROMPTMFILE
        INTEGER         RDTPROF

        EXTERNAL    CRLF, ENVINT, ENVYN, GETDATE, GETFLINE, GETNUM,
     &              INDEX1, MMDDYY, PROMPTFFILE, PROMPTMFILE, RDTPROF      
                                            
C.........  LOCAL PARAMETERS and their descriptions:

        CHARACTER*50, PARAMETER :: SCCSW = '@(#)$Id$'

C.........  Point sources work and output arrays
        REAL   , ALLOCATABLE :: EMIS ( :,: ) !  inventory emissions
        REAL   , ALLOCATABLE :: EMISV( :,: ) !  day-corrected emissions
        REAL   , ALLOCATABLE :: EMIST( :,: ) !  timestepped output emssions

C.........  Temporal allocation Matrix.  
        REAL, ALLOCATABLE :: TMAT( :, :, : ) ! temporal allocation factors

C.........  Array that contains the names of the inventory variables needed for
C           this program
        CHARACTER(LEN=IOVLEN3) IVARNAMS( MXINVARR )

C.........  Actual-SCC  table
        INTEGER                                NSCC
        CHARACTER(LEN=SCCLEN3), ALLOCATABLE :: SCCLIST( : )

C.........  Day-specific, hour-specific data, and elevated sources data. 
C.........  These need only to allow enough dimensions for one read per 
C           pollutant per time step.
        INTEGER                 NDAYPT       ! actual no. day-specific sources
        INTEGER, ALLOCATABLE :: INDXD( : )   ! SMOKE source IDs
        REAL   , ALLOCATABLE :: EMISD( : )   ! day-specific emissions

        INTEGER                 NHRPT        ! actual no. hour-specific sources
        INTEGER, ALLOCATABLE :: INDXH( : )   ! SMOKE source IDs
        REAL   , ALLOCATABLE :: EMISH( : )   ! hour-specific emissions

        INTEGER                 NPELV        ! optional elevated source-count
        INTEGER, ALLOCATABLE :: INDXE( : )   ! SMOKE source IDs
        REAL   , ALLOCATABLE :: EMISE( : )   ! elevated source emissions

C...........   Arrays for flexible application of day- and hour-specific data
        INTEGER                 NDAYPOL     ! no. of pols in day-spec file
        INTEGER                 NHRPOL      ! no. of pols in hr-spec file

        LOGICAL, ALLOCATABLE :: LDSPOL( : ) ! flags for day-spec pollutants
        LOGICAL, ALLOCATABLE :: LHSPOL( : ) ! flags for hr-spec pollutants

        CHARACTER(LEN=IOVLEN3), ALLOCATABLE :: DAYPNAM( : ) ! dy-spec pol names
        CHARACTER(LEN=IOVLEN3), ALLOCATABLE :: HRPNAM ( : ) ! hr-spec pol names

C.........  Reshaped inventory pollutants and associated variables
        INTEGER         NGRP                ! number of pollutant groups 
        INTEGER         NGSZ                ! number of pollutants per group 
        CHARACTER(LEN=IOVLEN3), ALLOCATABLE:: EINAM2D( :,: ) 

C...........   Logical names and unit numbers

        INTEGER         LDEV    !  unit number for log file
        INTEGER         RDEV    !  unit number for temporal profile file
        INTEGER         SDEV    !  unit number  for ASCII inventory file
        INTEGER         UDEV    !  unit no. for optional input elevated srcs
        INTEGER         XDEV    !  unit no. for cross-reference file

        CHARACTER*16    ANAME   !  logical name for ASCII inventory input file
        CHARACTER*16    ENAME   !  logical name for I/O API inventory input file
        CHARACTER*16 :: DNAME = 'NONE' !  " day-specific  input file, or "NONE"
        CHARACTER*16 :: HNAME = 'NONE' !  " hour-specific input file, or "NONE"
        CHARACTER*16    TNAME   !  " timestepped (low-level)   output file
        CHARACTER*16    UNAME   !  " (upper-level) output file, or "NONE"

C...........   Other local variables

        LOGICAL         DFLAG   !  day-specific  file available
        LOGICAL      :: EFLAG = .FALSE.  !  error-flag
        LOGICAL      :: EFLAG2= .FALSE.  !  error-flag (2)
        LOGICAL         HFLAG   !  hour-specific file available
        LOGICAL         PROMPTF !  iff PROMPTFLAG E.V. is true or not defined
        LOGICAL         UFLAG   !  generating upper-level output file

        INTEGER        I, J, K, L, L1, L2, N, S, T

        INTEGER         IOS, IOS1, IOS2, IOS3, IOS4 ! i/o status
        INTEGER         IOS5, IOS6, IOS7, IOS8      ! i/o status
        INTEGER         EDATE, ETIME        ! ending Julian date and time
        INTEGER         ENLEN               ! length of ENAME string
        INTEGER         JDATE, JTIME        ! Julian date and time
        INTEGER         JRUNLEN             ! models-3 "run length" in HHMMSS
        INTEGER         NLINE               ! tmp number of lines in ASCII file 
        INTEGER         NSRC                ! general number of sources
        INTEGER         NSTEPS              ! number of output time steps
        INTEGER         SDATE, STIME        ! starting Julian date and time
        INTEGER         TNLEN               ! length of TNAME string
        INTEGER         TSTEP               ! output time step
        INTEGER         TZONE               ! output-file time zone
        INTEGER         UNLEN               ! length of UNAME string

        CHARACTER*8     TREFFMT !  temporal cross-ref format (SOURCE|STANDARD)
        CHARACTER*14    DTBUF   !  buffer for MMDDYY
        CHARACTER*300   MESG    !  buffer for M3EXIT() messages
        CHARACTER(LEN=IOVLEN3) CBUF ! pollutant name temporary buffer 

        CHARACTER*16 :: PROGNAME = 'TMPPOINT' ! program name

C***********************************************************************
C   begin body of program TMPPOINT

        LDEV = INIT3()

C.........  Write out copywrite, version, web address, header info, and prompt
C           to continue running the program.
        CALL INITEM( LDEV, SCCSW, PROGNAME )

C.........  Get environment variables that control program behavior
C.........  Only abort if Models-3 environment variables are not set
        PROMPTF = ENVYN ( 'PROMPTFLAG', 'Prompt for inputs or not',
     &                    .FALSE., IOS )

        TZONE = ENVINT( 'OUTZONE', 'Output time zone', 0, IOS )

        DFLAG = ENVYN ( 'DAY_SPECIFIC_YN', 'Use day-specific data',
     &                   .FALSE., IOS )

        HFLAG = ENVYN ( 'HOUR_SPECIFIC_YN', 'Use hour-specific data',
     &                   .FALSE., IOS )

        UFLAG = ENVYN ( 'ELEVATED_PT_YN', 
     &                  'Generate elevated point-source file',
     &                   .FALSE., IOS )

C.........  Get Models-3-set episode settings
        SDATE  = ENVINT( 'G_STDATE', 'Start date (YYYYDDD)', 0, IOS )
        IF( IOS .NE. 0 ) THEN
            EFLAG = .TRUE.
            MESG = 'G_STDATE is not set properly to starting date!'
            CALL M3MSG2( MESG )
        ENDIF 

        STIME  = ENVINT( 'G_STTIME' , 'Start time (HHMMSS)', 0 , IOS )
        IF( IOS .NE. 0 ) THEN
            EFLAG = .TRUE.
            MESG = 'G_STTIME is not set properly to starting time!'
            CALL M3MSG2( MESG )
        ENDIF 

        TSTEP = ENVINT( 'G_TSTEP', 'Time step (HHMMSS)', 0, IOS )
        IF( IOS .NE. 0 ) THEN
            EFLAG = .TRUE.
            MESG = 'G_TSTEP is not set properly to time step!'
            CALL M3MSG2( MESG )
        ENDIF 

        JRUNLEN = ENVINT( 'G_RUNLEN', 'Duration (HHMMSS)', 0, IOS )
        IF( IOS .NE. 0 ) THEN
            EFLAG = .TRUE.
            MESG = 'G_RUNLEN is not set properly to duration!'
            CALL M3MSG2( MESG )
        ENDIF 

        IF( EFLAG ) THEN
            MESG = 'Bad environment variable setting(s)'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        ENDIF

C.........  Ensure that episode settings are consistent with SMOKE
        IF( .NOT. CHKEMEPI( SDATE, STIME, TSTEP, JRUNLEN ) ) THEN
            MESG = 'Invalid episode settings from environment'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        ENDIF

C.........   Set more time information based on environment variable inputs
        NSTEPS = JRUNLEN / 10000

C.........  Set source category based on environment variable setting
        CALL GETCTGRY

C.........  Get inventory file names given source category
        CALL GETINAME( CATEGORY, ENAME, ANAME )

C.........  Prompt for and open input I/O API and ASCII files
        ENAME = PROMPTMFILE( 
     &          'Enter logical name for the I/O API INVENTORY file',
     &          FSREAD3, ENAME, PROGNAME )
        ENLEN = LEN_TRIM( ENAME )

        SDEV = PROMPTFFILE( 
     &           'Enter logical name for the ASCII INVENTORY file',
     &           .TRUE., .TRUE., ANAME, PROGNAME )

        IF( DFLAG ) DNAME = PROMPTMFILE( 
     &          'Enter logical name for DAY-SPECIFIC file',
     &          FSREAD3, CRL // 'DAY', PROGNAME )

        IF( HFLAG ) HNAME = PROMPTMFILE( 
     &          'Enter logical name for HOUR-SPECIFIC file',
     &          FSREAD3, CRL // 'HOUR', PROGNAME )

        XDEV = PROMPTFFILE( 
     &           'Enter logical name for TEMPORAL XREF file',
     &           .TRUE., .TRUE., 'GTREF', PROGNAME )

        RDEV = PROMPTFFILE( 
     &           'Enter logical name for TEMPORAL PROFILES file',
     &           .TRUE., .TRUE., 'GTPRO', PROGNAME )

        IF( UFLAG ) UDEV = PROMPTFFILE( 
     &           'Enter logical name for ELEVATED SOURCES file',
     &            .TRUE., .TRUE., 'PELV', PROGNAME )

C.........  Get header description of inventory file 
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

C.............  Store non-category-specific header information
            NSRC = NROWS3D

        ENDIF

C.........  For day-specific data input...
        IF( DFLAG ) THEN

C.............  Get header description of day-specific input file
            IF( .NOT. DESC3( DNAME ) ) THEN
                CALL M3EXIT( PROGNAME, 0, 0, 
     &                       'Could not get description of file "' 
     &                       // DNAME( 1:LEN_TRIM( DNAME ) ) // '"', 2 )
            ENDIF

C.............  Allocate memory for pollutant pointer
            ALLOCATE( DAYPNAM( NVARS3D ), STAT=IOS )
            CALL CHECKMEM( IOS, 'DAYPNAM', PROGNAME )
           
C.............  Set day-specific file dates, check dates, and report problems
            CALL PDSETUP( DNAME, SDATE, STIME, EDATE, ETIME, NIPOL,  
     &                    EINAM, NDAYPOL, NDAYPT, EFLAG, DAYPNAM )

C.............  Allocate memory for reading day-specific emissions data
            ALLOCATE( INDXD( NDAYPT ), STAT=IOS )
            CALL CHECKMEM( IOS, 'INDXD', PROGNAME )
            ALLOCATE( EMISD( NDAYPT ), STAT=IOS )
            CALL CHECKMEM( IOS, 'EMISD', PROGNAME )

        ENDIF

C.........  For hour-specific data input...
        IF( HFLAG ) THEN

C............. Get header description of hour-specific input file
            IF( .NOT. DESC3( HNAME ) ) THEN
                CALL M3EXIT( PROGNAME, 0, 0, 
     &                       'Could not get description of file "' 
     &                       // HNAME( 1:LEN_TRIM( HNAME ) ) // '"', 2 )
            ENDIF

C.............  Allocate memory for pollutant pointer
            ALLOCATE( HRPNAM( NVARS3D ), STAT=IOS )
            CALL CHECKMEM( IOS, 'HRPNAM', PROGNAME )
           
C.............  Set day-specific file dates, check dates, and report problems
            CALL PDSETUP( HNAME, SDATE, STIME, EDATE, EDATE, NIPOL,  
     &                    EINAM, NHRPOL, NHRPT, EFLAG2, HRPNAM )

C.............  Allocate memory for reading hour-specific emissions data
            ALLOCATE( INDXH( NHRPT ), STAT=IOS )
            CALL CHECKMEM( IOS, 'INDXH', PROGNAME )
            ALLOCATE( EMISH( NHRPT ), STAT=IOS )
            CALL CHECKMEM( IOS, 'EMISH', PROGNAME )

        ENDIF

        IF( EFLAG .OR. EFLAG2 ) THEN
            MESG = 'Problem with day- or hour-specific inputs'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        ENDIF

C.........  Prompt for episode data only if PROMPTFLAG has been set to true
        IF( PROMPTF ) THEN

            TZONE = GETNUM( -12, 12, TZONE, 
     &             'Enter time zone (0 for GMT, 5 for EST, 4 for EDT)' )

            SDATE = GETDATE( SDATE,
     &           'Enter simulation starting date (YYYYDDD)|(YYYYMMDD)' )

            STIME = GETNUM( 0, 235959, STIME, 
     &                      'Enter simulation starting time (HHMMSS)' )
  
            NSTEPS= GETNUM( 1, 999999, NSTEPS,
     &                      'Enter output duration (hours)' )
        ENDIF

C.........  Report episode information that is being used for the output file(s)
        DTBUF = MMDDYY( SDATE )
        WRITE( MESG,94050 )
     &  'Output Time Zone:', TZONE,           CRLF() // BLANK5 //
     &  '      Start Date:', DTBUF( 1:LEN_TRIM( DTBUF ) ) //
     &                                        CRLF() // BLANK5 //
     &  '      Start Time:', STIME,'HHMMSS'// CRLF() // BLANK5 //
     &  '       Time Step:', 1    ,'hour'  // CRLF() // BLANK5 //
     &  '        Duration:', NSTEPS, 'hours'
 
        CALL M3MSG2( MESG )

        CALL M3MSG2( 'Reading source data from inventory file...' )

C.........  Set inventory variables to read for all source categories
        IVARNAMS( 1 ) = 'TZONES'
        IVARNAMS( 2 ) = 'TPFLAG'
        IVARNAMS( 3 ) = 'CSCC'

C.........  Allocate memory for and read required inventory characteristics
        IF( CATEGORY .EQ. 'AREA' ) THEN

        ELSE IF( CATEGORY .EQ. 'MOBILE' ) THEN

        ELSE IF( CATEGORY .EQ. 'POINT' ) THEN

            IVARNAMS( 4 ) = 'CSOURC'

            CALL RPNTSCHR( ENAME, SDEV, NSRC, 4, IVARNAMS, NCHARS )

        END IF

C.........  Build unique lists of SCCs per SIC from the inventory arrays
        CALL GENUSLST

C.........  Read elevated sources, set index to source list, and get actual 
C           number of valid entries in PELV file. Note that RDPELV will drop
C           records if they do not match inventory.
C.........  Allocate memory for EMISE
        IF( UFLAG ) THEN

            CALL M3MSG2( 'Reading elevated sources file...' )

            NLINE = GETFLINE( UDEV, 'Elevated sources file' )

            ALLOCATE( INDXE( NLINE ), STAT=IOS )
            CALL CHECKMEM( IOS, 'INDXE', PROGNAME )

            CALL RDPELV( UDEV, NSRC, NLINE, CSOURC, NPELV, INDXE )

            ALLOCATE( EMISE( NPELV ), STAT=IOS )
            CALL CHECKMEM( IOS, 'EMISE', PROGNAME )

C.............  Update UFLAG in case NPELV is zero
            IF( NPELV .EQ. 0 ) THEN
                UFLAG = .FALSE.
                MESG = 'No elevated sources matched the inventory ' //
     &                 'so no elevated file' // CRLF() // BLANK5 //
     &                 'will be written'
                CALL M3WARN( PROGNAME, 0, 0, MESG )

            ENDIF

        ENDIF

C.........  Read temporal-profile cross-reference file and put into tables
C.........  Only read entries for pollutants that are in the inventory.
        CALL RDTREF( XDEV, TREFFMT )

C.........  Read temporal-profiles file:  4 parts (monthly, weekly, 
C           weekday diurnal, and weekend diurnal)
        CALL M3MSG2( 'Reading temporal profiles file...' )

        NMON = RDTPROF( RDEV, 'MONTHLY' )
        NWEK = RDTPROF( RDEV, 'WEEKLY'  )
        NDIU = RDTPROF( RDEV, 'WEEKDAY' )
        NEND = RDTPROF( RDEV, 'WEEKEND' )

C.........  Adjust temporal profiles for use in generating temporal emissions
C.........  NOTE: All variables are passed by modules.
        CALL NORMTPRO

C.........  It is important that all major arrays must be allocated by this 
C           point because the next memory allocation step is going to pick a
C           data structure that will fit within the limits of the host.

C.........  Allocate memory, but allow flexibility in memory allocation
C           for second dimension.
C.........  The second dimension (the number of pollutants) can be different 
C           depending on the memory available.
C.........  To determine the approproate size, first attempt to allocate memory
C           for all pollutants to start, and if this fails, divide pollutants
C           into even groups and try again.

        NGSZ = NIPOL   ! Number of pollutant in each group
        NGRP = 1       ! Number of groups
        DO

            ALLOCATE( TMAT ( NSRC, NGSZ, 24 ), STAT=IOS1 )
            ALLOCATE( MDEX ( NSRC, NGSZ )    , STAT=IOS2 )
            ALLOCATE( WDEX ( NSRC, NGSZ )    , STAT=IOS3 )
            ALLOCATE( DDEX ( NSRC, NGSZ )    , STAT=IOS4 )
            ALLOCATE( EDEX ( NSRC, NGSZ )    , STAT=IOS5 )
            ALLOCATE( EMIS ( NSRC, NGSZ )    , STAT=IOS6 )
            ALLOCATE( EMIST( NSRC, NGSZ )    , STAT=IOS7 )
            ALLOCATE( EMISV( NSRC, NGSZ )    , STAT=IOS8 )

            IF( IOS1 .GT. 0 .OR. IOS2 .GT. 0 .OR. IOS3 .GT. 0 .OR.
     &          IOS4 .GT. 0 .OR. IOS5 .GT. 0 .OR. IOS6 .GT. 0 .OR.
     &          IOS7 .GT. 0 .OR. IOS8 .GT. 0 ) THEN

                IF( NGSZ .EQ. 1 ) THEN
                    J = 8 * NSRC * 31    ! Assume 8-byte reals
                    WRITE( MESG,94010 ) 
     &                'Insufficient memory to run program.' //
     &                CRLF() // BLANK5 // 'Could not allocate ' // 
     &                'pollutant-dependent block of', J, 'bytes.'
                    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                ENDIF

                NGRP = NGRP + 1
                NGSZ = NGSZ / NGRP + ( NIPOL - NGSZ * NGRP )

                DEALLOCATE( TMAT, MDEX, WDEX, DDEX, EDEX, EMIS, 
     &                      EMIST, EMISV )

            ELSE
                EXIT

            ENDIF

        ENDDO
        
C.........  Allocate a few small arrays based on the size of the groups
C.........  NOTE: this has a small potential for a problem if these little
C           arrays exceed the total memory limit.
        ALLOCATE( EINAM2D( NGSZ, NGRP ), STAT=IOS )
        CALL CHECKMEM( IOS, 'EINAM2D', PROGNAME )
        ALLOCATE( LDSPOL( NGSZ ), STAT=IOS )
        CALL CHECKMEM( IOS, 'LDSPOL', PROGNAME )
        ALLOCATE( LHSPOL( NGSZ ), STAT=IOS )
        CALL CHECKMEM( IOS, 'LHSPOL', PROGNAME )

C.........  Create 2-d array for storing pollutant names in groups
        EINAM2D = ' '
        EINAM2D = RESHAPE( EINAM, (/ NGSZ, NGRP /) )

C.........  Set up and open I/O API output file(s) header(s) variables...
        CALL OPENPTMP( ENAME, UFLAG, SDATE, STIME, TSTEP, TZONE, NPELV,
     &                 NIPOL, EINAM, TNAME, UNAME )

        TNLEN = LEN_TRIM( TNAME )
        UNLEN = LEN_TRIM( UNAME )

C.........  Loop through pollutant groups
        DO N = 1, NGRP

C.............  Write message stating the pollutants being processed
            CALL POLMESG( NGSZ, EINAM2D( 1,N ) )

C.............  Set up logical arrays that indicate which pollutants are day-
C               specific and which are hour-specific
            LDSPOL = .FALSE.   ! array
            DO I = 1, NDAYPOL
                J = INDEX1( DAYPNAM( I ), NGSZ, EINAM2D( 1,N ) )
                LDSPOL( J ) = .TRUE.
            END DO

            LHSPOL = .FALSE.   ! array
            DO I = 1, NHRPOL
                J = INDEX1( HRPNAM( I ), NGSZ, EINAM2D( 1,N ) )
                LHSPOL( J ) = .TRUE.
            END DO

C.............  Initialize emissions and other arrays for this pollutant group
            TMAT  = 0.
            MDEX  = IMISS3
            WDEX  = IMISS3
            DDEX  = IMISS3
            EDEX  = IMISS3
            EMIS  = 0.
            EMIST = 0.
            EMISV = 0.

C.............  Assign temporal profiles by source and pollutant
            CALL M3MSG2( 'Assigning temporal profiles to sources...' )

            CALL ASGNTPRO( NGSZ, EINAM2D( 1,N ), TREFFMT )

C.............  Read in pollutant emissions from inventory for current group
            DO I = 1, NGSZ

                CBUF = EINAM2D( I,N )

                IF( .NOT. READ3( ENAME, CBUF, ALLAYS3, 0, 0, 
     &                           EMIS( 1,I )                ) ) THEN
                    EFLAG = .TRUE.
                    L1 = LEN_TRIM( CBUF )
                    MESG = 'Error reading "' // CBUF( 1:L1 ) //
     &                     '" from file "' // ENAME( 1:ENLEN ) // '."'
                    CALL M3MSG2( MESG )

                ENDIF
            ENDDO

C NOTE: Note: add a control for processing emissions using days of the week
C   in the base year or the year of the G_SDATE (which might be a future year).
C   Also, need to get the base year from the FDESC3D packet, and add different
C   time variables (one for selecting day-of-week and one for output)

C.............  For each time step and pollutant in current group, generate
C               hourly emissions, write elevated emissions file (if any), and
C               write layer-1 emissions file (or all data).
            JDATE = SDATE
            JTIME = STIME
            DO T = 1, NSTEPS

C.................  Generate hourly emissions for current hour
                CALL GENHEMIS( CATEGORY, NSRC, NGSZ, NDAYPT, NHRPT, 
     &                         JDATE, JTIME, TZONE, DNAME, HNAME, 
     &                         LDSPOL, LHSPOL, EINAM2D( 1,N ), TZONES, 
     &                         TPFLAG, INDXD, EMISD, INDXH, EMISH, 
     &                         EMIS, EMISV, TMAT, EMIST )

C.................  Write index to elevated emissions file. Note that this
C                   is written for every time step, even though it currently
C                   does not depend on time.
                IF( UFLAG ) THEN 
                    IF ( .NOT. WRITE3( UNAME, 'INDXE', 
     &                                 JDATE, JTIME, INDXE ) ) THEN

                        MESG = 'Could not write "INDXE" to file "' //
     &                         UNAME( 1:UNLEN ) // '."'
                        CALL M3EXIT( PROGNAME, JDATE, JTIME, MESG, 2 )

                    ENDIF
                ENDIF

C.................  Loop through pollutants in this group
                DO I = 1, NGSZ

                    CBUF = EINAM2D( I,N )

C.....................  If there are elevated sources...
                    IF( UFLAG ) THEN 

C.........................  Updated compressed format elevated emissions for 
C                           the current pollutant I in group N
                        DO J = 1, NPELV
                            S = INDXE( J )
                            EMISE( J ) = EMIST( S,I )
                            EMIST( S,I ) = 0.0
                        ENDDO

C.........................  Write elevated emissions to I/O API NetCDF file
                        CBUF = EINAM2D( I,N )
                        IF( .NOT. WRITE3( UNAME, CBUF, 
     &                                    JDATE, JTIME, EMISE ) ) THEN

                            L = LEN_TRIM( CBUF )
                            MESG = 'Could not write "' // CBUF( 1:L ) //
     &                             '" to file "'//UNAME( 1:UNLEN )//'."'
                            CALL M3EXIT( PROGNAME,JDATE,JTIME,MESG,2 )

                        ENDIF

                    ENDIF       ! if writing elevated sources file

C.....................  Write hourly pollutant emissions to I/O API NetCDF file
                    IF( .NOT. WRITE3( TNAME, CBUF, JDATE, JTIME, 
     &                                EMIST( 1,I )              ) ) THEN

                        L = LEN_TRIM( CBUF )
                        MESG = 'Could not write "' // CBUF( 1:L ) // 
     &                         '" to file "' // UNAME( 1:UNLEN ) // '."'
                        CALL M3EXIT( PROGNAME, JDATE, JTIME, MESG, 2 )

                    ENDIF

                ENDDO  ! End loop on pollutants I in this group

C.................  Advance the date/time by one time step
                CALL NEXTIME( JDATE, JTIME, TSTEP )

            ENDDO      ! End loop on time steps T

        ENDDO          ! End loop on pollutant groups N

C.............  If not list format PTREF inputs, write out file of temporal 
C               profiles used per source

! NOTE: future

C.........  Exit program with normal completion
        CALL M3EXIT( PROGNAME, 0, 0, ' ', 0 )

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats.............94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

94050   FORMAT( A, 1X, I2.2, A, 1X, A, 1X, I6.6, 1X,
     &          A, 1X, I3.3, 1X, A, 1X, I3.3, 1X, A   )

        END PROGRAM TMPPOINT

