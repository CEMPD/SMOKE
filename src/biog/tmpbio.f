
        PROGRAM TMPBIO

C***********************************************************************
C  program body starts at line 219
C
C  DESCRIPTION:
C       Computes time stepped gridded biogenic emissions in terms of 
C       normalized gridded emissions from RAWBIO and postprocessed MM5
C       meteorology.
C
C  PRECONDITIONS REQUIRED:
C       Postprocessed MM5 meteorology that contains temperature and/or
C       solar radiation/cloud data. 
C       Normalized gridded emissions BGRD from RAWBIO or GRDBIO
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C       HRBIO, PREBMET
C
C  REVISION  HISTORY:
C       Prototype 11/99 by JMV from version 4.9 of TMPBIO SMOKE prototype
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
C***********************************************************************

C...........   Modules for public variables
C...........   This module contains the speciation profile tables
        USE MODSPRO, ONLY: MXSPEC, SPCNAMES

C...........   This module contains the global variables for the 3-d grid
        USE MODGRID, ONLY: NCOLS, NROWS, XOFF_A, YOFF_A, XOFF, YOFF

        IMPLICIT NONE

C...........   INCLUDES:

        INCLUDE 'PARMS3.EXT'      ! I/O API constants
        INCLUDE 'FDESC3.EXT'      ! I/O API file description data structure
        INCLUDE 'IODECL3.EXT'     ! I/O API function declarations
        INCLUDE 'EMCNST3.EXT'     !
        INCLUDE 'BIODIMS3.EXT'    ! biogenic-related constants


C...........   PARAMETERS and their descriptions:

        CHARACTER(50), PARAMETER :: CVSW = '$Name$' ! CVS release tag
     
C...........   EXTERNAL FUNCTIONS and their descriptions:

        CHARACTER(2)    CRLF   
        INTEGER         ENVINT 
        LOGICAL         ENVYN
        CHARACTER(50)   GETCFDSC
        INTEGER         GETDATE
        INTEGER         GETFLINE
        INTEGER         GETNUM
        LOGICAL         GETYN
        CHARACTER(10)   HHMMSS
        INTEGER         INDEX1
        CHARACTER(14)   MMDDYY
        CHARACTER(16)   PROMPTMFILE
        INTEGER         PROMPTFFILE
        INTEGER         SECSDIFF
        CHARACTER(16)   VERCHAR
        LOGICAL         SETENVVAR
        
        EXTERNAL        CRLF, ENVINT, ENVYN, GETDATE, GETFLINE, GETNUM, 
     &                  GETYN, HHMMSS, INDEX1, MMDDYY, PROMPTMFILE, 
     &                  PROMPTFFILE, SECSDIFF, VERCHAR, SETENVVAR

C.........  Meteorology data settings
        INTEGER, ALLOCATABLE :: METCHECK ( : )  ! dimension: no. time steps
                                                ! value indicates valid temperature file 
        CHARACTER(256), ALLOCATABLE :: METLIST( : ) ! list of temperature file names        
        INTEGER, ALLOCATABLE :: RADCHECK ( : )  ! valid radiation file for each time step
        CHARACTER(256), ALLOCATABLE :: RADLIST( : ) ! list of radiation file names

C.........  Gridded meteorology data
                
        REAL, ALLOCATABLE :: LAT  ( :, : )    !  grid lat (deg) -90 <= LAT <= 90
        REAL, ALLOCATABLE :: LON  ( :, : )    !  grid lon (deg) -180 <= LON <= 180 
        REAL, ALLOCATABLE :: TASFC ( :, : )    !  level-1 air  temperature (K)
        REAL, ALLOCATABLE :: PRSFC  ( :, : )    !  pressure (Pa)
        REAL, ALLOCATABLE :: TSOLAR ( :, :)     !  Photosynthetic Active Radiation (PAR)
        INTEGER, ALLOCATABLE :: SEASON( :, : )  !  Seasonal switch

C.......   Gridded normalized emissions to use in hrbio(s)

        REAL, ALLOCATABLE ::  PINE( :, :, : )         !   for pine
        REAL, ALLOCATABLE ::  DECD( :, :, : )         !   for deciduous forest
        REAL, ALLOCATABLE ::  CONF( :, :, : )         !   for coniferous forest
        REAL, ALLOCATABLE ::  AGRC( :, :, : )         !   for agricultural land
        REAL, ALLOCATABLE ::  LEAF( :, :, : )         !   for leaf area
        REAL, ALLOCATABLE ::  OTHR( :, :, : )         !   for other land

        REAL, ALLOCATABLE ::  AVLAI( :, : )           !  average LAI
        REAL, ALLOCATABLE ::  NORNO( :, :, : )        !  normalized NO emissions

C.......   Gridded winter normalized emissions from BGRDW file 

        REAL, ALLOCATABLE ::  PINEW( :, :, : )     !   for pine
        REAL, ALLOCATABLE ::  DECDW( :, :, : )     !   for deciduous forest
        REAL, ALLOCATABLE ::  CONFW( :, :, : )     !   for coniferous forest
        REAL, ALLOCATABLE ::  AGRCW( :, :, : )     !   for agricultural land
        REAL, ALLOCATABLE ::  LEAFW( :, :, : )     !   for leaf area
        REAL, ALLOCATABLE ::  OTHRW( :, :, : )     !   for other land

        REAL, ALLOCATABLE ::  AVLAIW( :, : )      !  average LAI
        REAL, ALLOCATABLE ::  NORNOW( :, :, :)    !  normalized NO emissions

C........  normalized emissions after seasonal switch 

        REAL, ALLOCATABLE ::  PINES( :, :, : )     !   for pine
        REAL, ALLOCATABLE ::  DECDS( :, :, : )     !   for deciduous forest
        REAL, ALLOCATABLE ::  CONFS( :, :, : )     !   for coniferous forest
        REAL, ALLOCATABLE ::  AGRCS( :, :, : )     !   for agricultural land
        REAL, ALLOCATABLE ::  LEAFS( :, :, : )     !   for leaf area
        REAL, ALLOCATABLE ::  OTHRS( :, :, : )     !   for other land

        REAL, ALLOCATABLE ::  AVLAIS( :, : )      !  average LAI
        REAL, ALLOCATABLE ::  NORNOS( :, :, :)    !  normalized NO emissions

 
        REAL, ALLOCATABLE ::  MLFAC( :, : )           !  mole factors 
        REAL, ALLOCATABLE ::  MSFAC( :, : )           !  mass factors (tons/hour)

C.......   BEIS2 internal, output species

        REAL, ALLOCATABLE :: EMPOL( :, :, : )         ! emissions of biogenic categories
        REAL, ALLOCATABLE :: EMISL( :, :, : )         ! emissions in moles/hour
        REAL, ALLOCATABLE :: EMISS( :, :, : )         ! emissions in tons/hour


        CHARACTER(5)     CTZONE     ! string of time zone
        CHARACTER(16)    RADNAM     !  string for shortwave radiation reaching ground
        CHARACTER(16)    TMPRNAM    !  string for temperature 
        CHARACTER(16)    PRESNAM    !  string for surface pressure
        CHARACTER(50) :: METSCEN    !  temporary string for met scenario name
        CHARACTER(50) :: CLOUDSHM   !  temporary string for cloud scheme name
        CHARACTER(50) :: LUSE       !  temporary string for land use description
        CHARACTER(50) :: LUSE2      !  temporary string for 2nd land use desc.

C.......   Name tables for file NNAME

        CHARACTER(16), ALLOCATABLE ::  NORMV( : )   ! names for VOC vbles
        CHARACTER(16), ALLOCATABLE ::  NORMN( : )   ! names for  NO-emission vbles
        CHARACTER(16), ALLOCATABLE ::  EMSPC( : )   ! names of emitting species 

        CHARACTER(SPNLEN3)       SPPRO        ! speciation profile to use

        CHARACTER(72)    PARMENU( 5 )            ! Methods to calc. PAR
        DATA     PARMENU
     &           / 'Use MM5 generated RGND or RSD',
     &             'Use KUO cloud attenuation',
     &             'Use KF  cloud attenuation',   
     &             'Use No deep convection param'  ,
     &             'Assume Clear Skies' /

C...........   Logical names and unit numbers

        INTEGER         LDEV    !  unit number for log device
        INTEGER         RDEV    !  unit number for speciation profiles file
        INTEGER         MDEV    !  unit number for temperature list file
        INTEGER         DDEV    !  unit number for radiation list file
            
        CHARACTER(16)   ENAME   !  logical name for emissions output (moles)
        CHARACTER(16)   SNAME   !  logical name for emissions output (mass)
        CHARACTER(16)   NNAME   !  logical name for normalized-emissions input
        CHARACTER(16)   NNAME2  !  logical name for 2nd norm emissions input
        CHARACTER(16)   BNAME   !  logical name for frost switch input
        CHARACTER(16)   GNAME   !  logical name for GRID_CRO_2D
        CHARACTER(16)   MNAME   !  logical name for gridded temperature file
        CHARACTER(16)   RNAME   !  logical name for gridded radiation file

C...........   Other variables and their descriptions:

        INTEGER         B, M    !  counters for biogenic, model species
        INTEGER         I, II, III, J, JJ, JJJ, K, L, C, N, R  !  loop counters and subscripts
        INTEGER         HR      !  current simulation hour
        INTEGER         HRPOS   !  current position in file checking arrays

        INTEGER         IOS     !  temporay IO status
        INTEGER         JDATE   !  current simulation date (YYYYDDD) in output time zone
        INTEGER         JTIME   !  current simulation time (HHMMSS) in output time zone
        INTEGER         LDATE   !  previous simulation date
        INTEGER         MDATE   !  meteorology start date in GMT
        INTEGER         MTIME   !  meteorology start time in GMT
        INTEGER         METNCOLS!  no. met file columns
        INTEGER         METNROWS!  no. met file rows
        INTEGER         METSTART  ! starting position in METCHECK array
        INTEGER         MSPCS   !  no. of emitting species
        INTEGER         NLINES  !  no. of lines in GSPRO speciation profiles file 
        INTEGER         NMETLINES ! no. of lines in temperature file list
        INTEGER         NRADLINES ! no. of lines in radiation file list 
        INTEGER      :: NSTEPS = 1!  duration of episode
        INTEGER         NMETSTEPS ! no. of time steps in temperature files
        INTEGER         NRADSTEPS ! no. of time steps in radiation files
        INTEGER         NEWSTEPS  ! temporary no. of time steps in simulation
        INTEGER         PARTYPE !  method number to calculate PAR
        INTEGER         RADSTART  ! starting position in RADCHECK array
        INTEGER      :: SWNCOLS = 0   !  bioseason 
        INTEGER      :: SWNROWS = 0   !  bioseason 
        INTEGER      :: SWXOFF  = 0   !  bioseason x offset from met grid
        INTEGER      :: SWYOFF  = 0   !  bioseason y offset from met grid
        INTEGER         TSTEP   !  Time step from environment
        INTEGER         TZONE   !  output-file time zone

        LOGICAL ::      EFLAG    = .FALSE.  ! error flag
        LOGICAL ::      SAMEFILE = .TRUE.   ! radiation/cld and tmpr data in same file 
        LOGICAL ::      SWITCH_FILE = .TRUE.  ! use frost switch file
        LOGICAL ::      DUPWARN  = .TRUE.   ! warn if duplicate data found
        LOGICAL ::      METOPEN = .FALSE.  ! true: a met I/O API file is open
        LOGICAL ::      RADOPEN = .FALSE.  ! true: a rad I/O API file is open

        CHARACTER(256)  METFILE  !  full gridded temperature file name
        CHARACTER(256)  RADFILE  !  full gridded radiation/cloud file name
        CHARACTER(256)  DUPFILE  !  duplicate file name
        CHARACTER(256)  NEXTFILE !  next file name to be opened
        
        CHARACTER(300)  MESG    !  message buffer for M3EXIT()

        CHARACTER(14)   DTBUF        ! date buffer
        
        CHARACTER(16) :: PROGNAME = 'TMPBIO'   !  program name

C***********************************************************************
C   begin body of program TMPBIO

        LDEV = INIT3()
 
C.........  Write out copywrite, version, web address, header info, and prompt
C           to continue running the program.

        CALL INITEM( LDEV, CVSW, PROGNAME )

C.........  Evaluate the environment variables...

C.........  Get the time zone for output of the emissions
        TZONE = ENVINT( 'OUTZONE', 'Output time zone', 0, IOS )

C.......   Check to see if frost date switch file to be used
        MESG = 'Using a frost date switch file?'
        SWITCH_FILE = ENVYN ( 'BIOSW_YN', MESG, .TRUE., IOS )

C.........  Write time zone to character string
        WRITE( CTZONE,94000 ) TZONE

        MESG = 'Speciation profile to use for biogenics'
        CALL ENVSTR( 'BIOG_SPRO', MESG, '0000', SPPRO, IOS )

C........  Are the rad/cld data in the same file as tmpr data

        MESG = 'Radiation/cloud in same file as temperature data?'
        SAMEFILE = ENVYN ( 'BIOMET_SAME', MESG, .FALSE., IOS )

        RDEV = PROMPTFFILE(
     &           'Enter logical name for SPECIATION PROFILES file',
     &           .TRUE., .TRUE., 'GSPRO', PROGNAME )


C.........  Scan speciation profiles file to get all of the pollutant-species
C           combinations that are valid for the pollutants in the inventory.
C.........  The species names are sorted in ABC order for each pollutant, and
C           and the pollutants are in the same order as BIOSPC.
C.........  Also retrieve the maximum number of species per pollutant and
C           maximum number of profile entries per pollutant.

        CALL DSCSPROF( RDEV, BSPCS, BIOSPC )

        NLINES = GETFLINE( RDEV, 'Speciation profile file' )

        ALLOCATE( EMSPC( NLINES ), STAT=IOS )
        CALL CHECKMEM( IOS, 'EMSPC', PROGNAME )

        MSPCS = 0

C.......... Find emitting species names

        DO I = 1, BSPCS
          DO J = 1, MXSPEC
            IF ( SPCNAMES( J, I ) .EQ. ' ' ) CYCLE

            K = INDEX1 ( SPCNAMES( J, I ), MSPCS, EMSPC ) 

            IF ( K .EQ. 0 ) THEN
               MSPCS = MSPCS + 1
               EMSPC( MSPCS ) = SPCNAMES( J, I )
            END IF
          END DO
        END DO

        ALLOCATE( MLFAC ( MSPCS, BSPCS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'LAT', PROGNAME )

        ALLOCATE( MSFAC ( MSPCS, BSPCS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'LAT', PROGNAME )
 
 
C.............  Read speciation profiles file
        MESG = BLANK5 // 'Reading biogenic speciation profile...'
        CALL M3MSG2( MESG )

        CALL RDBPRO( RDEV, SPPRO, BSPCS, BIOSPC, MSPCS, EMSPC, 
     &               MLFAC, MSFAC ) 

C.......   Get the method to calculate PAR

        PARTYPE = ENVINT( 'BG_CLOUD_TYPE', 
     &                    'How PAR will be calculated', 1, IOS )

        MESG =  'PAR calculation will ' // PARMENU( PARTYPE )  
        WRITE( LDEV, 92000 ) MESG

C.......   Get episode date and time from environment
        CALL GETM3EPI( TZONE, JDATE, JTIME, TSTEP, NSTEPS )

C........  Convert start date and time to output time zone.
        MDATE = JDATE
        MTIME = JTIME
        CALL NEXTIME( MDATE, MTIME, TZONE*10000 )

C......    Get name of temperature variable to use
        MESG = 'Variable name for temperature'
        CALL ENVSTR( 'TMPR_VAR', MESG, 'TA', TMPRNAM, IOS )

C.......   Get list of temperature files
        MDEV = PROMPTFFILE(
     &         'Enter name for list of gridded temperature input ' //
     &         'files (or "NONE")',
     &         .TRUE., .TRUE., 'METLIST', PROGNAME )

        MNAME = 'MET_FILE1'

C.......   If we don't have a list, get location of MET_FILE1 to use
        IF( MDEV < 0 ) THEN
            CALL ENVSTR( MNAME, 'Gridded temperature input file',
     &                   ' ', METFILE, IOS )
            IF( IOS /= 0 ) THEN
                MESG = 'No setting for gridded temperature file'
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF
            
            NMETLINES = 1  ! set one line in file
            
C.......   Otherwise, get number of lines in met list file
        ELSE
            NMETLINES = GETFLINE( MDEV, 'METLIST file' )
            
        END IF

C.......   Allocate and initialize arrays for storing met list and met check        
        ALLOCATE( METLIST( NMETLINES ), STAT=IOS )
        CALL CHECKMEM( IOS, 'METLIST', PROGNAME )
        ALLOCATE( METCHECK( NSTEPS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'METCHECK', PROGNAME )
        
        METLIST = ' '   ! array
        METCHECK = 0    ! array

C.......   Read lines from met list file   
        IF( MDEV > 0 ) THEN     
            CALL RDLINES( MDEV, 'METLIST file', NMETLINES, METLIST )
        ELSE
            METLIST( 1 ) = METFILE
        END IF
        
        IF( MDEV .GT. 0 ) CLOSE( MDEV )

        MESG = 'Checking temperature files...'
        CALL M3MSG2( MESG )

C.......   Loop through met files and check time period
        DO N = 1, NMETLINES
        
C.............  Close previous file if needed
            IF( METOPEN ) THEN
                IF( .NOT. CLOSE3( MNAME ) ) THEN
                    MESG = 'Could not close temperature file ' // 
     &                     METFILE( 1:LEN_TRIM( METFILE ) )
                    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                ELSE
                    METOPEN = .FALSE.
                END IF
            END IF

C.............  Get current file name
            METFILE = METLIST( N )

C.............  Reset duplicate warning - prints error once per met file
            DUPWARN = .TRUE.    
            
C.............  Skip any blank lines
            IF( METFILE == ' ' ) CYCLE

C.............  Set logical file name
            IF( .NOT. SETENVVAR( MNAME, METFILE ) ) THEN
            	EFLAG = .TRUE.
            	MESG = 'INTERNAL ERROR: Could not set logical file ' //
     &                 'name for file ' // TRIM( METFILE )
                CALL M3MSG2( MESG )
                CYCLE
            END IF       	

C.............  Try to open file            
            IF( .NOT. OPEN3( MNAME, FSREAD3, PROGNAME ) ) THEN
                EFLAG = .TRUE.
                MESG = 'ERROR: Could not open temperature file ' //
     &                 METFILE( 1:LEN_TRIM( METFILE ) )
                CALL M3MSG2( MESG )
                CYCLE
            ELSE
            	METOPEN = .TRUE.
            END IF

C.............  Read description of file            
            IF( .NOT. DESC3( MNAME ) ) THEN
                EFLAG = .TRUE.
                MESG = 'ERROR: Could not get description of file ' //
     &                 METFILE( 1:LEN_TRIM( METFILE ) )
                CALL M3MESG( MESG )
                CYCLE
            END IF

C.............  Check that requested variable is in file
            J = INDEX1( TMPRNAM, NVARS3D, VNAME3D )
            IF( J <= 0 ) THEN
            	EFLAG = .TRUE.
                MESG = 'ERROR: Could not find ' // TMPRNAM // 
     &                 'in file ' //
     &                 METFILE( 1:LEN_TRIM( METFILE ) )
                CALL M3MESG( MESG )
                CYCLE
            END IF

C.............  Check that all met files have the same grid.
            CALL CHKGRID( MNAME, 'GRID' , 0 , EFLAG )

C.............  If first file, store met grid settings and descriptions
            IF( N == 1 ) THEN
                METNCOLS = NCOLS
                METNROWS = NROWS

                METSCEN  = GETCFDSC( FDESC3D,'/MET SCENARIO/', .FALSE. )
                CLOUDSHM = GETCFDSC( FDESC3D,'/CLOUD SCHEME/', .FALSE. )
            END IF

C.............  Use start date, start time, and no. time steps to fill in METCHECK array
            HRPOS = SECSDIFF( MDATE, MTIME, SDATE3D, STIME3D )
            HRPOS = HRPOS / 3600

C.............  If met data starts before episode, skip ahead to start of episode
            IF( HRPOS <= 0 ) THEN
                NMETSTEPS = MXREC3D + HRPOS
                HRPOS = 1
            ELSE
            	NMETSTEPS = MXREC3D
                HRPOS = HRPOS + 1
            END IF

C.............  Loop through time steps in met data and fill in METCHECK array            
            DO J = HRPOS, NMETSTEPS
                IF( J > NSTEPS ) EXIT

C.................  Check if met data overlaps previous file and print warning                
                IF( METCHECK( J ) /= 0 .AND. DUPWARN ) THEN
                    DUPFILE = METLIST( METCHECK( J ) )
                    MESG = 'WARNING: Time period of temperature file '//
     &                     METFILE( 1:LEN_TRIM( METFILE ) ) //
     &                     ' overlaps that of file ' //
     &                     DUPFILE( 1:LEN_TRIM( DUPFILE ) ) // '.' //
     &                     CRLF() // BLANK10 // 'Data from ' // 
     &                     METFILE( 1:LEN_TRIM( METFILE ) ) //
     &                     ' will be used.' 
                    CALL M3MESG( MESG )
                    DUPWARN = .FALSE.
                END IF
                
                METCHECK( J ) = N
            END DO
        
        END DO

C.........  Exit if there was a problem with the temperature files
        IF( EFLAG ) THEN
            MESG = 'Problems checking temperature files'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

C.........  Check for missing temperature data
        MESG = 'Checking for missing temperature data...'
        CALL M3MSG2( MESG )

        NEWSTEPS = 0
        METSTART = 0

        DO N = 1, NSTEPS

            IF( METCHECK( N ) /= 0 ) THEN
                IF( METSTART == 0 ) THEN
                    NEWSTEPS = 1
                    CALL NEXTIME( JDATE, JTIME, (N-1)*10000 )
                    CALL NEXTIME( MDATE, MTIME, (N-1)*10000 )
                    METSTART = N
                ELSE
                    NEWSTEPS = NEWSTEPS + 1
                END IF
            ELSE
                IF( METSTART /= 0 ) EXIT
            END IF
        
        END DO

        IF( NEWSTEPS == 0 ) THEN
            MESG = 'Temperature data does not cover ' //
     &             'requested time period'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

C.........  Reset number of steps and print messages
        IF( NEWSTEPS /= NSTEPS ) THEN
            NSTEPS = NEWSTEPS
            MESG = 'WARNING: Resetting output start date, start ' //
     &             'time, or duration to match available temperature' //
     &             ' data.'
            CALL M3MSG2( MESG )

C.............  Get date buffer field
            DTBUF = MMDDYY( JDATE )
            L = LEN_TRIM( DTBUF )
            WRITE( MESG,94052 )
     &        '      Start Date:', DTBUF( 1:L ) //  CRLF() // BLANK5 //
     &        '      Start Time:', JTIME,'HHMMSS'// CRLF() // BLANK5 //
     &        '        Duration:', NSTEPS, 'hours'
            CALL M3MSG2( MESG )
        END IF

        IF ( PARTYPE .NE. 5 ) THEN

            IF ( PARTYPE .EQ. 1 ) THEN

C.................  Get name of radiation variable to use
                MESG = 'Variable name for radiation'
                CALL ENVSTR( 'RAD_VAR', MESG, 'RGRND', RADNAM, IOS )

            END IF

C............  Open second met file if needed

            IF ( .NOT. SAMEFILE ) THEN

                DDEV = PROMPTFFILE( 'Enter name for list of gridded ' //
     &              'radiation/cloud input files (or "NONE")', 
     &              .TRUE., .TRUE., 'RADLIST', PROGNAME )

                RNAME = 'MET_FILE2'

C.................  If we don't have a list, get location of MET_FILE2 to use
                IF( DDEV < 0 ) THEN
                    CALL ENVSTR( RNAME, 'Gridded radiation/' //
     &                           'cloud input file', ' ', RADFILE, IOS )
                    IF( IOS /= 0 ) THEN
                        MESG = 'No setting for gridded ' //
     &                         'radiation/cloud file'
                        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                    END IF
                    
                    NRADLINES = 1

C.................  Otherwise, get number of lines in met list file                    
                ELSE
                    NRADLINES = GETFLINE( DDEV, 'RADLIST file' )
                    
                END IF

C.................  Allocate and initialize arrays for storing rad list and rad check        
                ALLOCATE( RADLIST( NRADLINES ), STAT=IOS )
                CALL CHECKMEM( IOS, 'RADLIST', PROGNAME )
                ALLOCATE( RADCHECK( NSTEPS ), STAT=IOS )
                CALL CHECKMEM( IOS, 'RADCHECK', PROGNAME )
        
                RADLIST = ' '   ! array
                RADCHECK = 0    ! array

C.................  Read lines from rad list file   
                IF( DDEV > 0 ) THEN     
                    CALL RDLINES( DDEV, 'RADLIST file', NRADLINES, 
     &                            RADLIST )
                ELSE
                    RADLIST( 1 ) = RADFILE
                END IF

                IF( DDEV .GT. 0 ) CLOSE( DDEV )

                MESG = 'Checking radiation/cloud files...'
                CALL M3MSG2( MESG )

C.................  Loop through radiation files and check time period
                DO N = 1, NRADLINES

C.....................  Close previous file if needed
                    IF( RADOPEN ) THEN
                        IF( .NOT. CLOSE3( RNAME ) ) THEN
                            MESG = 'Could not close radiation file ' // 
     &                             RADFILE( 1:LEN_TRIM( RADFILE ) )
                            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                        ELSE
                            RADOPEN = .FALSE.
                        END IF
                    END IF

C.....................  Get current file name
                    RADFILE = RADLIST( N )

C.....................  Reset duplicate warning - prints error once per met file
                    DUPWARN = .TRUE.    
            
C.....................  Skip any blank lines
                    IF( RADFILE == ' ' ) CYCLE
                    
C.....................  Set env variable
                    IF( .NOT. SETENVVAR( RNAME, RADFILE ) ) THEN
                    	EFLAG = .TRUE.
                    	MESG = 'INTERNAL ERROR: Could not set ' //
     &               	       'logical file name for file ' //
     &                         TRIM( RADFILE )
                        CALL M3MSG2( MESG )
                        CYCLE
                    END IF
                    
C.....................  Try to open file   
                    IF( .NOT. OPEN3( RNAME, FSREAD3, 
     &                                  PROGNAME ) ) THEN
                        EFLAG = .TRUE.
                        MESG = 'ERROR: Could not open radiation ' // 
     &                         'file ' // 
     &                         RADFILE( 1:LEN_TRIM( RADFILE ) )
                        CALL M3MESG( MESG )
                        CYCLE
                    ELSE
                    	RADOPEN = .TRUE.
                    END IF

C.....................  Read description of file            
                    IF( .NOT. DESC3( RNAME ) ) THEN
                        EFLAG = .TRUE.
                        MESG = 'ERROR: Could not get description ' //
     &                         'of file ' //
     &                         RADFILE( 1:LEN_TRIM( RADFILE ) )
                        CALL M3MESG( MESG )
                        CYCLE
                    END IF

                    IF( PARTYPE .EQ. 1 ) THEN

C.........................  Check that requested variable is in file
                        J = INDEX1( RADNAM, NVARS3D, VNAME3D )
                        IF( J <= 0 ) THEN
                        	EFLAG = .TRUE.
                            MESG = 'ERROR: Could not find ' // RADNAM // 
     &                             'in file ' //
     &                             RADFILE( 1:LEN_TRIM( RADFILE ) )
                            CALL M3MESG( MESG )
                            CYCLE
                        END IF
                    END IF

C.....................  Check that all met files have the same grid.
                    CALL CHKGRID( RNAME, 'GRID' , 0 , EFLAG )

C.....................  If first file, store description
                    IF( N == 1 ) THEN
                        CLOUDSHM = GETCFDSC( FDESC3D, '/CLOUD SCHEME/', 
     &                                      .FALSE. )
                    END IF

C.....................  Use start date, start time, and no. time steps 
C                       to fill in RADCHECK array
                    HRPOS = SECSDIFF( MDATE, MTIME, SDATE3D, STIME3D )
                    HRPOS = HRPOS / 3600

C.....................  If met data starts before episode, 
C                       skip ahead to start of episode
                    IF( HRPOS <= 0 ) THEN
                        NRADSTEPS = MXREC3D + HRPOS
                        HRPOS = 1
                    ELSE
                    	NRADSTEPS = MXREC3D
                        HRPOS = HRPOS + 1
                    END IF

C.....................  Loop through time steps in met data and fill in RADCHECK array            
                    DO J = HRPOS, NRADSTEPS
                        IF( J > NSTEPS ) EXIT

C.........................  Check if met data overlaps previous file and print warning                
                        IF( RADCHECK( J ) /= 0 .AND. DUPWARN ) THEN
                            DUPFILE = RADLIST( RADCHECK( J ) )
                            MESG = 'WARNING: Time period of ' //
     &                             'radiation file '//
     &                             RADFILE( 1:LEN_TRIM( RADFILE ) ) //
     &                             ' overlaps that of file ' //
     &                             DUPFILE( 1:LEN_TRIM( DUPFILE ) ) // 
     &                             '.' //
     &                             CRLF() // BLANK10 // 'Data from ' // 
     &                             RADFILE( 1:LEN_TRIM( RADFILE ) ) //
     &                             ' will be used.' 
                            CALL M3MESG( MESG )
                            DUPWARN = .FALSE.
                        END IF
                        
                        RADCHECK( J ) = N
                    END DO
        
                END DO

C................  Exit if there was a problem with the radiation files
                IF( EFLAG ) THEN
                    MESG = 'Problems checking radiation/cloud files'
                    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                END IF

C.................  Check for missing radiation data
                MESG = 'Checking for missing radiation data...'
                CALL M3MSG2( MESG )
                
                NEWSTEPS = 0
                RADSTART = 0

                DO N = 1, NSTEPS

                    IF( RADCHECK( N ) /= 0 ) THEN
                        IF( RADSTART == 0 ) THEN
                            NEWSTEPS = 1
                            CALL NEXTIME( JDATE, JTIME, (N-1)*10000 )
                            CALL NEXTIME( MDATE, MTIME, (N-1)*10000 )
                            RADSTART = N
                        ELSE
                            NEWSTEPS = NEWSTEPS + 1
                        END IF
                    ELSE
                        IF( RADSTART /= 0 ) EXIT
                    END IF
        
                END DO

                IF( NEWSTEPS == 0 ) THEN
                    MESG = 'Radiation data does not cover ' //
     &                     'requested time period'
                    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                END IF

C................  Reset METSTART
                METSTART = METSTART + ( RADSTART - 1 )

C................  Reset number of steps and print messages
                IF( NEWSTEPS /= NSTEPS ) THEN
                    NSTEPS = NEWSTEPS
                    MESG = 'WARNING: Resetting output start date, ' //
     &                     'start time, or duration to match ' //
     &                     'available radiation data.'
                    CALL M3MSG2( MESG )

C....................  Get date buffer field
                    DTBUF = MMDDYY( JDATE )
                    L = LEN_TRIM( DTBUF )
                    WRITE( MESG,94052 )
     &        '      Start Date:', DTBUF( 1:L ) //  CRLF() // BLANK5 //
     &        '      Start Time:', JTIME,'HHMMSS'// CRLF() // BLANK5 //
     &        '        Duration:', NSTEPS, 'hours'
                    CALL M3MSG2( MESG )
                END IF
        
            ELSE       
                ALLOCATE( RADLIST( NMETLINES ), STAT=IOS )
                CALL CHECKMEM( IOS, 'RADLIST', PROGNAME )
                ALLOCATE( RADCHECK( NSTEPS ), STAT=IOS )
                CALL CHECKMEM( IOS, 'RADCHECK', PROGNAME )

                RADLIST  = METLIST
                RADCHECK = METCHECK
 
            END IF
        
        END IF

C........  if solar zenith angle calculation needed then get lat-lon
C........  coordinates from GRID_CRO_2D file

        IF ( PARTYPE .GT. 1  ) THEN

            GNAME = PROMPTMFILE( 
     &              'Enter name for 2D GRID PARAMETERS input file',
     &              FSREAD3, 'GRID_CRO_2D', 'TMPBIO' )

            IF ( .NOT. DESC3( GNAME ) ) THEN
                CALL M3EXIT( 'TMPBIO', 0, 0,
     &                      'Could not get description of file "'
     &                      // GNAME( 1:LEN_TRIM( GNAME ) ) // '"', 2 )

            END IF 

C.............  Check that all met files have the same grid.
            CALL CHKGRID( GNAME, 'GRID' , 0 , EFLAG ) 

        END IF

C.......    Get bioseason switch file, BIOSEASON
        IF ( SWITCH_FILE ) THEN

           BNAME = PROMPTMFILE( 
     &          'Enter name for season switch input file',
     &          FSREAD3, 'BIOSEASON', 'TMPBIO' )
           
C......    Read description of switch file

           IF ( .NOT. DESC3( BNAME ) ) THEN

              MESG = 'Could not get description of file "' //
     &             NNAME( 1:LEN_TRIM( BNAME ) ) // '"'
              CALL M3EXIT( 'TMPBIO', 0, 0, MESG, 2 )

           END IF

C............  Compare grid definition and call with temporary subgrid flag
C              (since this subgrid can be different from final subgrid)
           CALL CHKGRID( BNAME, 'GRID' , 2 , EFLAG )
           SWXOFF = XOFF_A
           SWYOFF = YOFF_A

        END IF

        SWNCOLS = NCOLS3D
        SWNROWS = NROWS3D

C.......   Get normalized emissions file, BGRD

        NNAME = PROMPTMFILE( 
     &          'Enter name for NORMALIZED EMISSIONS input file',
     &          FSREAD3, 'BGRD', 'TMPBIO' )

C......    Read description of normalized emissions file

        IF ( .NOT. DESC3( NNAME ) ) THEN

            MESG = 'Could not get description of file "' //
     &             NNAME( 1:LEN_TRIM( NNAME ) ) // '"'
            CALL M3EXIT( 'TMPBIO', 0, 0, MESG, 2 )

        END IF

C.........  Final grid definition 
        CALL CHKGRID( NNAME, 'GRID' , 1, EFLAG )

        LUSE  = GETCFDSC( FDESC3D, '/LANDUSE/', .FALSE. )

        IF ( SWITCH_FILE ) THEN

C.......   Get winter normalized emissions file, BGRDW
C.......   Note that second normalized file is assumed to be
C.......   the winter file

           NNAME2 = PROMPTMFILE( 
     &          'Enter name for winter NORMALIZED EMISSIONS input file',
     &          FSREAD3, 'BGRDW', 'TMPBIO' )

C......    Read description of second normalized emissions file

           IF ( .NOT. DESC3( NNAME2 ) ) THEN

              MESG = 'Could not get description of file "' //
     &             NNAME( 1:LEN_TRIM( NNAME2 ) ) // '"'
              CALL M3EXIT( 'TMPBIO', 0, 0, MESG, 2 )

           END IF

C............  Check normalized emissions are consistent with the subgrid.
           CALL CHKGRID( NNAME2, 'GRID' , 1 , EFLAG )

           LUSE2  = GETCFDSC( FDESC3D, '/LANDUSE/', .FALSE. )

        END IF

C........  If grid definitions do not match properly
        IF ( EFLAG ) THEN
          MESG = 'Problems opening input files. See ERROR(S) above.'
          CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

C.......   Build description for, and create/open output file

        SDATE3D = JDATE
        STIME3D = JTIME
        MXREC3D = NSTEPS
        NVARS3D = MSPCS
        NLAYS3D = 1
        TSTEP3D = 10000    ! only 1-hour time step supported

        DO  M = 1, MSPCS
            VNAME3D( M ) = EMSPC( M )
            UNITS3D( M ) = 'moles/hr'
            VDESC3D( M ) = 'biogenic emissions of the indicated species'
            VTYPE3D( M ) = M3REAL
        END DO

        FDESC3D = ' '   ! array

        FDESC3D( 1 ) = 'Gridded biogenic emissions from SMOKE-BEIS2'
        FDESC3D( 2 ) = '/FROM/ '    // PROGNAME
        FDESC3D( 3 ) = '/VERSION/ ' // VERCHAR( CVSW )
        FDESC3D( 4 ) = '/TZONE/ '   // CTZONE
        IF ( SWITCH_FILE ) THEN
           FDESC3D( 5 ) = '/LANDUSE/ ' // LUSE // LUSE2
        ELSE
           FDESC3D( 5 ) = '/LANDUSE/ ' // LUSE
        END IF

        FDESC3D( 6 ) = '/MET SCENARIO/ ' // METSCEN
        FDESC3D( 7 ) = '/CLOUD SCHEME/ ' // CLOUDSHM

        ENAME = PROMPTMFILE(
     &          'Enter name for BGTS output file - moles',
     &          FSUNKN3, 'BGTS_L', 'TMPBIO' )

        DO M = 1, MSPCS
            UNITS3D( M ) = 'tons/hr'
        END DO

        SNAME = PROMPTMFILE(
     &          'Enter name for BGTS output file - mass',
     &          FSUNKN3, 'BGTS_S', 'TMPBIO' )


C........  if solar zenith angle calculation needed then get lat-lon
C........  coordinates from GRID_CRO_2D file
        IF ( PARTYPE .GT. 1  ) THEN

            ALLOCATE( LAT ( METNCOLS, METNROWS ), STAT=IOS )
            CALL CHECKMEM( IOS, 'LAT', PROGNAME )

            ALLOCATE( LON ( METNCOLS, METNROWS ), STAT=IOS )
            CALL CHECKMEM( IOS, 'LON', PROGNAME )

            IF ( .NOT. READ3( GNAME, 'LAT', 1, 0, 0, LAT ) ) THEN
              MESG = 'Could not read LAT from file "' //
     &                GNAME( 1:LEN_TRIM( GNAME ) ) // '"'
              CALL M3EXIT( 'TMPBIO', 0, 0, MESG, 2 )
            END IF

            IF ( .NOT. READ3( GNAME, 'LON', 1, 0, 0, LON ) ) THEN
              MESG = 'Could not read LON from file "' //
     &                GNAME( 1:LEN_TRIM( GNAME ) ) // '"'
              CALL M3EXIT( 'TMPBIO', 0, 0, MESG, 2 )
            END IF

        END IF    ! Use cloud to calculate PAR

C.......   Build name table for variables in normalized emissions file"

        ALLOCATE( NORMV( BTYPES * ( BSPCS - 1 ) ), STAT=IOS )
        CALL CHECKMEM( IOS, 'NORMV', PROGNAME )

        ALLOCATE( NORMN( LUSES ), STAT=IOS )
        CALL CHECKMEM( IOS, 'NORMN', PROGNAME )

        I = 0

        DO  B = 1, BSPCS - 1
          DO  K = 1, BTYPES

            I = I + 1
            NORMV( I ) = BIOLTYPE( K ) // BIOSPC( B )

          END DO
        END DO

        DO  L = 1, LUSES

            NORMN( L ) = BIOLUSE( L )( 1:LEN_TRIM( BIOLUSE( L )))//'NO'

        END DO


C.......   Allocate memory for normalized emissions
        ALLOCATE( PINE( NCOLS, NROWS, BSPCS-1  ), STAT=IOS )
        CALL CHECKMEM( IOS, 'PINE', PROGNAME )

        ALLOCATE( DECD( NCOLS, NROWS, BSPCS-1  ), STAT=IOS )
        CALL CHECKMEM( IOS, 'DECD', PROGNAME )
        ALLOCATE( CONF( NCOLS, NROWS, BSPCS-1  ), STAT=IOS )
        CALL CHECKMEM( IOS, 'CONF', PROGNAME )
        ALLOCATE( AGRC( NCOLS, NROWS, BSPCS-1  ), STAT=IOS )
        CALL CHECKMEM( IOS, 'AGRC', PROGNAME )
        ALLOCATE( LEAF( NCOLS, NROWS, BSPCS-1  ), STAT=IOS )
        CALL CHECKMEM( IOS, 'LEAF', PROGNAME )
        ALLOCATE( OTHR( NCOLS, NROWS, BSPCS-1  ), STAT=IOS )
        CALL CHECKMEM( IOS, 'OTHR', PROGNAME )

        ALLOCATE( AVLAI( NCOLS, NROWS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'AVLAI', PROGNAME )

        ALLOCATE( NORNO( NCOLS, NROWS, LUSES ), STAT=IOS )
        CALL CHECKMEM( IOS, 'NORNO', PROGNAME )

        IF ( SWITCH_FILE ) THEN
           ALLOCATE( PINEW( NCOLS, NROWS, BSPCS-1 ), STAT=IOS )
           CALL CHECKMEM( IOS, 'PINEW', PROGNAME )
           ALLOCATE( DECDW( NCOLS, NROWS, BSPCS-1 ), STAT=IOS )
           CALL CHECKMEM( IOS, 'DECDW', PROGNAME )
           ALLOCATE( CONFW( NCOLS, NROWS, BSPCS-1 ), STAT=IOS )
           CALL CHECKMEM( IOS, 'CONFW', PROGNAME )
           ALLOCATE( AGRCW( NCOLS, NROWS, BSPCS-1 ), STAT=IOS )
           CALL CHECKMEM( IOS, 'AGRCW', PROGNAME )
           ALLOCATE( LEAFW( NCOLS, NROWS, BSPCS-1 ), STAT=IOS )
           CALL CHECKMEM( IOS, 'LEAFW', PROGNAME )
           ALLOCATE( OTHRW( NCOLS, NROWS, BSPCS-1 ), STAT=IOS )
           CALL CHECKMEM( IOS, 'OTHRW', PROGNAME )
           ALLOCATE( AVLAIW( NCOLS, NROWS ), STAT=IOS )
           CALL CHECKMEM( IOS, 'AVLAIW', PROGNAME )
           ALLOCATE( NORNOW( NCOLS, NROWS, LUSES ), STAT=IOS )
           CALL CHECKMEM( IOS, 'NORNOW', PROGNAME )

           ALLOCATE( PINES( NCOLS, NROWS, BSPCS-1 ), STAT=IOS )
           CALL CHECKMEM( IOS, 'PINEW', PROGNAME )
           ALLOCATE( DECDS( NCOLS, NROWS, BSPCS-1 ), STAT=IOS )
           CALL CHECKMEM( IOS, 'DECDW', PROGNAME )
           ALLOCATE( CONFS( NCOLS, NROWS, BSPCS-1 ), STAT=IOS )
           CALL CHECKMEM( IOS, 'CONFW', PROGNAME )
           ALLOCATE( AGRCS( NCOLS, NROWS, BSPCS-1 ), STAT=IOS )
           CALL CHECKMEM( IOS, 'AGRCW', PROGNAME )
           ALLOCATE( LEAFS( NCOLS, NROWS, BSPCS-1 ), STAT=IOS )
           CALL CHECKMEM( IOS, 'LEAFW', PROGNAME )
           ALLOCATE( OTHRS( NCOLS, NROWS, BSPCS-1 ), STAT=IOS )
           CALL CHECKMEM( IOS, 'OTHRW', PROGNAME )
           ALLOCATE( AVLAIS( NCOLS, NROWS ), STAT=IOS )
           CALL CHECKMEM( IOS, 'AVLAIW', PROGNAME )
           ALLOCATE( NORNOS( NCOLS, NROWS, LUSES ), STAT=IOS )
           CALL CHECKMEM( IOS, 'NORNOW', PROGNAME )

           ALLOCATE( SEASON( SWNCOLS, SWNROWS), STAT=IOS )
           CALL CHECKMEM( IOS, 'SEASON', PROGNAME )
           SEASON = 0   ! array

        END IF

C.......   Loops reading the various categories of normalized emissions:
        I = 0

        DO  M = 1, BSPCS - 1

            I = I + 1

            IF ( .NOT. READ3( NNAME, NORMV( I ), 1, 0, 0, 
     &                      PINE( 1, 1, M ) ) ) THEN
                 MESG = 'Could not read "' // 
     &                  NORMV( I )( 1 : LEN_TRIM( NORMV( I ) ) ) //
     &                  '" from file "' //
     &                  NNAME( 1:LEN_TRIM( NNAME ) ) // '"'
                CALL M3EXIT( 'TMPBIO', 0, 0, MESG, 2 )
            END IF
            IF ( SWITCH_FILE ) THEN
               IF ( .NOT. READ3( NNAME2, NORMV( I ), 1, 0, 0, 
     &                      PINEW( 1, 1, M ) ) ) THEN
                  MESG = 'Could not read "' // 
     &                  NORMV( I )( 1 : LEN_TRIM( NORMV( I ) ) ) //
     &                  '" from file "' //
     &                  NNAME2( 1:LEN_TRIM( NNAME2 ) ) // '"'
                  CALL M3EXIT( 'TMPBIO', 0, 0, MESG, 2 )
               END IF
            END IF

            I = I + 1

            IF ( .NOT. READ3( NNAME, NORMV( I ), 1, 0, 0,
     &                      DECD( 1, 1, M ) ) ) THEN
                 MESG = 'Could not read "' //
     &                  NORMV( I )( 1 : LEN_TRIM( NORMV( I ) ) ) //
     &                  '" from file "' //
     &                  NNAME( 1:LEN_TRIM( NNAME ) ) // '"'
                CALL M3EXIT( 'TMPBIO', 0, 0, MESG, 2 )
            END IF
            IF ( SWITCH_FILE ) THEN
               IF ( .NOT. READ3( NNAME2, NORMV( I ), 1, 0, 0, 
     &                      DECDW( 1, 1, M ) ) ) THEN
                  MESG = 'Could not read "' // 
     &                  NORMV( I )( 1 : LEN_TRIM( NORMV( I ) ) ) //
     &                  '" from file "' //
     &                  NNAME2( 1:LEN_TRIM( NNAME2 ) ) // '"'
                  CALL M3EXIT( 'TMPBIO', 0, 0, MESG, 2 )
               END IF
            END IF

            I = I + 1

            IF ( .NOT. READ3( NNAME, NORMV( I ), 1, 0, 0,
     &                      CONF( 1, 1, M ) ) ) THEN
                 MESG = 'Could not read "' //
     &                  NORMV( I )( 1 : LEN_TRIM( NORMV( I ) ) ) //
     &                  '" from file "' //
     &                  NNAME( 1:LEN_TRIM( NNAME ) ) // '"'
                CALL M3EXIT( 'TMPBIO', 0, 0, MESG, 2 )
            END IF
            IF ( SWITCH_FILE ) THEN
               IF ( .NOT. READ3( NNAME2, NORMV( I ), 1, 0, 0, 
     &                      CONFW( 1, 1, M ) ) ) THEN
                  MESG = 'Could not read "' // 
     &                  NORMV( I )( 1 : LEN_TRIM( NORMV( I ) ) ) //
     &                  '" from file "' //
     &                  NNAME2( 1:LEN_TRIM( NNAME2 ) ) // '"'
                  CALL M3EXIT( 'TMPBIO', 0, 0, MESG, 2 )
               END IF
            END IF

            I = I + 1

            IF ( .NOT. READ3( NNAME, NORMV( I ), 1, 0, 0,
     &                      AGRC( 1, 1, M ) ) ) THEN
                 MESG = 'Could not read "' //
     &                  NORMV( I )( 1 : LEN_TRIM( NORMV( I ) ) ) //
     &                  '" from file "' //
     &                  NNAME( 1:LEN_TRIM( NNAME ) ) // '"'
                CALL M3EXIT( 'TMPBIO', 0, 0, MESG, 2 )
            END IF
            IF ( SWITCH_FILE ) THEN
               IF ( .NOT. READ3( NNAME2, NORMV( I ), 1, 0, 0, 
     &                      AGRCW( 1, 1, M ) ) ) THEN
                  MESG = 'Could not read "' // 
     &                  NORMV( I )( 1 : LEN_TRIM( NORMV( I ) ) ) //
     &                  '" from file "' //
     &                  NNAME2( 1:LEN_TRIM( NNAME2 ) ) // '"'
                  CALL M3EXIT( 'TMPBIO', 0, 0, MESG, 2 )
               END IF
            END IF 

            I = I + 1

            IF ( .NOT. READ3( NNAME, NORMV( I ), 1, 0, 0,
     &                      LEAF( 1, 1, M ) ) ) THEN
                 MESG = 'Could not read "' //
     &                  NORMV( I )( 1 : LEN_TRIM( NORMV( I ) ) ) //
     &                  '" from file "' //
     &                  NNAME( 1:LEN_TRIM( NNAME ) ) // '"'
                CALL M3EXIT( 'TMPBIO', 0, 0, MESG, 2 )
            END IF
            IF ( SWITCH_FILE ) THEN
               IF ( .NOT. READ3( NNAME2, NORMV( I ), 1, 0, 0, 
     &                      LEAFW( 1, 1, M ) ) ) THEN
                  MESG = 'Could not read "' // 
     &                  NORMV( I )( 1 : LEN_TRIM( NORMV( I ) ) ) //
     &                  '" from file "' //
     &                  NNAME2( 1:LEN_TRIM( NNAME2 ) ) // '"'
                  CALL M3EXIT( 'TMPBIO', 0, 0, MESG, 2 )
               END IF
            END IF

            I = I + 1

            IF ( .NOT. READ3( NNAME, NORMV( I ), 1, 0, 0,
     &                      OTHR( 1, 1, M ) ) ) THEN
                 MESG = 'Could not read "' //
     &                  NORMV( I )( 1 : LEN_TRIM( NORMV( I ) ) ) //
     &                  '" from file "' //
     &                  NNAME( 1:LEN_TRIM( NNAME ) ) // '"'
                CALL M3EXIT( 'TMPBIO', 0, 0, MESG, 2 )
            END IF
            IF ( SWITCH_FILE ) THEN
               IF ( .NOT. READ3( NNAME2, NORMV( I ), 1, 0, 0, 
     &                      OTHRW( 1, 1, M ) ) ) THEN
                  MESG = 'Could not read "' // 
     &                  NORMV( I )( 1 : LEN_TRIM( NORMV( I ) ) ) //
     &                  '" from file "' //
     &                  NNAME2( 1:LEN_TRIM( NNAME2 ) ) // '"'
                  CALL M3EXIT( 'TMPBIO', 0, 0, MESG, 2 )
               END IF
            END IF

        END DO        !  end loop to read normalized emissions 

        IF ( .NOT. READ3( NNAME, 'AVLAI', 1, 0, 0, AVLAI ) ) THEN
             MESG = 'Could not read AVLAI from file "' //
     &              NNAME( 1:LEN_TRIM( NNAME ) ) // '"'
            CALL M3EXIT( 'TMPBIO', 0, 0, MESG, 2 )
        END IF
        IF ( SWITCH_FILE ) THEN
           IF ( .NOT. READ3( NNAME2, 'AVLAI', 1, 0, 0, 
     &                      AVLAIW ) ) THEN
              MESG = 'Could not read AVLAI from file "' //
     &              NNAME2( 1:LEN_TRIM( NNAME ) ) // '"'
              CALL M3EXIT( 'TMPBIO', 0, 0, MESG, 2 )
           END IF
        END IF

        I = 0
        DO  J = 1, BSPCS

            I = I + 1
            IF ( .NOT. READ3( NNAME, NORMN( I ), 1, 0, 0, 
     &                      NORNO( 1,1,J ) ) ) THEN
                 MESG = 'Could not read "' // 
     &                  NORMN( I )( 1 : LEN_TRIM( NORMN( I ) ) ) //
     &                  '" from file "' //
     &                  NNAME( 1:LEN_TRIM( NNAME ) ) // '"'
                CALL M3EXIT( 'TMPBIO', 0, 0, MESG, 2 )
            END IF

            IF ( SWITCH_FILE ) THEN
               IF ( .NOT. READ3( NNAME2, NORMN( I ), 1, 0, 0, 
     &                      NORNOW( 1,1,J ) ) ) THEN
                  MESG = 'Could not read "' // 
     &                  NORMN( I )( 1 : LEN_TRIM( NORMN( I ) ) ) //
     &                  '" from file "' //
     &                  NNAME2( 1:LEN_TRIM( NNAME2 ) ) // '"'
                  CALL M3EXIT( 'TMPBIO', 0, 0, MESG, 2 )
               END IF
            END IF

        END DO        !  end loop reading normalized NO's


C.......   Allocate memory for met
        ALLOCATE( TASFC( METNCOLS, METNROWS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'TASFC', PROGNAME )
        ALLOCATE( TSOLAR( METNCOLS, METNROWS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'TSOLAR', PROGNAME )
        ALLOCATE( PRSFC( METNCOLS, METNROWS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'PRSFC', PROGNAME )

C.......   Allocate memory for emissions 
        ALLOCATE( EMPOL( NCOLS, NROWS, BSPCS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'EMPOL', PROGNAME )
        ALLOCATE( EMISL( NCOLS, NROWS, MSPCS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'EMISL', PROGNAME )
        ALLOCATE( EMISS( NCOLS, NROWS, MSPCS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'EMISS', PROGNAME )

        LDATE = 0 
            
        METFILE = ' '
        RADFILE = ' '
        
C.........  Loop thru the number of time steps (hourly)
 
        DO  HR = 1, NSTEPS

           EMISL = 0        ! array
           EMISS = 0        ! array

           IF( JDATE .NE. LDATE ) THEN

               CALL WRDAYMSG( JDATE, MESG )               

C..........  If new date, read season switch 

               IF ( SWITCH_FILE ) THEN
                 
                  MESG = 'Reading gridded season switch data..'  
                  CALL M3MESG( MESG ) 
                  IF ( .NOT. READ3( BNAME, 'SEASON', 1, 
     &                JDATE, 0, SEASON ) ) THEN
                     MESG = 'Could not read SEASON from file ' //
     &               BNAME( 1:LEN_TRIM( BNAME ) )
                     CALL M3EXIT( 'TMPBIO', JDATE, 0, MESG, 2 )
                  END IF

                  MESG = 'Applying gridded season switch data..' 
                  CALL M3MESG( MESG )


                  DO I = 1, METNCOLS
                    DO J = 1, METNROWS

                        II = I - XOFF
                        JJ = J - YOFF

                        IF( II .LE. 0 .OR. II .GT. NCOLS .OR.
     &                      JJ .LE. 0 .OR. JJ .GT. NROWS       ) CYCLE

C........................  If switch equal to 0 use winter normalized emissions
C........................  Allow for subgrid in switch
C........................  If grid is larger subgrid, then extra cells will
C                          behave as the nearest available cell on the boundary.

                       III = MAX( MIN( I - SWXOFF, SWNCOLS ), 1 )
                       JJJ = MAX( MIN( J - SWYOFF, SWNROWS ), 1 )
                       IF ( SEASON ( III, JJJ ) .EQ. 0 ) THEN
                          DO M = 1, BSPCS - 1
                             PINES( II, JJ, M ) = PINEW( II, JJ, M ) 
                             DECDS( II, JJ, M ) = DECDW( II, JJ, M )
                             CONFS( II, JJ, M ) = CONFW( II, JJ, M )
                             LEAFS( II, JJ, M ) = LEAFW( II, JJ, M )
                             AGRCS( II, JJ, M ) = AGRCW( II, JJ, M )
                             OTHRS( II, JJ, M ) = OTHRW( II, JJ, M )
                          END DO

                          AVLAIS( II, JJ ) = AVLAIW ( II, JJ )

                          DO M = 1, BSPCS
                             NORNOS( II, JJ, M ) = NORNOW( II, JJ, M )
                          END DO

                       ELSE

                          DO M = 1, BSPCS - 1
                             PINES( II, JJ, M ) = PINE( II, JJ, M )
                             DECDS( II, JJ, M ) = DECD( II, JJ, M )
                             CONFS( II, JJ, M ) = CONF( II, JJ, M )
                             LEAFS( II, JJ, M ) = LEAF( II, JJ, M )
                             AGRCS( II, JJ, M ) = AGRC( II, JJ, M )
                             OTHRS( II, JJ, M ) = OTHR( II, JJ, M )
                          END DO

                          AVLAIS( II, JJ ) = AVLAI ( II, JJ )

                          DO M = 1, BSPCS
                             NORNOS( II, JJ, M ) = NORNO( II, JJ, M )
                          END DO

                       END IF

                    END DO
                  END DO
 
               END IF

           END IF

C.............  Write to screen because WRITE3 only writes to LDEV
           WRITE( *, 94030 ) HHMMSS( JTIME )

C.............  Open temperature file
           NEXTFILE = METLIST( METCHECK( METSTART + ( HR-1 ) ) )
           IF( NEXTFILE /= METFILE ) THEN           	
               IF( METOPEN ) THEN
                   IF( .NOT. CLOSE3( MNAME ) ) THEN
                       MESG = 'Could not close temperature file ' // 
     &                        METFILE( 1:LEN_TRIM( METFILE ) )
                       CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                   ELSE
                       METOPEN = .FALSE.
                   END IF
               END IF
               
               METFILE = NEXTFILE

C................  Set env variable
               IF( .NOT. SETENVVAR( MNAME, METFILE ) ) THEN
               	   MESG = 'Could not set logical file name ' //
     &          	  'for file ' // TRIM( METFILE )
                   CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
               END IF
               
C................  Try to open file               
               IF( .NOT. OPEN3( MNAME, FSREAD3, PROGNAME ) ) THEN
     	           MESG = 'Could not open temperature file ' // 
     &                    METFILE( 1:LEN_TRIM( METFILE ) )
     	           CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
     	       ELSE
     	       	   METOPEN = .TRUE.
     	       END IF
     	   END IF

C.............  Read temperature data
           IF ( .NOT. READ3( MNAME, TMPRNAM, 1, 
     &          MDATE, MTIME, TASFC ) ) THEN
              MESG = 'Could not read ' // TMPRNAM // 'from file ' //
     &                METFILE( 1:LEN_TRIM( METFILE ) )
              CALL M3EXIT( 'TMPBIO', MDATE, MTIME, MESG, 2 )
           END IF

C............  If necessary read solar radiation

           IF ( PARTYPE .EQ. 1 ) THEN      ! Use RADNAM

C...............  Open radiation file
              NEXTFILE = RADLIST( RADCHECK( RADSTART + ( HR-1 ) ) )
              IF( NEXTFILE /= RADFILE ) THEN
                  IF( RADOPEN ) THEN
                      IF( .NOT. CLOSE3( RNAME ) ) THEN
                          MESG = 'Could not close radiation file ' //
     &                           RADFILE( 1:LEN_TRIM( RADFILE ) )
                          CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                      ELSE
                      	  RADOPEN = .FALSE.
                      END IF
                  END IF 
                  
                  RADFILE = NEXTFILE
                  
                  IF( .NOT. SETENVVAR( RNAME, RADFILE ) ) THEN
               	      MESG = 'Could not set logical file name ' //
     &          	     'for file ' // TRIM( RADFILE )
                      CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                  END IF
                  
                  IF( .NOT. OPEN3( RNAME, FSREAD3, PROGNAME ) ) THEN
     	              MESG = 'Could not open radiation file ' // 
     &         	              RADFILE( 1:LEN_TRIM( RADFILE ) )
     	              CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
     	          ELSE
     	              RADOPEN = .TRUE.
     	          END IF
     	      END IF

              IF ( .NOT. READ3( RNAME, RADNAM, ALLAYS3, MDATE,
     &             MTIME, TSOLAR(1,1) ) ) THEN
                 MESG = 'Could not read ' // RADNAM // 'from file ' //
     &                RADFILE( 1:LEN_TRIM( RADFILE ) ) 
                 CALL M3EXIT( 'TMPBIO', MDATE, MTIME, MESG, 2 )
              END IF

C...............  Convert shortwave radiation to PAR  

              DO I = 1, METNCOLS
               DO J = 1, METNROWS
                    TSOLAR( I, J ) = TSOLAR( I, J ) * SOL2PAR ! Calc. PAR
               END DO
              END DO

C...............  Calculate non-speciated emissions
C...............     must pass met date and time here

              IF ( SWITCH_FILE ) THEN
  
                 CALL HRBIOS( MDATE, MTIME, METNCOLS, METNROWS, 
     &                   PINES, DECDS, CONFS, AGRCS, LEAFS, OTHRS,
     &                   AVLAIS, NORNOS, TASFC, TSOLAR, EMPOL )
          
              ELSE 
 
                 CALL HRBIOS( MDATE, MTIME, METNCOLS, METNROWS,
     &                   PINE, DECD, CONF, AGRC, LEAF, OTHR, AVLAI,
     &                   NORNO, TASFC, TSOLAR, EMPOL )

              END IF

           ELSE        ! Use clouds or clear skies to calculate PAR

C..............  Get name of surface pressure variable to use

              MESG = 'Variable name for surface pressure'
              CALL ENVSTR( 'PRES_VAR', MESG, 'PRSFC', PRESNAM, IOS )

C..............  Read surface pressure data from temperature file

              IF ( .NOT. READ3( MNAME, PRESNAM, 1, 
     &                        MDATE, MTIME, PRSFC ) ) THEN
                  MESG = 'Could not read '// PRESNAM // 'from file "'//
     &                   METFILE( 1:LEN_TRIM( METFILE ) ) // '"'
                  CALL M3EXIT( 'TMPBIO', MDATE, MTIME, MESG, 2 )
              END IF

C...............  convert to millibars

              DO  C = 1, METNCOLS
              DO  R = 1, METNROWS
                       PRSFC( C, R ) = PRSFC( C, R ) * 0.010  ! Pa to mb
              END DO
              END DO

C............. Calculate non-speciated emissions
C............. must pass met date and time here

              IF ( SWITCH_FILE ) THEN

                 CALL HRBIO( MDATE, MTIME, METNCOLS, METNROWS, LAT, LON,  
     &                  PRSFC, TASFC, PINES, DECDS, CONFS, AGRCS, 
     &                  LEAFS, OTHRS, AVLAIS, NORNOS, PARTYPE, EMPOL )
              ELSE
 
                 CALL HRBIO( MDATE, MTIME, METNCOLS, METNROWS, LAT, LON,
     &                  PRSFC, TASFC, PINE, DECD, CONF, AGRC, LEAF,
     &                  OTHR, AVLAI, NORNO, PARTYPE, EMPOL )

              END IF

           END IF 

C............. Speciate emissions

           DO I = 1, NCOLS
             DO J = 1, NROWS
               DO L = 1, MSPCS
                 DO K = 1, BSPCS
                   EMISL( I, J, L  ) = EMISL( I ,J, L ) +
     &                             EMPOL( I, J, K ) * MLFAC( L, K )
                   EMISS( I, J, L  ) = EMISS( I, J, L ) +
     &                             EMPOL( I, J, K ) * MSFAC( L, K ) 
                 END DO
               END DO
             END DO
           END DO

C.............  Write out speciated emissions    

           IF ( .NOT. WRITE3( ENAME, 'ALL', 
     &                        JDATE, JTIME, EMISL ) ) THEN
               CALL M3EXIT( 'TMPBIO', JDATE, JTIME, 
     &                      'Error writing BIO OUTPUT file' , 2 )
           END IF                              !  if write3 failed


           IF ( .NOT. WRITE3( SNAME, 'ALL',
     &                        JDATE, JTIME, EMISS ) ) THEN
               CALL M3EXIT( 'TMPBIO', JDATE, JTIME,
     &                      'Error writing BIO OUTPUT file' , 2 )
           END IF                              !  if write3 failed


C.............. Next time step

           LDATE = JDATE
           CALL NEXTIME( JDATE, JTIME, 10000 )

           CALL NEXTIME( MDATE, MTIME, 10000 ) 

      END DO                !  end loop on hours HR


C.........   End of program:

        CALL M3EXIT( 'TMPBIO', 0, 0, ' ', 0 )

C******************  FORMAT  STATEMENTS   ******************************

C...........   Informational (LOG) message formats... 92xxx

92000   FORMAT ( 5X , A )

C...........   Internal buffering formats............ 94xxx

94000   FORMAT( I2.2 )

94030   FORMAT( 8X, 'at time ', A8 )

94052   FORMAT( A, 1X, A, 1X, I6.6, 1X,
     &          A, 1X, I3.3, 1X, A, 1X, I3.3, 1X, A   )

        END PROGRAM TMPBIO

