
        PROGRAM TEMPORAL

C***********************************************************************
C  program body starts at line 213
C
C  DESCRIPTION:
C    This program computes the hourly emissions data from inventory emissions 
C    and/or activity and emission factor data. It can read average-inventory,
C    day-specific and hour-specific emissions and activity data.
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
C*************************************************************************

C.........  MODULES for public variables
C.........  This module contains the inventory arrays
        USE MODSOURC

C.........  This module contains the temporal cross-reference tables
        USE MODXREF

C.........  This module contains the temporal profile tables
        USE MODTMPRL

C.........  This module contains emission factor tables and related
        USE MODEMFAC

C.........  This module contains data for day- and hour-specific data
        USE MODDAYHR

C...........   This module is the derived meteorology data for emission factors
        USE MODMET

C.........  This module contains the lists of unique source characteristics
        USE MODLISTS

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
        INTEGER         FINDC
        INTEGER         GETDATE
        INTEGER         GETFLINE
        INTEGER         GETNUM
        INTEGER         INDEX1
        CHARACTER*14    MMDDYY
        INTEGER         RDTPROF

        EXTERNAL    CHKEMEPI, CRLF, ENVINT, ENVYN, FINDC, GETDATE, 
     &              GETFLINE, GETNUM, INDEX1, MMDDYY, RDTPROF
                        
C.........  LOCAL PARAMETERS and their descriptions:

        CHARACTER*50, PARAMETER :: CVSW = '$Name$'  ! CVS revision tag

C.........  Point sources work and output arrays
        REAL   , ALLOCATABLE :: EMAC ( :,: ) !  inven emissions or activities
        REAL   , ALLOCATABLE :: EMACV( :,: ) !  day-adjst emis or activities
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

        INTEGER                 NPELV        ! optional elevated source-count
        INTEGER, ALLOCATABLE :: INDXE( : )   ! SMOKE source IDs
        REAL   , ALLOCATABLE :: EMISE( : )   ! elevated source emissions

C...........   Gridded meteorology data
        REAL   , ALLOCATABLE :: TA( : )     ! one layer of temperature

C...........   Ungridding Matrix
        INTEGER, ALLOCATABLE :: UMAT( : )   ! contiguous ungridding matrix

C.........  Names of pollutants and activities associated with output variables
        CHARACTER(LEN=IOVLEN3), ALLOCATABLE:: ALLIN( : ) 

C.........  Reshaped input variables and output variables
        INTEGER         NGRP                ! no. of pol/emis-types groups 
        INTEGER         NGSZ                ! no. of pols/emis-types per group 
        CHARACTER(LEN=IOVLEN3), ALLOCATABLE:: ALLIN2D( :,: ) 
        CHARACTER(LEN=IOVLEN3), ALLOCATABLE:: EANAM2D( :,: ) 

C...........   Logical names and unit numbers

        INTEGER      :: CDEV = 0!  unit number for region codes file
        INTEGER      :: HDEV = 0!  unit number for holidays file
        INTEGER         LDEV    !  unit number for log file
        INTEGER      :: MDEV = 0!  unit number for mobile codes file
        INTEGER         RDEV    !  unit number for temporal profile file
        INTEGER         SDEV    !  unit number for ASCII inventory file
        INTEGER         TDEV    !  unit number for emission processes file
        INTEGER         FDEV    !  unit number for EF xref file
        INTEGER         XDEV    !  unit no. for cross-reference file

        CHARACTER*16 :: ANAME = ' '    !  logical name for ASCII inven input 
        CHARACTER*16 :: DNAME = 'NONE' !  day-specific  input file, or "NONE"
        CHARACTER*16 :: ENAME = ' '    !  logical name for I/O API inven input
        CHARACTER*16 :: FNAME = ' '    !  non-diurnal EF file
        CHARACTER*16 :: GNAME = ' '    !  ungridding matrix
        CHARACTER*16 :: HNAME = 'NONE' !  hour-specific input file, or "NONE"
        CHARACTER*16 :: MNAME = ' '    !  surface temperature file
        CHARACTER*16 :: NNAME = ' '    !  diurnal EF file
c        CHARACTER*16 :: SNAME = ' '    !  PSIs per source file
        CHARACTER*16 :: TNAME = ' '    !  timestepped (low-level) output file
        CHARACTER*16 :: WNAME = ' '    !  ungridded min/max temperatures

C...........   Other local variables

        INTEGER         I, J, K, L, L1, L2, N, S, T

        INTEGER         IOS, IOS1, IOS2, IOS3, IOS4 ! i/o status
        INTEGER         IOS6, IOS7, IOS8      ! i/o status
        INTEGER         EDATE, ETIME        ! ending Julian date and time
        INTEGER         ENLEN               ! length of ENAME string
        INTEGER         JDATE, JTIME        ! Julian date and time
        INTEGER         MDATE, MTIME        ! current met date/time
        INTEGER      :: MLDATE = 0          ! gridded met previous date
        INTEGER      :: MSDATE = 0          ! gridded met start date
        INTEGER      :: MSTIME = 0          ! gridded met start time
        INTEGER         NGRID               ! no. grid cells
        INTEGER         NINVARR             ! no. inventory variables to read
        INTEGER         NLINE               ! tmp number of lines in ASCII file 
        INTEGER         NMAJOR              ! no. major sources
        INTEGER         NMATX               ! size of ungridding matrix
        INTEGER         NPING               ! no. ping sources
        INTEGER         NSTEPS              ! number of output time steps
        INTEGER      :: PYEAR = 0           ! projected year
        INTEGER         SDATE, STIME        ! starting Julian date and time
        INTEGER         TNLEN               ! length of TNAME string
        INTEGER         TSTEP               ! output time step
        INTEGER         TZONE               ! output-file time zone
        INTEGER      :: WDATE               ! source min/max tmpr current date
        INTEGER      :: WEDATE = 0          ! source min/max tmpr end date
        INTEGER      :: WSDATE = 0          ! source min/max tmpr start date
        INTEGER      :: WSTIME = 0          ! source min/max tmpr start time
        INTEGER      :: WTIME               ! source min/max tmpr current time

        REAL            RTMP                ! tmp float

        LOGICAL         DFLAG   !  true: day-specific  file available
        LOGICAL      :: EFLAG = .FALSE.  !  error-flag
        LOGICAL      :: EFLAG2= .FALSE.  !  error-flag (2)
        LOGICAL         HFLAG   !  true: hour-specific file available
        LOGICAL         MFLAG   !  true: mobile codes file available
        LOGICAL         NFLAG   !  true: use all uniform temporal profiles
        LOGICAL         WFLAG   !  true: write QA on current time step
        LOGICAL         WTINDP  !  true: MINMAXT file is time-independent

        CHARACTER*8              TREFFMT ! tmprl x-ref format (SOURCE|STANDARD)
        CHARACTER*14             DTBUF   ! buffer for MMDDYY
        CHARACTER*20             MODELNAM! emission factor model name
        CHARACTER*300            MESG    ! buffer for M3EXIT() messages
        CHARACTER(LEN=IOVLEN3)   CBUF    ! pollutant name temporary buffer 
        CHARACTER(LEN=IOVLEN3):: TVARNAME = ' ' ! temperature variable name

        CHARACTER*16 :: PROGNAME = 'TEMPORAL' ! program name

C***********************************************************************
C   begin body of program TEMPORAL

        LDEV = INIT3()

C.........  Write out copywrite, version, web address, header info, and prompt
C           to continue running the program.
        CALL INITEM( LDEV, CVSW, PROGNAME )

C.........  Obtain settings from the environment...

C.........  Get the time zone for output of the emissions
        TZONE = ENVINT( 'OUTZONE', 'Output time zone', 0, IOS )

C.........  Get environment variable that overrides temporal profiles and 
C               uses only uniform profiles.
        NFLAG = ENVYN( 'UNIFORM_TPROF_YN', MESG, .FALSE., IOS )

C.........  Get environment variable that indicates if the MINMAXT file is
C           time independent (true) or not
        MESG   = 'MINMAXT is time-independent file? '  
        WTINDP = ENVYN( 'MINMAX_TINDP_YN', MESG, .FALSE., IOS )

C.........  Set source category based on environment variable setting
        CALL GETCTGRY

C.........  Get the name of the emission factor model to use for one run
        IF ( CATEGORY .EQ. 'MOBILE' ) THEN
            MESG = 'Emission factor model'
            CALL ENVSTR( 'SMK_EF_MODEL', MESG, 'MOBILE5', MODELNAM, IOS)
        ELSE
            MODELNAM = ' '
        END IF

C.........  Get inventory file names given source category
        CALL GETINAME( CATEGORY, ENAME, ANAME )

C.........  Prompt for and open input files
C.........  Also, store source-category specific information in the MODINFO 
C           module.
C.........  Also, determine name of temperature variable in the temperature
C           files (if any)
C.........  Also, compare min/max temperature settings from any files that have
C           them and populate the valid temperature arrays
        CALL OPENTMPIN( MODELNAM, NFLAG, ENAME, ANAME, DNAME, HNAME, 
     &                  FNAME, NNAME, MNAME, GNAME, WNAME, TVARNAME, 
     &                  SDEV, XDEV, RDEV, FDEV, CDEV, HDEV, TDEV, 
     &                  MDEV, PYEAR )

C.........  Determine status of some files for program control purposes
        DFLAG = ( DNAME .NE. 'NONE' )  ! Day-specific emissions
        HFLAG = ( HNAME .NE. 'NONE' )  ! Hour-specific emissions
        MFLAG = ( MDEV .NE. 0 )        ! Use mobile codes file

C.........  Get length of inventory file name
        ENLEN = LEN_TRIM( ENAME )

C.........  Get episode settings from the Models-3 environment variables
        SDATE  = 0
        STIME  = 0
        NSTEPS = 1
        TSTEP  = 10000  
        CALL GETM3EPI( TZONE, SDATE, STIME, NSTEPS )

C.........  Compare base year with episode and warn if not consistent
        IF( SDATE / 1000 .NE. BYEAR ) THEN

            WRITE( MESG,94010 ) 'WARNING: Inventory base year ', BYEAR, 
     &             'is inconsistent with year ' // CRLF() // BLANK10 //
     &             'of episode start date', SDATE/1000
            CALL M3MSG2( MESG )

        ENDIF

C.........  Give a note if running for a projected year
        IF( PYEAR .GT. 0 ) THEN

            WRITE( MESG,94010 ) 'NOTE: Emissions based on projected '//
     &             'year', PYEAR
            CALL M3MSG2( MESG )

        END IF

C.........  Calculate the ending date and time
        EDATE = SDATE
        ETIME = STIME
        CALL NEXTIME( EDATE, ETIME, NSTEPS * 10000 )

C.........  Set up gridded and min/max temperature files dates and times and
C           retrieve the number of grid cells
        IF( MNAME .NE. ' ' ) THEN
            CALL PRETMPR( MNAME, WNAME, TZONE, TSTEP, SDATE, STIME, 
     &                    NSTEPS, MSDATE, MSTIME, WSDATE, WSTIME, 
     &                    WEDATE, NGRID)
        END IF

C.........  For day-specific data input...
        IF( DFLAG ) THEN

C.............  Get header description of day-specific input file
            IF( .NOT. DESC3( DNAME ) ) THEN
                CALL M3EXIT( PROGNAME, 0, 0, 
     &                       'Could not get description of file "' 
     &                       // DNAME( 1:LEN_TRIM( DNAME ) ) // '"', 2 )
            END IF

C.............  Allocate memory for pollutant pointer
            ALLOCATE( DYPNAM( NVARS3D ), STAT=IOS )
            CALL CHECKMEM( IOS, 'DYPNAM', PROGNAME )
            ALLOCATE( DYPDSC( NVARS3D ), STAT=IOS )
            CALL CHECKMEM( IOS, 'DYPDSC', PROGNAME )
            DYPNAM = ' '  ! array
            DYPDSC = ' '  ! array

C.............  Set day-specific file dates, check dates, and report problems
            CALL PDSETUP( DNAME, SDATE, STIME, EDATE, ETIME, TZONE,  
     &                    NIPPA, EANAM, NDYPOA, NDYSRC, EFLAG, DYPNAM,
     &                    DYPDSC )

        ENDIF

C.........  Allocate memory for reading day-specific emissions data
C.........  NDYSRC is initialized to zero in case DFLAG is false
        ALLOCATE( INDXD( NDYSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'INDXD', PROGNAME )
        ALLOCATE( EMACD( NDYSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'EMACD', PROGNAME )

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
            ALLOCATE( HRPDSC( NVARS3D ), STAT=IOS )
            CALL CHECKMEM( IOS, 'HRPDSC', PROGNAME )
            HRPNAM = ' '  ! array
            HRPDSC = ' '  ! array

C.............  Set day-specific file dates, check dates, and report problems
            CALL PDSETUP( HNAME, SDATE, STIME, EDATE, ETIME, TZONE,  
     &                    NIPPA, EANAM, NHRPOA, NHRSRC, EFLAG2, HRPNAM,
     &                    HRPDSC )

        ENDIF

C.........  Allocate memory for reading hour-specific emissions data
C.........  NHRSRC is initialized to 0 in case HFLAG is false
        ALLOCATE( INDXH( NHRSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'INDXH', PROGNAME )
        ALLOCATE( EMACH( NHRSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'EMACH', PROGNAME )

        IF( EFLAG .OR. EFLAG2 ) THEN
            MESG = 'Problem with day- or hour-specific inputs'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        ENDIF

C.........  Set inventory variables to read for all source categories
        IVARNAMS( 1 ) = 'IFIP'
        IVARNAMS( 2 ) = 'TZONES'
        IVARNAMS( 3 ) = 'TPFLAG'
        IVARNAMS( 4 ) = 'CSCC'
        IVARNAMS( 5 ) = 'CSOURC'

C.........  Set inventory variables to read for specific source categories
        IF( CATEGORY .EQ. 'AREA' ) THEN
            NINVARR = 5

        ELSE IF( CATEGORY .EQ. 'MOBILE' ) THEN
            NINVARR = 9
            IVARNAMS( 6 ) = 'IRCLAS'
            IVARNAMS( 7 ) = 'IVTYPE'
            IVARNAMS( 8 ) = 'CLINK'
            IVARNAMS( 9 ) = 'CVTYPE'

        ELSE IF( CATEGORY .EQ. 'POINT' ) THEN
            NINVARR = 5

        END IF

C.........  Allocate memory for and read in required inventory characteristics
        CALL RDINVCHR( CATEGORY, ENAME, SDEV, NSRC, NINVARR, IVARNAMS )

C.........  Build unique lists of SCCs per SIC from the inventory arrays
        CALL GENUSLST

C.........  Read special files...

C.........  Read region codes file
        CALL RDSTCY( CDEV, NINVIFIP, INVIFIP )

C.........  Populate filter for sources that use daylight time
        CALL SETDAYLT

C.........  Read holidays file
        CALL RDHDAYS( HDEV, SDATE, EDATE )

C.........  When mobile codes file is being used read mobile codes file
        IF( MFLAG ) CALL RDMVINFO( MDEV )

C.........  Perform steps needed for using activities and emission factors
C.........  NOTE - the CHRT* variables that are part of the MODXREF module 
C           will be used by RDEFXREF and ASGNPSI and done with before these
C           variables are again needed by RDTREF and ASGNTPRO.

        IF( NIACT .GT. 0 ) THEN

C.............  Read header of ungridding matrix...
            IF( .NOT. DESC3( GNAME ) ) THEN
        	MESG = 'Could not get description for file ' // GNAME
        	CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF

C.............  Store number of ungridding factors
            NMATX = NCOLS3D

C.............  Allocate memory for ungridding matrix, gridded temperature, 
C               and PSI index matching
            ALLOCATE( UMAT( NSRC + 2*NMATX ), STAT=IOS )
            CALL CHECKMEM( IOS, 'UMAT', PROGNAME )
            ALLOCATE( TA( NGRID ), STAT=IOS )
            CALL CHECKMEM( IOS, 'TA', PROGNAME )
            ALLOCATE( TASRC( NSRC ), STAT=IOS )
            CALL CHECKMEM( IOS, 'TASRC', PROGNAME )
	    ALLOCATE( TKMIN( NSRC ), STAT=IOS )
	    CALL CHECKMEM( IOS, 'TKMIN', PROGNAME )
	    ALLOCATE( TKMAX( NSRC ), STAT=IOS )
	    CALL CHECKMEM( IOS, 'TKMAX', PROGNAME )
	    ALLOCATE( METIDX( NSRC,4 ), STAT=IOS )
	    CALL CHECKMEM( IOS, 'METIDX', PROGNAME )
 	    ALLOCATE( EFSIDX( NSRC,NIACT ), STAT=IOS )
	    CALL CHECKMEM( IOS, 'EFSIDX', PROGNAME )

C.............  Read ungridding matrix
            CALL RDUMAT( GNAME, NSRC, NMATX, NMATX, UMAT( 1 ),
     &                   UMAT( NSRC+1 ), UMAT( NSRC+NMATX+1 )  )

C.............  Read the cross-reference for assigning the PSIs
C               to the sources. This will populate parts of the MODXREF module.
C.............  Also create the list of unique parameter scheme indices.
            CALL RDEFXREF( FDEV, .TRUE. )

C.............  Read emission processes file.  Populate array in MODEMFAC.
            CALL RDEPROC( TDEV )

C.............   Map parameter scheme indexes (PSIs) onto sources for all actvty
            CALL ASGNPSI( NIACT, ACTVTY, NETYPE )

C.............  Loop through activities and...
C.............  NOTE - this is not fully implemented for multiple activities. 
C               To do this, the data structures and RDEFACS will need to be 
C               updated. Also, the variable names in the emission factor file
C               are not truly supporting 16-character pollutant and 
C               emission process names, because it is only set up for MOBILE5
            DO I = 1, NIACT

C.................  Skip activities that do not have emissions types
                IF( NETYPE( I ) .LE. 0 ) CYCLE            

C.................  Read emission factor tables for each activity and get the 
C                   emission factor variable names.
        	CALL RDEFACS( MODELNAM, FNAME, NNAME, ACTVTY( I ), 
     &                        UMAT(1) )

            END DO

        END IF

C.........  Determine all of the variables to be output based on the 
C           activities and input pollutants.  
C.........  NOTE - Uses NETYPE, EMTACTV, and EMTNAM, and sets units, and units
C           conversion factors for all pollutants and activities
        CALL TMNAMUNT

C.........  Reset the number of all output variables as the number of pollutants
C           and emission types, instead of the number of pollutants and 
C           activities
        NIPPA = NIPOL
        DO I = 1, NIACT
            NIPPA = NIPPA + NETYPE( I )
        END DO

C.........  Allocate memory for I/O pol names, activities, & emission types
C.........  Will be resetting EANAM to include the emission types instead
C           of the activities
        DEALLOCATE( EANAM )
        ALLOCATE( EANAM( NIPPA ), STAT=IOS )
        CALL CHECKMEM( IOS, 'EANAM', PROGNAME )
        ALLOCATE( ALLIN( NIPPA ), STAT=IOS )
        CALL CHECKMEM( IOS, 'ALLIN', PROGNAME )

C.........  Create 1-d arrays of I/O pol names
C.........  If pollutant is created by one or more activities, then skip it.
        N = 0
        DO I = 1, NIPOL

C.............  Look for pollutant in list of pollutants created by
C               activities
            J = INDEX1( EINAM( I ), NEPOL, EMTPOL )

C.............  If pollutant created by activity, skip from this list
            IF( J .LE. 0 ) THEN
                N = N + 1
                ALLIN( N ) = EAREAD( I )
                EANAM( N ) = EINAM ( I )
            END IF

        END DO

C.........  Add activities, & emission types to read and output lists
        J = NIPOL + 1
        DO I = 1, NIACT

            K = NETYPE( I )  ! Number of emission types

C.............  If any emissions types associated with this activity, store them
            IF ( K .GT. 0 ) THEN
                N = N + K
                ALLIN( J:J+K-1 ) = ACTVTY( I )
                EANAM( J:J+K-1 ) = EMTNAM( 1:K, I )
            END IF
            J = J + K

        END DO

C.........  Reset number of pollutants and emission types based on those used
        NIPPA = N

C.........  Read temporal-profile cross-reference file and put into tables
C.........  Only read entries for pollutants that are in the inventory.
C.........  Only read if not using uniform temporal profiles
        IF( .NOT. NFLAG ) CALL RDTREF( XDEV, TREFFMT )

C.........  Read temporal-profiles file:  4 parts (monthly, weekly, 
C           weekday diurnal, and weekend diurnal)
        CALL M3MSG2( 'Reading temporal profiles file...' )

        NMON = RDTPROF( RDEV, 'MONTHLY', NFLAG )
        NWEK = RDTPROF( RDEV, 'WEEKLY' , NFLAG )
        NHRL = RDTPROF( RDEV, 'DIURNAL', NFLAG )

C.........  Adjust temporal profiles for use in generating temporal emissions
C.........  NOTE - All variables are passed by modules.
        CALL NORMTPRO

C.........  It is important that all major arrays must be allocated by this 
C           point because the next memory allocation step is going to pick a
C           data structure that will fit within the limits of the host.

C.........  Allocate memory, but allow flexibility in memory allocation
C           for second dimension.
C.........  The second dimension (the number of pollutants and emission types)
C            can be different depending on the memory available.
C.........  To determine the approproate size, first attempt to allocate memory
C           for all pollutants & emission types to start, and if this fails,
C           divide pollutants into even groups and try again.

        NGSZ = NIPPA            ! No. of pollutant & emis types in each group
        NGRP = 1               ! Number of groups
        DO

            ALLOCATE( TMAT ( NSRC, NGSZ, 24 ), STAT=IOS1 )
            ALLOCATE( MDEX ( NSRC, NGSZ )    , STAT=IOS2 )
            ALLOCATE( WDEX ( NSRC, NGSZ )    , STAT=IOS3 )
            ALLOCATE( DDEX ( NSRC, NGSZ )    , STAT=IOS4 )
            ALLOCATE( EMAC ( NSRC, NGSZ )    , STAT=IOS6 )
            ALLOCATE( EMACV( NSRC, NGSZ )    , STAT=IOS7 )
            ALLOCATE( EMIST( NSRC, NGSZ )    , STAT=IOS8 )

            IF( IOS1 .GT. 0 .OR. IOS2 .GT. 0 .OR. IOS3 .GT. 0 .OR.
     &          IOS4 .GT. 0 .OR. IOS6 .GT. 0 .OR.
     &          IOS7 .GT. 0 .OR. IOS8 .GT. 0 ) THEN

                IF( NGSZ .EQ. 1 ) THEN
                    J = 8 * NSRC * 31    ! Assume 8-byte reals
                    WRITE( MESG,94010 ) 
     &                'Insufficient memory to run program.' //
     &                CRLF() // BLANK5 // 'Could not allocate ' // 
     &                'pollutant-dependent block of', J, 'bytes.'
                    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                END IF

                NGRP = NGRP + 1
                NGSZ = NGSZ / NGRP + ( NIPPA - NGSZ * NGRP )

                DEALLOCATE( TMAT, MDEX, WDEX, DDEX, EMAC, 
     &                      EMACV, EMIST )

            ELSE
                EXIT

            END IF

        END DO
        
C.........  Allocate a few small arrays based on the size of the groups
C.........  NOTE that this has a small potential for a problem if these little
C           arrays exceed the total memory limit.
        ALLOCATE( ALLIN2D( NGSZ, NGRP ), STAT=IOS )
        CALL CHECKMEM( IOS, 'ALLIN2D', PROGNAME )
        ALLOCATE( EANAM2D( NGSZ, NGRP ), STAT=IOS )
        CALL CHECKMEM( IOS, 'EANAM2D', PROGNAME )
        ALLOCATE( LDSPOA( NGSZ ), STAT=IOS )
        CALL CHECKMEM( IOS, 'LDSPOA', PROGNAME )
        ALLOCATE( LHSPOA( NGSZ ), STAT=IOS )
        CALL CHECKMEM( IOS, 'LHSPOA', PROGNAME )
        ALLOCATE( LHPROF( NGSZ ), STAT=IOS )
        CALL CHECKMEM( IOS, 'LHPROF', PROGNAME )

C.........  Create 2-d arrays of I/O pol names, activities, & emission types 
        ALLIN2D  = ' '
        ALLIN2D  = RESHAPE( ALLIN, (/ NGSZ, NGRP /) )
        EANAM2D = ' '
        EANAM2D = RESHAPE( EANAM, (/ NGSZ, NGRP /) )

C.........  Set up and open I/O API output file(s) ...
        CALL OPENTMP( ENAME, SDATE, STIME, TSTEP, TZONE, NPELV,
     &                TNAME )

        TNLEN = LEN_TRIM( TNAME )

C.........  Loop through pollutant/emission-type groups
        DO N = 1, NGRP

C.............  Write message stating the pols/emission-types being processed
            CALL POLMESG( NGSZ, EANAM2D( 1,N ) )

C.............  Set up logical arrays that indicate which pollutants/activities
C               are day-specific and which are hour-specific.
C.............  Also set flag for which hour-specific pollutants/activities
C               are actually diurnal profiles instead of emissions
            LDSPOA = .FALSE.   ! array
            DO I = 1, NDYPOA
                J = INDEX1( DYPNAM( I ), NGSZ,  EANAM2D( 1,N ) )
                LDSPOA( J ) = .TRUE.
            END DO

            LHSPOA = .FALSE.   ! array
            LHPROF = .FALSE.   ! array
            DO I = 1, NHRPOA
                J = INDEX1( HRPNAM( I ), NGSZ,  EANAM2D( 1,N ) )
                LHSPOA( J ) = .TRUE.

                CALL UPCASE( HRPDSC( I ) )
                K = INDEX( HRPDSC( I ), 'PROFILE' )
                IF( K .GT. 0 ) LHPROF( J ) = .TRUE.
            END DO

C.............  Initialize emissions, activities, and other arrays for this
C               pollutant/emission-type group
            TMAT  = 0.
            MDEX  = IMISS3
            WDEX  = IMISS3
            DDEX  = IMISS3
            EMAC  = 0.
            EMACV = 0.
            EMIST = 0.

C.............  Assign temporal profiles by source and pollutant
            CALL M3MSG2( 'Assigning temporal profiles to sources...' )

C.............  If using uniform profiles, set all temporal profile number
C               to 1; otherwise, assign profiles with cross-reference info
            IF( NFLAG ) THEN
                MDEX = 1
        	WDEX = 1
        	DDEX = 1

            ELSE
                CALL ASGNTPRO( NGSZ, EANAM2D( 1,N ), TREFFMT )

            END IF

C.............  Read in pollutant emissions or activities from inventory for 
C               current group
            DO I = 1, NGSZ

                CBUF = ALLIN2D( I,N )
                L1   = LEN_TRIM( CBUF )

                IF( .NOT. READ3( ENAME, CBUF, ALLAYS3, 0, 0, 
     &                           EMAC( 1,I )                ) ) THEN
                    EFLAG = .TRUE.
                    MESG = 'Error reading "' // CBUF( 1:L1 ) //
     &                     '" from file "' // ENAME( 1:ENLEN ) // '."'
                    CALL M3MSG2( MESG )

                END IF

C.................  If there are any missing values in the data, give an
C                   error to avoid problems in genhemis routine
                RTMP = MINVAL( EMAC( 1:NSRC,I ) )
                IF( RTMP .LT. 0 ) THEN
                    EFLAG = .TRUE.
                    MESG = 'ERROR: Missing or zero emission value(s) '//
     &                     'in inventory for "' // CBUF( 1:L1 ) // '".'
                    CALL M3MSG2( MESG )
                END IF

C.................  If pollutant name is ozone-season-based, remove the
C                   prefix from the input pollutant name
                K = INDEX1( CBUF, NIPPA, EAREAD )
                J = INDEX( CBUF, OZNSEART )
                IF( J .GT. 0 ) THEN
                    CBUF = CBUF( CPRTLEN3+1:L1 )
                    ALLIN2D( I,N ) = CBUF
                    EAREAD ( K )   = CBUF
                END IF

            END DO

C.............  Abort if error found
            IF( EFLAG ) THEN
                MESG = 'Problem with input data.'
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF

C.............  For each time step and pollutant or emission type in current 
C               group, generate hourly emissions, write elevated emissions 
C               file (if any), and write layer-1 emissions file (or all data).
            JDATE = SDATE
            JTIME = STIME
            MDATE = MSDATE
            MTIME = MSTIME

            IF ( WTINDP ) THEN

               IF( .NOT. READ3( WNAME, 'TKMIN', 1,
     &                          0, 0, TKMIN ) ) THEN
                 MESG = 'Could not read TKMIN from ' // WNAME
                 CALL M3EXIT( PROGNAME,0,0,MESG,2 )
               END IF

              IF( .NOT. READ3( WNAME, 'TKMAX', 1,
     &                         0, 0, TKMAX ) ) THEN
                MESG = 'Could not read TKMAX from ' // WNAME
                CALL M3EXIT( PROGNAME,0,0,MESG,2 )
              END IF

              IF( .NOT. READ3( WNAME, 'TMMI', ALLAYS3,
     &                         0, 0, METIDX ) ) THEN
                MESG = 'Could not read TMMI from ' // WNAME
                CALL M3EXIT( PROGNAME,0 ,0, MESG,2 )
              END IF

            ELSE

              WDATE = WSDATE
              WTIME = WSTIME

            ENDIF

            DO T = 1, NSTEPS

C.................  When there are activity data, assume that there is a 
C                   temperature dependence, as with Models-3
                IF( NIACT .GT. 0 ) THEN

C.....................  Read gridded temperature data and convert to
C                       source-based temperatures using the ungridding matrix
                    IF ( .NOT. READ3( MNAME, TVARNAME, 1, 
     &                                MDATE, MTIME, TA    ) ) THEN

                	L = LEN_TRIM( TVARNAME )
                	MESG = 'Could not read ' // TVARNAME( 1:L ) //
     &                         'from ' // TNAME 
                        CALL M3EXIT( PROGNAME, MDATE, MTIME, MESG, 2 )

                    END IF

                    CALL APPLUMAT( NSRC, NMATX, TA, UMAT(1), 
     &                             UMAT( NSRC+1 ), UMAT( NSRC+NMATX+1 ),
     &                             TASRC )

C.....................  Read source-based min/max temperatures for the current
C                       day when there is a new day in GMT (met ) time zone
                    IF( MDATE .NE. MLDATE .AND. .NOT. WTINDP ) THEN

                	IF( .NOT. READ3( WNAME, 'TKMIN', 1, 
     &                                   WDATE, WTIME, TKMIN ) ) THEN
                	    MESG = 'Could not read TKMIN from ' // WNAME
                            CALL M3EXIT( PROGNAME,WDATE,WTIME,MESG,2 )
                	END IF

                	IF( .NOT. READ3( WNAME, 'TKMAX', 1, 
     &                                   WDATE, WTIME, TKMAX ) ) THEN
                	    MESG = 'Could not read TKMAX from ' // WNAME
                            CALL M3EXIT( PROGNAME,WDATE,WTIME,MESG,2 )
                	END IF

                	IF( .NOT. READ3( WNAME, 'TMMI', ALLAYS3, 
     &                                   WDATE, WTIME, METIDX ) ) THEN
                	    MESG = 'Could not read TMMI from ' // WNAME
                            CALL M3EXIT( PROGNAME,WDATE,WTIME,MESG,2 )
                	END IF

                    END IF

                END IF

C.................  Generate hourly emissions for current hour
                CALL GENHEMIS( NGSZ, JDATE, JTIME, TZONE, DNAME, HNAME, 
     &                         ALLIN2D( 1,N ), EANAM2D( 1,N ), 
     &                         EMAC, EMACV, TMAT, EMIST )

C.................  Loop through pollutants/emission-types in this group
                DO I = 1, NGSZ

                    CBUF = EANAM2D( I,N )

C.....................  Write hourly emissions to I/O API NetCDF file
                    IF( .NOT. WRITE3( TNAME, CBUF, JDATE, JTIME, 
     &                                EMIST( 1,I )              ) ) THEN

                        L = LEN_TRIM( CBUF )
                        MESG = 'Could not write "' // CBUF( 1:L ) // 
     &                         '" to file "' // TNAME( 1:TNLEN ) // '."'
                        CALL M3EXIT( PROGNAME, JDATE, JTIME, MESG, 2 )

                    END IF

                END DO  ! End loop on pollutants/emission-types I in this group

C.................  Advance the output date/time by one time step
                CALL NEXTIME( JDATE, JTIME, TSTEP )

C.................  Advance the met date/time by one time step, saving the
C                   old date
                MLDATE = MDATE
                CALL NEXTIME( MDATE, MTIME, TSTEP )

C.................  Advance the min/max date/time by one day when the GMT date
C                   of the meterology has changed
C.................  Only advance time if the min/max file date is less than
C                   the ending date. Pretmpr has already given a warning about
C                   the mismatch.
                IF( MDATE .NE. MLDATE .AND. 
     &              WDATE .LT. WEDATE      ) THEN

                    CALL NEXTIME( WDATE, WTIME, 240000 )

                END IF

C.................  Call QA report routine
c               WFLAG = ( T .EQ. NSTEPS )
c               CALL QATMPR( LDEV, NGSZ, T, JDATE, JTIME, WFLAG, 
c    &                       EANAM2D( 1,N ), EMAC )

            END DO      ! End loop on time steps T

        END DO          ! End loop on pollutant groups N

C.........  Exit program with normal completion
        CALL M3EXIT( PROGNAME, 0, 0, ' ', 0 )

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats.............94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

94050   FORMAT( A, 1X, I2.2, A, 1X, A, 1X, I6.6, 1X,
     &          A, 1X, I3.3, 1X, A, 1X, I3.3, 1X, A   )

        END PROGRAM TEMPORAL

