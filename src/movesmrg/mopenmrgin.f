
        SUBROUTINE OPENMRGIN

C***********************************************************************
C  subroutine OPENMRGIN body starts at line
C
C  DESCRIPTION:
C      The purpose of this subroutine is to open all of the necessary
C      files for the merge routine and set the episode information 
C      for the calling program.
C
C  PRECONDITIONS REQUIRED:  
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C       Created 2/99 by M. Houyoux
C       09/2025 by HT UNC-IE:  Use M3UTILIO
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
C****************************************************************************
        USE M3UTILIO

C.........  MODULES for public variables
C.........  This module contains the major data structure and control flags
        USE MODMERGE, ONLY: 
     &          MENAME, MSDEV, CFDEV,
     &          NMSRC, MPRJFLAG, MFLAG_BD, MTNAME, MSDATE,
     &          MNIPPA, MEANAM, 
     &          MGNAME, MNGMAT,
     &          PDEV, CDEV, TZONE, SDATE, 
     &          STIME, TSTEP, NSTEPS, EDATE, ETIME, BYEAR, PYEAR,
     &          VARFLAG, SRCGRPFLAG, SGDEV, SUBSECFLAG

C.........  This module contains data structures and flags specific to Movesmrg
        USE MODMVSMRG, ONLY: RPDFLAG, RPVFLAG, RPPFLAG, RPHFLAG, ONIFLAG, RPSFLAG,
     &          TVARNAME, METNAME, XDEV, MDEV, FDEV, CFFLAG, SPDISTFLAG,
     &          SPDPROFLAG, MSNAME_L, MSNAME_S, MNSMATV_L, MNSMATV_S,
     &          MSVDESC_L, MSVDESC_S, MSVUNIT_L, MSVUNIT_S, ETABLEFLAG

C...........  This module contains the information about the source category
        USE MODINFO, ONLY: NMAP, MAPNAM, MAPFIL, NIACT, NSRC, CATEGORY,
     &          ACTVTY

C.........  This module contains the inventory arrays
        USE MODSOURC, ONLY: SPEED, VPOP, CIFIP

C...........   This module contains emission factor information
        USE MODEMFAC, ONLY: MXETYPE, EMTNAM

C.........  This module contains the global variables for the 3-d grid
        USE MODGRID, ONLY: GRDNM, NCOLS, NROWS, VGTYP, VGTOP, VGLVS

C.........  This module is required for the FileSetAPI
        USE MODFILESET

        IMPLICIT NONE

C.........  INCLUDES:
        
        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
        INCLUDE 'SETDECL.EXT'   !  FileSetAPI variables and functions
c       INCLUDE 'IODECL3.EXT'   !  I/O API function declarations

C.........  EXTERNAL FUNCTIONS and their descriptions:
        
c       CHARACTER(2)    CRLF
c       LOGICAL         DSCM3GRD
c       CHARACTER(50)   GETCFDSC  
c       INTEGER         GETIFDSC  
c       INTEGER         INDEX1
c       INTEGER         PROMPTFFILE  
c       CHARACTER(16)   PROMPTMFILE  
c       INTEGER         SECSDIFF  
c       LOGICAL         SETENVVAR

c       EXTERNAL  CRLF, INDEX1, GETCFDSC, GETIFDSC, PROMPTFFILE, 
c    &            PROMPTMFILE, SECSDIFF, SETENVVAR
        LOGICAL      , EXTERNAL :: DSCM3GRD
        CHARACTER(50), EXTERNAL :: GETCFDSC
        INTEGER      , EXTERNAL :: GETIFDSC

C.........   LOCAL VARIABLES and their descriptions:

C.........  Array that contains the names of the inventory variables to read
        CHARACTER(IOVLEN3) IVARNAMS( MXINVARR )

C.........  Other local variables

        INTEGER         I, J, K, M, N, V     ! counters and indices

        INTEGER         IDEV          ! tmp unit number if ENAME is map file
        INTEGER         TDEV          ! unit number for MEPROC file
        INTEGER         SPDEV         ! unit number for SPDPRO file
        INTEGER         SDDEV         ! unit number for SPDIST file
        INTEGER         IOS           ! tmp I/O status
        INTEGER         ISECS         ! tmp duration in seconds
        INTEGER         NPACT         ! no. variables per activity
        INTEGER         NPPOL         ! no. variables per pollutant
        INTEGER         NDIM          ! tmp dimensioning variable 
        INTEGER         NVAR          ! tmp no. variables 
        INTEGER         NINVARR       ! number inventory variables to read

        LOGICAL      :: EFLAG = .FALSE.  ! true: error in routine
        LOGICAL      :: IFLAG = .FALSE.  ! true: episode settings have been init
        LOGICAL      :: OFLAG = .FALSE.  ! true: met info has been init
        LOGICAL      :: YFLAG = .FALSE.  ! true: year/projection info been init
        LOGICAL      :: ZFLAG = .FALSE.  ! true: time zone has been init

        CHARACTER(16)   DUMNAME      ! tmp file name
        CHARACTER(16)   INAME        ! tmp name for inven file of unknown fmt
        CHARACTER(50)   METSCENR     ! met scenario name
        CHARACTER(50)   METCLOUD     ! met cloud scheme name
        CHARACTER(50)   METTMP       ! temporary buffer for met info
        CHARACTER(80)   GDESC        ! grid description
        CHARACTER(256)  MESG         ! message buffer
        CHARACTER(IOVLEN3) COORD3D   ! coordinate system name 
        CHARACTER(IOVLEN3) COORUN3D  ! coordinate system projection units
        CHARACTER(IOVLEN3) PROJTYPE  ! projection type
        CHARACTER(IOVLEN3) OUTGRDNM  ! output grid name

        CHARACTER(16) :: PROGNAME = 'OPENMRGIN' ! program name

C***********************************************************************
C   begin body of subroutine OPENMRGIN

C.........  Initialize gridded information with grid description file
        IF( .NOT. DSCM3GRD( GDNAM3D, GDESC, COORD3D, GDTYP3D, COORUN3D,
     &                      P_ALP3D, P_BET3D, P_GAM3D, XCENT3D, YCENT3D,
     &                      XORIG3D, YORIG3D, XCELL3D, YCELL3D,
     &                      NCOLS3D, NROWS3D, NTHIK3D ) ) THEN

            MESG = 'Could not get Models-3 grid description.'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF
        OUTGRDNM = GDNAM3D

C.........  Check or initialize the grid; do not allow subgrids 
C           when using a variable grid
        IF( VARFLAG ) THEN
            CALL CHKGRID( 'general', 'GRIDDESC', 0, EFLAG )
        ELSE
            CALL CHKGRID( 'general', 'GRIDDESC', 1, EFLAG )
        END IF

C.........  Get inventory file names given source category
        CALL GETINAME( 'MOBILE', MENAME, DUMNAME )

C.........  Prompt for and open inventory file 
        MESG= 'Enter logical name for the MAP ' //
     &        'MOBILE INVENTORY file'
        INAME = MENAME
        IDEV = PROMPTFFILE( MESG, .TRUE., .TRUE., INAME, PROGNAME )

C.........  Read map-formatted inventory file
        CALL RDINVMAP( INAME, IDEV, MENAME, DUMNAME, MSDEV )

C.........  Store source-category-specific header information, 
C           including the inventory pollutants in the file (if any).  Note that 
C           the I/O API header info is passed by include file and the
C           results are stored in module MODINFO.
        CALL GETSINFO( MENAME )

C.........  Ensure that there is at least one activity in the inventory 
C           file, or else this program does not need to be run
        IF( NIACT == 0 ) THEN
            MESG = 'No activities are found in the ' //
     &             'inventory file!  Program cannot be used.'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

C.........  Set inventory variables to read
        IVARNAMS( 1 ) = 'CIFIP'
        IVARNAMS( 2 ) = 'CSCC'
        IVARNAMS( 3 ) = 'TZONES'
        NINVARR = 3

C.........  Allocate memory for and read required inventory characteristics
        CALL RDINVCHR( CATEGORY, MENAME, MSDEV, NSRC, NINVARR, IVARNAMS )

C.........  Read speed and vehicle population data from the inventory
        IF( RPDFLAG ) THEN
            ALLOCATE( SPEED( NSRC ), STAT=IOS )
            CALL CHECKMEM( IOS, 'SPEED', PROGNAME )
            SPEED = 0.0

            IF( .NOT. (SPDPROFLAG .OR. SPDISTFLAG) ) THEN
                M = INDEX1( 'SPEED', NMAP, MAPNAM )
                IF( M <= 0 ) THEN
                    MESG = 'Mobile inventory does not include speed data'
                    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                END IF

                CALL RDMAPPOL( NSRC, 1, 1, 'SPEED', SPEED )
            END IF

C.............  Make sure inventory has VMT as activity (won't be using this
C               data but it needs to be there to make emission processes work)
            M = INDEX1( 'VMT', NMAP, MAPNAM )
            IF( M <= 0 ) THEN
                MESG = 'Mobile inventory does not include VMT data'
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF
        END IF

C.........  Read hotelling data from the inventory
        IF( RPHFLAG ) THEN
            M = INDEX1( 'HOTELLING', NMAP, MAPNAM )
            IF( M <= 0 ) THEN
                MESG = 'Mobile inventory does not include hotelling data'
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF
        END IF

        IF( RPSFLAG ) THEN
            M = INDEX1( 'STARTS', NMAP, MAPNAM )
            IF( M <= 0 ) THEN
                MESG = 'Mobile inventory does not include engine start data'
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF
        END IF

        IF( ONIFLAG ) THEN
            M = INDEX1( 'IDLING', NMAP, MAPNAM )
            IF( M <= 0 ) THEN
                MESG = 'Mobile inventory does not include idling data'
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF
        END IF
 
        IF( RPVFLAG .OR. RPPFLAG ) THEN
            M = INDEX1( 'VPOP', NMAP, MAPNAM )
            IF( M <= 0 ) THEN
                MESG = 'Mobile inventory does not include vehicle ' //
     &                 'population data'
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF
            
            ALLOCATE( VPOP( NSRC ), STAT=IOS )
            CALL CHECKMEM( IOS, 'VPOP', PROGNAME )
            CALL RDMAPPOL( NSRC, 1, 1, 'VPOP', VPOP )
        END IF

C.........  Build unique lists of SCCs and country/state/county codes
C           from the inventory arrays
        CALL GENUSLST

C.........  Get number of sources from MODINFO and store in MODMERGE variable
        NMSRC = NSRC

C.........  Determine the year and projection status of the inventory
        CALL CHECK_INVYEAR( MENAME, MPRJFLAG, FDESC3D )

        IF( RPDFLAG .OR. RPHFLAG .OR. RPSFLAG .OR. ONIFLAG  ) THEN

C.............  Open all temporal files for either by-day or standard
C               processing. 
C.............  Compare headers to make sure files are consistent.
            CALL OPEN_TMP_FILES( 'MOBILE', MFLAG_BD, MTNAME, MSDATE)

C.............  Determine the year and projection status of the hourly
            CALL CHECK_INVYEAR( MTNAME( 1 ), MPRJFLAG, FDESC3D )

        END IF

C.........  Open gridding matrix, compare number of sources, and
C           compare grid information
        MGNAME = PROMPTMFILE(
     &        'Enter logical name for the MOBILE GRIDDING MATRIX',
     &        FSREAD3, 'MGMAT', PROGNAME )
     
        IF( .NOT. DESC3( MGNAME ) ) THEN
            MESG = 'Could not get description of file "' //
     &             TRIM( MGNAME ) // '" '
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF
        
        IF( VARFLAG ) THEN
            DUMNAME = GETCFDSC( FDESC3D, '/VARIABLE GRID/', .TRUE. )
        END IF
        
        CALL CHKSRCNO( 'mobile', 'MGMAT', NTHIK3D, NMSRC, EFLAG )

C.........  Check the grid definition; do not allow subgrids if using
C           a variable grid
        IF( VARFLAG ) THEN
            CALL CHKGRID( 'mobile', 'GMAT', 0, EFLAG )
        ELSE
            CALL CHKGRID( 'mobile', 'GMAT', 1, EFLAG )
        END IF
        
        MNGMAT = NCOLS3D

C.........  Open mole-based speciation matrix, compare number of sources, and store
C           speciation variable descriptions.
        MSNAME_L = PROMPTSET( 
     &           'Enter logical name for the MOLE-BASED SPECIATION MATRIX',
     &           FSREAD3, 'MSMAT_L', PROGNAME )

        IF ( .NOT. DESCSET( MSNAME_L, ALLFILES ) ) THEN
            MESG = 'Could not get description of file set "' //
     &             MSNAME_L( 1:LEN_TRIM( MSNAME_L ) ) // '"'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        ENDIF

        CALL CHKSRCNO( 'mobile', 'MSMAT_L', NROWS3D, NMSRC, EFLAG)
        MNSMATV_L = NVARSET
        ALLOCATE( MSVDESC_L( MNSMATV_L ), STAT=IOS )
        CALL CHECKMEM( IOS, 'MSVDESC_L', PROGNAME )
        ALLOCATE( MSVUNIT_L( MNSMATV_L ), STAT=IOS )
        CALL CHECKMEM( IOS, 'MSVUNIT_L', PROGNAME )
        CALL STORE_VDESCS( 1, 1, MNSMATV_L, .TRUE., MSVDESC_L )
        CALL STORE_VUNITS( 1, 1, MNSMATV_L, .TRUE., MSVUNIT_L )

C.........  Open mass-based speciation matrix, compare number of sources, and store
C           speciation variable descriptions.
        MSNAME_S = PROMPTSET( 
     &           'Enter logical name for the MASS-BASED SPECIATION MATRIX',
     &           FSREAD3, 'MSMAT_S', PROGNAME )

        IF ( .NOT. DESCSET( MSNAME_S, ALLFILES ) ) THEN
            MESG = 'Could not get description of file set "' //
     &             MSNAME_S( 1:LEN_TRIM( MSNAME_S ) ) // '"'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        ENDIF

        CALL CHKSRCNO( 'mobile', 'MSMAT_S', NROWS3D, NMSRC, EFLAG)
        MNSMATV_S = NVARSET
        ALLOCATE( MSVDESC_S( MNSMATV_S ), STAT=IOS )
        CALL CHECKMEM( IOS, 'MSVDESC_S', PROGNAME )
        ALLOCATE( MSVUNIT_S( MNSMATV_S ), STAT=IOS )
        CALL CHECKMEM( IOS, 'MSVUNIT_S', PROGNAME )
        CALL STORE_VDESCS( 1, 1, MNSMATV_S, .TRUE., MSVDESC_S )
        CALL STORE_VUNITS( 1, 1, MNSMATV_S, .TRUE., MSVUNIT_S )

C.........  Check that variables in mole and mass speciation matrices match
        IF( MNSMATV_L .NE. MNSMATV_S ) THEN
            WRITE( MESG,94010 )
     &         'ERROR: Mole-based speciation matrix contains ', 
     &         MNSMATV_L, 'variables but mass-based speciation ' // 
     &         'matrix contains ', MNSMATV_S, 'variables.'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF
        
        DO I = 1, MNSMATV_L
            IF( MSVDESC_L( I ) .NE. MSVDESC_S( I ) ) THEN
                MESG = 'ERROR: Variable descriptions are not ' //
     &            'consistent between mole- and mass-based ' //
     &            'speciation matrices.'
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF
        END DO

C.........  Open meteorology file
        IF( RPDFLAG .OR. RPVFLAG .OR. RPHFLAG .OR. RPSFLAG .OR. ONIFLAG ) THEN
            IF( .NOT. ETABLEFLAG ) THEN
                METNAME = PROMPTMFILE(
     &              'Enter logical name for the METCRO2D meteorology file', 
     &              FSREAD3, 'MET_CRO_2D', PROGNAME )

                IF( .NOT. DESC3( METNAME ) ) THEN
                    MESG = 'Could not get description of file "' //
     &                  METNAME( 1:LEN_TRIM( METNAME ) ) // '" '
                    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                END IF

                CALL CHKGRID( 'mobile', 'GRID', 0, EFLAG )

C.................  Check modeling period
                CALL UPDATE_TIME_INFO( METNAME, .FALSE. )

C.................  Make sure met file contains requested temperature variable
                J = INDEX1( TVARNAME, NVARS3D, VNAME3D )
                IF( J <= 0 ) THEN
                    MESG = 'ERROR: Could not find "' // TRIM( TVARNAME ) //
     &                   '" in file "' // TRIM( METNAME )
                    CALL M3MESG( MESG )
                END IF

            END IF

        ELSE
            METNAME = PROMPTMFILE(
     &           'Enter logical name for the METMOVES meteorology file', 
     &           FSREAD3, 'METMOVES', PROGNAME )

            IF( .NOT. DESC3( METNAME ) ) THEN
                MESG = 'Could not get description of file "' //
     &                  METNAME( 1:LEN_TRIM( METNAME ) ) // '" '
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF

            CALL CHKGRID( 'mobile', 'GRID', 0, EFLAG )
        END IF
 
C.........  Get file name for inventory pollutants codes/names
        MESG = 'Enter logical name for INVENTORY DATA TABLE file'
        PDEV = PROMPTFFILE( MESG, .TRUE., .TRUE., 'INVTABLE',
     &                      PROGNAME )

C.........  Get country, state, and county names no matter what, because it is
C           needed to allocate memory for the state and county totals, even
C           when they aren't going to be output
        CDEV = PROMPTFFILE( 
     &             'Enter logical name for COUNTRY, STATE, AND ' //
     &             'COUNTY file', .TRUE., .TRUE., 'COSTCY', PROGNAME )

C.........  Open source groups file if needed
        IF( SRCGRPFLAG ) THEN
            MESG = 'Enter logical name for SOURCE GROUPS file'
            SGDEV = PROMPTFFILE( MESG, .TRUE., .TRUE., 
     &                           'SOURCE_GROUPS', PROGNAME )
        END IF

C.........  Open sub-sector source groups file if needed
        IF( SUBSECFLAG ) THEN
            MESG = 'Enter logical name for SUB-SECTOR SOURCE GROUPS file'
            SGDEV = PROMPTFFILE( MESG, .TRUE., .TRUE.,
     &                           'SUB_SEC_SOURCES', PROGNAME )
        END IF

C.........  Get emission processes file name
        TDEV = PROMPTFFILE( 
     &           'Enter logical name for EMISSION PROCESSES file',
     &           .TRUE., .TRUE., 'MEPROC', PROGNAME )

        CALL RDEPROC( TDEV )

C.........  Store process/pollutants combinations for correct activity
        IF( RPDFLAG ) THEN
            M = INDEX1( 'VMT', NIACT, ACTVTY )
        END IF

        IF( RPHFLAG ) THEN
            M = INDEX1( 'HOTELLING', NIACT, ACTVTY )
        END IF

        IF( RPSFLAG ) THEN
            M = INDEX1( 'STARTS', NIACT, ACTVTY )
        END IF

        IF( ONIFLAG ) THEN
            M = INDEX1( 'IDLING', NIACT, ACTVTY )
        END IF
 
        IF( RPVFLAG .OR. RPPFLAG ) THEN
            M = INDEX1( 'VPOP', NIACT, ACTVTY )
        END IF
        
        IF( M <= 0 ) THEN
            MESG = 'INTERNAL ERROR: Could not find expected activity'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

        MNIPPA = 0
        DO I = 1, MXETYPE
            IF( EMTNAM( I,M ) .NE. ' ' ) THEN
                MNIPPA = MNIPPA + 1
            END IF
        END DO

        ALLOCATE( MEANAM( MNIPPA ), STAT=IOS )
        CALL CHECKMEM( IOS, 'MEANAM', PROGNAME )

        DO I = 1, MNIPPA
            MEANAM( I ) = EMTNAM( I,M )
        END DO

C.........  Get county cross-reference file
        XDEV = PROMPTFFILE( 
     &           'Enter logical name for MCXREF cross-reference file',
     &           .TRUE., .TRUE., 'MCXREF', PROGNAME )

C.........  Get county fuel month file
        MDEV = PROMPTFFILE(
     &           'Enter logical name for fuel month reference file',
     &           .TRUE., .TRUE., 'MFMREF', PROGNAME )

C.........  Get reference county emission factors file list
        FDEV = PROMPTFFILE(
     &           'Enter logical name for reference county file list',
     &           .TRUE., .TRUE., 'MRCLIST', PROGNAME )

C.........  Open and read avereage speed distribution data
        IF( RPDFLAG .AND. SPDISTFLAG ) THEN
            SDDEV = PROMPTFFILE(
     &              'Enter logical name for average speed distribution file',
     &              .TRUE., .TRUE., 'SPDIST', PROGNAME )
            CALL RDSPDIST( SDDEV )
            SPDPROFLAG = .FALSE.
        END IF

C.........  Open and read hourly speed data
        IF( RPDFLAG .AND. SPDPROFLAG ) THEN
            SPDEV = PROMPTFFILE(
     &              'Enter logical name for speed profiles file',
     &              .TRUE., .TRUE., 'SPDPRO', PROGNAME )
            CALL RDSPDPRO( SPDEV )
        END IF

C.........  Get control factor file 
        IF( CFFLAG ) THEN
            CFDEV = PROMPTFFILE(
     &              'Enter logical name for control factor file',
     &              .TRUE., .TRUE., 'CFPRO', PROGNAME )
        END IF

C.........  If there were any errors inputing files or while comparing
C           with one another, then abort
        IF( EFLAG ) THEN
           MESG = 'Problems opening input files. See ERROR(S) above.'
           CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

C.........  If we are using temporalized emissions, then update date/time and
C           duration using environment variable settings, then prompt.
        CALL GETM3EPI( TZONE, SDATE, STIME, TSTEP, NSTEPS )
        TSTEP = 10000   ! only 1-hour time steps supported
        EDATE = SDATE
        ETIME = STIME
        CALL NEXTIME( EDATE, ETIME, ( NSTEPS-1 ) * TSTEP )

C.........  Compare base year with episode and warn if not consistent
        IF( BYEAR .NE. 0 .AND. SDATE / 1000 .NE. BYEAR ) THEN

            WRITE( MESG,94010 ) 'WARNING: Inventory base year ', BYEAR, 
     &             'is inconsistent with year ' // CRLF() // BLANK10 //
     &             'of episode start date', SDATE/1000
            CALL M3MSG2( MESG )

        ENDIF

C.........  Give a note if running for a projected year
        IF( PYEAR .NE. BYEAR ) THEN

            WRITE( MESG,94010 ) 'NOTE: Emissions based on projected '//
     &             'year', PYEAR
            CALL M3MSG2( MESG )

        END IF

C.........  Reset output grid name in case meteorology files use different name
        GRDNM = OUTGRDNM

C.........  Write message stating grid name and description
        N = LEN_TRIM( GRDNM )
        MESG = 'NOTE: Output grid "' // GRDNM( 1:N ) // 
     &         '" set; described as' // CRLF() // BLANK10 // GDESC
        CALL M3MSG2( MESG )

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats.............94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

C******************  INTERNAL SUBPROGRAMS  *****************************
 
        CONTAINS

C.............  This subprogram updates the time (episode) information
C               and compares to the existing information, if it has been
C               previously set.
            SUBROUTINE UPDATE_TIME_INFO( FILNAM, CHKTZONE )

C.............  Subprogram arguments
            CHARACTER(*) FILNAM
            LOGICAL      CHKTZONE

C.............  Local variables
            INTEGER ISECS   ! number of seconds different between dates/times
            INTEGER ED      ! tmp ending date
            INTEGER ET      ! tmp ending time
            INTEGER LOCZONE ! tmp time zone

C----------------------------------------------------------------------

C.............  If time information has already been initialized...
            IF( IFLAG ) THEN
                ISECS = SECSDIFF( SDATE, STIME, SDATE3D, STIME3D )

                IF( ISECS .GT. 0 ) THEN  ! SDATE3D/STIME3D are later
                    SDATE = SDATE3D
                    STIME = STIME3D
                END IF

                ED = SDATE3D
                ET = STIME3D
                CALL NEXTIME( ED, ET, ( MXREC3D-1 ) * TSTEP3D )
        
                ISECS = SECSDIFF( EDATE, ETIME, ED, ET )

                IF( ISECS .LT. 0 ) THEN  ! ED/ET are earlier
                    EDATE = ED
                    ETIME = ET
                END IF

                NSTEPS = 1+ SECSDIFF( SDATE, STIME, EDATE, ETIME )/ 3600

                IF( NSTEPS .LE. 0 ) THEN
                    MESG = 'Because of file ' // FILNAM // 
     &                     ', dates and times do not overlap at all!'
                    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

                END IF

C.............  If time information needs to be initialized...
            ELSE
                SDATE  = SDATE3D
                STIME  = STIME3D
                NSTEPS = MXREC3D

                EDATE  = SDATE
                ETIME  = STIME
                CALL NEXTIME( EDATE, ETIME, ( NSTEPS-1 ) * TSTEP3D )

                IFLAG = .TRUE.

            END IF

C.............  Make sure that time step is one hour
            IF( TSTEP3D .NE. 10000 ) THEN

                EFLAG = .TRUE.
                MESG = 'ERROR: Time step is not one hour in ' // 
     &                 FILNAM // ' file!'
                CALL M3MSG2( MESG )

            END IF

C.............  Retrieve and compare time zone
            IF( CHKTZONE ) THEN

                LOCZONE = GETIFDSC( FDESC3D, '/TZONE/', .TRUE. )

                IF( ZFLAG .AND. LOCZONE .NE. TZONE ) THEN
                    EFLAG = .TRUE.
                    WRITE( MESG,94010 )
     &                 'Time zone ', LOCZONE, 'in ' // FILNAM // 
     &                 ' hourly emissions file is not consistent ' //
     &                 'with initialized value of', TZONE
                    CALL M3MSG2( MESG )

                ELSE IF( .NOT. ZFLAG ) THEN
                    ZFLAG = .TRUE.
                    TZONE = LOCZONE

                    MESG = 'NOTE: Time zone initialized using ' // 
     &                     FILNAM // ' hourly emissions file.'

                    CALL M3MSG2( MESG )
                END IF
            END IF

C------------------  FORMAT  STATEMENTS   -----------------------------

C...........   Internal buffering formats.............94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

            END SUBROUTINE UPDATE_TIME_INFO

C----------------------------------------------------------------------
C----------------------------------------------------------------------
C.............  This subprogram opens the temporal emissions files. If their 
C               are multiple files, it compares the files to make sure that they
C               are consistent with each other.  The number of sources
C               are compared to the master number of sources.
            SUBROUTINE OPEN_TMP_FILES( LOCCAT, LBDSTAT, FNAME, SDATE )

C.............  Subprogram arguments
            CHARACTER(*), INTENT (IN) :: LOCCAT
            LOGICAL     , INTENT (IN) :: LBDSTAT
            CHARACTER(*), INTENT(OUT) :: FNAME( 7 )
            INTEGER     , INTENT(OUT) :: SDATE( 7 )

C.............  Local parameters
            CHARACTER(3), PARAMETER :: SUFFIX( 7 ) = 
     &                                ( / 'MON', 'TUE', 'WED', 'THU', 
     &                                    'FRI', 'SAT', 'SUN'        / )

C.............  Local allocatable arrays
            CHARACTER(IOVLEN3), ALLOCATABLE :: LOCVNAM ( : )
            CHARACTER(IOULEN3), ALLOCATABLE :: LOCVUNIT( : )

C.............  Local arrays
            INTEGER        IDX( 7 )     ! index for per-file arrays

C.............  Local variables
            INTEGER        D, L, N      ! counters and indices

            INTEGER        LOCZONE   ! tmp time zone
            INTEGER        LOCNVAR   ! tmp local number of variables in file 
            INTEGER        NFILE     ! no. hourly emission files

            LOGICAL     :: NFLAG = .FALSE.  ! true: no. vars inconsistent
            LOGICAL     :: VFLAG = .FALSE.  ! true: var names inconsistent
            LOGICAL     :: UFLAG = .FALSE.  ! true: var units inconsistent

            CHARACTER      CRL      ! 1-letter src category indicator
            CHARACTER(16)  TMPNAM   ! temporary logical file name
            CHARACTER(300) MESG     ! message buffer

C----------------------------------------------------------------------

            IF( LOCCAT .EQ. 'MOBILE' ) CRL = 'M'

C.............  Set the number of files and open the files...
C.............  For by-day processing...
            IF( LBDSTAT ) THEN
                NFILE = 7

                DO D = 1, NFILE

                    MESG = 'Enter logical name for the ' // SUFFIX( D )
     &                     // ' ' // LOCCAT // ' HOURLY EMISSIONS file'
                    TMPNAM = CRL // 'TMP_' // SUFFIX( D )

                    FNAME( D ) = PROMPTSET( MESG,FSREAD3,
     &                                        TMPNAM,PROGNAME )
                    IDX( D ) = D
                END DO

C.............  For standard processing...
            ELSE
                NFILE = 1

                MESG = 'Enter logical name for the ' // LOCCAT // 
     &                 ' HOURLY EMISSIONS file'
                TMPNAM = CRL // 'TMP'

                FNAME = PROMPTSET( MESG,FSREAD3,TMPNAM,PROGNAME ) ! array
                IDX( NFILE ) = 1

            END IF

C.............  Loop through each file and ensure they are consistent
            DO D = 1, NFILE

                TMPNAM = FNAME( IDX( D ) )

C.................  Get header and compare source number and time range
                IF ( .NOT. DESCSET( TMPNAM, ALLFILES ) ) THEN
                    MESG = 'Could not get description of file set "' //
     &                     TMPNAM( 1:LEN_TRIM( TMPNAM ) ) // '"'
                    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                ENDIF

C.................  Store the starting date
                SDATE( IDX( D ) ) = SDATE3D

C.................  Check the number of sources
                CALL CHKSRCNO( 'mobile', TMPNAM, NROWS3D, 
     &                         NMSRC, EFLAG )

C.................  For standard processing, compare time info to master
                IF( .NOT. LBDSTAT .AND. D .EQ. 1 ) THEN
                    CALL UPDATE_TIME_INFO( TMPNAM, .TRUE. )
                END IF

C.................  For by-day files, make sure that the file starts at hour 0
                IF( LBDSTAT .AND. STIME3D .NE. 0 ) THEN
                    EFLAG = .TRUE.
                    L = LEN_TRIM( TMPNAM )
                    WRITE( MESG,94010 ) 'ERROR: Start time of', STIME3D,
     &                     'in file "'// TMPNAM( 1:L ) // 
     &                     '" is invalid.' // CRLF() // BLANK10 //
     &                     'Only start time of 000000 is valid for' //
     &                     'processing by day.'
                    CALL M3MSG2( MESG )

                END IF

C.................  Make sure that the file has at least 24 hours 
                IF( LBDSTAT .AND. MXREC3D .LT. 24 ) THEN
                    EFLAG = .TRUE.
                    L = LEN_TRIM( TMPNAM )
                    WRITE( MESG,94010 ) 'ERROR: Number of hours', 
     &                     MXREC3D, 'in file "'// TMPNAM( 1:L ) // 
     &                     '" is invalid.' // CRLF() // BLANK10 //
     &                     'Minimum number of 24 hours is needed for' //
     &                     'processing by day.'
                    CALL M3MSG2( MESG )

                END IF

                LOCZONE = GETIFDSC( FDESC3D, '/TZONE/', .TRUE. )

                IF( ZFLAG .AND. LOCZONE .NE. TZONE ) THEN
                    EFLAG = .TRUE.
                    WRITE( MESG,94010 )
     &                 'Time zone ', LOCZONE, 'in ' // TMPNAM // 
     &                 ' hourly emissions file is not consistent ' //
     &                 'with initialized value of', TZONE
                    CALL M3MSG2( MESG )

                ELSE IF( .NOT. ZFLAG ) THEN
                    ZFLAG = .TRUE.
                    TZONE = LOCZONE

                    MESG = 'NOTE: Time zone initialized using ' // 
     &                     TMPNAM // ' hourly emissions file.'

                    CALL M3MSG2( MESG )
                END IF

C.................  For first file, store the pollutant names and units for
C                   making comparisons with other files.
                IF( D .EQ. 1 ) THEN

                    LOCNVAR = NVARS3D
                    ALLOCATE( LOCVNAM( LOCNVAR ), STAT=IOS )
                    CALL CHECKMEM( IOS, 'LOCVNAM', PROGNAME )
                    ALLOCATE( LOCVUNIT( LOCNVAR ), STAT=IOS )
                    CALL CHECKMEM( IOS, 'LOCVUNIT', PROGNAME )

                    LOCVNAM ( 1:LOCNVAR ) = VNAMESET( 1:LOCNVAR )
                    LOCVUNIT( 1:LOCNVAR ) = VUNITSET( 1:LOCNVAR )

C.................  Compare the pollutant names and units
                ELSE

C.....................  Check to make sure the number is consistent first
                    IF( NVARSET .NE. LOCNVAR ) NFLAG = .TRUE.

C.....................  Make sure no overflows                    
                    N = MIN( NVARSET, LOCNVAR )

C.....................  compare variable names and units among files
                    DO V = 1, N
                        IF( LOCVNAM( V ) .NE. VNAMESET( V ) ) THEN
                            VFLAG = .TRUE.
                        END IF

                        IF( LOCVUNIT( V ) .NE. VUNITSET( V ) ) THEN
                            UFLAG = .TRUE.
                        END IF
                    END DO

                END IF

            END DO

C.............  Write message and set error if any inconsistencies
            IF( NFLAG ) THEN
c bbh               EFLAG = .TRUE.  ! removed to prevent false errer of odd nubmer
c                                     of species for tmp files
                MESG = 'WARNING: ' // LOCCAT // ' source hourly ' //
     &                 'emission files have inconsistent ' //
     &                 CRLF() // BLANK10 // 'number of variables.'
                CALL M3MSG2( MESG )
            END IF

            IF( VFLAG ) THEN
                EFLAG = .TRUE.
                MESG = 'ERROR: ' // LOCCAT // ' source hourly ' //
     &                 'emission files have inconsistent ' //
     &                 CRLF() // BLANK10 // 'variable names.'
                CALL M3MSG2( MESG )
            END IF

            IF( UFLAG ) THEN
                EFLAG = .TRUE.
                MESG = 'ERROR: ' // LOCCAT // ' source hourly ' //
     &                 'emission files have inconsistent ' //
     &                 CRLF() // BLANK10 // 'variable units.'
                CALL M3MSG2( MESG )
            END IF

C.............  Deallocate local memory
            DEALLOCATE( LOCVNAM, LOCVUNIT )

            RETURN

C------------------  FORMAT  STATEMENTS   -----------------------------

C...........   Internal buffering formats.............94xxx

94010       FORMAT( 10( A, :, I8, :, 1X ) )

            END SUBROUTINE OPEN_TMP_FILES

C----------------------------------------------------------------------
C----------------------------------------------------------------------
C.............  This subprogram initializes and checks the inventory year
C               of the emissions and the projection status
            SUBROUTINE CHECK_INVYEAR( FNAME, PRJFLAG, IODESC )

C.............  Subprogram arguments
            CHARACTER(*), INTENT (IN)     :: FNAME
            LOGICAL     , INTENT (IN OUT) :: PRJFLAG
            CHARACTER(*), INTENT (IN)     :: IODESC( * )

C.............  Local variables
            INTEGER           L
            INTEGER           YY      ! tmp year
            LOGICAL           STRICT  ! flag for strict checks or not
            CHARACTER(20)     BUFFER  ! program name buffer
            INTEGER,  SAVE :: FLEN    ! name length of savnam
            CHARACTER(IOVLEN3), SAVE :: SAVNAM  ! name of file used to init

C----------------------------------------------------------------------

            STRICT = .TRUE.

C.............  First determine whether to abort when projected year does not
C               match.  This is used for reactivity matrices, which will
C               always have a projection year, even if the inventory isn't
C               projected.
            IF( .NOT. PRJFLAG ) THEN
                BUFFER = GETCFDSC( FDESC3D, '/FROM/', .FALSE. )
                IF( BUFFER .EQ. 'OPENRMAT' ) STRICT = .FALSE.
            END IF

C.............  If time information has already been initialized...
            IF( YFLAG ) THEN

                YY = GETIFDSC( IODESC, '/PROJECTED YEAR/', .FALSE. )
                IF( YY .LE. 0 ) THEN

                    YY = GETIFDSC( IODESC, '/BASE YEAR/', .FALSE. ) 
                    IF( YY .NE. BYEAR ) THEN
                        WRITE( MESG,94010 ) 
     &                        'Base year of ' // FNAME // ' file:', YY,
     &                        CRLF() // BLANK10 //
     &                        ', does not equal emissions year of ' //
     &                        SAVNAM( 1:FLEN ) // ' file:', BYEAR

C.........................  If there is projection, abort
                        IF ( PRJFLAG ) THEN
                            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

C.........................  Otherwise, make it a warning
                        ELSE
                            L = LEN_TRIM( MESG )
                            MESG = 'WARNING: ' // MESG( 1:L )
                            CALL M3MSG2( MESG )
                        END IF

                    END IF

                ELSE IF ( STRICT            .AND. 
     &                    YY     .GT. 0     .AND. 
     &                    YY     .NE. PYEAR      ) THEN

                    WRITE( MESG,94010 ) 
     &                    'Projected year of ' // FNAME // ' file:', YY,
     &                    CRLF() // BLANK10 //
     &                    ', does not equal emissions year of ' //
     &                    SAVNAM( 1:FLEN ) // ' file:', PYEAR
                    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

                END IF

C.............  If year information needs to be initialized...
            ELSE
                
                BYEAR = GETIFDSC( IODESC, '/BASE YEAR/', .FALSE. ) 
                PYEAR = GETIFDSC( IODESC, '/PROJECTED YEAR/', .FALSE. )

                IF( PYEAR .GT. 0 ) THEN
                    PRJFLAG = .TRUE.
                ELSE
                    PYEAR = BYEAR
                END IF

                SAVNAM = FNAME
                FLEN   = LEN_TRIM( SAVNAM )
                YFLAG  = .TRUE.

            END IF

C------------------  FORMAT  STATEMENTS   -----------------------------

C...........   Internal buffering formats.............94xxx

94010       FORMAT( 10( A, :, I8, :, 1X ) )

            END SUBROUTINE CHECK_INVYEAR

C----------------------------------------------------------------------
C----------------------------------------------------------------------
C.............  This subprogram updates the met information and compares to 
C               the existing information, if it has been previously set.
            SUBROUTINE CHECK_MET_INFO( CATDESC )

C.............  Subprogram arguments
            CHARACTER(*) CATDESC  ! category descriptions

C.............  Local variables
            INTEGER       L, L1, L2  ! length of strings
            CHARACTER(30) FILDESC    ! description of input file

C----------------------------------------------------------------------

C.............  Set tmp rows, columns, and total cells depending on file type
            IF( CATDESC .EQ. 'biogenics' ) THEN
                FILDESC = 'gridded emissions file'

            ELSEIF( CATDESC .EQ. 'mobile' ) THEN
                FILDESC = 'hourly emissions file'

            ELSEIF( CATDESC .EQ. 'point' ) THEN
                FILDESC = 'layer fractions file'

            ELSE
                MESG= 'INTERNAL ERROR: Category description "' // 
     &                CATDESC// '" not known in call to CHECK_MET_INFO!'
                CALL M3MSG2( MESG )
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

            ENDIF

            L = LEN_TRIM( FILDESC )

C.............  If met information has already been initialized, then compare
C               existing to this file.
            IF( OFLAG ) THEN

                METTMP = GETCFDSC( FDESC3D, '/MET SCENARIO/', .TRUE. )
                IF ( METTMP .NE. METSCENR ) THEN

                    L1 = LEN_TRIM( METTMP )
                    L2 = LEN_TRIM( METSCENR )

                    EFLAG = .TRUE.
                    MESG = 'ERROR: Meteorology scenario name "' // 
     &                     METTMP( 1:L1 ) // '" in ' // CATDESC //
     &                     FILDESC( 1:L ) // ' is inconsistent with '//
     &                     'initialized value "'// METSCENR(1:L2)// '"'
                    CALL M3MSG2( MESG )

                END IF

                METTMP = GETCFDSC( FDESC3D, '/CLOUD SCHEME/', .TRUE. )
                IF ( METTMP .NE. METCLOUD ) THEN

                    L1 = LEN_TRIM( METTMP )
                    L2 = LEN_TRIM( METCLOUD )

                    EFLAG = .TRUE.
                    MESG = 'ERROR: Meteorology cloud scheme "' // 
     &                     METTMP( 1:L1 ) // '" in ' // CATDESC //
     &                     FILDESC( 1:L ) // ' is inconsistent with '//
     &                     'initialized value "'// METCLOUD(1:L2)// '"'
                    CALL M3MSG2( MESG )

                END IF

C.............  Initialize meteorology information
            ELSE

                OFLAG    = .TRUE.
                METSCENR = GETCFDSC( FDESC3D, '/MET SCENARIO/', .TRUE. )
                METCLOUD = GETCFDSC( FDESC3D, '/CLOUD SCHEME/', .TRUE. )

                MESG = 'NOTE: Meteorology description initialized '//
     &                 'using '// CATDESC// ' '// FILDESC( 1:L )// '.'
                CALL M3MSG2( MESG )

            ENDIF

            END SUBROUTINE CHECK_MET_INFO

C----------------------------------------------------------------------
C----------------------------------------------------------------------
C.............  This subprogram stores I/O API NetCDF variable descriptions into
C               a local array based on indices in subprogram call.
            SUBROUTINE STORE_VDESCS( ISTART,INCRMT,NDESC,LFSET,DESCS )

            INCLUDE 'SETDECL.EXT'   !  FileSetAPI variables and functions

C.............  Subprogram arguments
            INTEGER     , INTENT (IN) :: ISTART   ! starting position in VDESCS of names
            INTEGER     , INTENT (IN) :: INCRMT   ! increment of VDESCS for names
            INTEGER     , INTENT (IN) :: NDESC    ! number of descriptions
            LOGICAL     , INTENT (IN) :: LFSET    ! number of descriptions
            CHARACTER(*), INTENT(OUT) :: DESCS( NDESC )! stored variable descriptions

C.............  Local variables
            INTEGER  I, J

C----------------------------------------------------------------------

            DESCS = ' '

            J = ISTART
            DO I = 1, NDESC

                IF( LFSET ) THEN  ! From FileSetAPI
                    DESCS( I ) = TRIM( VDESCSET( J ) )
                ELSE              ! From standard I/O API
                    DESCS( I ) = TRIM( VDESC3D( J ) )
                END IF

                J = J + INCRMT

            END DO
 
            END SUBROUTINE STORE_VDESCS

C----------------------------------------------------------------------
C----------------------------------------------------------------------
C.............  This subprogram stores I/O API NetCDF variable units into
C               a local array based on indices in subprogram call.
            SUBROUTINE STORE_VUNITS( ISTART,INCRMT,NUNIT,LFSET,UNITS )

            INCLUDE 'SETDECL.EXT'   !  FileSetAPI variables and functions

C.............  Subprogram arguments
            INTEGER     , INTENT (IN) :: ISTART        ! starting position in VDESCS of names
            INTEGER     , INTENT (IN) :: INCRMT        ! increment of VDESCS for names
            INTEGER     , INTENT (IN) :: NUNIT         ! number of units
            LOGICAL     , INTENT (IN) :: LFSET         ! number of descriptions
            CHARACTER(*), INTENT(OUT) :: UNITS( NUNIT )! stored variable units

C.............  Local variables
            INTEGER  I, J, L

C----------------------------------------------------------------------

            UNITS = ' '

            J = ISTART
            DO I = 1, NUNIT

                IF( LFSET ) THEN  ! From FileSetAPI
                    UNITS( I ) = TRIM( VUNITSET( J ) )
                ELSE              ! From standard I/O API
                    UNITS( I ) = TRIM( UNITS3D( J ) )
                END IF 

                J = J + INCRMT

            END DO
 
            END SUBROUTINE STORE_VUNITS

        END SUBROUTINE OPENMRGIN
