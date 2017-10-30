
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

C.........  MODULES for public variables
C.........  This module contains the major data structure and control flags
        USE MODMERGE, ONLY: 
     &          MENAME, MSDEV, NMSRC, MNIPPA, MEANAM, PDEV, CDEV, 
     &          TZONE, SDATE, STIME, TSTEP, NSTEPS, EDATE, ETIME,
     &          BYEAR, PYEAR, VARFLAG, APFLAG, APRT_ELEV, APRT_CODE,
     &          NAPRT, MPRJFLAG

C.........  This module contains data structures and flags specific to Movesmrg
        USE MODMVSMRG, ONLY: 
     &          METNAME, MGRNAME,
     &          MSNAME_L, MSNAME_S, MNSMATV_L, MNSMATV_S,
     &          MSVDESC_L, MSVDESC_S, MSVUNIT_L, MSVUNIT_S

C...........  This module contains the information about the source category
        USE MODINFO, ONLY: NMAP, MAPNAM, MAPFIL, NSRC, CATEGORY

C.........  This module contains the inventory arrays
        USE MODSOURC, ONLY: CIFIP, CSCC, TZONES, CLINK, CDPTID, CARRID

C.........  This module contains the global variables for the 3-d grid
        USE MODGRID, ONLY: GRDNM, NCOLS, NROWS, VGTYP, VGTOP, VGLVS, NLAYS

C.........  This module is required for the FileSetAPI
        USE MODFILESET

        IMPLICIT NONE

C.........  INCLUDES:
        
        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
        INCLUDE 'SETDECL.EXT'   !  FileSetAPI variables and functions
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations

C.........  EXTERNAL FUNCTIONS and their descriptions:
        
        CHARACTER(2)    CRLF
        LOGICAL         DSCM3GRD
        CHARACTER(50)   GETCFDSC  
        INTEGER         GETIFDSC  
        INTEGER         INDEX1
        INTEGER         PROMPTFFILE, GETFLINE  
        CHARACTER(16)   PROMPTMFILE  
        INTEGER         SECSDIFF  
        LOGICAL         SETENVVAR, BLKORCMT, USEEXPGEO
        REAL            STR2REAL

        EXTERNAL  CRLF, INDEX1, GETCFDSC, GETIFDSC, PROMPTFFILE, 
     &            PROMPTMFILE, SECSDIFF, SETENVVAR, BLKORCMT,
     &            GETFLINE, USEEXPGEO, STR2REAL

C.........   LOCAL VARIABLES and their descriptions:

C.........  Array that contains the names of the inventory variables to read
        CHARACTER(IOVLEN3) IVARNAMS( MXINVARR )

C.........  Other local variables

        INTEGER         I, J, K, M, N, V     ! counters and indices

        INTEGER         IDEV          ! tmp unit number if ENAME is map file
        INTEGER         IOS           ! tmp I/O status
        INTEGER         ISECS         ! tmp duration in seconds
        INTEGER         NPACT         ! no. variables per activity
        INTEGER         NPPOL         ! no. variables per pollutant
        INTEGER         NVAR          ! tmp no. variables 
        INTEGER         NINVARR       ! number inventory variables to read

        INTEGER         APDEV         ! channel no of aiprort height input file (option)

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

C.........  Set inventory variables to read
        IVARNAMS( 1 ) = 'CIFIP'
        IVARNAMS( 2 ) = 'CSCC'
        IVARNAMS( 3 ) = 'CLINK'
        IVARNAMS( 4 ) = 'TZONES'
        IVARNAMS( 5 ) = 'CDPTID'
        IVARNAMS( 6 ) = 'CARRID'
        IVARNAMS( 7 ) = 'CSOURC'
        NINVARR = 7

C.........  Allocate memory for and read required inventory characteristics
        CALL RDINVCHR( CATEGORY, MENAME, MSDEV, NSRC, NINVARR, IVARNAMS )

C.........  Build unique lists of SCCs and country/state/county codes
C           from the inventory arrays
        CALL GENUSLST

C.........  Get number of sources from MODINFO and store in MODMERGE variable
        NMSRC = NSRC

C.........  stores HAP factors for Turbine engines
        IF( APFLAG ) THEN
            MESG = 'Enter logical name for a file of airport elevation input file'
            APDEV = PROMPTFFILE( MESG, .TRUE., .TRUE., 'APRT_ELEVATION',PROGNAME )
        END IF

C.........  Determine the year and projection status of the inventory
        CALL CHECK_INVYEAR( MENAME, MPRJFLAG, FDESC3D )

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
        METNAME = PROMPTMFILE(
     &      'Enter logical name for the METCRO3D meteorology file', 
     &      FSREAD3, 'MET_CRO_3D', PROGNAME )
 
        IF( .NOT. DESC3( METNAME ) ) THEN
            MESG = 'Could not get description of file "' //
     &             METNAME( 1:LEN_TRIM( METNAME ) ) // '" '
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

C.........  Check the grid definition
        CALL CHKGRID( 'mobile', 'GRID', 0, EFLAG )
    
C.........  Get file name for inventory pollutants codes/names
        MESG = 'Enter logical name for INVENTORY DATA TABLE file'
        PDEV = PROMPTFFILE( MESG, .TRUE., .TRUE., 'INVTABLE',
     &                      PROGNAME )

C.........  Get country, state, and county names no matter what, because it is
C           needed to allocate memory for the state and county totals, even
C           when they aren't going to be output
        IF( .NOT. USEEXPGEO() ) THEN
            CDEV = PROMPTFFILE(
     &             'Enter logical name for COUNTRY, STATE, AND ' //
     &             'COUNTY file', .TRUE., .TRUE., 'COSTCY', PROGNAME )
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

C.........  Open Met input files
        MESG = 'Enter name for CROSS-POINT GRID file'
        MGRNAME = PROMPTMFILE( MESG, FSREAD3, 'GRD_CRO_2D', PROGNAME )

        MESG = 'Enter name for CROSS-POINT LAYERED MET file'
        METNAME = PROMPTMFILE( MESG, FSREAD3, 'MET_CRO_3D', PROGNAME )

C.........  Read description of 3d file for defining layer structure
        IF( .NOT. DESC3( METNAME ) ) THEN
             MESG = 'Could not get description of file "' //
     &               TRIM( METNAME ) // '"'
             CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        ENDIF

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

C******************  INTERNAL SUBPROGRAMS  *****************************

C  Read and store airport elevations to compute AGL from MSL height in unit of feet

            SUBROUTINE READ_AIRPORT_ELEVATION

C............   Local variables
            INTEGER         I, N, IREC            ! indices and counters

            INTEGER      :: NLINES = 0            ! number of lines in input file

            CHARACTER(256)  LINE                  ! Read buffer for a line
            CHARACTER(300)  MESG                  ! Message buffer
            CHARACTER(16 ):: SEGMENT( 3 ) = ' '   ! line parsing array

C......................................................................

C.............  Get the number of lines
            NLINES = GETFLINE( APDEV, 'Airport elevation input file' )

C..............  Count no of sources
            N = 0
            IREC  = 0

            DO I = 1, NLINES

                READ ( APDEV, 93000, IOSTAT=IOS ) LINE
                IREC = IREC + 1

                IF ( IOS .GT. 0 ) THEN
                    WRITE( MESG, 94010)
     &                   'I/O error', IOS, 'reading APRT_ELEVATION '//
     &                   'description file at line', IREC
                    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                END IF

C.................  Left adjust line
                LINE = TRIM( LINE )

C.................  Skip blank and comment lines
                IF( BLKORCMT( LINE ) ) CYCLE

                N = N + 1

            END DO    ! end of loop

            IF( N == 0 ) THEN
                MESG = 'ERROR: No entries of APRT_ELEVATION'
                CALL M3MSG2( MESG )
            END IF

            NAPRT = N

            REWIND( APDEV )

C...............  Determine number of lines in filelist; this will be the maximum
C                 number of airport sources
            ALLOCATE( APRT_CODE( NAPRT ), STAT=IOS )
            CALL CHECKMEM( IOS, 'APRT_CODE', PROGNAME )
            ALLOCATE( APRT_ELEV( NAPRT ), STAT=IOS )
            CALL CHECKMEM( IOS, 'APRT_ELEV', PROGNAME )

            APRT_CODE = ' '
            APRT_ELEV = 0.0

            N = 0
            IREC  = 0

            DO I = 1, NLINES

                READ ( APDEV, 93000, IOSTAT=IOS ) LINE
                IREC = IREC + 1

                IF ( IOS .GT. 0 ) THEN
                    WRITE( MESG, 94010)
     &                'I/O error', IOS, 'reading APRT_ELEVATION '//
     &                'description file at line', IREC
                    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                END IF

C.................  Left adjust line
                LINE = TRIM( LINE )

C.................  Skip blank and comment lines
                IF( BLKORCMT( LINE ) ) CYCLE

C.................  Get line
                CALL PARSLINE( LINE, 3, SEGMENT )

                N = N + 1
                APRT_CODE( N ) = SEGMENT( 1 )
                APRT_ELEV( N ) = STR2REAL( SEGMENT( 3 ) )  * FT2M

             END DO    ! end of loop

             NAPRT = N

             CLOSE( APDEV )

C...................  FORMAT  STATEMENTS   ............................

C.............  Formatted file I/O formats...... 93xxx
93000        FORMAT( A )

C..............  Internal buffering formats...... 94xxx
94010        FORMAT( 10 ( A, :, I8, :, 2X  ) )

             END SUBROUTINE READ_AIRPORT_ELEVATION


        END SUBROUTINE OPENMRGIN
