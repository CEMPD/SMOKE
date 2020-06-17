        PROGRAM MRGELEV

C***********************************************************************
C  program body starts at line 218
C
C  DESCRIPTION:
C       This program combines ASCII elevated files produced by 
C       Smkmerge. It can optionally use a corresponding list of elevated
C       point source files to flag combined PinG sources. The time period
C       of the output file is adjusted based on the latest starting file
C       and earliest ending file, unless MRG_DIFF_DAY is set in which case
C       the time period is based on the standard environment variables.
C
C  PRECONDITIONS REQUIRED:
C       ASCII elevated files created by Smkmerge
C       Point source elevated files created by Elevpoint
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C       INITEM
C       GETM3EPI
C       M3MESG, M3MSG2, M3EXIT
C       CHECKMEM
C       SORTIC, SORTR1
C       SETOUTDATE
C       NEXTIME
C       WRDAYMSG
C
C  REVISION HISTORY:
C       Created 4/2005 by C. Seppanen
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
C***********************************************************************

        IMPLICIT NONE
        
C.........  INCLUDES
        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
        
C.........  EXTERNAL FUNCTIONS
        LOGICAL       BLKORCMT
        CHARACTER(2)  CRLF
        LOGICAL       ENVYN
        INTEGER       FINDC
        INTEGER       GETFLINE
        CHARACTER(10) HHMMSS
        INTEGER       JUNIT
        INTEGER       PROMPTFFILE
        INTEGER       SECSDIFF
        INTEGER       SEC2TIME

        EXTERNAL      BLKORCMT, CRLF, ENVYN, FINDC, GETFLINE, HHMMSS,
     &                JUNIT, PROMPTFFILE, SECSDIFF, SEC2TIME

C.........  LOCAL PARAMETERS
        CHARACTER(50), PARAMETER :: CVSW = '$Name SMOKEv4.8_Jun2020$'  ! CVS release tag
        INTEGER,       PARAMETER :: MXPING = 300        ! final number of PinG sources

C.........  LOCAL VARIABLES

C.........  Allocatable arrays
        INTEGER, ALLOCATABLE :: FILEDEV( :,: )      ! input file device numbers
        INTEGER, ALLOCATABLE :: FILELNS( : )        ! num. lines per PELV file
        INTEGER, ALLOCATABLE :: NSRCS( : )          ! num. srcs per file
        INTEGER, ALLOCATABLE :: SDATE( : )          ! start date by file
        INTEGER, ALLOCATABLE :: STIME( : )          ! start time by file
        INTEGER, ALLOCATABLE :: NSTEP( : )          ! no. time steps by file
        INTEGER, ALLOCATABLE :: SRCIDX( : )         ! index for sorting sources
        INTEGER, ALLOCATABLE :: EMISIDX( : )        ! index for sorting emissions
        INTEGER, ALLOCATABLE :: SAVLNUM( : )        ! current line number for each file
        
        CHARACTER(256), ALLOCATABLE :: FILENAM( : ) ! input file names
        CHARACTER(10),  ALLOCATABLE :: SPCNAM( : )  ! model species names
        CHARACTER(FPLLEN3+CHRLEN3), ALLOCATABLE :: PELVSRCA(:) ! unsorted source list
        CHARACTER(FPLLEN3+CHRLEN3), ALLOCATABLE :: PELVSRC(:)  ! sorted source list
        
        REAL, ALLOCATABLE :: PELVEMISA( : )      ! unsorted source emissions
        REAL, ALLOCATABLE :: PELVEMIS( : )       ! sorted source emissions
        
        LOGICAL, ALLOCATABLE :: PINGFND( : )     ! true: PinG source was found
        LOGICAL, ALLOCATABLE :: ISPING( : )      ! true: source remains PinG
        LOGICAL, ALLOCATABLE :: USEFIRST( : )    ! true: use first time step of file

C.........  File units and logical names
        INTEGER          LDEV       ! unit for log file
        INTEGER          IDEV       ! unit for list of elevated files
        INTEGER          PDEV       ! unit for list of PELV files 
        INTEGER          ODEV       ! unit for ASCII output file
        INTEGER          BDEV       ! unit for binary output file
        INTEGER          RPTDEV     ! unit for merge report file
        INTEGER          RDEV       ! temporary unit for reading input
        INTEGER          TDEV       ! temporary unit for input files

C.........  ASCII elevated file header variables
        CHARACTER(IOULEN3) :: GRDENV,  TMPGRDENV  ! gridded ouput units
        CHARACTER(60)      :: FNOTE,   TMPNOTE    ! file header note
        INTEGER            :: NMSPC,   TMPNMSPC   ! number of model species
        INTEGER            :: TMPNOUT             ! number of output sources
        INTEGER            :: NPARAM,  TMPNPARAM  ! number of control pkt parameters
        INTEGER            :: PDEVOUT, TMPPDEVOUT ! PTSRCE output unit number
        INTEGER            :: GSWITCH, TMPGSWITCH ! output grid switch
        INTEGER            :: USWITCH, TMPUSWITCH ! units table switch
        INTEGER            :: LSWITCH, TMPLSWITCH ! source locs table switch
        INTEGER            :: MSWITCH, TMPMSWITCH ! methods table switch
        INTEGER            :: VSWITCH, TMPVSWITCH ! values table switch
        INTEGER            :: ESWITCH, TMPESWITCH ! vertical methods table switch
        INTEGER            :: DDEVOUT, TMPDDEVOUT ! diffbreak file unit number
        INTEGER            :: RDEVOUT, TMPRDEVOUT ! regiontop file unit number
        INTEGER            :: TDEVOUT, TMPTDEVOUT ! temperature file unit number
        INTEGER            :: MDEVOUT, TMPMDEVOUT ! metscalars file unit number
        INTEGER            :: WDEVOUT, TMPWDEVOUT ! wind file unit number
        CHARACTER(10)      :: TMPNAM              ! model species name
        
        INTEGER            :: P_ALPHA, TMPALPHA   ! alpha value of projection
        REAL               :: XORIG,   TMPXORIG   ! grid x-origin
        REAL               :: YORIG,   TMPYORIG   ! grid y-origin
        REAL               :: XCELL,   TMPXCELL   ! x-cell size
        REAL               :: YCELL,   TMPYCELL   ! y-cell size
        INTEGER            :: NCOLS,   TMPNCOLS   ! number of columns
        INTEGER            :: NROWS,   TMPNROWS   ! number of rows
        INTEGER            :: NULAYS,  TMPNULAYS  ! number of model layers
        INTEGER            :: NZLOWR,  TMPNZLOWR  ! layers below diffbreak
        INTEGER            :: NZUPPR,  TMPNZUPPR  ! layers above diffbreak
        REAL               :: HTSUR,   TMPHTSUR   ! height of surface layer
        REAL               :: HTLOWR,  TMPHTLOWR  ! min cell ht b/w sfc and diffbr
        REAL               :: HTUPPR,  TMPHTUPPR  ! min cell ht b/w diffbr and top

        CHARACTER(10)      :: VERTTYPE            ! vertical method type

C.........  Binary elevated file variables
        CHARACTER(10)   :: HDRKEY = 'PTSOURCE'
        INTEGER(4)      :: INTBIN(20)         ! storage for integers
        REAL(4)         :: REALBIN(20)        ! storage for reals
        CHARACTER(4)       SPECID(10)         ! species name
        
        REAL(4),    ALLOCATABLE :: SRCDEFS( :,: ) ! point source definitions
        REAL(4),    ALLOCATABLE :: SRCFLOW( : )   ! flow rate for point sources
        REAL(4),    ALLOCATABLE :: SRCHTS ( : )   ! plume height for sources
        REAL(4),    ALLOCATABLE :: SRCEMIS( : )   ! emissions for each source
        
        INTEGER(4), ALLOCATABLE :: SRCCOLS( : )   ! column for point sources
        INTEGER(4), ALLOCATABLE :: SRCROWS( : )   ! row for point sources
        INTEGER(4), ALLOCATABLE :: SRCLAYS( : )   ! layer for point sources
        
C.........  Report file variables
C..         Reused for each time step and species
        REAL, ALLOCATABLE :: RPTBYFILE ( : )   ! emissions by file
        REAL                 RPTASCII          ! emissions in ASCII file
        REAL                 RPTBINARY         ! emissions in binary file
C..         Summed over time steps
        REAL, ALLOCATABLE :: RPTALLFILE( :,: ) ! emissions by file and species
        REAL, ALLOCATABLE :: RPTALLASC ( : )   ! emissions in ASCII file by species
        REAL, ALLOCATABLE :: RPTALLBIN ( : )   ! emissions in binary file by species
        REAL                 RPTFILESUM        ! summed emissions across input files

C.........  Source information variables
        REAL             XCOORD               ! x-coordinate
        REAL             YCOORD               ! y-coordinate
        CHARACTER(10)    TFIP                 ! FIPS code
        CHARACTER(10)    FCID                 ! facility ID
        CHARACTER(10)    SKID                 ! stack ID

        REAL             STKHT                ! stack height
        REAL             STKDM                ! stack diameter
        REAL             STKTK                ! stack temperature
        REAL             STKVE                ! stack velocity

        CHARACTER(FIPLEN3)         FIP        ! FIPS code
        CHARACTER(PLTLEN3)         PLT        ! plant/facility code
        CHARACTER(CHRLEN3)         PNT        ! point/stack code
        CHARACTER(FPLLEN3+CHRLEN3) CSRC       ! FIP // PLT // CHAR1 string

C.........  Hourly emissions variables
        INTEGER          NUM                  ! source number
        CHARACTER(10)    VNAME                ! species name
        REAL             EMIS                 ! hourly emissions

C.........  Output time variables
        INTEGER ::       G_SDATE = 0          ! start date
        INTEGER ::       G_STIME = 0          ! start time
        INTEGER ::       G_NSTEPS = 1         ! number of time steps
        INTEGER ::       G_TSTEP = 0          ! time step

C.........  PinG source count variables
        INTEGER          NPINGSRC             ! no. sources in PELV files
        INTEGER          NPINGASCII           ! no. sources in ASCII elevated files
        INTEGER          NMISSPELV            ! no. sources missing from PELV files
        INTEGER          NMISSASCII           ! no. sources missing from ASCII files
        INTEGER          NPINGOUT             ! no. sources in output ASCII elevated file
        
C.........  Other local variables
        INTEGER          I, J, K, L           ! indexes and counters
        INTEGER          IDUM, IDUM2          ! dummy integers
        INTEGER          IOS                  ! i/o status
        INTEGER          JDATE, JTIME         ! date and time for looping
        INTEGER          LDATE                ! last date processed
        INTEGER          LNUM                 ! current line number
        INTEGER          MXFILES              ! maximum number of input files
        INTEGER          MXPFILES             ! maximum number of PELV files
        INTEGER          MXLINES              ! maximum number of lines
        INTEGER          NFILES               ! number of input files
        INTEGER          NTOTSRCS             ! total number of output sources
        INTEGER          IBD, IBT             ! start date and time in elev. format
        INTEGER          IED, IET             ! end date and time in elev. format
        INTEGER          TMPBD, TMPBT         ! tmp. start date/time in elev. format
        INTEGER          TMPED, TMPET         ! tmp. end date/time in elev. format
        INTEGER          RDATE                ! date to read data for
        INTEGER          SECS                 ! number of seconds
        INTEGER          TMPDATE              ! temporary date
        INTEGER          TMPTIME              ! temporary time
        INTEGER          TMPSTEP              ! temporary number of time steps
        INTEGER          BASENUM              ! base source number for hourly emissions
        INTEGER          ENUM, PNUM, GNUM     ! elevated, PinG, and group numbers
        INTEGER          ROW, COL             ! calculated row and column for source

        REAL             FDUM, FDUM2          ! dummy reals

        LOGICAL ::       EFLAG = .FALSE.      ! true: an error happened
        LOGICAL          FOUND                ! true: matching source was found
        LOGICAL          MRGDIFF              ! true: merge files from different days
        LOGICAL          MSGPRINT             ! true: print warning message
        LOGICAL          PINGFLAG             ! true: process PinG sources
        LOGICAL          LLGRID               ! true: grid is lat-lon
        LOGICAL          OUTFLAG              ! true: output merged ascii elev point src file

        CHARACTER(256)   DUMMY                ! dummy character string
        CHARACTER(30)    FMT                  ! output format for emissions
        CHARACTER(100)   RPTFMT               ! output format for report
        CHARACTER(100)   RPTALLFMT            ! output format for all time step report
        CHARACTER(256)   LINE                 ! input line
        CHARACTER(512)   MESG                 ! message buffer

        CHARACTER(16) :: PROGNAME = 'MRGELEV' ! program name

C***********************************************************************
C   begin body of program MRGELEV

        LDEV = INIT3()

C.........  Write out copyright, version, web address, header info, and 
C           prompt to continue running the program
        CALL INITEM( LDEV, CVSW, PROGNAME )

C.........  Open lists of input files
        MESG = 'Enter logical name for ASCII ELEVATED inputs list'
        IDEV = PROMPTFFILE( MESG, .TRUE., .TRUE., 'FILELIST', PROGNAME )

        MESG = 'Enter logical name for ELEVATED POINT SOURCE ' //
     &         'inputs list (or "NONE")'
        PDEV = PROMPTFFILE( MESG, .TRUE., .TRUE., 'PELVLIST', PROGNAME )
        IF( PDEV == -2 ) THEN
            PINGFLAG = .FALSE.
        ELSE
            PINGFLAG = .TRUE.
        END IF
        
C.........  Get environment variables
        MESG = 'Merge files from different days into single file'
        MRGDIFF = ENVYN( 'MRG_DIFF_DAYS', MESG, .FALSE., IOS )

        MESG = 'Output merged ASCII elevated point sources file or not'
        OUTFLAG = ENVYN( 'SMK_ASCIIELEV_YN', MESG, .TRUE., IOS )

C.........  Get date and time settings from environment        
        IF( MRGDIFF ) THEN
            CALL GETM3EPI( -1, G_SDATE, G_STIME, G_TSTEP, G_NSTEPS )
        END IF

C.........  Determine maximum number of input files in lists
        MXFILES  = GETFLINE( IDEV, 'List of files to merge' )
        
        IF( PINGFLAG ) THEN
            MXPFILES = GETFLINE( PDEV, 'List of PELV input files' )
        
            IF( MXFILES /= MXPFILES ) THEN
                MESG = 'ERROR: Number of ASCII elevated files does ' //
     &            'not match number of' // 
     &            CRLF() // BLANK10 // 'elevated point source files'
                CALL M3MSG2( MESG )
                
                MESG = 'Problem with input file list'
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF
        END IF
        
C.........  Allocate memory to store file info
        ALLOCATE( FILEDEV( MXFILES, 2 ), STAT=IOS )
        CALL CHECKMEM( IOS, 'FILEDEV', PROGNAME )
        ALLOCATE( FILENAM( MXFILES ), STAT=IOS )
        CALL CHECKMEM( IOS, 'FILENAM', PROGNAME )
        ALLOCATE( FILELNS( MXFILES ), STAT=IOS )
        CALL CHECKMEM( IOS, 'FILELNS', PROGNAME )

C.........  Read lists of input files
        RDEV = IDEV
        DO J = 1, 2
            IF( J == 1 ) THEN
                MESG = 'Reading list of ASCII elevated files...'
            ELSE
                IF( .NOT. PINGFLAG ) EXIT
                MESG = 'Reading list of elevated point source files...'
                MXLINES = 0
            END IF
            CALL M3MSG2( MESG )
            
            K = 0
            DO I = 1, MXFILES
            
                READ( RDEV, 93000, IOSTAT=IOS ) LINE
                
                IF( IOS /= 0 ) THEN
                    EFLAG = .TRUE.
                    WRITE( MESG,94010 ) 'ERROR: I/O error ', IOS,
     &                  'reading list of input files at line ', I
                    CALL M3MESG( MESG )
                    CYCLE
                END IF
            
C.................  Skip blank or comment lines
                IF( BLKORCMT( LINE ) ) CYCLE

C.................  When reading PELV list, allow file name to be NONE
                IF( J == 2 .AND. TRIM( LINE ) == 'NONE' ) THEN
                    K = K + 1
                    FILEDEV( K,J ) = 0
                    FILELNS( K ) = 0
                    CYCLE
                END IF

C.................  Open input file
                TDEV = JUNIT()
                OPEN( UNIT=TDEV, FILE=TRIM( LINE ), STATUS='OLD',
     &                ACTION='READ', IOSTAT=IOS )
            
                IF( IOS /= 0 ) THEN
                    EFLAG = .TRUE.
                    MESG = 'ERROR: Could not open input file ' // 
     &                     CRLF() // BLANK10 // TRIM( LINE )
                    CALL M3MESG( MESG )
                    CYCLE
                ELSE
                    K = K + 1
                    FILEDEV( K,J ) = TDEV
                    
                    IF( J == 1 ) THEN
                        FILENAM( K ) = TRIM( LINE )
                    END IF
                END IF
                
C.................  For PELV files, count number of lines
                IF( J == 2 ) THEN
                    FILELNS( K ) = 
     &                  GETFLINE( TDEV, 'Elevated point source file' )
                END IF
            END DO
            
            RDEV = PDEV        
        END DO

        IF( EFLAG ) THEN
            MESG = 'Problem opening input files'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

C.........  Store actual number of files opened
        NFILES = K

C.........  Build list of PinG sources from PELV files
        IF( PINGFLAG ) THEN
            MESG = 'Building list of PinG sources...'
            CALL M3MSG2( MESG )
        
            MXLINES = 0
            DO I = 1, NFILES
                MXLINES = MXLINES + FILELNS( I )
            END DO
    
            ALLOCATE( PELVSRCA( MXLINES ), STAT=IOS )
            CALL CHECKMEM( IOS, 'PELVSRCA', PROGNAME )
            ALLOCATE( PELVEMISA( MXLINES ), STAT=IOS )
            CALL CHECKMEM( IOS, 'PELVEMISA', PROGNAME )
    
            L = 0
            DO I = 1, NFILES
            
                TDEV = FILEDEV( I,2 )
                IF( TDEV == 0 ) CYCLE

C.................  Loop through lines in file
                DO J = 1, FILELNS( I )

C.....................  Read line as if all fields are present
                    READ( TDEV, 93095 ) ENUM, PNUM, GNUM,
     &                  FIP, PLT, PNT, EMIS

C.....................  For non-PinG lines, the elevated source
C                       number will be zero
                    IF( ENUM /= 0 ) CYCLE
            
                    PLT = ADJUSTL( PLT )
                    PNT = ADJUSTL( PNT )

C.....................  Build source string to match ASCII elevated file;
C                       plant name and characteristic are truncated to 10
                    CSRC = FIP // ADJUSTR( PLT( 1:10 ) ) // 
     &                     ADJUSTR( PNT( 1:10 ) )

C.....................  Check if we already have this source
                    FOUND = .FALSE.
                    DO K = 1, L
                        IF( PELVSRCA( K ) == CSRC ) THEN
                            PELVEMISA( K ) = PELVEMISA( K ) + EMIS
                            FOUND = .TRUE.
                            EXIT
                        END IF
                    END DO
                
                    IF( .NOT. FOUND ) THEN
                        L = L + 1
                        PELVSRCA( L ) = CSRC
                        PELVEMISA( L ) = EMIS
                    END IF
                
                END DO
            
            END DO

            NPINGSRC = L

C.............  Sort PinG list based on source information
            ALLOCATE( SRCIDX( NPINGSRC ), STAT=IOS )
            CALL CHECKMEM( IOS, 'SRCIDX', PROGNAME )
            ALLOCATE( EMISIDX( NPINGSRC ), STAT=IOS )
            CALL CHECKMEM( IOS, 'EMISIDX', PROGNAME )
            ALLOCATE( PELVSRC( NPINGSRC ), STAT=IOS )
            CALL CHECKMEM( IOS, 'PELVSRC', PROGNAME )
            ALLOCATE( PELVEMIS( NPINGSRC ), STAT=IOS )
            CALL CHECKMEM( IOS, 'PELVEMIS', PROGNAME )
            ALLOCATE( PINGFND( NPINGSRC ), STAT=IOS )
            CALL CHECKMEM( IOS, 'PINGFND', PROGNAME )
            ALLOCATE( ISPING( NPINGSRC ), STAT=IOS )
            CALL CHECKMEM( IOS, 'ISPING', PROGNAME )

            DO I = 1, NPINGSRC
                SRCIDX( I ) = I
                EMISIDX( I ) = I
            END DO
    
            CALL SORTIC( NPINGSRC, SRCIDX, PELVSRCA )
            
            DO I = 1, NPINGSRC
                PELVSRC( I ) = PELVSRCA( SRCIDX( I ) )
                PELVEMIS( I ) = PELVEMISA( SRCIDX( I ) )
            END DO
            
            DEALLOCATE( SRCIDX, PELVSRCA, PELVEMISA )
            
C.............  Create second sort index for emissions; this will be
C               used for matching the top 300 sources
            CALL SORTR1( NPINGSRC, EMISIDX, PELVEMIS )
            
            ISPING = .FALSE.
            DO I = NPINGSRC, MAX( NPINGSRC - MXPING + 1,1 ), -1
                ISPING( EMISIDX( I ) ) = .TRUE.
            END DO
        END IF

C.........  Allocate memory for file header information
        ALLOCATE( NSRCS( NFILES ), STAT=IOS )
        CALL CHECKMEM( IOS, 'NSRCS', PROGNAME )
        ALLOCATE( SDATE( NFILES ), STAT=IOS )
        CALL CHECKMEM( IOS, 'SDATE', PROGNAME )
        ALLOCATE( STIME( NFILES ), STAT=IOS )
        CALL CHECKMEM( IOS, 'STIME', PROGNAME )
        ALLOCATE( NSTEP( NFILES ), STAT=IOS )
        CALL CHECKMEM( IOS, 'EDATE', PROGNAME )
        ALLOCATE( SAVLNUM( NFILES ), STAT=IOS )
        CALL CHECKMEM( IOS, 'SAVLNUM', PROGNAME )
        
        IF( MRGDIFF ) THEN
            ALLOCATE( USEFIRST( NFILES ), STAT=IOS )
            CALL CHECKMEM( IOS, 'USEFIRST', PROGNAME )
        END IF

C.........  Check the header of each input file
        MESG = 'Checking headers of input files...'
        CALL M3MSG2( MESG )

        GRDENV = ''
        DO I = 1, NFILES
        
            TDEV = FILEDEV( I,1 )
            LNUM = 0
            
            READ( TDEV, 93010 ) DUMMY, TMPGRDENV
            LNUM = LNUM + 1
            CALL CHECK_HEADER( DUMMY, 'CONTROL', FILENAM( I ), LNUM )

            READ( TDEV, 93000 ) DUMMY
            LNUM = LNUM + 1
            CALL CHECK_HEADER( DUMMY, 'PTSOURCE', FILENAM( I ), LNUM )
            
            READ( TDEV, 93000 ) TMPNOTE
            LNUM = LNUM + 1
            
            READ( TDEV, 93015 ) TMPNMSPC, IDUM, TMPNOUT, 
     &                          IDUM2, TMPNPARAM
            LNUM = LNUM + 1
            IF( IDUM /= 0 .OR. IDUM2 /= 1 ) THEN
                CALL WRITE_HEADER_ERROR( FILENAM( I ), LNUM )
            END IF
            
            READ( TDEV, 93015 ) TMPPDEVOUT, IDUM, TMPGSWITCH
            LNUM = LNUM + 1
            IF( IDUM /= 0 ) THEN
                CALL WRITE_HEADER_ERROR( FILENAM( I ), LNUM )
            END IF
            
            READ( TDEV, 93015 ) TMPUSWITCH, TMPLSWITCH, IDUM, 
     &                          TMPMSWITCH, TMPVSWITCH
            LNUM = LNUM + 1
            IF( IDUM /= 0 ) THEN
                CALL WRITE_HEADER_ERROR( FILENAM( I ), LNUM )
            END IF
            
            READ( TDEV, 93015 ) IDUM, IDUM2, TMPESWITCH
            LNUM = LNUM + 1
            IF( IDUM /= 1 .OR. IDUM2 /= 0 ) THEN
                CALL WRITE_HEADER_ERROR( FILENAM( I ), LNUM )
            END IF
            
            READ( TDEV, 93015 ) TMPDDEVOUT, TMPRDEVOUT, IDUM, 
     &                          TMPTDEVOUT, TMPMDEVOUT, TMPWDEVOUT
            LNUM = LNUM + 1
            IF( IDUM /= 0 ) THEN
                CALL WRITE_HEADER_ERROR( FILENAM( I ), LNUM )
            END IF

C.............  Store variables the first time, otherwise compare to previous            
            IF( GRDENV == '' ) THEN
                GRDENV  = TMPGRDENV
                FNOTE   = TMPNOTE
                NMSPC   = TMPNMSPC
                NPARAM  = TMPNPARAM
                PDEVOUT = TMPPDEVOUT
                GSWITCH = TMPGSWITCH
                USWITCH = TMPUSWITCH
                LSWITCH = TMPLSWITCH
                MSWITCH = TMPMSWITCH
                VSWITCH = TMPVSWITCH
                ESWITCH = TMPESWITCH
                DDEVOUT = TMPDDEVOUT
                RDEVOUT = TMPRDEVOUT
                TDEVOUT = TMPTDEVOUT
                MDEVOUT = TMPMDEVOUT
                WDEVOUT = TMPWDEVOUT
            ELSE
                IF( TMPGRDENV  /= GRDENV  .OR.
     &              TMPNOTE    /= FNOTE   .OR.
     &              TMPNMSPC   /= NMSPC   .OR.
     &              TMPNPARAM  /= NPARAM  .OR.
     &              TMPPDEVOUT /= PDEVOUT .OR.
     &              TMPGSWITCH /= GSWITCH .OR.
     &              TMPUSWITCH /= USWITCH .OR.
     &              TMPLSWITCH /= LSWITCH .OR.
     &              TMPMSWITCH /= MSWITCH .OR.
     &              TMPVSWITCH /= VSWITCH .OR.
     &              TMPESWITCH /= ESWITCH .OR.
     &              TMPDDEVOUT /= DDEVOUT .OR.
     &              TMPRDEVOUT /= RDEVOUT .OR.
     &              TMPTDEVOUT /= TDEVOUT .OR.
     &              TMPMDEVOUT /= MDEVOUT .OR.
     &              TMPWDEVOUT /= WDEVOUT      ) THEN
                    CALL WRITE_HEADER_DIFF( FILENAM( I ) )
                END IF
            END IF
            
            NSRCS( I ) = TMPNOUT
            SAVLNUM( I ) = LNUM

        END DO     

C.........  Read and check species names for each file
        ALLOCATE( SPCNAM( NMSPC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'SPCNAM', PROGNAME )
        SPCNAM = ' '

        DO I = 1, NFILES
        
            TDEV = FILEDEV( I,1 )
            LNUM = SAVLNUM( I )
        
            DO J = 1, NMSPC
                READ( TDEV, 93020 ) TMPNAM
                LNUM = LNUM + 1
                
                IF( I == 1 ) THEN
                    SPCNAM( J ) = TMPNAM
                ELSE
                    IF( TMPNAM /= SPCNAM( J ) ) THEN
                        CALL WRITE_HEADER_DIFF( FILENAM( I ) )
                    END IF
                END IF
            END DO
            
            SAVLNUM( I ) = LNUM
        END DO

C.........  Continue checking file headers
        P_ALPHA = 0
        DO I = 1, NFILES
        
            TDEV = FILEDEV( I,1 )
            LNUM = SAVLNUM( I )

C.............  Read file start and end dates and times
            READ( TDEV, 93030 ) IBD, IBT, IED, IET
            LNUM = LNUM + 1
        
            IF( IBD > 70000 ) THEN
                SDATE( I ) = IBD + 1900000
            ELSE
                SDATE( I ) = IBD + 2000000
            END IF
            
            STIME( I ) = IBT * 100
            
            IF( IED > 70000 ) THEN
                TMPDATE = IED + 1900000
            ELSE
                TMPDATE = IED + 2000000
            END IF
        
            TMPTIME = IET * 100
            
            SECS = SECSDIFF( SDATE( I ), STIME( I ), TMPDATE, TMPTIME )
            NSTEP( I ) = 1 + (SECS / 3600)
            
            READ( TDEV, 93000 ) DUMMY
            LNUM = LNUM + 1
            CALL CHECK_HEADER( DUMMY, 'END', FILENAM( I ), LNUM )

C.............  Read grid and layer information            
            READ( TDEV, 93000 ) DUMMY
            LNUM = LNUM + 1
            CALL CHECK_HEADER( DUMMY, 'REGION', FILENAM( I ), LNUM )
            
            READ( TDEV, 93040 ) FDUM, FDUM2, TMPALPHA
            LNUM = LNUM + 1
            IF( FDUM /= 0. .OR. FDUM2 /= 0. ) THEN
                CALL WRITE_HEADER_ERROR( FILENAM( I ), LNUM )
            END IF
            
            READ( TDEV, 93040 ) TMPXORIG, TMPYORIG
            LNUM = LNUM + 1

C.............  Format for x- and y-cell sizes changes depending on if grid
C               is lat-lon; we don't have to worry about the difference on input
C               because Fortran ignores number of decimal places when input
C               has a decimal point
            READ( TDEV, 93040 ) TMPXCELL, TMPYCELL
            LNUM = LNUM + 1
            
            READ( TDEV, 93050 ) TMPNCOLS, TMPNROWS, TMPNULAYS
            LNUM = LNUM + 1
            READ( TDEV, 93060 ) TMPNZLOWR, TMPNZUPPR, TMPHTSUR, 
     &                          TMPHTLOWR, TMPHTUPPR
            LNUM = LNUM + 1
            
            IF( P_ALPHA == 0 ) THEN
                P_ALPHA = TMPALPHA
                XORIG   = TMPXORIG
                YORIG   = TMPYORIG
                XCELL   = TMPXCELL
                YCELL   = TMPYCELL
                NCOLS   = TMPNCOLS
                NROWS   = TMPNROWS
                NULAYS  = TMPNULAYS
                NZLOWR  = TMPNZLOWR
                NZUPPR  = TMPNZUPPR
                HTSUR   = TMPHTSUR
                HTLOWR  = TMPHTLOWR
                HTUPPR  = TMPHTUPPR
            ELSE
                IF( TMPALPHA  /= P_ALPHA .OR.
     &              TMPXORIG  /= XORIG   .OR.
     &              TMPYORIG  /= YORIG   .OR.
     &              TMPXCELL  /= XCELL   .OR.
     &              TMPYCELL  /= YCELL   .OR.
     &              TMPNCOLS  /= NCOLS   .OR.
     &              TMPNROWS  /= NROWS   .OR.
     &              TMPNULAYS /= NULAYS  .OR.
     &              TMPNZLOWR /= NZLOWR  .OR.
     &              TMPNZUPPR /= NZUPPR  .OR.
     &              TMPHTSUR  /= HTSUR   .OR.
     &              TMPHTLOWR /= HTLOWR  .OR.
     &              TMPHTUPPR /= HTUPPR       ) THEN
                    CALL WRITE_HEADER_DIFF( FILENAM( I ) )
                END IF
            END IF

            READ( TDEV, 93000 ) DUMMY
            LNUM = LNUM + 1
            CALL CHECK_HEADER( DUMMY, 'END', FILENAM( I ), LNUM )
            
            READ( TDEV, 93000 ) DUMMY
            LNUM = LNUM + 1
            CALL CHECK_HEADER( DUMMY, 'POINT SOURCES', FILENAM( I ), LNUM )

            SAVLNUM( I ) = LNUM
        END DO

C.........  Check if grid is lat-lon; if not lat-lon, then x-cell size won't
C           fit into F10.7 format
        LLGRID = .TRUE.
        WRITE( MESG,'( F10.7 )' ) XCELL
        IF( MESG( 1:1 ) == '*' ) THEN
            LLGRID = .FALSE.
        END IF

C.........  Calculate total number of output sources
        NTOTSRCS = 0
        DO I = 1, NFILES
            NTOTSRCS = NTOTSRCS + NSRCS( I )
        END DO

C.........  Check that environment settings are consistent with files
        IF( MRGDIFF ) THEN
            IF( G_TSTEP /= 10000 ) THEN
                MESG = 'Output time step must be 10000 for hourly data'
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF
        END IF

C.........  Determine output date, time, and number of time steps
        CALL SETOUTDATE( G_SDATE, G_STIME, G_NSTEPS, NFILES, SDATE,
     &                   STIME, NSTEP, FILENAM, MRGDIFF, USEFIRST )                         

        TMPDATE = G_SDATE
        TMPTIME = G_STIME
        CALL NEXTIME( TMPDATE, TMPTIME, ( G_NSTEPS - 1 ) * 10000 )
        
        IBD = REMOVE_4DIGIT_YEAR( G_SDATE )
        IED = REMOVE_4DIGIT_YEAR( TMPDATE )
        
        IBT = G_STIME / 100
        IET = TMPTIME / 100

C.........  Open output files
        IF( OUTFLAG ) THEN
            MESG = 'Enter logical name for output ' //
     &             'ASCII ELEVATED SOURCES file'
            ODEV = PROMPTFFILE( MESG, .FALSE., .TRUE., 'OUTFILE', PROGNAME )
        END IF

        MESG = 'Enter logical name for output ' //
     &         'BINARY ELEVATED SOURCES file'
        BDEV = PROMPTFFILE( MESG, .FALSE., .FALSE., 'PTSOURCE', PROGNAME )

C.........  Open report file
        MESG = 'Enter logical name for the MRGELEV REPORT file'
        RPTDEV = PROMPTFFILE( MESG, .FALSE., .TRUE., 'REPMRGELEV', 
     &                        PROGNAME )

C.........  Write header to output files
        IF( OUTFLAG ) CALL WRITE_ASCII_HEADER
        CALL WRITE_BINARY_HEADER
        CALL WRITE_REPORT_HEADER

C.........  Build report format string
        RPTFMT = "(I5, ';', 1X, I4, ';', 1X, A10, ';', "
        WRITE( RPTFMT, '(A,I2)' ) TRIM( RPTFMT ), NFILES + 3
        RPTFMT = TRIM( RPTFMT ) // "(1X, E12.5, ';'))"
        
        RPTALLFMT = "('All  ; All ; ', A10, ';',"
        WRITE( RPTALLFMT, '(A,I2)' ) TRIM( RPTALLFMT ), NFILES + 3
        RPTALLFMT = TRIM( RPTALLFMT ) // "(1X, E12.5, ';'))"

C.........  Allocate space for source information
        ALLOCATE( SRCDEFS( 6,NTOTSRCS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'SRCDEFS', PROGNAME )
        ALLOCATE( SRCCOLS( NTOTSRCS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'SRCCOLS', PROGNAME )
        ALLOCATE( SRCROWS( NTOTSRCS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'SRCROWS', PROGNAME )
        ALLOCATE( SRCLAYS( NTOTSRCS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'SRCLAYS', PROGNAME )
        ALLOCATE( SRCFLOW( NTOTSRCS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'SRCFLOW', PROGNAME )
        ALLOCATE( SRCHTS( NTOTSRCS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'SRCHTS', PROGNAME )
        ALLOCATE( SRCEMIS( NTOTSRCS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'SRCEMIS', PROGNAME )
        
        SRCEMIS = 0.    ! array

C.........  Read and output source information for each file
        K = 0
        NPINGASCII = 0
        NMISSPELV = 0
        NPINGOUT = 0
        MSGPRINT = .TRUE.
        DO I = 1, NFILES
        
            TDEV = FILEDEV( I,1 )
            LNUM = SAVLNUM( I )
        
            DO J = 1, NSRCS( I )
                READ( TDEV, 93070 ) IDUM, DUMMY, XCOORD, YCOORD,
     &              FCID, SKID, TFIP
                LNUM = LNUM + 1
                CALL CHECK_HEADER( DUMMY, 'STD', FILENAM( I ), LNUM )

C.................  Read stack parameters                    
                READ( TDEV, 93080 ) STKHT, STKDM, STKTK, STKVE
                LNUM = LNUM + 1
                
C.................  Check for source in PinG list
                IF( STKDM < 0 .AND. PINGFLAG ) THEN
                    NPINGASCII = NPINGASCII + 1
                    
                    TFIP = ADJUSTR( TFIP )
                
                    CSRC = TFIP( 5:10 ) // FCID // SKID
                
                    L = FINDC( CSRC, NPINGSRC, PELVSRC )
                    
                    IF( L < 1 ) THEN
                        NMISSPELV = NMISSPELV + 1
                        
                        IF( MSGPRINT ) THEN
                            MESG = 'WARNING: The following PinG ' //
     &                        'sources are in the ASCII elevated ' //
     &                        'files but' // CRLF() // BLANK10 // 
     &                        'are not in the corresponding PELV ' //
     &                        'files; the sources will remain as ' //
     &                        CRLF() // BLANK10 // 'PinG sources in ' //
     &                        'the output file.'
                            CALL M3MESG( MESG )
                            MSGPRINT = .FALSE.
                        END IF
                        
                        MESG = 'FIPS: ' // TFIP( 5:10 ) //
     &                         '    Facility: ' // FCID //
     &                         '    Stack: ' // SKID
                        CALL M3MESG( MESG )
                    ELSE
                        PINGFND( L ) = .TRUE.

C.........................  Check if source remains as PinG source
                        IF( .NOT. ISPING( L ) ) THEN
                            STKDM = -STKDM
                        ELSE
                            NPINGOUT = NPINGOUT + 1
                        END IF
                    END IF
                END IF
                
                K = K + 1

C.................  Write source to ASCII file
                IF( OUTFLAG ) THEN
                    IF( LLGRID ) THEN
                        WRITE( ODEV, 93065 ) K, 'STD       ', XCOORD,
     &                      YCOORD, FCID, SKID, TFIP
                    ELSE
                        WRITE( ODEV, 93070 ) K, 'STD       ', XCOORD,
     &                      YCOORD, FCID, SKID, TFIP
                    END IF
                    WRITE( ODEV, 93080 ) STKHT, STKDM, STKTK, STKVE
                END IF

C.................  Save source to write to binary file
                SRCDEFS( 1,K ) = XCOORD     ! x-coord of source
                SRCDEFS( 2,K ) = YCOORD     ! y-coord of source
                SRCDEFS( 3,K ) = STKHT      ! stack height
                SRCDEFS( 4,K ) = STKDM      ! stack diameter
                SRCDEFS( 5,K ) = STKTK      ! stack temperature
                SRCDEFS( 6,K ) = STKVE      ! stack exit velocity
                
                COL = 1 + INT( ( XCOORD - XORIG ) / XCELL ) ! src column
                IF( COL > NCOLS .OR. COL <= 0 ) THEN
                    SRCCOLS( K ) = 0
                ELSE
                    SRCCOLS( K ) = COL
                END IF
                
                ROW = 1 + INT( ( YCOORD - YORIG ) / YCELL ) ! src row
                IF( ROW > NROWS .OR. ROW <= 0 ) THEN
                    SRCROWS( K ) = 0
                ELSE
                    SRCROWS( K ) = ROW
                END IF
                
                SRCLAYS( K ) = 1            ! layer of source
                SRCFLOW( K ) = -9.          ! flow rate
                SRCHTS ( K ) = STKHT        ! plume height
            END DO
            
            READ( TDEV, 93000 ) DUMMY
            LNUM = LNUM + 1
            CALL CHECK_HEADER( DUMMY, 'END', FILENAM( I ), LNUM )

            SAVLNUM( I ) = LNUM

        END DO
        
        IF( OUTFLAG ) WRITE( ODEV, 93000 ) 'END'

C.........  Write point source definition record to binary file
        WRITE( BDEV ) SRCDEFS

C.........  Allocate space for report information
        ALLOCATE( RPTBYFILE( NFILES ), STAT=IOS )
        CALL CHECKMEM( IOS, 'RPTBYFILE', PROGNAME )
        ALLOCATE( RPTALLFILE( NMSPC, NFILES ), STAT=IOS )
        CALL CHECKMEM( IOS, 'RPTALLFILE', PROGNAME )
        ALLOCATE( RPTALLASC( NMSPC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'RPTALLASC', PROGNAME )
        ALLOCATE( RPTALLBIN( NMSPC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'RPTALLBIN', PROGNAME )

C.........  Read and output hourly emissions for each file
        LDATE = 0
        JDATE = G_SDATE
        JTIME = G_STIME
        
        RPTALLFILE = 0.   ! array
        RPTALLASC  = 0.   ! array
        RPTALLBIN  = 0.   ! array

        DO L = 1, G_NSTEPS

C.............  Write message for new day        
            IF( JDATE /= LDATE ) THEN
                CALL WRDAYMSG( JDATE, MESG )
            END IF

C.............  Write message for each hour
            WRITE( MESG,94020 ) HHMMSS( JTIME )
            CALL M3MSG2( MESG )
        
            DO I = 1, NFILES
            
                TDEV = FILEDEV( I,1 )
                LNUM = SAVLNUM( I )
                
C.................  Set read date
                IF( MRGDIFF .AND. USEFIRST( I ) ) THEN
                    IDUM = 0
                    TMPSTEP = SEC2TIME( 
     &                          SECSDIFF( 
     &                            G_SDATE, IDUM, JDATE, IDUM ) )
                    RDATE = SDATE( I )
                    CALL NEXTIME( RDATE, IDUM, TMPSTEP )
                ELSE
                    RDATE = JDATE
                END IF

C.................  Convert dates and times to elevated format
                TMPBD = REMOVE_4DIGIT_YEAR( RDATE )
                TMPBT = JTIME / 100

C.................  Calculate end date and time                
                TMPDATE = RDATE
                TMPTIME = JTIME
                CALL NEXTIME( TMPDATE, TMPTIME, 10000 )
                
                TMPED = REMOVE_4DIGIT_YEAR( TMPDATE )
                TMPET = TMPTIME / 100
                                
C.................  If not the first time step, read end of previous emissions section                
                IF( L /= 1 ) THEN
                    READ( TDEV, 93000 ) DUMMY
                    LNUM = LNUM + 1
                    CALL CHECK_HEADER( DUMMY, 'END', 
     &                                 FILENAM( I ), LNUM )
                    
                    READ( TDEV, 93000 ) DUMMY
                    LNUM = LNUM + 1
                    CALL CHECK_HEADER( DUMMY, 'ENDTIME', 
     &                                 FILENAM( I ), LNUM )
                END IF
                
                READ( TDEV, 93000 ) DUMMY
                LNUM = LNUM + 1
                CALL CHECK_HEADER( DUMMY, 'TIME INTERVAL', 
     &                             FILENAM( I ), LNUM )
                
                READ( TDEV, 93050 ) IBD, IBT, IED, IET
                LNUM = LNUM + 1

C.................  During first time step, skip ahead as needed to correct
C                   place in file; after first time step, all reads should be
C                   sequential so no skipping is needed
                IF( L == 1 ) THEN
                    DO WHILE( IBD /= TMPBD .AND. IBT /= TMPBT )

C.........................  Skip through data for wrong time steps
                        DO J = 1, 8
                            READ( TDEV, * ) DUMMY
                            LNUM = LNUM + 1
                        END DO

                        DO
                            READ( TDEV, * ) DUMMY
                            LNUM = LNUM + 1
                            IF( TRIM( DUMMY ) == 'END' ) EXIT
                        END DO
                        
                        DO J = 1, 2
                            READ( TDEV, * ) DUMMY
                            LNUM = LNUM + 1
                        END DO

                        READ( TDEV, 93050 ) IBD, IBT, IED, IET
                        LNUM = LNUM + 1
                    END DO
                ELSE
                    IF( IBD /= TMPBD .OR. IBT /= TMPBT .OR.
     &                  IED /= TMPED .OR. IET /= TMPET      ) THEN
                        CALL WRITE_HEADER_ERROR( FILENAM( I ), LNUM )
                    END IF
                END IF

C.................  Set output start date, start time, end date, and end time
                IBD = REMOVE_4DIGIT_YEAR( JDATE )
                IBT = JTIME / 100
                
                TMPDATE = JDATE
                TMPTIME = JTIME
                CALL NEXTIME( TMPDATE, TMPTIME, 10000 )
                
                IED = REMOVE_4DIGIT_YEAR( TMPDATE )
                IET = TMPTIME / 100
                
                READ( TDEV, 93000 ) DUMMY
                LNUM = LNUM + 1
                CALL CHECK_HEADER( DUMMY, 'METHOD', FILENAM( I ), LNUM )
                
                READ( TDEV, 93000 ) DUMMY
                LNUM = LNUM + 1
                CALL CHECK_HEADER( DUMMY, 'STD       ALL       ' //
     &              'EMVALUES  0.        50000.', FILENAM( I ), LNUM )
                
                READ( TDEV, 93000 ) DUMMY
                LNUM = LNUM + 1
                CALL CHECK_HEADER( DUMMY, 'END', FILENAM( I ), LNUM )
                
                READ( TDEV, 93000 ) DUMMY
                LNUM = LNUM + 1
                CALL CHECK_HEADER( DUMMY, 'VERTICAL METHOD', 
     &                             FILENAM( I ), LNUM )
     
                READ( TDEV, 93000 ) DUMMY
                LNUM = LNUM + 1
                CALL CHECK_HEADER( DUMMY( 1:20 ), 'STD       ' //
     &              'ALL       ', FILENAM( I ), LNUM )
                CALL CHECK_HEADER( DUMMY( 31:46 ), ' 0.       10000.',
     &                             FILENAM( I ), LNUM )
                VERTTYPE = DUMMY( 21:30 )
    
                READ( TDEV, 93000 ) DUMMY
                LNUM = LNUM + 1
                CALL CHECK_HEADER( DUMMY, 'END', FILENAM( I ), LNUM )
                
                READ( TDEV, 93000 ) DUMMY
                LNUM = LNUM + 1
                CALL CHECK_HEADER( DUMMY, 'EMISSIONS VALUES', 
     &                             FILENAM( I ), LNUM )
                
                READ( TDEV, 93000 ) DUMMY
                LNUM = LNUM + 1
                CALL CHECK_HEADER( DUMMY, 'ALL       ALL' //
     &              '            0.000', FILENAM( I ), LNUM )

                SAVLNUM( I ) = LNUM

            END DO

C.............  Write time step information to ASCII file
            IF( OUTFLAG ) THEN
                WRITE( ODEV, 93000 ) 'TIME INTERVAL'
                WRITE( ODEV, 93050 ) IBD, IBT, IED, IET
                WRITE( ODEV, 93000 ) 'METHOD'
                WRITE( ODEV, 93000 ) 'STD       ALL       ' //
     &              'EMVALUES  0.        50000.'
                WRITE( ODEV, 93000 ) 'END'
                WRITE( ODEV, 93000 ) 'VERTICAL METHOD'
                WRITE( ODEV, 93000 ) 'STD       ALL       ' //
     &              VERTTYPE // ' 0.       10000.'
                WRITE( ODEV, 93000 ) 'END'
                WRITE( ODEV, 93000 ) 'EMISSIONS VALUES'
                WRITE( ODEV, 93000 ) 'ALL       ALL' //
     &              '            0.000'
            END IF

C.............  Write time interval record to binary file
            INTBIN (1) = IBD        ! beginning date
            REALBIN(1) = IBT / 100  ! beginning time
            INTBIN (2) = IED        ! ending date
            REALBIN(2) = IET / 100  ! ending time
            
            WRITE( BDEV ) INTBIN(1), REALBIN(1), INTBIN(2), REALBIN(2)

C.............  Write counter record to binary file
            INTBIN (1) = 1        ! segment number
            INTBIN (2) = NTOTSRCS ! number of point sources in segment
            
            WRITE( BDEV ) INTBIN(1), INTBIN(2)

C.............  Write point source location record to binary file
            WRITE( BDEV ) ( SRCCOLS( I ), SRCROWS( I ), SRCLAYS( I ),
     &                      SRCFLOW( I ), SRCHTS( I ), I = 1, NTOTSRCS )

C.............  Loop through emissions by species
            DO I = 1, NMSPC
                BASENUM = 0
                RPTBYFILE = 0.  ! array
                RPTASCII  = 0.
                
                DO J = 1, NFILES
                    TDEV = FILEDEV( J,1 )
                    LNUM = SAVLNUM( J )
                    
                    DO
                        READ( UNIT=TDEV, FMT=93075, IOSTAT=IOS ) 
     &                      NUM, VNAME, EMIS
                        LNUM = LNUM + 1

C.........................  Check if we've reached the END line or
C                           if we've started a new species
                        IF( IOS /= 0 .OR. VNAME /= SPCNAM( I ) ) THEN
                            BACKSPACE( TDEV )
                            LNUM = LNUM - 1
                            EXIT
                        END IF

C.........................  Set up output format for emissions
                        IF( EMIS > 999999999. ) THEN
                            FMT = '( I10, A10, E10.3 )'
                        ELSE IF( EMIS > 99999999. ) THEN
                            FMT = '( I10, A10, F10.0 )'
                        ELSE IF( EMIS > 9999999. ) THEN
                            FMT = '( I10, A10, F10.1 )'
                        ELSE IF( EMIS > 999999. ) THEN
                            FMT = '( I10, A10, F10.2 )'
                        ELSE
                            FMT = '( I10, A10, F10.3 )'
                        END IF

C.........................  Write emissions to ASCII file     
                        IF( OUTFLAG ) THEN
                            WRITE( ODEV, FMT ) BASENUM + NUM, VNAME, EMIS
                        END IF

C.........................  Store emissions to write to binary file
                        SRCEMIS( BASENUM+NUM ) = EMIS

C.........................  Store reporting information
                        RPTBYFILE( J ) = RPTBYFILE( J ) + EMIS
                        RPTASCII = RPTASCII + EMIS
                        
                    END DO  ! loop over lines in file
            
                    BASENUM = BASENUM + NSRCS( J )
                    
                    SAVLNUM( J ) = LNUM
                    
                END DO  ! loop over files

C.................  Write point source emissions record to binary file
                INTBIN(1) = 1
                DO J = 1, 10
                    SPECID( J ) = SPCNAM( I )( J:J )
                END DO
                
                WRITE( BDEV ) INTBIN(1), SPECID, SRCEMIS
                
C.................  Calculate emissions total for report
                RPTBINARY = 0.
                DO J = 1, NTOTSRCS
                    RPTBINARY = RPTBINARY + SRCEMIS( J )
                END DO

C.................  Write report information for this time step
                RPTFILESUM = 0.
                DO K = 1, NFILES
                    RPTFILESUM = RPTFILESUM + RPTBYFILE( K )
                END DO

                WRITE( RPTDEV,RPTFMT ) IBD, IBT, SPCNAM( I ), 
     &              ( RPTBYFILE( K ), K = 1, NFILES ), RPTFILESUM, 
     &              RPTASCII, RPTBINARY

C.................  Store overall report information
                DO K = 1, NFILES
                    RPTALLFILE( I,K ) = RPTALLFILE( I,K ) + RPTBYFILE( K )
                END DO
                
                RPTALLASC( I ) = RPTALLASC( I ) + RPTASCII
                RPTALLBIN( I ) = RPTALLBIN( I ) + RPTBINARY
                
                SRCEMIS = 0.    ! array

            END DO  ! loop over species

            IF( OUTFLAG ) THEN
                WRITE( ODEV, 93000 ) 'END'
                WRITE( ODEV, 93000 ) 'ENDTIME'
            END IF

            LDATE = JDATE
            CALL NEXTIME( JDATE, JTIME, 10000 )
        
        END DO ! loop over time steps
        
        IF( EFLAG ) THEN
            MESG = 'Problem reading hourly emissions'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

C.........  Write report totals for all time steps
        DO I = 1, NMSPC
            RPTFILESUM = 0.
            DO K = 1, NFILES
                RPTFILESUM = RPTFILESUM + RPTALLFILE( I,K )
            END DO
        
            WRITE( RPTDEV,RPTALLFMT ) SPCNAM( I ), 
     &          ( RPTALLFILE( I,K ), K = 1, NFILES ), RPTFILESUM,
     &          RPTALLASC( I ), RPTALLBIN( I )
        END DO

C.........  Check for any sources that were in PELV but not in ASCII elevated file
        IF( PINGFLAG ) THEN

            NMISSASCII = 0
            MSGPRINT = .TRUE.        
            DO I = 1, NPINGSRC
            
                IF( .NOT. PINGFND( I ) ) THEN
                    NMISSASCII = NMISSASCII + 1
                    
                    IF( MSGPRINT ) THEN
                        MESG = 'WARNING: The following PinG ' //
     &                         'sources are in the PELV files but ' //
     &                         'are not' // CRLF() // BLANK10 //
     &                         'in the corresponding ASCII elevated ' //
     &                         'files.'
                        CALL M3MESG( MESG )
                        MSGPRINT = .FALSE.
                    END IF
                    
                    CSRC = PELVSRC( I )
                    MESG = 'FIPS: ' // CSRC( 1:6 ) //
     &                     '    Facility: ' // CSRC( 7:16 ) //
     &                     '    Stack: ' // CSRC( 17:26 )
                    CALL M3MESG( MESG )
                END IF
            END DO
        
            MESG = 'Number of PinG sources -'
            CALL M3MSG2( MESG )
            
            WRITE( MESG,94010 ) BLANK5 // 'In input ASCII elevated ' //
     &          'files:                       ', NPINGASCII
            CALL M3MSG2( MESG )
            WRITE( MESG,94010 ) BLANK5 // 'In PELV files:          ' //
     &          '                             ', NPINGSRC
            CALL M3MSG2( MESG )
            WRITE( MESG,94010 ) BLANK5 // 'In input ASCII elevated ' //
     &          'files but not in PELV files: ', NMISSPELV
            CALL M3MSG2( MESG )
            WRITE( MESG,94010 ) BLANK5 // 'In PELV files but not in' //
     &          ' input ASCII elevated files: ', NMISSASCII
            CALL M3MSG2( MESG )
            WRITE( MESG,94010 ) BLANK5 // 'In output ASCII elevated' //
     &          ' file:                       ', NPINGOUT
            CALL M3MSG2( MESG )
        
        END IF

C.........  End program successfully        
        CALL M3EXIT( PROGNAME, 0, 0, ' ', 0 )

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

93065   FORMAT( I10, A10, F10.5, F10.5, 2A10, A10 )

93070   FORMAT( I10, A10, F10.0, F10.0, 2A10, A10 )

93075   FORMAT( I10, A10, F10.0 )

93080   FORMAT( F10.1, F10.2, F10.1, F10.0 )

93090   FORMAT( 3(I8,1X) )

93095   FORMAT( 3(I8,1X), A6, 1X, 2(A20,1X), F10.3 )

C...........   Internal buffering formats.............94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

94020   FORMAT( 8X, 'at time ', A8 )

C******************  INTERNAL SUBPROGRAMS  *****************************
        
        CONTAINS

C.............  This internal subroutine writes an error message
C               when a bad header is read.
            SUBROUTINE WRITE_HEADER_ERROR( FILENAM, LNUM )
            
C.............  Subroutine arguments
            CHARACTER(*), INTENT(IN) :: FILENAM   ! file name
            INTEGER,      INTENT(IN) :: LNUM      ! line number

C.............  Local subroutine variables            
            CHARACTER(512) MESG
            
C.............................................................................            
                
            WRITE( MESG,94010 ) 'ERROR: Header incorrect in file ' //
     &             CRLF() // BLANK10 // TRIM( FILENAM ) //
     &             CRLF() // BLANK10 // 'at line ', LNUM
            CALL M3MESG( MESG )
            
            MESG = 'Problem with file header information'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )            

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

93065   FORMAT( I10, A10, F10.5, F10.5, 2A10, A10 )

93070   FORMAT( I10, A10, F10.0, F10.0, 2A10, A10 )

93075   FORMAT( I10, A10, F10.0 )

93080   FORMAT( F10.1, F10.2, F10.1, F10.0 )

93090   FORMAT( 3(I8,1X) )

93095   FORMAT( 3(I8,1X), A6, 1X, 2(A20,1X), F10.3 )

C...........   Internal buffering formats.............94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

94020   FORMAT( 8X, 'at time ', A8 )
                
            END SUBROUTINE WRITE_HEADER_ERROR

C-----------------------------------------------------------------------------
C-----------------------------------------------------------------------------

C.............  This internal subroutine checks the value of a line
C               against the expected value.
            SUBROUTINE CHECK_HEADER( STRING, VALUE, FILENAM, LNUM )
            
C.............  Subroutine arguments
            CHARACTER(*), INTENT(IN) :: STRING  ! string that was read
            CHARACTER(*), INTENT(IN) :: VALUE   ! expected value
            CHARACTER(*), INTENT(IN) :: FILENAM ! file name
            INTEGER,      INTENT(IN) :: LNUM    ! line number
            
C.............  Local subroutine variables            
            CHARACTER(512) MESG
            
C.............................................................................            
            
            IF( TRIM( STRING ) /= TRIM( VALUE ) ) THEN
                WRITE( MESG,94010 ) 'ERROR: Header incorrect ' //
     &                 'in file' // 
     &                 CRLF() // BLANK10 // TRIM( FILENAM ) //
     &                 CRLF() // BLANK10 // 'at line ', LNUM
                CALL M3MESG( MESG )
                MESG = BLANK10 // 'Expected string "' // 
     &                 TRIM( VALUE ) // '" but read "' //
     &                 TRIM( STRING ) // '"'
                CALL M3MESG( MESG )
                
                MESG = 'Problem with file header information'
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF

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

93065   FORMAT( I10, A10, F10.5, F10.5, 2A10, A10 )

93070   FORMAT( I10, A10, F10.0, F10.0, 2A10, A10 )

93075   FORMAT( I10, A10, F10.0 )

93080   FORMAT( F10.1, F10.2, F10.1, F10.0 )

93090   FORMAT( 3(I8,1X) )

93095   FORMAT( 3(I8,1X), A6, 1X, 2(A20,1X), F10.3 )

C...........   Internal buffering formats.............94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

94020   FORMAT( 8X, 'at time ', A8 )
            
            END SUBROUTINE CHECK_HEADER

C-----------------------------------------------------------------------------
C-----------------------------------------------------------------------------

C.............  This internal subroutine writes an error message
C               when a header value is inconsistent with previous values.
            SUBROUTINE WRITE_HEADER_DIFF( FILENAM )
            
C.............  Subroutine arguments
            CHARACTER(*), INTENT(IN) :: FILENAM   ! file name

C.............  Local subroutine variables            
            CHARACTER(512) MESG
            
C.............................................................................            
                
            MESG = 'ERROR: Header value in file' //
     &             CRLF() // BLANK10 // TRIM( FILENAM ) //
     &             CRLF() // BLANK10 // 'does not match previous values'
            CALL M3MESG( MESG )
            
            MESG = 'Problem with file header information'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            
            END SUBROUTINE WRITE_HEADER_DIFF

C-----------------------------------------------------------------------------
C-----------------------------------------------------------------------------
        
C.............  This internal subroutine removes the first two digits of the
C               year from a 7-digit Julian date.
            INTEGER FUNCTION REMOVE_4DIGIT_YEAR( JDATE )
            
C.............  Subroutine arguments
            INTEGER, INTENT(IN) :: JDATE
            
C.............  Local subroutine variables
            INTEGER YRREMOVE

C.............................................................................            

            YRREMOVE = ( JDATE / 100000 ) * 100000   ! integer math
            REMOVE_4DIGIT_YEAR = JDATE - YRREMOVE
            
            RETURN
            
            END FUNCTION REMOVE_4DIGIT_YEAR

C-----------------------------------------------------------------------------
C-----------------------------------------------------------------------------

C.............  This internal subroutine writes the header information of
C               an ASCII elevated file
            SUBROUTINE WRITE_ASCII_HEADER

            WRITE( ODEV, 93010 ) 'CONTROL', GRDENV
            WRITE( ODEV, 93000 ) 'PTSOURCE'
            WRITE( ODEV, 93000 ) FNOTE
            WRITE( ODEV, 93015 ) NMSPC, 0, NTOTSRCS, 1, NPARAM
            WRITE( ODEV, 93015 ) PDEVOUT, 0, GSWITCH
            WRITE( ODEV, 93015 ) USWITCH, LSWITCH, 0, MSWITCH, VSWITCH
            WRITE( ODEV, 93015 ) 1, 0, ESWITCH
            WRITE( ODEV, 93015 ) DDEVOUT, RDEVOUT, 0, TDEVOUT, 
     &                           MDEVOUT, WDEVOUT
         
            DO I = 1, NMSPC
                WRITE( ODEV, 93020 ) SPCNAM( I )
            END DO
            
            WRITE( ODEV, 93030 ) IBD, IBT, IED, IET
            WRITE( ODEV, 93000 ) 'END'
            WRITE( ODEV, 93000 ) 'REGION'
            WRITE( ODEV, 93040 ) 0., 0., P_ALPHA
            WRITE( ODEV, 93040 ) XORIG, YORIG
            
            IF( LLGRID ) THEN
                WRITE( ODEV, 93045 ) XCELL, YCELL
            ELSE
                WRITE( ODEV, 93040 ) XCELL, YCELL
            END IF
            
            WRITE( ODEV, 93050 ) NCOLS, NROWS, NULAYS
            WRITE( ODEV, 93060 ) NZLOWR, NZUPPR, HTSUR, HTLOWR, HTUPPR
            WRITE( ODEV, 93000 ) 'END'
            
            WRITE( ODEV, 93000 ) 'POINT SOURCES'

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

93065   FORMAT( I10, A10, F10.5, F10.5, 2A10, A10 )

93070   FORMAT( I10, A10, F10.0, F10.0, 2A10, A10 )

93075   FORMAT( I10, A10, F10.0 )

93080   FORMAT( F10.1, F10.2, F10.1, F10.0 )

93090   FORMAT( 3(I8,1X) )

93095   FORMAT( 3(I8,1X), A6, 1X, 2(A20,1X), F10.3 )

C...........   Internal buffering formats.............94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

94020   FORMAT( 8X, 'at time ', A8 )
            
            END SUBROUTINE WRITE_ASCII_HEADER

C-----------------------------------------------------------------------------
C-----------------------------------------------------------------------------

C.............  This internal subroutine writes the header information of
C               a binary PTSOURCE file
            SUBROUTINE WRITE_BINARY_HEADER

C.............  Local subroutine variables
            CHARACTER(4)       FILENAME(10)   ! file name
            CHARACTER(4)       FILEID(60)     ! file identifier
            CHARACTER(4), ALLOCATABLE :: SPECIDS(:,:)  ! species names

C.............................................................................            
            
C.............  Build and write file description header record            
            FILENAME = ' '
            DO I = 1, LEN_TRIM( HDRKEY )
                FILENAME( I ) = HDRKEY( I:I )
            END DO
            
            FILEID   = ' '
            DO I = 1, LEN_TRIM( FNOTE )
                FILEID( I ) = FNOTE( I:I )
            END DO
            
            INTBIN (1) = 1          ! number of segments
            INTBIN (2) = NMSPC      ! number of chemical species
            INTBIN (3) = IBD        ! beginning date of file
            REALBIN(1) = IBT / 100  ! beginning time of file
            INTBIN (4) = IED        ! ending date of file
            REALBIN(2) = IET / 100  ! ending time of file
            
            WRITE( BDEV ) FILENAME, FILEID, INTBIN(1), INTBIN(2), 
     &                    INTBIN(3), REALBIN(1), INTBIN(4), REALBIN(2)

C.............  Build and write region description header record
            REALBIN(1) = 0.       ! x-coord of reference origin
            REALBIN(2) = 0.       ! y-coord of reference origin
            INTBIN (1) = P_ALPHA  ! UTM zone
            REALBIN(3) = XORIG    ! x-coord of modeling region origin
            REALBIN(4) = YORIG    ! y-coord of modeling region origin
            REALBIN(5) = XCELL    ! x-dir cell size
            REALBIN(6) = YCELL    ! y-dir cell size
            INTBIN (2) = NCOLS    ! number of cells in x-dir
            INTBIN (3) = NROWS    ! number of cells in y-dir
            INTBIN (4) = NULAYS   ! number of cells in z-dir
            INTBIN (5) = NZLOWR   ! number of cells between surface and diff break
            INTBIN (6) = NZUPPR   ! number of cells between diff break and top
            REALBIN(7) = HTSUR    ! height of surface layer
            REALBIN(8) = HTLOWR   ! height of cells between surface and diff break
            REALBIN(9) = HTUPPR   ! height of cells between diff break and top
            
            WRITE( BDEV ) REALBIN(1), REALBIN(2), INTBIN(1), REALBIN(3),
     &                    REALBIN(4), REALBIN(5), REALBIN(6), INTBIN(2),
     &                    INTBIN(3), INTBIN(4), INTBIN(5), INTBIN(6),
     &                    REALBIN(7), REALBIN(8), REALBIN(9)

C.............  Build and write segment description header record
            INTBIN (1) = 0        ! x-location of segment origin
            INTBIN (2) = 0        ! y-location of segment origin
            INTBIN (3) = NCOLS    ! number of cells in x-dir in segment
            INTBIN (4) = NROWS    ! number of cells in y-dir in segment
            
            WRITE( BDEV ) INTBIN(1), INTBIN(2), INTBIN(3), INTBIN(4)

C.............  Build and write species description header records
            ALLOCATE( SPECIDS( 10, NMSPC ), STAT=IOS )
            CALL CHECKMEM( IOS, 'SPECIDS', PROGNAME )
            SPECIDS = ' '

            DO I = 1, NMSPC
                DO J = 1, 10
                    SPECIDS( J,I ) = SPCNAM( I )( J:J )
                END DO
            END DO
                
            WRITE( BDEV ) SPECIDS
            
            DEALLOCATE( SPECIDS )

C.............  Build and write time-invariant counter record
            INTBIN (1) = 1        ! segment number
            INTBIN (2) = NTOTSRCS ! number of point sources in segment
            
            WRITE( BDEV ) INTBIN(1), INTBIN(2)
            
            END SUBROUTINE WRITE_BINARY_HEADER

C-----------------------------------------------------------------------------
C-----------------------------------------------------------------------------

C.............  This internal subroutine writes the header of the report file
            SUBROUTINE WRITE_REPORT_HEADER

C.............  Local subroutine variables
            INTEGER LEN             ! current length of buffer

            CHARACTER(300) BUFFER

C.............................................................................            

            BUFFER = 'Date ; Time; Species   ; '
            LEN = 25
            
            DO I = 1, NFILES
                WRITE( BUFFER, '(A,I2,A)' ) 
     &              BUFFER( 1:LEN ) // 'File ', I, '     ; '
                LEN = LEN + 14
            END DO
            
            BUFFER = BUFFER( 1:LEN ) // 'Sum         ; ' //
     &               'ASCII       ; Binary      ;'

            WRITE( RPTDEV,'(A)' ) TRIM( BUFFER )

            END SUBROUTINE WRITE_REPORT_HEADER

        END PROGRAM MRGELEV
