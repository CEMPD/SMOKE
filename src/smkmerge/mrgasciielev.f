        PROGRAM MRGELEV

C***********************************************************************
C  program body starts at line
C
C  DESCRIPTION:
C       This program combines ASCII elevated files produced by 
C       Smkmerge. It can optionally use a corresponding list of elevated
C       point source files to flag combined PinG sources.
C
C  PRECONDITIONS REQUIRED:
C       ASCII elevated files created by Smkmerge
C       Point source elevated files created by Elevpoint
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C       
C
C  REVISION  HISTORY:
C       Created 4/2005 by C. Seppanen
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

        IMPLICIT NONE
        
C.........  INCLUDES
        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
        
C.........  EXTERNAL FUNCTIONS
        CHARACTER(2)  CRLF
        INTEGER       FINDC
        INTEGER       GETFLINE
        INTEGER       JUNIT
        INTEGER       PROMPTFFILE

        EXTERNAL      CRLF, FINDC, GETFLINE, JUNIT, PROMPTFFILE

C.........  LOCAL PARAMETERS
        CHARACTER(50), PARAMETER :: CVSW = '$Name$'  ! CVS release tag
        INTEGER,       PARAMETER :: MXPING = 300        ! final number of PinG sources

C.........  LOCAL VARIABLES

C.........  Allocatable arrays
        INTEGER, ALLOCATABLE :: FILEDEV( :,: )      ! input file device numbers
        INTEGER, ALLOCATABLE :: FILELNS( : )        ! num. lines per PELV file
        INTEGER, ALLOCATABLE :: NSRCS( : )          ! num. srcs per file
        INTEGER, ALLOCATABLE :: SDATE( : )          ! start date by file
        INTEGER, ALLOCATABLE :: STIME( : )          ! start time by file
        INTEGER, ALLOCATABLE :: EDATE( : )          ! end date by file
        INTEGER, ALLOCATABLE :: ETIME( : )          ! end time by file
        INTEGER, ALLOCATABLE :: SRCIDX( : )         ! index for sorting sources
        INTEGER, ALLOCATABLE :: EMISIDX( : )        ! index for sorting emissions
        
        CHARACTER(256), ALLOCATABLE :: FILENAM( :,: ) ! input file names
        CHARACTER(10),  ALLOCATABLE :: SPCNAM( : )  ! model species names
        CHARACTER(FPLLEN3+CHRLEN3), ALLOCATABLE :: PELVSRCA(:) ! unsorted source list
        CHARACTER(FPLLEN3+CHRLEN3), ALLOCATABLE :: PELVSRC(:)  ! sorted source list
        
        REAL, ALLOCATABLE :: PELVEMISA( : )      ! unsorted source emissions
        REAL, ALLOCATABLE :: PELVEMIS( : )       ! sorted source emissions
        
        LOGICAL, ALLOCATABLE :: PINGFND( : )     ! true: PinG source was found
        LOGICAL, ALLOCATABLE :: ISPING( : )      ! true: source remains PinG

C.........  File units and logical names
        INTEGER          LDEV       ! unit for log file
        INTEGER          IDEV       ! unit for list of elevated files
        INTEGER          PDEV       ! unit for list of PELV files 
        INTEGER          ODEV       ! unit for output file
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
        
C.........  Other local variables
        INTEGER          I, J, K, L           ! indexes and counters
        INTEGER          IDUM, IDUM2          ! dummy integers
        INTEGER          IOS                  ! i/o status
        INTEGER          MXFILES              ! maximum number of input files
        INTEGER          MXPFILES             ! maximum number of PELV files
        INTEGER          MXLINES              ! maximum number of lines
        INTEGER          NFILES               ! number of input files
        INTEGER          NPINGSRC             ! number of PinG sources
        INTEGER          NTOTSRCS             ! total number of output sources
        INTEGER          IBD, IBT             ! start date and time
        INTEGER          IED, IET             ! end date and time
        INTEGER          BASENUM              ! base source number for hourly emissions
        INTEGER          ENUM, PNUM, GNUM     ! elevated, PinG, and group numbers

        REAL             FDUM, FDUM2          ! dummy reals

        LOGICAL ::       EFLAG = .FALSE.      ! true: an error happened
        LOGICAL          FOUND                ! true: matching source was found
        LOGICAL          PINGFLAG             ! true: process PinG sources
        CHARACTER(256)   DUMMY                ! dummy character string
        CHARACTER(256)   LINE                 ! input line
        CHARACTER(300)   MESG                 ! message buffer

        CHARACTER(16) :: PROGNAME = 'MRGELEV' ! program name

C***********************************************************************
C   begin body of program MRGELEV

        LDEV = INIT3()

C.........  Write out copywrite, version, web address, header info, and 
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

C.........  Open output file
        MESG = 'Enter logical name for output ' //
     &         'ASCII ELEVATED SOURCES file'
        ODEV = PROMPTFFILE( MESG, .FALSE., .TRUE., 'OUTFILE', PROGNAME )

C.........  Determine maximum number of input files in lists
        MXFILES  = GETFLINE( IDEV, 'List of files to merge' )
        
        IF( PINGFLAG ) THEN
            MXPFILES = GETFLINE( PDEV, 'List of PELV input files' )
        
            IF( MXFILES /= MXPFILES ) THEN
                MESG = 'ERROR: Number of ASCII elevated files does ' //
     &            'not match number of elevated point source files'
                CALL M3MSG2( MESG )
                
                MESG = 'Problem with input files'
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF
        END IF
        
C.........  Allocate memory to store file device units
        ALLOCATE( FILEDEV( MXFILES, 2 ), STAT=IOS )
        CALL CHECKMEM( IOS, 'FILEDEV', PROGNAME )
        ALLOCATE( FILENAM( MXFILES, 2 ), STAT=IOS )
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
                    WRITE( MESG, 94010 ) 'I/O error', IOS,
     &                  'reading list of input files at line', I
                    CALL M3MESG( MESG )
                    CYCLE
                END IF
            
C.................  Skip any blank lines
                IF( LINE == ' ' ) CYCLE

C.................  Open input file
                TDEV = JUNIT()
                OPEN( UNIT=TDEV, FILE=TRIM( LINE ), STATUS='OLD',
     &                ACTION='READ', IOSTAT=IOS )
            
                IF( IOS /= 0 ) THEN
                    EFLAG = .TRUE.
                    MESG = 'Could not open file ' // CRLF() // BLANK5 // 
     &                     TRIM( LINE )
                    CALL M3MESG( MESG )
                    CYCLE
                ELSE
                    K = K + 1
                    FILEDEV( K,J ) = TDEV
                    FILENAM( K,J ) = TRIM( LINE )
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
            
            DO I = 1, NPINGSRC
                WRITE( *,* ) PELVEMIS( EMISIDX( I ) )
            END DO
            
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
        ALLOCATE( EDATE( NFILES ), STAT=IOS )
        CALL CHECKMEM( IOS, 'EDATE', PROGNAME )
        ALLOCATE( ETIME( NFILES ), STAT=IOS )
        CALL CHECKMEM( IOS, 'ETIME', PROGNAME )

C.........  Check the header of each input file
        MESG = 'Checking headers of input files...'
        CALL M3MSG2( MESG )

        GRDENV = ''
        DO I = 1, NFILES
        
            TDEV = FILEDEV( I,1 )
            
            READ( TDEV, 93010 ) DUMMY, TMPGRDENV
            IF( TRIM( DUMMY ) /= 'CONTROL' ) THEN
                EFLAG = .TRUE.
                CALL WRITE_HEADER_ERROR( FILENAM( I ) )
                CYCLE
            END IF
            
            READ( TDEV, 93000 ) DUMMY
            IF( TRIM( DUMMY ) /= 'PTSOURCE' ) THEN
                EFLAG = .TRUE.
                CALL WRITE_HEADER_ERROR( FILENAM( I ) )
                CYCLE
            END IF
            
            READ( TDEV, 93000 ) TMPNOTE
            
            READ( TDEV, 93015 ) TMPNMSPC, IDUM, TMPNOUT, 
     &                          IDUM2, TMPNPARAM
            IF( IDUM /= 0 .OR. IDUM2 /= 1 ) THEN
                EFLAG = .TRUE.
                CALL WRITE_HEADER_ERROR( FILENAM( I ) )
                CYCLE
            END IF
            
            READ( TDEV, 93015 ) TMPPDEVOUT, IDUM, TMPGSWITCH
            IF( IDUM /= 0 ) THEN
                EFLAG = .TRUE.
                CALL WRITE_HEADER_ERROR( FILENAM( I ) )
                CYCLE
            END IF
            
            READ( TDEV, 93015 ) TMPUSWITCH, TMPLSWITCH, IDUM, 
     &                          TMPMSWITCH, TMPVSWITCH
            IF( IDUM /= 0 ) THEN
                EFLAG = .TRUE.
                CALL WRITE_HEADER_ERROR( FILENAM( I ) )
                CYCLE
            END IF
            
            READ( TDEV, 93015 ) IDUM, IDUM2, TMPESWITCH
            IF( IDUM /= 1 .OR. IDUM2 /= 0 ) THEN
                EFLAG = .TRUE.
                CALL WRITE_HEADER_ERROR( FILENAM( I ) )
                CYCLE
            END IF
            
            READ( TDEV, 93015 ) TMPDDEVOUT, TMPRDEVOUT, IDUM, 
     &                          TMPTDEVOUT, TMPMDEVOUT, TMPWDEVOUT
            IF( IDUM /= 0 ) THEN
                EFLAG = .TRUE.
                CALL WRITE_HEADER_ERROR( FILENAM( I ) )
                CYCLE
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
                    EFLAG = .TRUE.
                    CALL WRITE_HEADER_DIFF( FILENAM( I ) )
                    CYCLE
                END IF
            END IF
            
            NSRCS( I ) = TMPNOUT

        END DO
        
        IF( EFLAG ) THEN
            MESG = 'Problem with file header information'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF        

C.........  Read and check species names for each file
        ALLOCATE( SPCNAM( NMSPC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'SPCNAM', PROGNAME )

        DO I = 1, NFILES
        
            TDEV = FILEDEV( I,1 )
        
            DO J = 1, NMSPC
                READ( TDEV, 93020 ) TMPNAM
                
                IF( I == 1 ) THEN
                    SPCNAM( J ) = TMPNAM
                ELSE
                    IF( TMPNAM /= SPCNAM( J ) ) THEN
                        EFLAG = .TRUE.
                        CALL WRITE_HEADER_DIFF( FILENAM( I ) )
                        EXIT
                    END IF
                END IF
            END DO
        END DO

        IF( EFLAG ) THEN
            MESG = 'Problem with file model species names'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

C.........  Continue checking file headers
        P_ALPHA = 0
        DO I = 1, NFILES
        
            TDEV = FILEDEV( I,1 )

C.............  Read file start and end dates and times
            READ( TDEV, 93030 ) SDATE( I ), STIME( I ),
     &                          EDATE( I ), ETIME( I )
        
            READ( TDEV, 93000 ) DUMMY
            IF( TRIM( DUMMY ) /= 'END' ) THEN
                EFLAG = .TRUE.
                CALL WRITE_HEADER_ERROR( FILENAM( I ) )
                CYCLE
            END IF

C.............  Read grid and layer information            
            READ( TDEV, 93000 ) DUMMY
            IF( TRIM( DUMMY ) /= 'REGION' ) THEN
                EFLAG = .TRUE.
                CALL WRITE_HEADER_ERROR( FILENAM( I ) )
                CYCLE
            END IF
            
            READ( TDEV, 93040 ) FDUM, FDUM2, TMPALPHA
            IF( FDUM /= 0. .OR. FDUM2 /= 0. ) THEN
                EFLAG = .TRUE.
                CALL WRITE_HEADER_ERROR( FILENAM( I ) )
                CYCLE
            END IF
            
            READ( TDEV, 93040 ) TMPXORIG, TMPYORIG
            
            READ( TDEV, 93045 ) TMPXCELL, TMPYCELL  ! FORMAT PROBLEM
            
            READ( TDEV, 93050 ) TMPNCOLS, TMPNROWS, TMPNULAYS
            READ( TDEV, 93060 ) TMPNZLOWR, TMPNZUPPR, TMPHTSUR, 
     &                          TMPHTLOWR, TMPHTUPPR
            
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
                    EFLAG = .TRUE.
                    CALL WRITE_HEADER_DIFF( FILENAM( I ) )
                    CYCLE
                END IF
            END IF

            READ( TDEV, 93000 ) DUMMY
            IF( TRIM( DUMMY ) /= 'END' ) THEN
                EFLAG = .TRUE.
                CALL WRITE_HEADER_ERROR( FILENAM( I ) )
                CYCLE
            END IF
            
            READ( TDEV, 93000 ) DUMMY
            IF( TRIM( DUMMY ) /= 'POINT SOURCES' ) THEN
                EFLAG = .TRUE.
                CALL WRITE_HEADER_ERROR( FILENAM( I ) )
                CYCLE
            END IF

        END DO

        IF( EFLAG ) THEN
            MESG = 'Problem with file region information'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

C.........  Calculate total number of output sources
        NTOTSRCS = 0
        DO I = 1, NFILES
            NTOTSRCS = NTOTSRCS + NSRCS( I )
        END DO

C.........  Determine start and end date/time
        IBD = SDATE( 1 )
        IBT = STIME( 1 )
        IED = EDATE( 1 )
        IET = ETIME( 1 )

C.........  Write header to output file
        WRITE( ODEV, 93010 ) 'CONTROL', GRDENV
        WRITE( ODEV, 93000 ) 'PTSOURCE'
        WRITE( ODEV, 93000 ) FNOTE
        WRITE( ODEV, 93015 ) NMSPC, 0, NTOTSRCS, 1, NPARAM
        WRITE( ODEV, 93015 ) PDEVOUT, 0, GSWITCH
        WRITE( ODEV, 93015 ) USWITCH, LSWITCH, 0, MSWITCH, VSWITCH
        WRITE( ODEV, 93015 ) 1, 0, ESWITCH
        WRITE( ODEV, 93015 ) DDEVOUT, RDEVOUT, 0, TDEVOUT, 
     &                       MDEVOUT, WDEVOUT
     
        DO I = 1, NMSPC
            WRITE( ODEV, 93020 ) SPCNAM( I )
        END DO
        
        WRITE( ODEV, 93030 ) IBD, IBT, IED, IET
        WRITE( ODEV, 93000 ) 'END'
        WRITE( ODEV, 93000 ) 'REGION'
        WRITE( ODEV, 93040 ) 0., 0., P_ALPHA
        WRITE( ODEV, 93040 ) XORIG, YORIG
        WRITE( ODEV, 93045 ) XCELL, YCELL   ! FORMAT PROBLEM
        WRITE( ODEV, 93050 ) NCOLS, NROWS, NULAYS
        WRITE( ODEV, 93060 ) NZLOWR, NZUPPR, HTSUR, HTLOWR, HTUPPR
        WRITE( ODEV, 93000 ) 'END'
        
        WRITE( ODEV, 93000 ) 'POINT SOURCES'

C.........  Read and output source information for each file
        K = 0
        DO I = 1, NFILES
        
            TDEV = FILEDEV( I,1 )
        
            DO J = 1, NSRCS( I )
                READ( TDEV, 93065 ) IDUM, DUMMY, XCOORD, YCOORD, ! FORMAT PROBLEM
     &              FCID, SKID, TFIP

                IF( TRIM( DUMMY ) /= 'STD' ) THEN
                    MESG = 'Problem reading source information'
                    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                END IF

C.................  Read stack parameters                    
                READ( TDEV, 93080 ) STKHT, STKDM, STKTK, STKVE
                
C.................  Check for source in PinG list
                IF( STKDM < 0 .AND. PINGFLAG ) THEN
                    TFIP = ADJUSTR( TFIP )
                
                    CSRC = TFIP( 5:10 ) // FCID // SKID
                
                    L = FINDC( CSRC, NPINGSRC, PELVSRC )
                    
                    IF( L < 1 ) THEN
                        MESG = 'WARNING: PinG source ' //
     &                      TRIM( CSRC ) // ' is in an ASCII ' //
     &                      'elevated file' // CRLF() // BLANK10 // 
     &                      'but not in the corresponding PELV ' //
     &                      'file; retaining PinG status for source'
                        CALL M3MESG( MESG )
                    ELSE
                        PINGFND( L ) = .TRUE.

C.........................  Check if source remains as PinG source
                        IF( .NOT. ISPING( L ) ) THEN
                            STKDM = -STKDM
                        END IF
                        WRITE( *,* ) ISPING( L ), PELVEMIS( L ), STKDM
                    END IF
                END IF
                
                K = K + 1
                WRITE( ODEV, 93065 ) K, 'STD       ', XCOORD, ! FORMAT PROBLEM
     &              YCOORD, FCID, SKID, TFIP                
                WRITE( ODEV, 93080 ) STKHT, STKDM, STKTK, STKVE
            END DO
            
            READ( TDEV, 93000 ) DUMMY
            IF( TRIM( DUMMY ) /= 'END' ) THEN
                MESG = 'Problem reading source information'
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF
        END DO

C.........  Check that we read the right number of sources
        IF( K /= NTOTSRCS ) THEN
            MESG = 'Wrong number of sources'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2)
        END IF
        
        WRITE( ODEV, 93000 ) 'END'

C.........  Read and output hourly emissions for each file
        DO L = 1, 25 ! loop over output hours somehow
        
            DO I = 1, NFILES
            
                TDEV = FILEDEV( I,1 )
                
                IF( L /= 1 ) THEN
                    READ( TDEV, 93000 ) DUMMY
                    IF( TRIM( DUMMY ) /= 'END' ) THEN
                        EFLAG = .TRUE.
                        CALL WRITE_HEADER_ERROR( FILENAM( I ) )
                        CYCLE
                    END IF
                    
                    READ( TDEV, 93000 ) DUMMY
                    IF( TRIM( DUMMY ) /= 'ENDTIME' ) THEN
                        EFLAG = .TRUE.
                        CALL WRITE_HEADER_ERROR( FILENAM( I ) )
                        CYCLE
                    END IF
                END IF
                
                READ( TDEV, 93000 ) DUMMY
                IF( TRIM( DUMMY ) /= 'TIME INTERVAL' ) THEN
                    EFLAG = .TRUE.
                    CALL WRITE_HEADER_ERROR( FILENAM( I ) )
                    CYCLE
                END IF
                
                READ( TDEV, 93050 ) IBD, IBT, IED, IET
                
                READ( TDEV, 93000 ) DUMMY
                IF( TRIM( DUMMY ) /= 'METHOD' ) THEN
                    EFLAG = .TRUE.
                    CALL WRITE_HEADER_ERROR( FILENAM( I ) )
                    CYCLE
                END IF
                
                READ( TDEV, 93000 ) DUMMY
                IF( TRIM( DUMMY ) /= 'STD       ALL       ' //
     &              'EMVALUES  0.        50000.' ) THEN
                    EFLAG = .TRUE.
                    CALL WRITE_HEADER_ERROR( FILENAM( I ) )
                    CYCLE
                END IF
                
                READ( TDEV, 93000 ) DUMMY
                IF( TRIM( DUMMY ) /= 'END' ) THEN
                    EFLAG = .TRUE.
                    CALL WRITE_HEADER_ERROR( FILENAM( I ) )
                    CYCLE
                END IF
                
                READ( TDEV, 93000 ) DUMMY
                IF( TRIM( DUMMY ) /= 'VERTICAL METHOD' ) THEN
                    EFLAG = .TRUE.
                    CALL WRITE_HEADER_ERROR( FILENAM( I ) )
                    CYCLE
                END IF
           
                READ( TDEV, 93000 ) DUMMY
                IF( DUMMY( 1:20 ) /= 'STD       ALL       ' ) THEN
                    EFLAG = .TRUE.
                    CALL WRITE_HEADER_ERROR( FILENAM( I ) )
                    CYCLE
                END IF
                IF( DUMMY( 31:46 ) /= ' 0.       10000.' ) THEN
                    EFLAG = .TRUE.
                    CALL WRITE_HEADER_ERROR( FILENAM( I ) )
                    CYCLE
                END IF
                VERTTYPE = DUMMY( 21:30 )
    
                READ( TDEV, 93000 ) DUMMY
                IF( TRIM( DUMMY ) /= 'END' ) THEN
                    EFLAG = .TRUE.
                    CALL WRITE_HEADER_ERROR( FILENAM( I ) )
                    CYCLE
                END IF
                
                READ( TDEV, 93000 ) DUMMY
                IF( TRIM( DUMMY ) /= 'EMISSIONS VALUES' ) THEN
                    EFLAG = .TRUE.
                    CALL WRITE_HEADER_ERROR( FILENAM( I ) )
                    CYCLE
                END IF
                
                READ( TDEV, 93000 ) DUMMY
                IF( TRIM( DUMMY ) /= 
     &              'ALL       ALL            0.000' ) THEN
                    EFLAG = .TRUE.
                    CALL WRITE_HEADER_ERROR( FILENAM( I ) )
                    CYCLE
                END IF
            
                IF( EFLAG ) THEN
                    MESG = 'Problem with emissions header'
                    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                END IF
            
                IF( I == 1 ) THEN
                    WRITE( ODEV, 93000 ) 'TIME INTERVAL'
                    WRITE( ODEV, 93050 ) IBD, IBT, IED, IET
                    WRITE( ODEV, 93000 ) 'METHOD'
                    WRITE( ODEV, 93000 ) 'STD       ALL       ' //
     &                  'EMVALUES  0.        50000.'
                    WRITE( ODEV, 93000 ) 'END'
                    WRITE( ODEV, 93000 ) 'VERTICAL METHOD'
                    WRITE( ODEV, 93000 ) 'STD       ALL       ' //
     &                  VERTTYPE // ' 0.       10000.'
                    WRITE( ODEV, 93000 ) 'END'
                    WRITE( ODEV, 93000 ) 'EMISSIONS VALUES'
                    WRITE( ODEV, 93000 ) 'ALL       ALL' //
     &                  '            0.000'
                END IF

            END DO

            DO I = 1, NMSPC
                BASENUM = 0
                
                DO J = 1, NFILES
                    TDEV = FILEDEV( J,1 )
                    
                    DO
                        READ( UNIT=TDEV, FMT=93070, IOSTAT=IOS ) 
     &                      NUM, VNAME, EMIS ! FORMAT PROBLEM

C.........................  Check if we've reached the END line or
C                           if we've started a new species
                        IF( IOS /= 0 .OR. VNAME /= SPCNAM( I ) ) THEN
                            BACKSPACE( TDEV )
                            EXIT
                        END IF
     
                        WRITE( ODEV, 93070 ) BASENUM + NUM, VNAME, EMIS ! FORMAT
                    END DO
            
                    BASENUM = BASENUM + NSRCS( J )
                END DO
            END DO

            WRITE( ODEV, 93000 ) 'END'
            WRITE( ODEV, 93000 ) 'ENDTIME'
        
        END DO ! loop over time steps
        
        IF( EFLAG ) THEN
            MESG = 'Problem reading hourly emissions'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

C.........  Check for any sources that were in PELV but not in ASCII elevated file
        IF( PINGFLAG ) THEN
        
            DO I = 1, NPINGSRC
            
                IF( .NOT. PINGFND( I ) ) THEN
                    CSRC = PELVSRC( I )
                    MESG = 'WARNING: PinG source ' // TRIM( CSRC ) //
     &                ' was in a PELV file' // CRLF() // BLANK10 // 
     &                'but not in the corresponding ASCII elevated file'
                    CALL M3MESG( MESG )
                END IF
            END DO
        
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

93070   FORMAT( I10, A10, F10.3 )

93080   FORMAT( F10.1, F10.2, F10.1, F10.0 )

93090   FORMAT( 3(I8,1X) )

93095   FORMAT( 3(I8,1X), A6, 1X, 2(A15,1X), F10.3 )

C...........   Internal buffering formats.............94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

C******************  INTERNAL SUBPROGRAMS  *****************************
        
        CONTAINS

C.............  This internal subroutine writes an error message
C               when a bad header is read.
            SUBROUTINE WRITE_HEADER_ERROR( FILENAM )
            
C.............  Subroutine arguments
            CHARACTER(*), INTENT(IN) :: FILENAM   ! file name

C.............  Local subroutine variables            
            CHARACTER(300) MESG
            
C.............................................................................            
                
            MESG = 'ERROR: header incorrect in file ' //
     &             CRLF() // BLANK5 // TRIM( FILENAM )
            CALL M3MESG( MESG )
                
            END SUBROUTINE WRITE_HEADER_ERROR

C.............  This internal subroutine writes an error message
C               when a header value is inconsistent with previous values.
            SUBROUTINE WRITE_HEADER_DIFF( FILENAM )
            
C.............  Subroutine arguments
            CHARACTER(*), INTENT(IN) :: FILENAM   ! file name

C.............  Local subroutine variables            
            CHARACTER(300) MESG
            
C.............................................................................            
                
            MESG = 'ERROR: header value in file ' //
     &             CRLF() // BLANK5 // TRIM( FILENAM ) //
     &             ' does not match previous values'
            CALL M3MESG( MESG )
                
            END SUBROUTINE WRITE_HEADER_DIFF
        
        END PROGRAM MRGELEV