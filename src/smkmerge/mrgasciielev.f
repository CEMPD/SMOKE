        PROGRAM MRGELEV

C***********************************************************************
C  program body starts at line
C
C  DESCRIPTION:
C       This program combines ASCII elevated files produced by 
C       Smkmerge.
C
C  PRECONDITIONS REQUIRED:
C       ASCII elevated files created by Smkmerge
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
        INTEGER       GETFLINE
        INTEGER       JUNIT
        INTEGER       PROMPTFFILE

        EXTERNAL      CRLF, GETFLINE, JUNIT, PROMPTFFILE

C.........  LOCAL PARAMETERS
        CHARACTER(50), PARAMETER :: CVSW = '$Name$'  ! CVS release tag

C.........  LOCAL VARIABLES

C.........  Allocatable arrays
        INTEGER, ALLOCATABLE :: FILEDEV( : )        ! input file device numbers
        INTEGER, ALLOCATABLE :: NSRCS( : )          ! num. srcs per file
        INTEGER, ALLOCATABLE :: SDATE( : )          ! start date by file
        INTEGER, ALLOCATABLE :: STIME( : )          ! start time by file
        INTEGER, ALLOCATABLE :: EDATE( : )          ! end date by file
        INTEGER, ALLOCATABLE :: ETIME( : )          ! end time by file
        
        CHARACTER(256), ALLOCATABLE :: FILENAM( : ) ! input file names
        CHARACTER(10),  ALLOCATABLE :: SPCNAM( : )  ! model species names

C.........  File units and logical names
        INTEGER          LDEV       ! unit for log file
        INTEGER          IDEV       ! unit for logical names list for elevated files
        INTEGER          ODEV       ! unit for output file
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
        INTEGER          FIP                  ! FIPS code
        CHARACTER(10)    FCID                 ! facility ID
        CHARACTER(10)    SKID                 ! stack ID

C.........  Hourly emissions variables
        INTEGER          NUM                  ! source number
        CHARACTER(10)    VNAME                ! species name
        REAL             EMIS                 ! hourly emissions
        
C.........  Other local variables
        INTEGER          I, J, K, L           ! indexes and counters
        INTEGER          IDUM, IDUM2          ! dummy integers
        INTEGER          IOS                  ! i/o status
        INTEGER          MXFILES              ! maximum number of input files
        INTEGER          NFILES               ! number of input files
        INTEGER          NTOTSRCS             ! total number of output sources
        INTEGER          IBD, IBT             ! start date and time
        INTEGER          IED, IET             ! end date and time
        INTEGER          BASENUM              ! base source number for hourly emissions

        REAL             FDUM, FDUM2          ! dummy reals

        LOGICAL ::       EFLAG = .FALSE.      ! true: an error happened

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

C.........  Open list of input files
        MESG = 'Enter logical name for ASCII ELEVATED inputs list'
        IDEV = PROMPTFFILE( MESG, .TRUE., .TRUE., 'FILELIST', PROGNAME )

C.........  Open output file
        MESG = 'Enter logical name for output ' //
     &         'ASCII ELEVATED SOURCES file'
        ODEV = PROMPTFFILE( MESG, .FALSE., .TRUE., 'OUTFILE', PROGNAME )

C.........  Determine maximum number of input files in list
        MXFILES = GETFLINE( IDEV, 'List of files to merge' )
        
C.........  Allocate memory to store file device units
        ALLOCATE( FILEDEV( MXFILES ), STAT=IOS )
        CALL CHECKMEM( IOS, 'FILEDEV', PROGNAME )
        ALLOCATE( FILENAM( MXFILES ), STAT=IOS )
        CALL CHECKMEM( IOS, 'FILENAM', PROGNAME )
        
        J = 0
        DO I = 1, MXFILES
        
            READ( IDEV, 93000, IOSTAT=IOS ) LINE
            
            IF( IOS /= 0 ) THEN
                EFLAG = .TRUE.
                WRITE( MESG, 94010 ) 'I/O error', IOS,
     &              'reading file list at line', I
                CALL M3MESG( MESG )
                CYCLE
            END IF
            
C.............  Skip any blank lines
            IF( LINE == ' ' ) CYCLE
            
C.............  Open input file
            TDEV = JUNIT()
            OPEN( UNIT=TDEV, FILE=TRIM( LINE ), STATUS='OLD',
     &            ACTION='READ', IOSTAT=IOS )
            
            IF( IOS /= 0 ) THEN
                EFLAG = .TRUE.
                MESG = 'Could not open file ' // CRLF() // BLANK5 // 
     &                 TRIM( LINE )
                CALL M3MESG( MESG )
                CYCLE
            ELSE
                J = J + 1
                FILEDEV( J ) = TDEV
                FILENAM( J ) = TRIM( LINE )
            END IF
        
        END DO

        IF( EFLAG ) THEN
            MESG = 'Problem opening input files'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

        NFILES = J

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
        GRDENV = ''
        DO I = 1, NFILES
        
            TDEV = FILEDEV( I )
            
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
        
            TDEV = FILEDEV( I )
        
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
        
            TDEV = FILEDEV( I )

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
        
            TDEV = FILEDEV( I )
        
            DO J = 1, NSRCS( I )
                READ( TDEV, 93065 ) IDUM, DUMMY, XCOORD, YCOORD, ! FORMAT PROBLEM
     &              FCID, SKID, FIP

                IF( TRIM( DUMMY ) /= 'STD' ) THEN
                    MESG = 'Problem reading source information'
                    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                END IF
                
                K = K + 1
                WRITE( ODEV, 93065 ) K, 'STD       ', XCOORD, ! FORMAT PROBLEM
     &              YCOORD, FCID, SKID, FIP

C................. Read stack parameters; this will change when we set PinG sources                    
                READ( TDEV, 93000 ) LINE
                WRITE( ODEV, 93000 ) TRIM( LINE )
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
            
                TDEV = FILEDEV( I )
                
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
                    TDEV = FILEDEV( J )
                    
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

93065   FORMAT( I10, A10, F10.5, F10.5, 2A10, I10.5 )

93070   FORMAT( I10, A10, F10.3 )

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