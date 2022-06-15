
        PROGRAM MRGGRID

C***********************************************************************
C  program body starts at line 
C
C  DESCRIPTION:
C    Program MRGGRID reads 2-D and 3-D I/O API files and merges them
C    into a single 2-D or 3-D file (depending on the inputs)
C    The time period merged is adjusted based on the latest
C    starting file and earliest ending file, unless MRG_DIFF_DAY is
C    set in which case the time period is based on the standard 
C    environment variables. All variables are merged, even if different 
C    variables are in each file.
C
C  PRECONDITIONS REQUIRED:  
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C    Original by M. Houyoux 4/98
C    Modified by B.H. Baek  4/08
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
 
C...........   MODULES for public variables
C.........  This module contains the global variables for the 3-d grid
        USE MODGRID, ONLY: NGRID, NCOLS, NROWS, NLAYS, 
     &                     VGLVS, VGTYP, VGTOP

C.........  This module is required for the FileSetAPI
        USE MODFILESET, ONLY : FILE_INFO, RNAMES, NVARSET, VNAMESET, 
     &                         VUNITSET, VDESCSET

        IMPLICIT NONE
 
C...........   INCLUDES:
        INCLUDE 'EMCNST3.EXT'
        INCLUDE 'PARMS3.EXT'
        INCLUDE 'IODECL3.EXT'
        INCLUDE 'FDESC3.EXT'
        INCLUDE 'SETDECL.EXT'   !  FileSetAPI variables and functions

C...........   EXTERNAL FUNCTIONS
        CHARACTER(2)  CRLF
        CHARACTER(14) MMDDYY
        LOGICAL       BLKORCMT
        LOGICAL       ENVYN, GETYN
        INTEGER       GETFLINE, GETEFILE
        INTEGER       INDEX1
        INTEGER       LBLANK
        INTEGER       PROMPTFFILE
        CHARACTER(16) PROMPTMFILE
        INTEGER       SEC2TIME
        INTEGER       SECSDIFF
        REAL          STR2REAL
        LOGICAL       CHKREAL
        LOGICAL       SETENVVAR

        EXTERNAL CRLF, ENVYN, GETFLINE, GETYN, INDEX1, LBLANK,
     &           PROMPTFFILE, PROMPTMFILE, SEC2TIME, SECSDIFF,
     &           BLKORCMT, STR2REAL, CHKREAL, MMDDYY, SETENVVAR

C.........  LOCAL PARAMETERS and their descriptions:

        CHARACTER(50), PARAMETER :: 
     &  CVSW = '$Name SMOKEv4.9_Jun2022$' ! CVS release tag

C...........   LOCAL VARIABLES and their descriptions:

C...........   Emissions arrays
        REAL, ALLOCATABLE :: E2D ( : )        ! 2-d emissions
        REAL, ALLOCATABLE :: EOUT( :,: )      ! output emissions
        REAL, ALLOCATABLE :: BEFORE_ADJ( : )  ! emissions before factors applied
        REAL, ALLOCATABLE :: AFTER_ADJ ( : )  ! emissions after factors applied
        REAL, ALLOCATABLE :: BEFORE_SPC( : )  ! emissions before factors applied
        REAL, ALLOCATABLE :: AFTER_SPC ( : )  ! emissions after factors applied

C...........   Input file descriptors
        INTEGER,       ALLOCATABLE :: DURATA( : ) ! no. time steps
        INTEGER,       ALLOCATABLE :: NCOLSA( : ) ! no. columns
        INTEGER,       ALLOCATABLE :: NROWSA( : ) ! no. rows
        INTEGER,       ALLOCATABLE :: NVARSA( : ) ! no. variables
        INTEGER,       ALLOCATABLE :: SDATEA( : ) ! start date
        INTEGER,       ALLOCATABLE :: STIMEA( : ) ! start time
        INTEGER,       ALLOCATABLE :: NLAYSA( : ) ! number of layers in the file
        INTEGER,       ALLOCATABLE :: NFILES( : ) ! number of files in each fileset
        CHARACTER(16), ALLOCATABLE :: IONAME ( : ) ! IOAPI 16chr 2-d input file names
        CHARACTER(32), ALLOCATABLE :: FNAME ( : ) ! 2-d input file names
        LOGICAL,       ALLOCATABLE :: USEFIRST(:) ! true: use first time step of file
        LOGICAL,       ALLOCATABLE :: LVOUTA( :,: ) ! iff out var in input file
        CHARACTER(16), ALLOCATABLE :: VNAMEA( :,: ) ! variable names
        CHARACTER(16), ALLOCATABLE :: VUNITA( :,: ) ! variable units
        CHARACTER(80), ALLOCATABLE :: VDESCA( :,: ) ! var descrip
        REAL,          ALLOCATABLE :: ADJ_FACTOR( : ) ! adjustment factors
        CHARACTER(32), ALLOCATABLE :: ADJ_LFN( : )    ! Species name
        CHARACTER(16), ALLOCATABLE :: ADJ_SPC( : )    ! logicalFileName
        CHARACTER(49), ALLOCATABLE :: ADJ_LFNSPC( : ) ! concatenated {logicalFileName}_{Species}
        CHARACTER(32), ALLOCATABLE :: TAG_LFN( : )    ! Tagging species name
        CHARACTER(16), ALLOCATABLE :: TAG_SPC( :    ) ! Tagging logicalFileName
        CHARACTER(49), ALLOCATABLE :: TAG_LFNSPC( : ) ! Tagging concatenated {logicalFileName}_{Species}
        CHARACTER(16), ALLOCATABLE :: TAG_APPEND( : )  ! Tagging index
        CHARACTER(66), ALLOCATABLE :: TAG_LFNSPCTAG( : )  ! Tagging index(logicalfile+species+tag)
        CHARACTER(16)                 VNAMEP( MXVARS3 ) ! pt variable names
        CHARACTER(16)                 VUNITP( MXVARS3 ) ! pt variable units
        CHARACTER(80)                 VDESCP( MXVARS3 ) ! pt var descrip

C...........   Intermediate output variable arrays
        INTEGER       INDXN ( MXVARS3 ) ! sorting index for OUTIDX
        INTEGER       OUTIDX( MXVARS3 ) ! index to master model species list

        CHARACTER(16) OUTNAM( MXVARS3 ) ! unsorted output variable names
        CHARACTER(16) VUNITU( MXVARS3 ) ! unsorted output variable units
        CHARACTER(80) VDESCU( MXVARS3 ) ! unsorted output variable descriptions

        LOGICAL       LVOUTP( MXVARS3 ) ! iff output var exists in point input

C...........   Logical names and unit numbers

        INTEGER       ADEV            ! unit for logical names list for SEG
        INTEGER       IDEV            ! unit for logical names list for 2d files
        INTEGER       LDEV            ! unit for log file
        INTEGER       RDEV            ! unit for merge report file
        INTEGER       ODEV            ! unit for QA report file
        INTEGER       SDEV            ! unit for overall QA report file
        INTEGER       TDEV            ! unit for taggin input file
        INTEGER       GDEV            ! unit for taggin species QA file
        CHARACTER(16) ONAME           ! Merged output file name
        CHARACTER(16) PNAME           ! Point source input file name 

C...........   Other local variables 
        INTEGER       C, DD, F, I, J, K, L, L1, L2, N, NL, V, T ! pointers and counters

        INTEGER       ADJ                        ! tmp adjustment factor main index
        INTEGER       ADJ1                       ! tmp adjustment factor index 1
        INTEGER       ADJ2                       ! tmp adjustment factor index 2
        INTEGER       TAG                        ! tmp tagging species index
        INTEGER       DUMMY                      ! dummy value for use with I/O API functions
        INTEGER       EDATE                      ! ending julian date
        INTEGER       ETIME                      ! ending time HHMMSS
        INTEGER    :: G_SDATE = 0                ! start date from environment
        INTEGER    :: G_STIME = 0                ! start time from environment
        INTEGER    :: G_NSTEPS = 1               ! number of time steps from environment
        INTEGER    :: G_TSTEP = 0                ! time step from environment
        INTEGER       ICNTFIL                    ! tmp count of fileset file count  
        INTEGER       IOS                        ! i/o status
        INTEGER       IREC                       ! line number count
        INTEGER       JDATE                      ! iterative julian date
        INTEGER       JTIME                      ! iterative time HHMMSS
        INTEGER       LB                         ! leading blanks counter
        INTEGER       LE                         ! location of end of string
        INTEGER       MXNF                       ! tmp no. of 2-d input files
        INTEGER       MXNFIL                     ! max no. of 2-d input files
        INTEGER       MXNFAC                     ! max no. of adjustment factors
        INTEGER       MXNTAG                     ! max no. of tagging species
        INTEGER    :: NADJ = 0                   ! no. of adjustment factors
        INTEGER    :: NTAG = 0                   ! no. of tagging species
        INTEGER       NFILE                      ! no. of 2-d input files
        INTEGER       NSTEPS                     ! no. of output time steps
        INTEGER       NVOUT                      ! no. of output variables
        INTEGER       RDATE                      ! reference date
        INTEGER       SAVLAYS                    ! number of layers
        INTEGER       SDATE                      ! starting julian date
        INTEGER       SECS                       ! tmp seconds
        INTEGER       SECSMAX                    ! seconds maximum
        INTEGER       SECSMIN                    ! seconds minimum
        INTEGER       STIME                      ! starting time HHMMSS
        INTEGER       STEPS                      ! tmp number of steps
        INTEGER       TIMET                      ! tmp time from seconds
        INTEGER       TSTEP                      ! time step
        INTEGER       VLB                        ! VGLVS3D lower bound 

        REAL       :: FACS = 1.0                 ! adjustment factor 
        REAL          RATIO                      ! ratio 

        CHARACTER(16)  FDESC                     ! tmp file description
        CHARACTER(16)  MRGFDESC                  ! name for file description EV
        CHARACTER(128) METADESC                  ! output meta description from MRGFDESC
        CHARACTER(16)  IO_NAM                    ! tmp 16 chr logical file name
        CHARACTER(32)  NAM                       ! tmp logical file name
        CHARACTER(32)  LNAM                      ! tmp previous file name
        CHARACTER(16)  VNM                       ! tmp variable name
        CHARACTER(16)  TVNM                      ! tmp2 variable name
        CHARACTER(49)  LFNSPC                    ! tmp spec and file name
        CHARACTER(66)  LFNSPCTAG                 ! tmp speC, file, and tag name
        CHARACTER(256) LINE                      ! input buffer
        CHARACTER(256) NAMBUF                    ! tmp buffer for logical file name
        CHARACTER(256) MESG                      ! message field
        CHARACTER(80)  NAME1                     ! tmp file name component
        CHARACTER(15)  RPTCOL                    ! single column in report line
        CHARACTER(10)  EFMT                      ! output emissions foamat
        CHARACTER(100) REPFMT                    ! output emissions foamat
        CHARACTER(300) REPFILE                   ! name of report file
        CHARACTER(300) RPTLINE                   ! line of report file
        CHARACTER(16)  SPCTMP                    ! tmp species name
        CHARACTER(32)  LFNTMP                    ! tmp file name

        LOGICAL    :: EFLAG   = .FALSE.   ! error flag
        LOGICAL    :: FIRST3D = .TRUE.    ! true: first 3-d file not yet input
        LOGICAL    :: LFLAG   = .FALSE.   ! true  if 3-d file input
        LOGICAL    :: TFLAG   = .FALSE.   ! true: grid didn't match
        LOGICAL       MRGDIFF             ! true: merge files from different days

        CHARACTER(16) :: PROGNAME = 'MRGGRID' ! program name
C***********************************************************************
C   begin body of program MRGGRID
 
        LDEV = INIT3()
 
C.........  Write out copyright, version, web address, header info, and prompt
C           to continue running the program.
        CALL INITEM( LDEV, CVSW, PROGNAME )

C.........  Read names of input files and open files
        MESG = 'Enter logical name for 2-D AND 3-D GRIDDED INPUTS list'

        IDEV = PROMPTFFILE( MESG, .TRUE., .TRUE.,
     &                      'FILELIST', PROGNAME   )

C.........  Get environment variables
        MESG = 'Merge files from different days into single file'
        MRGDIFF = ENVYN( 'MRG_DIFF_DAYS', MESG, .FALSE., IOS )

        IF( MRGDIFF ) THEN        
C.............  Get date and time settings from environment
            CALL GETM3EPI( -1, G_SDATE, G_STIME, G_TSTEP, G_NSTEPS )        
        END IF

C.........  Determine maximum number of input files in file
        MXNFIL = GETFLINE( IDEV, 'List of files to merge' )

C.........  Write message out about MXVARS3
        WRITE( MESG,94010 ) 'Mrggrid compiled with I/O API MXVARS3 =',
     &                      MXVARS3
        CALL M3MSG2( MESG )

C.........  Allocate memory for arrays that just depend on the maximum number
C           of files
        ALLOCATE( NFILES( MXNFIL ), STAT=IOS )
        CALL CHECKMEM( IOS, 'NFILES', PROGNAME )
        ALLOCATE( DURATA( MXNFIL ), STAT=IOS )
        CALL CHECKMEM( IOS, 'DURATA', PROGNAME )
        ALLOCATE( NCOLSA( MXNFIL ), STAT=IOS )
        CALL CHECKMEM( IOS, 'NCOLSA', PROGNAME )
        ALLOCATE( NROWSA( MXNFIL ), STAT=IOS )
        CALL CHECKMEM( IOS, 'NROWSA', PROGNAME )
        ALLOCATE( NLAYSA( MXNFIL ), STAT=IOS )
        CALL CHECKMEM( IOS, 'NLAYSA', PROGNAME )
        ALLOCATE( NVARSA( MXNFIL ), STAT=IOS )
        CALL CHECKMEM( IOS, 'NVARSA', PROGNAME )
        ALLOCATE( SDATEA( MXNFIL ), STAT=IOS )
        CALL CHECKMEM( IOS, 'SDATEA', PROGNAME )
        ALLOCATE( STIMEA( MXNFIL ), STAT=IOS )
        CALL CHECKMEM( IOS, 'STIMEA', PROGNAME )
        ALLOCATE( FNAME( MXNFIL ), STAT=IOS )
        CALL CHECKMEM( IOS, 'FNAME', PROGNAME )
        ALLOCATE( IONAME( MXNFIL ), STAT=IOS )
        CALL CHECKMEM( IOS, 'IONAME', PROGNAME )
        ALLOCATE( LVOUTA( MXVARS3,MXNFIL ), STAT=IOS )
        CALL CHECKMEM( IOS, 'LVOUTA', PROGNAME )
        ALLOCATE( VNAMEA( MXVARS3,MXNFIL ), STAT=IOS )
        CALL CHECKMEM( IOS, 'VNAMEA', PROGNAME )
        ALLOCATE( VUNITA( MXVARS3,MXNFIL ), STAT=IOS )
        CALL CHECKMEM( IOS, 'VUNITA', PROGNAME )
        ALLOCATE( VDESCA( MXVARS3,MXNFIL ), STAT=IOS )
        CALL CHECKMEM( IOS, 'VDESCA', PROGNAME )
        
        IF( MRGDIFF ) THEN
            ALLOCATE( USEFIRST( MXNFIL ), STAT=IOS )
            CALL CHECKMEM( IOS, 'USEFIRST', PROGNAME )
        END IF

C.........  Allocate output layer structure
        ALLOCATE( VGLVS( 0:MXLAYS3 ), STAT=IOS )
        CALL CHECKMEM( IOS, 'VGLVS', PROGNAME )
        VGLVS = 0.

C.........  Loop through input files and open them
        F = 0
        IREC = 0
        DO

C.............  Read file names - exit if read is at end of file
            READ( IDEV, 93000, END=27, IOSTAT=IOS ) LINE
            IREC = IREC + 1

            IF ( IOS .NE. 0 ) THEN
                EFLAG = .TRUE.
                WRITE( MESG,94010 ) 
     &              'I/O error', IOS, 
     &              'reading file list at line', IREC
                CALL M3MESG( MESG )
                CYCLE
            END IF

C.............  Skip blank and comment lines
            IF ( BLKORCMT( LINE ) ) CYCLE

            F = F + 1

            IF( F .LE. MXNFIL ) THEN

                LB = LBLANK ( LINE )
                LE = LEN_TRIM( LINE )
                FNAME( F ) = LINE( LB+1:LE )

C.................  Re-set logical file name
                IF( LE > 16 ) THEN
                    CALL GETENV( FNAME( F ), NAMBUF )
                    
                    WRITE( IO_NAM,93030 ) FNAME( F )( 1:13 ) // '_', F
                    
                    IF( .NOT. SETENVVAR( IO_NAM, NAMBUF ) ) THEN
                        MESG = 'Could not set logical file name for ' //
     &                         'file ' // TRIM( NAMBUF )
                        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                    END IF

                    IONAME( F ) = IO_NAM

                ELSE

                    IO_NAM = TRIM( FNAME( F ) )
                    IONAME( F ) = IO_NAM

                ENDIF

                IF ( .NOT. OPENSET( IONAME(F), FSREAD3, PROGNAME )) THEN
 
                    MESG = 'Could not open file "' //
     &                     FNAME( F )( 1:LEN_TRIM( FNAME(F) ) )// '".'
                    CALL M3MSG2( MESG )
                    MESG = 'Ending program "MRGGRID".'
                    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

                END IF      !  if open3() failed

C.................  Store whether it's a fileset file or not
                I = INDEX1( FNAME(F), MXFILE3, RNAMES )
                NFILES( F ) = SIZE( FILE_INFO( I )%LNAMES )

            END IF

        END DO
27      CONTINUE

        NFILE = F

        IF( NFILE .GT. MXNFIL ) THEN
            WRITE( MESG,94010 )
     &        'INTERNAL ERROR: Dimension mismatch.  Input file count:',
     &        NFILE, 'program allows (MXNFIL):', MXNFIL
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

        ELSEIF( NFILE .EQ. 0 ) THEN
            MESG = 'No input files in list!'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

        ENDIF

C.........  Get environment variable settings for adjustment factor input file
        CALL ENVSTR( 'ADJ_FACS', MESG, ' ', NAME1 , IOS )

C.........  Determine maximum number of input files in file
        IF( IOS < 0 ) THEN     !  failure to open
            ADEV = IOS
            MXNFAC = 1
            MESG = 'NOTE : No adjustment factors were applied because'//
     &             ' there is no ADJ_FACS environment variable defined' 
            CALL M3MSG2( MESG )

        ELSE
            MESG = 'Enter logical name for a list of adjustment factors'
            ADEV = PROMPTFFILE( MESG,.TRUE.,.TRUE.,'ADJ_FACS',PROGNAME )
            MXNFAC = GETFLINE( ADEV, 'List of adjustment factos' )

C............  Write summary of sector specific factor adjustment output
            ODEV = PROMPTFFILE(
     &         'Enter logical name for the MRGGRID QA REPORT file',
     &         .FALSE., .TRUE., 'REPMERGE_ADJ', PROGNAME )

            MXNF = 0
            DO         ! head of report file
                READ( ODEV, 93000, END=222 ) LINE
                MXNF = MXNF + 1
            ENDDO
222         CONTINUE

C.............  Write header line to report     
            IF( MXNF == 0 ) THEN
              WRITE( ODEV,93000 ) '#MRGGRID logical file QA Report'
              WRITE( ODEV,93000 ) '#COLUMN_TYPES=Int(4)|Varchar(32)|' // 
     &                     'Varchar(16)|Real(8)|Real(8)|Real(8)|Real(8)'
              WRITE( ODEV,93000 ) 'DATE,FileName,Species,Factor,'//
     &                          'Before,After,Ratio'
            END IF

C............  Write summary of overall factor adjustment output by species
            SDEV = PROMPTFFILE(
     &         'Enter logical name for the MRGGRID Overall REPORT file',
     &         .FALSE., .TRUE., 'REPMERGE_SUM', PROGNAME ) 

            MXNF = 0
            DO         ! head of report file
                READ( SDEV, 93000, END=333 ) LINE
                MXNF = MXNF + 1
            ENDDO
333         CONTINUE

C.............  Write header line to report     
            IF( MXNF == 0 ) THEN
              WRITE( SDEV,93000 ) '#MRGGRID Overall Summary by Species'
              WRITE( SDEV,93000 ) '#COLUMN_TYPES=Int(4)|Varchar(16)|' //
     &                          'Real(8)|Real(8)|Real(8)'
              WRITE( SDEV,93000 ) 'DATE,Species,Before,After,Ratio'
            END IF

        END IF

C.........  Allocate memory for arrays that just depend on the maximum number
C           of adjustment factors in ADJ_FACS input file.
        ALLOCATE( ADJ_LFN( MXNFAC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'ADJ_LFN', PROGNAME )
        ALLOCATE( ADJ_SPC( MXNFAC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'ADJ_SPC', PROGNAME )
        ALLOCATE( ADJ_LFNSPC( MXNFAC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'ADJ_LFNSPC', PROGNAME )
        ALLOCATE( ADJ_FACTOR( MXNFAC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'ADJ_FACTOR', PROGNAME )

        ADJ_SPC = ' '
        ADJ_LFN = ' '
        ADJ_LFNSPC = ' '
        ADJ_FACTOR = 0.0

C.........  Define a number of adjustment factors
        IF( ADEV < 0 ) THEN 
            NADJ = 1

        ELSE
C.............  Store a list of adjustment factors
            CALL READ_ADJ_FACS( NADJ )

        END IF

C.........  Duplicate Check of ADJ_FACS file
        DO  F = 1, NADJ
            LFNSPC = ADJ_LFNSPC( F )
            DD = 0
            DO I = 1, NADJ
                IF( LFNSPC == ADJ_LFNSPC( I ) ) DD = DD + 1
            END DO
	    
            IF( DD > 1 ) THEN
                MESG = 'ERROR: Duplicate entries of '// TRIM(ADJ_SPC(F))
     &               // ' species from the ' // TRIM( ADJ_LFN(F) ) // 
     &              ' file in the ADJ_FACS file.' // LFNSPC
                CALL M3MSG2( MESG )
                EFLAG = .TRUE.
            END IF
        ENDDO

C.........  Give error message and end program unsuccessfully
        IF( EFLAG ) THEN
            MESG = 'ERROR: Duplicate entries in the ADJ_FACS file'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

C.........  Allocate arrays that will store sector-specific daily/gridded total emissinos
        ALLOCATE( BEFORE_ADJ( NADJ ), STAT=IOS )
        CALL CHECKMEM( IOS, 'BEFORE_ADJ', PROGNAME )
        ALLOCATE( AFTER_ADJ( NADJ ), STAT=IOS )
        CALL CHECKMEM( IOS, 'AFTER_ADJ', PROGNAME )
        BEFORE_ADJ = 0.0
        AFTER_ADJ  = 0.0

C.........  Get environment variable settings for tagging species input file
        CALL ENVSTR( 'TAG_SPECIES', MESG, ' ', NAME1 , IOS )

C.........  Determine maximum number of input files in file
        IF( IOS < 0 ) THEN     !  failure to open
            TDEV = IOS
            MESG = 'NOTE : No tagging species were available because'//
     &             ' there is no TAG_SPECIES environment variable defined' 
            CALL M3MSG2( MESG )

        ELSE

            MESG = 'NOTE : Tagging species based upon TAG_SPECIES file' 
            CALL M3MSG2( MESG )

C.............  Store a list of tagging species
            CALL READ_TAG_SPECIES( NTAG )

C............  Write summary of sector specific factor adjustment output
            GDEV = PROMPTFFILE(
     &         'Enter logical name for the MRGGRID Tagging REPORT file',
     &         .FALSE., .TRUE., 'REPMERGE_TAG', PROGNAME )

C.............  Write header line to report     
            WRITE( GDEV,93000 ) '#MRGGRID Tagging species Report'
            WRITE( GDEV,93000 ) '#COLUMN_TYPES=Varchar(32)|' // 
     &                          'Varchar(16)|Varchar(16)'
            WRITE( GDEV,93000 ) 'FileName,OriginalSpecies,TaggedSpecies'

        END IF

C.........  Determine I/O API layer storage lower bound
        VLB = LBOUND( VGLVS3D,1 )

C.........  Get file descriptions and store for all input files
C.........  Loop through 2D input files
        NLAYS = 1
        DO F = 1, NFILE

            NAM    = FNAME( F )
            IO_NAM = IONAME( F )   ! retrieve 16 char ioapi local file name

            ICNTFIL = ALLFILES
            IF( NFILES( F ) .EQ. 1 ) ICNTFIL = 1   ! send ALLFILES if more than one file, send 1 otherwise
            IF ( .NOT. DESCSET( IO_NAM, ICNTFIL ) ) THEN
                MESG = 'Could not get description of file "'  //
     &                  NAM( 1:LEN_TRIM( NAM ) ) // '"'
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            ELSE
                NROWSA( F ) = NROWS3D
                NCOLSA( F ) = NCOLS3D
                NLAYSA( F ) = NLAYS3D
                NVARSA( F ) = NVARSET
                SDATEA( F ) = SDATE3D
                STIMEA( F ) = STIME3D
                DURATA( F ) = MXREC3D
                
                IF( F == 1 ) TSTEP = TSTEP3D
                
                DO V = 1, NVARSET
                    VNAMEA( V,F ) = VNAMESET( V )
                    VUNITA( V,F ) = VUNITSET( V )
                    VDESCA( V,F ) = VDESCSET( V )

                    IF( TDEV > 0 ) THEN
C.........................  Set tmp variables
                        VNM = VNAMESET( V )

C..........................  Search tagged species for the current file
                        LFNSPC = TRIM( NAM ) // '~' // TRIM( VNM )

                        TAG = INDEX1( LFNSPC, NTAG, TAG_LFNSPC )
                        ADJ1= INDEX1( LFNSPC, NADJ, ADJ_LFNSPC )

C.........................  Assign tagged species for the current species
                        IF( TAG > 0 ) THEN
                            TVNM = TRIM( VNM ) // '_' // 
     &                             TRIM( TAG_APPEND( TAG ) )

                            MESG = 'NOTE : Appending a tag (' // 
     &                          TRIM( TAG_APPEND(TAG) ) // ') to the '
     &                          //'species ' //TRIM( VNM )// ' from the '
     &                          //TRIM( NAM )// ' file'
                            CALL M3MSG2( MESG )

C.............................  Replace a species name with a tagged one
                            VNAMEA( V,F ) = TVNM

C.............................  Write the changes to tagging summary report
                            WRITE( GDEV, 92000 ) NAM, VNM, TVNM

C.............................  Error and Warning messages for the tagged species
C                               before you apply the adjustment factors if necessary

                            LFNSPC = TRIM( NAM ) // '~' // TRIM( TVNM )

                            ADJ2 = INDEX1( LFNSPC, NADJ, ADJ_LFNSPC )

                            IF( ADJ1 > 0 .AND. ADJ2 < 1 ) THEN
                                MESG ='WARNING : Adjustment factor ' //
     &                              ' for the species ' // TRIM( VNM )
     &                              // ' from file ' // TRIM( NAM ) // 
     &                             ' will be skipped due to ' //
     &                              'the change of species name to ' //
     &                              TRIM( TVNM )
                                CALL M3MSG2( MESG )
                                
                            END IF

                        END IF
                        
                    END IF

                END DO

            END IF

C.............  Search for tag species in the logical file 
            IF( TDEV > 0 ) THEN
                DO N = 1, NTAG
                    LFNTMP = TAG_LFN( N )  ! retriev logical file name from TAG_SPECIES

                    IF( LFNTMP == NAM ) THEN

                        SPCTMP = TAG_SPC( N ) ! retrieve spcieces name from TAG_SPECIES 

                        K = INDEX1( SPCTMP, NVARSET, VNAMESET )

                        IF( K <= 0 ) THEN
                            EFLAG = .TRUE.
                            MESG = 'ERROR: The species ' //
     &                          TRIM( TAG_SPC(N) )//' you want to tag '//
     &                          'is not available from file '//TRIM(NAM)
                            CALL M3MSG2( MESG )
                        END IF

                    END IF

                END DO

            END IF

C.............  Search for adj factor species in the logical file
            DO J = 1, NADJ

                LFNTMP = ADJ_LFN( J ) ! retriev logical file from ADJ_FACS

                IF( LFNTMP == NAM ) THEN
                    SPCTMP = ADJ_SPC( J ) ! retrieve spcieces name from ADJ_FACS

                    IF( TDEV > 0 ) THEN

C..........................  Search tagged species for the current file
                        LFNSPC = TRIM( LFNTMP ) // '~' // TRIM( SPCTMP )
                        TAG = INDEX1( LFNSPC, NTAG, TAG_LFNSPCTAG )
     
C.........................  Assign adjustment factor for the current species
                        IF( TAG > 0 ) SPCTMP = TAG_SPC( TAG )
                        
                    END IF

                    K = INDEX1( SPCTMP, NVARSET, VNAMESET )
                    IF( K <= 0 ) THEN
                        EFLAG = .TRUE.
                        MESG = 'ERROR: The species '//TRIM(ADJ_SPC(J))// 
     &                  ' you want to adjust is not available in the '//
     &                  TRIM(NAM) // ' file'
                        CALL M3MSG2( MESG )
                    END IF

                END IF

            END DO

C.............  Compare all other time steps back to first file.
C.............  They must match exactly.
            IF( TSTEP3D /= TSTEP ) THEN
                EFLAG = .TRUE.
                WRITE( MESG,94010 ) 'ERROR: Time step', TSTEP3D,
     &                 'in file "' // TRIM( NAM ) //
     &                 '" is inconsistent with first file value of',
     &                 TSTEP
                CALL M3MSG2( MESG )
            END IF

C.............  Compare all other grids back to first grid.
C.............  They must match exactly.
            WRITE( FDESC, '(A,I3.3)' ) 'FILE', F
            TFLAG = .FALSE.
            CALL CHKGRID( FDESC, 'GRID', 0, TFLAG )

            IF( TFLAG ) THEN
                EFLAG = .TRUE.
                L = LEN_TRIM( NAM )
                WRITE( MESG,94010 ) 'ERROR: File "' // NAM( 1:L ) //
     &            '" (NX,NY)  : (', NCOLSA( F ), ',', NROWSA( F ), ')'//
     &            CRLF() // BLANK10 // 'is inconsistent with first ' //
     &            'file (NX,NY) : (', NCOLS, ',', NROWS, ')'
                CALL M3MSG2( MESG )
            END IF

C.............  Compare layer structures for 3-d files. The number of layers do 
C               not need to match, but the layer structures do need to match.
            NLAYS = MAX( NLAYS, NLAYSA( F ) )
            IF( NLAYSA( F ) .GT. 1 ) THEN
                LFLAG = .TRUE.

C.................  For the first file that is 3-d, initialize output layer
C                   structure       
                IF ( FIRST3D ) THEN

                    NLAYS = NLAYSA( F )
                    VGTYP = VGTYP3D
                    VGTOP = VGTOP3D
                    VGLVS( 0:NLAYS ) = VGLVS3D( 0+VLB:NLAYS+VLB )   ! array
                    FIRST3D = .FALSE.

C.................  For additional 3-d files, compare the layer structures
                ELSE

C.....................  Check vertical type
                    IF( VGTYP3D .NE. VGTYP ) THEN
                        EFLAG = .TRUE.
                        L = LEN_TRIM( NAM )
                        WRITE( MESG, 94010 ) 'ERROR: Vertical ' //
     &                         'coordinate type', VGTYP3D, 
     &                         'in file "'// NAM(1:L) //
     &                         '" is inconsistent with first 3-d'//
     &                         'file value of', VGTYP
                        CALL M3MSG2( MESG )
                    END IF

C.....................  Check vertical top
                    IF( VGTOP3D .NE. VGTOP ) THEN
                        EFLAG = .TRUE.
                        L = LEN_TRIM( NAM )
                        WRITE( MESG, 94010 ) 'ERROR: Vertical ' //
     &                         'top value', VGTOP3D, 
     &                         'in file "'// NAM(1:L) //
     &                         '" is inconsistent with first 3-d'//
     &                         'file value of', VGTOP
                        CALL M3MSG2( MESG )
                    END IF

C.....................  Loop through layers of current file F
                    DO NL = 0, NLAYSA( F )

C.........................  For layers that are common to this file and previous
C                           files
                        IF( NL .LE. NLAYS ) THEN

                            IF( VGLVS3D( NL+VLB ) .NE. VGLVS( NL )) THEN
                                EFLAG = .TRUE.
                                L = LEN_TRIM( NAM )
                                WRITE( MESG, 94020 ) 'ERROR: Layer', NL,
     &                            'in file "'// NAM( 1:L ) // 
     &                            '" with level value', VGLVS3D(NL+VLB), 
     &                            CRLF()//BLANK10//'is inconsistent '//
     &                            'with first file value of', VGLVS(NL)
                                CALL M3MSG2( MESG )
                            END IF

C.........................  Add additional layers from current file to output 
C                           layer structure
                        ELSE
                            VGLVS( NL ) = VGLVS3D( NL+VLB )
                        END IF

                    END DO    ! End checking layers

C.....................  Reset the global number of layers as the maximum between
C                       the current file and all previous files
                    NLAYS = MAX( NLAYS, NLAYSA( F ) )

                END IF        ! End first 3-d file or not
            END IF            ! End 3-d files

        END DO                ! End loop through files

C.........  Give error message and end program unsuccessfully
        IF( EFLAG ) THEN
            MESG = 'Inconsistent time step, grid, species, or layers '//
     &              'among the files!'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

C.........  Check that environment settings are consistent with files
        IF( MRGDIFF ) THEN
            IF( TSTEP /= G_TSTEP ) THEN
                WRITE( MESG,94010 ) 'ERROR: Value for G_TSTEP ',
     &              G_TSTEP, 'is inconsistent with the time step' //
     &              CRLF() // BLANK10 // 'of the input files', TSTEP
                CALL M3MSG2( MESG )
                
                MESG = 'Inconsistent environment settings'
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF
        END IF

C.........  Deterimine output date, time, and number of time steps
        SDATE = G_SDATE
        STIME = G_STIME
        NSTEPS = G_NSTEPS
        CALL SETOUTDATE( SDATE, STIME, NSTEPS, NFILE, SDATEA,
     &                   STIMEA, DURATA, FNAME, MRGDIFF, USEFIRST )

C.........  Build master output variables list
        NVOUT = 0

C.........  Loop through input files and build an output variable list
        DO F = 1, NFILE

C.............  Loop through variables in the files
            DO V = 1, NVARSA( F )

                VNM = VNAMEA( V,F )

C.................  Look for variable name in output list
                K = INDEX1( VNM, NVOUT, OUTNAM  )  ! look in output list

C.................  If its not in the output list, add it
                IF( K .LE. 0 ) THEN
                    NVOUT = NVOUT + 1
                    INDXN ( NVOUT ) = NVOUT
                    OUTNAM( NVOUT ) = VNM
                    VDESCU( NVOUT ) = VDESCA( V,F )
                    VUNITU( NVOUT ) = VUNITA( V,F )

C.................  If variable is in the output list, check the units
                ELSE
                    IF ( VUNITA( V,F ) .NE. VUNITU( K ) ) THEN
                        EFLAG = .TRUE.
                        L  = LEN_TRIM( VNM )
                        L1 = LEN_TRIM( VUNITA( V,F ) )
                        L2 = LEN_TRIM( VUNITU( K )   )
                        WRITE( MESG,94010 ) 'ERROR: Variable "' //
     &                         VNM( 1:L ) // '" in file', F,
     &                         'has units "'// VUNITA( V,F )( 1:L1 ) //
     &                         '"' // CRLF() // BLANK10 //
     &                         'that are inconsistent with a '//
     &                         'previous file that had units "' //
     &                         VUNITU(K)( 1:L2 )// '" for this variable'
                        CALL M3MSG2( MESG )

                    END IF  ! End check of units

                END IF      ! End variable in output list already or not

            END DO          ! End loop through variables in this file

        END DO              ! End loop through files.

C.........  Give error message and end program unsuccessfully
        IF( EFLAG ) THEN
            MESG = 'Inconsistent units for common variables among '//
     &             'the files!'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

C.........  Sort output variables into alphabetical order
        CALL SORTIC( NVOUT, INDXN, OUTNAM )

C.........  Set up for opening output file...
C.........  Get grid information
        ICNTFIL = ALLFILES
        IF( NFILES( 1 ) .EQ. 1 ) ICNTFIL = 1   ! send ALLFILES if more than one file, send 1 otherwise

        IF( .NOT. DESCSET( IONAME( 1 ), ICNTFIL ) ) THEN
            MESG = 'Could not get description of file "'  //
     &              FNAME( 1 )( 1:LEN_TRIM( FNAME(1) ) ) // '"'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        ENDIF

        SDATE3D = SDATE
        STIME3D = STIME
        NVARS3D = NVOUT

C.........  Set up layer structure for output file
        NLAYS3D = NLAYS
        VGTOP3D = VGTOP
        VGTYP3D = VGTYP
        VGLVS3D = 0.     ! initialize array
        DO NL = 0, NLAYS
            VGLVS3D( NL+VLB ) = VGLVS( NL )
        END DO

C.........  Set the EV name for met data description
        MESG = 'Setting for the environment variable name for meta ' //
     &           'file description for output file'
        CALL ENVSTR( 'MRG_FILEDESC', MESG, ' ', MRGFDESC, IOS  )
        FDESC3D( 1 ) = 'Merged emissions output file from Mrggrid'

        IF( IOS >= 0 ) THEN
            MESG = 'Use this meta file description for output file'
            CALL ENVSTR( MRGFDESC, MESG, ' ', METADESC, IOS )
            FDESC3D( 1 ) = METADESC
        END IF

C........  Update variable names in sorted order, and also 
C........  set up logical arrays for which files have which species
        DO V = 1, NVOUT
            VNM = OUTNAM( INDXN( V ) )

            VNAME3D( V ) = VNM                  ! store sorted output vars, etc.
            VDESC3D( V ) = VDESCU( INDXN( V ) )
            UNITS3D( V ) = VUNITU( INDXN( V ) )
            VTYPE3D( V ) = M3REAL

            DO F = 1, NFILE
                LVOUTA( V,F ) = .FALSE.

                J = INDEX1( VNM, NVARSA( F ), VNAMEA( 1,F ) )
                IF( J .GT. 0 ) LVOUTA( V,F ) = .TRUE.
                                
            END DO 
        END DO

C.........  Allocate memory for the number of grid cells and layers
        NGRID = NROWS * NCOLS
        ALLOCATE( E2D( NGRID ), STAT=IOS )
        CALL CHECKMEM( IOS, 'E2D', PROGNAME )
        ALLOCATE( EOUT( NGRID, NLAYS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'EOUT', PROGNAME )

C.........  Prompt for and open output file
        ONAME = PROMPTMFILE( 
     &          'Enter logical name for MERGED GRIDDED OUTPUT file',
     &          FSUNKN3, 'OUTFILE', PROGNAME )

C.........  Prompt for and open report file
        IF( MRGDIFF ) THEN
            RDEV = PROMPTFFILE(
     &             'Enter logical name for the MRGGRID REPORT file',
     &             .FALSE., .TRUE., 'REPMRGGRID', PROGNAME ) 

C.............  Write header line to report     
            WRITE( RPTLINE,93010 ) 'Output date'
            WRITE( RPTCOL,93010 ) 'Output time'
            RPTLINE = TRIM( RPTLINE ) // RPTCOL
            
            DO F = 1, NFILE
                NAM = FNAME( F )
                WRITE( RPTCOL,93010 ) TRIM( NAM ) // ' date'
                RPTLINE = TRIM( RPTLINE ) // RPTCOL
            END DO
            
            WRITE( RDEV,93000 ) TRIM( RPTLINE )
        END IF

C.........  Warning missing logical file names from the ADJ_FACS list 
        LNAM = ' '
        DO I = 1, NADJ
             NAM = ADJ_LFN( I )
             L = INDEX1( NAM, NFILE, FNAME )
             IF( L <= 0 ) THEN
                 MESG = 'WARNING: The logical file '//TRIM(NAM) //
     &               ' in the adjustment file (ADJ_FACS) is not '// 
     &               'found in the FILELIST on DATE : '//
     &               MMDDYY(SDATE)
                 IF( LNAM /= NAM ) THEN
                     CALL M3MSG2( MESG )
                     LNAM = NAM
                 END IF
              END IF
        END DO

C.........  Warning missing logical file names from the TAG_SPECIES list 
        LNAM = ' '
        IF( TDEV > 0 ) THEN
            DO I = 1, NTAG
                NAM = TAG_LFN( I )
                L = INDEX1( NAM, NFILE, FNAME )
                IF( L <= 0 ) THEN
                     MESG = 'WARNING: The logical file '//TRIM(NAM) //
     &                   ' in the tagging file (TAG_SPECIES) is not '// 
     &                   'found in the FILELIST on DATE : '//
     &                   MMDDYY(SDATE)
                     IF( LNAM /= NAM ) THEN
                         CALL M3MSG2( MESG )
                         LNAM = NAM
                     END IF
                END IF
            END DO
        END IF

C.........  Allocate arrays that will store overall daily/gridded total emissinos by species
        ALLOCATE( BEFORE_SPC( NVOUT ), STAT=IOS )
        CALL CHECKMEM( IOS, 'BEFORE_SPC', PROGNAME )
        ALLOCATE( AFTER_SPC( NVOUT ), STAT=IOS )
        CALL CHECKMEM( IOS, 'AFTER_SPC', PROGNAME )
        BEFORE_SPC = 0.0
        AFTER_SPC  = 0.0

C.........  Loop through hours
        JDATE = SDATE
        JTIME = STIME
        FACS  = 1.0
        DO T = 1, NSTEPS

C.............  Loop through species
            DO V = 1, NVOUT

                VNM = VNAME3D( V ) 

C.................  Output array
                EOUT = 0.   ! array

                DO F = 1, NFILE

C.....................  Set read date
                    IF( MRGDIFF ) THEN
                      IF( USEFIRST( F ) ) THEN
                        DUMMY = 0
                        STEPS = SEC2TIME( 
     &                            SECSDIFF( 
     &                              SDATE, DUMMY, JDATE, DUMMY ) )
                        RDATE = SDATEA( F )
                        CALL NEXTIME( RDATE, DUMMY, STEPS )
                      END IF 
                    ELSE
                        RDATE = JDATE
                    END IF

C.....................  Set tmp variables
                    NAM    = FNAME ( F )   ! input file name
                    IO_NAM = IONAME( F )   ! retrieve 16 char ioapi input file name

                    NL  = NLAYSA( F )       ! number of layers

C.....................  Search adjustment factor for the current file
                    LFNSPC = TRIM( NAM ) // '~' // TRIM( VNM )

C.....................  Assign adjustment factor for the current species
                    ADJ = INDEX1( LFNSPC, NADJ, ADJ_LFNSPC )
                    IF( ADJ > 0 ) THEN
                        FACS = ADJ_FACTOR( ADJ )

                        WRITE( MESG,93011 )'Apply adjustment factor' ,
     &                      FACS, ' to the '  // TRIM( VNM ) //
     &                      ' species from the '//TRIM( NAM )// ' file'
                        CALL M3MSG2( MESG )
                    ELSE
                        FACS = 1.0
                       
                    END IF

C.....................  Search tagged species for the current file
                    IF( TDEV > 0 ) THEN
                        TAG = INDEX1( LFNSPC, NTAG, TAG_LFNSPCTAG )
                        IF( TAG > 0 ) THEN
                            TVNM = TAG_SPC( TAG )
                        ELSE
                            TVNM = VNM
                        END IF
                    ELSE
                        TVNM = VNM

                    END IF

C.....................  If file has species, read (do this for all files)...
                    IF( LVOUTA( V,F ) ) THEN

                        ICNTFIL = ALLFILES
                        IF( NFILES( F ) .EQ. 1 ) ICNTFIL = 1   ! send ALLFILES if more than one file, send 1 otherwise

C.........................  If 2-d input file, read, and add
                        IF( NL .EQ. 1 ) THEN
                            IF( .NOT. 
     &                           READSET( IO_NAM, TVNM, 1, ICNTFIL,
     &                                    RDATE, JTIME, E2D     )) THEN

                                MESG = 'Could not read "' // VNM //
     &                                 '" from file "' //
     &                                 NAM( 1:LEN_TRIM( NAM ) )// '".'
                                CALL M3EXIT( PROGNAME, RDATE, JTIME, 
     &                                       MESG, 2 )
                            ENDIF

C.............................  Logical file specific summary
                            IF( ADJ > 0 ) THEN
                                BEFORE_ADJ( ADJ ) = BEFORE_ADJ( ADJ ) + 
     &                                          SUM( E2D(1:NGRID) )

                                AFTER_ADJ ( ADJ ) = AFTER_ADJ ( ADJ ) + 
     &                                          SUM( E2D(1:NGRID)*FACS )

C.............................  Overall summary by species
                                BEFORE_SPC( V )  = BEFORE_SPC( V ) + 
     &                                          SUM( E2D(1:NGRID) )

                                AFTER_SPC ( V ) = AFTER_SPC ( V ) + 
     &                                          SUM( E2D(1:NGRID)*FACS )
                            END IF

                            EOUT( 1:NGRID,1 ) = EOUT( 1:NGRID,1 ) + 
     &                                          E2D( 1:NGRID) * FACS

C.........................  If 3-d input file, allocate memory, read, and add
                        ELSE

                            DO K = 1, NL
                                IF( .NOT. 
     &                               READSET( IO_NAM,TVNM,K,ICNTFIL,
     &                                        RDATE, JTIME, E2D  )) THEN

                                    MESG = 'Could not read "' // VNM //
     &                                     '" from file "' //
     &                                   NAM( 1:LEN_TRIM( NAM ) )// '".'
                                    CALL M3EXIT( PROGNAME, RDATE, JTIME,
     &                                           MESG, 2 )
                                END IF

C.................................  Logical file specific summary
                                IF( ADJ > 0 ) THEN

                                  BEFORE_ADJ( ADJ ) = BEFORE_ADJ( ADJ )+ 
     &                                              SUM( E2D(1:NGRID) )

                                  AFTER_ADJ ( ADJ ) = AFTER_ADJ ( ADJ )+ 
     &                                            SUM(E2D(1:NGRID)*FACS)

C.................................  Overall summary by species
                                  BEFORE_SPC( V )  = BEFORE_SPC( V ) + 
     &                                           SUM( E2D(1:NGRID) )

                                  AFTER_SPC ( V ) = AFTER_SPC ( V ) + 
     &                                           SUM(E2D(1:NGRID)*FACS)

                                END IF

                                EOUT( 1:NGRID,K )= EOUT( 1:NGRID,K ) + 
     &                                             E2D( 1:NGRID )*FACS
                            END DO

                        END IF  ! if 2-d or 3-d
                    END IF      ! if pollutant is in this file

C.....................  Build report line if needed
                    IF( MRGDIFF .AND. V == 1 ) THEN
                        IF( F == 1 ) THEN
                            WRITE( RPTLINE,93020 ) JDATE
                            WRITE( RPTCOL,93020 ) JTIME
                            RPTLINE = TRIM( RPTLINE ) // RPTCOL
                        END IF
                        
                        WRITE( RPTCOL,93020 ) RDATE
                        RPTLINE = TRIM( RPTLINE ) // RPTCOL
                    END IF

                END DO          ! loop over input files

C.................  Write species/hour to output file
                IF( .NOT. WRITE3( ONAME, VNM, JDATE, JTIME, EOUT )) THEN

                    MESG = 'Could not write "'// VNM// '" to file "'// 
     &                      ONAME( 1:LEN_TRIM( ONAME ) ) // '".'
     &                        
                    CALL M3EXIT( PROGNAME, JDATE, JTIME, MESG, 2 )

                END IF

            END DO   ! loop through variables

C.............  Write this time step to report
            IF( MRGDIFF ) THEN
                WRITE( RDEV,93000 ) TRIM( RPTLINE )
            END IF

            CALL NEXTIME( JDATE, JTIME, TSTEP )
      
        END DO       ! loop through timesteps

C........  Write summary of sector specific factor adjustment output
C          Columns: Date, Sector, Species, value before, value after, ratio of before/after
C          Later we can add the total amount of the adjusted species summed accross all of the input files
C          and the total amount of the adjusted species in the output file.

C.........  Write header line to report     
        DO F = 1, NADJ

            VNM = ADJ_SPC( F )     ! species name
            NAM = ADJ_LFN( F )     ! logical file name
            FACS   = ADJ_FACTOR( F )   ! adjustment factor

            IF( BEFORE_ADJ( F ) == 0.0 ) CYCLE

            RATIO = ( AFTER_ADJ( F ) / BEFORE_ADJ( F ) )

            REPFMT = "(I8,2(',',A),',',F10.6,',',"

C.............  Define the format of real values
            CALL GET_FORMAT( VNM, BEFORE_ADJ( F ), EFMT )
            REPFMT = TRIM( REPFMT ) // TRIM( EFMT )

            CALL GET_FORMAT( VNM, AFTER_ADJ( F ), EFMT )
            REPFMT = TRIM( REPFMT ) // TRIM( EFMT )
            REPFMT = TRIM( REPFMT ) // "F10.6)"

            WRITE( RPTLINE,REPFMT ) SDATE, NAM, VNM, FACS,
     &                            BEFORE_ADJ( F ), AFTER_ADJ( F ), RATIO

            IF( ADEV > 0 .AND. RATIO /= 1.0 ) THEN
                WRITE( ODEV,93000 ) TRIM( RPTLINE )
            ENDIF
        END DO

C.........  Write header line to overall summary report     
        DO V = 1, NVOUT

            VNM   = VNAME3D( V )     ! species name

            IF( BEFORE_SPC( V ) == 0.0 ) CYCLE

            RATIO = ( AFTER_SPC( V ) / BEFORE_SPC( V ) )

            REPFMT = "( I8,',',A,',',"

C.............  Define the format of real values
            CALL GET_FORMAT( VNM, BEFORE_SPC( V ), EFMT )
            REPFMT = TRIM( REPFMT ) // TRIM( EFMT )

            CALL GET_FORMAT( VNM, AFTER_SPC( V ), EFMT )
            REPFMT = TRIM( REPFMT ) // TRIM( EFMT )
            REPFMT = TRIM( REPFMT ) // "F10.6)"

            WRITE( RPTLINE,REPFMT ) SDATE, VNM, BEFORE_SPC(V),
     &                              AFTER_SPC(V),RATIO

            IF( SDEV > 0 .AND.  RATIO /= 1.0 ) THEN
                WRITE( SDEV,93000 ) TRIM( RPTLINE )
            END IF
        END DO

C......... Normal Completion
        CALL M3EXIT( PROGNAME, 0, 0, ' ', 0)
    
C******************  FORMAT  STATEMENTS   ******************************

C...........   Formatted file I/O formats............ 93xxx

92000   FORMAT(  A,',',A,',',A  )

93000   FORMAT(  A )

93010   FORMAT( A15 )

93011   FORMAT(  A, F8.5, A )

93020   FORMAT( I15 )

93030   FORMAT( A,I2.2 )

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I7, :, 1X ) )

94020   FORMAT( A, :, I3, :, 1X, 10 ( A, :, F8.5, :, 1X ) )


C*****************  INTERNAL SUBPROGRAMS  ******************************

        CONTAINS

C----------------------------------------------------------------------

C.............  This internal subprogram determines the format to output
C               emission values
            SUBROUTINE GET_FORMAT( VBUF, VAL, FMT )

C.............  Subroutine arguments
            CHARACTER(*), INTENT (IN) :: VBUF
            REAL        , INTENT (IN) :: VAL
            CHARACTER(*), INTENT(OUT) :: FMT

C----------------------------------------------------------------------

C.............  Value is too large for 
            IF( VAL .GT. 999999999. ) THEN
                FMT = "E10.3,',',"

                L = LEN_TRIM( VBUF )
                WRITE( MESG,95020 ) 'WARNING: "' // VBUF( 1:L ) // 
     &              '"Emissions value of', VAL, CRLF()// BLANK10// 
     &              '" is too large for file format, so writing ' //
     &              'in scientific notation for source'

                CALL M3MESG( MESG )

            ELSE IF( VAL .GT. 99999999. ) THEN
                FMT = "F10.0,',',"

            ELSE IF( VAL .GT. 9999999. ) THEN
                FMT = "F10.1,',',"

            ELSE IF( VAL .GT. 999999. ) THEN
                FMT = "F10.2,',',"

            ELSE IF( VAL .GT. 0. .AND. VAL .LT. 1. ) THEN
                FMT = "E10.4,',',"

            ELSE
                FMT = "F10.3,',',"

            END IF

            RETURN

95020       FORMAT( 10( A, :, E12.5, :, 1X ) )

            END SUBROUTINE GET_FORMAT

C*****************  INTERNAL SUBPROGRAMS  ******************************

C----------------------------------------------------------------------
C.............  This internal subprogram determines to store a list of 
C               adjustment factos  from ADJ_FACS input file.
            SUBROUTINE READ_ADJ_FACS( F )

C.............  Subroutine arguments
            INTEGER,     INTENT( OUT ) :: F

            CHARACTER( 32 ) SEGMENT( 6 )          ! line parsing array

C.......................................................................

C.............  Loop through input files and open them
            IREC = 0
            F = 0
            DO
        
C.................  Read file names - exit if read is at end of file
                READ( ADEV, 93000, END = 30, IOSTAT=IOS ) LINE
                IREC = IREC + 1

                IF ( IOS .NE. 0 ) THEN
                    WRITE( MESG,94010 ) 
     &                  'I/O error', IOS, 
     &                  'reading adustment factor file at line', IREC
                    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                    CYCLE
                END IF

C..................  Skip blank and comment lines
                IF ( BLKORCMT( LINE ) ) CYCLE

C..................  Get line
                CALL PARSLINE( LINE, 3, SEGMENT )

                CALL UPCASE( SEGMENT( 1 ) )   ! species name
                CALL UPCASE( SEGMENT( 2 ) )   ! logical file name

C.................  Search adjustment factor for the current file
                NAM = TRIM( SEGMENT( 2 ) )
            
                L = INDEX1( NAM, NFILE, FNAME )

C.................  Skip EMF-specific header line
                IF( L <= 0 .AND. .NOT. CHKREAL( SEGMENT( 3 ) ) ) CYCLE

                F = F + 1

                ADJ_SPC( F ) = TRIM( SEGMENT( 1 ) )
                ADJ_LFN( F ) = TRIM( SEGMENT( 2 ) ) 
                ADJ_LFNSPC( F ) = TRIM( SEGMENT( 2 ) ) // '~' // 
     &                            TRIM( SEGMENT( 1 ) )
                ADJ_FACTOR( F ) = STR2REAL( SEGMENT( 3 ) )

                IF( ADJ_FACTOR( F ) < 0 ) THEN
                    MESG = 'WARNING: ' // TRIM( ADJ_SPC(F) ) // 
     &                 ' emissions from the ' //TRIM(NAM)// ' file' //
     &                 ' will be substracted due to a negative' //
     &                 ' adjustment factor' 
                    CALL M3MSG2( MESG )

                ELSE IF( ADJ_FACTOR( F ) == 0 ) THEN
                    MESG = 'WARNING: ' // TRIM( ADJ_SPC(F) ) // 
     &                 ' emissions from the ' //TRIM(NAM)// ' file' //
     &                 ' will be zero due to a zero adjustment factor' 
                    CALL M3MSG2( MESG )

                END IF

            END DO

30          CONTINUE

            RETURN

93000       FORMAT(  A )

94010       FORMAT( 10( A, :, I7, :, 1X ) )

            END SUBROUTINE READ_ADJ_FACS

C*****************  INTERNAL SUBPROGRAMS  ******************************

C----------------------------------------------------------------------
C.............  This internal subprogram determines to store a list of 
C               adjustment factos  from ADJ_FACS input file.
            SUBROUTINE READ_TAG_SPECIES( F )

C.............  Subroutine arguments
            INTEGER,     INTENT( OUT ) :: F

            CHARACTER( 32 ) SEGMENT( 6 )          ! line parsing array

C.......................................................................

            MESG = 'Enter logical name for a list of tagging species'
            TDEV = PROMPTFFILE( MESG,.TRUE.,.TRUE.,'TAG_SPECIES',PROGNAME )
            MXNTAG = GETFLINE( TDEV, 'List of tagging species' )

C.............  Allocate memory for arrays that just depend on the maximum number
C               of adjustment factors in ADJ_FACS input file.
            ALLOCATE( TAG_LFN( MXNTAG ), STAT=IOS )
            CALL CHECKMEM( IOS, 'TAG_LFN', PROGNAME )
            ALLOCATE( TAG_SPC( MXNTAG ), STAT=IOS )
            CALL CHECKMEM( IOS, 'TAG_SPC', PROGNAME )
            ALLOCATE( TAG_LFNSPC( MXNTAG ), STAT=IOS )
            CALL CHECKMEM( IOS, 'TAG_LFNSPC', PROGNAME )
            ALLOCATE( TAG_APPEND( MXNTAG ), STAT=IOS )
            CALL CHECKMEM( IOS, 'TAG_APPEND', PROGNAME )
            ALLOCATE( TAG_LFNSPCTAG( MXNTAG ), STAT=IOS )
            CALL CHECKMEM( IOS, 'TAG_LFNSPCTAG', PROGNAME )

            TAG_SPC = ' '
            TAG_LFN = ' '
            TAG_LFNSPC = ' '
            TAG_APPEND = ' '
            TAG_LFNSPCTAG = ' '

C.............  Loop through input files and open them
            IREC = 0
            F = 0
            DO
        
C.................  Read file names - exit if read is at end of file
                READ( TDEV, 93000, END = 40, IOSTAT=IOS ) LINE
                IREC = IREC + 1

                IF ( IOS .NE. 0 ) THEN
                    WRITE( MESG,94010 ) 
     &                  'I/O error', IOS, 
     &                  'reading APPEND_TAGS file at line', IREC
                    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                    CYCLE
                END IF

C..................  Skip blank and comment lines
                IF ( BLKORCMT( LINE ) ) CYCLE

C..................  Get line
                CALL PARSLINE( LINE, 3, SEGMENT )
                CALL UPCASE( SEGMENT( 1 ) )   ! logical file name
                CALL UPCASE( SEGMENT( 2 ) )   ! species name
                CALL UPCASE( SEGMENT( 3 ) )   ! tagging name

C.................  Skip EMF-specific header line
                IF( TRIM( SEGMENT( 1 ) ) == 'SECTOR' ) THEN
                IF( TRIM( SEGMENT( 2 ) ) == 'SPECIES' ) THEN
                IF( TRIM( SEGMENT( 3 ) ) == 'TAG' ) THEN
                    CYCLE
                END IF
                END IF
                END IF

C.................  store tagging species
                F = F + 1

                TAG_LFN( F ) = TRIM( SEGMENT( 1 ) ) 
                TAG_SPC( F ) = TRIM( SEGMENT( 2 ) )
                TAG_LFNSPC( F ) = TRIM( SEGMENT( 1 ) ) // '~' // 
     &                            TRIM( SEGMENT( 2 ) )
                TAG_APPEND( F ) = TRIM( SEGMENT( 3 ) )
                TAG_LFNSPCTAG(F)= TRIM( SEGMENT( 1 ) ) // '~' // 
     &                            TRIM( SEGMENT( 2 ) ) // '_' // 
     &                            TRIM( SEGMENT( 3 ) )

                L1 = LEN_TRIM( SEGMENT( 2 ) )
                L2 = LEN_TRIM( SEGMENT( 3 ) )
                L  = L1 + L2 + 1

                IF( L > 16 ) THEN
                    MESG = 'ERROR: Can not append a tag to the species '
     &                  // TRIM( TAG_SPC( F ) ) // ' from the ' // 
     &                  TRIM( NAM ) // ' file due to exceeding ' //
     &                  'max size(16-char) of tagged species name.'
                    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                END IF

            END DO

40          CONTINUE

            RETURN

93000       FORMAT(  A )

94010       FORMAT( 10( A, :, I7, :, 1X ) )

            END SUBROUTINE READ_TAG_SPECIES

        END PROGRAM MRGGRID
