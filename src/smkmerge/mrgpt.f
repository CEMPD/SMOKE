           PROGRAM MRGPT

C***********************************************************************
C  program body starts at line 
C
C  DESCRIPTION:
C    Program MRGPT reads STACK_GROUPS and INLN I/O API files and merges them
C    into a single STACK_GROUP/INLN file pair(depending on the inputs)
C
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
C    Original by G. Pouliot 11/30/2007
C    Revised by M. Omary 08/05/2010
C    Revised by H. Tran 09/25/2018
C    Revised by H. Tran 03/28/2024
C***********************************************************************
C
C Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
C                System
C File: @(#)$Id: mrggrid.f,v 1.19 2007/07/11 19:30:43 bbaek Exp $
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
C Pathname: $Source: /afs/isis/depts/cep/emc/apps/archive/smoke/smoke/src/smkmerge/mrggrid.f,v $
C Last updated: $Date: 2007/07/11 19:30:43 $ 
C

C       Updated with USE M3UTILIO by Huy Tran UNC-IE on 2026-01
C***********************************************************************
 
C...........   MODULES for public variables
C.........  This module contains the global variables for the 3-d grid

        USE M3UTILIO

        IMPLICIT NONE

        INTEGER :: NCOLS, NROWS, NLAYS, NGRID
 
C...........   INCLUDES:
        INCLUDE 'EMCNST3.EXT'
C        INCLUDE 'PARMS3.EXT'
C        INCLUDE 'IODECL3.EXT'
C        INCLUDE 'FDESC3.EXT'
      
C...........   EXTERNAL FUNCTIONS
C       CHARACTER(2)  CRLF
C       LOGICAL       ENVYN
        INTEGER       GETFLINE
C       LOGICAL       GETYN       
C       INTEGER       INDEX1
C       INTEGER       LBLANK
C       INTEGER       PROMPTFFILE
C       CHARACTER(16) PROMPTMFILE
C       INTEGER       SEC2TIME
C       INTEGER       SECSDIFF
        LOGICAL       BLKORCMT
        LOGICAL       CHKREAL
C       REAL          STR2REAL
C       CHARACTER(14) MMDDYY

C        EXTERNAL CRLF, ENVYN, GETFLINE, GETYN, INDEX1, LBLANK,
C     &           PROMPTFFILE, PROMPTMFILE, SEC2TIME, SECSDIFF
C     &           BLKORCMT, CHKREAL, STR2REAL, MMDDYY
C     
CC.........  LOCAL PARAMETERS and their descriptions:
        EXTERNAL     GETFLINE, BLKORCMT, CHKREAL, C, LOCAL, PARAMETERS, 
     &               AND, THEIR, DESCRIPTIONS

C       CHARACTER(50), PARAMETER :: 
C    &  CVSW = '$Name: SMOKEv5.1_Jul2024$' ! CVS release tag

C...........   LOCAL VARIABLES and their descriptions:

C...........   Emissions arrays
        REAL, ALLOCATABLE     :: VAR_R(:,:,:)
        REAL, ALLOCATABLE    :: STACK_PARAM_R_IN ( :,: , :)   ! 

        INTEGER, ALLOCATABLE  :: VAR_I(:,:,:)
        INTEGER, ALLOCATABLE :: STACK_PARAM_I_IN ( :, :, :)

        REAL, ALLOCATABLE    :: STACK_PARAM_R_OUT ( :,: )   ! 
        INTEGER, ALLOCATABLE :: STACK_PARAM_I_OUT ( :, :)
                
        REAL, ALLOCATABLE :: EIN ( : ,: , :)
        REAL, ALLOCATABLE :: EOUT( : ,:)
        REAL, ALLOCATABLE :: EOUT1(:, : ,:)

        REAL, ALLOCATABLE :: EMIS_TOTALS_BA(: ,:, :)  !Total emissins Before Adjustment
        REAL, ALLOCATABLE :: EMIS_TOTALS_AA(: ,:, :)  !Total emissins After Adjustment

C...........   Input file descriptors
        INTEGER,       ALLOCATABLE :: DURATA( : , : ) ! no. time steps
        INTEGER,       ALLOCATABLE :: NCOLSA( : , : ) ! no. columns
        INTEGER,       ALLOCATABLE :: NROWSA( : , :) ! no. rows
        INTEGER,       ALLOCATABLE :: NVARSA( : , :) ! no. variables
        INTEGER,       ALLOCATABLE :: SDATEA( : , :) ! start date
        INTEGER,       ALLOCATABLE :: STIMEA( : , :) ! start time
        INTEGER,       ALLOCATABLE :: NLAYSA( : , :) ! number of layers in the file
        CHARACTER(16), ALLOCATABLE :: FNAME ( : , :) ! 2-d input file names
        LOGICAL,       ALLOCATABLE :: USEFIRST(:) ! true: use first time step of file emission file only

        LOGICAL,       ALLOCATABLE :: LVOUTA( :,: , :) ! iff out var in input file
        CHARACTER(16), ALLOCATABLE :: VNAMEA( :,: , :) ! variable names
        CHARACTER(16), ALLOCATABLE :: VUNITA( :,: , :) ! variable units
        CHARACTER(80), ALLOCATABLE :: VDESCA( :,: , :) ! var descrip
        INTEGER      , ALLOCATABLE :: VTYPEA( :,: , :) ! var type
        INTEGER      , ALLOCATABLE :: NSRC(:,:)

C...........   Adujustemnt related arrays
        REAL,          ALLOCATABLE :: ADJ_FACTOR( : ) ! adjustment factors
        CHARACTER(32), ALLOCATABLE :: ADJ_LFN( : )    ! Species name
        CHARACTER(16), ALLOCATABLE :: ADJ_SPC( : )    ! logicalFileName
        CHARACTER(80)  NAME1                     ! tmp file name component
        CHARACTER(49), ALLOCATABLE :: ADJ_LFNSPC( : ) ! concatenated {logicalFileName}_{Species}
        CHARACTER(49)  LFNSPC                    ! tmp spec and file name
        CHARACTER(16)  SPCTMP                    ! tmp species name
        CHARACTER(32)  LFNTMP                    ! tmp file name
        CHARACTER(32)  LNAM                      ! tmp previous file name
        CHARACTER(16)  NAM1                      ! tmp file name

        REAL, ALLOCATABLE :: BEFORE_ADJ( : )  ! emissions before factors applied
        REAL, ALLOCATABLE :: AFTER_ADJ ( : )  ! emissions after factors applied
        REAL, ALLOCATABLE :: BEFORE_SPC( : )  ! emissions before factors applied
        REAL, ALLOCATABLE :: AFTER_SPC ( : )  ! emissions after factors applied

        REAL       :: FACS = 1.0                 ! adjustment factor 
        REAL        :: A

        INTEGER       ADJ                        ! tmp adjustment factor main index
        INTEGER       ADJ1                       ! tmp adjustment factor index 1
        INTEGER       ADJ2                       ! tmp adjustment factor index 2
        INTEGER       MXNF                       ! tmp no. of 2-d input files
        INTEGER       MXNFAC                     ! max no. of adjustment factors
        INTEGER       SDEV                       ! unit for overall QA report file
        INTEGER       EDEV
        INTEGER    :: NADJ = 0                   ! no. of adjustment factors

C...........   Intermediate output variable arrays
        INTEGER       INDXN ( MXVARS3,2 ) ! sorting index for OUTIDX
        INTEGER       OUTIDX( MXVARS3,2 ) ! index to master model species list

        CHARACTER(16) OUTNAM( MXVARS3,2 ) ! unsorted output variable names
        CHARACTER(16) VUNITU( MXVARS3,2 ) ! unsorted output variable units
        CHARACTER(80) VDESCU( MXVARS3,2 ) ! unsorted output variable descriptions
        INTEGER       VTYPEU( MXVARS3,2 ) ! unsorted output variable type

C...........   Logical names and unit numbers

        INTEGER       IDEV            ! unit for logical names list for stack_groups files
        INTEGER       JDEV            ! unit for logical names list for inline_pt_emis files
        INTEGER       LDEV            ! unit for log file
        INTEGER       RDEV            ! unit for merge report file
        INTEGER       ADEV            ! unit for logical names list for SEG
        INTEGER       ODEV            ! unit for QA report file

        CHARACTER(16) ONAME_S           ! Merged output stack groups file name
        CHARACTER(16) ONAME_E           ! merge output emission file name

C...........   Other local variables 
        INTEGER       I, C,DD, F,F1, J,J1, K, L, L1, L2, NL, V, T, S ! pointers and counters

        INTEGER       ISTART, IEND, ICT
        INTEGER       NVAR_INT                   !Number if Integer Variables
        INTEGER       NVAR_REAL                  !Number if REAL Variables
        INTEGER       J_LOOP
        INTEGER       DUMMY                      ! dummy value for use with I/O API functions
        INTEGER       EDATE                      ! ending julian date
        INTEGER       ETIME                      ! ending time HHMMSS
        INTEGER    :: G_SDATE = 0                ! start date from environment
        INTEGER    :: G_STIME = 0                ! start time from environment
        INTEGER    :: G_NSTEPS = 1               ! number of time steps from environment
        INTEGER    :: G_TSTEP = 0                ! time step from environment
        INTEGER       IOS                        ! i/o status
        INTEGER       IREC                       ! line number count
        INTEGER       JDATE                      ! iterative julian date
        INTEGER       JTIME                      ! iterative time HHMMSS
        INTEGER       LB                         ! leading blanks counter
        INTEGER       LE                         ! location of end of string
        INTEGER       MXNFIL_1                   ! max no. of stack groups files
        INTEGER       MXNFIL_2                   ! max no. of emission files
        INTEGER       MXNFIL                     ! max no. of both
        INTEGER       NFILE                      ! no. of 2-d input files
        INTEGER       NSTEPS                     ! no. of output time steps
        INTEGER       NVOUT(2)                   ! no. of output variables
        INTEGER       NVOUT_I(2)                 ! no. of integer type output variables        
        INTEGER       NVOUT_R(2)                 ! no. of real type output variables        
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
        INTEGER       VTYPE                      
        INTEGER, PARAMETER :: MAXFIP = 120       ! Maximum number of FIPS; HTmd            
C       INTEGER       NSTATES(2)                 ! no. of states   
        INTEGER, ALLOCATABLE :: NSTATES(:)       ! no. of states; HTmd
C       INTEGER       STATEFIPS(120,2)           ! State fips 
        INTEGER, ALLOCATABLE :: STATEFIPS(:,:)   ! State fips in SSCCC format; HTmd
        INTEGER       FIPS                       ! Temporay State fips   
        INTEGER, PARAMETER :: NLIST = 2          ! Number of list files to be read   
                           
        INTEGER       TMAX                       ! max time counter
        
        INTEGER       IDX_FIP, IFIP1
        INTEGER    :: MAXSTATES = 299
        
        CHARACTER(16)  FDESC                     ! tmp file description
        CHARACTER(16)  NAM                       ! tmp file name
        CHARACTER(16)  VNM                       ! tmp variable name
        CHARACTER(256) LINE                      ! input buffer
        CHARACTER(256) MESG                      ! message field
        CHARACTER(15)  RPTCOL                    ! single column in report line
        CHARACTER(300) RPTLINE                   ! line of report file
        CHARACTER(6)   FIP_CODE                  ! 

        LOGICAL    :: EFLAG   = .FALSE.   ! error flag
        LOGICAL    :: FIRST3D = .TRUE.    ! true: first 3-d file not yet input
        LOGICAL    :: LFLAG   = .FALSE.   ! true iff 3-d file input
        LOGICAL    :: TFLAG   = .FALSE.   ! true: grid didn't match
        LOGICAL       CREATE_REPORT       ! true: create report by state
        LOGICAL       MRGDIFF             ! true: merge files from different days
   
        CHARACTER(16) :: PROGNAME = 'MRGPT' ! program name
        CHARACTER(5)  :: IFIP2
C       CHARACTER(5)  :: STATEFIP2(120,NLIST)
        CHARACTER(5), ALLOCATABLE :: STATEFIP2(:,:)  ! State code in SS format; HTmd
        CHARACTER(17) :: UNDERLINE(120) 
C***********************************************************************
C   begin body of program MRGGRID
 
        LDEV = INIT3()
 
C.........  Write out copyright, version, web address, header info, and prompt
C           to continue running the program.
        CALL INITEM( LDEV, CVSW, PROGNAME )

C.........  Read names of input files and open files
        MESG = 'Enter logical name for STACK_GROUP INPUTS list'

        IDEV = PROMPTFFILE( MESG, .TRUE., .TRUE.,
     &                      'FILELIST_STACK', PROGNAME )

        MESG = 'Enter logical name for INLINEPT INPUTS list'
        JDEV = PROMPTFFILE( MESG, .TRUE., .TRUE.,
     &                      'FILELIST_INLN', PROGNAME )
     
C.........  Determine maximum number of input files in file
        MXNFIL_1 = GETFLINE( IDEV, 'List of files to merge' )

        MXNFIL_2 = GETFLINE( JDEV, 'List of files to merge' )
        
        IF (MXNFIL_1 .NE. MXNFIL_2) THEN
           MESG = 'Inconsistent number of input files '
           CALL M3MSG2( MESG )
           MESG = 'Ending program "MRGGRID".'
           CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )                            
        
        ENDIF
        MXNFIL = MXNFIL_1

C.........  Allocate memory for arrays that just depend on the maximum number
C           of files
        ALLOCATE( DURATA( MXNFIL ,NLIST), STAT=IOS )
        CALL CHECKMEM( IOS, 'DURATA', PROGNAME )
        ALLOCATE( NCOLSA( MXNFIL ,NLIST), STAT=IOS )
        CALL CHECKMEM( IOS, 'NCOLSA', PROGNAME )
        ALLOCATE( NROWSA( MXNFIL ,NLIST), STAT=IOS )
        CALL CHECKMEM( IOS, 'NROWSA', PROGNAME )
        ALLOCATE( NLAYSA( MXNFIL ,NLIST), STAT=IOS )
        CALL CHECKMEM( IOS, 'NLAYSA', PROGNAME )
        ALLOCATE( NVARSA( MXNFIL ,NLIST), STAT=IOS )
        CALL CHECKMEM( IOS, 'NVARSA', PROGNAME )
        ALLOCATE( SDATEA( MXNFIL ,NLIST), STAT=IOS )
        CALL CHECKMEM( IOS, 'SDATEA', PROGNAME )
        ALLOCATE( STIMEA( MXNFIL ,NLIST), STAT=IOS )
        CALL CHECKMEM( IOS, 'STIMEA', PROGNAME )
        ALLOCATE( FNAME( MXNFIL,NLIST ), STAT=IOS )
        CALL CHECKMEM( IOS, 'FNAME', PROGNAME )
               
        ALLOCATE( LVOUTA( MXVARS3,MXNFIL,NLIST ), STAT=IOS )
        CALL CHECKMEM( IOS, 'LVOUTA', PROGNAME )
        ALLOCATE( VNAMEA( MXVARS3,MXNFIL,NLIST ), STAT=IOS )
        CALL CHECKMEM( IOS, 'VNAMEA', PROGNAME )
        ALLOCATE( VUNITA( MXVARS3,MXNFIL,NLIST ), STAT=IOS )
        CALL CHECKMEM( IOS, 'VUNITA', PROGNAME )
        ALLOCATE( VDESCA( MXVARS3,MXNFIL,NLIST ), STAT=IOS )
        CALL CHECKMEM( IOS, 'VDESCA', PROGNAME )
        ALLOCATE( VTYPEA( MXVARS3,MXNFIL,NLIST ), STAT=IOS )
        CALL CHECKMEM( IOS, 'VTYPEA', PROGNAME )        

        ALLOCATE( NSRC( NLIST,MXNFIL+1), STAT=IOS )
        CALL CHECKMEM( IOS, 'NSRC', PROGNAME )
 
C.........  Loop through input files (1) and open them
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

            IF( LINE .EQ. ' ' ) CYCLE

            F = F + 1

            IF( F .LE. MXNFIL ) THEN

                LB = LBLANK ( LINE )
                LE = LEN_TRIM( LINE )
                FNAME( F,1 ) = LINE( LB+1:LE )

                IF ( .NOT. OPEN3( FNAME(F,1), FSREAD3, PROGNAME )) THEN
 
                    MESG = 'Could not open file "' //
     &                 FNAME( F,1 )( 1 : LEN_TRIM( FNAME(F,1) ) )// '".'
                    CALL M3MSG2( MESG )
                    MESG = 'Ending program "MRGGRID".'
                    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                END IF      !  if open3() failed
                
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
      
C.........  Loop through input files (2) and open them
        F = 0
        IREC = 0
        DO

C.............  Read file names - exit if read is at end of file
            READ( JDEV, 93000, END=28, IOSTAT=IOS ) LINE
            IREC = IREC + 1

            IF ( IOS .NE. 0 ) THEN
               EFLAG = .TRUE.
                WRITE( MESG,94010 ) 
     &              'I/O error', IOS, 
     &              'reading file list at line', IREC
                CALL M3MESG( MESG )
                CYCLE
            END IF

            IF( LINE .EQ. ' ' ) CYCLE

            F = F + 1

            IF( F .LE. MXNFIL ) THEN

                LB = LBLANK ( LINE )
                LE = LEN_TRIM( LINE )
                FNAME( F,2 ) = LINE( LB+1:LE )

                IF ( .NOT. OPEN3( FNAME(F,2), FSREAD3, PROGNAME )) THEN
 
                    MESG = 'Could not open file "' //
     &                 FNAME( F,2 )( 1 : LEN_TRIM( FNAME(F,2) ) )// '".'
                    CALL M3MSG2( MESG )
                    MESG = 'Ending program "MRGGRID".'
                    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

                END IF      !  if open3() failed
                
            END IF

        END DO
28      CONTINUE

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
      WRITE(*,*)'   '
      WRITE(*,*)'########  Finished Reading the List Files and ########'
      WRITE(*,*)'######## Opennig the input Files              ########'
      WRITE(*,*)'######################################################'      

c----HTmd: Allocate memory for some variables with the defined NFILE----
      ALLOCATE( NSTATES( NFILE ) )
      CALL CHECKMEM( IOS, 'STATES', PROGNAME )
      ALLOCATE( STATEFIPS( MAXFIP, NFILE ) )
      CALL CHECKMEM( IOS, 'STATEFIPS', PROGNAME )
      ALLOCATE( STATEFIP2( MAXFIP, NFILE ) )
      CALL CHECKMEM( IOS, 'STATEFIP2', PROGNAME )
      ALLOCATE ( USEFIRST ( NFILE ) ) 
      CALL CHECKMEM( IOS, 'USEFIRST', PROGNAME )
c----End HTmd

c#########################################################################
C.........  Get file descriptions and store for all input files
C.........  Loop through 2 sets of input files
        NLAYS = 1
        DO S = 1,NLIST !2
         DO F = 1, NFILE

            NAM = FNAME( F,S )
            
            IF ( .NOT. DESC3( NAM ) ) THEN
                MESG = 'Could not get description of file "'  //
     &                  NAM( 1:LEN_TRIM( NAM ) ) // '"'
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            ELSE
                NROWSA( F ,S) = NROWS3D
                NCOLSA( F ,S) = NCOLS3D
                NLAYSA( F ,S) = NLAYS3D
                NVARSA( F ,S) = NVARS3D
                SDATEA( F ,S) = SDATE3D
                STIMEA( F ,S) = STIME3D
                DURATA( F ,S) = MXREC3D
                
                IF( F == 1 ) TSTEP = TSTEP3D
                
                DO V = 1, NVARS3D
                    VNAMEA( V,F,S ) = VNAME3D( V )
                    VUNITA( V,F,S ) = UNITS3D( V )
                    VDESCA( V,F,S ) = VDESC3D( V )
                    VTYPEA( V,F,S ) = VTYPE3D( V )
                END DO
            END IF

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

         END DO      ! End loop through files
        END DO     ! s_loop

c.............  Check if the # of recrords (NROWS) are tha same in the Stack group file and the Emmissions file

       DO F = 1, NFILE
         IF( NROWSA( F ,1) .NE. NROWSA( F ,2) ) THEN
             WRITE(*,*)'Number of sources in the stack group file ',FNAME( F,1 ),NROWSA( F ,1)
             WRITE(*,*)'Number of sources in the emissions file   ',FNAME( F,2 ),NROWSA( F ,2)
 
            MESG = 'Inconsistent number of sources between:' //CRLF() //
     &              BLANK10 //'stack group file '//TRIM(FNAME( F,1 ))//' '//
     &              'and Emissions file '//TRIM(FNAME( F,2 ))//CRLF() //
     &              BLANK10 //'The order of the input files in the FILELIST_STACK should be'//
     &              CRLF() //BLANK10 //'the same as in the FILELIST_INLN  '               
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
         END IF
       END DO
c#####################################################################
c#####################################################################
c#####################################################################
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

c#####################################################################
c#####################################################################
c#####################################################################

C.........  Give error message and end program unsuccessfully
        IF( EFLAG ) THEN
            MESG = 'Inconsistent time step, grid, or layers ' //
     &              'among the files!'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

C.........  Deterimine output date, time, and number of time steps
        SDATE = G_SDATE
        STIME = G_STIME
        NSTEPS = G_NSTEPS
        CALL SETOUTDATE( SDATE, STIME, NSTEPS, NFILE, SDATEA(1,2),
     &                   STIMEA(1,2), DURATA(1,2), FNAME(1,2), MRGDIFF,
     &                   USEFIRST )
        WRITE(*,*)'======================================================'

c#####################################################################
c#####################################################################
c#####################################################################

      DO S = 1,NLIST   ! THIS loop is to read in the input files only

C.........  Build master output variables list

          NVOUT(S) = 0
          NSRC(S,1:NFILE+1) = 0
          NVOUT_R(S) = 0
          NVOUT_I(S) = 0

C.........  Loop through input files and build an output variable list
          WRITE(*,*)'===>>>       List#(S)     File#(F)       Records# '

         DO F = 1, NFILE

	   NSRC(S,F) = NSRC(S,F) + NCOLSA(F,S)*NROWSA(F,S)
           WRITE(*,*)'===>>> ',S,F,NSRC(S,F),'  ',FNAME( F,S )	
C.............  Loop through variables in the files
            DO V = 1, NVARSA( F,S )

                VNM = VNAMEA( V,F,S )       

C.................  Look for variable name in output list
                K = INDEX1( VNM, NVOUT(S), OUTNAM(1,S)  )  ! look in output list

C.................  If its not in the output list, add it
                IF( K .LE. 0 ) THEN
                    NVOUT(S) = NVOUT(S) + 1
                    INDXN ( NVOUT(S),S ) = NVOUT(S)
                    OUTNAM( NVOUT(S),S ) = VNM
                    VDESCU( NVOUT(S),S ) = VDESCA( V,F,S )
                    VUNITU( NVOUT(S),S ) = VUNITA( V,F,S )
                    VTYPEU( NVOUT(S),S ) = VTYPEA( V,F,S )

                    IF (VTYPEA(V,F,S) .eq. M3REAL) THEN
                       NVOUT_R(S) = NVOUT_R(S) + 1
                    END IF

                    IF (VTYPEA(V,F,S) .eq. M3INT) THEN
                       NVOUT_I(S) = NVOUT_I(S) + 1
                    END IF

C.................  If variable is in the output list, check the units
                ELSE	
                    IF (VNM(2:5) .ne. 'LOCA') THEN
                      ! don't check units on XLOCA and YLOCA

                    IF ( VUNITA( V,F,S ) .NE. VUNITU( K,S ) ) THEN
                        EFLAG = .TRUE.
                        L  = LEN_TRIM( VNM )
                        L1 = LEN_TRIM( VUNITA( V,F,S ) )
                        L2 = LEN_TRIM( VUNITU( K,S )   )
                        WRITE( MESG,94010 ) 'ERROR: Variable "' //
     &                         VNM( 1:L ) // '" in file', F,
     &                         'has units "'// VUNITA( V,F,S )( 1:L1 ) //
     &                         '"' // CRLF() // BLANK10 //
     &                         'that are inconsistent with a '//
     &                         'previous file that had units "' //
     &                         VUNITU(K,S)(1:L2)//'" for this variable'
                        CALL M3MSG2( MESG )

                    END IF  ! End check of units
                    END IF

                END IF      ! End variable in output list already or not

            END DO          ! End loop through variables in this file
         END DO              ! End loop through files.

        NSRC(S,NFILE+1) = SUM(NSRC(S,1:NFILE))
       
C.........  Give error message and end program unsuccessfully
        IF( EFLAG ) THEN
            MESG = 'Inconsistent units for common variables among '//
     &             'the files!'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF


C.........  Sort output variables into alphabetical order
        CALL SORTIC( NVOUT(S), INDXN(1,S), OUTNAM(1,S) )

        DO V = 1, NVOUT(S)
            VNM = OUTNAM( INDXN( V,S ),S )

            DO F = 1, NFILE
                LVOUTA( V,F,S ) = .FALSE.

                J = INDEX1( VNM, NVARSA( F,S ), VNAMEA( 1,F,S ) )
                IF( J .GT. 0 ) LVOUTA( V,F,S ) = .TRUE.
            END DO 
        END DO
        
         WRITE (*,*) S, 'TOTAL NUMBER OF SOURCES = ', NSRC(S,NFILE+1)

         IF (S .EQ. 1) THEN

            ALLOCATE (STACK_PARAM_I_IN(NSRC(S,NFILE+1),NVOUT_I(S),
     &                    NFILE),STAT=IOS)
            CALL CHECKMEM( IOS, 'STACK_PARAM_I_IN', PROGNAME ) 

            ALLOCATE (STACK_PARAM_R_IN(NVOUT_R(S),NSRC(S,NFILE+1),
     &                    NFILE),STAT=IOS)
            CALL CHECKMEM( IOS, 'STACK_PARAM_R_IN', PROGNAME )
            
            ALLOCATE (STACK_PARAM_I_OUT(NSRC(S,NFILE+1),NVOUT_I(S)),
     &                STAT=IOS)
            CALL CHECKMEM( IOS, 'STACK_PARAM_I_OUT', PROGNAME )       

            ALLOCATE (STACK_PARAM_R_OUT(NSRC(S,NFILE+1),NVOUT_R(S)),
     &                STAT=IOS)
            CALL CHECKMEM( IOS, 'STACK_PARAM_R_OUT', PROGNAME )
                         
         ENDIF

         IF ( S .EQ. 2) THEN

             ALLOCATE (EIN( NSRC(S,NFILE+1), NVOUT(S),NFILE),STAT=IOS)
             CALL CHECKMEM( IOS, 'EIN', PROGNAME )
         
            ALLOCATE (EOUT( NSRC(S,NFILE+1), NVOUT(S)),STAT=IOS)
            CALL CHECKMEM( IOS, 'EOUT', PROGNAME )
    
            ALLOCATE (EOUT1( NSTEPS,NSRC(S,NFILE+1), NVOUT(S)),STAT=IOS)
            CALL CHECKMEM( IOS, 'EOUT1', PROGNAME )
           
            ALLOCATE (EMIS_TOTALS_BA(MAXSTATES,NVOUT(S),NFILE),STAT=IOS)
            CALL CHECKMEM( IOS, 'EMIS_TOTALS_BA', PROGNAME )
            ALLOCATE (EMIS_TOTALS_AA(MAXSTATES,NVOUT(S),NFILE),STAT=IOS)
            CALL CHECKMEM( IOS, 'EMIS_TOTALS_AA', PROGNAME )
            
            EMIS_TOTALS_BA(1:MAXSTATES,1:NVOUT(S),1:NFILE) = 0.0
	    EMIS_TOTALS_AA(1:MAXSTATES,1:NVOUT(S),1:NFILE) = 0.0
         ENDIF
	
        WRITE(*,*)'======================================================'
    
C.........  Loop through hours
          JDATE  = SDATE
          JTIME  = STIME
         IF (S .EQ. 1) THEN
            TMAX = 1
          ELSE
            TMAX = NSTEPS 
          ENDIF
          NSTATES = 0 
	  STATEFIPS = 0
          STATEFIP2 = ' '
         DO T = 1,TMAX

C.............  Loop through species
            NVAR_INT = 0
            NVAR_REAL = 0
          DO V = 1, NVOUT( S )
      
           VNM = OUTNAM( INDXN( V,S ),S ) 
	   
C.................  Output array
!                EOUT = 0.   ! array
                ICT = 0

           DO F = 1, NFILE

           NAM = FNAME(F,S)
           IF ( ALLOCATED(VAR_I) ) DEALLOCATE(VAR_I) !HTmd
           IF ( ALLOCATED(VAR_R) ) DEALLOCATE(VAR_R) !HTmd
           ALLOCATE(VAR_I  (NCOLSA(F,S),NROWSA(F,S),NLAYSA(F,S)))
	   ALLOCATE(VAR_R  (NCOLSA(F,S),NROWSA(F,S),NLAYSA(F,S)))
           VAR_R = 0.0
	   VAR_I = 0
	 
C.....................  Set read date
		        IF(S .EQ. 1) THEN
			RDATE = SDATEA( F, S )
			ELSE
                        RDATE =  JDATE
			END IF

C.....................  If file has species, read (do this for all files)...
          IF( LVOUTA( V,F,S ) ) THEN   !If Pollutant is in the file

	   VTYPE = VTYPEA(INDXN( V,S ),F,S)

	     IF (VTYPE .EQ.M3REAL) THEN
	       
	       IF(.NOT. READ3( NAM, VNM, ALLAYS3 , RDATE, 
     &            JTIME, VAR_R) )THEN
                  MESG = 'Could not read "' // VNM //
     &                    '" from file "' //
     &                     NAM( 1:LEN_TRIM( NAM ) )// '".'
                  CALL M3EXIT( PROGNAME, RDATE, JTIME, 
     &                                       MESG, 2 )
               END IF
	       
                    IF (S .EQ. 1) THEN
                     IF ( F .EQ. 1) NVAR_REAL = NVAR_REAL + 1
                      ISTART = ICT + 1
                      IEND   = ISTART + NSRC(S,F) -1
		      STACK_PARAM_R_OUT(ISTART:IEND,NVAR_REAL) = VAR_R(1,1:NSRC(S,F),1)	
		    END IF
		  
                    IF ( (S .EQ. 2)) THEN
                        ISTART = ICT + 1
                        IEND   = ISTART + NSRC(S,F) - 1

C.............  Search for adj factor species in the logical file
                      IF(ADEV .GT. 0) THEN
                        DO J = 1, NADJ
                          LFNTMP = ADJ_LFN( J ) ! retriev logical file from ADJ_FACS

                          IF( LFNTMP == NAM ) THEN
                             SPCTMP = ADJ_SPC( J ) ! retrieve spcieces name from ADJ_FACS

                           K = INDEX1( SPCTMP, NVARSA(F,S), VNAMEA(1,F,S) ) !VNAMEA( V,F,S ) )

                           IF( K <= 0 ) THEN
                               EFLAG = .TRUE.
                               MESG = 'WARNING: The species '//TRIM(ADJ_SPC(J))// 
     &                         ' you want to adjust is not available in the '//
     &                         TRIM(NAM) // ' file'
                               CALL M3MSG2( MESG )
                           END IF
                          END IF
                        END DO   !END NADJ

C.........  Warning missing emissions logical file names from the ADJ_FACS list 
                         LNAM = ' '
                         DO I = 1, NADJ
                            NAM1 = ADJ_LFN( I )
                            L = INDEX1( NAM1, NFILE, FNAME(1,2) )
                           IF( L <= 0 ) THEN
                             MESG = 'WARNING: The logical file '//TRIM(NAM1) //
     &                       ' in the adjustment file (ADJ_FACS) is not '// 
     &                       'found in the FILELIST_INLN on DATE : '//
     &                       MMDDYY(SDATE)
                              IF( LNAM /= NAM1 ) THEN
                                 CALL M3MSG2( MESG )
                                 LNAM = NAM1
                              END IF
                           END IF
                         END DO
C.....................  Search adjustment factor for the current emissions file
                         LFNSPC = TRIM( NAM ) // '~' // TRIM( VNM )

C.....................  Assign adjustment factor for the current species

                         ADJ = INDEX1( LFNSPC, NADJ, ADJ_LFNSPC )
                         IF( ADJ > 0 ) THEN
                             FACS = ADJ_FACTOR( ADJ )

                             WRITE( MESG,93011 )'Apply adjustment factor' ,
     &                           FACS, ' to the '  // TRIM( VNM ) //
     &                           ' species from the '//TRIM( NAM )// ' file at time ',JTIME
                             CALL M3MSG2( MESG )
                         ELSE
                             FACS = 1.0
                         END IF
			    
	              ELSE		    
		        FACS = 1.0  
                      END IF   !ADEV > 0   

		        DO J = 1,NSRC(S,F)
			   IFIP1 =STACK_PARAM_I_IN(J,IDX_FIP,F)/1000
                           WRITE(IFIP2,'(I3)')IFIP1
C.................  Look for state fips in output list
                K = INDEX1( IFIP2, NSTATES(F), STATEFIP2(1,F)  )  ! look in output list

C.................  If its not in the output list, add it
                         IF( K .LE. 0 ) THEN
                           NSTATES(F) = NSTATES(F) + 1
		           STATEFIP2(NSTATES(F),F ) = IFIP2		    
                           STATEFIPS( NSTATES(F),F ) = IFIP1
                         END IF			   

                         EMIS_TOTALS_BA(IFIP1,V,F) =
     &                        EMIS_TOTALS_BA(IFIP1,V,F) + VAR_R(1,J,1)
                        END DO


cxx                        EOUT(ISTART:IEND,V) = VAR_R(1,1:NSRC(S,F),1)*FACS
                        EOUT1(T,ISTART:IEND,V) = VAR_R(1,1:NSRC(S,F),1)*FACS
			
		         DO J = 1,NSRC(S,F)
			   IFIP1 =STACK_PARAM_I_IN(J,IDX_FIP,F)/1000
                           EMIS_TOTALS_AA(IFIP1,V,F) =
     &                           EMIS_TOTALS_AA(IFIP1,V,F) + VAR_R(1,J,1)*FACS
                         END DO

                    END IF  ! S = 2

	     ELSEIF (VTYPE .EQ. M3INT) THEN
	        IF(.NOT. READ3( NAM, VNM, ALLAYS3 , RDATE, 
     &            JTIME, VAR_I) )THEN
                  MESG = 'Could not read.... "' // VNM //
     &                   '" from file "' //
     &                    NAM( 1:LEN_TRIM( NAM ) )// '".'
                  CALL M3EXIT( PROGNAME, RDATE, JTIME, 
     &                                       MESG, 2 )
                END IF	  

                   IF (S .EQ. 1) THEN
                     IF ( F .EQ. 1) NVAR_INT = NVAR_INT + 1
                        ISTART = ICT + 1
                        IEND   = ISTART + NSRC(S,F) -1
                        IF (VNM .NE. 'ISTACK') THEN        
                            IF (VNM .EQ. 'IFIP') IDX_FIP = NVAR_INT
			STACK_PARAM_I_IN(1:NSRC(S,F),NVAR_INT ,F) = 
     &                                           VAR_I(1,1:NSRC(S,F),1)
        
                        STACK_PARAM_I_OUT(ISTART:IEND,NVAR_INT)  = VAR_I(1,1:NSRC(S,F),1)
                        ELSE
                           DO J_LOOP = ISTART, IEND
                           STACK_PARAM_I_OUT(J_LOOP,NVAR_INT) = J_LOOP
                           ENDDO
                        ENDIF
                   ENDIF

	     END IF ! IF (VTYPE
          ELSE
	    IF (T .EQ. 1) THEN	    
	    WRITE(*,*)'Variable ...',VNM,' is not in the file ..',NAM
            END IF   
          END IF !LVOUTA  If Pollutant is in the file 

C.....................  Build report line if needed
	   
	     ICT = ICT + NSRC(S,F) 
        END DO       ! End loop through files

        END DO      ! End loop through Variable

C.............  Write this time step to report
	      IF(S .EQ.2) THEN
                  CALL NEXTIME( JDATE, JTIME, TSTEP3D)
             END IF

       END DO    ! End loop throught T   

      END DO     ! s_loop

c############################################################################
c############################################################################
c############################################################################

       DO S = 1,NLIST  ! this S_loop is for writing the output files and doing reports

C.........  Set up for opening output file...
C.........  Get grid information
        IF( .NOT. DESC3( FNAME( 1,S ) ) ) THEN
            MESG = 'Could not get description of file "'  //
     &              FNAME( 1,S )( 1:LEN_TRIM( FNAME(1,S) ) ) // '"'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        ENDIF

        SDATE3D = SDATE
        STIME3D = STIME
        IF (S .EQ. 1) TSTEP3D = 0
        NCOLS3D = 1
        NROWS3D = NSRC(S,NFILE+1)
        
        NVARS3D = NVOUT(S)

C.........  Set up layer structure for output file
        NLAYS3D = 1
!        VGTOP3D = VGTOP
!        VGTYP3D = VGTYP
        VGLVS3D = 0.     ! initialize array


C........  Update variable names in sorted order, and also 
C........  set up logical arrays for which files have which species
        DO V = 1, NVOUT(S)
            VNM = OUTNAM( INDXN( V,S ),S )

            VNAME3D( V ) = VNM                  ! store sorted output vars, etc.
            VDESC3D( V ) = VDESCU( INDXN( V,S ),S )
            UNITS3D( V ) = VUNITU( INDXN( V,S ),S )
            VTYPE3D( V ) = VTYPEU( INDXN( V,S ),S )
        END DO

        IF (S .EQ. 1) THEN
C.........  Prompt for and open output file
        ONAME_S = PROMPTMFILE( 
     &          'Enter logical name for MERGED STACK GROUPS file',
     &          FSUNKN3, 'OUTFILE_S', PROGNAME )
        ENDIF
        IF (S .EQ. 2) THEN
C.........  Prompt for and open output file
        ONAME_E = PROMPTMFILE( 
     &          'Enter logical name for MERGED PT EMISSION FILE file',
     &          FSUNKN3, 'OUTFILE_E', PROGNAME )
        ENDIF

C.........  Loop through hours
         IF ( S .EQ. 1) THEN
            TMAX = 1
         ELSE
            TMAX = NSTEPS
         ENDIF

            JDATE = SDATE
            JTIME = STIME            
              
        DO T = 1, TMAX

C.............  Loop through species
            NVAR_REAL = 0
            NVAR_INT = 0
            DO V = 1, NVOUT(S)

                VNM = OUTNAM( INDXN( V,S ),S )

C.................  Output array
                 IF (S .EQ. 2 ) THEN
		 EOUT(1:NSRC(S,NFILE+1),V) = EOUT1(T,1:NSRC(S,NFILE+1),V)
                   IF( .NOT. WRITE3( ONAME_E, VNM, JDATE, JTIME,
     &               EOUT(1,V) )) THEN		   
                     MESG = 'Could not write "'// VNM// '" to file "'// 
     &                      ONAME_E( 1:LEN_TRIM( ONAME_E ) ) // '".'
                     CALL M3EXIT( PROGNAME, JDATE, JTIME, MESG, 2 )
                   END IF
                 END IF

                 IF ((S .EQ. 1) .AND. (VTYPE3D(V) .EQ. M3REAL)) THEN
                   NVAR_REAL = NVAR_REAL + 1
                    IF( .NOT. WRITE3( ONAME_S, VNM, JDATE, JTIME, 
     &                 STACK_PARAM_R_OUT(1,NVAR_REAL) )) THEN
                       MESG = 'Could not write "'// VNM// '" to file "'// 
     &                      ONAME_S( 1:LEN_TRIM( ONAME_S ) ) // '".'
     &                        
                    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
		   END IF
                 END IF
 
                 IF ((S .EQ. 1) .AND. (VTYPE3D(V) .EQ. M3INT)) THEN
                  NVAR_INT = NVAR_INT + 1
                   IF( .NOT. WRITE3( ONAME_S, VNM, JDATE, JTIME, 
     &               STACK_PARAM_I_OUT(1,NVAR_INT) )) THEN
                     MESG = 'Could not write "'// VNM// '" to file "'// 
     &                      ONAME_S( 1:LEN_TRIM( ONAME_S ) ) // '".'
                     CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                   END IF
                 ENDIF
            END DO   ! loop through variables
            
C.............  Write this time step to report

            IF (S .EQ. 2) CALL NEXTIME( JDATE, JTIME, TSTEP )
       END DO       ! loop through timesteps

       END DO  ! end of S loop

c.... writing emissions report before adjustemnt
        MESG = 'Create report for total emissions by state'
        CREATE_REPORT = ENVYN( 'CREATE_REPORT', MESG, .FALSE., IOS )
	
C.........  Create report file

        IF( CREATE_REPORT ) THEN
            EDEV = PROMPTFFILE(
     &             'Enter logical name for the MRGGRID REPORT file',
     &             .FALSE., .TRUE., 'REPORT_BY_STATE', PROGNAME ) 
     
	        UNDERLINE = '-----------------'
		
	         DO F = 1, NFILE
		    WRITE(EDEV,93000)'Emissions Before Adjustment'
		    WRITE(EDEV,*)'Emissions from ......',FNAME(F,2)
                    WRITE(EDEV,950)'STATE FIP',(OUTNAM( INDXN( V,2 ),2 ),V =1,NVOUT(2))		
                    WRITE(EDEV,950)'         ',(VUNITU( INDXN( V,2 ),2 ),V =1,NVOUT(2))
		    WRITE(EDEV,960)'-----------',(UNDERLINE(K),K=1,NVOUT(2))
		  DO J = 1,NSTATES(F)
		    I = STATEFIPS(J,F)
		    FIPS = STATEFIPS(J,F)*1000
		    WRITE(FIP_CODE,'(I6)')FIPS
		      IF (FIPS .LT. 9999) THEN
		        WRITE(FIP_CODE,'(I4)')FIPS
			FIP_CODE = '00'//TRIM(FIP_CODE)
		      ELSE IF (FIPS .GT. 9999 .AND. FIPS .LT. 99999) THEN
		        WRITE(FIP_CODE,'(I5)')FIPS
			FIP_CODE = '0'//TRIM(FIP_CODE)
		      END IF
		    WRITE(EDEV,900)TRIM(FIP_CODE),(EMIS_TOTALS_BA(I,V,F),V=1,NVOUT(2))
                  END DO
		    WRITE(EDEV,93000)'  '
		 END DO

c.... writing emissions report after adjustemnt, if there is any

               IF(ADEV .GT. 0) THEN
		    WRITE(EDEV,93000)'  '
		    WRITE(EDEV,93000)'  '
		    WRITE(EDEV,93000)'Emissions After Adjustment'
	  
		 DO F = 1, NFILE
		   DO K =  1,NADJ
		      J1 = 0
                     IF(ADJ_LFN( K )  == FNAME(F,2) ) THEN
		       J1 = J1 + 1
                      WRITE( MESG,93012 )'Adjustment factor',ADJ_FACTOR(K),
     &                   ' applied to the '  // TRIM( ADJ_SPC(K)) //
     &                           ' specie from the '//TRIM(FNAME(F,2)) // ' file'
cxx                             CALL M3MSG2( MESG )
		      WRITE(EDEV,93000)MESG
		     END IF
                   END DO   !END NADJ
		     IF(J1 .EQ. 0) THEN
		     WRITE(EDEV,93000)'No Adjustment factors for '//FNAME(F,2)
		     END IF
		    		
		    WRITE(EDEV,*)'Emissions from ......',FNAME(F,2)
                    WRITE(EDEV,950)'STATE FIP',(OUTNAM( INDXN( V,2 ),2 ),V =1,NVOUT(2))		
                    WRITE(EDEV,950)'         ',(VUNITU( INDXN( V,2 ),2 ),V =1,NVOUT(2))
		    WRITE(EDEV,960)'-----------',(UNDERLINE(K),K=1,NVOUT(2))
		  DO J = 1,NSTATES(F)
		    I = STATEFIPS(J,F)
		    FIPS = STATEFIPS(J,F)*1000
		    WRITE(FIP_CODE,'(I6)')FIPS
		      IF (FIPS .LT. 9999) THEN
		        WRITE(FIP_CODE,'(I4)')FIPS
			FIP_CODE = '00'//TRIM(FIP_CODE)
		      ELSE IF (FIPS .GT. 9999 .AND. FIPS .LT. 99999) THEN
		        WRITE(FIP_CODE,'(I5)')FIPS
			FIP_CODE = '0'//TRIM(FIP_CODE)
		      END IF
		    WRITE(EDEV,900)FIP_CODE,(EMIS_TOTALS_AA(I,V,F),V=1,NVOUT(2))
                  END DO
		    WRITE(EDEV,93000)'  '
		 END DO
               ELSE
		    WRITE(EDEV,93000)'  '
		    WRITE(EDEV,93000)'  '
		    WRITE(EDEV,93000)'  No Adjustment Factors were applied'
               END IF

        ELSE
	WRITE(*,*)'########################################'
	WRITE(*,*)'No Report was created ....'
	WRITE(*,*)'########################################'
        END IF


C......... Normal Completion
        CALL M3EXIT( PROGNAME, 0, 0, ' ', 0)

C******************  FORMAT  STATEMENTS   ******************************

C...........   Informational (LOG) message formats... 92xxx
900     FORMAT(A10,'|',100(F16.5,'|'))
950     FORMAT(A10,'|',100(A16,'|'))
960     FORMAT(A11,100A17)
cx901     FORMAT(A4,3I6,3F13.5)
cx90000   FORMAT(A20,I4,2A16,I10)

92000   FORMAT( 5X, A )
 
C...........   Formatted file I/O formats............ 93xxx

93000   FORMAT(  A )

93010   FORMAT( A15 )

93020   FORMAT( I15 )

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I7, :, 1X ) )

93011   FORMAT(  A, F7.3, A,I8 )
93012   FORMAT(  A, F8.5, A )
94020   FORMAT( A, :, I3, :, 1X, 10 ( A, :, F8.5, :, 1X ) )

cxx	           END PROGRAM MRGPT
	   
C*****************  INTERNAL SUBPROGRAMS  ******************************
        CONTAINS
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
!                IF ( BLKORCMT( LINE ) ) CYCLE
                 L1 = LEN_TRIM(LINE)             !MO add
		 IF(L1 .EQ. 0) CYCLE             !MO add

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







	           END PROGRAM MRGPT
