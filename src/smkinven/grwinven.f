
         PROGRAM GRWINVEN

C***************************************************************************
C  program body starts at line 168
C
C  DESCRIPTION:
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C
C***************************************************************************
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
C***************************************************************************

C...........   MODULES for public variables
C...........   This module is the inventory arrays
        USE MODSOURC, ONLY: INVYR, FUGHGT, FUGWID, FUGLEN, FUGANG

C...........   This module contains the information about the source category
        USE MODINFO, ONLY: CATEGORY, CRL, CATDESC, BYEAR, 
     &                     NIPOL, NPPOL, NIACT, NPACT, NIPPA,
     &                     NSRC, NMAP, MAPNAM, MAPFIL, EANAM, ACTVTY

        USE MODFILESET

        IMPLICIT NONE

C...........   INCLUDES:

        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
        INCLUDE 'SETDECL.EXT'   !  FileSetAPI variables and functions

C...........   EXTERNAL FUNCTIONS and their descriptions:
        
        CHARACTER(2)            CRLF
        INTEGER                 ENVINT
        LOGICAL                 ENVYN
        INTEGER                 GETIFDSC
        INTEGER                 INDEX1
        INTEGER                 PROMPTFFILE
        INTEGER                 STR2INT

        EXTERNAL CRLF, ENVINT, ENVYN, GETIFDSC,
     &           INDEX1, PROMPTFFILE, STR2INT

C...........  LOCAL PARAMETERS and their descriptions:

        CHARACTER(50), PARAMETER :: 
     &  CVSW = '$Name SMOKEv4.8_Jun2020$' ! CVS release tag

C...........  LOCAL VARIABLES and their descriptions:
C...........  Local parameters
        INTEGER, PARAMETER :: MXVARS = 500     ! max no of pollutants in control/projection matrices
 
C.........  Array that contains the names of the inventory variables needed 
C           for this program
        CHARACTER(IOVLEN3) IVARNAMS( MXINVARR )

C...........  ALLOCATABLE VARIABLES and their decriptions...

C...........  Variables for reading and writing inventory pollutants
        INTEGER, ALLOCATABLE :: SRCID  ( : ) ! source numbers
        INTEGER, ALLOCATABLE :: INVYR_BASE( : ) ! input inventory year
        REAL   , ALLOCATABLE :: DATAVAR( :,: ) ! emissions and associated vars
        REAL   , ALLOCATABLE :: SRCDAT ( :,: ) ! emissions and associated vars

C...........  Variables for each control/projection matrix
        INTEGER, ALLOCATABLE :: CCNTRA ( : ) ! unsorted original position
        INTEGER, ALLOCATABLE :: CINDXA ( : ) ! sorting index
        INTEGER, ALLOCATABLE :: CTYPEA ( : ) ! unsorted control matrix type code
        INTEGER, ALLOCATABLE :: CTYPE  ( : ) ! sorted control matrix type code
        INTEGER, ALLOCATABLE :: IDXALL ( : ) ! index of matrixes with "All"
        INTEGER, ALLOCATABLE :: NCPVARS( : ) ! sorted no. variables

        CHARACTER(NAMLEN3), ALLOCATABLE :: CNAMEA ( : )  ! unsort mtx names
        CHARACTER(NAMLEN3), ALLOCATABLE :: CNAME  ( : )  ! sorted mtx names
        CHARACTER(NAMLEN3), ALLOCATABLE :: CPVNAMS( :,: )! mtx var names

C...........  Control matrix factors
        REAL, ALLOCATABLE :: CFAC   ( : )  ! by-variable factors
        REAL, ALLOCATABLE :: CFACALL( :,: )  ! all-variable factors
        REAL, ALLOCATABLE :: CEFF   ( : )  ! by-variable control efficiency
        REAL, ALLOCATABLE :: REFF   ( : )  ! by-variable rule effectivness
        REAL, ALLOCATABLE :: RPEN   ( : )  ! by-variable rule penetration

C...........  Temporary IDA pollutant files unit nos.
        INTEGER, ALLOCATABLE :: TDEV( : )

C...........  FIXED DIMENSION VARIABLES and their descriptions...

C...........  Inventory file variable names
        CHARACTER(IOVLEN3), ALLOCATABLE :: IVNAMES( : )

C...........   File units and logical/physical names

        INTEGER  :: DDEV = 0 !  unit no. for output IDA emissions file
        INTEGER  :: IDEV = 0 !  unit no. for inventory table
        INTEGER  :: VDEV = 0 !  unit no. for output IDA activity file
        INTEGER  :: LDEV = 0 !  log-device
        INTEGER  :: MDEV = 0 !  tmp unit number if ENAME is map file
        INTEGER  :: ODEV = 0 !  for output map inventory file
        INTEGER  :: RDEV = 0 !  for output ORL emissions file
        INTEGER  :: SDEV = 0 !  for ASCII input inventory file
        INTEGER  :: ZDEV = 0 !  for country/state/county file

        CHARACTER(NAMLEN3) :: ANAME = ' '! inven ASCII input logical name
        CHARACTER(NAMLEN3) :: ENAME = ' '! emis input inven logical name
        CHARACTER(NAMLEN3) :: INAME      ! tmp name for inven file of unknown fmt
        CHARACTER(NAMLEN3)    MNAME      ! tmp control/proj matrix name
        CHARACTER(NAMLEN3) :: ONAME      ! output main i/o api
        CHARACTER(NAMLEN3) :: TNAME = 'IOAPI_DAT'  ! input name for pol/act

C...........   Other local variables
                                
        INTEGER         C, I, J, K, K1, K2, L, L1, L2, 
     &                  LM, M, N, S, V                  !  counters and indices
        INTEGER         IOS        ! i/o status
        INTEGER         IYEAR      ! inventory year
        INTEGER         MXPVAR     ! max of NPPOL and NPACT
        INTEGER      :: NCALL  = 0 ! number of variables applying to "all"
        INTEGER         NCMAT      ! no. control matrices
        INTEGER         NINVARR    ! no. of inventory characteristics
        INTEGER         NNPVAR     ! no. non-pollutant inventory variables
        INTEGER         NPVAR      ! no. variables per data variable
        INTEGER         NREC       ! no. records per inventory pol/act
        INTEGER      :: PYEAR  = 0 ! projected inventory year
        INTEGER      :: PBYEAR = 0 ! projection matrix base year
        INTEGER      :: PPYEAR = 0 ! projection matrix destination year

        REAL            ALLFAC     ! tmp "all" control factor
        REAL            PSFAC      ! tmp pollutant-source control factor

        LOGICAL      :: EFLAG   = .FALSE.  ! true: error found
        LOGICAL      :: LAR2PT  = .FALSE.  ! true: area inven has xloc/yloc
        LOGICAL      :: ORLFLAG = .FALSE.  ! true: ORL output file
        LOGICAL      :: PFLAG   = .FALSE.  ! true: project matrix encountered
        LOGICAL      :: SFLAG   = .FALSE.  ! true: inventory read error
        LOGICAL      :: SMKFLAG = .FALSE.  ! true: SMOKE intermediate output file(s)
        LOGICAL      :: ADJVAL  = .FALSE.  ! true: adjust data value to 0-100 scale

        CHARACTER(16)   EVNAME  !  tmp environment variable name
        CHARACTER(80)   NAME1   ! tmp file name component
        CHARACTER(80)   NAME2   ! tmp file name component
        CHARACTER(100)  BUFFER  !  text buffer
        CHARACTER(256)  MESG    !  message buffer
        CHARACTER(IOVLEN3) :: VARBUF   ! tmp variable name
        CHARACTER(PHYLEN3) :: BPATH    ! base path
        CHARACTER(PHYLEN3) :: DATPATH  ! path for pol/act output files
        CHARACTER(PHYLEN3) :: RPATH    ! relative path
        CHARACTER(PHYLEN3) :: RFNAME   ! relative physical file name

        CHARACTER(16) :: PROGNAME= 'GRWINVEN' !  program name

C***********************************************************************
C   begin body of program GRWINVEN

        LDEV = INIT3()

C.........  Write out copyright, version, web address, header info, and prompt
C           to continue running the program.
        CALL INITEM( LDEV, CVSW, PROGNAME )

C.........  Get environment variables that control this program
        EVNAME = 'SMK_NUM_CTLMAT'
        BUFFER = 'Number of control and projection matrices'
        NCMAT = ENVINT( EVNAME, BUFFER, 1, IOS )

        IF( IOS .NE. 0 ) THEN
            WRITE( MESG,94010 ) 'WARNING: Environment variable ' //
     &             'SMK_NUM_CTLMAT is not set or is invalid.' //
     &             CRLF() // BLANK16 // 'A default value of', NCMAT,
     &             'will be used.'
            CALL M3MSG2( MESG )
        END IF

        MESG = 'First name of output inventory files'
        CALL ENVSTR( 'INVNAME1', MESG, ' ', NAME1, IOS )
        IF( IOS .NE. 0 ) THEN
            EFLAG = .TRUE.
            MESG = 'ERROR: INVNAME1 environment variable is not' //
     &             'defined for output file name'
            CALL M3MSG2( MESG )
        END IF

        MESG = 'Second name of output inventory files'
        CALL ENVSTR( 'INVNAME2', MESG, ' ', NAME2, IOS )
        IF( IOS .NE. 0 ) THEN
            EFLAG = .TRUE.
            MESG = 'ERROR: INVNAME2 environment variable is not' //
     &             'defined for output file name'
            CALL M3MSG2( MESG )
        END IF

        EVNAME = 'SMK_GRWSMKOUT_YN'
        MESG   = 'Output SMOKE inventory file'
        SMKFLAG = ENVYN( EVNAME, MESG, .TRUE., IOS )
        IF( IOS .GT. 0 ) EFLAG = .TRUE.
        
        EVNAME = 'SMK_GRWORLOUT_YN'
        MESG   = 'Output ORL inventory file'
        ORLFLAG = ENVYN( EVNAME, MESG, .FALSE., IOS )
        IF( IOS .GT. 0 ) EFLAG = .TRUE.

        IF( EFLAG ) THEN
            MESG = 'Problem with input environment variables'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

C.........  Set source category based on environment variable setting
        CALL GETCTGRY

C.........  Get inventory file names given source category
        CALL GETINAME( CATEGORY, ENAME, ANAME )

C.........  Prompt for and open inventory file 
        INAME = ENAME
        MESG = 'Enter logical name for the MAP INVENTORY file'
        MDEV = PROMPTFFILE( MESG, .TRUE., .TRUE., INAME, PROGNAME )

C.........  Open and read map file
        CALL RDINVMAP( INAME, MDEV, ENAME, ANAME, SDEV )

C.........  Open country, state, and county file
        ZDEV = PROMPTFFILE(
     &         'Enter logical name for COUNTRY, STATE, AND ' //
     &         'COUNTY file', .TRUE., .TRUE., 'COSTCY', PROGNAME )

C.........  If using ORL output format, open and read inventory table
        IF( ORLFLAG ) THEN
            IDEV = PROMPTFFILE(
     &             'Enter logical name for INVENTORY TABLE',
     &             .TRUE., .TRUE., 'INVTABLE', PROGNAME )
            CALL RDCODNAM( IDEV )
        END IF

C.........  Store source-category-specific header information, 
C           including the inventory pollutants in the file (if any).  Note that 
C           the I/O API head info is passed by include file and the
C           results are stored in module MODINFO.
        CALL GETSINFO( ENAME )

C.........  Store varible names
        NINVARR = NVARSET
        ALLOCATE( IVNAMES( NINVARR ), STAT=IOS )
        CALL CHECKMEM( IOS, 'IVNAMES', PROGNAME )
        IVNAMES( 1:NINVARR ) = VNAMESET( 1:NINVARR )

C.........  Check to see if the file has already been projected, and if
C           so, update the inventory year and print a warning
C.........  If BYEAR was not set in GETSINFO, this will be discovered later,
C           once a projection matrix is needed. It is not important if
C           a projection matrix is not being used.

        PYEAR   = GETIFDSC( FDESC3D, '/PROJECTED YEAR/', .FALSE. )
        IF( PYEAR .GT. 0 ) THEN
            IYEAR = PYEAR
            MESG = 'WARNING: Inventory file has already been ' //
     &                 'projected to a future year.'
            CALL M3MSG2( MESG )

        ELSE
            IYEAR = BYEAR

        END IF 

C.........  Read country, state, and county file for country codes
        CALL RDSTCY( ZDEV, 1, 0 )

C.........  Set number of non-pollutant inventory variables
        NNPVAR = NVARSET - NIPOL * NPPOL - NIACT * NPACT

C.........  Allocate memory based on number of control/projection matrices
        ALLOCATE( CNAMEA( NCMAT ), STAT=IOS )
        CALL CHECKMEM( IOS, 'CNAMEA', PROGNAME )
        ALLOCATE( CTYPEA( NCMAT ), STAT=IOS )
        CALL CHECKMEM( IOS, 'CTYPEA', PROGNAME )
        ALLOCATE( CINDXA( NCMAT ), STAT=IOS )
        CALL CHECKMEM( IOS, 'CINDXA', PROGNAME )
        ALLOCATE( CCNTRA( NCMAT ), STAT=IOS )
        CALL CHECKMEM( IOS, 'CCNTRA', PROGNAME )
        ALLOCATE( CNAME( NCMAT ), STAT=IOS )
        CALL CHECKMEM( IOS, 'CNAME', PROGNAME )
        ALLOCATE( CTYPE( NCMAT ), STAT=IOS )
        CALL CHECKMEM( IOS, 'CTYPE', PROGNAME )
        ALLOCATE( NCPVARS( NCMAT ), STAT=IOS )
        CALL CHECKMEM( IOS, 'NCPVARS', PROGNAME )
        ALLOCATE( IDXALL( NCMAT ), STAT=IOS )
        CALL CHECKMEM( IOS, 'IDXALL', PROGNAME )  
        ALLOCATE( CPVNAMS( MXVARS, NCMAT ), STAT=IOS )
        CALL CHECKMEM( IOS, 'CPVNAMS', PROGNAME )

C.........  Allocate memory for storing input inventory year by source
        ALLOCATE( INVYR_BASE( NSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'INVYR_BASE', PROGNAME )

C.............  Initialize arrays
        CNAMEA  = ' '  ! array
        CTYPEA  = 0    ! array
        CINDXA  = 0    ! array
        CCNTRA  = 0    ! array
        CNAME   = ' '  ! array
        CTYPE   = 0    ! array
        NCPVARS = 0    ! array
        IDXALL  = 0    ! array
        CPVNAMS = ' '  ! array

C.........  Loop through potential control and projection matrices
        M = 0
        DO I = 1, NCMAT 

C.............  Generate default name
            WRITE( MNAME, '(A5,I2.2)' ) CRL // 'CMAT', I

C.............  Prompt for name 
            WRITE( MESG, 94100 ) 'Enter logical name for ' //
     &             ' PROJECTION or CONTROL MATRIX file', I

            MNAME = PROMPTSET( MESG, FSREAD3, MNAME, PROGNAME )

            LM = LEN_TRIM( MNAME )

            IF( .NOT. DESCSET( MNAME,ALLFILES ) ) THEN
                EFLAG = .TRUE.
                MESG = 'Could not read description for "' //
     &                 MNAME( 1:LM ) // '"'
                CALL M3MSG2( MESG )
                CYCLE               ! to head of files loop
            END IF

            CNAMEA( I ) = MNAME  ! Store name in unsorted list
            CTYPEA( I ) = GETIFDSC( FDESC3D, '/CTYPE/', .TRUE. )
            CINDXA( I ) = I
            CCNTRA( I ) = I

C.............  Give error if reactivity matrix is provided
            IF( CTYPEA( I ) .EQ. CTYPREAC ) THEN
                EFLAG = .TRUE.
                MESG = 'ERROR: Grwinven does not support reactivity ' //
     &                 'matrix, as provided with file '//CRLF()//
     &                 BLANK10//TRIM(CNAMEA(I))
                CALL M3MSG2( MESG )
                CYCLE
            END IF

C............  Give error if a matrix contains more than max no of pollutatns(>500)
            IF( NVARSET > 500 ) THEN
                MESG = 'ERROR: Grwinven can not process '// 
     &                 TRIM( CNAMEA(I) ) // ' file which has '//
     &                 'more than 500 pollutants'
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF

C.............  Compare number of sources in matrix to NSRC
            CALL CHKSRCNO( CATDESC, MNAME, NROWS3D, NSRC, EFLAG )

C.............  Store number of variables in current matrix
            NCPVARS( I ) = NVARSET

C.............  For projection and control matrices, interpret variable 
C               names and compare to pollutant list. Determine whether 
C               matrix applies to the inventory or not.
            N = 0
            DO V = 1, NVARSET

                VARBUF = VNAMESET( V )
                CPVNAMS( V,I ) = VARBUF

                K = INDEX1( VARBUF, NIPPA, EANAM )

C.................  Count variables that apply to all pollutants and acts
                IF( VARBUF .EQ. 'all'  .OR. 
     &              VARBUF .EQ. 'pfac'      ) THEN
                    M = M + 1

                ELSE IF( K .LE. 0 ) THEN   ! Report names don't match with EANAM
                    N = N + 1

                    WRITE( MESG,94010 ) 'WARNING: variable "' //
     &                    TRIM( VARBUF ) // '" in matrix', I, 
     &                    CRLF() // BLANK10 // 'does not apply to ' //
     &                    'any inventory pollutants.' 
                    CALL M3MSG2( MESG )

                END IF

            END DO

C............  Give error if a matrix doesn't match with any pollutants
            IF( N .EQ. NVARSET ) THEN
                WRITE( MESG,94010 ) 'No variables in matrix', 
     &                              I, 'apply to the inventory'
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF

C.............  When a projection matrix is encountered...
            IF( CTYPEA( I ) .EQ. CTYPPROJ ) THEN

C.................  Check if there is more than one projection matrix (this
C                   is not allowed)
                IF( PFLAG ) THEN

                    EFLAG = .TRUE.
                    MESG = 'ERROR: Multiple projection matrices ' //
     &                     'encountered at matrix "' // 
     &                     MNAME( 1:LM ) // '"'
                    CALL M3MSG2( MESG )
                    CYCLE

                ELSE
                    PFLAG = .TRUE.
              
                END IF

C.................  Check if projection matrix base year is consistent with 
C                   inventory
                PBYEAR = GETIFDSC( FDESC3D, '/BASE YEAR/', .TRUE. )
                PPYEAR = GETIFDSC( FDESC3D, '/PROJECTED YEAR/', .TRUE. )

                IF( PBYEAR .NE. IYEAR ) THEN
                    EFLAG = .TRUE.
                    WRITE( MESG,94010 ) 'ERROR: Projection matrix ' //
     &                     'base year is', PBYEAR, CRLF() // BLANK10 //
     &                     'but inventory year is', IYEAR
                    CALL M3MSG2( MESG )

                ELSE
                    WRITE( MESG,94010 ) 'NOTE: Projecting from ', IYEAR,
     &                     'to', PPYEAR
                    CALL M3MSG2( MESG )

                END IF

            END IF

        END DO  ! End loop on matrices, I

        NCALL = M

C.........  Abort if error occurred while opening matrices
        IF( EFLAG ) THEN
            MESG = 'Problem with input matrices'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

C.........  Give warning if no control or projection matrices entered
C.........  This is valid if the program is only being used for format
C           conversion to IDA format
        ELSE IF( NCMAT .EQ. 0 ) THEN
            MESG = 'WARNING: No control or projection matrices. '
            CALL M3MSG2( MESG )

        END IF
        
C.........  Sort control matrices in order of precedence, and sort sorted list.
C           The sort must make sure that when two matrices of the same type
C           are present, the input order is maintained.
        CALL SORTI2( NCMAT, CINDXA, CTYPEA, CCNTRA )

        DO I = 1, NCMAT
            J = CINDXA( I )
            CNAME( I ) = CNAMEA( J )
            CTYPE( I ) = CTYPEA( J )
        END DO

C.........  Allocate memory for control factors that apply to all pollutants
C           and/or activities
        
        ALLOCATE( CFACALL( NSRC, NCALL ), STAT=IOS )
        CALL CHECKMEM( IOS, 'CFACALL', PROGNAME )
        ALLOCATE( CFAC( NSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'CFAC', PROGNAME )
        ALLOCATE( CEFF( NSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'CEFF', PROGNAME )
        ALLOCATE( REFF( NSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'REFF', PROGNAME )
        ALLOCATE( RPEN( NSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'RPEN', PROGNAME )
        CEFF = 0.
        REFF = 0.
        RPEN = 0.

C.........  Read in control matrix variables that are "all".  First, store
C           position of data in storage array, then read.
C.........  NOTE - lowercase variable names are used to permit pollutants
C           named "ALL" and "PFAC"
        N = 0
        DO I = 1, NCMAT

            J = CINDXA( I )

            K1 = INDEX1( 'all' , NCPVARS( J ), CPVNAMS( 1,J ) )
            K2 = INDEX1( 'pfac', NCPVARS( J ), CPVNAMS( 1,J ) )
            IF( K1 .GT. 0 .OR. K2 .GT. 0 ) THEN

                N = N + 1
                IDXALL( I ) = N

                IF( K1 .GT. 0 ) VARBUF = 'all'
                IF( K2 .GT. 0 ) VARBUF = 'pfac'
                IF( .NOT. READ3( CNAME( I ), VARBUF, ALLAYS3, 
     &                           0, 0, CFACALL( 1,N )         ) ) THEN

                    L = LEN_TRIM( VARBUF )
                    MESG = 'ERROR: Could not read variable "' // 
     &                     VARBUF( 1:L ) // '" from file.'
                    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

                END IF

            ELSE
                IDXALL( I ) = 0

            END IF

        END DO ! End loop on matrices

c note: the "all" feature is not documented in cntlmat because it has not
C    n: been implemented

        IVARNAMS( 1 ) = 'CIFIP'
        IVARNAMS( 2 ) = 'TZONES'
        IVARNAMS( 3 ) = 'TPFLAG'
        IVARNAMS( 4 ) = 'INVYR'

        SELECT CASE ( CATEGORY )

        CASE ( 'AREA' )
            NINVARR = 11
            IVARNAMS( 5 )  = 'CELLID'
            IVARNAMS( 6 )  = 'CISIC'
            IVARNAMS( 7 )  = 'CSCC'
            IVARNAMS( 8 )  = 'CSOURC'
            IVARNAMS( 9 )  = 'CMACT'
            IVARNAMS( 10 ) = 'CNAICS'
            IVARNAMS( 11 ) = 'CSRCTYP'

        CASE ( 'MOBILE' )
            NINVARR = 14
            IVARNAMS( 5  ) = 'IRCLAS'
            IVARNAMS( 6  ) = 'IVTYPE'
            IVARNAMS( 7  ) = 'XLOC1'
            IVARNAMS( 8  ) = 'YLOC1'
            IVARNAMS( 9  ) = 'XLOC2'
            IVARNAMS( 10 ) = 'YLOC2'
            IVARNAMS( 11 ) = 'CSCC'
            IVARNAMS( 12 ) = 'CLINK'
            IVARNAMS( 13 ) = 'CVTYPE'
            IVARNAMS( 14 ) = 'CSOURC'

        CASE ( 'POINT' )
            NINVARR = 22
            IVARNAMS( 5  ) = 'CISIC'
            IVARNAMS( 6  ) = 'XLOCA'
            IVARNAMS( 7  ) = 'YLOCA'
            IVARNAMS( 8  ) = 'STKHT'
            IVARNAMS( 9  ) = 'STKDM'
            IVARNAMS( 10 ) = 'STKTK'
            IVARNAMS( 11 ) = 'STKVE'
            IVARNAMS( 12 ) = 'CSCC'
            IVARNAMS( 13 ) = 'CORIS'
            IVARNAMS( 14 ) = 'CBLRID'
            IVARNAMS( 15 ) = 'CPDESC'
            IVARNAMS( 16 ) = 'CSOURC'
            IVARNAMS( 17 ) = 'CMACT'
            IVARNAMS( 18 ) = 'CNAICS'
            IVARNAMS( 19 ) = 'CSRCTYP'
            IVARNAMS( 20 ) = 'CERPTYP'
            IVARNAMS( 21 ) = 'CNEIUID'
            IVARNAMS( 22 ) = 'CEXTORL'

            ALLOCATE( FUGHGT( NSRC ), STAT=IOS )
            CALL CHECKMEM( IOS, 'FUGHGT', PROGNAME )
            ALLOCATE( FUGWID( NSRC ), STAT=IOS )
            CALL CHECKMEM( IOS, 'FUGWID', PROGNAME )
            ALLOCATE( FUGLEN( NSRC ), STAT=IOS )
            CALL CHECKMEM( IOS, 'FUGLEN', PROGNAME )
            ALLOCATE( FUGANG( NSRC ), STAT=IOS )
            CALL CHECKMEM( IOS, 'FUGANG', PROGNAME )

            FUGHGT   = 0.        ! array
            FUGWID   = 0.        ! array
            FUGLEN   = 0.        ! array
            FUGANG   = 0.        ! array

        END SELECT

C.........  If area sources and XLOC and YLOC are being used, then
C           set flag
        IF( CATEGORY  .EQ. 'AREA' ) THEN
            J = INDEX1( 'XLOCA', NINVARR, IVNAMES )
            IF( J .GT. 0 ) THEN
                NINVARR = 13
                IVARNAMS( 12 ) = 'XLOCA'
                IVARNAMS( 13 ) = 'YLOCA'
                LAR2PT = .TRUE.
            END IF
        END IF

C.........  Allocate memory for and read in required inventory characteristics
        CALL RDINVCHR( CATEGORY, ENAME, SDEV, NSRC, NINVARR, IVARNAMS )

C.........  Open output file(s)
        CALL OPENGRWOUT( ENAME, PPYEAR, NAME1, SMKFLAG,
     &                   ORLFLAG, ODEV, DDEV, VDEV, RDEV, 
     &                   ONAME, DATPATH )

C.........  Separate the pol/act file path into two parts to be
C           able to build relative file names for the map file.
        L = LEN_TRIM( DATPATH )
        DO N = L, 1, -1

            IF( DATPATH( N:N ) .EQ. '/' .OR.
     &          DATPATH( N:N ) .EQ. '\'      ) THEN
                BPATH = DATPATH( 1:N )
                RPATH = DATPATH( N+1:L )
                EXIT
            END IF

        END DO
        
C.........  Save the input inventory year by source for application
C           of the projection factors
C.........  Update the year of the data if the projection year is non-zero
        IF( PPYEAR .NE. 0 ) THEN
            DO S = 1, NSRC
                INVYR_BASE( S ) = INVYR( S )
                IF( INVYR( S ) .EQ. PBYEAR ) INVYR( S ) = PPYEAR
            END DO
        END IF

C.........  Write out I/O API file source characteristics, if the file has been
C           requested.  Do not write out the companion ASCII file as it would
C           be identical to that of the original inventory.
        IF( ONAME .NE. 'NONE' ) THEN
            CALL WRINVCHR( ONAME, 0, LAR2PT, .FALSE. )
        END IF

C.........  Write out the map-formatted intermediate inven header
        IF( SMKFLAG ) THEN
            WRITE( ODEV, '(A)' ) CATDESC
            WRITE( ODEV, '(A)' ) '/IOAPI/ '// TRIM( NAME1 )// '.ncf'
            WRITE( ODEV, '(A)' ) '/TEXT/ '// TRIM( NAME2 )// '.txt'
            WRITE( ODEV, '(A,I8)' ) '/NDAT/ ', NMAP
            WRITE( ODEV, '(A)' ) '/DATMAP/'
        END IF

        MESG = 'Processing inventory data...'
        CALL M3MSG2( MESG )

C.........  Loop through inventory pollutants and activities
        DO V = 1, NIPPA

            VARBUF = EANAM( V )

            K1 = INDEX1( VARBUF, NIACT, ACTVTY )

C.............  Determine name index and count depending on pol or activity
            NPVAR = NPPOL
            IF( K1 .GT. 0 ) NPVAR = NPACT

            M = INDEX1( EANAM( V ), NMAP, MAPNAM )
            CALL OPENPHYS( PROGNAME, TNAME, FSREAD3, MAPFIL(M), EFLAG )
            NREC = NROWS3D

C.............  Deallocate and allocate memory for input arrays
            IF( ALLOCATED( SRCID ) ) DEALLOCATE( SRCID, DATAVAR )
            ALLOCATE( SRCID( NREC ), STAT=IOS )
            CALL CHECKMEM( IOS, 'SRCID', PROGNAME )
            ALLOCATE( DATAVAR( NREC,NPVAR ), STAT=IOS )
            CALL CHECKMEM( IOS, 'DATAVAR', PROGNAME )
            
            DATAVAR = 0.0
            
            CALL RDINVPOL( TNAME, NREC, NPVAR, VNAMESET( 2:NPVAR+1 ), 
     &                    VTYPESET( 2:NPVAR+1 ), SRCID, DATAVAR, SFLAG )

C............  Close output file for this variable
            IF( .NOT. CLOSESET( TNAME ) ) THEN
                MESG = 'Could not close file:'//CRLF()//BLANK10//
     &                     TRIM( MAPFIL(M) )
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF

C.............  If there was a read error, then go to next variable
            IF( SFLAG ) THEN
                EFLAG = .TRUE.
                CYCLE
            END IF

C.............  Loop through matrices, if any
            DO I = 1, NCMAT

                J = CINDXA( I )
                
C.................  Search for pollutant name in list for this matrix
                K = INDEX1( VARBUF, NCPVARS( J ), CPVNAMS( 1,J ) )

C.................  Read in pollutant-specific array for multiplicative
C                   and projection matrices
                
                IF( K .GT. 0 ) THEN

                   IF( .NOT. READSET( CNAME( I ), VARBUF, ALLAYS3,
     &                                ALLFILES, 0, 0, CFAC ) ) THEN
                       MESG = 'ERROR: Could not read "' //
     &                        TRIM( VARBUF ) //'" from file "' // 
     &                        TRIM( CNAME( I ) ) // '"'
                       CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                   END IF

                   IF( CTYPE(I) .NE. CTYPPROJ ) THEN  ! not for PROJECTION packet/matrix 

                      L = LEN_TRIM( VARBUF )
                      IF( .NOT. READSET( CNAME( I ), 'CE_'//VARBUF(1:L),
     &                          ALLAYS3, ALLFILES, 0, 0, CEFF ) ) THEN
                          MESG = 'ERROR: Could not read "' //
     &                           TRIM( VARBUF ) //'" from file "' // 
     &                           TRIM( CNAME( I ) ) // '"'
                          CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                      END IF

                      IF( .NOT. READSET( CNAME( I ), 'RP_'//VARBUF(1:L),
     &                          ALLAYS3, ALLFILES, 0, 0, RPEN ) ) THEN
                          MESG = 'ERROR: Could not read "' //
     &                           TRIM( VARBUF ) //'" from file "' // 
     &                           TRIM( CNAME( I ) ) // '"'
                          CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                      END IF

                    END IF

                END IF

C.................  Set default factors
                ALLFAC = 1.
                PSFAC  = 1.

C.................  Change operation based on type of matrix, then
C                   loop through sources, applying factors as needed
                IF( CTYPE(I) .EQ. CTYPPROJ ) THEN! Projection

                    N = IDXALL( I )
                    DO C = 1, NREC

                        S = SRCID( C )   ! Get source from input index

                        IF( N .GT. 0 ) THEN
                            ALLFAC = CFACALL( S,N )
                        ELSE IF ( K .GT. 0 ) THEN
                            PSFAC  = CFAC( S ) 
                        END IF

C.........................  Skip if the base year is not the same as the
C                           base year of the projection matrix...
                        IF( INVYR_BASE( S ) .NE. PBYEAR ) CYCLE

C.........................  Adjust annual emissions or activity
                        DATAVAR( C,1 ) = 
     &                       MAX( DATAVAR( C,1 ) * ALLFAC* PSFAC, 0. )

C.........................  Adjust average-day emissions (OS_* variable)
                        IF( NPVAR .GT. 1 ) DATAVAR( C,2 ) =
     &                       MAX( DATAVAR( C,2 ) * ALLFAC * PSFAC,0. )

                    END DO

                ELSEIF( CTYPE(I) .EQ. CTYPMULT ) THEN ! Multiply (standard)

                    N = IDXALL( I )
                    DO C = 1, NREC
                        S = SRCID( C )   ! Get source from input index

                        IF( N .GT. 0 ) THEN
                            ALLFAC = CFACALL( S,N )
                        ELSE IF( K .GT. 0 ) THEN
                            PSFAC  = CFAC( S )
                        END IF

                        DATAVAR( C,1 ) = 
     &                           MAX( DATAVAR( C,1 )*ALLFAC*PSFAC, 0. )

C.........................  Adjust average-day emissions (OS_* variable)
                        IF( NPVAR .GT. 1 ) THEN
                            IF( DATAVAR( C,2 ) .GT. 0. ) 
     &                          DATAVAR(C,2)=
     &                            MAX( DATAVAR(C,2)*ALLFAC*PSFAC, 0. )

C.........................  Set the control efficiency, rule effectiveness
C                           and rule penetration to the values from the
C                           control matrix.

                            IF( CATEGORY .EQ. 'AREA' ) THEN
                              DATAVAR( C,4 ) = CEFF( S )
                              DATAVAR( C,5 ) = REFF( S )
                              DATAVAR( C,6 ) = RPEN( S )
                            ELSE IF( CATEGORY .EQ. 'POINT' ) THEN
                              DATAVAR( C,3 ) = CEFF( S )
                              DATAVAR( C,4 ) = REFF( S )
                            END IF
                        END IF

                    END DO

                END IF

            END DO  ! End loop on control/projection matrices

C.............  Write out pollutant-based variables to SMOKE file
            IF( SMKFLAG ) THEN
      
C.................  The FileSet headers are still correct from
C                   opening this pollutant's file, so no need
C                   to reset the header before opening.
                RFNAME = TRIM( RPATH )// '/'// TRIM( VARBUF )// '.ncf'     
                CALL WRINVPOL( CATEGORY, BPATH, RFNAME, NREC,
     &                         NPVAR, VARBUF, SRCID, DATAVAR, EFLAG )

C.................  Write the map inventory file entry
                WRITE( ODEV, '(A)' ) MAPNAM( M )// ' '//TRIM( RFNAME )

            END IF

C.............  Write out ORL format for current pollutant
            IF( ORLFLAG ) THEN
                CALL WRORLOUT( RDEV, VARBUF, NREC, NPVAR, SRCID,
     &                         DATAVAR, IOS )
            END IF

            IF( IOS .GT. 0 ) EFLAG = .TRUE.

        END DO  ! End loop on inventory variables


C.........  Write out the map-formatted intermediate inven header
        IF( SMKFLAG ) THEN
            WRITE( ODEV, '(A)' ) '/END/'
        END IF

        IF( EFLAG ) THEN
            MESG = 'ERROR: Could not read and write all pollutant-' //
     &             'specific variables'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

C.........  Exit program with normal completion
        CALL M3EXIT( PROGNAME, 0, 0, ' ', 0 )

C******************  FORMAT  STATEMENTS   ******************************

C...........   Informational (LOG) message formats... 92xxx

92000   FORMAT( 5X, A )

92010   FORMAT( 5X, A, :, I10 )


C...........   Formatted file I/O formats............ 93xxx

93000   FORMAT( A )

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

94040   FORMAT( A, I2.2 )

94060   FORMAT( 10( A, :, E10.3, :, 1X ) )

94100   FORMAT( 10( A, :, I3, :, 1X ) )

        END PROGRAM GRWINVEN
