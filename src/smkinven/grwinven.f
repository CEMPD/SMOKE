
         PROGRAM GRWINVEN

C***************************************************************************
C  program body starts at line 
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
C***************************************************************************

C...........   MODULES for public variables
C...........   This module is the inventory arrays
        USE MODSOURC

C...........   This module contains the information about the source category
        USE MODINFO

        IMPLICIT NONE

C...........   INCLUDES:

        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
        INCLUDE 'FDESC3.EXT'    !  I/O API file description data structures.

C...........   EXTERNAL FUNCTIONS and their descriptions:
        
        CHARACTER*2             CRLF
        INTEGER                 ENVINT
        INTEGER                 GETIFDSC
        INTEGER                 INDEX1
        INTEGER                 PROMPTFFILE
        CHARACTER(LEN=NAMLEN3)  PROMPTMFILE
        INTEGER                 STR2INT

        EXTERNAL CRLF, ENVINT, GETIFDSC,
     &           INDEX1, PROMPTFFILE, PROMPTMFILE, STR2INT

C...........  LOCAL PARAMETERS and their descriptions:

        CHARACTER*50, PARAMETER :: SCCSW = '@(#)$Id$'

C...........  LOCAL VARIABLES and their descriptions:

C.........  Array that contains the names of the inventory variables needed 
C           for this program
        CHARACTER(LEN=IOVLEN3) IVARNAMS( MXINVARR )

C...........  ALLOCATABLE VARIABLES and their decriptions...

C...........  Variables for each control/projection matrix
        INTEGER, ALLOCATABLE :: CCNTRA ( : ) ! unsorted original position
        INTEGER, ALLOCATABLE :: CINDXA ( : ) ! sorting index
        INTEGER, ALLOCATABLE :: CTYPEA ( : ) ! unsorted control matrix type code
        INTEGER, ALLOCATABLE :: CTYPE  ( : ) ! sorted control matrix type code
        INTEGER, ALLOCATABLE :: IDXALL ( : ) ! index of matrixes with "All"
        INTEGER, ALLOCATABLE :: NCPVARS( : ) ! sorted no. variables

        CHARACTER(LEN=NAMLEN3), ALLOCATABLE :: CNAMEA ( : )  ! unsort mtx names
        CHARACTER(LEN=NAMLEN3), ALLOCATABLE :: CNAME  ( : )  ! sorted mtx names
        CHARACTER(LEN=NAMLEN3), ALLOCATABLE :: CPVNAMS( :,: )! mtx var names

C...........  Control matrix factors
        REAL, ALLOCATABLE :: CFAC   ( : )  ! by-variable factors
        REAL, ALLOCATABLE :: CFACALL( :,: )  ! all-variable factors

C...........  Temporary IDA pollutant files unit nos.
        INTEGER, ALLOCATABLE :: TDEV( : )

C...........  FIXED DIMENSION VARIABLES and their descriptions...

C...........  Inventory file variable names
        CHARACTER(LEN=IOVLEN3) :: IVNAMES( MXVARS3 )

C...........  Matrix variables not matching any inventory pollutants
        CHARACTER(LEN=IOVLEN3) CVODDLST( MXVARS3 ) 

C...........   File units and logical/physical names

        INTEGER  :: DDEV = 0 !  unit no. for output IDA emissions file
        INTEGER  :: VDEV = 0 !  unit no. for output IDA activity file
        INTEGER  :: LDEV = 0 !  log-device
        INTEGER  :: SDEV = 0 !  for ASCII input inventory file
        INTEGER  :: ZDEV = 0 !  for country/state/county file

        CHARACTER(LEN=NAMLEN3) :: ANAME = ' '! inven ASCII input logical name
        CHARACTER(LEN=NAMLEN3) :: ENAME = ' '!  emis input inven logical name
        CHARACTER(LEN=NAMLEN3)    MNAME      !  tmp control/proj matrix name
        CHARACTER(LEN=NAMLEN3) :: ONAME = ' '!  emis output inven logical name

C...........   Other local variables
                                
        INTEGER         I, IS, J, K, K1, K2, L, L1, L2, 
     &                  LM, M, N, S, V                  !  counters and indices

        INTEGER         IOS        ! i/o status
        INTEGER         IYEAR      ! inventory year
        INTEGER      :: NCALL  = 0 ! number of variables applying to "all"
        INTEGER         NCMAT      ! no. control matrices
        INTEGER         NINVARR    ! no. of inventory characteristics
        INTEGER         NNPVAR     ! no. non-pollutant inventory variables
        INTEGER      :: PYEAR  = 0 ! projected inventory year
        INTEGER      :: PBYEAR = 0 ! projection matrix base year
        INTEGER      :: PPYEAR = 0 ! projection matrix destination year

        REAL            ALLFAC     ! tmp "all" control factor
        REAL            PSFAC      ! tmp pollutant-source control factor

        LOGICAL      :: EFLAG   = .FALSE.  ! true: error found
        LOGICAL      :: IDAFLAG = .FALSE.  ! true: IDA output file(s) 
        LOGICAL      :: PFLAG   = .FALSE.  ! true: project matrix encountered

        CHARACTER*16    EVNAME  !  tmp environment variable name
        CHARACTER*100   BUFFER  !  text buffer
        CHARACTER*300   MESG    !  message buffer

        CHARACTER(LEN=IOVLEN3) :: VARBUF  ! tmp variable name

        CHARACTER*16 :: PROGNAME= 'GRWINVEN' !  program name

C***********************************************************************
C   begin body of program GRWINVEN

        LDEV = INIT3()

C.........  Write out copywrite, version, web address, header info, and prompt
C           to continue running the program.
        CALL INITEM( LDEV, SCCSW, PROGNAME )

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

C.........  Set source category based on environment variable setting
        CALL GETCTGRY

C.........  Get inventory file names given source category
        CALL GETINAME( CATEGORY, ENAME, ANAME )

C.........  Prompt for and open I/O API inventory input file
        ENAME = PROMPTMFILE( 
     &          'Enter logical name for the I/O API INVENTORY file',
     &          FSREAD3, ENAME, PROGNAME )

C.........  Prompt for and open ASCII inventory input file
        SDEV = PROMPTFFILE( 
     &         'Enter logical name for the ASCII INVENTORY file',
     &         .TRUE., .TRUE., ANAME, PROGNAME )

        ZDEV = PROMPTFFILE(
     &         'Enter logical name for COUNTRY, STATE, AND ' //
     &         'COUNTY file', .TRUE., .TRUE., 'COSTCY', PROGNAME )

C.........  Get header description of inventory file, error if problem
        IF( .NOT. DESC3( ENAME ) ) THEN
            L2 = LEN_TRIM( ENAME )
            MESG = 'Could not get description of file "' //
     &             ENAME( 1:L2 ) // '"'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

C.........  Otherwise, store source-category-specific header information, 
C           including the inventory pollutants in the file (if any).  Note that 
C           the I/O API head info is passed by include file and the
C           results are stored in module MODINFO.
        ELSE

            CALL GETSINFO

C.............  Store varible names
            IVNAMES( 1:NVARS3D ) = VNAME3D( 1:NVARS3D )

C.............  Check to see if the file has already been projected, and if
C               so, update the inventory year and print a warning
C.............  If BYEAR was not set in GETSINFO, this will be discovered later,
C               once a projection matrix is needed. It is not important if
C               a projection matrix is not being used.

            PYEAR   = GETIFDSC( FDESC3D, '/PROJECTED YEAR/', .FALSE. )
            IF( PYEAR .GT. 0 ) THEN
                IYEAR = PYEAR
                MESG = 'WARNING: Inventory file has already been ' //
     &                 'projected to a future year.'
                CALL M3MSG2( MESG )

            ELSE
                IYEAR = BYEAR

            END IF 
         END IF

C.........  Read country, state, and county file for country codes
        CALL RDSTCY( ZDEV, 1, 0 )

C.........  Set number of non-pollutant inventory variables
        NNPVAR = NVARS3D - NIPOL * NPPOL - NIACT * NPACT

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
        ALLOCATE( CPVNAMS( MXVARS3, NCMAT ), STAT=IOS )
        CALL CHECKMEM( IOS, 'CPVNAMS', PROGNAME )

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

            MNAME = PROMPTMFILE( MESG, FSREAD3, MNAME, PROGNAME )

            LM = LEN_TRIM( MNAME )

            IF( .NOT. DESC3( MNAME ) ) THEN
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

C.............  Compare number of sources in matrix to NSRC 
            CALL CHKSRCNO( CATDESC, MNAME, NROWS3D, NSRC, EFLAG )

C.............  Store number of variables in current matrix
            NCPVARS( I ) = NVARS3D

C.............  Interpret variable names and compare to pollutant list.  Keep
C               track of those that are not in the inventory file.
            N = 0
            DO V = 1, NVARS3D
 
                VARBUF = VNAME3D( V )
                CPVNAMS( V,I ) = VARBUF

                K = INDEX1( VARBUF, NIPPA, EANAM )

C.................  Count variables that apply to all pollutants and activities
                IF( VARBUF .EQ. 'all'  .OR. 
     &              VARBUF .EQ. 'pfac'      ) THEN
                    M = M + 1

                ELSE IF( K .EQ. -1 ) THEN   ! Store names don't match with EANAM
                    N = N + 1
                    CVODDLST( N ) = VARBUF

                    L2 = LEN_TRIM( VARBUF )
                    WRITE( MESG,94010 ) 'WARNING: factors in matrix', I, 
     &                 'for variable "' // VARBUF( 1:L2 ) //
     &                 '"' // CRLF() // BLANK10 // 'do not apply to ' //
     &                 'inventory.' 
                    CALL M3MSG2( MESG )

                END IF

            END DO

C.............  Report matrices that don't match with any pollutants or 
C               activities
            IF( N .EQ. NVARS3D ) THEN
                WRITE( MESG,94010 ) 'WARNING: No variables in matrix', 
     &                              I, 'apply to the inventory'
                CALL M3MSG2( MESG )

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
            MESG = 'Problem with control or projection matrices'
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

C.........  Read in control matrix variables that are "all".  First, store
C           position of data in storage array, then read.
C.........  NOTE - lowercase variable names are used to permit pollutants
C           named "ALL" and "PFAC"
        J = 0
        DO I = 1, NCMAT

            K1 = INDEX1( 'all' , NCPVARS( I ), CPVNAMS( 1,I ) )
            K2 = INDEX1( 'pfac', NCPVARS( I ), CPVNAMS( 1,I ) )
            IF( K1 .GT. 0 .OR. K2 .GT. 0 ) THEN

                J = J + 1
                IDXALL( I ) = J

                IF( K1 .GT. 0 ) VARBUF = 'all'
                IF( K2 .GT. 0 ) VARBUF = 'pfac'
                IF( .NOT. READ3( CNAME( I ), VARBUF, ALLAYS3, 
     &                           0, 0, CFACALL( 1,J )         ) ) THEN

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

C.........  Set inventory variables to read depending on source category
        IVARNAMS( 1 ) = 'IFIP'
        IVARNAMS( 2 ) = 'TZONES'
        IVARNAMS( 3 ) = 'TPFLAG'
        IVARNAMS( 4 ) = 'INVYR'

        SELECT CASE ( CATEGORY )

        CASE ( 'AREA' )
            NINVARR = 7
            IVARNAMS( 5 ) = 'CELLID'
            IVARNAMS( 6 ) = 'CSCC'
            IVARNAMS( 7 ) = 'CSOURC'

        CASE ( 'MOBILE' )
            NINVARR = 15
            IVARNAMS( 5  ) = 'IRCLAS'
            IVARNAMS( 6  ) = 'IVTYPE'
            IVARNAMS( 7  ) = 'XLOC1'
            IVARNAMS( 8  ) = 'YLOC1'
            IVARNAMS( 9  ) = 'XLOC2'
            IVARNAMS( 10 ) = 'YLOC2'
            IVARNAMS( 11 ) = 'SPEED'
            IVARNAMS( 12 ) = 'CSCC'
            IVARNAMS( 13 ) = 'CLINK'
            IVARNAMS( 14 ) = 'CVTYPE'
            IVARNAMS( 15 ) = 'CSOURC'

        CASE ( 'POINT' )
            NINVARR = 16
            IVARNAMS( 5  ) = 'ISIC'
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

        END SELECT

C.........  Allocate memory for and read in required inventory characteristics
        CALL RDINVCHR( CATEGORY, ENAME, SDEV, NSRC, NINVARR, IVARNAMS )

C.........  Open output file(s)
        CALL OPENGRWOUT( ENAME, PPYEAR, ONAME, DDEV, VDEV )

        IDAFLAG = ( DDEV .GT. 0 .OR. VDEV .GT. 0 )

C.........  Update the year of the data if the projection year is non-zero
        IF( PPYEAR .NE. 0 ) THEN
            INVYR = PPYEAR       ! array
        END IF

C.........  Write out I/O API file source characteristics, if the file has been
C           requested.  Do not write out the companion ASCII file as it would
C           be identical to that of the original inventory.
        IF( ONAME .NE. 'NONE' ) THEN

            CALL WRINVCHR( ONAME, 0 )

        END IF

C.........  Processing when IDA output is needed...
        IF( IDAFLAG ) THEN

C.............  Allocate memory for unit numbers for IDA temporary files
            ALLOCATE( TDEV( NIPPA ), STAT=IOS )
            CALL CHECKMEM( IOS, 'TDEV', PROGNAME )
            TDEV = 0   ! array
 
C.........  If no IDA output, deallocate memory for source characteristics
        ELSE

            CALL SRCMEM( CATEGORY, 'SORTED', .FALSE., .FALSE., 1, 1, 1 )

        END IF

C.........  Allocate memory for pollutant- and activity-specific data (one
C           pollutant or activity at a time)
        CALL SRCMEM( CATEGORY, 'SORTED', .TRUE., .TRUE., NSRC, 
     &               NSRC, NPPOL )

        MESG = 'Processing inventory data...'
        CALL M3MSG2( MESG )

C.........  Loop through inventory pollutants and activities
        DO V = 1, NIPPA

            VARBUF = EANAM( V )

C.............  Read inventory pollutant or activity from I/O API file
            IS = NNPVAR + ( V-1 ) * NPPOL + 1
            CALL RDINVPOL( ENAME, NSRC, NPPOL, 0, 0, IVNAMES( IS ),
     &                     POLVAL, IOS )

C.............  If there was a read error, then go to next variable
            IF( IOS .GT. 0 ) THEN
                EFLAG = .TRUE.
                CYCLE
            END IF

C.............  Loop through matrices, if any
            DO I = 1, NCMAT

                J = CINDXA( I )

C.................  Search for pollutant name in list for this matrix
                K = INDEX1( VARBUF, NCPVARS( J ), CPVNAMS( 1,J ) )

C.................  Read in pollutant-specific array 
                IF( K .GT. 0 ) THEN
                    IF( .NOT. READ3( CNAME( I ), VARBUF, 
     &                               ALLAYS3, 0, 0, CFAC ) ) THEN

                        L1 = LEN_TRIM( VARBUF )
                        L2 = LEN_TRIM( CNAME( I ) )
                        MESG = 'ERROR: Could not read "' //
     &                         VARBUF( 1:L1 ) //'" from file "' // 
     &                         CNAME( I )( 1:L2 )
                        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

                    END IF
                END IF

C.................  Set default factors
                ALLFAC = 1.
                PSFAC  = 1.

C.................  Change operation based on type of matrix, then
C                   loop through sources, applying factors as needed
                IF( CTYPE(I) .EQ. CTYPPROJ ) THEN! Projection

                    N = IDXALL( I )
                    DO S = 1, NSRC
                        IF( N .GT. 0 ) ALLFAC = CFACALL( S,N )
                        IF( K .GT. 0 ) PSFAC  = CFAC( S )
                        POLVAL( S,1 ) = POLVAL( S,1 ) * ALLFAC* PSFAC
                    END DO

                ELSE IF( CTYPE(I) .EQ. CTYPREAC ) THEN ! Reactivity 

                    DO S = 1, NSRC
                        POLVAL( S,1 ) = POLVAL( S,1 ) * CFAC( S )
                    END DO

                ELSEIF( CTYPE(I) .EQ. CTYPADD .AND. K .GT. 0 ) THEN ! Add

                    DO S = 1, NSRC
                        POLVAL( S,1 ) = POLVAL( S,1 ) + CFAC( S )
                    END DO

                ELSEIF( CTYPE(I) .EQ. CTYPMULT ) THEN ! Multiply (standard)

                    N = IDXALL( I )
                    DO S = 1, NSRC
                        IF( N .GT. 0 ) ALLFAC = CFACALL( S,N )
                        IF( K .GT. 0 ) PSFAC  = CFAC( S )
                        POLVAL( S,1 ) = POLVAL( S,1 )* ALLFAC* PSFAC
                    END DO

                END IF

            END DO  ! End loop on control/projection matrices

C.............  Write out pollutant-based variables to SMOKE file
            IF( ONAME .NE. 'NONE' ) THEN
                CALL WRINVPOL( ONAME, CATEGORY, NSRC, 1, NPPOL, 
     &                         IVNAMES( IS ), POLVAL )

            END IF

            IF( IOS .GT. 0 ) EFLAG = .TRUE.

C.............  Write out pollutant-based variables to temporary files for IDA
            IF( IDAFLAG ) THEN
                    
                CALL WRIDAPOL( CATEGORY, VARBUF, NSRC, NPPOL, POLVAL, 
     &                         TDEV( V ), IOS )

            END IF

            IF( IOS .GT. 0 ) EFLAG = .TRUE.

        END DO  ! End loop on inventory variables
            
        IF( EFLAG ) THEN
            MESG = 'ERROR: Could not read and write all pollutant-' //
     &             'specific variables'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

C.........  Deallocate input pollutant arrays
        CALL SRCMEM( CATEGORY, 'SORTED', .FALSE., .TRUE., 
     &               NSRC, NSRC, NPPOL )

C.........  For IDA output
        IF( IDAFLAG ) THEN

            MESG = 'Writing IDA-formatted outputs...'
            CALL M3MSG2( MESG )

            CALL WRIDAOUT( DDEV, VDEV, TDEV, IOS )

            IF( IOS .GT. 0 ) EFLAG = .TRUE.

        END IF
      
        IF( EFLAG ) THEN
            MESG = 'ERROR: Could not read and write all source ' //
     &             'characteristics for IDA output.'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

C.........  End program successfully
        MESG = ' '
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 0 )

C******************  FORMAT  STATEMENTS   ******************************

C...........   Informational (LOG) message formats... 92xxx

92000   FORMAT( 5X, A )

92010   FORMAT( 5X, A, :, I10 )


C...........   Formatted file I/O formats............ 93xxx

93000   FORMAT( A )

c93500   FORMAT( A248, <NIPOL>(A52) )


C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

94020   FORMAT( A, 1X, I5.5, 1X, A, 1X, I8.8, 1X,
     &          A, I6, X, A, I6, X, A, :, I6 )

94040   FORMAT( A, I2.2 )

94060   FORMAT( 10( A, :, E10.3, :, 1X ) )

94080   FORMAT( '************  ', A, I7, ' ,  ' , A, I12 )
 
94100   FORMAT( 10( A, :, I3, :, 1X ) )

        END

