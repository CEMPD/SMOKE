
        SUBROUTINE WRINVEMIS( IDEV, DATPATH, VAR_FORMULA )

C***********************************************************************
C  subroutine body starts at line 135
C
C  DESCRIPTION:
C      This subroutine writes the average inventory emissions to the inventory
C      files
C
C  PRECONDITIONS REQUIRED:
C      Logical name of output file defined
C      Emission arrays populated
C      MODINFO values assigned
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C      Subroutines: I/O API subroutine
C
C  REVISION  HISTORY:
C      Created 12/99 by M. Houyoux
C
C*************************************************************************
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

C.........  MODULES for public variables
C...........   This module is the inventory arrays
        USE MODSOURC, ONLY: CSOURC, NPCNT, IPOSCOD, POLVAL

C.........  This module contains the lists of unique inventory information
        USE MODLISTS, ONLY: MXIDAT, INVDNAM, INVDUNT, FIREFLAG

C.........  This module contains the information about the source category
        USE MODINFO, ONLY: CATEGORY, CATDESC, NSRC, NMAP,
     &                     NIPOL, NIACT, NIPPA, NPPOL, EANAM, ACTVTY,
     &                     MAPNAM, MAPFIL, EINAM, EIIDX, NPACT, AVIDX,
     &                     NDY, NC1, NC2, NCE, NRE, NRP

C.........  This module is required by the FileSetAPI
        USE MODFILESET
        
       IMPLICIT NONE

C...........   INCLUDES

        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
        INCLUDE 'SETDECL.EXT'   !  FileSetAPI variables and functions

C.........  EXTERNAL FUNCTIONS
        CHARACTER(2)  CRLF
        LOGICAL       ENVYN
        INTEGER       INDEX1
        LOGICAL       SETENVVAR

        EXTERNAL      CRLF, ENVYN, INDEX1, SETENVVAR

C.........  SUBROUTINE ARGUMENTS
        INTEGER     , INTENT (IN) :: IDEV           ! unit number for map file
        CHARACTER(PHYLEN3), INTENT (IN) :: DATPATH ! path for pol/act files
        CHARACTER(*), INTENT (IN) :: VAR_FORMULA    ! formula string

C...........   LOCAL PARAMETERS
        CHARACTER(16), PARAMETER :: FORMEVNM = 'SMKINVEN_FORMULA'

C.........  Inventory temporay arrays
        INTEGER, ALLOCATABLE:: IPPTR ( : ) ! position in POLVAL sparse array
        INTEGER, ALLOCATABLE:: IPMAX ( : ) ! max IPPTR by source
        INTEGER, ALLOCATABLE:: SRCID ( : ) ! source index
        REAL   , ALLOCATABLE:: SRCPOL( :,: )  ! data-spec values by source
        REAL   , ALLOCATABLE:: COMPUTED( :,: )! computed data-spec values by src

C...........   Names, Units, types, & descriptions for pollutant-specific 
C              output variables.  NOTE - second dimension will work so long
C              as NPPOL > NPACT, which is expected to always be the case

        CHARACTER(IOVLEN3), ALLOCATABLE :: EONAMES( :,: ) ! Names 
        INTEGER           , ALLOCATABLE :: EOTYPES( :,: ) ! Types (Real|Int)
        CHARACTER(IOULEN3), ALLOCATABLE :: EOUNITS( :,: ) ! Units  
        CHARACTER(IODLEN3), ALLOCATABLE :: EODESCS( :,: ) ! Dscriptions  

C...........   Other local allocatable arrays
        CHARACTER(IOVLEN3), ALLOCATABLE :: SAVEANAM( : )  ! tmp variables

C...........   Other local variables
        INTEGER         I, S, L, L2, N, V1, V2, VA, VB     ! counters and indices

        INTEGER         IOS       ! i/o status
        INTEGER         LEQU      ! position of '=' in formula
        INTEGER         LDIV      ! position of '-' or '+' in formula
        INTEGER         LMNS      ! position of '-' in formula
        INTEGER         LPLS      ! position of '+' in formula
        INTEGER         MCNT      ! count of actual mapped pol/act files

        INTEGER         RIMISS3          ! real value of integer missing

        LOGICAL      :: CHKPLUS  = .FALSE. ! true: formula uses a + sign
        LOGICAL      :: CHKMINUS = .FALSE. ! true: formula uses a - sign
        LOGICAL      :: EFLAG    = .FALSE. ! true: error found
        LOGICAL      :: FFLAG    = .FALSE. ! true: formula in use
        LOGICAL         ZFLAG              ! true: write zeros to output file

        CHARACTER(80)   NAME1            ! tmp file name component
        CHARACTER(80)   NAME2            ! tmp file name component
        CHARACTER(128)  BUFFER           ! message buffer
        CHARACTER(256)  MESG             ! message buffer

        CHARACTER(IOVLEN3 ) VIN_A    ! first variable in equation
        CHARACTER(IOVLEN3 ) VIN_B    ! second variable in equation
        CHARACTER(IOVLEN3 ) VNAME    ! computed variable name
        CHARACTER(PHYLEN3) :: BPATH  ! base path
        CHARACTER(PHYLEN3) :: RPATH  ! relative path

        CHARACTER(16) :: PROGNAME = 'WRINVEMIS' !  program name

C***********************************************************************
C   begin body of program WRINVEMIS

C..........  Get environment variables
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

        IF( EFLAG ) THEN
            MESG = 'Problem with input environment variables'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

        MESG = 'Write zero values to annual emissions inventory'
        ZFLAG = ENVYN( 'WRITE_ANN_ZERO', MESG, .FALSE., IOS )

C.........  Compute real value of integer missing
        RIMISS3 = REAL( IMISS3 )

C.........  Allocate memory for local source-specific arrays used for output   
C.........  Allocate memory for indices IPPTR & IPMAX for pointing to position
C           in sparsely stored data array.  
        ALLOCATE( IPPTR( NSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'IPPTR', PROGNAME )
        ALLOCATE( IPMAX( NSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'IPMAX', PROGNAME )
        ALLOCATE( SRCPOL( NSRC, NPTPPOL3 ), STAT=IOS )
        CALL CHECKMEM( IOS, 'SRCPOL', PROGNAME )
        ALLOCATE( SRCID( NSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'SRCID', PROGNAME )
       
C.........  Initialize global index based on count of number of sources per 
C           pollutant
C.........  Since the emission values are already sorted in the output order of
C           the pollutants/activities, can use the IPPTR and IPMAX indices 
C           to keep track of what pollutant/activity we are on for each source
C           (sources can have different numbers and types of output pollutants
C           and activities).
        IPPTR( 1 ) = 1
        IPMAX( 1 ) = NPCNT( 1 )
        DO S = 2, NSRC
            IPPTR( S ) = IPPTR( S-1 ) + NPCNT( S-1 )
            IPMAX( S ) = IPPTR( S )   + NPCNT( S ) - 1
        END DO

C.........  Set maximum number of map variables for map-formatted outputs
        NMAP = NIPOL + NIACT
        MCNT = 0              ! initialize actual pol/act file count for later

C.........  If there is a computed output variable, get set up for that
        L = LEN_TRIM( VAR_FORMULA )
        IF( L .GT. 0 ) THEN

            FFLAG = .TRUE.
            NMAP = NMAP + 1  ! one more variable to map

C.............  Make sure formula makes sense
            LEQU = INDEX( VAR_FORMULA, '=' )
            LPLS = INDEX( VAR_FORMULA, '+' )
            LMNS = INDEX( VAR_FORMULA, '-' )

            CHKPLUS  = ( LPLS .GT. 0 )
            CHKMINUS = ( LMNS .GT. 0 )

            LDIV = LPLS
            IF( CHKMINUS ) LDIV = LMNS

            IF( LEQU .LE. 0 .OR. 
     &        ( .NOT. CHKPLUS .AND. .NOT. CHKMINUS ) ) THEN

                L = LEN_TRIM( FORMEVNM )
                MESG = 'Could not interpret formula for extra ' //
     &                 'pollutant from environment variable ' //
     &                 CRLF() // BLANK10 // '"' // FORMEVNM( 1:L ) //
     &                 '": ' // VAR_FORMULA
                EFLAG = .TRUE.
            END IF

C.............  Extract formula variable names
            L     = LEN_TRIM( VAR_FORMULA )
            VNAME = ADJUSTL ( VAR_FORMULA(      1:LEQU-1 ) )
            VIN_A = ADJUSTL ( VAR_FORMULA( LEQU+1:LDIV-1 ) )
            VIN_B = ADJUSTL ( VAR_FORMULA( LDIV+1:L      ) )

C.............  Find formula inputs in existing variable list
            VA = INDEX1( VIN_A, NIPPA, EANAM )
            VB = INDEX1( VIN_B, NIPPA, EANAM )

            IF( VA .LE. 0 ) THEN
                EFLAG = .TRUE.
                L = LEN_TRIM( VIN_A )
                MESG = 'Variable "'// VIN_A( 1:L ) // 
     &                 '" from formula was not found in inventory.'
                CALL M3MSG2( MESG )
            END IF

            IF( VB .LE. 0 ) THEN
                EFLAG = .TRUE.
                L = LEN_TRIM( VIN_B )
                MESG = 'Variable "'// VIN_B( 1:L ) // 
     &                 '" from formula was not found in inventory.'
                CALL M3MSG2( MESG )
            END IF

            V1 = INDEX1( VIN_A, NIACT, ACTVTY )
            V2 = INDEX1( VIN_B, NIACT, ACTVTY )

            IF( V1 .GT. 0 ) THEN
                EFLAG = .TRUE.
                L = LEN_TRIM( VIN_A )
                MESG = 'ERROR: Variable "'// VIN_A( 1:L )// '" is an'//
     &                 'activity, which is not allowed in a formula.'
                CALL M3MSG2( MESG )
            END IF

            IF( V2 .GT. 0 ) THEN
                EFLAG = .TRUE.
                L = LEN_TRIM( VIN_B )
                MESG = 'ERROR: Variable "'// VIN_B( 1:L )// '" is an'//
     &                 'activity, which is not allowed in a formula.'
                CALL M3MSG2( MESG )
            END IF

            IF( EFLAG ) THEN
                MESG = 'Problem processing formula.'
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF

C.............  Allocate memory for computed variable
            ALLOCATE( COMPUTED( NSRC, NPTPPOL3 ), STAT=IOS )
            CALL CHECKMEM( IOS, 'COMPUTED', PROGNAME )
            COMPUTED = 0.   ! array

        END IF

C.........  Allocate memory for map file arrays
        ALLOCATE( MAPNAM( NMAP ), STAT=IOS )
        CALL CHECKMEM( IOS, 'MAPNAM', PROGNAME )
        ALLOCATE( MAPFIL( NMAP ), STAT=IOS )
        CALL CHECKMEM( IOS, 'MAPFIL', PROGNAME )
        MAPNAM = ' '
        MAPFIL = ' '

C.........  Allocate memory for temporary variable names etc.
        ALLOCATE( EONAMES( NIPOL,NPPOL ), STAT=IOS )
        CALL CHECKMEM( IOS, 'EONAMES', PROGNAME )
        ALLOCATE( EOUNITS( NIPOL,NPPOL ), STAT=IOS )
        CALL CHECKMEM( IOS, 'EOUNITS', PROGNAME )
        ALLOCATE( EOTYPES( NIPOL,NPPOL ), STAT=IOS )
        CALL CHECKMEM( IOS, 'EOTYPES', PROGNAME )
        ALLOCATE( EODESCS( NIPOL,NPPOL ), STAT=IOS )
        CALL CHECKMEM( IOS, 'EODESCS', PROGNAME )

C.........  Get names, units, etc. of output pollutant-specific records
        CALL BLDENAMS( CATEGORY, NIPOL, NPPOL, EINAM, 
     &                 EONAMES, EOUNITS, EOTYPES, EODESCS )

C.........  Set up for opening I/O API sparse pollutant output files
        CALL HDRMISS3  ! Initialize for emissions

C.........  Set number of variables and allocate file description arrays
        NVARSET = 1 + NPPOL   ! Additional 1 for the SRCID variable
        WRITE( FDESC3D( 1 ), '(A,1X,I8)' ) '/NSRC/', NSRC

        IF( ALLOCATED( VTYPESET ) ) 
     &      DEALLOCATE( VTYPESET, VNAMESET, VUNITSET, VDESCSET )
        ALLOCATE( VTYPESET( NVARSET ), STAT=IOS )
        CALL CHECKMEM( IOS, 'VTYPESET', PROGNAME )
        ALLOCATE( VNAMESET( NVARSET ), STAT=IOS )
        CALL CHECKMEM( IOS, 'VNAMESET', PROGNAME )
        ALLOCATE( VUNITSET( NVARSET ), STAT=IOS )
        CALL CHECKMEM( IOS, 'VUNITSET', PROGNAME )
        ALLOCATE( VDESCSET( NVARSET ), STAT=IOS )
        CALL CHECKMEM( IOS, 'VDESCSET', PROGNAME )
        
        IF( ALLOCATED( VARS_PER_FILE ) ) DEALLOCATE( VARS_PER_FILE )
        
        VTYPESET = 0    ! array initialization
        VNAMESET = ' '  ! array initialization
        VUNITSET = ' '  ! array initialization
        VDESCSET = ' '  ! array initialization

C.........  Set up source ID variable information, which is the same for
C           all pollutant files
        VTYPESET( 1 ) = M3INT
        VNAMESET( 1 ) = 'SRCID'
        VUNITSET( 1 ) = 'n/a'
        VDESCSET( 1 ) = 'Source ID number'

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
C.........  Loop through pollutants, store, and write to inventory file
        IF( NIPOL .GT. 0 ) 
     &      CALL LOOP_FOR_OUTPUT( NIPOL, NPPOL, IDEV, EIIDX, 
     &                            EINAM, MCNT )

C.........  If computed variable, update EANAM, in case needed for day- 
C           and hour-specific data
        IF( FFLAG ) THEN

            ALLOCATE( SAVEANAM( NIPPA ), STAT=IOS )
            CALL CHECKMEM( IOS, 'SAVEANAM', PROGNAME )
            SAVEANAM = EANAM       ! array

            DEALLOCATE( EANAM )
            ALLOCATE( EANAM( NIPPA+1 ), STAT=IOS )
            CALL CHECKMEM( IOS, 'EANAM', PROGNAME )

            EANAM( 1:NIPOL ) = SAVEANAM( 1:NIPOL )
            EANAM( NIPOL+1 ) = VNAME
            EANAM( NIPOL+2:NIPPA+1 ) = SAVEANAM( NIPOL+1:NIPPA )
            NIPPA = NIPPA + 1
            DEALLOCATE( SAVEANAM )

        END IF

C.........  Output for activities...

C.........  Deallocate arrays for variable names
        DEALLOCATE( EONAMES, EOUNITS, EOTYPES, EODESCS )

C.........  Allocate memory for temporary variable names etc.
        ALLOCATE( EONAMES( NIACT,NPACT ), STAT=IOS )
        CALL CHECKMEM( IOS, 'EONAMES', PROGNAME )
        ALLOCATE( EOUNITS( NIACT,NPACT ), STAT=IOS )
        CALL CHECKMEM( IOS, 'EOUNITS', PROGNAME )
        ALLOCATE( EOTYPES( NIACT,NPACT ), STAT=IOS )
        CALL CHECKMEM( IOS, 'EOTYPES', PROGNAME )
        ALLOCATE( EODESCS( NIACT,NPACT ), STAT=IOS )
        CALL CHECKMEM( IOS, 'EODESCS', PROGNAME )
        
C.........  Get names, units, etc. of output activity-specific records
        CALL BLDENAMS( CATEGORY, NIACT, NPACT, ACTVTY, 
     &                 EONAMES, EOUNITS, EOTYPES, EODESCS )

C.........  Set up for opening I/O API sparse pollutant output files
        CALL HDRMISS3  ! Initialize for emissions

C.........  Set number of variables and allocate file description arrays
        NVARSET = 1 + NPACT
        WRITE( FDESC3D( 1 ), '(A,1X,I8)' ) '/NSRC/', NSRC
        
        IF( ALLOCATED( VARS_PER_FILE ) ) DEALLOCATE( VARS_PER_FILE )

C.........  Loop through activity data, store, and write to inventory file
        IF( NIACT .GT. 0 ) 
     &      CALL LOOP_FOR_OUTPUT( NIACT, NPACT, IDEV, AVIDX, 
     &                            ACTVTY, MCNT )

C.........  Deallocate local arrays
        IF( ALLOCATED( IPPTR ) )    DEALLOCATE( IPPTR )
        IF( ALLOCATED( IPMAX ) )    DEALLOCATE( IPMAX )
        IF( ALLOCATED( SRCPOL ) )   DEALLOCATE( SRCPOL )
        IF( ALLOCATED( SRCID ) )    DEALLOCATE( SRCID )
        IF( ALLOCATED( COMPUTED ) ) DEALLOCATE( COMPUTED )

        DEALLOCATE( EONAMES, EOTYPES, EOUNITS, EODESCS )

C.........  Reset the number of map records, in case some had no data
        NMAP = MCNT

C.........  Write the map inventory file
        WRITE( IDEV, '(A)' ) CATDESC
        WRITE( IDEV, '(A)' ) '/IOAPI/ ' // TRIM( NAME1 ) // '.ncf'
        WRITE( IDEV, '(A)' ) '/TEXT/ ' // TRIM( NAME2 ) // '.txt'
        WRITE( IDEV, '(A,I8)' ) '/NDAT/ ', NMAP
        WRITE( IDEV, '(A)' ) '/DATMAP/'

        DO I = 1, NMAP
C.............  Write the relative file name
            WRITE( IDEV, '(A)' ) MAPNAM( I ) // ' ' // TRIM( MAPFIL(I) )

C.............  Reset with the full path, in case needed for hour-specific 
C               or day-specific processing emissions reading
            MAPFIL( I ) = TRIM( BPATH ) // TRIM( MAPFIL(I) ) 

        END DO

        WRITE( IDEV, '(A)' ) '/END/'

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

C******************  INTERNAL SUBPROGRAMS  *****************************

        CONTAINS

C.............  This internal subprogram is for writing out the inventory
C               data, whether it is the pollutant data or the activity data.
C.............  Most variables are defined from the main subroutine
            SUBROUTINE LOOP_FOR_OUTPUT( NOUT, NPVAR, IDEV, INDX, 
     &                                  NAMES, MCNT )

C.............  Subroutine arguments 
            INTEGER     , INTENT (IN) :: NOUT          ! no. pols/act for output
            INTEGER     , INTENT (IN) :: NPVAR         ! no. vars per data
            INTEGER     , INTENT (IN) :: IDEV          ! unit no. for map file
            INTEGER     , INTENT (IN) :: INDX ( NOUT ) ! index to master list
            CHARACTER(*), INTENT (IN) :: NAMES( NOUT ) ! names of pols/act
            INTEGER , INTENT (IN OUT) :: MCNT          ! map counter

C.............  Local allocatable arrays...
C.............  Array for storing temporary variable names for computed var
            CHARACTER(IOVLEN3), ALLOCATABLE :: VNAMFORM( : )

C.............  Array for output of pol/act and assoc data
            REAL , ALLOCATABLE :: OUTVAL( :,: )

C.............  Local variables
            INTEGER    I, J, K, L, M, N, S

            INTEGER          NREC      ! number of output records for each pollutant

            REAL           MINANN      ! min negative annual value
            REAL           MINAVD      ! min negative seasonal/ave day value

            LOGICAL     :: FFLAGA = .FALSE.  ! true: first part of formula processed
            LOGICAL     :: FFLAGB = .FALSE.  ! true: 2nd part of formula processed

            CHARACTER(PHYLEN3) :: RFNAME ! relative physical file name

C----------------------------------------------------------------------

            DO I = 1, NOUT

C.................  Initialize SRCPOL pollutant/activity-specific data array
C                   as missing
                SRCPOL( :,1 ) = 0.      ! array
                IF( NDY .GT. 0 ) SRCPOL( :,NDY ) = 0.      ! array
                IF( NC1 .GT. 0 ) SRCPOL( :,NC1 ) = RIMISS3 ! array
                IF( NC2 .GT. 0 ) SRCPOL( :,NC2 ) = RIMISS3 ! array
                IF( NCE .GT. 0 ) SRCPOL( :,NCE ) = 0.      ! array
                IF( NRE .GT. 0 ) SRCPOL( :,NRE ) = 100.    ! array
                IF( NRP .GT. 0 ) SRCPOL( :,NRP ) = 100.    ! array

C.................  Transfer emissions or activity data to output SRCPOL array
                K = 0
                N = 0
                DO S = 1, NSRC

C.....................  Set position in sparse array by comparing pointer to max
                    K = MIN( IPPTR( S ), IPMAX( S ) ) 

C.....................  Retrieve emissions from the pollutant array if the 
C                       source pollutant ID equals the pollutant ID of the 
C                       pollutant of iteration I.
                    IF( IPOSCOD( K ) .EQ. INDX( I ) ) THEN

                        IPPTR( S ) = K + 1  ! pointer for source S

C......................  Skip records with zero data or missing data
                        IF( .NOT. ZFLAG ) THEN
                            IF( POLVAL( K,1 ) .LE. 0 . AND.
     &                          POLVAL( K,2 ) .LE. 0        ) CYCLE
                        END IF

                        N = N + 1
                        SRCID( N ) = S

                        DO J = 1, NPVAR     ! rearrange pollutant-specific info
                            SRCPOL( N,J ) = POLVAL( K,J )
                        END DO

C......................  If current data variable is the first variable in the
C                        formula, then store data in formula arrays
                        IF( NAMES( I ) .EQ. VIN_A ) THEN
                            FFLAGA = .TRUE.
                            COMPUTED( S,1 ) = COMPUTED( S,1 ) +         ! Annual value
     &                                        MAX( 0., SRCPOL( N,1 ) )

                            IF( NPVAR .GT. 1 ) COMPUTED(S,2)=           ! Average day value
     &                          COMPUTED(S,2) + MAX( 0., SRCPOL(N,2) )

                            IF( NPVAR .GT. 2 )                          ! Other info
     &                          COMPUTED(S,3:NPVAR) = SRCPOL(N,3:NPVAR) 

C...........................  Store variable names for computed variable
                            IF( .NOT. ALLOCATED( VNAMFORM ) ) THEN
                                ALLOCATE( VNAMFORM( NPVAR ), STAT=IOS )
                                CALL CHECKMEM( IOS,'VNAMFORM',PROGNAME )
                                VNAMFORM = ' '
                                DO J = 1, NPVAR
                                    L= LEN_TRIM( NAMES( I ) )
                                    L= INDEX(EONAMES(I,J),NAMES(I)(1:L))
                                    IF ( L .GT. 1 ) THEN
                                        VNAMFORM(J)=EONAMES(I,J)(1:L-1)
     &                                              // TRIM( VNAME )
                                    ELSE
                                        VNAMFORM( J ) = VNAME
                                    END IF
                                END DO
                            END IF  ! need to store variable name info

                        END IF  ! end if pollutant I is first var in formula

C......................  If current data variable is the second variable in the
C                        formula, then use data in formula to compute output value
                        IF( NAMES( I ) .EQ. VIN_B ) THEN

                            FFLAGB = .TRUE.
                            IF( CHKPLUS ) THEN
                                COMPUTED( S,1 )= COMPUTED( S,1 ) + 
     &                                           MAX( 0., SRCPOL(N,1) )

                                IF ( NPVAR .GT. 1 ) COMPUTED( S,2 ) = 
     &                               COMPUTED(S,2) + MAX(0.,SRCPOL(N,2))

                            ELSE IF( CHKMINUS ) THEN
                                COMPUTED( S,1 )= COMPUTED( S,1 ) - 
     &                                           MAX( 0., SRCPOL(N,1) )
                                IF ( NPVAR .GT. 1 ) COMPUTED( S,2 ) = 
     &                               COMPUTED(S,2) - MAX(0.,SRCPOL(N,2))

                            END IF  ! If formula uses plus or minus
                        END IF      ! If pol I is 2nd var in formula

                   END IF           ! If pol/act available for current source

                END DO  ! end of loop through sources
                NREC = N

C.................  Give warning if no emissions data available
                IF( NREC .EQ. 0 ) THEN                    
                    MESG = 'WARNING: no sources with data for output '//
     &                     'of data variable "' // TRIM( NAMES(I) )//'"'
     &                     //CRLF()//BLANK10// 
     &                     'File will not be written.'
                    CALL M3MSG2( MESG )
                    CYCLE
                END IF

C.................  Give warning about maximum value of negative numbers

C.................  Set up to open output file...
C.................  Set the pollutant-specific items for the I/O API header
                NROWS3D = NREC

                M = 1
                DO J = 1, NPVAR
                    M = M + 1
                    VNAMESET( M ) = EONAMES( I,J )
                    VTYPESET( M ) = EOTYPES( I,J )
                    VUNITSET( M ) = EOUNITS( I,J )
                    VDESCSET( M ) = EODESCS( I,J )
                END DO    

C.................  Reset units for the primary (annual) data value
                K = INDEX1( EONAMES( I,1 ), MXIDAT, INVDNAM )
                VUNITSET( 2 ) = INVDUNT( K )

C.................  Allocate memory for output arrays
                ALLOCATE( OUTVAL( NREC, NPVAR ), STAT=IOS )
                CALL CHECKMEM( IOS, 'OUTVAL', PROGNAME )

                DO J = 1, NPVAR
                    OUTVAL( 1:NREC,J ) = SRCPOL( 1:NREC, J )
                END DO

C.................  Build name and output pollutant or activity file
                RFNAME = TRIM( RPATH )// '/'// TRIM( NAMES(I) )// '.ncf'
                CALL WRINVPOL( CATEGORY, BPATH, RFNAME, NREC, NPVAR, 
     &                         NAMES(I), SRCID(1:NREC), OUTVAL, EFLAG )

                DEALLOCATE( OUTVAL )

C...............  Update map file counter and map file names
                MCNT = MCNT + 1
                MAPNAM( MCNT ) = NAMES( I )
                MAPFIL( MCNT ) = TRIM( RFNAME )

            END DO  ! end loop I through output variables

C...........  If formula has been finished on this iteration, then
C             also output the emissions from the formula
            IF( FFLAGA .AND. FFLAGB ) THEN

C...............  Reset flags so that the information won't be written again
                FFLAGA = .FALSE.
                FFLAGB = .FALSE.

C...............  Count up the number of non-zero records for the computed
C                 emissions.  No need to reformat SRCPOL or SRCID, because it is 
C                 sparse storage.  Populate sparse arrays.
                NREC = 0
                MINANN = 0.
                MINAVD = 0.
                DO S = 1, NSRC
                
                    IF( FIREFLAG ) ZFLAG = .TRUE.    ! disregard zero values on PMC(wildfire only)

C.....................  Skip records with zero data
                    IF( .NOT. ZFLAG ) THEN
                        IF( COMPUTED( S,1 ) == 0 .AND.
     &                      COMPUTED( S,2 ) == 0       ) CYCLE
                    END IF

                    NREC = NREC + 1
                    SRCID ( NREC ) = S

C......................  Check for negative values for annual value
                    IF( COMPUTED( S,1 ) .LT. 0 ) THEN
                        L = LEN_TRIM( VNAMFORM( 1 ) )
                        CALL FMTCSRC( CSOURC( S ), 7, BUFFER, L2 )
                        WRITE( MESG,94020 ) 
     &                    'WARNING: Resetting negative value of "'//
     &                    VNAMFORM(1)(1:L)//'" from', COMPUTED(S,1),
     &                    'to 0. for source:'// CRLF() // BLANK10// 
     &                    BUFFER( 1:L2 )
                        CALL M3MESG( MESG )
                        MINANN = MIN( MINANN, COMPUTED( S,1 ) )
                        COMPUTED( S,1 ) = 0.
                    END IF

C......................  Check for negative values for average day value
                    IF( COMPUTED( S,2 ) .LT. 0 ) THEN
                        L = LEN_TRIM( VNAMFORM( 2 ) )
                        CALL FMTCSRC( CSOURC( S ), 7, BUFFER, L2 )
                        WRITE( MESG,94020 ) 
     &                    'WARNING: Resetting negative value of "'//
     &                    VNAMFORM(2)(1:L)//'" from', COMPUTED(S,2),
     &                    'to 0. for source:'// CRLF() // BLANK10// 
     &                    BUFFER( 1:L2 )
                        CALL M3MESG( MESG )
                        MINAVD = MIN( MINAVD, COMPUTED( S,2 ) )
                        COMPUTED( S,2 ) = 0.
                    END IF
                    
                    SRCPOL( NREC,1:NPVAR ) = COMPUTED( S,1:NPVAR )

                END DO          ! End loop through sources for sparse storage

C...............  Give warning for that includes negative annual value
                IF( MINANN .LT. 0. ) THEN
                    WRITE( MESG,94020 ) 'WARNING: Largest negative '//
     &                     'value for "'// TRIM( VNAMFORM(1) ) //
     &                     '" was', MINANN, 'and was reset to 0.'
                    CALL M3MSG2( MESG )
                END IF

C...............  Give warning for that includes negative seasonal value
                IF( MINAVD .LT. 0. ) THEN
                    WRITE( MESG,94020 ) 'WARNING: Largest negative '//
     &                     'value for "'// TRIM( VNAMFORM(2) ) //
     &                     '" was', MINAVD, 'and was reset to 0.'
                    CALL M3MSG2( MESG )
                END IF

C...............  Set up output arrays in sparse format....
C...............  Reset I/O API header information
                NROWS3D = NREC

C...............  For header variable info, rely on the fact that all pollutants
C                 have the SRCID in position 1 and the same types and units.
                VNAMESET( 2:NPVAR+1 ) = VNAMFORM( 1:NPVAR )
          
C...............  Build name and output pollutant file
                RFNAME = TRIM( RPATH )// '/'// TRIM( VNAME )// '.ncf'
                CALL WRINVPOL( CATEGORY, BPATH, RFNAME, NSRC, NPVAR, 
     &                         VNAME, SRCID, SRCPOL, EFLAG )

C...............  Update map file counter and map relative file names
                MCNT = MCNT + 1
                MAPNAM( MCNT ) = VNAME
                MAPFIL( MCNT ) = TRIM( RFNAME )

            END IF  ! if formula is ready to be written

            IF( EFLAG ) THEN
                MESG = 'Problem writing data.'
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF

            RETURN

C----------------------  FORMAT  STATEMENTS   --------------------------

C...........   Internal buffering formats............ 94xxx

94020       FORMAT( 10( A, :, E10.2, :, 1X ) )

            END SUBROUTINE LOOP_FOR_OUTPUT
  
        END SUBROUTINE WRINVEMIS
