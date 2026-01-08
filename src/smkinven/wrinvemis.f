
        SUBROUTINE WRINVEMIS( IDEV, DATPATH )

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
C     09/2025 by HT UNC-IE:  Use M3UTILIO
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
        USE M3UTILIO

C.........  MODULES for public variables
C...........   This module is the inventory arrays
        USE MODSOURC, ONLY: CSOURC, NPCNT, IPOSCOD, POLVAL

C.........  This module contains the lists of unique inventory information
        USE MODLISTS, ONLY: MXIDAT, INVDNAM, INVDUNT, FIREFLAG, NINVTBL,
     &                      ITNAMA, ITCASA, FF10FLAG, MEDSFLAG

C.........  This module contains the information about the source category
        USE MODINFO, ONLY: CATEGORY, CATDESC, NSRC, NMAP,
     &                     NIPOL, NIACT, NIPPA, NPPOL, EANAM, ACTVTY,
     &                     MAPNAM, MAPFIL, EINAM, EIIDX, NPACT, AVIDX,
     &                     NDY, NC1, NC2, NCE, NRE, NRP, VAR_FORMULA,
     &                     NCOMP, CHKPLUS, CHKMINUS, FORMULAS, VIN_A,
     &                     VIN_B, VNAME

C.........  This module is required by the FileSetAPI
        USE MODFILESET
        
       IMPLICIT NONE

C...........   INCLUDES

        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
c       INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
        INCLUDE 'SETDECL.EXT'   !  FileSetAPI variables and functions

C.........  EXTERNAL FUNCTIONS
c       CHARACTER(2)  CRLF
c       INTEGER       ENVINT
c       LOGICAL       ENVYN
c       INTEGER       INDEX1
c       LOGICAL       SETENVVAR

c       EXTERNAL      CRLF, ENVINT, ENVYN, INDEX1, SETENVVAR

C.........  SUBROUTINE ARGUMENTS
        INTEGER     , INTENT (IN) :: IDEV          ! unit number for map file
        CHARACTER(PHYLEN3), INTENT (IN) :: DATPATH ! path for pol/act files

C...........   LOCAL PARAMETERS
        CHARACTER(16), PARAMETER :: FORMEVNM = 'SMKINVEN_FORMULA'

C.........  Inventory temporary arrays
        INTEGER, ALLOCATABLE:: IPPTR ( : ) ! position in POLVAL sparse array (for input pols/actv)
        INTEGER, ALLOCATABLE:: IPPTR2( : ) ! position in POLVAL sparse array (for computed values)
        INTEGER, ALLOCATABLE:: IPMAX ( : ) ! max IPPTR by source
        INTEGER, ALLOCATABLE:: SRCID ( : ) ! source index
        REAL   , ALLOCATABLE:: SRCPOL( :,: )  ! data-spec values by source

C...........   Names, Units, types, & descriptions for pollutant-specific 
C              output variables.  NOTE - second dimension will work so long
C              as NPPOL > NPACT, which is expected to always be the case

        CHARACTER(IOVLEN3), ALLOCATABLE :: EONAMES( :,: ) ! Names 
        INTEGER           , ALLOCATABLE :: EOTYPES( :,: ) ! Types (Real|Int)
        CHARACTER(IOULEN3), ALLOCATABLE :: EOUNITS( :,: ) ! Units  
        CHARACTER(IODLEN3), ALLOCATABLE :: EODESCS( :,: ) ! Dscriptions  

C...........   Other local allocatable arrays
        CHARACTER(IOVLEN3), ALLOCATABLE :: SAVEANAM( : ) ! tmp variables

C...........   Other local variables
        INTEGER         F, I, J, S, L, L2, N     ! counters and indices

        INTEGER         IOS       ! i/o status

        INTEGER         MCNT      ! count of actual mapped pol/act files
        INTEGER         MXWARN    ! maximum number of warnings of each type to write

        INTEGER         RIMISS3          ! real value of integer missing

        LOGICAL      :: EFLAG    = .FALSE. ! true: error found
        LOGICAL      :: FFLAG    = .FALSE. ! true: formula in use
        LOGICAL      :: NEGOK    = .FALSE. ! true: okay to output negative emission values
        LOGICAL         ZFLAG              ! true: write zeros to output file

        CHARACTER(80)   NAME1            ! tmp file name component
        CHARACTER(80)   NAME2            ! tmp file name component
        CHARACTER(128)  BUFFER           ! message buffer
        CHARACTER(256)  MESG             ! message buffer

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

        MXWARN  = ENVINT( WARNSET  , ' ', 100, IOS )

        NEGOK = ENVYN( 'ALLOW_NEGATIVE',
     &                 'Allow negative output data',
     &                 .FALSE., IOS )

        IF( EFLAG ) THEN
            MESG = 'Problem with input environment variables'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

C.........  If processing fires or FF10 formats, then must write out zero values, since
C           zero annual will be filled with daliy/hourly fire/F10 inventories 
        IF ( FIREFLAG .OR. FF10FLAG .OR. MEDSFLAG ) THEN
            ZFLAG = .TRUE.

C.........  Otherwise, check environment to see if user wants to output
C           zero values
        ELSE
            MESG = 'Write zero values to annual emissions inventory'
            ZFLAG = ENVYN( 'WRITE_ANN_ZERO', MESG, .FALSE., IOS )
        END IF

C.........  Compute real value of integer missing
        RIMISS3 = REAL( IMISS3 )

C.........  Set maximum number of map variables for map-formatted outputs
        NMAP = NIPOL + NIACT
        MCNT = 0              ! initialize actual pol/act file count for later

C.........  If there is one or more computed output variable, get set up
        IF( LEN_TRIM( VAR_FORMULA ) .GT. 0 ) THEN
            CALL FORMLIST
            FFLAG = .TRUE.
            NMAP = NMAP + NCOMP  ! NCOMP more variable(s) to map
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

C.........  Allocate memory for indices IPPTR & IPMAX for pointing to position
C           in sparsely stored data array.  
        ALLOCATE( IPPTR( NSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'IPPTR', PROGNAME )
        ALLOCATE( IPPTR2( NSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'IPPTR2', PROGNAME )
        ALLOCATE( IPMAX( NSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'IPMAX', PROGNAME )

C.........  Allocate memory for storing and writing emissions
        ALLOCATE( SRCPOL( NSRC, NPPOL ), STAT=IOS )
        CALL CHECKMEM( IOS, 'SRCPOL', PROGNAME )
        ALLOCATE( SRCID( NSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'SRCID', PROGNAME )

C.........  Loop through non-formula pollutants, store, and write to inventory file
        IF( NIPOL .GT. 0 ) 
     &      CALL LOOP_FOR_OUTPUT( NIPOL, NPPOL, IDEV, EIIDX, 
     &                            EINAM, MCNT )

C.........  Loop through formula-based pollutants, compute, and write to inventory file
        IF( NCOMP .GT. 0 ) 
     &      CALL OUTPUT_FORMULAS( NIPOL, NPPOL, IDEV, EIIDX, 
     &                            EINAM, MCNT )

C.........  If computed variable, update EANAM, in case needed for day- 
C           and hour-specific data
        IF( FFLAG ) THEN

            ALLOCATE( SAVEANAM( NIPPA ), STAT=IOS )
            CALL CHECKMEM( IOS, 'SAVEANAM', PROGNAME )
            SAVEANAM = EANAM       ! array

            DEALLOCATE( EANAM )
            ALLOCATE( EANAM( NIPPA+NCOMP ), STAT=IOS )
            CALL CHECKMEM( IOS, 'EANAM', PROGNAME )

C.............  Save names of existing pollutants
            EANAM( 1:NIPOL ) = SAVEANAM( 1:NIPOL )

C.............  Add names of computed pollutants
            DO F = 1, NCOMP
                EANAM( NIPOL+F ) = VNAME( F )
            END DO

C.............  Save names of existing activities
            EANAM( NIPOL+NCOMP+1:NIPPA+NCOMP )=SAVEANAM( NIPOL+1:NIPPA )

C.............  Update pollutant/activity counts
            NIPOL = NIPOL + NCOMP
            NIPPA = NIPPA + NCOMP

C.............  Release memory for temporary array
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
        IF ( CATEGORY .EQ. 'POINT' )  NVARSET = NVARSET + 4     !!  FUG_* variables...
        WRITE( FDESC3D( 1 ), '(A,1X,I8)' ) '/NSRC/', NSRC
        
        IF( ALLOCATED( VARS_PER_FILE ) ) DEALLOCATE( VARS_PER_FILE )

C.........  Loop through activity data, store, and write to inventory file
        IF( NIACT .GT. 0 ) 
     &      CALL LOOP_FOR_OUTPUT( NIACT, NPACT, IDEV, AVIDX, 
     &                            ACTVTY, MCNT )
C.........  Deallocate local arrays
        IF( ALLOCATED( IPPTR ) )    DEALLOCATE( IPPTR )
        IF( ALLOCATED( IPPTR2 ) )   DEALLOCATE( IPPTR2 )
        IF( ALLOCATED( IPMAX ) )    DEALLOCATE( IPMAX )
        IF( ALLOCATED( SRCPOL ) )   DEALLOCATE( SRCPOL )
        IF( ALLOCATED( SRCID ) )    DEALLOCATE( SRCID )

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

C.............  Local variables
            INTEGER    I, J, K, L, M

            INTEGER          NREC      ! number of output records for each pollutant

            CHARACTER(PHYLEN3) :: RFNAME ! relative physical file name

C----------------------------------------------------------------------

C.............  Loop through output variables (pollutants or activities)
            DO I = 1, NOUT

                NREC = 0
                CALL FILL_OUTPUT_DATA( INDX(I), NPVAR, 'ORDERED', NREC )

C.................  Give warning if no emissions data available
                IF( NREC .EQ. 0 ) THEN                    
                    MESG = 'WARNING: no sources with data for output '//
     &                     'of data variable "' // TRIM( NAMES(I) )//'"'
     &                     //CRLF()//BLANK10// 
     &                     'File will not be written.'
                    CALL M3MSG2( MESG )
                    CYCLE
                END IF

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
                IF( K .LE. 0 ) THEN
                    MESG = 'ERROR: Cannot find pollutant "'//
     &                     TRIM( EONAMES( I,1 ) ) // '" from '//
     &                     'inventory in INVTABLE file.'
                    CALL M3MSG2( MESG )
                    EFLAG = .TRUE.
                    CYCLE
                END IF

                VUNITSET( 2 ) = INVDUNT( K )

C.................  Build name and output pollutant or activity file
                RFNAME = TRIM( RPATH )// '/'// TRIM( NAMES(I) )// '.ncf'

                CALL WRINVPOL( CATEGORY, BPATH, RFNAME, NSRC, NPVAR, 
     &                         NAMES(I), SRCID, SRCPOL, EFLAG )

C...............  Update map file counter and map file names
                MCNT = MCNT + 1
                MAPNAM( MCNT ) = NAMES( I )
                MAPFIL( MCNT ) = TRIM( RFNAME )

            END DO  ! end loop I through output variables

C.............  Exit program
            IF( EFLAG ) THEN
                MESG = 'Problem writing data.'
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF

            RETURN

C----------------------  FORMAT  STATEMENTS   --------------------------

C...........   Internal buffering formats............ 94xxx

94020       FORMAT( 10( A, :, E10.2, :, 1X ) )

            END SUBROUTINE LOOP_FOR_OUTPUT

C-----------------------------------------------------------------------
C-----------------------------------------------------------------------

C.............  This internal subprogram is for writing out the computed formulas' 
C               inventory data.
C.............  Most variables are defined from the main subroutine
            SUBROUTINE OUTPUT_FORMULAS( NOUT, NPVAR, IDEV, INDX, 
     &                                  NAMES, MCNT )

C.............  Subroutine arguments 
            INTEGER     , INTENT (IN) :: NOUT          ! no. pols/act for output
            INTEGER     , INTENT (IN) :: NPVAR         ! no. vars per data
            INTEGER     , INTENT (IN) :: IDEV          ! unit no. for map file
            INTEGER     , INTENT (IN) :: INDX ( NOUT ) ! index to master list
            CHARACTER(*), INTENT (IN) :: NAMES( NOUT ) ! names of pols/act
            INTEGER , INTENT (IN OUT) :: MCNT          ! map counter

C.............  Arrays for emissions or activities output
            REAL   , ALLOCATABLE:: COMPUTED( :,: ) ! computed data-spec values by src

C.............  Local variables
            INTEGER    F, J, K, L, M, N, S

            INTEGER          IDXA      ! position of first formula variable in list
            INTEGER          IDXB      ! position of 2nd formula variable in list
            INTEGER          NREC      ! number of output records for each pollutant
            INTEGER          WARNCNT_A ! number of times output warning A 
            INTEGER          WARNCNT_B ! number of times output warning B 

            REAL           MINANN      ! min negative annual value
            REAL           MINAVD      ! min negative seasonal/ave day value

            LOGICAL     :: FFLAGA = .FALSE.  ! true: first part of formula processed
            LOGICAL     :: FFLAGB = .FALSE.  ! true: 2nd part of formula processed

            CHARACTER*256      :: BUFFER
            CHARACTER(PHYLEN3) :: RFNAME ! relative physical file name


C----------------------------------------------------------------------

C.............  Write message saying that formula variables are being computed
            WRITE( MESG,94010 ) 'NOTE: Computing variables for',
     &                           NCOMP, 'formulas.'
            CALL M3MSG2( MESG )

C.............  Allocate full-source sized array for performing emissions computations
            ALLOCATE( COMPUTED( NSRC, NPTPPOL3 ), STAT=IOS )
            CALL CHECKMEM( IOS, 'COMPUTED', PROGNAME )

C.............  Deallocate arrays for variable names for output to SMOKE intermediate files
            DEALLOCATE( EONAMES, EOUNITS, EOTYPES, EODESCS )

C.............  Allocate memory for temporary variable names etc.
            ALLOCATE( EONAMES( NCOMP,NPVAR ), STAT=IOS )
            CALL CHECKMEM( IOS, 'EONAMES', PROGNAME )
            ALLOCATE( EOUNITS( NCOMP,NPVAR ), STAT=IOS )
            CALL CHECKMEM( IOS, 'EOUNITS', PROGNAME )
            ALLOCATE( EOTYPES( NCOMP,NPVAR ), STAT=IOS )
            CALL CHECKMEM( IOS, 'EOTYPES', PROGNAME )
            ALLOCATE( EODESCS( NCOMP,NPVAR ), STAT=IOS )
            CALL CHECKMEM( IOS, 'EODESCS', PROGNAME )
        
C.............  Get names, units, etc. of output activity-specific records
            CALL BLDENAMS( CATEGORY, NCOMP, NPVAR, VNAME, 
     &                     EONAMES, EOUNITS, EOTYPES, EODESCS )

C.............  Loop through formulas. Note that the variables have already
C               been checked in the main routine to make sure the names match if
C               the code has gotten this far, so no additional checks on the INDEX1
C               calls included here.
            DO F = 1, NCOMP

C.................  Write message about current formula
                MESG = 'NOTE: Computing variable "'//TRIM( VNAME(F) )// 
     &                 '" using formula "'// TRIM( FORMULAS(F) )// '"'
                CALL M3MSG2( MESG )

C.................  Initialize values for minimum negative values
                MINANN = 0.
                MINAVD = 0.

C.................  Initialize values for warning counts (once per formula)
                WARNCNT_A = 0
                WARNCNT_B = 0

C.................  Initialize array to account for sparse storage when filling in
                COMPUTED = 0.  ! array

C.................  Find first variable in formula in list of pollutants
                IDXA = INDEX1( VIN_A( F ), NOUT, NAMES )

C.................  Populate SRCPOL array from first input variable
                NREC = 0
                CALL FILL_OUTPUT_DATA( INDX(IDXA),NPVAR,'SEARCH',NREC )

C.................  Give error if no emissions data available
                IF( NREC .EQ. 0 ) THEN                    
                    EFLAG = .TRUE.
                    MESG = 'ERROR: no sources with data for output '//
     &                     'of data variable "'// TRIM(NAMES(IDXA))//'"'
     &                     //CRLF()//BLANK10// 
     &                     'Variable "'//TRIM(VNAME(F))//'" will '//
     &                     'not be computed.'
                    CALL M3MSG2( MESG )
                    CYCLE
                END IF

C.................  Initialize full-source computational array with first variable
                DO N = 1, NREC
                    S = SRCID( N )
                    COMPUTED( S,1 ) = MAX( 0., SRCPOL( N,1 ) )

                    IF( NPVAR .GT. 1 ) COMPUTED( S,2 ) = 
     &                                 MAX( 0., SRCPOL( N,2 ) ) ! Average day value

                    IF( NPVAR .GT. 2 )                                ! Other info
     &                  COMPUTED( S,3:NPVAR ) = SRCPOL( N,3:NPVAR ) 

                END DO

C.................  Find second variable in formula in list of pollutants
                IDXB = INDEX1( VIN_B( F ), NOUT, NAMES )

C.................  Populate SRCPOL array with second input variable
                NREC = 0
                CALL FILL_OUTPUT_DATA( INDX(IDXB),NPVAR,'SEARCH',NREC )
               
C.................  Give error if no emissions data available
                IF( NREC .EQ. 0 ) THEN                    
                    EFLAG = .TRUE.
                    MESG = 'ERROR: no sources with data for output '//
     &                     'of data variable "'// TRIM(NAMES(IDXB))//'"'
     &                     //CRLF()//BLANK10// 
     &                     'Variable "'//TRIM(VNAME(F))//'" will '//
     &                     'not be computed.'
                    CALL M3MSG2( MESG )
                    CYCLE
                END IF

C.................  Compute variable by adding or subtracting second variable
                DO N = 1, NREC
                    S = SRCID( N )

                    IF( CHKPLUS( F ) ) THEN     ! Formula is an addition

                        COMPUTED( S,1 ) = COMPUTED( S,1 ) + 
     &                                    MAX( 0., SRCPOL( N,1 ) )

                        IF ( NPVAR .GT. 1 ) 
     &                       COMPUTED( S,2 ) = 
     &                       COMPUTED( S,2 ) + MAX( 0.,SRCPOL( N,2 ) )

                    ELSE IF( CHKMINUS( F ) ) THEN    ! Formula is a subtraction

                        COMPUTED( S,1 ) = COMPUTED( S,1 ) - 
     &                                    MAX( 0., SRCPOL( N,1 ) )

                        IF ( NPVAR .GT. 1 ) 
     &                       COMPUTED( S,2 ) = 
     &                       COMPUTED( S,2 ) - MAX( 0.,SRCPOL( N,2 ) )

C.........................  Check to see if the computed value is now negative, and
C                           if so, reset to zero.
                        IF( COMPUTED( S,1 ) .LT. 0 ) THEN
                            WARNCNT_A = WARNCNT_A + 1

C.............................  If warning count is less than max and not
C                               fires.
                            IF ( WARNCNT_A .LE. MXWARN .AND.
     &                           .NOT. FIREFLAG ) THEN

                              CALL FMTCSRC( CSOURC( S ), 7, BUFFER, L2 )

                              IF( NEGOK ) THEN

                                WRITE( MESG,94020 ) 'WARNING: '//
     &                          'Retaining negative value of annual "'//
     &                          TRIM(VNAME(F))//'" of', COMPUTED(S,1),
     &                          ' for source:'//CRLF()//BLANK10// 
     &                          BUFFER( 1:L2 )

                              ELSE

                                WRITE( MESG,94020 ) 'WARNING: '//
     &                          'Resetting negative value of annual "'//
     &                          TRIM(VNAME(F))//'" from', COMPUTED(S,1),
     &                          'to 0. for source:'//CRLF()//BLANK10// 
     &                          BUFFER( 1:L2 )

                              END IF

                              CALL M3MESG( MESG )
                            END IF

                            MINANN = MIN( MINANN, COMPUTED( S,1 ) )
                            IF( .NOT. NEGOK ) COMPUTED( S,1 ) = 0.

                        END IF

C..........................  Check for negative values for average day value
                        IF( COMPUTED( S,2 ) .LT. 0 ) THEN
                            WARNCNT_B = WARNCNT_B + 1

                            IF ( WARNCNT_B .LE. MXWARN ) THEN
                              CALL FMTCSRC( CSOURC( S ), 7, BUFFER, L2 )

                              IF( NEGOK ) THEN
                                WRITE( MESG,94020 ) 'WARNING: '//
     &                           'Retaining negative value of "'//
     &                           'average-day "'//TRIM(VNAME(F))// 
     &                           '" of',COMPUTED(S,2),'for source:'//
     &                           CRLF()//BLANK10//BUFFER(1:L2)
                              ELSE
                                WRITE( MESG,94020 ) 'WARNING: '//
     &                           'Resetting negative value of "'//
     &                           'average-day "'//TRIM(VNAME(F))// 
     &                           '" from',COMPUTED(S,2),'to 0. for '//
     &                          'source:'//CRLF()//BLANK10//BUFFER(1:L2)
                              END IF 
                              CALL M3MESG( MESG )
                            END IF

                            MINAVD = MIN( MINAVD, COMPUTED( S,2 ) )
                            IF ( .NOT. NEGOK ) COMPUTED( S,2 ) = 0.

                        END IF

                    END IF   ! If formula is a plus or a minus

                END DO       ! End loop over values in second variable

C................   Copy data to output structure (no condensing output structure)
                IF ( ZFLAG ) THEN
                    DO S = 1, NSRC
                        SRCID( S ) = S
                    END DO

                    SRCPOL( 1:NSRC,1:NPVAR ) = COMPUTED( 1:NSRC,1:NPVAR)   ! array
                    NREC = NSRC

C................   Condense data to output structure by excluding zeros
                ELSE
                    N = 0
                    DO S = 1, NSRC

C.........................  Skip records with zero data, unless option (ZFLAG) says not to
                        IF( COMPUTED( S,1 ) == 0 .AND.
     &                      COMPUTED( S,2 ) == 0       ) CYCLE

                        N = N + 1
                        SRCID ( N ) = S                
                        SRCPOL( N,1:NPVAR ) = COMPUTED( S,1:NPVAR )

                    END DO          ! End loop through sources for sparse storage
                    NREC = N

                END IF

C.................  Give an error when all computed values are zero
                IF( NREC == 0 ) THEN
                    MESG = 'WARNING: All computed values are zero. Nothing to output'
                    CALL M3MESG( MESG )
                    CYCLE
                END IF

C...............  Give warning for that includes negative annual value
                IF( MINANN .LT. 0. ) THEN
                    WRITE( MESG,94020 ) 'WARNING: Largest negative '//
     &                     'value for "'// TRIM( EONAMES(F,1) ) //
     &                     '" was', MINANN
                    IF( .NOT. NEGOK )
     &                  MESG = TRIM( MESG ) // ' and was reset to 0.'

                    CALL M3MSG2( MESG )
                END IF

C...............  Give warning for that includes negative seasonal value
                IF( MINAVD .LT. 0. ) THEN
                    WRITE( MESG,94020 ) 'WARNING: Largest negative '//
     &                     'value for "'// TRIM( EONAMES(F,2) ) //
     &                     '" was', MINAVD
                    IF( .NOT. NEGOK )
     &                  MESG = TRIM( MESG ) // ' and was reset to 0.'

                    CALL M3MSG2( MESG )
                END IF

C.................  Set up to open output file...
C.................  Set the pollutant-specific items for the I/O API header
                NROWS3D = NREC

                M = 1
                DO J = 1, NPVAR
                    M = M + 1
                    VNAMESET( M ) = EONAMES( F,J )
                    VTYPESET( M ) = EOTYPES( F,J )
                    VUNITSET( M ) = EOUNITS( F,J )
                    VDESCSET( M ) = EODESCS( F,J )
                END DO    

C.................  Reset units for the primary (annual) data value
                K = INDEX1( EONAMES( F,1 ), MXIDAT, INVDNAM )
                IF( K .LE. 0 ) THEN
                    MESG = 'ERROR: Cannot find pollutant "'//
     &                     TRIM( EONAMES( F,1 ) ) // '" from '//
     &                     'inventory in INVTABLE file.'
                    CALL M3MSG2( MESG )
                    EFLAG = .TRUE.
                    CYCLE
                END IF

                VUNITSET( 2 ) = INVDUNT( K )

C...............  Build name and output pollutant file
                RFNAME = TRIM( RPATH )// '/'// TRIM( VNAME(F) )// '.ncf'
                CALL WRINVPOL( CATEGORY, BPATH, RFNAME, NSRC, NPVAR, 
     &                         VNAME(F), SRCID, SRCPOL, EFLAG )

C...............  Update map file counter and map relative file names
                MCNT = MCNT + 1
                MAPNAM( MCNT ) = VNAME( F )
                MAPFIL( MCNT ) = TRIM( RFNAME )

            END DO   !  End of formulas, loop on F

C.............  Deallocate local memory
            IF( ALLOCATED( COMPUTED ) ) DEALLOCATE( COMPUTED )

C.............  Exit if error found (all formulas couldn't be computed)
            IF ( EFLAG ) THEN
                MESG = 'Could not output all formulas. ' //
     &                 'See previous error messages.'
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF

            RETURN

C......................  FORMAT  STATEMENTS   ..........................

C...............   Internal buffering formats............ 94xxx

94010       FORMAT( 10( A, :, I8, :, 1X ) )

94020       FORMAT( 10( A, :, E10.2, :, 1X ) )

            END SUBROUTINE OUTPUT_FORMULAS

C-----------------------------------------------------------------------
C-----------------------------------------------------------------------

C.............  This internal subprogram is for writing out the computed formulas' 
C               inventory data.
C.............  Most variables are defined from the main subroutine
            SUBROUTINE FILL_OUTPUT_DATA( POLINDEX, NPVAR, 
     &                                   CALLTYPE, NREC   )

C.............  Subroutine arguments 
            INTEGER     , INTENT(IN) :: POLINDEX ! index of pollutant or activity to fill
            INTEGER     , INTENT(IN) :: NPVAR    ! no. variables per pol/act
            CHARACTER(*), INTENT(IN) :: CALLTYPE ! ordered or search
            INTEGER     , INTENT(OUT):: NREC     ! number of non-zero records in array for current pollutant

C.............  Local variables
            INTEGER    I, J, K, N, S

            LOGICAL       :: FFLAG             ! pollutant has been found
            LOGICAL, SAVE :: FIRSTIME = .TRUE. ! first time called

            CHARACTER*256    MESG  ! error message buffer

C----------------------------------------------------------------------

C.............  Since the emission values are already sorted in the output order of
C               the pollutants/activities, can use the IPPTR and IPMAX indices 
C               to keep track of what pollutant/activity we are on for each source
C               (sources can have different numbers and types of output pollutants
C               and activities).
C.............  IPPTR2 has been added for use in "search" mode.  It is never
C               reset from its initial value, unlike IPPTR.  Need to have both
C               because IPPTR is used for the input pollutants and activities,
C               whereas IPPTR2 is used for just the computed pollutants. Since
C               the computed pollutants are inserted between the input pollutants
C               and the activities, it's necessary to have both arrays
            IF ( FIRSTIME ) THEN

                IPPTR ( 1 ) = 1
                IPPTR2( 1 ) = 1
                IPMAX( 1 ) = NPCNT( 1 )
                DO S = 2, NSRC
                    IPPTR ( S ) = IPPTR( S-1 ) + NPCNT( S-1 )
                    IPPTR2( S ) = IPPTR( S )
                    IPMAX ( S ) = IPPTR( S )   + NPCNT( S ) - 1
                END DO
                FIRSTIME = .FALSE.

            END IF

C.............  Initialize SRCPOL pollutant/activity-specific data array
C               as missing
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

                FFLAG = .FALSE.
                SELECT CASE( CALLTYPE )
                CASE ( 'ORDERED' )

C.....................  Set position in sparse array by comparing pointer to max
                    K = MIN( IPPTR( S ), IPMAX( S ) ) 
                    IF( IPOSCOD( K ) .EQ. POLINDEX ) THEN
                        IPPTR( S ) = K + 1  ! pointer for source S
                        FFLAG = .TRUE.
                    END IF

                CASE ( 'SEARCH' )

                    K = MIN( IPPTR2( S ), IPMAX( S ) ) 
                    DO J = 1, NPCNT( S )

                        IF( IPOSCOD( K ) .EQ. POLINDEX ) THEN
                            FFLAG = .TRUE.
                            EXIT            ! pollutant has been found, break out of loop
                        ELSE
                            K = MIN( K + 1, IPMAX( S ) ) 

                        END IF

                    END DO

                END SELECT

C.....................  Retrieve emissions from the pollutant array if the 
C                       source pollutant ID equals the pollutant ID of the 
C                       pollutant of iteration I.
                IF( FFLAG ) THEN

C......................  Skip records with zero data or missing data
                    IF( .NOT. ZFLAG ) THEN
                        IF( POLVAL( K,1 ) .LE. 0 . AND.
     &                      POLVAL( K,2 ) .LE. 0        ) CYCLE
                    END IF

                    N = N + 1
                    SRCID( N ) = S

                    DO J = 1, NPVAR     ! rearrange pollutant-specific info
                        SRCPOL( N,J ) = POLVAL( K,J )
                    END DO

                END IF           ! If pol/act available for current source

            END DO  ! end of loop through sources
            NREC = N

            RETURN

            END SUBROUTINE FILL_OUTPUT_DATA

        END SUBROUTINE WRINVEMIS
