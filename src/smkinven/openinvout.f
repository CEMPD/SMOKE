
        SUBROUTINE OPENINVOUT( GRDNM, ENAME, ANAME, SDEV, A2PFLAG )

C*************************************************************************
C  subroutine body starts at line 119
C
C  DESCRIPTION:
C      This subroutine sets up the header and variables for the I/O API 
C      inventory file, and opens the I/O API and ASCII files for the SMOKE
C      area, mobile, or point source inventory.
C
C  PRECONDITIONS REQUIRED:
C      Correct number of pollutants and names EANAM are set
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C      Subroutines: I/O API subroutines, BLDENAMS
C      Functions: I/O API functions, VERCHAR
C
C  REVISION  HISTORY:
C      Created 4/99 by M. Houyoux
C
C*************************************************************************
C
C Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
C                System
C File: @(#)$Id$
C
C COPYRIGHT (C) 2002, MCNC Environmental Modeling Center
C All Rights Reserved
C
C See file COPYRIGHT for conditions of use.
C
C Environmental Modeling Center
C MCNC
C P.O. Box 12889
C Research Triangle Park, NC  27709-2889
C
C smoke@emc.mcnc.org
C
C Pathname: $Source$
C Last updated: $Date$ 
C
C***************************************************************************

C.........  MODULES for public variables
C.........  This module contains the lists of unique source characteristics
        USE MODLISTS

C.........  This module contains the information about the source category
        USE MODINFO

C.........  This module is required by the FileSetAPI
        USE MODFILESET
        
        IMPLICIT NONE

C...........   INCLUDES
        INCLUDE 'EMCNST3.EXT'   !  emissions constat parameters
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
        INCLUDE 'SETDECL.EXT'   !  FileSetAPI function declarations

C...........   EXTERNAL FUNCTIONS and their descriptionsNRAWIN
        CHARACTER*2     CRLF
        INTEGER         ENVINT
        INTEGER         INDEX1
        INTEGER         PROMPTFFILE
        CHARACTER*16    VERCHAR

        EXTERNAL CRLF, ENVINT, INDEX1, PROMPTFFILE, VERCHAR

C...........   SUBROUTINE ARGUMENTS
        CHARACTER(*), INTENT(IN)  :: GRDNM   ! grid name if any gridded data
        CHARACTER(*), INTENT(OUT) :: ENAME   ! emis i/o api inven logical name
        CHARACTER(*), INTENT(OUT) :: ANAME   ! emis ASCII inven logical name
        INTEGER     , INTENT(OUT) :: SDEV    ! ascii output inven file unit no.
        LOGICAL     , INTENT(IN)  :: A2PFLAG ! true: using area-to-point processing

C...........   LOCAL PARAMETERS
        CHARACTER*16, PARAMETER :: FORMEVNM = 'SMKINVEN_FORMULA'
        CHARACTER*50, PARAMETER :: CVSW = '$Name$' ! CVS release tag

C...........   Names, Units, types, & descriptions for pollutant-specific 
C              output variables.  NOTE - second dimension will work so long
C              as NPPOL > NPACT, which is expected to always be the case

        CHARACTER(LEN=IOVLEN3), ALLOCATABLE :: EONAMES( :,: ) ! Names 
        INTEGER               , ALLOCATABLE :: EOTYPES( :,: ) ! Types (Real|Int)
        CHARACTER(LEN=IOULEN3), ALLOCATABLE :: EOUNITS( :,: ) ! Units  
        CHARACTER(LEN=IODLEN3), ALLOCATABLE :: EODESCS( :,: ) ! Dscriptions  

        CHARACTER(LEN=IOVLEN3), ALLOCATABLE :: AONAMES( :,: ) ! Names 
        INTEGER               , ALLOCATABLE :: AOTYPES( :,: ) ! Types (Real|Int)
        CHARACTER(LEN=IOULEN3), ALLOCATABLE :: AOUNITS( :,: ) ! Units  
        CHARACTER(LEN=IODLEN3), ALLOCATABLE :: AODESCS( :,: ) ! Dscriptions  

C...........   Other local allocatable arrays
        CHARACTER(LEN=IOVLEN3), ALLOCATABLE :: SAVEANAM( : )  ! tmp variables

C...........   Other local variables

        INTEGER       I, J, K, L, L2, V, V1, V2     ! counter and indices

        INTEGER    :: BASYR_OVR = 0     ! base year override
        INTEGER    :: FYEAR = 0         ! future year
        INTEGER       IOS       ! i/o status
        INTEGER       LEQU      ! position of '=' in formula
        INTEGER       LDIV      ! position of '-' or '+' in formula
        INTEGER       LMNS      ! position of '-' in formula
        INTEGER       LPLS      ! position of '+' in formula
        INTEGER       NIOVARS   ! Number of I/O API file non-emis variables
        INTEGER       NDATMAX   ! Max no of pols+activitys, based on I/O API
        INTEGER       NNPVAR    ! No. non-pollutant inventory variables
        INTEGER       YEAR      ! predominant year of inventory

        LOGICAL    :: CHKPLUS  = .FALSE. ! true: formula uses a + sign
        LOGICAL    :: CHKMINUS = .FALSE. ! true: formula uses a - sign
        LOGICAL    :: EFLAG    = .FALSE. ! true: error found

        CHARACTER*60  VAR_FORMULA
        CHARACTER*300 MESG      ! message buffer 

        CHARACTER(LEN=NAMLEN3)  NAMBUF ! file name buffer
        CHARACTER(LEN=IOVLEN3)  VIN_A
        CHARACTER(LEN=IOVLEN3)  VIN_B
        CHARACTER(LEN=IOVLEN3)  VNAME
        CHARACTER(LEN=IOULEN3)  UNITS  ! tmp units name

        CHARACTER*16 :: PROGNAME = 'OPENINVOUT' ! program name

C***********************************************************************
C   begin body of subroutine OPENINVOUT

C.........  Get environment variable settings
        CALL ENVSTR( FORMEVNM, MESG, ' ', VAR_FORMULA, IOS )

        BASYR_OVR = ENVINT( 'SMK_BASEYR_OVERRIDE', 
     &                      'Base year override', 0, IOS )

C.........  Get output inventory file names given source category
        CALL GETINAME( CATEGORY, ENAME, ANAME )

C.........  Depending on source category, set number of non-pollutant 
C           inventory variables
        SELECT CASE( CATEGORY )
        CASE( 'AREA' )
            IF( A2PFLAG ) THEN
                NNPVAR = NARVAR3 + 2
            ELSE
                NNPVAR = NARVAR3
            END IF
        CASE( 'MOBILE' )
            NNPVAR = NMBVAR3
        CASE( 'POINT' )
            NNPVAR = NPTVAR3
        END SELECT

C.........  Compute actual request output variables based on input data
        NIOVARS = NNPVAR + NIPOL * NPPOL + NIACT * NPACT

C.........  Compute conservative maximum number of pollutants and activities
C           (conservative b/c NPPOL > NPACT )
ccs        NDATMAX = INT( ( MXVARS3 - NNPVAR ) / NPPOL )

C.........  If there are too many output variables, reset NIPOL

ccs        IF( NIOVARS .GT. MXVARS3 ) THEN
ccs
ccs            WRITE( MESG,94010 ) 
ccs     &             'WARNING: Maximum number of pollutants or ' //
ccs     &             'activities that can be be'// CRLF()// BLANK10//
ccs     &             'written to the I/O API file is', NDATMAX, 
ccs     &             '. This limitation is caused by'// CRLF()// BLANK10//
ccs     &             'the I/O API variable limit of', MXVARS3, '.'
ccs            CALL M3MSG2( MESG )
ccs 
ccs            WRITE( MESG,94010 ) 
ccs     &             'WARNING: Reseting total number of output '//
ccs     &             'pollutants and activities to', NDATMAX
ccs            CALL M3MSG2( MESG )

C.............  If the number of pollutants alone are too much, then reset
C               that number.
ccs            IF( NIPOL .GT. NDATMAX ) THEN
ccs                NIPOL   = NDATMAX
ccs                NIACT   = 0
ccs                NIOVARS = NNPVAR + NPPOL * NIPOL                

C.............  Otherwise, reset the number of activities by subtracting the
C               number of pollutants from the maximum allowed number of 
C               variables
ccs            ELSE
ccs                NIACT = NDATMAX - NIPOL
ccs                NIOVARS = NNPVAR + NIPOL * NPPOL + NIACT * NPACT
ccs
ccs            END IF
ccs
ccs        ENDIF

C.........  Determine the predominant year of the inventory
        CALL GETBASYR( NSRC, YEAR )

C.........  Set the base year and past/future year, if needed
        IF( BASYR_OVR .GT. 0 .AND. YEAR .NE. BASYR_OVR ) THEN
            BYEAR = BASYR_OVR
            FYEAR = YEAR
        ELSE
            BYEAR = YEAR
        END IF

C.........  Set up for opening I/O API output file header

        CALL HDRMISS3  ! Initialize for emissions 

C.........  Set number of variables and allocate file description arrays
        NVARSET = NIOVARS + NPPOL  ! add extra space in case of formula based output
        
        ALLOCATE( VTYPESET( NVARSET ), STAT=IOS )
        CALL CHECKMEM( IOS, 'VTYPESET', PROGNAME )
        ALLOCATE( VNAMESET( NVARSET ), STAT=IOS )
        CALL CHECKMEM( IOS, 'VNAMESET', PROGNAME )
        ALLOCATE( VUNITSET( NVARSET ), STAT=IOS )
        CALL CHECKMEM( IOS, 'VUNITSET', PROGNAME )
        ALLOCATE( VDESCSET( NVARSET ), STAT=IOS )
        CALL CHECKMEM( IOS, 'VDESCSET', PROGNAME )
        
        VTYPESET = 0    ! array initialization
        VNAMESET = ' '  ! array initialization
        VUNITSET = ' '  ! array initialization
        VDESCSET = ' '  ! array initialization        
        
        NROWS3D = NSRC   !  number of rows = # of sources.

        FDESC3D( 1 ) = CATDESC // ' source inventory'
        FDESC3D( 2 ) = '/FROM/ ' // PROGNAME
        FDESC3D( 3 ) = '/VERSION/ ' // VERCHAR( CVSW )
        WRITE( FDESC3D( 4 ),94010 ) '/NON POLLUTANT/ ', NNPVAR

        IF( NIPOL .GT. 0 ) THEN
            WRITE( FDESC3D( 5 ),94010 ) '/POLLUTANTS/', NIPOL
            WRITE( FDESC3D( 6 ),94010 ) '/PER POLLUTANT/ ', NPPOL
        END IF

        IF( NIACT .GT. 0 ) THEN
            WRITE( FDESC3D( 7 ),94010 ) '/ACTIVITIES/', NIACT
            WRITE( FDESC3D( 8 ),94010 ) '/PER ACTIVITY/ ', NPACT
        END IF

        WRITE( FDESC3D( 9  ),94010 ) '/NUMBER CHARS/ ' , NCHARS 
        WRITE( FDESC3D( 10 ),94010 ) '/SCC POSITION/ ' , JSCC 
        WRITE( FDESC3D( 11 ),94010 ) '/STACK POSITION/ ' , JSTACK 
        WRITE( FDESC3D( 12 ),94010 ) '/BASE YEAR/ '    , BYEAR

        IF( GRDNM .NE. ' ' ) FDESC3D( 13 ) = '/GRIDNAME/ ' // GRDNM

        IF( FYEAR .GT. 0 ) THEN
             WRITE( FDESC3D( 14 ),94010 ) '/PROJECTED YEAR/ ', FYEAR
        END IF

C.........  Define source characteristic variables that are not strings

        J = 1
        VNAMESET( J ) = 'IFIP'
        VTYPESET( J ) = M3INT
        VUNITSET( J ) = 'n/a'
        VDESCSET( J ) = 'State and county FIPS code'
        J = J + 1

        VNAMESET( J ) = 'TZONES'
        VTYPESET( J ) = M3INT
        VUNITSET( J ) = 'n/a'
        VDESCSET( J ) = 'Time zone for site'
        J = J + 1

        VNAMESET( J ) = 'TPFLAG'
        VTYPESET( J ) = M3INT
        VUNITSET( J ) = 'T|2? T|3?'
        VDESCSET( J ) = 'Use week(2), month(3) temporal profiles or not'
        J = J + 1

        VNAMESET( J ) = 'INVYR'
        VTYPESET( J ) = M3INT
        VUNITSET( J ) = 'year AD'
        VDESCSET( J ) = 'Year of inventory for this record'
        J = J + 1

        SELECT CASE( CATEGORY )

        CASE( 'AREA' )
        
            IF( A2PFLAG ) THEN
                VNAMESET( J ) = 'XLOCA'
                VTYPESET( J ) = M3REAL
                VUNITSET( J ) = 'degrees'
                VDESCSET( J ) = 'longitude'
                J = J + 1
                
                VNAMESET( J ) = 'YLOCA'
                VTYPESET( J ) = M3REAL
                VUNITSET( J ) = 'degrees'
                VDESCSET( J ) = 'latitude'
                J = J + 1
            END IF

            VNAMESET( J ) = 'CELLID'
            VTYPESET( J ) = M3INT
            VUNITSET( J ) = 'n/a'
            VDESCSET( J ) = 'Cell number'
            J = J + 1

        CASE( 'MOBILE' )

            VNAMESET( J ) = 'IRCLAS'
            VTYPESET( J ) = M3INT
            VUNITSET( J ) = 'n/a'
            VDESCSET( J ) = 'Roadway type'
            J = J + 1

            VNAMESET( J ) = 'IVTYPE'
            VTYPESET( J ) = M3INT
            VUNITSET( J ) = 'n/a'
            VDESCSET( J ) = 'Vehicle type code'
            J = J + 1

            VNAMESET( J ) = 'XLOC1'
            VTYPESET( J ) = M3REAL
            VUNITSET( J ) = 'degrees'
            VDESCSET( J ) = 'Longitude at beginning of link'
            J = J + 1

            VNAMESET( J ) = 'YLOC1'
            VTYPESET( J ) = M3REAL
            VUNITSET( J ) = 'degrees'
            VDESCSET( J ) = 'Latitude at beginning of link'
            J = J + 1

            VNAMESET( J ) = 'XLOC2'
            VTYPESET( J ) = M3REAL
            VUNITSET( J ) = 'degrees'
            VDESCSET( J ) = 'Longitude at end of link'
            J = J + 1

            VNAMESET( J ) = 'YLOC2'
            VTYPESET( J ) = M3REAL
            VUNITSET( J ) = 'degrees'
            VDESCSET( J ) = 'Latitude at end of link'
            J = J + 1

        CASE( 'POINT' )

            VNAMESET( J ) = 'ISIC'
            VTYPESET( J ) = M3INT
            VUNITSET( J ) = 'n/a'
            VDESCSET( J ) = 'Source Industrial Code'
            J = J + 1

            VNAMESET( J ) = 'XLOCA'
            VTYPESET( J ) = M3REAL
            VUNITSET( J ) = 'degrees'
            VDESCSET( J ) = 'longitude'
            J = J + 1

            VNAMESET( J ) = 'YLOCA'
            VTYPESET( J ) = M3REAL
            VUNITSET( J ) = 'degrees'
            VDESCSET( J ) = 'latitude'
            J = J + 1

            VNAMESET( J ) = 'STKHT'
            VTYPESET( J ) = M3REAL
            VUNITSET( J ) = 'm'
            VDESCSET( J ) = 'Stack height'
            J = J + 1

            VNAMESET( J ) = 'STKDM'
            VTYPESET( J ) = M3REAL
            VUNITSET( J ) = 'm'
            VDESCSET( J ) = 'Stack diameter'
            J = J + 1

            VNAMESET( J ) = 'STKTK'
            VTYPESET( J ) = M3REAL
            VUNITSET( J ) = 'deg K'
            VDESCSET( J ) = 'Stack exhaust temperature'
            J = J + 1

            VNAMESET( J ) = 'STKVE'
            VTYPESET( J ) = M3REAL
            VUNITSET( J ) = 'm/s'
            VDESCSET( J ) = 'Stack exhaust velocity'
            J = J + 1

        END SELECT

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

C.........  Create output variables for pollutants
        DO V = 1 , NIPOL
            
            DO I = 1, NPPOL ! Loop through number of variables per pollutant

C.................  Set units for the primary data value
                IF( I .EQ. 1 ) THEN
                    K = INDEX1( EONAMES( V, 1 ), MXIDAT, INVDNAM )
                    UNITS = INVDUNT( K )

C.................  Set units for the other data values (per activity)
                ELSE
                    UNITS = EOUNITS( V, I )

                END IF

C.................  Store variable names and information
                VNAMESET( J ) = EONAMES( V, I )
                VTYPESET( J ) = EOTYPES( V, I )
                VUNITSET( J ) = UNITS
                VDESCSET( J ) = EODESCS( V, I )
                J = J + 1

            END DO    !  end loop on number of variables per pollutant

        END DO        !  end loop on inventory pollutants V

C.........  If there is a computed output variable, add it to the variable list.
C           Base the other variable settings (besides the name) on the first 
C           of the two variables in the formula.
        IF( VAR_FORMULA .NE. ' ' ) THEN

C.............  Make sure formula makes sense
            L2   = LEN_TRIM( VAR_FORMULA ) 
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
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF

C.............  Extract formula variable names
            VNAME = ADJUSTL ( VAR_FORMULA(      1:LEQU-1 ) )
            VIN_A = ADJUSTL ( VAR_FORMULA( LEQU+1:LDIV-1 ) )
            VIN_B = ADJUSTL ( VAR_FORMULA( LDIV+1:L2     ) )

C.............  Find formula inputs in existing variable list
            V1 = INDEX1( VIN_A, NIPOL, EINAM )
            V2 = INDEX1( VIN_B, NIPOL, EINAM )

            IF( V1 .LE. 0 ) THEN
                EFLAG = .TRUE.
                L = LEN_TRIM( VIN_A )
                MESG = 'ERROR: Variable "'// VIN_A( 1:L ) // '" from '//
     &                 'formula was not found in emissions inputs.'
                CALL M3MSG2( MESG )
            END IF

            IF( V2 .LE. 0 ) THEN
                EFLAG = .TRUE.
                L = LEN_TRIM( VIN_B )
                MESG = 'ERROR: Variable "'// VIN_B( 1:L ) // '" from '//
     &                 'formula was not found in emissions inputs.'
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

            IF ( EFLAG ) THEN
                MESG = 'ERROR: Problem processing formulas.'
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF

            V = INDEX1( VIN_A, NIPPA, EANAM )
            DO I = 1, NPPOL ! Loop through number of variables per pollutant

                L = LEN_TRIM( VIN_A )
                L = INDEX( EONAMES( V, I ), VIN_A(1:L) )
                IF ( L .GT. 1 ) THEN
        	    VNAMESET( J ) = EONAMES( V, I )( 1:L-1 ) // VNAME
                ELSE
                    VNAMESET( J ) = VNAME
                END IF
        	VTYPESET( J ) = EOTYPES( V, I )
        	VUNITSET( J ) = EOUNITS( V, I )
        	VDESCSET( J ) = EODESCS( V, I )
        	J = J + 1

            END DO    !  end loop on number of variables per pollutant

C.............  Update header settings
            WRITE( FDESC3D( 5 ),94010 ) '/POLLUTANTS/', NIPOL + 1

C.............  Update EANAM, in case needed for day- and hour-specific data
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

        ELSE

C.............  Decrease total number of variables since there isn't a formula based pollutant
            NVARSET = NVARSET - NPPOL

        END IF

C.........  Allocate names for activities
        ALLOCATE( AONAMES( NIACT,NPACT ), STAT=IOS )
        CALL CHECKMEM( IOS, 'AONAMES', PROGNAME )
        ALLOCATE( AOUNITS( NIACT,NPACT ), STAT=IOS )
        CALL CHECKMEM( IOS, 'AOUNITS', PROGNAME )
        ALLOCATE( AOTYPES( NIACT,NPACT ), STAT=IOS )
        CALL CHECKMEM( IOS, 'AOTYPES', PROGNAME )
        ALLOCATE( AODESCS( NIACT,NPACT ), STAT=IOS )
        CALL CHECKMEM( IOS, 'AODESCS', PROGNAME )

C.........  Get names, units, etc. of output activity-specific records
        CALL BLDENAMS( CATEGORY, NIACT, NPACT, ACTVTY, 
     &                 AONAMES, AOUNITS, AOTYPES, AODESCS )

C.........  Create output variables for activities
        DO V = 1 , NIACT
            
            DO I = 1, NPACT ! Loop through number of variables per activity

C.................  Set units for the primary data value
                IF( I .EQ. 1 ) THEN
                    K = INDEX1( AONAMES( V, 1 ), MXIDAT, INVDNAM )
                    UNITS = INVDUNT( K )

C.................  Set units for the other data values (per activity)
                ELSE
                    UNITS = AOUNITS( V, I )

                END IF

C.................  Store variable names and information
                VNAMESET( J ) = AONAMES( V, I )
                VTYPESET( J ) = AOTYPES( V, I )
                VUNITSET( J ) = UNITS
                VDESCSET( J ) = AODESCS( V, I )
                J = J + 1

            END DO    !  end loop on number of variables per activity

        END DO        !  end loop on inventory activities V

C.........  Prompt for and open I/O API output file
        NAMBUF= PROMPTSET( 
     &       'Enter logical name for the I/O API INVENTORY output file',
     &       FSUNKN3, ENAME, PROGNAME )
        ENAME = NAMBUF
        
C.........  Prompt for and open ASCII output file
        SDEV= PROMPTFFILE( 
     &      'Enter logical name for the ASCII INVENTORY output file',
     &      .FALSE., .TRUE., ANAME, PROGNAME )

C.........  Deallocate local memory
        DEALLOCATE( EONAMES, EOTYPES, EOUNITS, EODESCS )
        DEALLOCATE( AONAMES, AOTYPES, AOUNITS, AODESCS )

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )
 
        END SUBROUTINE OPENINVOUT

