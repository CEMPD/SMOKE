
        SUBROUTINE OPENINVOUT( A2PFLAG, GRDNM, ENAME, ANAME, IDEV, 
     &                        SDEV, ADEV, VARPATH )

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
C.........  This module contains the information about the source category
        USE MODINFO, ONLY: CATEGORY, CATDESC, NSRC, BYEAR,
     &                     NIPOL, NPPOL, NIACT, NPACT, NCHARS,
     &                     JSCC, JSTACK, NCOMP, VAR_FORMULA

C.........  This module is required by the FileSetAPI
        USE MODFILESET
        
        IMPLICIT NONE

C...........   INCLUDES
        INCLUDE 'EMCNST3.EXT'   !  emissions constat parameters
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
        INCLUDE 'SETDECL.EXT'   !  FileSetAPI variables and functions

C...........   EXTERNAL FUNCTIONS and their descriptionsNRAWIN
        CHARACTER(2)    CRLF
        INTEGER         ENVINT
        INTEGER         GETEFILE
        INTEGER         INDEX1
        INTEGER         PROMPTFFILE
        LOGICAL         SETENVVAR
        CHARACTER(16)   VERCHAR

        EXTERNAL    CRLF, ENVINT, GETEFILE, INDEX1, PROMPTFFILE, 
     &              SETENVVAR, VERCHAR

C...........   SUBROUTINE ARGUMENTS
        LOGICAL     , INTENT(IN)  :: A2PFLAG  ! true: using area-to-point processing
        CHARACTER(*), INTENT(IN)  :: GRDNM    ! grid name if any gridded data
        CHARACTER(*), INTENT(OUT) :: ENAME    ! emis i/o api inven logical name
        CHARACTER(*), INTENT(OUT) :: ANAME    ! emis ASCII inven logical name
        INTEGER     , INTENT(OUT) :: IDEV     ! map inventory file
        INTEGER     , INTENT(OUT) :: SDEV     ! ascii output inven file unit no.
        INTEGER     , INTENT(OUT) :: ADEV     ! REPINVEN output file unit no.
        CHARACTER(PHYLEN3), INTENT( OUT ) :: VARPATH ! path for pol/act output files

C...........   LOCAL PARAMETERS
        CHARACTER(16), PARAMETER :: FORMEVNM = 'SMKINVEN_FORMULA'
        CHARACTER(50), PARAMETER :: 
     &  CVSW = '$Name SMOKEv4.9_Jun2022$' ! CVS release tag

C...........   Other local variables

        INTEGER       I, J, K, L, L2, N, V, V1, V2   ! counter and indices

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

        CHARACTER(256) MESG      ! message buffer 

        CHARACTER(80)       NAME1  ! tmp file name component
        CHARACTER(80)       NAME2  ! tmp file name component
        CHARACTER(NAMLEN3)  NAMBUF ! file name buffer
        CHARACTER(IOULEN3)  UNITS  ! tmp units name
        CHARACTER(PHYLEN3)  APHYS  ! ASCII physical file name
        CHARACTER(PHYLEN3)  EPHYS  ! I/O API physical file name
        CHARACTER(PHYLEN3)  PATH   ! path name

        CHARACTER(16) :: PROGNAME = 'OPENINVOUT' ! program name

C***********************************************************************
C   begin body of subroutine OPENINVOUT

C.........  Get environment variable settings
        CALL ENVSTR( FORMEVNM, MESG, ' ', VAR_FORMULA, IOS )

        BASYR_OVR = ENVINT( 'SMK_BASEYR_OVERRIDE', 
     &                      'Base year override', 0, IOS )

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

C.........  Get output inventory file names given source category
        CALL GETINAME( CATEGORY, ENAME, ANAME )

C.........  Open map-formatted inventory file without prompting
        IDEV = GETEFILE( ENAME, .FALSE., .TRUE., PROGNAME )
        IF ( IDEV .LT. 0 ) THEN     !  failure to open

            MESG = 'Could not open INVENTORY MAP file:' // CRLF() // 
     &              BLANK10 // TRIM( ENAME ) // '.'
            CALL M3MSG2( MESG )

            MESG = 'Ending program.'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

        END IF      !  if getefile() failed

C.........  Evaluate physical file name of inventory map
        MESG = 'Inventory map file name'
        CALL ENVSTR( ENAME, MESG, BLANK16, APHYS, IOS )

        IF( IOS .NE. 0 ) THEN
            MESG = 'Unable to evaluate environment variable "' //
     &             TRIM( ENAME ) // '"'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

C.........  Determine path of inventory map
        L = LEN_TRIM( APHYS )
        DO N = L, 1, -1

            IF( APHYS( N:N ) .EQ. '/' .OR.
     &          APHYS( N:N ) .EQ. '\'      ) THEN
                PATH = APHYS( 1:N )
                EXIT
            END IF

        END DO

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

C.........  Determine the predominant year of the inventory
        CALL GETBASYR( NSRC, YEAR )

C.........  Set the base year and past/future year, if needed
        IF( BASYR_OVR .GT. 0 .AND. YEAR .NE. BASYR_OVR ) THEN
            BYEAR = BASYR_OVR
            FYEAR = YEAR
        ELSE
            BYEAR = YEAR
        END IF

C.........  Set up for opening I/O API part of inventory output
        CALL HDRMISS3  ! Initialize for emissions 

C.........  Set number of variables and allocate file description arrays
        NVARSET = NNPVAR
        
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

C........  Get setup up for later processing of a formula for adding a variable
        L = LEN_TRIM( VAR_FORMULA )
        IF( L .GT. 0 ) THEN

C.............  Figure out how many variables there are based on the
C               number of commas found in the string.
            NCOMP = 1
            DO I = 1, L
                IF( VAR_FORMULA( I:I ) == ',' ) NCOMP = NCOMP + 1
            ENDDO

            WRITE( FDESC3D( 5 ),94010 ) '/POLLUTANTS/', NIPOL + NCOMP
        END IF  

C.........  Define source characteristic variables that are not strings

        J = 1
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
                VTYPESET( J ) = M3DBLE
                VUNITSET( J ) = 'degrees'
                VDESCSET( J ) = 'longitude'
                J = J + 1
                
                VNAMESET( J ) = 'YLOCA'
                VTYPESET( J ) = M3DBLE
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
            VTYPESET( J ) = M3DBLE
            VUNITSET( J ) = 'degrees'
            VDESCSET( J ) = 'Longitude at beginning of link'
            J = J + 1

            VNAMESET( J ) = 'YLOC1'
            VTYPESET( J ) = M3DBLE
            VUNITSET( J ) = 'degrees'
            VDESCSET( J ) = 'Latitude at beginning of link'
            J = J + 1

            VNAMESET( J ) = 'XLOC2'
            VTYPESET( J ) = M3DBLE
            VUNITSET( J ) = 'degrees'
            VDESCSET( J ) = 'Longitude at end of link'
            J = J + 1

            VNAMESET( J ) = 'YLOC2'
            VTYPESET( J ) = M3DBLE
            VUNITSET( J ) = 'degrees'
            VDESCSET( J ) = 'Latitude at end of link'
            J = J + 1

        CASE( 'POINT' )

            VNAMESET( J ) = 'XLOCA'
            VTYPESET( J ) = M3DBLE
            VUNITSET( J ) = 'degrees'
            VDESCSET( J ) = 'longitude'
            J = J + 1

            VNAMESET( J ) = 'YLOCA'
            VTYPESET( J ) = M3DBLE
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

            VNAMESET( J ) = 'FUG_HEIGHT'
            VTYPESET( J ) = M3REAL
            VUNITSET( J ) = 'm'
            VDESCSET( J ) = 'Fugitive emissions height'
            J = J + 1

            VNAMESET( J ) = 'FUG_WIDTH'
            VTYPESET( J ) = M3REAL
            VUNITSET( J ) = 'm'
            VDESCSET( J ) = 'Fugitive emissions width (Y DIM)'
            J = J + 1

            VNAMESET( J ) = 'FUG_LENGTH'
            VTYPESET( J ) = M3REAL
            VUNITSET( J ) = 'm'
            VDESCSET( J ) = 'Fugitive emissions length (X DIM)'
            J = J + 1

            VNAMESET( J ) = 'FUG_ANGLE'
            VTYPESET( J ) = M3REAL
            VUNITSET( J ) = 'm'
            VDESCSET( J ) = 'Fugitive emissions angle'
            J = J + 1

        END SELECT

C.........  Build output file physical name
        ENAME = 'IOAPI_INV'
        EPHYS = TRIM( PATH ) // TRIM( NAME1 ) // '.ncf'

C.........  Set output logical file name
        IF( .NOT. SETENVVAR( ENAME, EPHYS ) ) THEN
            EFLAG = .TRUE.
            MESG = 'ERROR: Could not set logical file name for file:' //
     &              CRLF() // BLANK10 // TRIM( EPHYS )
            CALL M3MSG2( MESG )

C........  Open I/O API file
        ELSE IF( .NOT. OPENSET( ENAME, FSNEW3, PROGNAME ) ) THEN
            EFLAG = .TRUE.
            MESG = 'ERROR: Could not open I/O API inventory file ' //
     &             'for file name:' // CRLF() // BLANK10 // 
     &             TRIM( EPHYS )
            CALL M3MSG2( MESG )

        END IF

C.........  Build ASCII file physical name
        APHYS = TRIM( PATH ) // TRIM( NAME2 ) // '.txt'
        
C.........  Set output logical file name for ASCII inventory
        IF( .NOT. SETENVVAR( ANAME, APHYS ) ) THEN
            EFLAG = .TRUE.
            MESG = 'ERROR: Could not set logical file name for file:' //
     &              CRLF() // BLANK10 // TRIM( APHYS )
            CALL M3MSG2( MESG )
        END IF

C.........  Open ASCII output file
        SDEV = GETEFILE( ANAME, .FALSE., .TRUE., PROGNAME )
        IF ( SDEV .LT. 0 ) THEN     !  failure to open

            MESG = 'Could not open output file:' // CRLF() // 
     &              BLANK10 // TRIM( ANAME ) // '.'
            CALL M3MSG2( MESG )

            MESG = 'Ending program.'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

        END IF      !  if getefile() failed

C.........  Prompt for and open REPINVEN file
        ADEV = PROMPTFFILE(
     &      'Enter logical name for the REPINVEN file',
     &      .FALSE., .TRUE., 'REPINVEN', PROGNAME )

C.........  Abort if error
        IF( EFLAG ) THEN
            MESG = 'Problem opening inventory output files. ' //
     &             'See listed error(s) above.'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 ) 
        END IF

C.........  Provide variable path
        VARPATH = TRIM( PATH ) // TRIM( NAME1 ) // '_dat'

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

        END SUBROUTINE OPENINVOUT

