
        SUBROUTINE OPENSEF( SRCCT, DESC, SDATE, EDATE, EMISDIR, FNAME )

C...........   MODULES for public variables
C...........   This module contains the information about the source category
        USE MODINFO

C.........  This module contains emission factor tables and related
        USE MODEMFAC
        
        IMPLICIT NONE

C...........   INCLUDES:
        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
        INCLUDE 'PARMS3.EXT'     !  I/O API parameters
        INCLUDE 'IODECL3.EXT'    !  I/O API function declarations
        INCLUDE 'FDESC3.EXT'     !  I/O API file description data structures

C...........   EXTERNAL FUNCTIONS and their descriptions:
        CHARACTER(LEN=IODLEN3) GETCFDSC
        INTEGER                INDEX1
        CHARACTER*16           PROMPTMFILE
        CHARACTER*16           VERCHAR
        LOGICAL                SETENVVAR

        EXTERNAL     GETCFDSC, INDEX1, PROMPTMFILE, VERCHAR, SETENVVAR

C...........   SUBROUTINE ARGUMENTS
        INTEGER     , INTENT    (IN) :: SRCCT    ! total number of sources
        CHARACTER(*), INTENT    (IN) :: DESC     ! description of file type
        INTEGER     , INTENT    (IN) :: SDATE    ! julian start date
        INTEGER     , INTENT    (IN) :: EDATE    ! julian end date
        CHARACTER(*), INTENT    (IN) :: EMISDIR  ! directory of output files
        CHARACTER(*), INTENT(IN OUT) :: FNAME    ! name output emission factors file

C...........   LOCAL PARAMETERS
        CHARACTER*50, PARAMETER :: CVSW = '$Name$' ! CVS release tag

C...........   Other local variables
        INTEGER     J, K

        CHARACTER*300   MESG      ! message buffer
        CHARACTER*256   FULLNAME  ! full file name
        CHARACTER*16    CURRVNAME ! current variable name

        CHARACTER*16 :: PROGNAME = 'OPENSEF' ! program name

C***********************************************************************
C   begin body of subroutine OPENSEF

C.........  Initialize I/O API output file headers
        CALL HDRMISS3

        FDESC3D( 1 ) = CATEGORY( 1:LEN_TRIM( CATEGORY ) ) //
     &                 ' emission factors file'
        FDESC3D( 2 ) = '/FROM/ '    // PROGNAME
        FDESC3D( 3 ) = '/VERSION/ ' // VERCHAR( CVSW )
        FDESC3D( 4 ) = '/NOTE/ Time 000000 in file represents ' //
     &                 '6 AM in local time zone'
        WRITE( FDESC3D( 5 ), 94010 ) '/END DATE/ ', EDATE

C.........  Set header values that cannot be default

        SDATE3D = SDATE
        STIME3D = 0
        TSTEP3D = 10000
        NVARS3D = MXETYPE + 1
        NROWS3D = SRCCT
        NLAYS3D = 1
 
        J = 1
        VNAME3D( J ) = 'SOURCES'
        UNITS3D( J ) = 'n/a'
        VDESC3D( J ) = 'Source number'
        VTYPE3D( J ) = M3INT
        
C.........  Loop through emission process/pollutant combos 
C           (EMTNAM is created from MEPROC file and contains only the ones we want, 
C           EFS* contains all possible MOBILE6 outputs with units and descriptions)        
        DO J = 2, NVARS3D
            CURRVNAME = TRIM( EMTNAM( J-1,1 ) )
            K = INDEX1( CURRVNAME, NEFS, EFSNAM )
            VNAME3D( J ) = EFSNAM( K )
            UNITS3D( J ) = EFSUNT( K )
            VDESC3D( J ) = EFSDSC( K )
            VTYPE3D( J ) = M3REAL
        END DO

C.........  Create full file name
        WRITE( FULLNAME,94010 ) EMISDIR( 1:LEN_TRIM( EMISDIR ) ) //
     &                          '/emisfacs.' // 
     &                          DESC( 1:LEN_TRIM( DESC ) ) // '.', 
     &                          SDATE, '.ncf'

C.........  Set logical file name
        IF( .NOT. SETENVVAR( FNAME, FULLNAME ) ) THEN
            MESG = 'Could not set logical file name for file ' //
     &             TRIM( FULLNAME )
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF       

C.........  Open new file
        IF( .NOT. OPEN3( FNAME, FSUNKN3, PROGNAME ) ) THEN
            MESG = 'Could not create new output file ' // 
     &             TRIM( FULLNAME )
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( A, I7, A )

94020   FORMAT( A, I7, A1, I7, A )

        END SUBROUTINE OPENSEF


