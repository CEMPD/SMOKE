
        SUBROUTINE OPENSEF( SRCCT, DESC, SDATE, STIME, EMISDIR, FNAME )

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
        CHARACTER*16           PROMPTMFILE
        CHARACTER*16           VERCHAR
        LOGICAL                OPNFULL3

        EXTERNAL        GETCFDSC, PROMPTMFILE, VERCHAR, OPNFULL3

C...........   SUBROUTINE ARGUMENTS
        INTEGER     , INTENT    (IN) :: SRCCT    ! total number of sources
        CHARACTER(*), INTENT    (IN) :: DESC     ! description of file type
        INTEGER     , INTENT    (IN) :: SDATE    ! julian start date
        INTEGER     , INTENT    (IN) :: STIME    ! HHMMSS start time
        CHARACTER(*), INTENT    (IN) :: EMISDIR  ! directory of output files
        CHARACTER(*), INTENT(IN OUT) :: FNAME    ! name output emission factors file

C...........   LOCAL PARAMETERS
        CHARACTER*50, PARAMETER :: CVSW = '$Name$' ! CVS release tag

C...........   Other local variables
        INTEGER     J

        CHARACTER*300   MESG    !  message buffer
        CHARACTER*256   FULLNAME ! full file name

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

C.........  Set header values that cannot be default

        SDATE3D = SDATE
        STIME3D = 0
        TSTEP3D = 10000
        NVARS3D = NEFS + 1
        NROWS3D = SRCCT
        NLAYS3D = 1
 
        J = 1
        VNAME3D( J ) = 'SOURCES'
        UNITS3D( J ) = 'n/a'
        VDESC3D( J ) = 'Source number'
        VTYPE3D( J ) = M3INT
        
        DO J = 2, NVARS3D
            VNAME3D( J ) = EFSNAM( J - 1 )
            UNITS3D( J ) = EFSUNT( J - 1 )
            VDESC3D( J ) = EFSDSC( J - 1 )
            VTYPE3D( J ) = M3REAL
        END DO

C.........  Create full file name
        WRITE( FULLNAME,94010 ) EMISDIR( 1:LEN_TRIM( EMISDIR ) ) //
     &                          '/emisfacs.' // DESC // '.', SDATE,
     &                          '.ncf' 

C.........  Open new file
        IF( .NOT. OPNFULL3( FNAME, FSNEW3, FULLNAME, PROGNAME ) ) THEN
            MESG = 'Could not create new output file ' // FULLNAME
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( A, I7, A )

94030   FORMAT( A, F15.9, 1X, A )

        END SUBROUTINE OPENSEF


