
        SUBROUTINE OPENSEF( SRCCT, SDATE, STIME, FNAME )

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

        EXTERNAL        GETCFDSC, PROMPTMFILE, VERCHAR

C...........   SUBROUTINE ARGUMENTS
        INTEGER     , INTENT    (IN) :: SRCCT    ! total number of sources
        INTEGER     , INTENT    (IN) :: SDATE    ! julian start date
        INTEGER     , INTENT    (IN) :: STIME    ! HHMMSS start time
        CHARACTER(*), INTENT(IN OUT) :: FNAME    ! name output emission factors file

C...........   LOCAL PARAMETERS
        CHARACTER*50, PARAMETER :: CVSW = '$Name$' ! CVS release tag

C...........   Other local variables
        INTEGER     J

        CHARACTER*300   MESG    !  message buffer

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

C.........  Prompt for source-based output file
        FNAME = PROMPTMFILE(
     &          'Enter logical name for EMISSION FACTORS file',
     &          FSUNKN3, FNAME, PROGNAME )

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

94030   FORMAT( A, F15.9, 1X, A )

        END SUBROUTINE OPENSEF


