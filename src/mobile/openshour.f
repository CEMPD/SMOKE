
        SUBROUTINE OPENSHOUR( ENAME, DESC, SDATE, EDATE, TVARNAME, 
     &                        NCOUNTY, TEMPDIR, FNAME )

C***********************************************************************
C  subroutine body starts at line 81
C
C  DESCRIPTION:
C      This subroutine opens the derived meteorology files needed for 
C      processing activity data with emission factors.
C 
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C      Copied from opensmet 4/01 by C. Seppanen
C
C************************************************************************
C
C Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
C                System
C File: @(#)$Id$
C
C COPYRIGHT (C) 1999, MCNC--North Carolina Supercomputing Center
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
C...........   This module contains the information about the source category
        USE MODINFO

C...........   This module is the derived meteorology data for emission factors
        USE MODMET

        IMPLICIT NONE

C...........   INCLUDES:
        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
        INCLUDE 'PARMS3.EXT'     !  I/O API parameters
        INCLUDE 'IODECL3.EXT'    !  I/O API function declarations
        INCLUDE 'FDESC3.EXT'     !  I/O API file description data structures

C...........   EXTERNAL FUNCTIONS and their descriptions:
        CHARACTER(LEN=IODLEN3) GETCFDSC
        CHARACTER*16           VERCHAR
        LOGICAL                OPNFULL3
        INTEGER                WKDAY

        EXTERNAL        GETCFDSC, VERCHAR, OPNFULL3, WKDAY

C...........   SUBROUTINE ARGUMENTS
        CHARACTER(*), INTENT    (IN) :: ENAME    ! name of inventory file
        CHARACTER(*), INTENT    (IN) :: DESC     ! description of file type
        INTEGER     , INTENT    (IN) :: SDATE    ! julian start date
        INTEGER     , INTENT    (IN) :: EDATE    ! julian end date
        CHARACTER(*), INTENT    (IN) :: TVARNAME ! name of variable for tmpr
        INTEGER     , INTENT    (IN) :: NCOUNTY  ! no. of counties in file
        CHARACTER(*), INTENT    (IN) :: TEMPDIR  ! directory of output files
        CHARACTER(*), INTENT(IN OUT) :: FNAME    ! name output hourly tmpr file

C...........   LOCAL PARAMETERS
        CHARACTER*50, PARAMETER :: CVSW = '$Name$' ! CVS release tag

C...........   Other local variables
        INTEGER     J
        INTEGER     DUMMYTIME    ! dummy time variable
        INTEGER     DUMMYDAY     ! dummy day variable
        INTEGER     FILE_EDATE   ! ending date of file
        INTEGER     CURRMNTH     ! current month
        INTEGER     PREVMNTH     ! previous month

        LOGICAL         FEXIST   ! true: file exists

        CHARACTER*300   MESG     ! message buffer
        CHARACTER*256   FULLNAME ! full file name

        CHARACTER(LEN=IODLEN3)  IFDESC2, IFDESC3 ! fields 2 & 3 from INVEN FDESC

        CHARACTER*16 :: PROGNAME = 'OPENSHOUR' ! program name

C***********************************************************************
C   begin body of subroutine OPENSHOUR

C.........  Determine end date based on file type and episode end
        FILE_EDATE = SDATE
        DUMMYTIME = 0

        IF( DESC == 'monthly' ) THEN
            CALL DAYMON( FILE_EDATE, PREVMNTH, DUMMYDAY )
        END IF

        DO
        	
C.............  For daily files, the end date is the start date
            IF( DESC == 'daily' ) EXIT

C.............  If the proposed end date is later than the episode end date
C               or we're using episode length files, the end date is the 
C               episode end date
            IF( FILE_EDATE >= EDATE .OR. DESC == 'episode' ) THEN
                FILE_EDATE = EDATE
                EXIT
            END IF

C.............  For weekly files, the end date is the next Saturday        
            IF( DESC == 'weekly' ) THEN
                IF( WKDAY( FILE_EDATE ) == 6 ) EXIT
            END IF

C.............  For monthly files, the end date is the last day of the month            
            IF( DESC == 'monthly' ) THEN
            	CALL DAYMON( FILE_EDATE, CURRMNTH, DUMMYDAY )
            	IF( CURRMNTH /= PREVMNTH ) THEN
            	    CALL NEXTIME( FILE_EDATE, DUMMYTIME, 24*10000 )
            	    EXIT
            	END IF
            END IF

            CALL NEXTIME( FILE_EDATE, DUMMYTIME, 24*10000 )
        END DO

C.........  Get header from inventory file
        IF ( .NOT. DESC3( ENAME ) ) THEN
            MESG = 'Could not get description of file "' 
     &             // ENAME( 1:LEN_TRIM( ENAME ) ) // '".'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

        IFDESC2 = GETCFDSC( FDESC3D, '/FROM/', .TRUE. )
        IFDESC3 = GETCFDSC( FDESC3D, '/VERSION/', .TRUE. )

C.........  Initialize I/O API output file headers
        CALL HDRMISS3

        FDESC3D( 1 ) = CATEGORY( 1:LEN_TRIM( CATEGORY ) ) // DESC //
     &                 ' temperature profiles file'
        FDESC3D( 2 ) = '/FROM/ '    // PROGNAME
        FDESC3D( 3 ) = '/VERSION/ ' // VERCHAR( CVSW )
        WRITE( FDESC3D( 4 ), 94030 ) '/MINTEMP/', MINTEMP
        WRITE( FDESC3D( 5 ), 94030 ) '/MAXTEMP/', MAXTEMP
        WRITE( FDESC3D( 6 ), 94030 ) '/T_UNITS/ "deg K"'
        FDESC3D( 7 ) = '/T_VNAME/ ' // TVARNAME
        FDESC3D( 8 ) = '/NOTE/ Time 000000 in file represents ' //
     &                 '6 AM in local time zone'
        WRITE( FDESC3D( 9 ), 94010 ) '/END DATE/ ', FILE_EDATE

        FDESC3D( 21 ) = '/INVEN FROM/ ' // IFDESC2
        FDESC3D( 22 ) = '/INVEN VERSION/ ' // IFDESC3

C.........  Set header values that cannot be default

        SDATE3D = SDATE
        STIME3D = 0
        TSTEP3D = 10000
        NVARS3D = 2
        NROWS3D = NCOUNTY
        NLAYS3D = 1
 
        J = 1
        VNAME3D( J ) = 'COUNTIES'
        UNITS3D( J ) = 'n/a'
        VDESC3D( J ) = 'County FIPS code'
        VTYPE3D( J ) = M3INT
        
        J = 2
        VNAME3D( J ) = 'TKCOUNTY'
        UNITS3D( J ) = 'K'
        VDESC3D( J ) = 'Hourly source temperature by county'
        VTYPE3D( J ) = M3REAL

C.........  Create full file name
        WRITE( FULLNAME,94010 ) TEMPDIR( 1:LEN_TRIM( TEMPDIR ) ) //
     &                          '/' // DESC // '.', SDATE, '.ncf'

C.........  Open new file
        IF( .NOT. OPNFULL3( FNAME, FSUNKN3, FULLNAME, PROGNAME ) ) THEN
            MESG = 'Could not create new output file ' // FULLNAME
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( A, I7, A )

94030   FORMAT( A, F15.9, 1X, A )

        END SUBROUTINE OPENSHOUR


