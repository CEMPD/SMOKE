
        SUBROUTINE OPENSHOUR( ENAME, SDATE, STIME, TVARNAME, FNAME )

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
        CHARACTER*16           PROMPTMFILE
        CHARACTER*16           VERCHAR

        EXTERNAL        GETCFDSC, PROMPTMFILE, VERCHAR

C...........   SUBROUTINE ARGUMENTS
        CHARACTER(*), INTENT    (IN) :: ENAME    ! name of inventory file
        INTEGER     , INTENT    (IN) :: SDATE    ! julian start date
        INTEGER     , INTENT    (IN) :: STIME    ! HHMMSS start time
        CHARACTER(*), INTENT    (IN) :: TVARNAME ! name of variable for tmpr
        CHARACTER(*), INTENT(IN OUT) :: FNAME    ! name output hourly tmpr file

C...........   LOCAL PARAMETERS
        CHARACTER*50, PARAMETER :: CVSW = '$Name$' ! CVS release tag

C...........   Other local variables
        INTEGER     J

        CHARACTER*300   MESG    !  message buffer

        CHARACTER(LEN=IODLEN3)  IFDESC2, IFDESC3 ! fields 2 & 3 from INVEN FDESC

        CHARACTER*16 :: PROGNAME = 'OPENSHOUR' ! program name

C***********************************************************************
C   begin body of subroutine OPENSHOUR

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

        FDESC3D( 1 ) = CATEGORY( 1:LEN_TRIM( CATEGORY ) ) //
     &                 ' hourly temperature file'
        FDESC3D( 2 ) = '/FROM/ '    // PROGNAME
        FDESC3D( 3 ) = '/VERSION/ ' // VERCHAR( CVSW )
        WRITE( FDESC3D( 4 ), 94030 ) '/MINT_MIN/', MINT_MIN
        WRITE( FDESC3D( 5 ), 94030 ) '/MINT_MAX/', MINT_MAX
        WRITE( FDESC3D( 6 ), 94030 ) '/MAXT_MIN/', MAXT_MIN
        WRITE( FDESC3D( 7 ), 94030 ) '/MAXT_MAX/', MAXT_MAX
        WRITE( FDESC3D( 9 ), 94030 ) '/T_MAXINTVL/', TMXINVL
        WRITE( FDESC3D( 10 ), 94030 ) '/T_UNITS/ "deg F"'
        WRITE( FDESC3D( 11 ), 94030 ) '/T_VNAME/ ' // TVARNAME

        FDESC3D( 21 ) = '/INVEN FROM/ ' // IFDESC2
        FDESC3D( 22 ) = '/INVEN VERSION/ ' // IFDESC3

C.........  Set header values that cannot be default

        SDATE3D = SDATE
        STIME3D = 0
        TSTEP3D = 10000
        NVARS3D = 1
        NROWS3D = NSRC
        NLAYS3D = 1
 
        J = 1
        VNAME3D( J ) = 'TKHOUR'
        UNITS3D( J ) = 'K'
        VDESC3D( J ) = 'Hourly source temperature'
        VTYPE3D( J ) = M3REAL

C.........  Prompt for Source-based Output file
        FNAME = PROMPTMFILE(
     &          'Enter logical name for HOURLY TEMPERATURE file',
     &          FSUNKN3, FNAME, PROGNAME )

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

94030   FORMAT( A, F15.9, 1X, A )

        END SUBROUTINE OPENSHOUR


