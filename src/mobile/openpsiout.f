
        SUBROUTINE OPENPSIOUT( ENAME, TNAME, FNAME, FDEV, VNAME )

C***********************************************************************
C  subroutine body starts at line 
C
C  DESCRIPTION:
C      This subroutine opens the intermediate files for processing 
C      parameter scheme indices (PSIs). One file contains the min/max
C      temperature combinations used for each PSI and activity, and 
C      one file contains the PSIs to be used for each source for
C      24 hours.
C 
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C      Created 6/99 by M. Houyoux
C
C****************************************************************************/
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
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
        INCLUDE 'FDESC3.EXT'    !  I/O API file description data structures

C...........   EXTERNAL FUNCTIONS and their descriptions:
        CHARACTER(LEN=IODLEN3) GETCFDSC
        INTEGER                PROMPTFFILE
        CHARACTER*16           PROMPTMFILE
        CHARACTER*16           VERCHAR

        EXTERNAL        GETCFDSC, PROMPTFFILE, PROMPTMFILE, VERCHAR

C...........   SUBROUTINE ARGUMENTS
        CHARACTER(*), INTENT    (IN) :: ENAME  ! name of inventory file
        CHARACTER(*), INTENT(IN OUT) :: TNAME  ! name of output tmpr/PSIs file
        CHARACTER(*), INTENT(IN OUT) :: FNAME  ! name of output src/PSIs file
        INTEGER     , INTENT   (OUT) :: FDEV   ! Unit no. of output tmpr/PSIs
        CHARACTER(*), INTENT   (OUT) :: VNAME( NIACT ) ! variable names

C...........   LOCAL PARAMETERS
        CHARACTER*50, PARAMETER :: SCCSW = ! SCCS string with version no. at end
     &               '@(#)$Id$'

C...........   Other local variables
        INTEGER     J, L        ! counters and indices

        CHARACTER*2     NUM     ! tmp variable number
        CHARACTER*300   MESG    ! message buffer

        CHARACTER(LEN=IODLEN3)  IFDESC2, IFDESC3 ! fields 2 & 3 from INVEN FDESC

        CHARACTER*16 :: PROGNAME = 'OPENPSIOUT' ! program name

C***********************************************************************
C   begin body of subroutine OPENPSIOUT

C.........  Get header from inventory file
        IF ( .NOT. DESC3( ENAME ) ) THEN
            MESG = 'Could not get description of file "' 
     &             // ENAME( 1:LEN_TRIM( ENAME ) ) // '".'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

        IFDESC2 = GETCFDSC( FDESC3D, '/FROM/', .TRUE. )
        IFDESC3 = GETCFDSC( FDESC3D, '/VERSION/', .TRUE. )

C.........  Initialize source/PSI I/O API output file headers
        CALL HDRMISS3

        FDESC3D( 1 ) = CATEGORY( 1:LEN_TRIM( CATEGORY ) ) //
     &                 ' PSI by source'
        FDESC3D( 2 ) = '/FROM/ '    // PROGNAME
        FDESC3D( 3 ) = '/VERSION/ ' // VERCHAR( SCCSW )

        FDESC3D( 21 ) = '/INVEN FROM/ ' // IFDESC2
        FDESC3D( 22 ) = '/INVEN VERSION/ ' // IFDESC3

C.........  Set header values that cannot be default

        TSTEP3D = 240000
        NVARS3D = 1
        NROWS3D = NSRC  
        NCOLS3D = 24
 
        DO J = 1, NIACT

            L = LEN_TRIM( ACTVTY( J ) )

            WRITE( NUM, '(I2.2)' ) J
            VNAME  ( J ) = 'SRCPSI' // NUM
            VNAME3D( J ) = VNAME( J )
            UNITS3D( J ) = 'n/a' 
            VDESC3D( J ) = ACTVTY( J )( 1:L ) // ' PSI' // 
     &                     ', by source & 24 hours'          
            VTYPE3D( J ) = M3INT

        END DO

C.........  Prompt for source-based output file of PSIs
        FNAME = PROMPTMFILE(
     &          'Enter logical name for SOURCE PSIs file',
     &          FSUNKN3, FNAME, PROGNAME )

C.........  Prompt for file with min-max temperatures per PSI
        FDEV  = PROMPTFFILE(
     &          'Enter logical name for PSI/MIN-MAX-TEMPERATURE file',
     &          .FALSE., .TRUE., TNAME, PROGNAME )

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

94030   FORMAT( A, F10.3, 1X, A )

        END SUBROUTINE OPENPSIOUT


