
        SUBROUTINE OPENCMAT( ENAME, MATTYP, MNAME )

C***********************************************************************
C  subroutine body starts at line 88
C
C  DESCRIPTION:
C      This subroutine opens the file in which the control matrix
C      will be written.
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C     
C***************************************************************************
C
C Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
C                System
C File: @(#)$Id$
C
C COPYRIGHT (C) 2000, MCNC--North Carolina Supercomputing Center
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

C.........  MODULES for public variables
C.........  This module contains the control packet data and control matrices
        USE MODCNTRL

C.........  This module contains the information about the source category
        USE MODINFO

        IMPLICIT NONE

C...........   INCLUDES
        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
        INCLUDE 'FDESC3.EXT'    !  I/O API file description data structures.

C...........   EXTERNAL FUNCTIONS and their descriptions:
        CHARACTER*2            CRLF
        CHARACTER(LEN=IODLEN3) GETCFDSC
        INTEGER                GETIFDSC   
        CHARACTER(LEN=NAMLEN3) PROMPTMFILE
        CHARACTER*16           VERCHAR

        EXTERNAL     CRLF, GETCFDSC, GETIFDSC, PROMPTMFILE, VERCHAR

C...........   LOCAL PARAMETERS
        CHARACTER*50, PARAMETER :: CVSW = '$Name$' ! CVS release tag

C.........  SUBROUTINE ARGUMENTS
        CHARACTER(*), INTENT (IN) :: ENAME      ! emissions inven logical name
        CHARACTER(*), INTENT (IN) :: MATTYP     ! matrix type
        CHARACTER(*), INTENT(OUT) :: MNAME      ! controls file name

C.........  Other local variables
        INTEGER          J              !  counters and indices

        CHARACTER(LEN=NAMLEN3) NAMBUF   ! file name buffer
        CHARACTER*300          MESG     ! message buffer

        CHARACTER(LEN=IODLEN3) IFDESC2, IFDESC3 ! fields 2 & 3 from inven FDESC
        CHARACTER(LEN=IOVLEN3) UNITS    ! emissions units

        CHARACTER*16 :: PROGNAME = 'OPENCMAT' ! program name

C***********************************************************************
C   begin body of subroutine OPENCMAT

C.........  Get header information from inventory file

        IF ( .NOT. DESC3( ENAME ) ) THEN
            MESG = 'Could not get description of file "' 
     &             // ENAME( 1:LEN_TRIM( ENAME ) ) // '".'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

        IFDESC2 = GETCFDSC( FDESC3D, '/FROM/', .TRUE. )
        IFDESC3 = GETCFDSC( FDESC3D, '/VERSION/', .TRUE. )
        J       = GETIFDSC( FDESC3D, '/NON POLLUTANT/', .TRUE. )
        UNITS   = UNITS3D( J + 1 )

C.........  Initialize I/O API output file headers
        CALL HDRMISS3

C.........  Set I/O API header parms that need values
        NROWS3D = NSRC
        NVARS3D = NVCMULT

        FDESC3D( 1 ) = CATEGORY( 1:CATLEN ) // ' control matrix'
        FDESC3D( 2 ) = '/FROM/ '    // PROGNAME
        FDESC3D( 3 ) = '/VERSION/ ' // VERCHAR( CVSW )
        WRITE( FDESC3D( 4 ), '(A,I4)' ) '/CTYPE/ ', CTYPMULT

        FDESC3D( 11 ) = '/INVEN FROM/ ' // IFDESC2
        FDESC3D( 12 ) = '/INVEN VERSION/ ' // IFDESC3

C.........  Set up non-speciation variables
        DO J = 1,NVCMULT

           VNAME3D( J )= PNAMMULT( J )
           VTYPE3D( J )= M3REAL
           UNITS3D( J )= 'fraction'
           VDESC3D( J )= 'Multiplicative control factor for pollutant '
     &                   // PNAMMULT( J )

        END DO

C.........  Error if number of variables is passed maximum because we can't
C           store the names of the variables.
C.........  DO NOT end program here because it will be ended when the write
C           attempt is made for these extra variables.
        IF( J .GT. MXVARS3 ) THEN

            WRITE( MESG, 94010 ) 
     &             'Maximum I/O API variables exceeded:' //
     &             CRLF() // BLANK10 // 'Max: ', MXVARS3, 'Actual:', J
            CALL M3MSG2( MESG )

        ENDIF

        MESG = 'Enter logical name for control matrix...'
        CALL M3MSG2( MESG )

C.........  Open control matrix.
C.........  Using NAMBUF is needed for HP to ensure string length consistencies
        MESG = 'I/O API ' // MATTYP // ' CONTROL MATRIX'

        NAMBUF = PROMPTMFILE( MESG, FSUNKN3, CRL // 'CMAT', 
     &                        PROGNAME )
        MNAME = NAMBUF

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

        END SUBROUTINE OPENCMAT
