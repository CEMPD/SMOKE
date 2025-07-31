
        SUBROUTINE OPENPMAT( ENAME, BYEARIN, PYEAR, PNAME )

C***********************************************************************
C  subroutine body starts at line 87
C
C  DESCRIPTION:
C      Open the projection matrix.
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C     
C
C***********************************************************************
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
        USE MODINFO, ONLY: CATEGORY, CATLEN, CRL, NSRC

C.........  This module contains the control packet data and control matrices
        USE MODCNTRL, ONLY: POLSFLAG, NVPROJ, PNAMPROJ

C.........This module is required by the FileSetAPI
        USE MODFILESET

        IMPLICIT NONE

C...........   INCLUDES
        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
        INCLUDE 'SETDECL.EXT'   !  FileSetAPI variables and functions

C...........   EXTERNAL FUNCTIONS and their descriptions:
        CHARACTER(2)       CRLF
        CHARACTER(IODLEN3) GETCFDSC
        INTEGER            GETIFDSC   
        CHARACTER(16)      VERCHAR

        EXTERNAL     CRLF, GETCFDSC, GETIFDSC, VERCHAR

C...........   LOCAL PARAMETERS
        CHARACTER(50), PARAMETER :: 
     &  CVSW = '$Name SMOKEv5.2_Jul2025$' ! CVS release tag

C.........  SUBROUTINE ARGUMENTS
        CHARACTER(*), INTENT (IN) :: ENAME      ! emissions inven logical name
        INTEGER     , INTENT (IN) :: BYEARIN    ! base year of proj factors
        INTEGER     , INTENT (IN) :: PYEAR      ! projected year of proj factors
        CHARACTER(*), INTENT(OUT) :: PNAME      ! projection file name

      
C.........  Other local variables
        INTEGER          I, J           !  counters and indices
        INTEGER          IOS            !  i/o status

        CHARACTER(NAMLEN3) NAMBUF   ! file name buffer
        CHARACTER(300)     MESG     ! message buffer

        CHARACTER(IODLEN3) IFDESC2, IFDESC3 ! fields 2 & 3 from inven FDESC

        CHARACTER(16) :: PROGNAME = 'OPENPMAT' ! program name

C***********************************************************************
C   begin body of subroutine OPENPMAT

C.........  Get header information from inventory file

        IF ( .NOT. DESCSET( ENAME,-1 ) ) THEN
            MESG = 'Could not get description of file "' 
     &             // ENAME( 1:LEN_TRIM( ENAME ) ) // '".'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

        IFDESC2 = GETCFDSC( FDESC3D, '/FROM/', .TRUE. )
        IFDESC3 = GETCFDSC( FDESC3D, '/VERSION/', .TRUE. )

C.........  Initialize I/O API output file headers
        CALL HDRMISS3

C.........  Set I/O API header parms that need values
        NROWS3D = NSRC

        FDESC3D( 1 ) = CATEGORY( 1:CATLEN ) // ' projection matrix'
        FDESC3D( 2 ) = '/FROM/ '    // PROGNAME
        FDESC3D( 3 ) = '/VERSION/ ' // VERCHAR( CVSW )

        WRITE( FDESC3D( 4 ), '(A,I4)' ) '/CTYPE/ ', CTYPPROJ
        WRITE( FDESC3D( 5 ), '(A,I4)' ) '/BASE YEAR/ ', BYEARIN
        WRITE( FDESC3D( 6 ), '(A,I4)' ) '/PROJECTED YEAR/ ', PYEAR

        FDESC3D( 11 ) = '/INVEN FROM/ ' // IFDESC2
        FDESC3D( 12 ) = '/INVEN VERSION/ ' // IFDESC3

        IF( ALLOCATED( VTYPESET ) ) 
     &      DEALLOCATE( VTYPESET, VNAMESET, VUNITSET, VDESCSET )

C.........  Set number of variables in output file based on whether pollutant-specific
C           assignments are being used
        IF( POLSFLAG ) THEN
            NVARSET = NVPROJ
        ELSE
            NVARSET = 1
        END IF

        ALLOCATE( VTYPESET( NVARSET ), STAT=IOS )
        CALL CHECKMEM( IOS, 'VTYPESET', PROGNAME )
        ALLOCATE( VNAMESET( NVARSET ), STAT=IOS )
        CALL CHECKMEM( IOS, 'VNAMESET', PROGNAME )
        ALLOCATE( VUNITSET( NVARSET ), STAT=IOS )
        CALL CHECKMEM( IOS, 'VUNITSET', PROGNAME )
        ALLOCATE( VDESCSET( NVARSET ), STAT=IOS )
        CALL CHECKMEM( IOS, 'VDESCSET', PROGNAME )

C.........  Also deallocate the number of variables per file so
C           that this will be set automatically by openset
        DEALLOCATE( VARS_PER_FILE )

C.........  If pollutant-specific assignments, then set up projection matrix
C           with one variable for each pollutant being projected
        IF( POLSFLAG ) THEN
            DO J = 1, NVPROJ
                VNAMESET( J )= PNAMPROJ( J )  ! Lowercase used to permit inv data named "PFAC"
                VTYPESET( J )= M3REAL
                VUNITSET( J )= 'n/a'
                VDESCSET( J )= 'Projection factor for ' // PNAMPROJ( J )
            END DO

C.........  If no pollutant-specific assignments, then set up projection
C           matrix to reflect that all pollutants affected the 
        ELSE
            J = 1
            VNAMESET( J )= 'pfac'  ! Lowercase used to permit inv data named "PFAC"
            VTYPESET( J )= M3REAL
            VUNITSET( J )= 'n/a'
            VDESCSET( J )= 'Projection factor'
        END IF

        MESG = 'Enter logical name for projection matrix...'
        CALL M3MSG2( MESG )

C.........  Open projection matrix.
C.........  Using NAMBUF is needed for HP to ensure string length consistencies

        MESG = 'I/O API PROJECTION MATRIX'

        NAMBUF = PROMPTSET( MESG, FSUNKN3, CRL // 'PMAT', 
     &                      PROGNAME )
        PNAME = NAMBUF

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

        END SUBROUTINE OPENPMAT
