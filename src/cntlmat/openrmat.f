        SUBROUTINE OPENRMAT( ENAME, RPOL, SFLAG, LFLAG, 
     &                       BYEARIN, PYEAR, NSREAC, NMSPC, SPECIES, 
     &                       SDEV, SNAME, LNAME, SVNAMES, LVNAMES )

C***********************************************************************
C  subroutine body starts at line
C
C  DESCRIPTION:
C      Open the mass-based and mole-based reactivity matrices
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C     Created 3/99 by M. Houyoux
C
C****************************************************************************/
C
C Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
C                System
C File: @(#)$Id$
C
C COPYRIGHT (C) 1998, MCNC--North Carolina Supercomputing Center
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
        INTEGER                PROMPTFFILE
        CHARACTER*16           PROMPTMFILE
        CHARACTER*16           VERCHAR

        EXTERNAL     CRLF, GETCFDSC, GETIFDSC, PROMPTFFILE, PROMPTMFILE, 
     &               VERCHAR

C...........   LOCAL PARAMETERS
        INTEGER     , PARAMETER :: NBASVAR = 4
        CHARACTER*50, PARAMETER :: SCCSW   = '@(#)$Id$' ! SCCS string w/ vers no.

C.........  SUBROUTINE ARGUMENTS
        CHARACTER(*), INTENT (IN) :: ENAME      ! emissions inven logical name
        CHARACTER(*), INTENT (IN) :: RPOL       ! pollutant for matrices
        LOGICAL     , INTENT (IN) :: SFLAG      ! true: open mass-based file
        LOGICAL     , INTENT (IN) :: LFLAG      ! true: open mole-based file
        INTEGER     , INTENT (IN) :: BYEARIN    ! base year of proj factors
        INTEGER     , INTENT (IN) :: PYEAR      ! projected year of proj factors
        INTEGER     , INTENT (IN) :: NSREAC     ! number of reactvty sources
        INTEGER     , INTENT (IN) :: NMSPC      ! number of reactivity species
        CHARACTER(*), INTENT (IN) :: SPECIES( NMSPC ) ! model spec nams
        INTEGER     , INTENT(OUT) :: SDEV       ! unit number of supplement file
        CHARACTER(*), INTENT(OUT) :: SNAME      ! mass-based spec file name 
        CHARACTER(*), INTENT(OUT) :: LNAME      ! mole-based spec file name
        CHARACTER(*), INTENT(OUT) :: SVNAMES( NMSPC )   ! species mass out vars
        CHARACTER(*), INTENT(OUT) :: LVNAMES( NMSPC )   ! species mole out vars
      
C.........  Other local variables
        INTEGER          I, J           !  counters and indices

        CHARACTER*300          MESG     ! message buffer

        CHARACTER(LEN=IOVLEN3) CPOL     ! pollutant name buffer
        CHARACTER(LEN=IODLEN3) IFDESC2, IFDESC3 ! fields 2 & 3 from inven FDESC
        CHARACTER(LEN=IOVLEN3) UNITS    ! emissions units

        CHARACTER*16 :: PROGNAME = 'OPENRMAT' ! program name

C***********************************************************************
C   begin body of subroutine OPENRMAT

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

C.........  Initialize variable names for I/O API files
        SVNAMES = ' '  ! Array
        LVNAMES = ' '  ! Array

C.........  Initialize I/O API output file headers
        CALL HDRMISS3

C.........  Set I/O API header parms that need values
        NROWS3D = NSREAC
        NVARS3D = MIN( NBASVAR + NMSPC, MXVARS3 )
        NTHIK3D = NSRC

        FDESC3D( 1 ) = CATEGORY( 1:CATLEN ) // ' reactivity matrix'
        FDESC3D( 2 ) = '/FROM/ '    // PROGNAME
        FDESC3D( 3 ) = '/VERSION/ ' // VERCHAR( SCCSW )

        WRITE( FDESC3D( 5 ), '(A,I4)' ) '/BASE YEAR/ ', BYEARIN
        WRITE( FDESC3D( 6 ), '(A,I4)' ) '/PROJECTED YEAR/ ', PYEAR
        WRITE( FDESC3D( 7 ), '(A,I4)' ) '/SPECIES VARS/ ', NMSPC

        FDESC3D( 11 ) = '/INVEN FROM/ ' // IFDESC2
        FDESC3D( 12 ) = '/INVEN VERSION/ ' // IFDESC3

C.........  Set up non-speciation variables
        J = 1
        VNAME3D( J ) = 'SRCID'
        VTYPE3D( J ) = M3INT
        UNITS3D( J ) = 'n/a'
        VDESC3D( J ) = 'Inventory source ID (position)'
        J = J + 1

        VNAME3D( J )= 'REPEMIS'
        VTYPE3D( J )= M3REAL
        UNITS3D( J )= UNITS
        VDESC3D( J )= 'Reactivity base-year emissions' 
        J = J + 1

        VNAME3D( J )= 'PRJFAC'
        VTYPE3D( J )= M3REAL
        UNITS3D( J )= 'n/a'
        VDESC3D( J )= 'Reactivity projection factor'

        J = J + 1
        VNAME3D( J ) = 'MKTPEN'
        VTYPE3D( J ) = M3REAL
        UNITS3D( J ) = 'fraction'
        VDESC3D( J ) = 'Reactivity control market penetration'

C.........  Make sure program has not been modified improperly to cause a 
C           hard-to-detect error

        IF( J .NE. NBASVAR ) THEN
            MESG = 'INTERNAL ERROR: Number of variables NBASVAR is ' //
     &             'inconsistent with other subroutine modifications'
            CALL M3MSG2( MESG )
            CALL M3EXIT( PROGNAME, 0, 0, ' ', 2 )
        END IF

C.........  Set speciation variables for I/O API header
        CPOL = ADJUSTL( RPOL )
        DO I = 1, NMSPC

            J = J + 1

C.............  Set SVNAMES and LVNAMES here instead of VNAME3D, so that it 
C               will be easier to change these to have different name roots
C               if needed in the future
C.............  Check total number of output variables with I/O API max.  
            IF( J .LE. MXVARS3 ) THEN
                VTYPE3D( J ) = M3REAL
                VDESC3D( J ) = CPOL // SPJOIN // SPECIES( I )               
            ENDIF

            WRITE( SVNAMES( I ), '(A4,I3.3)' ) 'SVAR', I
            WRITE( LVNAMES( I ), '(A4,I3.3)' ) 'SVAR', I

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

        MESG = 'Enter logical name(s) for ...'
        CALL M3MSG2( MESG )

C.........  Set up variables specifically for mass-based file, and open it
        IF( SFLAG ) THEN

            FDESC3D( 4 ) = '/SMATTYPE/ ' // ' Mass'

            I = 0
            DO J = NBASVAR+1, NVARS3D

                I = I + 1
                VNAME3D( J ) =  SVNAMES( I ) 
                UNITS3D( J ) = 'gm/ton'

            END DO

            MESG = 'I/O API MASS-BASED REACTIVITY MATRIX for ' // RPOL

            SNAME = PROMPTMFILE( MESG, FSUNKN3, CRL // 'RMAT_S', 
     &                           PROGNAME )

        ENDIF

C.........  Set up variables specifically for mole-based file, and open it
        IF( LFLAG ) THEN

            FDESC3D( 4 ) = '/SMATTYPE/ ' // ' Mole'

            I = 0
            DO J = NBASVAR+1, NVARS3D

                I = I + 1
                VNAME3D( J ) =  LVNAMES( I ) 
                UNITS3D( J ) = 'mole/ton'

            END DO

            MESG = 'I/O API MOLE-BASED REACTIVITY MATRIX for ' // RPOL

            LNAME = PROMPTMFILE( MESG, FSUNKN3, CRL // 'RMAT_L', 
     &                           PROGNAME )

        ENDIF

C.........  Open the supplementary file (for SCCs and SPROFs)

        MESG = 'ASCII REACTIVITY MATRIX SUPPLEMENT file for ' // RPOL

        SDEV = PROMPTFFILE( MESG, .FALSE., .TRUE., CRL // 'RSUP', 
     &                      PROGNAME )

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

        END SUBROUTINE OPENRMAT
