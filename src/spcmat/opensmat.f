        SUBROUTINE OPENSMAT( CATEGORY, ENAME, SFLAG, LFLAG, 
     &                       MXSPEC, NIPOL, EINAM, SPCNAMES, 
     &                       SNAME, LNAME, SVNAMES, LVNAMES )

C***********************************************************************
C  subroutine body starts at line
C
C  DESCRIPTION:
C      Abridge full pollutants table to a data structure that can be used
C      in ASGNSPRO
C
C  PRECONDITIONS REQUIRED:
C      Expects cross-reference tables to be set to IMISS3 if not defined
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C     Created 2/99 by M. Houyoux
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

        IMPLICIT NONE

C...........   INCLUDES
        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
        INCLUDE 'FDESC3.EXT'    !  I/O API file description data structures.

C...........   EXTERNAL FUNCTIONS and their descriptions:
        CHARACTER*2            CRLF
        INTEGER                FINDC
        CHARACTER(LEN=IODLEN3) GETCFDSC
        CHARACTER*16           PROMPTMFILE
        CHARACTER*16           VERCHAR

        EXTERNAL        CRLF, FINDC, GETCFDSC, PROMPTMFILE, VERCHAR

C.........  SUBROUTINE ARGUMENTS
        CHARACTER(*), INTENT (IN) :: CATEGORY   ! source category
        CHARACTER(*), INTENT (IN) :: ENAME      ! emissions inven logical name
        LOGICAL     , INTENT (IN) :: SFLAG      ! true: open mass-based file
        LOGICAL     , INTENT (IN) :: LFLAG      ! true: open mole-based file
        INTEGER     , INTENT (IN) :: MXSPEC     ! max no. of spec per pol
        INTEGER     , INTENT (IN) :: NIPOL      ! number of inventory pols
        CHARACTER(*), INTENT (IN) :: EINAM( NIPOL ) ! names of actual pols
        CHARACTER(*), INTENT (IN) :: SPCNAMES( MXSPEC, NIPOL ) ! model spec nams
        CHARACTER(*), INTENT(OUT) :: SNAME           ! mass-based spec file name 
        CHARACTER(*), INTENT(OUT) :: LNAME           ! mole-based spec file name
        CHARACTER(*), INTENT(OUT) :: SVNAMES( MXSPEC, NIPOL )   ! mass out vars
        CHARACTER(*), INTENT(OUT) :: LVNAMES( MXSPEC, NIPOL )   ! mole out vars
      
C...........   LOCAL PARAMETERS
        CHARACTER*50  SCCSW          ! SCCS string with version number at end

        PARAMETER   ( SCCSW   = '@(#)$Id$'
     &              )

C.........  Count of species per inventory pollutant
        INTEGER    NSPEC( NIPOL )

C.........  Other local variables
        INTEGER          I, J, K, V     !  counters and indices

        INTEGER          ICNT     ! cntr for the total number of output vars
        INTEGER          NCNT     ! cntr for number of species per inv pol

        CHARACTER*300    MESG     ! message buffer

        CHARACTER(LEN=IODLEN3) IFDESC2, IFDESC3 ! fields 2 & 3 from inven FDESC
        CHARACTER(LEN=SPNLEN3) PCODE      ! current speciation profile code
        CHARACTER(LEN=SPNLEN3) PREVCODE   ! previous speciation profile code

        CHARACTER*16 :: PROGNAME = 'OPENSMAT' ! program name

C***********************************************************************
C   begin body of subroutine OPENSMAT

C.........  Initialize variable names
        SVNAMES = ' '  ! Array
        LVNAMES = ' '  ! Array

C.........  Create the output variable names and count them.  Handle the special
C           case where the total number of variables exceeds the I/O API max.
C.........  The names of the output variables have been set up so that it will 
C           be easy to make the mass-based and the mole-based ones different.
        ICNT = 0
        DO V = 1, NIPOL

            NCNT = 0
            DO J = 1, MXSPEC

C.................  End loop early if species is blank
                IF( SPCNAMES( J,V ) .EQ. ' ' ) EXIT

C.................  Count total number of output variables
                ICNT = ICNT + 1

C.................  Check total number of output variables with I/O API max
                IF( ICNT .LE. MXVARS3 ) THEN

                    NCNT = NCNT + 1
                    WRITE( SVNAMES( J,V ), '(A4,I3.3)' ) 'SVAR', ICNT
                    WRITE( LVNAMES( J,V ), '(A4,I3.3)' ) 'SVAR', ICNT

                ENDIF

            END DO

            NSPEC( V ) = NCNT

        END DO

C.........  Print Error if number of variables is passed maximum.
C.........  DO NOT end program here because it will be ended when the write
C           attempt is made for these extra variables.
        IF( ICNT .GT. MXVARS3 ) THEN
            WRITE( MESG, 94010 ) 
     &             'ERROR: maximum I/O API variables exceeded:' //
     &             CRLF() // BLANK10 // 'Max: ', MXVARS3, 'Actual:',ICNT
            CALL M3MSG2( MESG )

            ICNT = MXVARS3
        ENDIF

C.........  Set up file header(s) for opening I/O API output(s). Base this on
C           inventory header...

C.........  Get header information from inventory file

        IF ( .NOT. DESC3( ENAME ) ) THEN
            MESG = 'Could not get description of file "' 
     &             // ENAME( 1:LEN_TRIM( ENAME ) ) // '".'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

        IFDESC2 = GETCFDSC( FDESC3D, '/FROM/' )
        IFDESC3 = GETCFDSC( FDESC3D, '/VERSION/' )

        NVARS3D = ICNT

        FDESC3D = ' '   ! array

        FDESC3D( 1 ) = CATEGORY( 1:LEN_TRIM( CATEGORY ) ) //
     &                 ' speciation matrix'
        FDESC3D( 2 ) = '/FROM/ '    // PROGNAME
        FDESC3D( 3 ) = '/VERSION/ ' // VERCHAR( SCCSW )

        FDESC3D( 11 ) = '/PNTS FROM/ ' // IFDESC2
        FDESC3D( 12 ) = '/PNTS VERSION/ ' // IFDESC3


C.........  Set up variable descriptions that will be used to indicate the 
C           inventory pollutant and model species names

        VDESC3D = ' '  ! array initialization
        VTYPE3D = 0    ! array initialization

        I = 0
        DO V = 1, NIPOL
            DO J = 1, NSPEC( V )

                I = I + 1
                VDESC3D( I ) = EINAM( V ) // '_' // SPCNAMES( J,V )               
                VTYPE3D( I ) = M3REAL

            END DO
        END DO

C.........  Set up variables specifically for mass-based file, and open it
        IF( SFLAG ) THEN

            FDESC3D( 4 ) = '/SMATTYPE/ ' // ' Mass'

            I = 0
            DO V = 1, NIPOL
                DO J = 1, NSPEC( V )

                    I = I + 1
                    VNAME3D( I ) =  SVNAMES( J,V ) 
                    UNITS3D( I ) = 'gm/ton'

                END DO
            END DO

            SNAME = PROMPTMFILE( 
     &        'Enter logical name for MASS-BASED SPECIATION MATRIX',
     &        FSUNKN3, 'PSMAT_S', PROGNAME )

        ENDIF

C.........  Set up variables specifically for mole-based file, and open it
        IF( LFLAG ) THEN

            FDESC3D( 4 ) = '/SMATTYPE/ ' // ' Mole'

            I = 0
            DO V = 1, NIPOL
                DO J = 1, NSPEC( V )

                    I = I + 1
                    VNAME3D( I ) =  LVNAMES( J,V ) 
                    UNITS3D( I ) = 'mole/ton'

                END DO
            END DO

            LNAME = PROMPTMFILE( 
     &        'Enter logical name for MOLE-BASED SPECIATION MATRIX',
     &        FSUNKN3, 'PSMAT_L', PROGNAME )

        ENDIF

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

        END SUBROUTINE OPENSMAT
