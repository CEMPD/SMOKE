
        SUBROUTINE OPENSMAT( ENAME, SFLAG, LFLAG, NOPOL, MXSPEC, 
     &                       EALLOUT, EAIDX, SPCNAMES, MOLUNITS, 
     &                       SDEV, SNAME, LNAME, SVNAMES, LVNAMES )

C***********************************************************************
C  subroutine body starts at line 106
C
C  DESCRIPTION:
C      Open the mass-based and mole-based speciation matrices
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C     Created 2/99 by M. Houyoux
C
C***************************************************************************
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
        USE MODINFO, ONLY: CATEGORY, CATLEN, CRL, NIPPA

C.........  This module is required by the FileSetAPI
        USE MODFILESET

        IMPLICIT NONE

C...........   INCLUDES
        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
        INCLUDE 'SETDECL.EXT'   !  FileSetAPI function declarations

C...........   EXTERNAL FUNCTIONS and their descriptions:
        CHARACTER*2            CRLF
        INTEGER                FINDC
        CHARACTER(LEN=IODLEN3) GETCFDSC
        INTEGER                PROMPTFFILE
        CHARACTER*16           PROMPTMFILE
        CHARACTER*16           VERCHAR

        EXTERNAL        CRLF, FINDC, GETCFDSC, PROMPTFFILE, 
     &                  PROMPTMFILE, VERCHAR

C.........  SUBROUTINE ARGUMENTS
        CHARACTER(*), INTENT (IN) :: ENAME      ! emissions inven logical name
        LOGICAL     , INTENT (IN) :: SFLAG      ! true: open mass-based file
        LOGICAL     , INTENT (IN) :: LFLAG      ! true: open mole-based file
        INTEGER     , INTENT (IN) :: NOPOL      ! no. output pollutants
        INTEGER     , INTENT (IN) :: MXSPEC     ! max no. of spec per pol
        CHARACTER(*), INTENT (IN) :: EALLOUT ( NIPPA ) ! output pol/emistypes
        INTEGER     , INTENT (IN) :: EAIDX   ( NIPPA ) ! index to SPCNAMES
        CHARACTER(*), INTENT (IN) :: SPCNAMES( MXSPEC, NOPOL ) ! model spec nams
        CHARACTER(*), INTENT (IN) :: MOLUNITS( MXSPEC, NOPOL ) ! mole-based unts
        INTEGER     , INTENT(OUT) :: SDEV            ! suplmt file unit no.
        CHARACTER(*), INTENT(OUT) :: SNAME           ! mass-based spec file name 
        CHARACTER(*), INTENT(OUT) :: LNAME           ! mole-based spec file name
        CHARACTER(*), INTENT(OUT) :: SVNAMES( MXSPEC, NIPPA )   ! mass out vars
        CHARACTER(*), INTENT(OUT) :: LVNAMES( MXSPEC, NIPPA )   ! mole out vars
      
C...........   LOCAL PARAMETERS
        CHARACTER*50, PARAMETER :: CVSW = '$Name$'  ! CVS revision tag

C.........  Count of species per inventory pollutant/emission type
        INTEGER    NSPEC( NIPPA )

C.........  Other local variables
        INTEGER          I, J, K, V     !  counters and indices
        INTEGER          IOS            !  I/O status

        INTEGER          ICNT     ! cntr for the total number of output vars
        INTEGER          NCNT     ! cntr for number of species per inv pol

        CHARACTER*300    MESG     ! message buffer

        CHARACTER(LEN=NAMLEN3) NAMBUF     ! file name buffer
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
        DO K = 1, NIPPA

            V = EAIDX( K )

            NCNT = 0
            DO J = 1, MXSPEC

C.................  End inner loop if species is blank
                IF( SPCNAMES( J,V ) .EQ. ' ' ) EXIT

C.................  Count total number of output variables
                ICNT = ICNT + 1

                NCNT = NCNT + 1
                WRITE( SVNAMES( J,K ), '(A4,I3.3)' ) 'SVAR', ICNT
                WRITE( LVNAMES( J,K ), '(A4,I3.3)' ) 'SVAR', ICNT

            END DO

            NSPEC( K ) = NCNT

        END DO

C.........  Set up file header(s) for opening I/O API output(s). Base this on
C           inventory header...

C.........  Get header information from inventory file
        IF ( .NOT. DESCSET( ENAME,-1 ) ) THEN
            MESG = 'Could not get description of file "' 
     &             // ENAME( 1:LEN_TRIM( ENAME ) ) // '".'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

        IFDESC2 = GETCFDSC( FDESC3D, '/FROM/', .TRUE. )
        IFDESC3 = GETCFDSC( FDESC3D, '/VERSION/', .TRUE. )

        FDESC3D = ' '   ! array

        FDESC3D( 1 ) = CATEGORY( 1:CATLEN ) // ' speciation matrix'
        FDESC3D( 2 ) = '/FROM/ '    // PROGNAME
        FDESC3D( 3 ) = '/VERSION/ ' // VERCHAR( CVSW )

        FDESC3D( 11 ) = '/INVEN FROM/ ' // IFDESC2
        FDESC3D( 12 ) = '/INVEN VERSION/ ' // IFDESC3

C.........  Resset output number of variables
        NVARSET = ICNT

C.........  Deallocate, then allocate, output arrays
        DEALLOCATE( VNAMESET, VTYPESET, VUNITSET, VDESCSET )
        ALLOCATE( VNAMESET( NVARSET ), STAT=IOS )
        CALL CHECKMEM( IOS, 'VNAMESET', PROGNAME )
        ALLOCATE( VTYPESET( NVARSET ), STAT=IOS )
        CALL CHECKMEM( IOS, 'VTYPESET', PROGNAME )
        ALLOCATE( VUNITSET( NVARSET ), STAT=IOS )
        CALL CHECKMEM( IOS, 'VUNITSET', PROGNAME )
        ALLOCATE( VDESCSET( NVARSET ), STAT=IOS )
        CALL CHECKMEM( IOS, 'VDESCSET', PROGNAME )

C.........  Also deallocate the number of variables per file so
C           that this will be set automatically by openset
        DEALLOCATE( VARS_PER_FILE )

C.........  Set up variable descriptions that will be used to indicate the 
C           inventory pollutant and model species names

        VTYPESET = 0    ! array initialization
        VNAMESET = ' '  ! array initialization
        VUNITSET = ' '  ! array initialization
        VDESCSET = ' '  ! array initialization

        I = 0
        DO K = 1, NIPPA

            V = EAIDX( K )

            DO J = 1, NSPEC( K )

                I = I + 1
                VDESCSET( I ) = EALLOUT( K )// SPJOIN// SPCNAMES( J,V )
                VTYPESET( I ) = M3REAL

            END DO
        END DO

C.........  Set up variables specifically for mass-based file, and open it
        IF( SFLAG ) THEN

            FDESC3D( 4 ) = '/SMATTYPE/ ' // MASSSTR

            I = 0
            DO K = 1, NIPPA
                DO J = 1, NSPEC( K )

                    I = I + 1
                    VNAMESET( I ) = SVNAMES( J,K ) 
                    VUNITSET( I ) = SMASUNIT

                END DO
            END DO

C.............  Open with NAMBUF for HP
            NAMBUF = PROMPTSET( 
     &        'Enter logical name for MASS-BASED SPECIATION MATRIX',
     &        FSUNKN3, CRL // 'SMAT_S', PROGNAME )
            SNAME = NAMBUF

        END IF

C.........  Set up variables specifically for mole-based file, and open it
        IF( LFLAG ) THEN

            FDESC3D( 4 ) = '/SMATTYPE/ ' // MOLESTR

            I = 0
            DO K = 1, NIPPA

                V = EAIDX( K )

                DO J = 1, NSPEC( K )

                    I = I + 1
                    VNAMESET( I ) =  LVNAMES( J,K ) 
                    VUNITSET( I ) = MOLUNITS( J,V )

                END DO
            END DO

C.............  Open with NAMBUF for HP
            NAMBUF = PROMPTSET( 
     &        'Enter logical name for MOLE-BASED SPECIATION MATRIX',
     &        FSUNKN3, CRL // 'SMAT_L', PROGNAME )
            LNAME = NAMBUF

        END IF

C.........  Open supplemental speciation file
        MESG = 'Enter logical name for the SPECIATION SUPPLEMENTAL '//
     &         'file'
        SDEV = PROMPTFFILE( MESG, .FALSE., .TRUE., 
     &                      CRL // 'SSUP', PROGNAME )


        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

        END SUBROUTINE OPENSMAT
