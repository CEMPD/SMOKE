
        SUBROUTINE RDIUSRC( NSRC, UINAME )

C***********************************************************************
C
C  DESCRIPTION:
C     Subroutine Rdiusrc reads the uncertainty statistical information
C     regarding methods (distribution), parameters, ets.
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C     Created 3/02 by G. Cano
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

C.........  MODULES for public variables

C...........   This module contains the uncertainty arrays and variables
        USE MODUNCERT

        IMPLICIT NONE

C...........   INCLUDES
        INCLUDE 'EMCNST3.EXT'   ! emissions constants
        INCLUDE 'PARMS3.EXT'    ! I/O API parameters
        INCLUDE 'IODECL3.EXT'   ! I/O API function declarations
        INCLUDE 'FDESC3.EXT'    ! I/O API file description data structures

C...........   EXTERNAL FUNCTIONS and their descriptions:
        INTEGER     ENVINT
        INTEGER     FIND1
        INTEGER     INDEX1
        INTEGER     GETIFDSC
        INTEGER     SEC2TIME
        CHARACTER(LEN=NAMLEN3) PROMPTMFILE

        EXTERNAL    ENVINT, FIND1, INDEX1, SEC2TIME, PROMPTMFILE
     &              GETIFDSC
 
C...........   SUBROUTINE ARGUMENTS
        INTEGER,      INTENT(IN)        :: NSRC   ! number of sources
        CHARACTER(*), INTENT(IN OUT)    :: UINAME ! I/O API uncertainty inven input

C.........  Local variables
        INTEGER       F, H, I, J, K, L, L2, S, T, U  ! counters and indices
        INTEGER       IOS                     ! i/o status
        INTEGER, ALLOCATABLE :: INUNSRC( : )  ! array for uncertainty source IDs

        LOGICAL       LIOS               ! true: error found closing file
        LOGICAL    :: EFLAG = .FALSE.    ! true: error found
        LOGICAL    :: FIRSTIME = .TRUE.  ! true: first time the routine called

        CHARACTER*300 MESG                    ! message buffer
        CHARACTER*16 :: NAMBUF = ' '          ! temporary buffer
        CHARACTER*16 :: VARBUF = ' '          ! temporary buffer
        CHARACTER*16 :: PROGNAME = 'RDIUSRC' ! program name

C***********************************************************************
C   begin body of subroutine RDIUSRC

        IF ( ALLOCATED( SRCNUM  ) ) DEALLOCATE( SRCNUM  )
        IF ( ALLOCATED( USTAT   ) ) DEALLOCATE( USTAT   )
        IF ( ALLOCATED( UNAMES  ) ) DEALLOCATE( UNAMES  )

        ALLOCATE( SRCNUM( NSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'SRCNUM', PROGNAME )
        SRCNUM = IMISS3  ! Array

        ALLOCATE( USTAT( NSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'USTAT', PROGNAME )
        USTAT = .FALSE.  ! Array

C.........  Prompt for and open I/O API output file(s)...
        MESG = 'Reading uncertainty input file...'
        NAMBUF = PROMPTMFILE( MESG, FSREAD3, UCAT( 1:1 ) // 'UCOUT',
     &                        PROGNAME )
        UINAME = NAMBUF

C.........  Get header description of non-diurnal file to get date for reading
        IF( .NOT. DESC3( UINAME ) ) THEN
            MESG = 'Could not get description of file "' 
     &              // UINAME // '"'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

        UNSRC  = NROWS3D
        UNVAR  = MAX( 0,GETIFDSC( FDESC3D,'/NON POLLUTANT/',.FALSE.) )
        UNIPOL = MAX( 0,GETIFDSC( FDESC3D,'/UNCERT POLL/',.TRUE.) )
        UNPPOL = MAX( 0,GETIFDSC( FDESC3D,'/PER UNCERT POLL/',.FALSE.) )
        UNIACT = MAX( 0,GETIFDSC( FDESC3D,'/UNCERT ACT/',.FALSE.) )
        UNPACT = MAX( 0,GETIFDSC( FDESC3D,'/PER UNCERT ACT/',.FALSE.) )
        UNIPPA = UNIPOL + UNIACT

        UPINVAR = NVARS3D
        ALLOCATE( UNAMES( UPINVAR ), STAT=IOS )
        CALL CHECKMEM( IOS, 'UNAMES', PROGNAME )
        DO I = 1, UPINVAR
            UNAMES( I ) = VNAME3D( I )
        END DO

C.........  Allocate memory for emission factors and initialize to zero
        ALLOCATE( INSRC( UNSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'INSRC', PROGNAME )
        IF ( .NOT. READ3( UINAME, 'SRCNUM' , ALLAYS3, 0, 0,
     &                    INSRC ) ) THEN      
            WRITE( MESG,94010 ) 
     &             'ERROR: Could not read ' // UINAME // ' for SRCNUM'
            CALL M3MSG2( MESG )
            EFLAG = .TRUE.
        END IF

C.........  Make uncertainty sources 0 for output
        WRITE( MESG,94010 ) 'DEBUG: uncertainty rows read =', UNSRC
        WRITE( MESG,94010 ) 'DEBUG: sources read =', NSRC
        DO I = 1, UNSRC
            SRCNUM( INSRC( I ) ) = I
            USTAT( INSRC( I ) ) = .TRUE.
        END DO



        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C.........  Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

        END SUBROUTINE RDIUSRC
