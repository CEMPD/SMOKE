
        SUBROUTINE RDUSTATI( UINAME, UENAME, UPNAME )

C***********************************************************************
C
C  DESCRIPTION:
C     Subroutine Rdustati reads the uncertainty statistical information
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
C...........   This module contains the information about the source category
        USE MODINFO

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

        CHARACTER(*), INTENT(IN OUT) :: UINAME ! I/O API uncertainty inven input
        CHARACTER(*), INTENT(IN OUT) :: UENAME ! I/O API uncertainty empirical input
        CHARACTER(*), INTENT(IN OUT) :: UPNAME ! I/O API uncertainty parametric input

C.........  Local variables
        INTEGER       F, H, I, J, K, L, L2, S, T, U  ! counters and indices
        INTEGER       IOS                     ! i/o status

        LOGICAL    :: EFLAG = .FALSE.    ! true: error found
        LOGICAL    :: FIRSTIME = .TRUE.  ! true: first time the routine called

        CHARACTER*300 MESG                    ! message buffer
        CHARACTER*16 :: NAMBUF = ' '          ! temporary buffer
        CHARACTER*16 :: VARBUF = ' '          ! temporary buffer
        CHARACTER*16 :: PROGNAME = 'RDUSTATI' ! program name

C***********************************************************************
C uncertainty constructs
        INTEGER, ALLOCATABLE :: INUNSRC( : )  ! array for uncertainty source IDs

C***********************************************************************
C   begin body of subroutine RDEFACS

C.........  Get the number of simulations to draw from sampling
        MESG = 'Realizations to draw'
        DRAWRLZN = ENVINT( 'SMK_REALIZATIONS', MESG, 1 , IOS )

        IF ( ALLOCATED( APRCH   ) ) DEALLOCATE( APRCH   ) 
        IF ( ALLOCATED( EPTYP   ) ) DEALLOCATE( EPTYP   ) 
        IF ( ALLOCATED( EMFVAL  ) ) DEALLOCATE( EMFVAL  )
        IF ( ALLOCATED( PARMS   ) ) DEALLOCATE( PARMS   )
        IF ( ALLOCATED( PROBVAL ) ) DEALLOCATE( PROBVAL )
        IF ( ALLOCATED( METHOD  ) ) DEALLOCATE( METHOD  )  
        IF ( ALLOCATED( NUMEP   ) ) DEALLOCATE( NUMEP   )
        IF ( ALLOCATED( SRCNUM  ) ) DEALLOCATE( SRCNUM  )
        IF ( ALLOCATED( UNCIDX  ) ) DEALLOCATE( UNCIDX  )
        IF ( ALLOCATED( USTAT   ) ) DEALLOCATE( USTAT   )

        ALLOCATE( SRCNUM( NSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'SRCNUM', PROGNAME )
        SRCNUM = IMISS3  ! Array

        ALLOCATE( USTAT( NSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'USTAT', PROGNAME )
        USTAT = .FALSE.  ! Array

C.........  Prompt for and open I/O API output file(s)...
        MESG = 'Reading uncertainty input file...'
        NAMBUF = PROMPTMFILE( MESG, FSREAD3, CRL // 'UCOUT', PROGNAME )
        UINAME = NAMBUF

C.........  Get header description of non-diurnal file to get date for reading
        IF( .NOT. DESC3( UINAME ) ) THEN
            MESG = 'Could not get description of file "' 
     &              // UINAME // '"'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

        UNSRC = NROWS3D
        UNVAR  = MAX( 0,GETIFDSC( FDESC3D,'/NON POLLUTANT/',.FALSE.) )
        UNIPOL = MAX( 0,GETIFDSC( FDESC3D,'/UNCERT POLL/',.TRUE.) )
        UNPPOL = MAX( 0,GETIFDSC( FDESC3D,'/PER UNCERT POLL/',.FALSE.) )
        UNIACT = MAX( 0,GETIFDSC( FDESC3D,'/UNCERT ACT/',.FALSE.) )
        UNPACT = MAX( 0,GETIFDSC( FDESC3D,'/PER UNCERT ACT/',.FALSE.) )
        UNIPPA = UNIPOL + UNIACT

C.........  Allocate memory for emission factors and initialize to zero
        ALLOCATE( INUNSRC( UNSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'INUNSRC', PROGNAME )
        IF ( .NOT. READ3( UINAME, 'SRCNUM' , ALLAYS3, 0, 0,
     &                    INUNSRC ) ) THEN      
            WRITE( MESG,94010 ) 
     &             'ERROR: Could not read ' // UINAME // ' for SRCNUM'
            CALL M3MSG2( MESG )
            EFLAG = .TRUE.
        END IF

C.........  Make uncertainty sources 0 for output
        WRITE( MESG,94010 ) 'DEBUG: uncertainty rows read =', UNSRC
        WRITE( MESG,94010 ) 'DEBUG: sources read =', NSRC
        DO I = 1, UNSRC
            SRCNUM( INUNSRC( I ) ) = I
            USTAT( INUNSRC( I ) ) = .TRUE.
        END DO

C.........  Compute total number of output variables and allocate arrays
C           for uncertainty pollutants and activities
        ALLOCATE( UACTVTY( NIACT ), STAT=IOS )
        CALL CHECKMEM( IOS, 'UACTVTY', PROGNAME )
        UACTVTY = .FALSE.  ! array initialization

        ALLOCATE( UEINAM( NIPOL ), STAT=IOS )
        CALL CHECKMEM( IOS, 'UEINAM', PROGNAME )
        UEINAM = .FALSE.  ! array initialization

        ALLOCATE( UEFINAM( NIPOL ), STAT=IOS )
        CALL CHECKMEM( IOS, 'UEFINAM', PROGNAME )
        UEFINAM = .FALSE.  ! array initialization

        ALLOCATE( UEAREAD( NIPPA ), STAT=IOS )
        CALL CHECKMEM( IOS, 'UEAREAD', PROGNAME )
        UEAREAD = .FALSE.  ! array initialization

        U = 0
C.........  Load emission factors pollutants and activities
        DO I = 1, NIPOL
            L = LEN_TRIM( EINAM( I ) )
            VARBUF = 'MTH_EF_' // EINAM( I )( 1:L )
            K = LEN_TRIM( VARBUF )
            DO S = 1, NVARS3D
                IF ( VNAME3D( S ) .EQ. VARBUF ) THEN
                    U = U + 1
                    UEINAM( I ) = .TRUE.
                    UEAREAD( I ) = .TRUE.
                    EXIT
                END IF
            END DO
        END DO

        IF ( U .NE. UNIPOL) THEN
            MESG = 'ERROR: The number of uncertainty pollutants ' //
     &             'does not match the number in the uncertainty file' 
            CALL M3MSG2( MESG )
        END IF

        DO I = 1, NIACT
            L = LEN_TRIM( ACTVTY( I ) )
            VARBUF = 'MTH_' // ACTVTY( I )( 1:L )
            K = LEN_TRIM( VARBUF )
            DO S = 1, NVARS3D
                IF ( VNAME3D( S ) .EQ. VARBUF ) THEN
                    U = U + 1
                    UACTVTY( I ) = .TRUE.
                    UEAREAD( I + NIPOL ) = .TRUE.
                    EXIT
                END IF
            END DO
        END DO

        IF ( U .NE. UNIPPA) THEN
            MESG = 'ERROR: The number of uncertainty activities ' //
     &             'does not match the number in the uncertainty file' 
            CALL M3MSG2( MESG )
        END IF


c bug: will need to add a loop for uncertainty pollutants
c this is different from emissions factors but for the time being they are the same.
         UEFINAM = UEINAM
c
c        DO I = 1, NIPOL
c            L = LEN_TRIM( EINAM( I ) )
c            VARBUF = 'MTH_' // EINAM( I )( 1:L )
c            K = LEN_TRIM( VARBUF )
c            DO S = 1, NVARS3D
c                IF ( VNAME3D( S ) .EQ. VARBUF ) THEN
c                    U = U + 1
c                    UEINAM( I ) = .TRUE.
c                    UEAREAD( I ) = .TRUE.
c                    EXIT
c                END IF
c            END DO
c        END DO


C.........  Get method of sampling empirical(0) or parametric(1)
        ALLOCATE( METHOD( UNSRC, UNIPPA ), STAT=IOS )
        CALL CHECKMEM( IOS, 'METHOD', PROGNAME )
        METHOD = IMISS3 ! array initilization

        ALLOCATE( EPTYP( UNSRC, UNIPPA ), STAT=IOS )
        CALL CHECKMEM( IOS, 'EPTYP', PROGNAME )
        EPTYP = IMISS3 ! array initilization

        ALLOCATE( NUMEP( UNSRC, UNIPPA ), STAT=IOS )
        CALL CHECKMEM( IOS, 'NUMEP', PROGNAME )
        NUMEP = IMISS3 ! array initilization

        ALLOCATE( UNCIDX( UNSRC, UNIPPA ), STAT=IOS )
        CALL CHECKMEM( IOS, 'UNCIDX', PROGNAME )
        UNCIDX = IMISS3 ! array initilization

        ALLOCATE( APRCH( UNSRC, UNIPPA ), STAT=IOS )
        CALL CHECKMEM( IOS, 'APRCH', PROGNAME )
        APRCH = IMISS3 ! array initilization

        ALLOCATE( SAMPGEN( UNSRC, UNIPPA ), STAT=IOS )
        CALL CHECKMEM( IOS, 'SAMPGEN', PROGNAME )
        SAMPGEN = AMISS3 ! array initilization

C.........  Get emission factors information for pollutants
C BUG: will also need to add code for uncertainty pollutants and activities 
        S = 0
        DO I = 1, NIPPA ! BUG: may need to change to access pollutants array!

            IF ( .NOT. ( UEAREAD( I ) ) ) CYCLE
            S = S + 1

            L = LEN_TRIM( EINAM( I ) )
            VARBUF = 'MTH_EF_' // EINAM(I)( 1:L ) 
            L2 = LEN_TRIM( VARBUF )
            IF ( .NOT. READ3( UINAME, VARBUF( 1:L2 ), ALLAYS3, 0, 0,
     &                        METHOD( 1, S ) ) ) THEN
                MESG = 'ERROR: Could not read '//UINAME// ' for METHOD'
                CALL M3MSG2( MESG )
                EFLAG = .TRUE.
            END IF

            VARBUF = 'TYP_EF_' // EINAM(I)( 1:L )
            L2 = LEN_TRIM( VARBUF )
            IF ( .NOT. READ3( UINAME, VARBUF( 1:L2 ), ALLAYS3, 0, 0,
     &                        EPTYP( 1, S ) ) ) THEN
                MESG = 'ERROR: Could not read '//UINAME// ' for EPTYP'
                CALL M3MSG2( MESG )
                EFLAG = .TRUE.
            END IF

            VARBUF = 'NEP_EF_' // EINAM(I)( 1:L )
            L2 = LEN_TRIM( VARBUF )
            IF ( .NOT. READ3( UINAME, VARBUF( 1:L2 ), ALLAYS3, 0, 0,
     &                        NUMEP( 1, S ) ) ) THEN
                MESG = 'ERROR: Could not read '//UINAME// ' for NUMEP'
                CALL M3MSG2( MESG )
                EFLAG = .TRUE.
            END IF

            VARBUF = 'UIX_EF_' // EINAM(I)( 1:L )
            L2 = LEN_TRIM( VARBUF )
            IF ( .NOT. READ3( UINAME, VARBUF( 1:L2 ), ALLAYS3, 0, 0,
     &                        UNCIDX( 1, S ) ) ) THEN
                MESG = 'ERROR: Could not read '//UINAME// ' for UNCIDX'
                CALL M3MSG2( MESG )
                EFLAG = .TRUE.
            END IF

            VARBUF = 'APR_EF_' // EINAM(I)( 1:L )
            L2 = LEN_TRIM( VARBUF )
            IF ( .NOT. READ3( UINAME, VARBUF( 1:L2 ), ALLAYS3, 0, 0,
     &                        APRCH( 1, S ) ) ) THEN
                MESG = 'ERROR: Could not read '//UINAME// ' for APRCH'
                CALL M3MSG2( MESG )
                EFLAG = .TRUE.
            END IF

        END DO ! end loading METHOD, EPTYP, NUMEP, UNCIDX

C.........  Prompt for and open I/O API empirical file(s)...
        MESG = 'Reading uncertainty input empirical file...'
        NAMBUF = PROMPTMFILE( MESG, FSREAD3, CRL // 'UCEOUT', PROGNAME )
        UENAME = NAMBUF

C.........  Get header description of the uncertainty empirical statistics file
        IF ( .NOT. DESC3( UENAME ) ) THEN
            MESG = 'Could not get description of file "' 
     &              // UENAME // '"'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

C.........  Get uncertainty statistics information
C.........  Get empirical statistics information
        ALLOCATE( EMFVAL( NROWS3D, NCOLS3D ), STAT=IOS )
        CALL CHECKMEM( IOS, 'EMFVAL', PROGNAME )
        ALLOCATE( PROBVAL( NROWS3D, NCOLS3D ), STAT=IOS )
        CALL CHECKMEM( IOS, 'PROBVAL', PROGNAME )
        VARBUF = "EMFVAL"
        CALL GETDTPOS(UENAME, VARBUF, EMFVAL )
        VARBUF = "PROBVAL"
        CALL GETDTPOS(UENAME, VARBUF, PROBVAL )

C.........  Prompt for and open I/O API parametric file(s)...
        MESG = 'Reading uncertainty input parametric file...'
        NAMBUF = PROMPTMFILE( MESG, FSREAD3, CRL // 'UCPOUT', PROGNAME )
        UPNAME = NAMBUF

C.........  Get header description of the uncertainty empirical statistics file
        IF( .NOT. DESC3( UPNAME ) ) THEN
            MESG = 'Could not get description of file "' 
     &              // UPNAME // '"'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF
C.........  Get paramertric statistics information
        ALLOCATE( PARMS( NROWS3D, NCOLS3D ), STAT=IOS )
        CALL CHECKMEM( IOS, 'PARMS', PROGNAME )
        CALL GETDTPOS(UPNAME, "PARMS", PARMS )

C******************  FORMAT  STATEMENTS   ******************************

C.........  Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

C.........  test buffering formats............ 94xxx

98010   FORMAT( 10( A, :, F8.4, :, 1X ) )

C******************  INTERNAL SUBPROGRAMS  *****************************
        CONTAINS

            SUBROUTINE GETDTPOS( FNAME, FDATA, FVARIN )

            CHARACTER*16, INTENT(IN) :: FNAME       ! file name buffer
            CHARACTER(*), INTENT(IN) :: FDATA       ! file data to read
            REAL,        INTENT(OUT) :: FVARIN(:, :)! file variable buffer for input

            INTEGER               I, J, L           ! temporary work variables

            REAL, ALLOCATABLE :: TMPWORK( :, : )    ! temporary work matrix 

C*************************  BEGIN CODE  ********************************

            ALLOCATE( TMPWORK( NCOLS3D, NROWS3D ), STAT=IOS )
            CALL CHECKMEM( IOS, 'TMPWORK', PROGNAME )
  
            IF ( .NOT. READ3( FNAME, FDATA, ALLAYS3, 0, 0,
     &                        TMPWORK( 1, 1 ) ) ) THEN
                MESG = 'ERROR: Could not read ' //FNAME// ' for '//FDATA
                CALL M3MSG2( MESG )
                EFLAG = .TRUE.
            END IF

            IF ( NROWS3D .GT. 1 .AND. NCOLS3D .GT. 1 ) THEN
                DO I = 1, NROWS3D
                    DO J = 1, NCOLS3D
                        FVARIN( I, J ) = TMPWORK( J, I )
                    END DO
                END DO
            END IF

            DEALLOCATE( TMPWORK )

C.........  Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

            END SUBROUTINE GETDTPOS

        END SUBROUTINE RDUSTATI
