
        SUBROUTINE RDMETHS( UINAME, UENAME, UPNAME )

C***********************************************************************
C
C  DESCRIPTION:
C     Subroutine Rdmeths reads the uncertainty statistical information
C     regarding methods (distributions), parameters, etc.
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
        INTEGER     GETIFDSC
        INTEGER     INDEX1
        CHARACTER(LEN=NAMLEN3) PROMPTMFILE

        EXTERNAL     GETIFDSC, INDEX1, PROMPTMFILE
 
C...........   SUBROUTINE ARGUMENTS
        CHARACTER(*), INTENT(IN)     :: UINAME ! I/O API uncertainty info 
        CHARACTER(*), INTENT(OUT)    :: UENAME ! I/O API uncertainty empirical input
        CHARACTER(*), INTENT(OUT)    :: UPNAME ! I/O API uncertainty parametric input

C.........  Local variables
        INTEGER       F, H, I, J, K, L, L2, S, T, U  ! counters and indices
        INTEGER       IOS                     ! i/o status

        LOGICAL    :: EFLAG = .FALSE.    ! true: error found
        LOGICAL    :: FIRSTIME = .TRUE.  ! true: first time the routine called

        CHARACTER*300 MESG                    ! message buffer
        CHARACTER*16 :: NAMBUF = ' '          ! temporary buffer
        CHARACTER*16 :: VARBUF = ' '          ! temporary buffer
        CHARACTER*16 :: PROGNAME = 'RDMETHS'  ! program name

C***********************************************************************
C   begin body of subroutine RDMETHS

        IF ( ALLOCATED( APRCH   ) ) DEALLOCATE( APRCH   ) 
        IF ( ALLOCATED( EPTYP   ) ) DEALLOCATE( EPTYP   ) 
        IF ( ALLOCATED( EMFVAL  ) ) DEALLOCATE( EMFVAL  )
        IF ( ALLOCATED( PARMS   ) ) DEALLOCATE( PARMS   )
        IF ( ALLOCATED( PROBVAL ) ) DEALLOCATE( PROBVAL )
        IF ( ALLOCATED( METHOD  ) ) DEALLOCATE( METHOD  )  
        IF ( ALLOCATED( NUMEP   ) ) DEALLOCATE( NUMEP   )
        IF ( ALLOCATED( UNCIDX  ) ) DEALLOCATE( UNCIDX  )

C.....  Compute total number of output variables and allocate arrays
C       for uncertainty pollutants and activities
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
C.........  Load uncertainty emissions
        DO I = 1, NIPOL
            L = LEN_TRIM( EINAM( I ) )
            VARBUF = 'MTH_' // EINAM( I )( 1:L )
            K = INDEX1( VARBUF, UPINVAR, UNAMES )
            IF ( K .GT. 0 ) THEN
                S = LEN_TRIM( VARBUF )
                U = U + 1
                UEINAM( I ) = .TRUE.
                UEAREAD( I ) = .TRUE.
            END IF
        END DO
      
        IF ( U .GT. UNIPOL) THEN
            MESG = 'ERROR: The number of uncertainty emissions ' //
     &             'exceeded the number in the inventory file'
            CALL M3MSG2( MESG )
        END IF

c BUG may need to incorporate for uncertainty emission factors.
C.........  Load uncertainty emission factors
c        DO I = 1, NIPOL
c            L = LEN_TRIM( EINAM( I ) )
c            VARBUF = 'MTH_EF_' // EINAM( I )( 1:L )
c            K = INDEX1( VARBUF, UPINVAR, UNAMES )
c            IF ( K .GT. 0 ) THEN
c                S = LEN_TRIM( VARBUF )
c                U = U + 1
c                UEFINAM( I ) = .TRUE.
c                UEAREAD( I ) = .TRUE.
c            END IF
c        END DO
c        IF ( U .GT. UNIPOL) THEN
c            MESG = 'ERROR: The number of uncertainty ' //
c     &              'emission factors does not match ' //
c     &              'the number in the inventory file'
c            CALL M3MSG2( MESG )
c        END IF

c BUG may need to incorporate for uncertainty activities.
c        DO I = 1, NIACT
c            L = LEN_TRIM( ACTVTY( I ) )
c            VARBUF = 'MTH_' // ACTVTY( I )( 1:L )
c            K = INDEX1( VARBUF, UPINVAR, UNAMES )
c            IF( VNAME3D( S ) .EQ. VARBUF ) THEN
c                U = U + 1
c                UACTVTY( I ) = .TRUE.
c                UEAREAD( I + NIPOL ) = .TRUE.
c            END IF
c        END DO
c        IF ( U .NE. UNIPPA) THEN
c            MESG = 'ERROR: The number of uncertainty activities ' //
c     &             'does not match the number in the inventory file' 
c            CALL M3MSG2( MESG )
c        END IF
          
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
C BUG: will also need to add code for uncertainty emission factors 
        S = 0
        DO I = 1, NIPPA ! BUG: may need to change to access pollutants array!
      
            IF ( .NOT. ( UEAREAD( I ) ) ) CYCLE
            S = S + 1
      
            L = LEN_TRIM( EINAM( I ) )
            VARBUF = 'MTH_' // EINAM(I)( 1:L ) 
            L2 = LEN_TRIM( VARBUF )
            IF ( .NOT. READ3( UINAME, VARBUF( 1:L2 ), ALLAYS3, 0, 0,
     &                        METHOD( 1, S ) ) ) THEN
                MESG = 'ERROR: Could not read ' // UINAME
     &                 // ' for METHOD'
                CALL M3MSG2( MESG )
                EFLAG = .TRUE.
            END IF


c bug: will need to add a loop for uncertainty pollutants
c      this may also be handled by adding an ELSE to the above IF statement
c this is different from emissions factors but for the time being they are the same.
c
c            DO I = 1, NIPOL
c                L = LEN_TRIM( EINAM( I ) )
c                VARBUF = 'MTH_EF' // EINAM( I )( 1:L )
c                K = LEN_TRIM( VARBUF )
c                DO S = 1, NVARS3D  ! nvars3d from inventory file AUCOUT
c                    IF ( VNAME3D( S ) .EQ. VARBUF ) THEN
c                        U = U + 1
c                        UEINAM( I ) = .TRUE.
c                        UEAREAD( I ) = .TRUE.
c                        EXIT
c                    END IF
c                END DO
c            END DO

          
            VARBUF = 'TYP_' // EINAM(I)( 1:L )
            L2 = LEN_TRIM( VARBUF )
            IF ( .NOT. READ3( UINAME, VARBUF( 1:L2 ), ALLAYS3, 0, 0,
     &                        EPTYP( 1, S ) ) ) THEN
                MESG = 'ERROR: Could not read ' // UINAME
     &                 // ' for EPTYP'
                CALL M3MSG2( MESG )
                EFLAG = .TRUE.
            END IF
      
            VARBUF = 'APR_' // EINAM(I)( 1:L )
            L2 = LEN_TRIM( VARBUF )
            IF ( .NOT. READ3( UINAME, VARBUF( 1:L2 ), ALLAYS3, 0, 0,
     &                        APRCH( 1, S ) ) ) THEN
                MESG = 'ERROR: Could not read ' // UINAME
     &                 // ' for APRCH'
                CALL M3MSG2( MESG )
                EFLAG = .TRUE.
            END IF

            VARBUF = 'NEP_' // EINAM(I)( 1:L )
            L2 = LEN_TRIM( VARBUF )
            IF ( .NOT. READ3( UINAME, VARBUF( 1:L2 ), ALLAYS3, 0, 0,
     &                        NUMEP( 1, S ) ) ) THEN
                MESG = 'ERROR: Could not read ' // UINAME
     &                 // ' for NUMEP'
                CALL M3MSG2( MESG )
                EFLAG = .TRUE.
            END IF
      
            VARBUF = 'UIX_' // EINAM(I)( 1:L )
            L2 = LEN_TRIM( VARBUF )
            IF ( .NOT. READ3( UINAME, VARBUF( 1:L2 ), ALLAYS3, 0, 0,
     &                        UNCIDX( 1, S ) ) ) THEN
                MESG = 'ERROR: Could not read ' // UINAME
     &                 // ' for UNCIDX'
                CALL M3MSG2( MESG )
                EFLAG = .TRUE.
            END IF
            
        END DO ! end loading METHOD, EPTYP, NUMEP, UNCIDX
       
C.........  Prompt for and open I/O API empirical file(s)...
        MESG = 'Reading uncertainty input empirical file...'
        NAMBUF = PROMPTMFILE( MESG, FSREAD3, 
     &                        UCAT( 1:1 ) // 'UCEOUT', PROGNAME )
        UENAME = NAMBUF
      
C.....  Get header description of the uncertainty empirical statistics file
        IF ( .NOT. DESC3( UENAME ) ) THEN
            MESG = 'Could not get description of file "' 
     &              // UENAME // '"'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF
      
C.........  Get uncertainty statistics information
C.........  Get empirical statistics information
        ALLOCATE( EMFVAL( NCOLS3D, NROWS3D ), STAT=IOS )
        CALL CHECKMEM( IOS, 'EMFVAL', PROGNAME )
        ALLOCATE( PROBVAL( NCOLS3D, NROWS3D ), STAT=IOS )
        CALL CHECKMEM( IOS, 'PROBVAL', PROGNAME )

        IF ( .NOT. READ3( UENAME, 'EMFVAL', ALLAYS3, 0, 0,
     &                    EMFVAL( 1, 1 ) ) ) THEN
            MESG = 'ERROR: Could not read ' // UENAME
     &             // ' for empirical data values'
            CALL M3MSG2( MESG )
            EFLAG = .TRUE.
        END IF
        
        IF ( .NOT. READ3( UENAME, 'PROBVAL', ALLAYS3, 0, 0,
     &                    PROBVAL( 1, 1 ) ) ) THEN
            MESG = 'ERROR: Could not read ' // UENAME
     &             // ' for empirical probability data'
            CALL M3MSG2( MESG )
            EFLAG = .TRUE.
        END IF

c bug not currently used
C.........  Prompt for and open I/O API parametric file(s)...
c        MESG = 'Reading uncertainty input parametric file...'
c        NAMBUF = PROMPTMFILE( MESG, FSREAD3, 
c     &                        UCAT( 1:1 ) // 'UCPOUT', PROGNAME )
c        UPNAME = NAMBUF
C.........  Get header description of the uncertainty empirical statistics file
c        IF( .NOT. DESC3( UPNAME ) ) THEN
c            MESG = 'Could not get description of file "' 
c     &              // UPNAME // '"'
c            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
c        END IF
C.........  Get paramertric statistics information
c        ALLOCATE( PARMS( NROWS3D, NCOLS3D ), STAT=IOS )
c        CALL CHECKMEM( IOS, 'PARMS', PROGNAME )
c        DO I = 1, NCOLS3D       
c            IF ( .NOT. READ3( UPNAME, 'PARMS', ALLAYS3, 0, 0,
c     &                        PARMS( 1, I ) ) ) THEN
c                MESG = 'ERROR: Could not read ' // UPNAME
c     &                 // ' for parametric data values'
c                CALL M3MSG2( MESG )
c                EFLAG = .TRUE.
c            END IF
c        END DO


c        CALL GETDTPOS(UPNAME, "PARMS", PARMS )

      
        RETURN

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

        END SUBROUTINE RDMETHS 
