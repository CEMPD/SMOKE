
        INTEGER FUNCTION RDTPROF( FDEV, PROFTYPE )

C***********************************************************************
C  function body starts at line 
C
C  DESCRIPTION:
C      Reads the input temporal profiles and computes fractions based
C      on the input total without renormalization
C
C  PRECONDITIONS REQUIRED:
C       this segment of the input file sorted
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C       TRIMLEN, M3ERR
C
C  REVISION  HISTORY:
C       Copied from RDTPROF.F version 1.3 by M Houyoux 1/99
C
C***********************************************************************
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
C****************************************************************************

C.........  MODULES for public variables
C.........  For temporal profiles
        USE MODTPRO

        IMPLICIT NONE

C..........  EXTERNAL FUNCTIONS
        INTEGER         STR2INT         !  read (unjustified) int from str

        EXTERNAL        STR2INT

C.........  SUBROUTINE ARGUMENTS
        INTEGER     , INTENT (IN) :: FDEV     ! unit number for profiles file
        CHARACTER(*), INTENT (IN) :: PROFTYPE ! 'MONTHLY', 'WEEKLY', etc.

C.........  Unsorted temporal profiles

        INTEGER, ALLOCATABLE :: CODEA( : )    !  list of profile codes
        INTEGER, ALLOCATABLE :: INDXA( : )    !  list of subscripts
        INTEGER, ALLOCATABLE :: WT   ( : )    !  weight factors

        REAL   , ALLOCATABLE :: PFACA( :,: )  !  returned table of factors

C...........   SCRATCH LOCAL VARIABLES and their descriptions:

        INTEGER         I, J, K, N    !  counters and indices
        INTEGER         IOS           ! i/o status
        INTEGER         IREC          ! line counter
        INTEGER         ISUM          ! input row-total
        INTEGER         NFAC          ! tmp number of factors per profile
        INTEGER         NPROF         ! actual number of profiles
        INTEGER         NSKIP         ! number of lines before desired profiles
        INTEGER         W             ! temporary profile entry

        REAL            DIV             !  scratch divisor

        LOGICAL      :: EFLAG  = .FALSE.  !  input error flag
        LOGICAL      :: FOUND             !  true if profile type has been found

        CHARACTER*120   LINE    !  character-string buffer for profiles input
        CHARACTER*300   MESG    !  message buffer

        CHARACTER*16 :: PROGNAME = 'RDTPROF' ! program name

C***********************************************************************
C   begin body of subroutine  RDTPROF

C.........  Ensure that profile type name has not yet been found
        FOUND = .FALSE.

C.........  Determine the number of entries for the requested profile type
        I    = 0
        IREC = 0
        DO     ! head of first read loop

C.............  Read line of profile.  Use END in case the profile type
C               requested is not in the file - such as weekend
            READ( FDEV, 93000, END=111, IOSTAT=IOS ) LINE
            IREC = IREC + 1

C.............  Check read error status
            IF( IOS .GT. 0 ) THEN
                EFLAG = .TRUE.
                WRITE( MESG,94010 ) 'I/O error', IOS, 
     &                 'reading TEMPORAL PROFILE file at line', IREC

                CALL M3MESG( MESG )
                CYCLE
            ENDIF

C.............  Scan line for profile type (e.g., /MONTHLY/)
            IF( .NOT. FOUND ) THEN
                J = INDEX( LINE, PROFTYPE )
                IF( J .GT. 0 ) FOUND = .TRUE.
                NSKIP = IREC

C.............  Count records of input profile type, and look of end of section
            ELSE

                K = INDEX( LINE, '/END/' )

                IF( K .GT. 0 ) EXIT

                I = I + 1

            ENDIF

        ENDDO   ! End of first read though to determine memory needs

111     NPROF = I

        REWIND( FDEV )

        IF( EFLAG ) THEN
            MESG = 'Problem reading temporal profiles file.'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        ENDIF

C.........  Allocate memory for sorted arrays, depending on profile type. 
C.........  Initialize factors to 1.0
        SELECT CASE( PROFTYPE ) 

        CASE( 'WEEKDAY' )
            NFAC = 24
            ALLOCATE( DIUREF( NPROF ), STAT=IOS )
            CALL CHECKMEM( IOS, 'DIUREF', PROGNAME )
            ALLOCATE( DIUFAC( NFAC,NPROF ), STAT=IOS )
            CALL CHECKMEM( IOS, 'DIUFAC', PROGNAME )

            DIUFAC = 1.0

        CASE( 'WEEKEND' )
            NFAC = 24
            ALLOCATE( ENDREF( NPROF ), STAT=IOS )
            CALL CHECKMEM( IOS, 'ENDREF', PROGNAME )
            ALLOCATE( ENDFAC( NFAC,NPROF ), STAT=IOS )
            CALL CHECKMEM( IOS, 'ENDFAC', PROGNAME )

            ENDFAC = 1.0

        CASE( 'WEEKLY' )
            NFAC = 7
            ALLOCATE( WEKREF( NPROF ), STAT=IOS )
            CALL CHECKMEM( IOS, 'WEKREF', PROGNAME )
            ALLOCATE( WEKFAC( NFAC,NPROF ), STAT=IOS )
            CALL CHECKMEM( IOS, 'WEKFAC', PROGNAME )

            WEKFAC = 1.0

        CASE( 'MONTHLY' ) 
            NFAC = 12
            ALLOCATE( MONREF( NPROF ), STAT=IOS )
            CALL CHECKMEM( IOS, 'MONREF', PROGNAME )
            ALLOCATE( MONFAC( NFAC,NPROF ), STAT=IOS )
            CALL CHECKMEM( IOS, 'MONFAC', PROGNAME )

            MONFAC = 1.0

        CASE DEFAULT

            MESG = 'INTERNAL ERROR: Profile type "' // PROFTYPE // 
     &             '" not known in program ' // 
     &             PROGNAME( 1:LEN_TRIM( PROGNAME ) )
            CALL M3MSG2( MESG )
            CALL M3EXIT( PROGNAME, 0, 0, ' ', 2 )

        END SELECT

C.........  Check for profiles. Note that the arrays should have still been 
C           allocated with zero dimension.
        IF( NPROF .EQ. 0 ) THEN

            MESG = 'No temporal profiles of type "' // PROFTYPE //
     &             '" were found.'
            CALL M3WARN( PROGNAME, 0, 0, MESG )

            RDTPROF = NPROF
            RETURN

        ENDIF


C.........  Allocate memory for unsorted arrays
        ALLOCATE( INDXA( NPROF ), STAT=IOS )
        CALL CHECKMEM( IOS, 'INDXA', PROGNAME )
        ALLOCATE( CODEA( NPROF ), STAT=IOS )
        CALL CHECKMEM( IOS, 'CODEA', PROGNAME )
        ALLOCATE( WT( NFAC + 1 ), STAT=IOS )
        CALL CHECKMEM( IOS, 'WT', PROGNAME )
        ALLOCATE( PFACA( NFAC + 1, NPROF ), STAT=IOS )
        CALL CHECKMEM( IOS, 'PFACA', PROGNAME )

C.........  Skip irrelevant lines in input profiles file
        CALL SKIPL( FDEV, NSKIP )

C.........  Read unsorted entries of the requested profile type
        IREC = NSKIP
        DO I = 1, NPROF

            READ( FDEV, 93000, IOSTAT=IOS ) LINE
            IREC = IREC + 1

            IF( IOS .GT. 0 ) THEN
                EFLAG = .TRUE.
                WRITE( MESG,94010 ) 'I/O error', IOS, 
     &                 'reading TEMPORAL PROFILE file at line', IREC

                CALL M3MESG( MESG )
                CYCLE
            ENDIF

            CODEA( I ) = STR2INT( LINE( 1:5 ) )
            INDXA( I ) = I

C.............  Convert columns from ASCII to integer in groups of 4
            J    = 6
            K    = 9	!  j:k spans 4 characters
            ISUM = 0
            DO N = 1, NFAC
                WT( N ) = STR2INT( LINE( J:K ) )
                J       = J + 4
                K       = K + 4
            ENDDO

C.............  Final field is 1-character wider than others
            WT( NFAC + 1 ) = STR2INT( LINE( J:K+1 ) )

            IF ( WT( NFAC+1 ) .NE. 0 ) THEN
                DIV = 1.0 / FLOAT( WT( NFAC+1 ) )
            ELSE
                DIV = 0.0
            END IF

            DO N = 1, NFAC
                PFACA( N,I ) = DIV * FLOAT( WT( N ) )
            ENDDO
            
        ENDDO

        REWIND( FDEV )

        IF( EFLAG ) THEN
            MESG = 'Problem reading temporal profiles file.'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        ENDIF

C.........  Sort requested profile type
        CALL SORTI1( NPROF, INDXA, CODEA )

C.........  Store sorted profiles, depending on profile type
        SELECT CASE( PROFTYPE )

        CASE( 'WEEKDAY' )
            CALL STORE_TPROF( NPROF, NFAC, DIUREF, DIUFAC )

        CASE( 'WEEKEND' )
            CALL STORE_TPROF( NPROF, NFAC, ENDREF, ENDFAC )

        CASE( 'WEEKLY' )
            CALL STORE_TPROF( NPROF, NFAC, WEKREF, WEKFAC )

        CASE( 'MONTHLY' )
            CALL STORE_TPROF( NPROF, NFAC, MONREF, MONFAC )

        END SELECT

C.........  Deallocate unsorted arrays
        DEALLOCATE( INDXA, CODEA, WT, PFACA )

        RDTPROF = NPROF

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Formatted file I/O formats............ 93xxx

93000   FORMAT( A )

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10 ( A, :, I8, :, 1X ) )

C******************  INTERNAL SUBPROGRAMS  *****************************

        CONTAINS

C.............  This internal subprogram writes a warning message for 
C               duplicate entries in the cross-reference file.
            SUBROUTINE STORE_TPROF( N, M, CODES, PFACS )

C.............  Subroutine arguments
            INTEGER, INTENT (IN) :: N             ! number of profiles
            INTEGER, INTENT (IN) :: M             ! number of facs per profile
            INTEGER, INTENT(OUT) :: CODES( N )    ! profile codes
            REAL   , INTENT(OUT) :: PFACS( M, N ) ! profile factors

C.............  Local variables
            INTEGER  I, J, K

C----------------------------------------------------------------------

            DO I = 1, N

                K = INDXA( I )
                CODES( I ) = CODEA( K )

                DO J = 1, M

                    PFACS( J,I ) = PFACA( J,K )

                ENDDO

            ENDDO
        
            END SUBROUTINE STORE_TPROF

        END FUNCTION RDTPROF

