
        INTEGER FUNCTION RDTPROF( FDEV, PROFTYPE, UFLAG )

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
C****************************************************************************

C.........  MODULES for public variables
C.........  For temporal profiles
        USE MODTMPRL, ONLY: HRLFAC, WEKFAC, MONFAC, MONFAC_ORG,
     &                      HRLREF, WEKREF, MONREF, METPROFLAG

        IMPLICIT NONE

C...........   INCLUDES:

        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters

C..........  EXTERNAL FUNCTIONS
        LOGICAL         BLKORCMT
        CHARACTER(2)    CRLF
        INTEGER         STR2INT         !  read (unjustified) int from str

        EXTERNAL        BLKORCMT, CRLF, STR2INT

C.........  SUBROUTINE ARGUMENTS
        INTEGER     , INTENT (IN) :: FDEV     ! unit number for profiles file
        CHARACTER(*), INTENT (IN) :: PROFTYPE ! 'MONTHLY', 'WEEKLY', etc.
        LOGICAL     , INTENT (IN) :: UFLAG    ! true: use uniform profiles

C.........  Local parameters
        INTEGER, PARAMETER :: NDTYPE = 9      ! no. diurnal profile types

        CHARACTER(9), PARAMETER :: DIURTYPE( NDTYPE ) = 
     &                   ( / 'WEEKDAY  ', 'WEEKEND  ', 'MONDAY   ',
     &                       'TUESDAY  ', 'WEDNESDAY', 'THURSDAY ', 
     &                       'FRIDAY   ', 'SATURDAY ', 'SUNDAY   '   / )

C.........  Unsorted temporal profiles
        INTEGER, ALLOCATABLE :: CODEA( : )    !  list of profile codes
        INTEGER, ALLOCATABLE :: INDXA( : )    !  list of subscripts
        INTEGER, ALLOCATABLE :: WT   ( : )    !  weight factors

        REAL   , ALLOCATABLE :: PFACA( :,: )  !  returned table of factors

C.........  Local, sorted temporal profile codes
        INTEGER, ALLOCATABLE :: TMPREF( : ) 

C...........   Local arrays
        INTEGER         DNPROF( NDTYPE )  ! no. diurnal profiles of each type
        INTEGER         DNSKIP( NDTYPE )  ! no. lines to skip in input file

C...........   SCRATCH LOCAL VARIABLES and their descriptions:

        INTEGER         I, J, K, L, N    !  counters and indices

        INTEGER         IOS           ! i/o status
        INTEGER         NFAC          ! tmp number of factors per profile
        INTEGER         NPROF         ! actual number of profiles
        INTEGER         NSKIP         ! number of lines before desired profiles
        INTEGER         W             ! temporary profile entry

        LOGICAL      :: DFLAG  = .FALSE.   ! diurnal profiles initialized
        LOGICAL      :: EFLAG  = .FALSE.  !  input error flag

        CHARACTER(300)  MESG    !  message buffer

        CHARACTER(16) :: PROGNAME = 'RDTPROF' ! program name

C***********************************************************************
C   begin body of subroutine  RDTPROF

C.........  Initialize counts for each call of  function
        NPROF = 0
        NSKIP = 0

C.........  If using uniform profiles...
        IF( UFLAG ) THEN

            NPROF = 1

C.........  Otherwise, if not overriding with uniform profiles...
        ELSE

C.............  Determine the number of entries for the requested profile type
C.............  For monthly and weekly, send the requested type name directly...
            IF( PROFTYPE .EQ. 'MONTHLY' .OR.
     &          PROFTYPE .EQ. 'WEEKLY'       ) THEN

                CALL COUNT_TPROF( PROFTYPE, NSKIP, NPROF )

C.............  For diurnal, check for all profile types, and ensure that 
C               each type has the same number of profiles
            ELSE

C.................  Initialize count of diurnal profiles of each type
                DNPROF = 0   !  array

C.................  Search for the nine different possible profiles
                DO I = 1, 9

                    CALL COUNT_TPROF( DIURTYPE( I ), 
     &                                DNSKIP( I ), DNPROF( I ) )

C.....................  When profiles of current type are found...
                    IF( DNPROF( I ) .GT. 0 ) THEN

                        L = LEN_TRIM( DIURTYPE( I ) )

C.........................  Ensure that this total is consistent with all other
C                           diurnal temporal profiles
                        IF( NPROF .GT. 0          .AND. 
     &                      NPROF .NE. DNPROF( I )      ) THEN

                            EFLAG = .TRUE.
                            WRITE( MESG,94010 ) 'ERROR: Number of' //
     &                         'diurnal temporal profiles in previous'//
     &                         'packet'// CRLF() // BLANK10 // 'was',
     &                         NPROF, 'but "' // DIURTYPE( I )( 1:L ) //
     &                         '" packet has', DNPROF( I )
                            CALL M3MSG2( MESG )

C.........................  Otherwise, initialize number of diurnal profiles
                        ELSE
                            NPROF = DNPROF( I )

                        END IF

C.........................  Message about profiles being used
                        IF( DNPROF( I ) .GT. 0 ) THEN

                            MESG = 'NOTE: Storing '// DIURTYPE( I )(1:L)
     &                             // ' temporal profiles'
                            CALL M3MSG2( MESG )

                        END IF

                    END IF

                END DO

C.................  Ensure that appropriate diurnal temporal profiles are 
C                   available
C.................  If weekday diurnal is not provided and monday-friday is not
C                   either, then error.
                IF( DNPROF( 1 ) .EQ. 0 .AND. 
     &            ( DNPROF( 3 ) .EQ. 0 .OR. 
     &              DNPROF( 4 ) .EQ. 0 .OR. 
     &              DNPROF( 5 ) .EQ. 0 .OR. 
     &              DNPROF( 6 ) .EQ. 0 .OR. 
     &              DNPROF( 7 ) .EQ. 0      ) ) THEN

                    EFLAG = .TRUE.
                    MESG = 'ERROR: Weekday diurnal profiles are ' //
     &                     'needed when any of Monday through'// 
     &                     CRLF()// BLANK10// 'Friday profiles are ' //
     &                     'not provided.'
                    CALL M3MSG2( MESG )

                END IF

C.................  If weekend diurnal and weekday diurnal are not provided, 
C                   and saturday or sunday are not available, then error
                IF( DNPROF( 1 ) .EQ. 0 .AND. 
     &              DNPROF( 2 ) .EQ. 0 .AND.
     &            ( DNPROF( 8 ) .EQ. 0 .OR.
     &              DNPROF( 9 ) .EQ. 0      ) ) THEN

                    EFLAG = .TRUE.
                    MESG = 'ERROR: Weekday or weekend diurnal ' //
     &                     'profiles are needed when either'// CRLF() //
     &                     BLANK10// 'Saturday or Sunday profiles ' //
     &                     'are not provided.'
                    CALL M3MSG2( MESG )

                END IF

            END IF     ! End if Monthly/Weekly or Diurnal

            IF( EFLAG ) THEN
                MESG = 'Problem reading temporal profiles file.'
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            ENDIF

        END IF   !  End if uniform profiles or not

C.........  Allocate memory for sorted arrays, depending on profile type. 
C.........  Initialize factors to 1.0
        SELECT CASE( PROFTYPE ) 

        CASE( 'DIURNAL' )
            NFAC = 24
            ALLOCATE( HRLREF( NPROF ), STAT=IOS )
            CALL CHECKMEM( IOS, 'HRLREF', PROGNAME )
            ALLOCATE( TMPREF( NPROF ), STAT=IOS )
            CALL CHECKMEM( IOS, 'TMPREF', PROGNAME )
            ALLOCATE( HRLFAC( NFAC, NPROF, 7 ), STAT=IOS )
            CALL CHECKMEM( IOS, 'HRLFAC', PROGNAME )

            HRLFAC = 1.0

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
            ALLOCATE( MONFAC_ORG( NFAC,NPROF ), STAT=IOS )
            CALL CHECKMEM( IOS, 'MONFAC_ORG', PROGNAME )
            MONFAC_ORG = 1.0
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

            MESG = 'WARNING: No temporal profiles of type "' // 
     &             PROFTYPE // '" were found.'
            CALL M3MSG2( MESG )

            RDTPROF = NPROF
            RETURN

        END IF

C.........  Allocate memory for unsorted arrays
        ALLOCATE( INDXA( NPROF ), STAT=IOS )
        CALL CHECKMEM( IOS, 'INDXA', PROGNAME )
        ALLOCATE( CODEA( NPROF ), STAT=IOS )
        CALL CHECKMEM( IOS, 'CODEA', PROGNAME )
        ALLOCATE( WT( NFAC + 1 ), STAT=IOS )
        CALL CHECKMEM( IOS, 'WT', PROGNAME )
        ALLOCATE( PFACA( NFAC + 1, NPROF ), STAT=IOS )
        CALL CHECKMEM( IOS, 'PFACA', PROGNAME )

C.........  If using uniform profiles, set these
        IF( UFLAG ) THEN

            CODEA( 1 ) = 1
            INDXA( 1 ) = 1
            PFACA( 1:NFAC,1 ) = 1.
            PFACA( NFAC+1,1 ) = REAL( NFAC )

C.............  Store as sorted profiles, depending on profile type
            SELECT CASE( PROFTYPE )

            CASE( 'DIURNAL' )
                CALL STORE_TPROF( NPROF, NFAC, HRLREF, HRLFAC( 1,1,1 ) )
                CALL STORE_TPROF( NPROF, NFAC, HRLREF, HRLFAC( 1,1,2 ) )
                CALL STORE_TPROF( NPROF, NFAC, HRLREF, HRLFAC( 1,1,3 ) )
                CALL STORE_TPROF( NPROF, NFAC, HRLREF, HRLFAC( 1,1,4 ) )
                CALL STORE_TPROF( NPROF, NFAC, HRLREF, HRLFAC( 1,1,5 ) )
                CALL STORE_TPROF( NPROF, NFAC, HRLREF, HRLFAC( 1,1,6 ) )
                CALL STORE_TPROF( NPROF, NFAC, HRLREF, HRLFAC( 1,1,7 ) )

            CASE( 'WEEKLY' )
                CALL STORE_TPROF( NPROF, NFAC, WEKREF, WEKFAC )

            CASE( 'MONTHLY' )
                CALL STORE_TPROF( NPROF, NFAC, MONREF, MONFAC )

            END SELECT

C.........  If not using uniform profiles, read and store profiles
        ELSE

C.............  Depending on profile type....
            SELECT CASE( PROFTYPE )

            CASE( 'DIURNAL' )

C.................  Weekday diurnal
                IF( DNPROF( 1 ) .GT. 0 ) THEN

                    CALL READ_TPROF ( DNSKIP( 1 ), NPROF )
                    CALL STORE_TPROF( NPROF,NFAC,HRLREF,HRLFAC( 1,1,1 ))
                    CALL STORE_TPROF( NPROF,NFAC,HRLREF,HRLFAC( 1,1,2 ))
                    CALL STORE_TPROF( NPROF,NFAC,HRLREF,HRLFAC( 1,1,3 ))
                    CALL STORE_TPROF( NPROF,NFAC,HRLREF,HRLFAC( 1,1,4 ))
                    CALL STORE_TPROF( NPROF,NFAC,HRLREF,HRLFAC( 1,1,5 ))
                    CALL STORE_TPROF( NPROF,NFAC,HRLREF,HRLFAC( 1,1,6 ))
                    CALL STORE_TPROF( NPROF,NFAC,HRLREF,HRLFAC( 1,1,7 ))

                    DFLAG = .TRUE.

                END IF

C.................  Weekend diurnal. Read, store, and make sure consistent with
C                   previously read profile codes.
                IF( DNPROF( 2 ) .GT. 0 ) THEN

                    CALL READ_TPROF ( DNSKIP( 2 ), NPROF )
                    CALL STORE_TPROF( NPROF,NFAC,TMPREF,HRLFAC( 1,1,6 ))
                    CALL STORE_TPROF( NPROF,NFAC,TMPREF,HRLFAC( 1,1,7 ))

                    IF ( DFLAG ) THEN
                        CALL CHECK_PROFCODE( DIURTYPE(2), NPROF, 
     &                                       HRLREF, TMPREF      )
                    ELSE
                        HRLREF = TMPREF   ! array
                    END IF

                    DFLAG = .TRUE.

                END IF

C.................  Each day of the week diurnal
                DO I = 3, 9

                    IF( DNPROF( I ) .GT. 0 ) THEN

                        CALL READ_TPROF ( DNSKIP( I ), NPROF )
                        CALL STORE_TPROF( NPROF, NFAC, TMPREF, 
     &                                    HRLFAC( 1,1,I-2 )    )

                        IF ( DFLAG ) THEN
                            CALL CHECK_PROFCODE( DIURTYPE( I ), NPROF,
     &                                           HRLREF, TMPREF )
                        ELSE
                            HRLREF = TMPREF   ! array
                        END IF

                        DFLAG = .TRUE.

                    END IF

                END DO

C............  Weekly profiles
            CASE( 'WEEKLY' )
                CALL READ_TPROF ( NSKIP, NPROF ) 
                CALL STORE_TPROF( NPROF, NFAC, WEKREF, WEKFAC )

C............  Monthly profiles
            CASE( 'MONTHLY' )

                CALL READ_TPROF ( NSKIP, NPROF ) 
                CALL STORE_TPROF( NPROF, NFAC, MONREF, MONFAC )

            END SELECT

            IF( EFLAG ) THEN
                MESG = 'Problem reading temporal profiles file.'
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF

        END IF  ! End of uniform profiles or not


C.........  Deallocate unsorted arrays
        DEALLOCATE( INDXA, CODEA, WT, PFACA )

C.........  Deallocate other local memory, if needed
        IF( ALLOCATED( TMPREF ) ) DEALLOCATE( TMPREF )

        RDTPROF = NPROF

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10 ( A, :, I8, :, 1X ) )

C******************  INTERNAL SUBPROGRAMS  *****************************

        CONTAINS

C.............  This internal subprogram counts the number of temporal
C               profiles of a certain type

            SUBROUTINE COUNT_TPROF( PTYPE, NSKIP, NPROF )

C.............  Subroutine arguments
            CHARACTER(*), INTENT (IN) :: PTYPE    ! profile type
            INTEGER     , INTENT(OUT) :: NSKIP    ! number of lines to skip
            INTEGER     , INTENT(OUT) :: NPROF    ! number of profiles

C.............  Local variables
            INTEGER  I, J, K, L

            INTEGER         IREC            ! record counter
            INTEGER         IOS             ! i/o status

            LOGICAL      :: FOUND           !  true: profile type has been found
            
            CHARACTER(120)  LINE    !  char-string buffer for profiles input
            CHARACTER(300)  MESG    !  message buffer

C----------------------------------------------------------------------

C.............  Initialize output values
            NSKIP = 0
            NPROF = 0

C.............  Ensure that profile type name has not yet been found
            FOUND = .FALSE.

            I    = 0
            IREC = 0
            L    = LEN_TRIM( PTYPE )
C.............  Head of counting read loop
            DO

C.................  Read line of profile.  Use END in case the profile type
C                   requested is not in the file - such as weekend
                READ( FDEV, 93000, END=111, IOSTAT=IOS ) LINE
                IREC = IREC + 1

C.................  Check read error status
                IF( IOS .GT. 0 ) THEN
                    EFLAG = .TRUE.
                    WRITE( MESG,94010 ) 'I/O error', IOS, 
     &                     'reading TEMPORAL PROFILE file at line', IREC
 
                    CALL M3MESG( MESG )
                    CYCLE
                END IF
                
                IF( BLKORCMT( LINE ) ) CYCLE

C.................  Scan line for profile type (e.g., /MONTHLY/)
                IF( .NOT. FOUND ) THEN
                    J = INDEX( LINE, PTYPE( 1:L ) )
                    IF( J .GT. 0 ) FOUND = .TRUE.
                    NSKIP = IREC

C.................  Count records of input profile type, and look of end of 
C                   section
                ELSE

                    K = INDEX( LINE, '/END/' )

                    IF( K .GT. 0 ) EXIT

                    I = I + 1

                END IF
          
            END DO   ! End of first read though to determine memory needs

111         NPROF = I


            REWIND( FDEV )        
       
C-------------------    FORMAT  STATEMENTS   ---------------------------

C...........   Formatted file I/O formats............ 93xxx

93000   FORMAT( A )

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10 ( A, :, I8, :, 1X ) )

            END SUBROUTINE COUNT_TPROF

C----------------------------------------------------------------------
C----------------------------------------------------------------------

C.............  This internal subprogram reads the temporal profiles of
C               a given type, and stores in unsorted order

            SUBROUTINE READ_TPROF( NSKIP, NPROF )

C.............  Subroutine arguments
            INTEGER, INTENT (IN) :: NSKIP   ! number of lines before profiles
            INTEGER, INTENT (IN) :: NPROF   ! number of profiles expected

C.............  Local variables
            INTEGER    I, J, K

            INTEGER         IREC            ! record counter
            INTEGER         IOS             ! i/o status
            INTEGER         ISUM            ! input row-total

            REAL            DIV             !  scratch divisor

            CHARACTER(120)  LINE    !  char-string buffer for profiles input
            CHARACTER(300)  MESG    !  message buffer

C----------------------------------------------------------------------
        
C.............  Initialize profiles and factors
            CODEA = 0   ! array
            INDXA = 0   ! array
            PFACA = 0.  ! array

C.............  Skip irrelevant lines in input profiles file
            CALL SKIPL( FDEV, NSKIP )

C.............  Read unsorted entries of the requested profile type
            IREC = NSKIP
            DO I = 1, NPROF

C.................  Loop until non-blank or comment line
                DO
                    READ( FDEV, 93000, IOSTAT=IOS ) LINE
                    IREC = IREC + 1

                    IF( IOS .GT. 0 ) THEN
                        EFLAG = .TRUE.
                        WRITE( MESG,94010 ) 'I/O error', IOS, 
     &                    'reading TEMPORAL PROFILE file at line', IREC
                        CALL M3MESG( MESG )
                        CYCLE
                    END IF
                
                    IF( .NOT. BLKORCMT( LINE ) ) EXIT
                END DO 

                CODEA( I ) = STR2INT( LINE( 1:5 ) )
                INDXA( I ) = I

C.................  Check for bad cross-reference code
                IF( .NOT. METPROFLAG .AND. CODEA( I ) == 99999 ) THEN
                    WRITE( MESG, 94010 )
     &                  'ERROR: CAN NOT USE temporal profile code ',
     &                  CODEA(I), ' at line ', IREC,
     &                  'while processing Met-based profiles' // 
     &                  CRLF() // BLANK16 // 'MUST be less than 99999'
                    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                END IF

C.................  Convert columns from ASCII to integer in groups of 4
                J    = 6
                K    = 9         !  j:k spans 4 characters
                ISUM = 0
                DO N = 1, NFAC
                    WT( N ) = STR2INT( LINE( J:K ) )
                    J       = J + 4
                    K       = K + 4
                END DO

C.................  Final field is 1-character wider than others
                WT( NFAC + 1 ) = STR2INT( LINE( J:K+1 ) )

                IF ( WT( NFAC+1 ) .NE. 0 ) THEN
                    DIV = 1.0 / FLOAT( WT( NFAC+1 ) )
                ELSE
                    DIV = 0.0
                END IF

                DO N = 1, NFAC
                    PFACA( N,I ) = DIV * FLOAT( WT( N ) )
                END DO
            END DO

            REWIND( FDEV )

C-------------------    FORMAT  STATEMENTS   ---------------------------

C...........   Formatted file I/O formats............ 93xxx

93000   FORMAT( A )

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10 ( A, :, I8, :, 1X ) )

            END SUBROUTINE READ_TPROF

C----------------------------------------------------------------------
C----------------------------------------------------------------------

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

C.............  Sort requested profile type
            CALL SORTI1( N, INDXA, CODEA )

C.............  Store in sorted order
            DO I = 1, N

                K = INDXA( I )
                CODES( I ) = CODEA( K )

                DO J = 1, M

                    PFACS( J,I ) = PFACA( J,K )

                END DO

            END DO
        
            END SUBROUTINE STORE_TPROF

C----------------------------------------------------------------------
C----------------------------------------------------------------------

C.............  This internal subprogram unsures that the temporal profile
C               codes are consistent
            SUBROUTINE CHECK_PROFCODE( PTYPE, N, REFCODES, CODES )

C.............  Subroutine arguments
            CHARACTER(*), INTENT (IN) :: PTYPE         ! profile type
            INTEGER     , INTENT (IN) :: N             ! number of profiles
            INTEGER     , INTENT (IN) :: REFCODES( N ) ! reference profile codes
            INTEGER     , INTENT (IN) :: CODES   ( N ) ! profile codes to check

C.............  Local variables
            INTEGER  I

            CHARACTER(300) MESG

C----------------------------------------------------------------------

            DO I = 1, N

                IF( REFCODES( I ) .NE. CODES( I ) ) THEN

                    EFLAG = .TRUE.
                    MESG = 'ERROR: Temporal profile codes for ' //
     &                     PTYPE // 'profiles are not consistent' //
     &                     CRLF() // BLANK10 // 'with initialized ' //
     &                     'values.'
                    CALL M3MSG2( MESG )

                    EXIT   ! Exit from loop

                END IF

            END DO
        
            END SUBROUTINE CHECK_PROFCODE

        END FUNCTION RDTPROF

