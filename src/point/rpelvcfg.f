
        SUBROUTINE RPELVCFG( FDEV )

C***********************************************************************
C  subroutine body starts at line 
C
C  DESCRIPTION:
C      Subroutine RPELVCFG reads the PELVCONFIG file that is used to set
C      the elevated groups, the elevated sources, and the plume-in-grid
C      (PinG) sources.
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C     Created 7/2001 by M Houyoux
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
C***********************************************************************

C...........   MODULES for public variables
C.........  This module contains Smkreport-specific settings
        USE MODREPRT, ONLY: PKT_IDX, ELG_IDX, PNG_IDX, ELV_IDX,
     &                      SPCF_NAND, INSPCIFY, LIN_SPCIFY,
     &                      SPCF_NOR, NCRTSYBL, CRTSYBL

C.........  This module contains arrays for plume-in-grid and major sources
        USE MODELEV, ONLY: GRPVALS, GRPTYPES, 
     &                     ELVVALS, ELVTYPES, ELVCHRS,
     &                     PNGVALS, PNGTYPES, PNGCHRS,
     &                     EVPEMIDX, EVPESTAT, EVPPSTAT,
     &                     NGRPVAR, NEVPVAR, NEVPEMV, 
     &                     NGRPCRIT, NELVCRIT, NPNGCRIT,
     &                     MXGRPCHK, MXELVCHK, MXPNGCHK,
     &                     LCUTOFF, LPNGRNK, LELVRNK

C.........  This module contains the information about the source category
        USE MODINFO, ONLY: EINAM, NIPOL

        IMPLICIT NONE

C...........   INCLUDES
        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters

C...........   EXTERNAL FUNCTIONS and their descriptions:
        CHARACTER(2)  CRLF
        INTEGER       GETFLINE

        EXTERNAL   CRLF, GETFLINE

C...........   SUBROUTINE ARGUMENTS
        INTEGER     , INTENT (IN) :: FDEV       ! input file unit

C...........   Local allocatable arrays
        INTEGER, ALLOCATABLE :: TMPIDX( : )  ! tmp index for flagging pollutants

        LOGICAL, ALLOCATABLE :: EISTAT( : )  ! true: pollutant used as criteria

C...........   Other local variables
        INTEGER         I, V   ! indices and counters

        INTEGER         IOS                  ! i/o status
        INTEGER,SAVE :: NALLPOL = 0          ! number of pols used as criteria
        INTEGER         NLINES               ! number of lines

        LOGICAL,SAVE :: EFLAG = .FALSE.      ! true: error found

        CHARACTER(300)  MESG                 ! message buffer

        CHARACTER(16) :: PROGNAME = 'RPELVCFG' ! program name

C***********************************************************************
C   begin body of subroutine RPELVCFG

C.........  Write status message
        MESG = 'Reading elevated source configuration file...'
        CALL M3MSG2( MESG )

C.........  Initialize the number of variables for grouping, selecting PinG
C           sources, and selecting elevated sources. 
        NGRPVAR = 5     ! The number of stack parmeters
        NEVPVAR = MAX( HT_IDX, DM_IDX, TK_IDX, VE_IDX, FL_IDX,
     &                 SRC_IDX, FIP_IDX, PLT_IDX )

C.........  Get number of lines in file
        NLINES = GETFLINE( FDEV, 'Elevated configuration' )

C.........  Allocate for local temporary indices
        ALLOCATE( TMPIDX( NIPOL ), STAT=IOS )
        CALL CHECKMEM( IOS, 'TMPIDX', PROGNAME )
        ALLOCATE( EISTAT( NIPOL ), STAT=IOS )
        CALL CHECKMEM( IOS, 'EISTAT', PROGNAME )
        TMPIDX = 0
        EISTAT = .FALSE.

C.........  Read through file to determine how many ORs and ANDs for each
C           section of file
        CALL READ_PELVCONFIG( FDEV, 'COUNT', EFLAG )

C.........  Abort if error
        IF ( EFLAG ) THEN
            MESG = 'Problem scanning elevated configuration file'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

        NEVPVAR = NEVPVAR + NALLPOL

        NEVPEMV = 0
        DO V = 1, NIPOL
            IF ( EISTAT( V ) ) NEVPEMV = NEVPEMV + 1
        END DO  

C.........  Limit group criteria to 1 for later usage in ASGNGRPS
        NGRPCRIT = MAX( NGRPCRIT, 1 )
        MXGRPCHK = MAX( MXGRPCHK, 1 )

C.........  Allocate memory for criteria arrays
        ALLOCATE( GRPVALS( NGRPCRIT, MXGRPCHK, NGRPVAR ), STAT=IOS )
        CALL CHECKMEM( IOS, 'GRPVALS', PROGNAME )
        ALLOCATE( GRPTYPES( NGRPCRIT, MXGRPCHK, NGRPVAR ), STAT=IOS )
        CALL CHECKMEM( IOS, 'GRPTYPES', PROGNAME )
        GRPVALS  = 0.
        GRPTYPES = '<' ! bogus comparison type. By default, no groups get set.

        ALLOCATE( ELVVALS( NELVCRIT, MXELVCHK, NEVPVAR ), STAT=IOS )
        CALL CHECKMEM( IOS, 'ELVVALS', PROGNAME )
        ALLOCATE( ELVCHRS( NELVCRIT, MXELVCHK, NEVPVAR ), STAT=IOS )
        CALL CHECKMEM( IOS, 'ELVCHRS', PROGNAME )
        ALLOCATE( ELVTYPES( NELVCRIT, MXELVCHK, NEVPVAR ), STAT=IOS )
        CALL CHECKMEM( IOS, 'ELVTYPES', PROGNAME )
        ELVVALS  = 0.
        ELVCHRS  = ' '
        ELVTYPES = ' '

        ALLOCATE( PNGVALS( NPNGCRIT, MXPNGCHK, NEVPVAR ), STAT=IOS )
        CALL CHECKMEM( IOS, 'PNGVALS', PROGNAME )
        ALLOCATE( PNGCHRS( NPNGCRIT, MXPNGCHK, NEVPVAR ), STAT=IOS )
        CALL CHECKMEM( IOS, 'PNGCHRS', PROGNAME )
        ALLOCATE( PNGTYPES( NPNGCRIT, MXPNGCHK, NEVPVAR ), STAT=IOS )
        CALL CHECKMEM( IOS, 'PNGTYPES', PROGNAME )
        PNGVALS  = 0.
        PNGCHRS  = ' '
        PNGTYPES = ' '

        ALLOCATE( EVPEMIDX( NEVPEMV ), STAT=IOS )
        CALL CHECKMEM( IOS, 'EVPEMIDX', PROGNAME )
        ALLOCATE( EVPESTAT( NIPOL ), STAT=IOS )
        CALL CHECKMEM( IOS, 'EVPESTAT', PROGNAME )
        ALLOCATE( EVPPSTAT( NIPOL ), STAT=IOS )
        CALL CHECKMEM( IOS, 'EVPPSTAT', PROGNAME )
        EVPEMIDX = 0
        EVPESTAT = .FALSE.
        EVPPSTAT = .FALSE.

C.........  Store index from pollutants used as selection criteria to master
        I = 0
        DO V = 1, NIPOL
            IF ( EISTAT( V ) ) THEN
                I = I + 1
                EVPEMIDX( I ) = V
                TMPIDX  ( V ) = I
            END IF
        END DO  

C.........  Read through file to store criteria arrays
        CALL READ_PELVCONFIG( FDEV, 'STORE', EFLAG )

        IF ( EFLAG ) THEN
            MESG = 'Problem reading elevated configuration file'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

C.........  Deallocate local memory
        DEALLOCATE( TMPIDX, EISTAT )

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I10, :, 1X ) )

C******************  INTERNAL SUBPROGRAMS  *****************************
 
        CONTAINS

            SUBROUTINE READ_PELVCONFIG( FDEV, READTYPE, EFLAG )

C.............  External functions
            LOGICAL     BLKORCMT
            LOGICAL     CHKREAL
            CHARACTER(2) CRLF
            INTEGER     INDEX1
            INTEGER     GETNLIST
            REAL        STR2REAL

            EXTERNAL    BLKORCMT, CHKREAL, CRLF, INDEX1, GETNLIST, 
     &                  STR2REAL

C.............  Subprogram arguments
            INTEGER     , INTENT (IN) :: FDEV       ! input file unit
            CHARACTER(*), INTENT (IN) :: READTYPE   ! Reading type
            LOGICAL  , INTENT(IN OUT) :: EFLAG      ! true: error found

C.............  Subprogram local allocatable arrays
            CHARACTER(32), ALLOCATABLE :: SEGMENT( : )

C.............  Local variables
            INTEGER       I, I1, I2, I3, J, K, L, N, V  ! counters and indices

            INTEGER       IOS      ! i/o status
            INTEGER       IREC     ! record counter
            INTEGER       MX_IDX   ! maximum of stack parm indices
            INTEGER    :: NS = 1   ! no. segments in line
            INTEGER       RCNT     ! record count

            REAL          VAL      ! tmp value

            LOGICAL, SAVE :: FIRSTSET = .TRUE. ! true: 1st time group crit found

            CHARACTER(600) BUFFER   ! tmp line buffer as uppercase
            CHARACTER(600) LINE     ! tmp line buffer
            CHARACTER(300) MESG     ! mesg buffer

C----------------------------------------------------------------------

C.............  Rewind input file
            REWIND( FDEV )

C.............  Set local constant variables
            MX_IDX = MAX( HT_IDX, DM_IDX, TK_IDX, VE_IDX, FL_IDX,
     &                    RISE_IDX, SRC_IDX, FIP_IDX, PLT_IDX )

C.............  Read file with different steps depending on READTYPE arguments
            DO IREC = 1, NLINES

                READ( FDEV, 93000, END=999, IOSTAT=IOS ) LINE

                IF ( IOS .NE. 0 ) THEN
                    EFLAG = .TRUE.
                    WRITE( MESG,94010 ) 
     &                'I/O error', IOS, 
     &                'reading elevated source configuration file '//
     &                'at line', IREC
                    CALL M3MESG( MESG )
                    CYCLE
                END IF

C.................  Skip blank lines and comment lines
                IF( BLKORCMT( LINE ) ) CYCLE

C.................  Screen for appended comments and remove them
                CALL RMCOMMNT( '##', LINE )

C.................  Left-justify and convert line to upper case
                BUFFER = ADJUSTL( LINE )
                CALL UPCASE( BUFFER )

C.................  Deallocate segment from previous iteration
                IF ( ALLOCATED( SEGMENT ) ) DEALLOCATE( SEGMENT )

C.................  Allocate memory for segments and parse line into segments
                L = LEN_TRIM( BUFFER )
                NS = GETNLIST( L, BUFFER )
                ALLOCATE( SEGMENT( NS ), STAT=IOS )
                CALL CHECKMEM( IOS, 'SEGMENT', PROGNAME )

                CALL PARSLINE( BUFFER, NS, SEGMENT )

C.................  Interpret line of code.  Set global variables in MODREPRT.
                CALL PRCLINRC( IREC, NS, BUFFER, SEGMENT )

C.................  Skip line if it is the start or end of a packet
                IF ( .NOT. INSPCIFY .OR. LIN_SPCIFY ) CYCLE

C.................  Store the maximum number of ANDs, set the number of ORs (the
C                   number of records), and update the number of variables
                SELECT CASE( PKT_IDX )
                CASE( ELG_IDX )
                    NGRPCRIT = SPCF_NOR
                    MXGRPCHK = MAX( MXGRPCHK, SPCF_NAND )

C.................  Select PinG sources
                CASE( PNG_IDX )
                    NPNGCRIT = SPCF_NOR
                    MXPNGCHK = MAX( MXPNGCHK, SPCF_NAND )

C.....................  Search segments for a match with pollutant names
                    N = 0
                    DO I = 1, NS, 4
                        N = N + 1
                        I1 = ( N-1 )*4 + 1
                        I2 = ( N-1 )*4 + 2

                        J = INDEX1( SEGMENT( I1 ), NIPOL, EINAM )

C.........................  If pollutant found with ranking, increase pollutant
C                           ranking counter
                        IF ( J             .GT. 0     .AND. 
     &                       SEGMENT( I2 ) .EQ. 'TOP'       ) THEN
                            IF ( .NOT. EISTAT(J) ) NALLPOL = NALLPOL + 1
                            EISTAT( J ) = .TRUE.

C.........................  If pollutant found without ranking, increase
C                           pollutant value counter
                        ELSE IF ( J .GT. 0 ) THEN
                            IF ( .NOT. EISTAT(J) ) NALLPOL = NALLPOL + 1
                            EISTAT( J ) = .TRUE.

                        END IF
                    END DO

C.................  Select elevated sources
                CASE( ELV_IDX )
                    NELVCRIT = SPCF_NOR
                    MXELVCHK = MAX( MXELVCHK, SPCF_NAND )

C.....................  Search segments for a match with pollutant names, and
C                       increase counter if found
                    N = 0
                    DO I = 1, NS, 4
                        N = N + 1
                        I1 = ( N-1 )*4 + 1
                        I2 = ( N-1 )*4 + 2

                        J = INDEX1( SEGMENT( I1 ), NIPOL, EINAM )

C.........................  If pollutant found with ranking, increase pollutant
C                           ranking counter
                        IF ( J             .GT. 0     .AND. 
     &                       SEGMENT( I2 ) .EQ. 'TOP'       ) THEN
                            IF ( .NOT. EISTAT(J) ) NALLPOL = NALLPOL + 1
                            EISTAT( J ) = .TRUE.

C.........................  If pollutant found without ranking, increase
C                           pollutant value counter
                        ELSE IF ( J .GT. 0 ) THEN
                            IF ( .NOT. EISTAT(J) ) NALLPOL = NALLPOL + 1
                            EISTAT( J ) = .TRUE.

                        END IF
                    END DO

                END SELECT

C.................  End loop if just counting the number of records
                IF ( READTYPE .EQ. 'COUNT' ) CYCLE

C.................  Check fields
                N = 0
                DO I = 1, NS, 4
                    N  = N + 1
                    I1 = ( N-1 )*4 + 1
                    I2 = ( N-1 )*4 + 2
                    I3 = ( N-1 )*4 + 3

                    IF ( I2 .GT. NS .OR. I3 .GT. NS ) THEN
                        EFLAG = .TRUE.
                        WRITE( MESG,94010 ) 'ERROR: Bad use of ' //
     &                         'delimeters at line ', IREC, 
     &                         '. Line skipped.'
                        CALL M3MSG2( MESG )
                        CYCLE
                    END IF

C.....................  Check and store first part of each AND component 
                    SELECT CASE( SEGMENT( I1 ) )
                    CASE( 'HEIGHT', 'HT' )
                        K = HT_IDX
                    CASE( 'DIAMETER', 'DM' )
                        K = DM_IDX
                    CASE( 'TEMPERATURE', 'TK' )
                        K = TK_IDX
                    CASE( 'VELOCITY', 'VE' )
                        K = VE_IDX
                    CASE( 'FLOW', 'FL' )
                        K = FL_IDX
                    CASE( 'RISE' )
                        K = RISE_IDX

C.........................  Make sure RISE is not used for setting groups
                        IF( PKT_IDX .EQ. ELG_IDX ) THEN
                            WRITE( MESG,94010 ) 
     &                            'WARNING: Variable "RISE" '// 
     &                            'cannot be used to set elevated '//
     &                            'groups at line', IREC  
                            CALL M3MSG2( MESG )
                            CYCLE
                        END IF

C.........................  Flag file as using cutoff approach
                        LCUTOFF = .TRUE.

                    CASE( 'SOURCE' )
                        K = SRC_IDX

C.........................  Make sure SOURCE is not used for setting groups
                        IF( PKT_IDX .EQ. ELG_IDX ) THEN
                            WRITE( MESG,94010 ) 
     &                            'WARNING: Variable "SOURCE" '// 
     &                            'cannot be used to set elevated '//
     &                            'groups at line', IREC  
                            CALL M3MSG2( MESG )
                            CYCLE
                        END IF

                    CASE( 'FIPS' )
                        K = FIP_IDX

C.........................  Make sure FIPS is not used for setting groups
                        IF( PKT_IDX .EQ. ELG_IDX ) THEN
                            WRITE( MESG,94010 ) 
     &                            'WARNING: Variable "FIPS" '// 
     &                            'cannot be used to set elevated '//
     &                            'groups at line', IREC  
                            CALL M3MSG2( MESG )
                            CYCLE
                        END IF

                    CASE( 'PLANT' )
                        K = PLT_IDX

C.........................  Make sure PLANT is not used for setting groups
                        IF( PKT_IDX .EQ. ELG_IDX ) THEN
                            WRITE( MESG,94010 ) 
     &                            'WARNING: Variable "PLANT" '// 
     &                            'cannot be used to set elevated '//
     &                            'groups at line', IREC  
                            CALL M3MSG2( MESG )
                            CYCLE
                        END IF

C.....................  Make sure comparison type is consistent with 
C                       plant
                        IF( SEGMENT( I2 ) .NE. 'IS' ) THEN
                            EFLAG = .TRUE.
                            L = LEN_TRIM( SEGMENT( I2 ) )
                            WRITE( MESG,94010 ) 
     &                        'ERROR: Comparison type "'// 
     &                        SEGMENT(I2)(1:L) // '" at line', IREC, 
     &                        'is not allowed with PLANT criteria.'
                            CALL M3MSG2( MESG )
                        END IF

C.....................  Otherwise, determine if the field is a pollutant, and
C                       if so, store index accordingly
                    CASE DEFAULT

                        V = INDEX1( SEGMENT( I1 ), NIPOL, EINAM )
                        IF ( V .GT. 0 ) THEN

                            SELECT CASE( PKT_IDX )
                            CASE( ELG_IDX )
                                L = LEN_TRIM( EINAM( V ) )
                                WRITE( MESG,94010 ) 
     &                            'WARNING: Pollutant "'// EINAM(V)(1:L)
     &                            // '" cannot be used to set ' //
     &                            'elevated groups at line', IREC
                                CALL M3MSG2( MESG )
                                CYCLE       ! To next line of file

                            CASE( PNG_IDX )
                                 EVPPSTAT( V ) = .TRUE.

                            CASE( ELV_IDX )
                                 EVPESTAT( V ) = .TRUE.

                            END SELECT

                            K = MX_IDX + TMPIDX( V )

C..........................  Otherwise, error
                        ELSE
                            L = LEN_TRIM( SEGMENT( I1 ) )
                            EFLAG = .TRUE.
                            WRITE( MESG,94010 ) 
     &                        'ERROR: Criteria "'// SEGMENT(I1)(1:L) //
     &                        '" at line', IREC, 'is not recognized.'
                            CALL M3MSG2( MESG )

                        END IF

                    END SELECT

C.....................  Check and store second part of each AND component
C.....................  These must match what is recognized in Evalcrit routine
                    J = INDEX1( SEGMENT( I2 ), NCRTSYBL, CRTSYBL )
                    IF ( J .LE. 0 ) THEN
                        EFLAG = .TRUE.
                        L = LEN_TRIM( SEGMENT( I2 ) )
                        WRITE( MESG,94010 ) 
     &                         'ERROR: Comparison type "'// 
     &                         SEGMENT(I2)(1:L) //
     &                         '" at line', IREC, 'is not recognized.'
                        CALL M3MSG2( MESG )

C.....................  Store type
                    ELSE
                        SELECT CASE( PKT_IDX )
                        CASE( ELG_IDX )

C.............................  Make sure TOP is not used for setting groups
                            IF( SEGMENT( I2 ) .EQ. 'TOP' ) THEN
                                WRITE( MESG,94010 ) 
     &                            'WARNING: Comparison type "TOP" '// 
     &                            'ignored at line', IREC, 'for ' //
     &                            'stack group specification'
                                CALL M3MSG2( MESG )
                                CYCLE
C.............................  Otherwise store it
                            ELSE

                                IF( FIRSTSET ) THEN
                                    GRPTYPES = ' '
                                    FIRSTSET = .FALSE.
                                END IF
                                GRPTYPES(NGRPCRIT,N,K) = SEGMENT( I2 )

                            END IF

                        CASE( PNG_IDX )
                            PNGTYPES( NPNGCRIT, N, K ) = SEGMENT( I2 )
                            IF( SEGMENT(I2) .EQ. 'TOP' ) LPNGRNK= .TRUE.
                        CASE( ELV_IDX )
                            ELVTYPES( NELVCRIT, N, K ) = SEGMENT( I2 )
                            IF( SEGMENT(I2) .EQ. 'TOP' ) LELVRNK= .TRUE.
                        END SELECT

                    END IF

C.....................  Check value of each AND component...
C.....................  If plant is specified, then store string information...
                    IF ( K .EQ. PLT_IDX ) THEN

                        SELECT CASE( PKT_IDX )
                        CASE( ELG_IDX )
                            
                        CASE( PNG_IDX )
                            PNGCHRS( NPNGCRIT, N, K ) = SEGMENT( I3 )
                        CASE( ELV_IDX )
                            ELVCHRS( NELVCRIT, N, K ) = SEGMENT( I3 )
                        END SELECT

C.....................  If plant not specified, then check if integer or real
                    ELSE IF ( .NOT. CHKREAL( SEGMENT( I3 ) ) ) THEN
                        EFLAG = .TRUE.
                        L = LEN_TRIM( SEGMENT( I3 ) )
                        WRITE( MESG,94010 ) 
     &                         'ERROR: Value "'// SEGMENT(I3)(1:L) //
     &                         '" at line', IREC, 'is an invalid '//
     &                         'numeric value.'
                        CALL M3MSG2( MESG )

C.....................  Store value of each AND component as a real value
                    ELSE

                        VAL = STR2REAL( SEGMENT( I3 ) )

                        SELECT CASE( PKT_IDX )
                        CASE( ELG_IDX )
                            GRPVALS( NGRPCRIT, N, K ) = VAL
                        CASE( PNG_IDX )
                            PNGVALS( NPNGCRIT, N, K ) = VAL
                        CASE( ELV_IDX )
                            ELVVALS( NELVCRIT, N, K ) = VAL
                        END SELECT

                    END IF

                END DO  ! End loop over ANDs on line

            END DO      ! End loop over lines of file

C.............  Check if there are no group criteria and write warning
            IF( READTYPE .EQ. 'COUNT' .AND.
     &          NGRPCRIT .LE. 0             ) THEN
                MESG = 'WARNING: No grouping criteria in PELVCONFIG ' //
     &                 'will likely result in wasteful ' //
     &                 CRLF() // BLANK10 // ' separation of PinG '//
     &                 'and/or elevated records that could be grouped.'
                CALL M3MSG2( MESG )     
            END IF

            RETURN

C.............  Problem(s) reading input file...
999         WRITE( MESG,94010 ) 'INTERNAL ERROR: Unexpected end of ' //
     &             'file at line', IREC
            CALL M3MSG2( MESG )
            CALL M3EXIT( PROGNAME, 0, 0, ' ', 2 )

C......................  FORMAT  STATEMENTS   ..........................

C...........   Formatted file I/O formats............ 93xxx
93000       FORMAT( A )

C...............   Internal buffering formats............ 94xxx

94010       FORMAT( 10( A, :, I10, :, 1X ) )

            END SUBROUTINE READ_PELVCONFIG

        END SUBROUTINE RPELVCFG
