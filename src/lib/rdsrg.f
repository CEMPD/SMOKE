
        SUBROUTINE RDSRG( FDEV, SRGFMT, XCENT, YCENT, XORIG, YORIG,
     &                    XCELL, YCELL, NCOLS, NROWS ) 

C***********************************************************************
C  subroutine body starts at line
C
C  DESCRIPTION:
C      This subroutine allocates memory for the spatial surrogate
C      arrays, reads the the spatial surrogates file, and then fills
C      the spatial surrogate arrays.
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
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

C...........   Modules for public variables
C...........   This module contains the gridding surrogates tables
        USE MODSURG

        IMPLICIT NONE

C...........   EXTERNAL FUNCTIONS and their descriptions:

        INTEGER        FIND1
        INTEGER        STR2INT
        REAL           STR2REAL

        EXTERNAL       FIND1, STR2INT, STR2REAL

C...........   Subroutine arguments

        INTEGER      , INTENT  (IN) :: FDEV       ! File unit number
        CHARACTER(*) , INTENT  (IN) :: SRGFMT     ! Format of surrogates file
        REAL         , INTENT  (IN) :: XCENT      ! Center of coordinate system
        REAL         , INTENT  (IN) :: YCENT      ! Center of coordinate system
        REAL         , INTENT  (IN) :: XORIG      ! X origin
        REAL         , INTENT  (IN) :: YORIG      ! Y origin
        REAL         , INTENT  (IN) :: XCELL      ! Cell size, X direction
        REAL         , INTENT  (IN) :: YCELL      ! Cell size, Y direction
        INTEGER      , INTENT  (IN) :: NCOLS      ! # cells in X direction
        INTEGER      , INTENT  (IN) :: NROWS      ! # cells in Y direction

C...........   Local parameters

        INTEGER, PARAMETER :: MXSEG = 5           ! # of potential line segments

C...........   Other arrays

        CHARACTER*20 SEGMENT( MXSEG )             ! Segments of parsed lines

C...........   Local variables
        INTEGER         I, J, K               ! indices and counters

        INTEGER         C                     ! tmp cell number
        INTEGER         CELCNT                ! cell counter
        INTEGER         COL                   ! Temp grid column number (x-axis)
        INTEGER         FIP                   ! tmp country/state/county code
        INTEGER         IOS                   ! i/o status
        INTEGER         IREC                  ! Record counter
        INTEGER         LC, LR                ! length of COLRANGE & ROWRANGE
        INTEGER         LCEL                  ! cell ID from previous iteration
        INTEGER         LFIP                  ! county code from prev iteration
        INTEGER         LSSC                  ! srg ID from previous iteration
        INTEGER         NSRGALL               ! No. entries in surrgoates file
        INTEGER         MXCFIP                ! Max cells per county code
        INTEGER         ROW                   ! Temp grid row number (y-axis)
        INTEGER         SSC                   ! Temp spatial surrogate code

        REAL            RATIO                 ! Temp spatial surrogate Ratio

        LOGICAL      :: EFLAG = .FALSE.       ! true: error found
        LOGICAL      :: HFLAG = .FALSE.       ! true: header found
        LOGICAL      :: OFLAG = .FALSE.       ! true: overall warning flag
        LOGICAL         WFLAG                 ! true: per iteration warning flag

        CHARACTER*20    COLRANGE              ! buffer w/ column range
        CHARACTER*20    ROWRANGE              ! buffer w/ row range
        CHARACTER*80    LINE                  ! Read buffer for a line
        CHARACTER*300   MESG                  ! Message buffer

        CHARACTER*16 :: PROGNAME = 'RDSRG'    !  program name

C***********************************************************************
C   Begin body of subroutine RDSRG

        IREC    = 0
	J       = 0
        NSRGREC = 0

C......... Determine the number surrogate file entries

        REWIND( FDEV )

        SELECT CASE( SRGFMT )
	
        CASE( 'MODELS3' )

            HFLAG = .FALSE.
            DO

                READ ( FDEV, 93000, END=12, IOSTAT=IOS ) LINE

                IREC = IREC + 1
             
                IF ( IOS .GT. 0 ) THEN
                    WRITE( MESG, 94010)
     &                'I/O error', IOS, 'reading gridding surrogates '//
     &                'file at line', IREC
                    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                END IF

                IF ( .NOT. HFLAG ) THEN

                    CALL UPCASE( LINE )

                    I = INDEX( LINE, '#GRID' )! check current line for header

                    IF( I .GT. 0 ) THEN       ! skip if current line is header
                        HFLAG = .TRUE.
                        CYCLE
                    END IF

                ELSEIF ( LINE .EQ. ' ' ) THEN ! skip if current line is blank
                    CYCLE
                ELSE
                    NSRGREC = NSRGREC + 1
                END IF

            END DO

12          CONTINUE    ! end of read on input file

        END SELECT

C......... Allocate memory for surrogate arrays

        ALLOCATE( IDXSRGA( NSRGREC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'IDXSRGA', PROGNAME )

        ALLOCATE( IDXSRGB( NSRGREC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'IDXSRGB', PROGNAME )

        ALLOCATE( SCELLA( NSRGREC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'SCELLA', PROGNAME )

        ALLOCATE( SFIPSA( NSRGREC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'SFIPSA', PROGNAME )

        ALLOCATE( SSRGIDA( NSRGREC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'SSRGIDA', PROGNAME )

        ALLOCATE( SFRACA( NSRGREC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'SFRACA', PROGNAME )

C.........  Create message fields for errors
        WRITE( COLRANGE, '( "( ", I1, " to ", I4, " )" )' ) 1, NCOLS
        WRITE( ROWRANGE, '( "( ", I1, " to ", I4, " )" )' ) 1, NROWS

        LC = LEN_TRIM( COLRANGE )
        LR = LEN_TRIM( ROWRANGE )

C.........  Fill surrogate arrays

        REWIND( FDEV )
	IREC    = 0

        SELECT CASE( SRGFMT )
	
        CASE( 'MODELS3' )

            HFLAG = .FALSE.
            DO

                READ ( FDEV, 93000, END=24, IOSTAT=IOS ) LINE

                IREC = IREC + 1
             
                IF ( IOS .GT. 0 ) THEN
                    WRITE( MESG, 94010)
     &                'I/O error', IOS, 'reading gridding surrogates '//
     &                'file at line', IREC
                    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                END IF

                IF ( .NOT. HFLAG ) THEN

                    CALL UPCASE( LINE )

                    I = INDEX( LINE, '#GRID' )! check current line for header

                    IF( I .GT. 0 ) THEN   ! skip if current line is header
                        HFLAG = .TRUE.
                        CYCLE
                    END IF

                ELSEIF ( LINE .EQ. ' ' ) THEN ! skip if current line is blank
                    CYCLE

                END IF

C.................  Parse the line of data into segments based on the rules
C                   for "list-formatted" in fortran, but not requiring 
C                   quotes around the text strings

                CALL PARSLINE( LINE, MXSEG, SEGMENT )

                SSC    = STR2INT ( SEGMENT( 1 ) )
                FIP    = STR2INT ( SEGMENT( 2 ) )
                COL    = STR2INT ( SEGMENT( 3 ) )
                ROW    = STR2INT ( SEGMENT( 4 ) )
                RATIO  = STR2REAL( SEGMENT( 5 ) )

                WFLAG = .FALSE.
C.................  Check the value of the column number
                IF( COL .LT. 0 .OR.  COL .GT. NCOLS  .OR.
     &            ( ROW .EQ. 0 .AND. COL .NE. 0     )    ) THEN                    WFLAG = .TRUE.
                    WFLAG = .TRUE.
                    OFLAG = .TRUE.
                    WRITE( MESG,94010 ) 'WARNING: Column value ', COL,
     &                     'is outside range ' // COLRANGE( 1:LC ) // 
     &                     ' at line', IREC
                    CALL M3MESG( MESG )
                END IF

C.................  Check the value of the row number
                IF( ROW .LT. 0 .OR.  ROW .GT. NROWS  .OR.
     &            ( ROW .EQ. 0 .AND. COL .NE. 0     )    ) THEN
                    WFLAG = .TRUE.
                    OFLAG = .TRUE.
                    WRITE( MESG,94010 ) 'WARNING: Row value ', ROW,
     &                     'is outside range ' // ROWRANGE( 1:LR ) // 
     &                     ' at line', IREC
                    CALL M3MESG( MESG )                    

C.................  Special treatment for cell (0,0)
                ELSE IF( ROW .EQ. 0 .AND. COL. EQ. 0 ) THEN
                    CYCLE

                END IF

C.................  Check the value of the ratio value
                IF( RATIO .GT. 1. ) THEN
                    EFLAG = .TRUE.
                    WRITE( MESG,94010 )
     &                     'ERROR: surrogates ratio is greater than ' //
     &                     '1.0 at line', IREC
                    CALL M3MESG( MESG )
                END IF

C.................  Skip entry if rows and columns are out of range
                IF( WFLAG ) CYCLE

                J = J + 1 

                IDXSRGA( J ) = J
                IDXSRGB( J ) = J
                SCELLA ( J ) = (ROW-1)*NCOLS + COL
                SFIPSA ( J ) = FIP
                SSRGIDA( J ) = SSC
                SFRACA ( J ) = RATIO     

            END DO

24          CONTINUE    ! end of read on input file

            NSRGALL = J

C.............  Write out final warning for row/col out of range
            IF( OFLAG ) THEN
                MESG = 'Lines skipped in surrogate file because ' //
     &                 'rows or columns were out of range.'
                CALL M3WARN( PROGNAME, 0, 0, MESG )
            END IF

        CASE DEFAULT
            MESG = 'INTERNAL ERROR: Surrogate format "' // SRGFMT //
     &             '" not understood'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

        END SELECT

        IF( EFLAG ) THEN
            MESG = 'Problem reading gridding surrogates file'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

C.........  Now create the derived surrogates tables from the original data...

C.........  Sort surrogates by county code & cell & surrogate code
        CALL SORTI3( NSRGALL, IDXSRGA, SFIPSA, SCELLA, SSRGIDA )

C.........  Sort surrogates by surrogate code
        CALL SORTI1( NSRGALL, IDXSRGB, SSRGIDA )

C.........  Count county codes in surrogates file and maximum number of cells
C           per cy/st/co code.
        LFIP     = -1
        LCEL     = -1
        MXCFIP   = 0
        NSRGFIPS = 0
        CELCNT   = 0
        DO I = 1, NSRGALL

            J   = IDXSRGA( I )
            FIP = SFIPSA ( J )
            C   = SCELLA ( J )

            IF( FIP .NE. LFIP ) THEN

                IF( CELCNT .GT. MXCFIP ) MXCFIP = CELCNT

                NSRGFIPS = NSRGFIPS + 1  ! incrmt cntr for county in srg file
                CELCNT   = 0             ! init cell counter per county
                LFIP     = FIP
              
            END IF

            IF( C   .NE. LCEL .OR. CELCNT .EQ. 0 ) THEN
                CELCNT = CELCNT + 1      ! increment cell counter per county
                LCEL = C
            END IF

        END DO

        IF( CELCNT .GT. MXCFIP ) MXCFIP = CELCNT

C.........  Count surrogate codes in surrogates file
        LSSC  = -1
        NSRGS = 0
        DO I = 1, NSRGALL

            J   = IDXSRGB( I )
            SSC = SSRGIDA( J )

            IF( SSC .NE. LSSC ) THEN
                NSRGS = NSRGS + 1
                LSSC = SSC
            END IF

        END DO

C.........  Allocate memory for derived surrogates tables
        ALLOCATE( NCELLS( NSRGFIPS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'NCELLS', PROGNAME )

        ALLOCATE( SRGFIPS( NSRGFIPS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'SRGFIPS', PROGNAME )

        ALLOCATE( SRGLIST( NSRGS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'SRGLIST', PROGNAME )

        ALLOCATE( FIPCELL( MXCFIP, NSRGFIPS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'FIPCELL', PROGNAME )

        ALLOCATE( SRGFRAC( NSRGS, MXCFIP, NSRGFIPS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'SRGFRAC', PROGNAME )

        ALLOCATE( SRGCSUM( NSRGS, NSRGFIPS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'SRGCSUM', PROGNAME )

C.........  Initialize arrays
        NCELLS  = 0  ! array
        SRGFIPS = 0  ! array
        SRGLIST = 0  ! array
        FIPCELL = 0  ! array
        SRGFRAC = 0. ! array
        SRGCSUM = 0. ! array

C.........  Store derived surrogates tables...

C.........  Store the surrogate ID list
        LSSC  = -1
        NSRGS = 0
        DO I = 1, NSRGALL

            J    = IDXSRGB( I )
            SSC = SSRGIDA( J )

            IF( SSC .NE. LSSC ) THEN
                NSRGS = NSRGS + 1
                SRGLIST( NSRGS ) = SSC
                LSSC = SSC
            END IF

        END DO

C.........  Initialize arrays that might not be totally populated
        FIPCELL = 0   ! array
        SRGFRAC = 0   ! array

C.........  Store the surrogate fractions, FIPS codes, and cell numbers...
        LFIP     = -1
        LCEL     = -1
        NSRGFIPS = 0
        DO I = 1, NSRGALL

            J     = IDXSRGA( I )
            FIP   = SFIPSA ( J )
            C     = SCELLA ( J )
            SSC   = SSRGIDA( J )
            RATIO = SFRACA ( J )

            IF( FIP .NE. LFIP ) THEN

                NSRGFIPS = NSRGFIPS + 1  ! incrmt cntr for county in srg file
                SRGFIPS( NSRGFIPS ) = FIP
                CELCNT   = 0             ! init cell counter per county
                LFIP     = FIP
              
            END IF

            IF( C   .NE. LCEL .OR. CELCNT .EQ. 0 ) THEN
                CELCNT = CELCNT + 1      ! increment cell counter per county
                LCEL = C
            END IF

C.............  Store cell count at every iteration, and it will get 
C               overwritten until the maximum value for each county is stored
            NCELLS ( NSRGFIPS ) = CELCNT  

C.............  Store cell number  at every iteration, but it will be the same
C               each time
            FIPCELL( CELCNT, NSRGFIPS ) = C  

C.............  Find surrogate code in sorted list and use position to store
C               the surrogates fraction
            K = FIND1( SSC, NSRGS, SRGLIST )

            SRGFRAC( K, CELCNT, NSRGFIPS ) = RATIO
            SRGCSUM( K, NSRGFIPS ) = SRGCSUM( K, NSRGFIPS ) + RATIO

        END DO

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Formatted file I/O formats............ 93xxx

93000   FORMAT( A )

94010   FORMAT( 10( A, :, I8, :, 1X ) )

        END
