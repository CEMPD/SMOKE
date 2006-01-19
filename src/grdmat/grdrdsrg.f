
        SUBROUTINE GRDRDSRG( NT, TMPLINE, VFLAG ) 

C***********************************************************************
C  subroutine body starts at line 117
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
C**************************************************************************
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

C...........   Modules for public variables
C...........   This module contains the gridding surrogates tables
        USE MODSURG, ONLY: IDXSRGA, IDXSRGB, SCELLA, SFIPSA, SSRGIDA,
     &                     SFRACA, NSRGREC, NSRGFIPS, SRGFIPS, FIPCELL,
     &                     NCELLS, SRGFRAC, SRGCSUM, SRGFMT, SRGNCOLS,
     &                     SRGNROWS

C.........  This module contains the global variables for the 3-d grid
        USE MODGRID, ONLY: XOFF, YOFF, NCOLS, NROWS

        IMPLICIT NONE

C...........   INCLUDES
        INCLUDE 'EMCNST3.EXT'   !  emissions constat parameters

C...........   EXTERNAL FUNCTIONS and their descriptions:

        CHARACTER(2)   CRLF
        INTEGER        FIND1
        INTEGER        STR2INT
        REAL           STR2REAL
        LOGICAL        BLKORCMT

        EXTERNAL       CRLF, FIND1, STR2INT, STR2REAL, BLKORCMT

C...........   Subroutine arguments
        INTEGER      , INTENT  (IN) :: NT            ! no of surrogate entry
        CHARACTER(*) , INTENT  (IN) :: TMPLINE( NT ) ! tmp line buffer
        LOGICAL      , INTENT  (IN) :: VFLAG         ! true: using variable grid

C...........   Local parameters

        INTEGER, PARAMETER :: MXSEG = 5           ! # of potential line segments

C...........   Other arrays

        CHARACTER(20) SEGMENT( MXSEG )             ! Segments of parsed lines

C...........   Local variables
        INTEGER         I, F, J, K, L, N         ! indices and counters

        INTEGER         C                     ! tmp cell number
        INTEGER         CELCNT                ! cell counter
        INTEGER         CNTCHK                ! check for message overflow count
        INTEGER         COL                   ! Temp grid column number (x-axis)
        INTEGER         FIP                   ! tmp country/state/county code
        INTEGER         IOS                   ! i/o status
        INTEGER         LCEL                  ! cell ID from previous iteration
        INTEGER         LFIP                  ! county code from prev iteration
        INTEGER         LSSC                  ! srg ID from previous iteration
        INTEGER         NSRGALL               ! No. entries in surrgoates file
        INTEGER         MXCFIP                ! Max cells per county code
        INTEGER         ROW                   ! Tmp grid row number (y-axis)
        INTEGER         SSC                   ! Tmp spatial surrogate code

        LOGICAL, SAVE :: FIRSTIME = .TRUE.

        REAL            RATIO                 ! Temp spatial surrogate Ratio

        LOGICAL      :: GFLAG = .FALSE.       ! true: per iteration surg > 1
        LOGICAL      :: HFLAG = .FALSE.       ! true: header found
        LOGICAL      :: OFLAG = .FALSE.       ! true: overall warning flag
        LOGICAL      :: RFLAG = .FALSE.       ! true: renormalized surrogates
        LOGICAL         WFLAG                 ! true: per iteration warning flag

        CHARACTER(20)   COLRANGE              ! buffer w/ column range
        CHARACTER(20)   ROWRANGE              ! buffer w/ row range
        CHARACTER(80)   LINE                  ! Read buffer for a line
        CHARACTER(200)  MESG                  ! Message buffer

        CHARACTER(16) :: PROGNAME = 'GRDRDSRG'    !  program name

C***********************************************************************
C   Begin body of subroutine GRDRDSRG

        NSRGREC = NT
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

C.........  Initialize arrays
        IDXSRGA = 0  ! array
        IDXSRGB = 0  ! array
        SCELLA  = 0  ! array
        SFIPSA  = 0. ! array
        SSRGIDA = 0. ! array
        SFRACA  = 0. ! array

C.........  Fill surrogate arrays

        SELECT CASE( SRGFMT )

        CASE( 'MODELS3' )

            HFLAG = .FALSE.
            J = 0
            DO I = 1, NT

                LINE = TMPLINE( I )

                IF ( .NOT. HFLAG ) CALL UPCASE( LINE )

C............  Parse the line of data into segments based on the rules
C              for "list-formatted" in fortran, but not requiring 
C              quotes around the text strings

                CALL PARSLINE( LINE, MXSEG, SEGMENT )

                SSC    = STR2INT ( SEGMENT( 1 ) )
                FIP    = STR2INT ( SEGMENT( 2 ) )
                COL    = STR2INT ( SEGMENT( 3 ) )
                ROW    = STR2INT ( SEGMENT( 4 ) )
                RATIO  = STR2REAL( SEGMENT( 5 ) )
                
                WFLAG = .FALSE.
                
C.................  Adjust column and row for subgrid
                COL = COL - XOFF
                ROW = ROW - YOFF

C.................  Skip entry after subgrid adjustment
                IF( COL .LE. 0 .OR. COL .GT. NCOLS .OR.
     &              ROW .LE. 0 .OR. ROW .GT. NROWS       ) CYCLE

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
            
            CONTINUE

            NSRGALL = J

C.............  Write out final warning for row/col out of range
            IF( OFLAG ) THEN
                MESG = 'WARNING: Lines skipped in surrogate file ' //
     &                 'because rows or columns were out of range.'
                CALL M3MSG2( MESG )
            END IF

        CASE DEFAULT
            MESG = 'INTERNAL ERROR: Surrogate format "' // SRGFMT //
     &             '" not understood'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

        END SELECT

C.........  Now create the derived surrogates tables from the original data...

C.........  Sort surrogates by county code & cell & surrogate code
        CALL SORTI3( NSRGALL, IDXSRGA, SFIPSA, SCELLA, SSRGIDA )

C.........  Count county codes in surrogates file and maximum number of cells
C       per cy/st/co code.
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
      
        K = 1            ! single srg code per run
  
C.........  Allocate memory for derived surrogates tables
        ALLOCATE( NCELLS( NSRGFIPS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'NCELLS', PROGNAME )
      
        ALLOCATE( SRGFIPS( NSRGFIPS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'SRGFIPS', PROGNAME )

        ALLOCATE( FIPCELL( MXCFIP, NSRGFIPS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'FIPCELL', PROGNAME )
        
        ALLOCATE( SRGFRAC( K, MXCFIP, NSRGFIPS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'SRGFRAC', PROGNAME )
      
        ALLOCATE( SRGCSUM( K, NSRGFIPS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'SRGCSUM', PROGNAME )
      
C.........  Initialize arrays
        NCELLS  = 0  ! array
        SRGFIPS = 0  ! array
        FIPCELL = 0  ! array
        SRGFRAC = 0. ! array
        SRGCSUM = 0. ! array
      
C.............  Store the surrogate fractions, FIPS codes, and cell numbers...
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
      
            SRGFRAC( K, CELCNT, NSRGFIPS ) = RATIO
            SRGCSUM( K, NSRGFIPS ) = SRGCSUM( K, NSRGFIPS ) + RATIO
      
        END DO
      
C.........  Now check to make sure the surrogate totals for each county are 
C           less than or equal to 1.
      
        DO I = 1, NSRGFIPS
      
            GFLAG = .FALSE.
            CNTCHK = 0
      
C.........  Check if county total surrogates greater than 1
            IF( SRGCSUM( K,I ) .GT. 1.001 ) THEN
      
C.........  If first problem on this line
                IF( .NOT. GFLAG ) THEN
                    WRITE( MESG,94030 ) 'WARNING: County ' //
     &              'surrogate total greater than 1. for '//
     &              'county' // CRLF() // BLANK10, SRGFIPS( I ), 
     &              ', SSC(', K , '):', SRGCSUM( K,I )
                    GFLAG = .TRUE.
                    CNTCHK = CNTCHK + 1
      
C.............  If multiple problems on this line, but fewer than
C               the MESG length will permit, add to message
                ELSE IF( CNTCHK .LE. 27 ) THEN
                    L = LEN_TRIM( MESG )
                    WRITE( MESG,94031 ) MESG( 1:L )// ', SSC(', 
     &                       K , '):', SRGCSUM( K,I )
                    CNTCHK = CNTCHK + 1
      
                END IF
      
           END IF
      
C.............  Renormalize all surrogates with totals > 1
            IF( SRGCSUM( K,I ) .GT. 1. ) THEN
                RFLAG = .TRUE.  ! Set to give global warning
                DO C = 1, NCELLS( I )
                    SRGFRAC( K,C,I ) = SRGFRAC( K,C,I ) /
     &                                 SRGCSUM( K,I )
                END DO
            END IF
      
C.........  Give a warning message for significant counties
            IF( GFLAG ) CALL M3MESG( MESG )

        END DO
      
C.........  Reset number of surrogate records stored in the module with
C           the correct number after reading file and removing records that
C           are outside the subgrid (if any)

C.....  Deallocate local variables

        DEALLOCATE( IDXSRGA, IDXSRGB, SCELLA, SFIPSA, SSRGIDA, SFRACA )

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Formatted file I/O formats............ 93xxx

93000   FORMAT( A )


C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

94020   FORMAT( A, 1X, I8, 1X, A, 1X, F10.6, 1X, A )

94030   FORMAT( A, 1X, I6.6, A, I3.2, A, F8.4 )

94031   FORMAT( A, 1X, I3.2, A, F8.4 )

        END SUBROUTINE GRDRDSRG