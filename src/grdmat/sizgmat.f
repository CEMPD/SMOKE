
        SUBROUTINE SIZGMAT( CATEGORY, NSRC, VFLAG, DEFSRGID, SRGFLAG,
     &                      MXSCEL, MXCSRC, MXCCL, NMATX, NMATXU )

C***********************************************************************
C  subroutine body starts at line 102
C
C  DESCRIPTION:
C      This subroutine determines the sizes needed for creating the for the 
C      area and mobile gridding matrices
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
C       Updated with USE M3UTILIO by Huy Tran UNC-IE on 2026-01
C***************************************************************************

C...........   MODULES for public variables
C...........   This module is the source inventory arrays
        USE M3UTILIO

        USE MODSOURC, ONLY: XLOCA, YLOCA, CIFIP, CELLID, CLINK,
     &                      XLOC1, YLOC1, XLOC2, YLOC2

C...........   This module contains the cross-reference tables
        USE MODXREF, ONLY: ASRGID

C.........  This module contains the global variables for the 3-d grid
        USE MODGRID, ONLY: NGRID, NCOLS, NROWS, XOFF, YOFF

C...........   This module contains the gridding surrogates tables
        USE MODSURG, ONLY: NCELLS, FIPCELL, NSRGS, SRGLIST, NSRGFIPS,
     &                     SRGFIPS, NTSRGDSC, SRGFNAM, SRGFCOD, NTLINES,
     &                     MXCFIP, SRGNCOLS, SRGNROWS, SRGCSUM, SRGFRAC

        USE MODGRDLIB

        IMPLICIT NONE

C...........   INCLUDES:
        
        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
C        INCLUDE 'PARMS3.EXT'    !  I/O API parameters

C...........   EXTERNAL FUNCTIONS and their descriptions:
C       CHARACTER(2)    CRLF
C       INTEGER         ENVINT
C       INTEGER         FIND1
C       INTEGER         FINDC
        LOGICAL         BLKORCMT
C       LOGICAL         SETENVVAR
C       INTEGER         STR2INT
        INTEGER         GETFLINE
C       INTEGER         GETEFILE

C        EXTERNAL        CRLF, FIND1, FINDC, BLKORCMT, SETENVVAR,
C     &                  STR2INT, GETFLINE, ENVINT, GETEFILE
        EXTERNAL     BLKORCMT, GETFLINE


C...........   SUBROUTINE ARGUMENTS
        CHARACTER(*), INTENT (IN) :: CATEGORY  ! source category
        INTEGER     , INTENT (IN) :: NSRC      ! local number of sources
        LOGICAL     , INTENT (IN) :: VFLAG     ! true: using variable grid
        INTEGER     , INTENT (IN) :: DEFSRGID  ! default surrogate code
        LOGICAL     , INTENT (IN) :: SRGFLAG   ! true: using default fallback surrogate
        INTEGER     , INTENT(OUT) :: MXSCEL    ! max sources per cell   
        INTEGER     , INTENT(OUT) :: MXCSRC    ! max cells per source   
        INTEGER     , INTENT(OUT) :: MXCCL     ! max cells per county or link   
        INTEGER     , INTENT(OUT) :: NMATX     ! no. src-cell intersections   
        INTEGER     , INTENT(OUT) :: NMATXU    ! no. county-cell intrsctns for all sources

C...........   Local arrays dimensioned by subroutine arguments
C...........   Note that the NGRID dimension could conceivably be too small if 
C              a link winds through the whole domain, but this is a case that
C              is not worth going to extra trouble for since it is not realistic
        INTEGER         NX  ( NGRID )    ! number of srcs per cell  
        INTEGER         ACEL( NGRID )    ! number of cell intersections per src
        REAL            AFAC( NGRID )    ! fraction of link in cell

C...........   Other local variables
        INTEGER         C, F, J, JJ, K, KK, I, II, N, NT, S   ! counters and indices

        INTEGER         CCNT             ! counters for no. non-zero-surg cells
        INTEGER      :: CELLSRC = 0      ! cell number as source char
        INTEGER         COL              ! tmp column
        INTEGER         TCOL             ! tmp column
        INTEGER         GDEV             !  for surrogate coeff file
        INTEGER         ID1, ID2         ! primary and 2ndary surg codes
        INTEGER         LC, LR           ! length of COLRANGE & ROWRANGE
        INTEGER         ISIDX            ! tmp surrogate ID code index
        INTEGER         IOS              ! i/o status
        INTEGER         IREC             ! Record counter
        INTEGER         ISDEF            !  default surrogate ID code index
        INTEGER         NWARN            ! tmp number of WARNING msg
        INTEGER         NCEL             ! tmp number of cells
        INTEGER         ROW              ! tmp row
        INTEGER         TROW             ! tmp row
        INTEGER         NTL              ! max no. of line buffers
        INTEGER         TGTSRG           ! target surrogates code
        INTEGER         SSC              ! surrogates code
        INTEGER      :: NLINES = 0       ! number of lines in input file

        INTEGER, PARAMETER :: MXSEG = 10           ! # of potential line segments

C...........   Other arrays
        CHARACTER(20) SEGMENT( MXSEG )             ! Segments of parsed lines
        CHARACTER(60), ALLOCATABLE :: TMPLINE( : )   ! tmp line buffer

        REAL            ALEN        ! link length
        REAL*8          XX, YY

        LOGICAL      :: EFLAG = .FALSE. ! true: error flag
        LOGICAL      :: LFLAG = .FALSE. ! true: location data available
        LOGICAL      :: XYSET = .FALSE. ! true: X/Y available for src
        LOGICAL      :: CFLAG = .TRUE.  ! true: called by sizgmat, false: called by gen[a|m]gmat
        LOGICAL      :: WFLAG = .FALSE. ! true: per iteration warning flag

        CHARACTER(FIPLEN3)  CFIP             ! country/state/county code
        CHARACTER(FIPLEN3)  LFIP             ! tmp country/state/county code
        CHARACTER(200)      LINE             ! Read buffer for a line
        CHARACTER(300)      MESG             !  message buffer
        CHARACTER(256)      NAMBUF           !  surrogate file name buffer
        CHARACTER(256)      NAMBUFT          !  tmp surrogate file name buffer
        CHARACTER(256)      TSRGFNAM         !  tmp surrogate file name buffer
        CHARACTER(20)       COLRANGE         ! buffer w/ column range
        CHARACTER(20)       ROWRANGE         ! buffer w/ row range

        CHARACTER(LNKLEN3) :: CLNK = ' '   ! tmp link ID

        CHARACTER(16) :: PROGNAME = 'SIZGMAT' ! program name

C***********************************************************************
C   begin body of subroutine SIZGMAT

C.........  Print status message
        MESG = 'Computing gridding matrix size...'
        CALL M3MSG2( MESG )

C..........  Set flag to indicate that XLOCA/YLOCA are available
        LFLAG = ALLOCATED( XLOCA )

C.........  Initialize the count of sources per cell
        NX = 0   ! array

C.........  Loop through sources
        MXCSRC  = 0
        MXCCL   = 0
        NMATXU  = 0
        CELLSRC = 0

C.........  Allocate memory for indices to surrogates tables for each source
        ALLOCATE( NTLINES( NSRGS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'NTLINES', PROGNAME )

C.........  Create message fields for errors
        WRITE( COLRANGE, '( "( ", I1, " to ", I4, " )" )' ) 1, SRGNCOLS
        WRITE( ROWRANGE, '( "( ", I1, " to ", I4, " )" )' ) 1, SRGNROWS

        LC = LEN_TRIM( COLRANGE )
        LR = LEN_TRIM( ROWRANGE )
 
C.........  Count total line buffers to define memory size
        DO II = 1, NSRGS    ! loop through only the surrogate code assigned by sources

            NTL = 0
            TSRGFNAM = ' '

            TGTSRG = SRGLIST( II )
            KK = II

            ISDEF  = FIND1( DEFSRGID, NSRGS, SRGLIST )

C..............  default fallback surrogate will run for the last for gap filling.
            IF( II >= ISDEF ) THEN

C.................  default fallback surrogate will run at last after
C               re-assigned zero fraction surrogate
                IF( II == NSRGS ) THEN
                    TGTSRG = DEFSRGID
                    KK = ISDEF
                ELSE
                    TGTSRG = SRGLIST( II + 1 )
                    KK = II + 1
                END IF
              
            END IF
            
C.............   Count total line buffers to define memory size of TMPLINE
            DO I = 1, NTSRGDSC  ! Open all surrogate files using the same srg code
       
C..................  Prompt for and open I/O API output file(s)...
                IF( TGTSRG .NE. SRGFCOD( I ) ) CYCLE

                CALL GETENV( 'SRGPRO_PATH', NAMBUF )
                WRITE( NAMBUFT, '( 2A )' ) TRIM( NAMBUF ), SRGFNAM( I )
       
                IF( NAMBUFT .NE. TSRGFNAM  ) THEN
                    CALL OPEN_SRGFILE

                    IREC = 0
C.....................  Reading surrogate files
                    DO JJ = 1, NLINES
             
                       READ ( GDEV, 93000, END=111, IOSTAT=IOS ) LINE
                       IREC = IREC + 1
             
                       IF ( IOS .GT. 0 ) THEN
                            WRITE( MESG, 94010)
     &                      'I/O error', IOS, 'reading gridding ' //
     &                      'surrogates file at line', IREC
                            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                       END IF
             
                       IF ( BLKORCMT( LINE ) ) CYCLE

                       CALL PARSLINE( LINE, MXSEG, SEGMENT )
                       SSC    = STR2INT ( SEGMENT( 1 ) )
                       TCOL   = STR2INT ( SEGMENT( 3 ) )
                       TROW   = STR2INT ( SEGMENT( 4 ) )

C..........................  Check the value of the column number
                        IF( TCOL .LT. 0 .OR.  TCOL .GT. SRGNCOLS  .OR.
     &                    ( TROW .EQ. 0 .AND. TCOL .NE. 0 ) ) THEN
                            WFLAG = .TRUE.
                            WRITE( MESG,94010 ) 'WARNING: Column value',
     &                             TCOL, 'is outside range ' //
     &                             COLRANGE( 1:LC ) // ' from FIPS ' //
     &                             CFIP // ' and surrogate', SSC
                            CALL M3MESG( MESG )
                        END IF

C..........................  Check the value of the row number
                        IF( TROW .LT. 0 .OR.  TROW .GT. SRGNROWS  .OR.
     &                    ( TCOL .EQ. 0 .AND. TROW .NE. 0 ) ) THEN
                            WFLAG = .TRUE.
                            WRITE( MESG,94010 ) 'WARNING: Row value ', 
     &                             TROW, 'is outside range ' // 
     &                             ROWRANGE( 1:LR ) // ' from FIPS ' //
     &                             CFIP // ' and surrogate', SSC
                            CALL M3MESG( MESG )                    

C..........................  Special treatment for cell (0,0) (skip for now)
                        ELSE IF( TROW .EQ. 0 .AND. TCOL. EQ. 0 ) THEN
                            CYCLE

                        END IF

C..........................  Adjust column and row for subgrid
                        TCOL = TCOL - XOFF
                        TROW = TROW - YOFF

C.........................  Skip entry after subgrid adjustment
                        IF( TCOL .LE. 0 .OR. TCOL .GT. NCOLS .OR.
     &                      TROW .LE. 0 .OR. TROW .GT. NROWS ) CYCLE

C.......................... Skip entry if rows and columns are out of range
                        IF( WFLAG ) CYCLE

C.........................  Skip entry if SSC is not in the assigned SRGLIST by source
                        IF( SSC .EQ. TGTSRG ) NTL = NTL + 1
                        
111                 END DO

                    TSRGFNAM = NAMBUFT    ! store a previous surrogate file name buffer
                    CLOSE( GDEV )
       
                END IF      ! skip if surrogate file has the same srg file
            
            END DO       ! loop over all surrogate files in SRGDESC file
            
C............  Store no of line buffers of each surrogate
            NTLINES( KK ) = NTL

C.............  Write the status of reading surrogate files.
            WRITE( MESG,94010 ) 'Reading surrogate', TGTSRG,
     &                          'to define the size of gridding matrix'
            CALL M3MSG2( MESG )

C.............  Warning message when there are no surrogate available.
            IF( NTL .EQ. 0 ) THEN
                WRITE( MESG,94010 ) 'WARNING: The surrogate', TGTSRG,
     &             ' does not exist within a modeling domain'
                CALL M3MESG( MESG )
            END IF

C.............  Allocate memory for indices to surrogates tables for each source
            ALLOCATE( TMPLINE( NTL ), STAT=IOS )
            CALL CHECKMEM( IOS, 'TMPLINE', PROGNAME )

C............  If surrogates are needed, read and store the gridding surrogates, 
C              allocate memory for the surrogate assignments, and assign
C              surrogates to each source.
            NT = 0
            TMPLINE = ' '
            TSRGFNAM = ' '

            DO I = 1, NTSRGDSC  ! Open all surrogate files using the same srg code

C.................  Prompt for and open I/O API output file(s)...
                CALL GETENV( 'SRGPRO_PATH', NAMBUF )
                WRITE( NAMBUFT, '( 2A )' ) TRIM( NAMBUF ), SRGFNAM( I )

                IF( TGTSRG .NE. SRGFCOD( I ) ) CYCLE

                IF( NAMBUFT .NE. TSRGFNAM  ) THEN
                    CALL OPEN_SRGFILE

C.....................  Reading surrogate files
                    DO JJ = 1, NLINES

                       READ ( GDEV, 93000, END=113, IOSTAT=IOS )LINE

                       IF ( BLKORCMT( LINE ) ) CYCLE

                       CALL PARSLINE( LINE, MXSEG, SEGMENT )
                       SSC    = STR2INT ( SEGMENT( 1 ) )
                       TCOL   = STR2INT ( SEGMENT( 3 ) )
                       TROW   = STR2INT ( SEGMENT( 4 ) )

C........................  Check the value of the column number
                       IF( TCOL .LT. 0 .OR.  TCOL .GT. SRGNCOLS .OR.
     &                   ( TROW .EQ. 0 .AND. TCOL .NE. 0 ) ) THEN
                           WFLAG = .TRUE.
                       END IF
              
C........................  Check the value of the row number
                       IF( TROW .LT. 0 .OR.  TROW .GT. SRGNROWS .OR.
     &                   ( TCOL .EQ. 0 .AND. TROW .NE. 0 ) ) THEN
                           WFLAG = .TRUE.
                           CALL M3MESG( MESG )                    

C........................  Special treatment for cell (0,0) (skip for now)
                       ELSE IF( TROW .EQ. 0 .AND. TCOL. EQ. 0 ) THEN
                           CYCLE

                       END IF

C........................  Adjust column and row for subgrid
                       TCOL = TCOL - XOFF
                       TROW = TROW - YOFF

C........................  Skip entry after subgrid adjustment
                       IF( TCOL .LE. 0 .OR. TCOL .GT. NCOLS .OR.
     &                     TROW .LE. 0 .OR. TROW .GT. NROWS ) CYCLE

C........................  Skip entry if rows and columns are out of range
                       IF( WFLAG ) CYCLE

C........................  Skip entry if SSC is not in the assigned SRGLIST by source
                       IF( SSC .EQ. TGTSRG ) THEN
                           NT = NT + 1
                           TMPLINE( NT ) = LINE
                       END IF
                       
113                 END DO

                    TSRGFNAM = NAMBUFT    ! store a previous surrogate file name buffer
                    CLOSE( GDEV )

                END IF      ! skip if surrogate file has the same srg file

            END DO       ! loop over all surrogate files in SRGDESC file

            CALL RDSRG4GRD( NT, TMPLINE, CFLAG ) ! populating surrogates

            DEALLOCATE( TMPLINE )
            
            LFIP = ' '
            NWARN = 0

C.............  Loop over sources per each assigned surrogate
            DO S = 1, NSRC

                CFIP = CIFIP( S )
                SSC  = ASRGID( S )
                IF( CATEGORY .EQ. 'AREA' ) CELLSRC = CELLID( S )
                IF( CATEGORY .EQ. 'MOBILE' ) CLNK = CLINK( S )

                IF( SSC .NE. TGTSRG ) CYCLE
            
C.................  Determine if x/y location is available
                XYSET = .FALSE.
                IF( LFLAG ) XYSET = ( XLOCA( S ) .GT. AMISS3 )
            
C.................  If cell-specific source...
                IF ( CELLSRC .GT. 0 ) THEN
                    NCEL = 1
                    ACEL( 1 ) = CELLID( S )
                    AFAC( 1 ) = 1.
            
C.................  Check if source has been converted to point src
                ELSE IF( XYSET ) THEN
            
C....................  If source is in the domain....
                    XX = XLOCA( S )
                    YY = YLOCA( S )
                    IF( INGRID( XX, YY, NCOLS, NROWS, COL, ROW  ) ) THEN
            
C........................  Set as 1 cell and get the cell number
                        NCEL = 1
                        ACEL( 1 ) = ( ROW-1 ) * NCOLS + COL
                        AFAC( 1 ) = 1.
            
C.....................  Otherwise, skip this source because it's outside the grid
                    ELSE
                        NCEL = 0

                    END IF
            
C................  If area/non-link source...
                ELSE IF( CLNK .EQ. ' ' ) THEN
            
C.....................  Retrieve the index to the surrogates cy/st/co list
                    ISIDX = 1
                    F     = FINDC( CFIP, NSRGFIPS, SRGFIPS )

C.....................  Retrieve the cell intersection info from the
C                       surrogates tables from MODSURG
                    IF ( F .GT. 0 ) THEN
            
                        NCEL = NCELLS( F )
                        ACEL( 1:NCEL ) = FIPCELL( 1:NCEL, F )   ! arrays
            
                        DO K = 1, NCEL
                            CALL SETFRAC( S, ISIDX, TGTSRG, K, F, 1,  
     &                                  .FALSE.,' ', DEFSRGID, SRGFLAG,
     &                                   ID1, ID2, AFAC( K ), CFLAG )

C.........................  Re-assigning org assigned srg to default fallback srg
                            IF( ID2 .EQ. DEFSRGID .AND. SRGFLAG ) THEN
                                ASRGID( S ) = DEFSRGID
                                CYCLE
                            END IF

                        END DO

C.....................  Otherwise, skip this source because it's outside the grid
                    ELSE
                        IF( SRGFLAG ) THEN
                            ASRGID( S ) = DEFSRGID
                            IF( II < NSRGS ) CYCLE
                        END IF

C..........................  Critical warning message about zeroing emission 
C                            due to no surrogates for this co/st/cy
                        IF( LFIP .NE. CFIP ) THEN
                            WRITE( MESG,94010 ) 'WARNING: Causing ' //
     &                         'zeroing emissions due to missing '//
     &                         'surrogate', TGTSRG, ' for co/st/cy :: '//
     &                         CFIP
                            NWARN = NWARN + 1
                            IF( NWARN < 100 ) CALL M3MESG( MESG )
                            LFIP = CFIP
                        END IF

                        CYCLE

                    END IF


C.................  If link source, determine the number of cells for this source
                ELSE
                    CALL LNK2GRD( NGRID, XLOC1( S ), YLOC1( S ), 
     &                            XLOC2( S ), YLOC2( S ), NCEL, ACEL,
     &                            AFAC, ALEN, EFLAG)

C.....................  Make sure that there was enough storage 
                    IF( EFLAG ) THEN
                        WRITE( MESG,94010 )
     &                      'INTERNAL ERROR: Overflow for source', S
                        CALL M3MSG2( MESG )
                        CYCLE
                    END IF

                END IF
C           
C.................  Loop through the cells for this source and increment the number
C                   of sources per cell. 
                CCNT = 0
                DO N = 1, NCEL
            
                    IF( AFAC( N ) .GT. 0.0 ) THEN
                        C = ACEL( N )
                        NX( C ) = NX( C ) + 1
                        CCNT = CCNT + 1
                    END IF

                END DO    ! End loop on cells for this source
            
C.................  Update the maximum number of cells per source
                IF( CCNT .GT. MXCSRC ) MXCSRC = CCNT
            
C.................  Update the maximum number of cells per county or link
                IF( NCEL .GT. MXCCL ) MXCCL = NCEL
            
            END DO        ! End loop on sources

        END DO        ! End loop on assigned surrogates
        
C.........  Abort if error
        IF( EFLAG ) THEN
            MESG = 'Problem determining memory for gridding matrix.'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

C.........  Determine maximum number of sources per cell
C.........  And determine the total number of source-cell intersections
        MXSCEL = NX( 1 )
        NMATX  = NX( 1 )          
        DO C = 2, NGRID
            
            J = NX( C )
            IF( J .GT. MXSCEL ) THEN
                MXSCEL = J
            END IF

            NMATX = NMATX + J

        END DO        ! End loop on cells

C..........................................................
C.........  Estimating ungridding matrices sizes ..........
C..........................................................

C...... Loop over sources per each assigned surrogate
        DO S = 1, NSRC

            CFIP = CIFIP( S )
            SSC  = ASRGID( S )
            IF( CATEGORY .EQ. 'AREA' ) CELLSRC = CELLID( S )
            IF( CATEGORY .EQ. 'MOBILE' ) CLNK = CLINK( S )

C.........  Determine if x/y location is available
            XYSET = .FALSE.
            IF( LFLAG ) XYSET = ( XLOCA( S ) .GT. AMISS3 )

C.........  If cell-specific source...
            IF ( CELLSRC .GT. 0 ) THEN
                NCEL = 1
                ACEL( 1 ) = CELLID( S )
                AFAC( 1 ) = 1.
       
C.........  Check if source has been converted to point src
            ELSE IF( XYSET ) THEN
       
C..............  If source is in the domain....
                XX = XLOCA( S )
                YY = YLOCA( S )
                IF( INGRID( XX, YY, NCOLS, NROWS, COL, ROW  ) ) THEN

C....................  Set as 1 cell and get the cell number
                    NCEL = 1
                    ACEL( 1 ) = ( ROW-1 ) * NCOLS + COL
                    AFAC( 1 ) = 1.
       
C.................  Otherwise, skip this source because it's outside the grid
                ELSE
                    NCEL = 0

                END IF
       
C............  If area/non-link source...
            ELSE IF( CLNK .EQ. ' ' ) THEN
       
C.................  Retrieve the index to the surrogates cy/st/co list
                ISIDX = 1
                F     = FINDC( CFIP, NSRGFIPS, SRGFIPS )

C.............  Retrieve the cell intersection info from the
C               surrogates tables from MODSURG
                IF ( F .GT. 0 ) THEN
       
                    NCEL = NCELLS( F )

C.............  Otherwise, skip this source because it's outside the grid
                ELSE
                    NCEL = 1 

                END IF

C.............  If link source, determine the number of cells for this source
            ELSE
                CALL LNK2GRD( NGRID, XLOC1( S ), YLOC1( S ), 
     &                       XLOC2( S ), YLOC2( S ), NCEL, ACEL,
     &                       AFAC, ALEN, EFLAG)

C.............  Make sure that there was enough storage 
                IF( EFLAG ) THEN
                    WRITE( MESG,94010 )
     &                  'INTERNAL ERROR: Overflow for source', S
                    CALL M3MSG2( MESG )
                    CYCLE
                END IF

            END IF
C      
C.............  Count all county/cell intersections for all sources.  This
C               is needed for ungridding matrix.
            DO N = 1, NCEL

                NMATXU = NMATXU + 1

            END DO    ! End loop on cells for this source

        END DO        ! End loop on sources

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats............ 94xxx

93000   FORMAT( A )

94010   FORMAT( 10( A, :, I8, :, 1X ) )

C******************  INTERNAL SUBPROGRAMS  *****************************
        CONTAINS

C.........  This internal subprogram opens individual surrogate file 

            SUBROUTINE OPEN_SRGFILE
C----------------------------------------------------------------------
                 
C.........  Set logical file name
            IF( .NOT. SETENVVAR( 'SRGPRO_PATH', NAMBUFT )) THEN
                MESG = 'Could not set logical file ' //
     &                 'name of file ' // TRIM( NAMBUFT )
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF

C.........  Get the number of lines in the surrogate description file desription file
            GDEV = GETEFILE( 'SRGPRO_PATH',.TRUE., .TRUE., PROGNAME )

            IF( GDEV .LT. 0 ) THEN
                MESG = 'Could not open input surrogate file' // 
     &                 TRIM( NAMBUFT )
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF

            REWIND( GDEV )

            NLINES = GETFLINE( GDEV, 'Reading surrogate files' )
            
            IF( .NOT. SETENVVAR( 'SRGPRO_PATH', NAMBUF )) THEN
                MESG = 'Could not set logical file ' //
     &                 'name of file ' // TRIM( NAMBUF )
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF

            RETURN

C------------------- SUBPROGRAM FORMAT STATEMENTS ----------------------

            END SUBROUTINE OPEN_SRGFILE

        END SUBROUTINE SIZGMAT
