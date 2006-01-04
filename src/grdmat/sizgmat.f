
        SUBROUTINE SIZGMAT( CATEGORY, NSRC, VFLAG, OSDEF, FSGFLAG,
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
C      Created by M. Houyoux 5/99
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

C...........   MODULES for public variables
C...........   This module is the source inventory arrays
        USE MODSOURC, ONLY: XLOCA, YLOCA, IFIP, CELLID, CLINK,
     &                      XLOC1, YLOC1, XLOC2, YLOC2

C...........   This module contains the cross-reference tables
        USE MODXREF, ONLY: ASRGID

C.........  This module contains the global variables for the 3-d grid
        USE MODGRID, ONLY: NGRID, NCOLS, NROWS

C...........   This module contains the gridding surrogates tables
        USE MODSURG, ONLY: NCELLS, FIPCELL, NSRGS, SRGLIST, NSRGFIPS,
     &                     SRGFIPS, TMPLINE, NTSRGDSC, SRGFNAM, SRGFCOD,
     &                     NTLINES, MXCFIP
        IMPLICIT NONE

C...........   INCLUDES:
        
        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters

C...........   EXTERNAL FUNCTIONS and their descriptions:
        CHARACTER(2)    CRLF
        INTEGER         ENVINT
        INTEGER         FIND1
        LOGICAL         INGRID
        LOGICAL         BLKORCMT
        LOGICAL         SETENVVAR
        INTEGER         STR2INT
        INTEGER         GETFLINE
        INTEGER         PROMPTFFILE

        EXTERNAL        CRLF, FIND1, INGRID, BLKORCMT, SETENVVAR,
     &                  STR2INT, GETFLINE, PROMPTFFILE, ENVINT


C...........   SUBROUTINE ARGUMENTS
        CHARACTER(*), INTENT (IN) :: CATEGORY  ! source category
        INTEGER     , INTENT (IN) :: NSRC      ! local number of sources
        LOGICAL     , INTENT (IN) :: VFLAG     ! true: using variable grid
        INTEGER     , INTENT (IN) :: OSDEF     ! original index of default surrogate ID
        LOGICAL     , INTENT (IN) :: FSGFLAG   ! true: using default fallback surrogate
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
        INTEGER         DEFSRGID         !  default surrogate ID
        INTEGER         FIP              ! country/state/county code
        INTEGER         GDEV             !  for surrogate coeff file
        INTEGER         ID1, ID2         ! primary and 2ndary surg codes
        INTEGER         ISIDX            ! tmp surrogate ID code index
        INTEGER         IOS              ! i/o status
        INTEGER         IREC             ! Record counter
        INTEGER         ISDEF            !  default surrogate ID code index
        INTEGER         NCEL             ! tmp number of cells
        INTEGER         ROW              ! tmp row
        INTEGER         NTL              ! max no. of line buffers
        INTEGER         TGTSRG           ! target surrogates code
        INTEGER         SSC              ! surrogates code
        INTEGER      :: NLINES = 0       ! number of lines in input file

        INTEGER, PARAMETER :: MXSEG = 5           ! # of potential line segments

C...........   Other arrays
        CHARACTER(20) SEGMENT( MXSEG )             ! Segments of parsed lines

        REAL            ALEN        ! link length

        LOGICAL      :: EFLAG = .FALSE. ! true: error flag
        LOGICAL      :: LFLAG = .FALSE. ! true: location data available
        LOGICAL      :: XYSET = .FALSE. ! true: X/Y available for src

        CHARACTER(300)      MESG             !  message buffer
        CHARACTER(80)       LINE             ! Read buffer for a line
        CHARACTER(196)      NAMBUF           !  surrogate file name buffer
        CHARACTER(256)      NAMBUFT          !  tmp surrogate file name buffer
        CHARACTER(256)      TSRGFNAM         !  tmp surrogate file name buffer

        CHARACTER(LNKLEN3) :: CLNK = ' '   ! tmp link ID

        CHARACTER(16) :: PROGNAME = 'SIZGMAT' ! program name

C***********************************************************************
C   begin body of subroutine SIZGMAT

C.....  Print status message
        MESG = 'Computing gridding matrix size...'
        CALL M3MSG2( MESG )

C.....  Determine default surrogate number from the environment
C.....  Default surrogate code 8 is population
        DEFSRGID = ENVINT( 'SMK_DEFAULT_SRGID', 'Default surrogate',
     &                      8, IOS )

C.....  Set flag to indicate that XLOCA/YLOCA are available
        LFLAG = ALLOCATED( XLOCA )

C.....  Initialize the count of sources per cell
        NX = 0   ! array

C.....  Loop through sources
        MXCSRC  = 0
        MXCCL   = 0
        NMATXU  = 0
        CELLSRC = 0

C.....  Allocate memory for indices to surrogates tables for each source
        ALLOCATE( NTLINES( NSRGS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'NTLINES', PROGNAME )
 
C.....  Count total line buffers to define memory size
        DO II = 1, NSRGS + 1  ! loop through only the surrogate code assigned by sources

            NTL = 0
            TSRGFNAM = ' '

            ISDEF  = FIND1( DEFSRGID, NSRGS, SRGLIST )

C.........  default fallback surrogate will run at last after
C           re-assigned zero fraction surrogate
            IF( ISDEF .EQ. 1 ) THEN

                TGTSRG = SRGLIST( II )
                KK = II
                
                IF( II .EQ. NSRGS + 1 ) THEN
                    TGTSRG = DEFSRGID
                    KK = ISDEF
                END IF
                
            END IF

            IF( ISDEF .GT. 1 ) THEN
                
                IF( II .EQ. 1 ) THEN
                    TGTSRG = SRGLIST( ISDEF )
                    KK = ISDEF
                ELSE IF( II .GT. 1 .AND. II .LE. ISDEF ) THEN
                    TGTSRG = SRGLIST( II - 1 )
                    KK = II - 1
                ELSE IF( II .GT. ISDEF .AND. II .LT. NSRGS + 1 ) THEN
                    TGTSRG = SRGLIST( II )
                    KK = II
                ELSE IF( II .EQ. NSRGS + 1 ) THEN
                    TGTSRG = DEFSRGID
                    KK = ISDEF
                END IF
            
            END IF
                

C.........   Count total line buffers to define memory size
            DO I = 1, NTSRGDSC  ! Open all surrogate files using the same srg code
       
C.............  Prompt for and open I/O API output file(s)...
                CALL GETENV( "SRG_PATH", NAMBUF )
                WRITE( NAMBUFT, '( 2A )' ) TRIM( NAMBUF ), SRGFNAM( I )

                IF( TGTSRG .NE. SRGFCOD( I ) ) CYCLE
       
                IF( NAMBUFT .NE. TSRGFNAM  ) THEN
                    CALL OPEN_SRGFILE

                    IREC = 0
C.................  Reading surrogate files
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

C.................  Skip entry if SSC is not in the assigned SRGLIST by source
                       IF( SSC .EQ. TGTSRG ) NTL = NTL + 1

111                 END DO

                    TSRGFNAM = NAMBUFT    ! store a previous surrogate file name buffer
                    CLOSE( GDEV )
       
                END IF      ! skip if surrogate file has the same srg file
             
            END DO       ! loop over all surrogate files in SRGDESC file
            
C.........  Store no of line buffers of each surrogate
            IF( NTL .EQ. 0 ) CYCLE
            NTLINES( KK ) = NTL

C.........  skip default surrogate if FSGFLAG is false
            IF( .NOT. FSGFLAG .AND. OSDEF .LE. 0 ) NTLINES( ISDEF ) = 0

C.........  Allocate memory for indices to surrogates tables for each source
            ALLOCATE( TMPLINE( NTL ), STAT=IOS )
            CALL CHECKMEM( IOS, 'TMPLINE', PROGNAME )
        
C.........  If surrogates are needed, read and store the gridding surrogates, 
C           allocate memory for the surrogate assignments, and assign
C           surrogates to each source.
            NT = 0
            TMPLINE = ' '
            TSRGFNAM = ' '
      
            DO I = 1, NTSRGDSC  ! Open all surrogate files using the same srg code
       
C.............  Prompt for and open I/O API output file(s)...
                CALL GETENV( "SRG_PATH", NAMBUF )
                WRITE( NAMBUFT, '( 2A )' ) TRIM( NAMBUF ), SRGFNAM( I )
                
                IF( TGTSRG .NE. SRGFCOD( I ) ) CYCLE
       
                IF( NAMBUFT .NE. TSRGFNAM  ) THEN
                    CALL OPEN_SRGFILE

                    IREC = 0
C.................  Reading surrogate files
                    DO JJ = 1, NLINES

                       READ ( GDEV, 93000, END=113, IOSTAT=IOS ) LINE
                       IREC = IREC + 1

                       IF ( BLKORCMT( LINE ) ) CYCLE

                       CALL PARSLINE( LINE, MXSEG, SEGMENT )
                       SSC    = STR2INT ( SEGMENT( 1 ) )

C.................  Skip entry if SSC is not in the assigned SRGLIST by source
                       IF( SSC .EQ. TGTSRG ) THEN
                           NT = NT + 1
                           TMPLINE( NT ) = LINE
                       END IF
                       
113                 END DO

                    TSRGFNAM = NAMBUFT    ! store a previous surrogate file name buffer
                    CLOSE( GDEV )

                 END IF      ! skip if surrogate file has the same srg file

            END DO       ! loop over all surrogate files in SRGDESC file

            CALL SIZRDSRG( NT, VFLAG )

            DEALLOCATE( TMPLINE )
            
            IF( II == 1 ) CYCLE

C.........  Loop over sources per each assigned surrogate
            DO S = 1, NSRC

                FIP = IFIP( S )
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
                    IF( INGRID( XLOCA( S ), YLOCA( S ), 
     &                          NCOLS, NROWS, COL, ROW  ) ) THEN
            
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
                    F     = FIND1( FIP, NSRGFIPS, SRGFIPS )

C.....................  Retrieve the cell intersection info from the
C                       surrogates tables from MODSURG
                    IF ( F .GT. 0 ) THEN
            
                        NCEL = NCELLS( F )
                        ACEL( 1:NCEL ) = FIPCELL( 1:NCEL, F )      ! arrays
            
                        DO K = 1, NCEL
                            CALL SETFRAC( S, ISIDX, TGTSRG, K, F, 1,  
     &                               .FALSE.,' ', ID1, ID2, AFAC( K ) )
                        END DO
                        
                        
C.....................  Otherwise, skip this source because it's outside the grid
                    ELSE
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
            
                    IF( AFAC( N ) .GT. 0. ) THEN
                        C = ACEL( N )
                        NX( C ) = NX( C ) + 1
                        CCNT = CCNT + 1
                    END IF

C...................  Count all county/cell intersections for all sources.  This
C                     is needed for ungridding matrix.
                    NMATXU = NMATXU + 1

                END DO    ! End loop on cells for this source
            
C.................  Update the maximum number of cells per source
                IF( CCNT .GT. MXCSRC ) MXCSRC = CCNT
            
C.................  Update the maximum number of cells per county or link
                IF ( NCEL .GT. MXCCL ) MXCCL = NCEL
            
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

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats............ 94xxx

93000   FORMAT( A )

94010   FORMAT( 10( A, :, I8, :, 1X ) )

C******************  INTERNAL SUBPROGRAMS  *****************************
        CONTAINS

C.........  This internal subprogram open individual surrogate file 

            SUBROUTINE OPEN_SRGFILE
C----------------------------------------------------------------------
                 
C.........  Set logical file name
            IF( .NOT. SETENVVAR( "SRG_PATH", NAMBUFT )) THEN
                MESG = 'Could not set logical file ' //
     &                 'name of file ' // TRIM( NAMBUFT )
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF

C.........  Get the number of lines in the surrogate description file desription file
            GDEV = PROMPTFFILE( 'Reading surrogate files..',
     &             .TRUE., .TRUE., 'SRG_PATH', PROGNAME )
     
            REWIND( GDEV )

            NLINES = GETFLINE( GDEV, 'Reading srg files' )
            
            IF( .NOT. SETENVVAR( "SRG_PATH", NAMBUF )) THEN
                MESG = 'Could not set logical file ' //
     &                 'name of file ' // TRIM( NAMBUF )
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF

            RETURN

C------------------- SUBPROGRAM FORMAT STATEMENTS ----------------------

            END SUBROUTINE OPEN_SRGFILE

        END SUBROUTINE SIZGMAT
