
        SUBROUTINE GENMGMAT( ENAME, GNAME, UNAME, FDEV, MXSCEL, MXCSRC, 
     &                       MXCCL, NMSRC, NMATX, NMATXU, UFLAG, NX, 
     &                       IX, CX, NU, IU, CU, NCOEF, CMAX, CMIN, 
     &                       NCOEFU )

C***********************************************************************
C  subroutine body starts at line 
C
C  DESCRIPTION:
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
C COPYRIGHT (C) 2002, MCNC Environmental Modeling Center
C All Rights Reserved
C
C See file COPYRIGHT for conditions of use.
C
C Environmental Modeling Center
C MCNC
C P.O. Box 12889
C Research Triangle Park, NC  27709-2889
C
C smoke@emc.mcnc.org
C
C Pathname: $Source$
C Last updated: $Date$ 
C
C***************************************************************************

C...........   MODULES for public variables
C...........   This module is the source inventory arrays
        USE MODSOURC

C...........   This module contains the cross-reference tables
        USE MODXREF

C.........  This module contains the global variables for the 3-d grid
        USE MODGRID

C...........   This module contains the gridding surrogates tables
        USE MODSURG

C.........  This module contains the lists of unique source characteristics
        USE MODLISTS

C.........  This module contains the information about the source category
        USE MODINFO

        IMPLICIT NONE

C...........   INCLUDES

        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
        INCLUDE 'FDESC3.EXT'    !  I/O API file description data structures

C...........   EXTERNAL FUNCTIONS and their descriptions:
        CHARACTER*2     CRLF
        INTEGER         ENVINT
        INTEGER         FIND1
        INTEGER         FINDC
        LOGICAL         INGRID
        LOGICAL         DSCM3GRD

        EXTERNAL        CRLF, ENVINT, FIND1, FINDC, INGRID, DSCM3GRD

C...........   SUBROUTINE ARGUMENTS
        CHARACTER(*), INTENT (IN) :: ENAME         ! inventory file name
        CHARACTER(*), INTENT (IN) :: GNAME         ! gridding mtx logical name
        CHARACTER(*), INTENT (IN) :: UNAME         ! ungridding mtx logical name
        INTEGER     , INTENT (IN) :: FDEV          ! surg codes report file
        INTEGER     , INTENT (IN) :: MXSCEL        ! max sources per cell
        INTEGER     , INTENT (IN) :: MXCSRC        ! max cells per source
        INTEGER     , INTENT (IN) :: MXCCL         ! max cells per county or link
        INTEGER     , INTENT (IN) :: NMSRC         ! no. mobile sources
        INTEGER     , INTENT (IN) :: NMATX         ! no. source-cell intersects
        INTEGER     , INTENT (IN) :: NMATXU        ! county-cell intrscts for all sources
        LOGICAL     , INTENT (IN) :: UFLAG         ! true: create gridding matrix
        INTEGER     , INTENT(OUT) :: NX  ( NGRID ) ! no. srcs per cell
        INTEGER     , INTENT(OUT) :: IX  ( NMATX ) ! src IDs 
        REAL        , INTENT(OUT) :: CX  ( NMATX ) ! gridding coefficients
        INTEGER     , INTENT(OUT) :: NU  ( NMSRC ) ! no. cells per source
        INTEGER     , INTENT(OUT) :: IU  ( NMATXU ) ! cell numbers
        REAL        , INTENT(OUT) :: CU  ( NMATXU ) ! ungridding coefficients
        INTEGER     , INTENT(OUT) :: NCOEF         ! no. of gridding coeffs
        INTEGER     , INTENT(OUT) :: CMAX          ! max no. of sources per cell
        INTEGER     , INTENT(OUT) :: CMIN          ! min no. of sources per cell
        INTEGER     , INTENT(OUT) :: NCOEFU        ! no. of ungridding coeffs

C...........   Local allocatable arrays...

C...........   Scratch Gridding Matrix (subscripted by source-within-cell, cell)

        INTEGER, ALLOCATABLE :: IDXSRT ( : )! sorting index
        INTEGER, ALLOCATABLE :: IDXSRT2( : )! 2nd sorting index
        INTEGER, ALLOCATABLE :: IC     ( : )! cell IDs
        INTEGER, ALLOCATABLE :: IS     ( : )! source IDs
        INTEGER, ALLOCATABLE :: NCL    ( : )! position for county and/or link
        REAL   , ALLOCATABLE :: CS     ( : )! facs (no adjustments), for ungrid
        REAL   , ALLOCATABLE :: CSJ    ( : )! factors (with adjustments)

C...........   Scratch Ungridding Matrix information
        INTEGER, ALLOCATABLE :: CLIDX     ( : ) ! county or link index by source
        INTEGER, ALLOCATABLE :: CNT_CL    ( : ) ! cell count per county or link
        INTEGER, ALLOCATABLE :: VMT_CELL( :,: ) ! cell numbers for county or link
        REAL   , ALLOCATABLE :: VMT_FRAC( :,: ) ! VMT fraction for cell/county or cell/link
        REAL   , ALLOCATABLE :: VMT_CL_INV( : ) ! inverse of VMT by county or link
        INTEGER, ALLOCATABLE :: vmt_label ( : ) ! FIPS codes 

C...........   Temporary array for flagging sources that are outside the
C              domain
        LOGICAL, ALLOCATABLE :: INDOMAIN( : ) ! true: src is in the domain

C...........   Arrays for links intersecting with cells
C...........   Note that the NGRID dimension could conceivably be too small if
C              a link winds through the whole domain, but this is a case that
C              is not worth going to extra trouble for.
        INTEGER, ALLOCATABLE :: ACEL( : )    ! number of cell intersctns per src
        REAL   , ALLOCATABLE :: AFAC( : )    ! fraction of link in cell

        INTEGER, ALLOCATABLE :: FIPNOSRG( : )  ! cy/st/co codes w/o surrogates
        CHARACTER(LEN=SRCLEN3), ALLOCATABLE :: LKOGRD( : ) ! link srcs outside
                                                           !    the grid
C...........   Temporary arrays for storing surrogate codes to use
        INTEGER, ALLOCATABLE :: SURGID1( : ) ! primary surrogate code
        INTEGER, ALLOCATABLE :: SURGID2( : ) ! secondary surrogate code

C...........   Other local variables

        INTEGER         C, CNT, F, I, IDX, J, K, N, S !  indices and counters
        INTEGER         LC, LN              !  "previous" values of these

        INTEGER         CNTMAX   !  counter for storing correct max dimensions
        INTEGER         COL      !  tmp column
        INTEGER         FIP      !  tmp country/state/county code
        INTEGER         ID1, ID2 !  tmp primary and secondary surg codes
        INTEGER         IOS      !  i/o status
        INTEGER         ISIDX    !  surrogate ID code index
        INTEGER         L2       !  string length
        INTEGER         LFIP     !  cy/st/co code from previous iteration
        INTEGER         LNKEND   !  width of sourc info to end of link ID
        INTEGER         NCEL     !  tmp number of cells 
        INTEGER         NCOULNK  !  number of counties and links
        INTEGER         NLKOGRD  !  no. of link sources outside the grid
        INTEGER         NNOSRG   !  no. of cy/st/co codes with no surrogates
        INTEGER         ROW      ! tmp row
        INTEGER         RWT      !  tmp roadway type

        REAL            ADJ, ADJC   ! tmp adjustment factors
        REAL            ALEN        ! link length
        REAL            FRAC        ! tmp surrogate fraction
        REAL            XBEG, YBEG  ! tmp X and Y link start coordinates
        REAL            XEND, YEND  ! tmp X and Y link end   coordinates
        real            sum

        LOGICAL      :: EFLAG = .FALSE.  ! true: error detected
        LOGICAL      :: LFLAG = .FALSE.  ! true: location data available
        LOGICAL      :: XYSET = .FALSE. ! true: X/Y available for src

        CHARACTER*16    COORUNIT  !  coordinate system projection units
        CHARACTER*80    GDESC     !  grid description
        CHARACTER*256   BUFFER    !  source fields buffer
        CHARACTER*256   MESG      !  message buffer 

        CHARACTER(LEN=LNKLEN3)    CLNK   ! tmp link ID
        CHARACTER(LEN=LNKLEN3)    LLNK   ! previous link ID
        CHARACTER(LEN=SRCLEN3)    CSRC   ! tmp source chars string
        CHARACTER(LEN=SRCLEN3)    CSRC2  ! tmp truncated source chars string

        CHARACTER*16 :: PROGNAME = 'GENMGMAT' ! program name

C***********************************************************************
C   begin body of subroutine GENMGMAT

C.........  Get grid name from the environment and read grid parameters
        IF( .NOT. DSCM3GRD( GRDNM, GDESC, COORD, GDTYP3D, COORUNIT, 
     &                      P_ALP3D, P_BET3D, P_GAM3D, XCENT3D, YCENT3D,
     &                      XORIG3D, YORIG3D, XCELL3D, YCELL3D,
     &                      NCOLS3D, NROWS3D, NTHIK3D ) ) THEN

            MESG = 'Could not get Models-3 grid description'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

C.........  Store grid parameters for later processing
        ELSE

            IF( NCOLS3D * NROWS3D .NE. NGRID ) THEN
 
                MESG = 'INTERNAL ERROR: Number of cells in "' //
     &                 PROGNAME( 1:16 ) // '" are inconsistent '//
     &                 'with calling program.'
                CALL M3MSG2( MESG )
                CALL M3EXIT( PROGNAME, 0, 0, ' ', 2 )

            END IF

        END IF

C.........  If ungridding matrix is needed, allocate memory for VMT, 
C           and read the VMT data.
        IF ( UFLAG ) THEN
            ALLOCATE( VMT( NSRC ), STAT=IOS )
            CALL CHECKMEM( IOS, 'VMT', PROGNAME )
            CALL RDMAPPOL( NSRC, 1, 1, 'VMT', VMT )
        END IF

C.........  Initialize number of sources per cell counter
        NX = 0   ! Array

C.........  Allocate memory for temporary gridding matrix and other
        ALLOCATE( IDXSRT( NMATX ), STAT=IOS )
        CALL CHECKMEM( IOS, 'IDXSRT', PROGNAME )
        ALLOCATE( IC( NMATX ), STAT=IOS )
        CALL CHECKMEM( IOS, 'IC', PROGNAME )
        ALLOCATE( IS( NMATX ), STAT=IOS )
        CALL CHECKMEM( IOS, 'IS', PROGNAME )
        ALLOCATE( NCL( NMATX ), STAT=IOS )
        CALL CHECKMEM( IOS, 'NCL', PROGNAME )
        ALLOCATE( CS( NMATX ), STAT=IOS )
        CALL CHECKMEM( IOS, 'CS', PROGNAME )
        ALLOCATE( CSJ( NMATX ), STAT=IOS )
        CALL CHECKMEM( IOS, 'CSJ', PROGNAME )
        ALLOCATE( INDOMAIN( NMSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'INDOMAIN', PROGNAME )
        ALLOCATE( ACEL( NGRID ), STAT=IOS )
        CALL CHECKMEM( IOS, 'ACEL', PROGNAME )
        ALLOCATE( AFAC( NGRID ), STAT=IOS )
        CALL CHECKMEM( IOS, 'AFAC', PROGNAME )
        ALLOCATE( FIPNOSRG( NINVIFIP ), STAT=IOS )
        CALL CHECKMEM( IOS, 'FIPNOSRG', PROGNAME )
        ALLOCATE( LKOGRD( NMSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'LKOGRD', PROGNAME )
        ALLOCATE( SURGID1( NSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'SURGID1', PROGNAME )
        ALLOCATE( SURGID2( NSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'SURGID2', PROGNAME )
        IC  = 0
        IS  = 0
        NCL = 0
        CS  = 0.
        CSJ = 0.
        FIPNOSRG = 0
        SURGID1 = 0   ! array
        SURGID2 = 0   ! array

C.........  Set flag to indicate that XLOCA/YLOCA are available
        LFLAG = ALLOCATED( XLOCA )

C.......   Compute gridding matrix:
C.......       First case:   explicit link (ILINK > 0
C.......       Second case:  some LNKDEF entry applies
C.......       Third case:   FIP+roadtype cross-reference match
C.......       Fourth case:  state+roadtype cross-reference match
C.......       fifth case:   roadtype cross-reference match
C.......       sixth case:   fallback default

        MESG = 'Computing gridding matrix and statistics...'
        CALL M3MSG2( MESG )

        J       = 0
        LFIP    = 0
        LLNK    = ' '
        LNKEND  = VIDPOS3 - 1
        NNOSRG  = 0
        NLKOGRD = 0
        ADJ     = 1.    !  temporary
        ADJC    = 1.    !  temporary
        N       = 0

        DO S = 1, NSRC
            
            FIP  = IFIP  ( S )
            RWT  = IRCLAS( S )
            CLNK = CLINK ( S )
            CSRC = CSOURC( S )
            CSRC2= CSRC( 1:VIDPOS3-1 )  ! abridged for these purposes

C.............  Initialize sources as being in the domain
            INDOMAIN( S ) = .TRUE.

C.............  Count the number of county and links
            IF ( FIP .NE. LFIP .OR. CLNK .NE. LLNK ) N = N + 1            

C.............  Determine if x/y location is available
            XYSET = .FALSE.
            IF( LFLAG ) XYSET = ( XLOCA( S ) .GT. AMISS3 )

C.............  Special case for source has an x/y location
            IF( XYSET ) THEN

C................  If source is in the domain, get cell number and store
                IF( INGRID( XLOCA( S ), YLOCA( S ), 
     &                      NCOLS, NROWS, COL, ROW  ) ) THEN

                    C = ( ROW-1 ) * NCOLS + COL
                    NX( C ) = NX( C ) + 1
                    J = J + 1

C.....................  Check that the maximum number of sources per cell is ok
                    IF ( J .LE. NMATX ) THEN
                	IS ( J ) = S
                        IC ( J ) = C
                        NCL( J ) = N
                	CS ( J ) = 1.0
                    END IF

C................  Otherwise, mark source as being outside domain
                ELSE
                    INDOMAIN( S ) = .FALSE.
                END IF

                CYCLE           ! To head of loop over sources

            END IF   ! End if assigned point location or not

C.............  Find FIP/RoadWayType adjustment for FIP and RWT
c            ADJ = ADJMV( NADJ1, FIP, RWT, ADJFIP, ADJRWT, ADJFAC1 )

C.............  Process for link source...
            IF ( CLNK .NE. ' ' ) THEN

C.................  Convert coordinate system from UTM to required system

                XBEG = XLOC1 ( S )
                YBEG = YLOC1 ( S )
                XEND = XLOC2 ( S )
                YEND = YLOC2 ( S )

C.................  Compute the fractions of the link in each grid cell
                CALL LNK2GRD( NGRID, XBEG, YBEG, XEND, YEND,
     &                        NCEL, ACEL, AFAC, ALEN, EFLAG  )

C.................  Make sure that there was enough storage 
                IF ( EFLAG ) THEN
                    CALL FMTCSRC( CSRC, NCHARS, BUFFER, L2 )
                    WRITE( MESG,94010 )
     &                 'INTERNAL ERROR: Overflow for link to ' //
     &                 'grid conversion for:' // CRLF()// BLANK10//
     &                 BUFFER( 1:L2 )

                    CALL M3MESG( MESG )
                    CYCLE

C.................  If link is outside the grid, store its information and 
C                   go to next loop iteration
                ELSE IF( NCEL .EQ. 0 ) THEN

                    NLKOGRD = NLKOGRD + 1
                    LKOGRD( NLKOGRD ) = CSRC( 1:LNKEND )

                    CYCLE
                    
C.................  Write error if the link has no starting coordinates = ending
                ELSE IF( NCEL .EQ. -1 ) THEN

                    EFLAG = .TRUE.
                    CALL FMTCSRC( CSRC, NCHARS, BUFFER, L2 )
                    WRITE( MESG,94010 )
     &                  'Zero-length link for ' // BUFFER( 1:L2 )
                    CALL M3MESG( MESG )
                    CYCLE

                END IF

C.................  Loop through cells intersecting the current link
                DO K = 1, NCEL

                    C = ACEL( K )

C.....................  Find cell-based adjustment for cell C
c                    ADJC = ADJMV( NADJ2, C, 0, ADJCELL, ADJCELL, 
c     &                            ADJFAC2 )

C.....................  Warn when both types of adjustments are applied
c                    IF( ADJ .NE. 1. .AND. ADJC .NE. 1. ) THEN

c                        CALL FMTCSRC( CSRC, NCHARS, BUFFER, L2 )
c                        WRITE( MESG,94010 ) 'WARNING: Both FIP/'//
c     &                         'RoadWayType and cell-specific adjustments'
c     &                         // ' applied for ' // CRLF() // BLANK10//
c     &                         BUFFER( 1:L2 )
c                        CALL M3MSG2( MESG )

c                    END IF

C.....................  Make sure cell ID is in valid range.
                    IF( C .GT. NGRID .OR. C .LT. 0 ) THEN
                        EFLAG = .TRUE.
                        CALL FMTCSRC( CSRC, NCHARS, BUFFER, L2 )
                        WRITE( MESG,94010 ) 'Tried to assign link to'//
     &                         ' illegal cell =', C, 'for:' //
     &                         CRLF() // BLANK5 // BUFFER( 1:L2 )
                        CALL M3MESG( MESG )
                        CYCLE    ! to next cell in link

C.....................  Otherwise, increase the count of records and
C.....................  store the source count for this cell
                    ELSE
                        J = J + 1
                        NX( C ) = NX( C ) + 1

                    END IF

C.....................  Check that the maximum number of src-cell intersections
C                       is okay
                    IF ( J .LE. NMATX ) THEN

                        IDXSRT( J ) = J
                        IS    ( J ) = S
                        IC    ( J ) = C
                        NCL   ( J ) = N
                        CS    ( J ) = AFAC( K )
                        CSJ   ( J ) = ADJ * ADJC * AFAC( K )
                    END IF

                END DO  ! end loop over cells for this link

                CYCLE   ! to next source

            END IF      ! end of link-specific processing

C............. Start of non-link processing (will only get here if the
C              source is a non-link source)...

C............. Look for source in the grid-by-link table
c            S = FIND2( FIP, RWT, NUMLNK, LNKFIP, LNKRWT )  !  use LNKDEF factor
c            IF ( S .GT. 0 ) THEN
                
c                DO K = 1, LNKCNT( S )
c                    C = LNKCEL( J,S )
c                    J = NX( C ) + 1

C.....................  Find cell-based adjustment for cell C
c                    ADJC = ADJMV( NADJ2, C, 0, ADJCELL, ADJCELL, 
c     &                            ADJFAC2 )

C.....................  Warn when both types of adjustments are applied
c                    IF( ADJ .NE. 1. .AND. ADJC .NE. 1. ) THEN

c                        CALL FMTCSRC( CSRC, NCHARS, BUFFER, L2 )
c                        WRITE( MESG,94010 ) 'WARNING: Both FIP/'//
c     &                         'RoadWayType and cell-specific adjustments'
c     &                         // ' applied for ' // CRLF() // BLANK10//
c     &                         BUFFER( 1:L2 ) // ' Cell:', C
c                        CALL M3MSG2( MESG )

c                    END IF

C.....................  Check that the maximum number of sources per cell is ok
c                    IF ( J .LE. MXSCEL ) THEN
c                        IS ( J,C ) = S
c                        CS ( J,C ) = LNKFAC( K,S )
c                        CSJ( J,C ) = ADJ * ADJC * LNKFAC( K,S )
c                    END IF

C.....................  Store the count of sources for current cell
c                    NX( C )   = J

c                END DO   ! end of loop on cells for this source

c                CYCLE    ! to next source

c            END IF       ! end of grid-by-link processing

C.............  Process for non-link, non-grid-by-link sources...

C.............  Retrieve the indices to the surrogates tables
            ISIDX = SRGIDPOS( S )
            F     = SGFIPPOS( S )

C.............  Store the number and values of unfound cy/st/co codes
C.............  Keep track of sources that are outside the domain
            IF ( F .LE. 0 ) THEN

                IF( FIP .NE. LFIP ) THEN
                    NNOSRG = NNOSRG + 1
                    FIPNOSRG( NNOSRG ) = FIP
                    LFIP = FIP
                END IF
              
                INDOMAIN( S ) = .FALSE.
                CYCLE   ! To next source

            END IF

C.............  Loop through all of the cells intersecting this co/st/cy code. 
            DO K = 1, NCELLS( F )
            
                C = FIPCELL( K,F )   ! Retrieve cell number

C.................  Set the surrogate fraction
                CALL SETFRAC( S, ISIDX, K, F, 2, INDOMAIN( S ), 
     &                        CSRC2, ID1, ID2, FRAC )

C.................  Store surg IDs for reporting
                SURGID1( S ) = ID1
                SURGID2( S ) = ID2

                IF( FRAC .GT. 0 ) THEN

                    J = J + 1              ! increment no. src-cell intersection
                    NX( C ) = NX( C ) + 1  ! increment no. srcs per cell

C.....................  Find cell-based adjustment for cell C
c                    ADJC = ADJMV( NADJ2, C, 0, ADJCELL, ADJCELL, ADJFAC2 )

C.....................  Warn when both types of adjustments are applied
c                    IF( ADJ .NE. 1. .AND. ADJC .NE. 1. ) THEN

c                        CALL FMTCSRC( CSRC, NCHARS, BUFFER, L2 )
c                        WRITE( MESG,94010 ) 'WARNING: Both FIP/'//
c     &                         'RoadWayType and cell-specific adjustments'//
c     &                         ' applied for ' // CRLF() // BLANK10 //
c     &                         BUFFER( 1:L2 ) // ' Cell:', C
c                        CALL M3MSG2( MESG )

c                    END IF


C.....................  Check that the maximum number of src-cell intersections
C                       is okay
                    IF ( J .LE. NMATX ) THEN
                        IDXSRT( J ) = J
                        IS    ( J ) = S
                        IC    ( J ) = C
                        NCL   ( J ) = N
                        CS    ( J ) = FRAC
                        CSJ   ( J ) = ADJ * ADJC * FRAC
                    END IF

                END IF  ! if surrogate fraction > 0.

            END DO    !  end of loop on cells K for this FIP

            LFIP = FIP
            LLNK = CLNK

        END DO        !  end loop on sources S, computing gridding matrix.
        
C.........  Abort if error
        IF( EFLAG ) THEN
            MESG = 'Problem creating gridding matrix.'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

        IF( NLKOGRD .NE. 0 ) THEN

            WRITE( MESG,94010 ) 'NOTE: ', NLKOGRD, 
     &                          'link sources were outside grid.'
            CALL M3MSG2( MESG )

        END IF

        NCOEF = J
        NCOULNK = N

C.........  Abort if overflow occurred
        IF ( NCOEF .GT. NMATX ) THEN   
            
            WRITE( MESG,94010 )
     &        'INTERNAL ERROR: Gridding and ungridding matrices not ' //
     &        'written.' // CRLF() // BLANK10 //
     &        'Arrays would have overflowed.' 
     &        // CRLF() // BLANK10 // 
     &        'Current maximum cell-source intersections (NMATX) =', 
     &        NMATX, '.' // CRLF() // BLANK10 // 
     &        'Actual  maximum cell-source intersections         =', 
     &        NCOEF, '.'
            CALL M3MSG2( MESG )
            CALL M3EXIT( PROGNAME, 0, 0, ' ', 2 )

        END IF

C.........  Sort the scratch gridding matrix arrays to organize by cell and by
C           source.
        CALL SORTI2( NCOEF, IDXSRT, IC, IS )

C.........  Compress scratch gridding matrix into output representation
C.........  and compute statistics.  Already NX part of gridding matrix

        CMAX = 0
        CMIN = 99999999
        DO C = 1, NGRID
            J = NX( C )
            IF (      J .GT. CMAX ) THEN
                CMAX = J
            ELSE IF ( J .GT. 0 .AND. J .LT. CMIN ) THEN
                CMIN = J
            END IF
        END DO
 
        DO K = 1, NCOEF  ! Loop through cell-src intersections

            J = IDXSRT( K )                           
            IX( K ) = IS ( J )
            CX( K ) = CSJ( J )
                   
        END DO

C.........  Write gridding matrix
        MESG = 'Writing out GRIDDING MATRIX file...'
        CALL M3MSG2( MESG )

        IF( .NOT. WRITE3( GNAME, 'ALL', 0, 0, NX ) ) THEN
            MESG = 'Error writing GRIDDING MATRIX file.'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

C.........  Write output surrogates codes
        DO S = 1, NSRC
            WRITE( FDEV,93360 ) S, SURGID1( S ), SURGID2( S )
        END DO

C.........  Deallocate memory that is no longer needed
        DEALLOCATE( CSJ, IDXSRT )

C............................................................................
C..........   Generate ungridding matrix ....................................
C............................................................................

C........  If we need to create the ungridding matrix...
        IF( UFLAG ) THEN

C............  Create a new sorting index
            ALLOCATE( IDXSRT2( NMATX ), STAT=IOS )
            CALL CHECKMEM( IOS, 'IDXSRT2', PROGNAME )

            DO J = 1, NMATX
                IDXSRT2( J ) = J
            END DO

C.........  Sort the scratch gridding matrix arrays to organize by county/link
C           and then cell
            CALL SORTI2( NCOEF, IDXSRT2, NCL, IC )

C..........  Allocate memory for county/link VMT within grid
            ALLOCATE( CLIDX( NSRC ), STAT=IOS )
            CALL CHECKMEM( IOS, 'CLIDX', PROGNAME )
            ALLOCATE( CNT_CL( NCOULNK ), STAT=IOS )
            CALL CHECKMEM( IOS, 'CNT_CL', PROGNAME )
            ALLOCATE( VMT_CELL( NCOULNK,MXCCL ), STAT=IOS )
            CALL CHECKMEM( IOS, 'VMT_CELL', PROGNAME )
            ALLOCATE( VMT_FRAC( NCOULNK,MXCCL ), STAT=IOS )
            CALL CHECKMEM( IOS, 'VMT_FRAC', PROGNAME )
            ALLOCATE( VMT_CL_INV( NCOULNK ), STAT=IOS )
            CALL CHECKMEM( IOS, 'VMT_CL_INV', PROGNAME )
            ALLOCATE( vmt_label( NCOULNK ), STAT=IOS )
            CALL CHECKMEM( IOS, 'vmt_label', PROGNAME )
            CLIDX       = 0    ! array
            CNT_CL      = 0    ! array
            VMT_CELL    = 0    ! array
            VMT_FRAC    = 0.   ! array
            VMT_CL_INV  = 0.   ! array
            vmt_label = 0

C..........  Compute county/link VMT total within grid
            DO K = 1, NCOEF
                J = IDXSRT2( K )
                S = IS ( J )
                N = NCL( J )
                IF ( N .EQ. 0 ) CYCLE   ! skip records outside the grid

C...............  County and/or link total VMT
                if ( s .gt. 0 ) then
                    vmt_label( n ) = ifip( s )
                else
                    vmt_label( n ) = 0
                endif
                VMT_CL_INV( N ) = VMT_CL_INV( N ) + VMT(S) * CS(J)

            END DO       ! end loop to compute VMT totals

C..........  Inverse of county and/or link total VMT
            write(87,*) ' N   IFIP  County VMT in grid [miles/yr]'
            DO N = 1, NCOULNK
                write (87,*) n, vmt_label(n), vmt_cl_inv( n ) 
                IF( VMT_CL_INV( N ) .GT. 0. ) VMT_CL_INV( N ) = 
     &                                             1. / VMT_CL_INV( N )
            END DO

C..........  Loop through cell-source intersections and compute county total
C            VMT by cell/source over County total VMT
            LN = -1
            LC = -1
            CNTMAX = 0
            DO K = 1, NCOEF
   
                J   = IDXSRT2( K ) 
                C   = IC( J )
                S   = IS( J )
                N   = NCL( J )
                IF ( N .EQ. 0 ) CYCLE   ! skip records outside the grid

                CLIDX( S ) = N

                IF ( N .NE. LN ) CNT = 0
                IF ( C .NE. LC ) CNT = CNT + 1
                IF( CNT .LE. MXCCL ) THEN
                    CNT_CL  ( N )     = CNT
                    VMT_CELL( N,CNT ) = C
                    VMT_FRAC( N,CNT ) =  VMT_FRAC( N,CNT ) + 
     &                                   VMT(S) * CS(J) * VMT_CL_INV(N)
                END IF
                IF( CNT > CNTMAX ) THEN
                    CNTMAX = CNT
                END IF 
                LN = N
                LC = C
            END DO

        END IF

        IF( CNTMAX .GT. MXCCL ) THEN
            WRITE( MESG,94010 )
     &       'INTERNAL ERROR: Ungridding matrix not written'
     &       // CRLF() // BLANK10 //
     &       'Arrays would have overflowed.' 
     &       // CRLF() // BLANK10 // 
     &       'Current maximum cells per county/link (MXCCL) =', MXCCL,
     &       '.' // CRLF() // BLANK10 // 
     &       'Actual  maximum cells per county/link    =', CNTMAX  , '.'
            CALL M3MSG2( MESG )
            CALL M3EXIT( PROGNAME, 0, 0, ' ', 2 )
        END IF

C.........  Create ungridding matrix. The factors have been created so that any 
C           non-link sources within a county have the same factors used for
C           computing temperatures.  This means that a source (e.g., urban 
C           interstates) that appear in only some cells in the county will
C           use temperatures that are based on all sources in the county).
C           The link sources within a county have different factors than the
C           non-link sources, and these use only the cells that the link intersects.

        write(88,*) '  source     N Cnt  Col  Row UngridFactor'
        K  = 0
        ln = -1
        DO S = 1, NSRC
            N = CLIDX( S )

C............  Skip sources that are outside the grid
            IF( N .EQ. 0 ) THEN
                NU( S ) = 0
                CYCLE

C............  Store the number of cells per source as the number of cell 
C              intersections with the county or link
            ELSE
                NU( S ) = CNT_CL( N )
            END IF

            sum = 0
C............  Store the cell numbers and VMT fractions into the ungridding matrix
            DO I = 1, CNT_CL( N )
                K = K + 1
                IU( K ) = VMT_CELL( N,I )
                CU( K ) = VMT_FRAC( N,I )
                row = iu(k) / ncols + 1
                col = iu(k) - (row-1) * ncols
                if( cu(k) .ne. 0 ) 
     &              write (88,'(I8,I6,I4,I5,I5,E12.5)' ) 
     &                     s,n,i,col,row,cu(k)
                sum = sum + cu(k)
            END DO
            write( 88,'(40X, E12.5)' ) sum
            if ( sum .NE. 0. .AND.
     &           ( sum > 1.0001 .or. sum < 0.9999 ) ) then
                print *,'warning: sum=', sum
            endif
            ln = 0
            
        END DO
        NCOEFU = K

C.........  Check for overflow
        IF( NCOEFU .GT. NMATXU ) THEN

            WRITE( MESG,94010 )
     &       'INTERNAL ERROR: Ungridding matrix not written'
     &       // CRLF() // BLANK10 //
     &       'Arrays would have overflowed.' 
     &       // CRLF() // BLANK10 // 
     &       'Current cell-source intersections (NMATXU) =', NMATXU, '.'
     &       // CRLF() // BLANK10 // 
     &       'Actual  cell-source intersections          =', NCOEFU, '.'
            CALL M3MSG2( MESG )
            CALL M3EXIT( PROGNAME, 0, 0, ' ', 2 )

        END IF

C.........  Write out ungridding matrix
 
        MESG = 'Writing out UNGRIDDING MATRIX file...'
        CALL M3MSG2( MESG )

        IF ( .NOT. WRITE3( UNAME, 'ALL', 0, 0, NU ) ) THEN
            MESG = 'Problem writing UNGRIDDING MATRIX file.'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

C.........  Report FIPS that don't have surrogate data
C.........  Report links that are outside the grid
c        CALL RPSRCOUT( NNOSRG, NLKOGRD, FIPNOSRG, LKOGRD )

C.........  Dellallocate locally allocated memory
        DEALLOCATE( IS, CS, NCL, IDXSRT2, INDOMAIN )
        DEALLOCATE( ACEL, AFAC, FIPNOSRG, LKOGRD )

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Formatted file I/O formats............ 93xxx

93360   FORMAT( I8, 1X, I4, 1X, I4 )

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I9, :, 1X ) )
 
        END SUBROUTINE GENMGMAT

