
        SUBROUTINE GENMGMAT( GNAME, UNAME, FDEV, MXSCEL, MXCSRC, NMSRC, 
     &                       NGRID, NMATX, NX, IX, CX, NU, IU, CU,
     &                       NCOEF, CMAX, CMIN, NCOEFU, CMAXU, CMINU )

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
        LOGICAL         DSCM3GRD

        EXTERNAL        CRLF, ENVINT, FIND1, FINDC, DSCM3GRD

C...........   SUBROUTINE ARGUMENTS
        CHARACTER(*), INTENT (IN) :: GNAME         ! gridding mtx logical name
        CHARACTER(*), INTENT (IN) :: UNAME         ! ungridding mtx logical name
        INTEGER     , INTENT (IN) :: FDEV          ! surg codes report file
        INTEGER     , INTENT (IN) :: MXSCEL        ! max sources per cell
        INTEGER     , INTENT (IN) :: MXCSRC        ! max cells per source
        INTEGER     , INTENT (IN) :: NMSRC         ! no. mobile sources
        INTEGER     , INTENT (IN) :: NGRID         ! actual grid cell count
        INTEGER     , INTENT (IN) :: NMATX         ! no. source-cell intersects
        INTEGER     , INTENT(OUT) :: NX  ( NGRID ) ! no. srcs per cell
        INTEGER     , INTENT(OUT) :: IX  ( NMATX ) ! src IDs 
        REAL        , INTENT(OUT) :: CX  ( NMATX ) ! gridding coefficients
        INTEGER     , INTENT(OUT) :: NU  ( NMSRC ) ! no. cells per source
        INTEGER     , INTENT(OUT) :: IU  ( NMATX ) ! cell numbers
        REAL        , INTENT(OUT) :: CU  ( NMATX ) ! ungridding coefficients
        INTEGER     , INTENT(OUT) :: NCOEF         ! no. of gridding coeffs
        INTEGER     , INTENT(OUT) :: CMAX          ! max no. of sources per cell
        INTEGER     , INTENT(OUT) :: CMIN          ! min no. of sources per cell
        INTEGER     , INTENT(OUT) :: NCOEFU        ! no. of ungridding coeffs
        INTEGER     , INTENT(OUT) :: CMAXU         ! max no. of cells per source
        INTEGER     , INTENT(OUT) :: CMINU         ! min no. of cells per source

C...........   Local allocatable arrays...

C...........   Scratch Gridding Matrix (subscripted by source-within-cell, cell)

        INTEGER, ALLOCATABLE :: IDXSRT( : )! sorting index
        INTEGER, ALLOCATABLE :: IC    ( : )! cell IDs
        INTEGER, ALLOCATABLE :: IS    ( : )! source IDs
        REAL   , ALLOCATABLE :: CS    ( : )! facs (no adjustments), for ungrid
        REAL   , ALLOCATABLE :: CSJ   ( : )! factors (with adjustments)

C...........   Scratch Ungridding Matrix (subscripted by cell, source)
  
        INTEGER, ALLOCATABLE :: IT( :,: ) ! cell IDs for each source
        REAL   , ALLOCATABLE :: CT( :,: ) ! ungridding coefficients
        REAL   , ALLOCATABLE :: DN( : )   ! source normalization coefficient

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
C...........   Other local variables

        INTEGER         C, F, I, J, K, N, S !  indices and counters.

        INTEGER         FIP      !  tmp country/state/county code
        INTEGER         IOS      !  i/o status
        INTEGER         ISIDX    !  surrogate ID code index
        INTEGER         L2       !  string length
        INTEGER         LFIP     !  cy/st/co code from previous iteration
        INTEGER         LNKEND   !  width of sourc info to end of link ID
        INTEGER         NCEL     !  tmp number of cells 
        INTEGER         NLKOGRD  !  no. of link sources outside the grid
        INTEGER         NMAX     !  counter for storing correct max dimensions
        INTEGER         NNOSRG   !  no. of cy/st/co codes with no surrogates
        INTEGER         RWT      !  tmp roadway type

        REAL            ADJ, ADJC   ! tmp adjustment factors
        REAL            ALEN        ! link length
        REAL            FRAC        ! tmp surrogate fraction
        REAL            XBEG, YBEG  ! tmp X and Y link start coordinates
        REAL            XEND, YEND  ! tmp X and Y link end   coordinates

        LOGICAL      :: EFLAG = .FALSE.  !  true: error detected

        CHARACTER*16    COORD     !  coordinate system name
        CHARACTER*16    COORUNIT  !  coordinate system projection units
        CHARACTER*16    GRDNM     !  grid name
        CHARACTER*80    GDESC     !  grid description
        CHARACTER*300   BUFFER    !  source fields buffer
        CHARACTER*300   MESG      !  message buffer 

        CHARACTER(LEN=LNKLEN3)    CLNK   ! tmp link ID
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

C.........  Initialize number of sources per cell counter
        NX = 0   ! Array

C.........  Allocate memory for temporary gridding matrix and other
        ALLOCATE( IDXSRT( NMATX ), STAT=IOS )
        CALL CHECKMEM( IOS, 'IDXSRT', PROGNAME )
        ALLOCATE( IC( NMATX ), STAT=IOS )
        CALL CHECKMEM( IOS, 'IC', PROGNAME )
        ALLOCATE( IS( NMATX ), STAT=IOS )
        CALL CHECKMEM( IOS, 'IS', PROGNAME )
        ALLOCATE( CS( NMATX ), STAT=IOS )
        CALL CHECKMEM( IOS, 'CS', PROGNAME )
        ALLOCATE( CSJ( NMATX ), STAT=IOS )
        CALL CHECKMEM( IOS, 'CSJ', PROGNAME )
        ALLOCATE( IT( MXCSRC, NMSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'IT', PROGNAME )
        ALLOCATE( CT( MXCSRC, NMSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'CT', PROGNAME )
        ALLOCATE( DN( NMSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'DN', PROGNAME )
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
        LNKEND  = VIDPOS3 - 1
        NNOSRG  = 0
        NLKOGRD = 0
        ADJ     = 1.    !  temporary
        ADJC    = 1.    !  temporary

        DO S = 1, NSRC
            
            FIP  = IFIP  ( S )
            RWT  = IRCLAS( S )
            CLNK = CLINK ( S )
            CSRC = CSOURC( S )
            CSRC2= CSRC( 1:VIDPOS3-1 )  ! abridged for these purposes

C.............  Initialize sources as being in the domain
            INDOMAIN( S ) = .TRUE.

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
                CALL SETFRAC( FDEV, S, ISIDX, K, F, 2, INDOMAIN( S ), 
     &                        CSRC2, FRAC )

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
                        CS    ( J ) = FRAC
                        CSJ   ( J ) = ADJ * ADJC * FRAC
                    END IF

                END IF  ! if surrogate fraction > 0.

            END DO    !  end of loop on cells K for this FIP

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
        CMAX = MAXVAL( NX )
        CMIN = MINVAL( NX )
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

C............................................................................
C..........   Generate ungridding matrix ....................................
C............................................................................

C.........  Initialize ungridding matrix related arrays
        NU = 0   ! array
        DN = 0.  ! array

C.........  Create renormalization factors for each source
        DO K = 1, NCOEF  ! Loop through cell-src intersections

            J = IDXSRT( K ) 
            S = IS( J )                          
            DN( S ) = DN( S ) + CS( J )

        END DO

C.........  Pre-divide normalization factors - make sure not zero
C.........  Have already given notes and warnings about zero-surrogates
        DO S = 1, NMSRC

            IF( DN( S ) .GT. 0. ) DN( S ) = 1. / DN( S )

        END DO

C.........  Renormalize and Transpose ((IS,CS) into scratch arrays (IT,CT):
        NMAX = -1
        DO K = 1, NCOEF  ! Loop through cell-src intersections

            J = IDXSRT( K ) 
            C = IC( J )
            S = IS( J )

            N = NU( S ) + 1

C.............  Ensure that the number of cells per source is okay
            IF ( N .LE. MXCSRC ) THEN
                IT( N,S ) = C
                CT( N,S ) = DN( S ) * CS( J )

C.............  Keep track of the maximum sources per cell for reporting
            ELSE
                IF( N .GT. NMAX ) NMAX = N
            END IF

            NU( S ) = N

        END DO

        IF( NMAX .GT. MXCSRC ) THEN
            WRITE( MESG,94010 )
     &       'INTERNAL ERROR: Ungridding matrix not written'
     &       // CRLF() // BLANK10 //
     &       'Arrays would have overflowed.' 
     &       // CRLF() // BLANK10 // 
     &       'Current maximum cells per source (MXCSRC) =', MXCSRC, '.'
     &       // CRLF() // BLANK10 // 
     &       'Actual  maximum cells per source          =', NMAX  , '.'
            CALL M3MSG2( MESG )
            CALL M3EXIT( PROGNAME, 0, 0, ' ', 2 )
        END IF

C.........  Compress matrix into I/O representation from scratch 
C           representation and compute statistics.
        K  = 0
        DO I = 1, NU( 1 )
            K = K + 1
            IU( K ) = IT( I,1 )
            CU( K ) = CT( I,1 )
        END DO
                   
        CMAXU = NU( 1 )
        CMINU = CMAXU
        DO S = 2, NMSRC

            J = NU( S )
            IF ( J .GT. CMAXU ) THEN
                CMAXU = J
            ELSE IF ( J .LT. CMINU ) THEN
                CMINU = J
            END IF
                
            DO I = 1, J
                K = K + 1
                IF( K .LE. NMATX ) THEN
                    IU( K ) = IT( I,S )
                    CU( K ) = CT( I,S )
                END IF
            END DO

        END DO

        NCOEFU = K

C.........  Check for overflow
        IF( NCOEFU .GT. NMATX ) THEN

            WRITE( MESG,94010 )
     &       'INTERNAL ERROR: Ungridding matrix not written'
     &       // CRLF() // BLANK10 //
     &       'Arrays would have overflowed.' 
     &       // CRLF() // BLANK10 // 
     &       'Current cell-source intersections (NMATX) =', NMATX, '.'
     &       // CRLF() // BLANK10 // 
     &       'Actual  cell-source intersections         =', NCOEFU, '.'
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
        DEALLOCATE( IS, CS, CSJ, IT, CT, DN, INDOMAIN )
        DEALLOCATE( ACEL, AFAC, FIPNOSRG, LKOGRD )

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I9, :, 1X ) )
 
        END SUBROUTINE GENMGMAT

