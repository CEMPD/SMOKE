
C Version "@(#)$Id$ $Source$ $Date$ 

        PROGRAM GRDMMAT

C***********************************************************************
C  program body starts at line 207
C
C  DESCRIPTION:
C       Construct mobile source gridding and ungridding matrices from 
C       data contained in EPS-style surrogates file and/or links file.
C
C  PRECONDITIONS REQUIRED:
C       Sorted, cut-down input data for surrogate coeffs.
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C       Models-3 I/O
C       FIND1, FIND2, FIND3, GETYN, PROMPTFFILE, PROMPTMFILE, TRIMLEN
C
C  REVISION  HISTORY:
C       Prototype 1/96 by CJC 
C***********************************************************************

      IMPLICIT NONE

C...........   INCLUDES:
        
        INCLUDE 'GRDIMS3.EXT'   !  grid parameters
        INCLUDE 'MBDIMS3.EXT'   !  mobile-source dimensioning parameters
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
        INCLUDE 'FDESC3.EXT'    !  I/O API file description data structures.


C...........   EXTERNAL FUNCTIONS and their descriptions:
        
        LOGICAL         DSCGRID
        INTEGER         FIND1, FIND2, FIND3
        LOGICAL         GETYN
        INTEGER         PROMPTFFILE
        CHARACTER*16    PROMPTMFILE
        LOGICAL         SETLAM
        INTEGER         TRIMLEN
        LOGICAL         UTM2LAM
        
        EXTERNAL        DSCGRID, FIND1, FIND2, FIND3, GETYN, 
     &                  PROMPTFFILE, PROMPTMFILE, TRIMLEN, 
     &                  SETLAM, UTM2LAM


C...........   PARAMETER:

        INTEGER         NMMAX   !  max number of gridding-xref entries
        INTEGER         MXLNKD  !  max number of entries in link defs file

        PARAMETER     ( NMMAX  = NRCLAS*( 1 + MXSTA + MXFIP + MXLNK ),
     &                  MXLNKD = NRCLAS * MXFIP )

C...........   LOCAL VARIABLES and their descriptions:

        INTEGER     FIP
        INTEGER     SID
        INTEGER     RDT
        INTEGER     LNK
        INTEGER     LFIP, LRDT
        INTEGER     COL
        INTEGER     ROW
        INTEGER     ZONE, GZONE, UZONE
        INTEGER     CMAX, CMIN
        REAL        XDUM, YDUM
        REAL        XBEG, YBEG, XEND, YEND
        REAL        X00, Y00, XCELL, YCELL, XX, YY, DDX, DDY, DD
        REAL        FRAC( NMSRG )         !  (scratch) surrogate fractions
        REAL        DIV

C...........   LNKDEF tables:
        
        INTEGER     NUMLNK           !  number of entries in the link defs table
        INTEGER     LNKFIP( MXLNKD ) !  fip code
        INTEGER     LNKRDT( MXLNKD ) !  road type
        INTEGER     LNKCNT( MXLNKD ) !  number of cells for this fip, road type
        INTEGER     LNKCEL( MXCFIP, MXLNKD )
        REAL        LNKFAC( MXCFIP, MXLNKD )
        REAL        LNKTOT( MXLNKD )
        
        INTEGER     NCEL
        INTEGER     ACEL( MXCFIP )
        REAL        AFAC( MXCFIP )
        REAL        ALEN

C...........   XREF before sorting:                                 

        INTEGER     NUMSCA           !  current surrogates-table count
        INTEGER     INDEXA( NMMAX )  !  FIP associated with this surrogate
        INTEGER     SRFIPA( NMMAX )  !  FIP associated with this surrogate
        INTEGER     SRGIDA( NMMAX )  !  surrogate number (1...15)
        INTEGER     SRRDTA( NMMAX )  !  road class
        INTEGER     SRLNKA( NMMAX )  !  link ID
                
C...........   XREF after sorting and processing (four tables, one for: 
C.......   fallback, one for roadtype-only, one for state+roadtype
C.......   one for FIP+roadtype

        INTEGER     SRGID0           !  fallback surrogate number

        INTEGER     NUMSCR           !  RDT:  current surrogates-table count
        INTEGER     SRRDTR( NRCLAS ) !  road type
        INTEGER     SRGIDR( NRCLAS ) !  surrogate number (1...15)

        INTEGER     NUMSCS           !  State-RDT table:
        INTEGER     SRFIPS( MXSREF ) !  state associated with this surrogate
        INTEGER     SRRDTS( MXSREF ) !  road type
        INTEGER     SRGIDS( MXSREF ) !  surrogate number (1...15)

        INTEGER     NUMSCF           !  FIP-RDT table:
        INTEGER     SRFIPF( MXFREF ) !  FIP associated with this surrogate
        INTEGER     SRRDTF( MXFREF ) !  road type
        INTEGER     SRGIDF( MXFREF ) !  surrogate number (1...15)

C...........   Surrogate-cell::FIPS table
        
        INTEGER     NFIPLIST                         !  # of FIPS
        INTEGER     NCELLS (                MXFIP )  !  # cells per FIP
        INTEGER     FIPLIST(                MXFIP )  !  FIP-list
        INTEGER     FIPCELL(        MXCFIP, MXFIP )  !  cty codes::cell number
        REAL        SRGFRAC( NMSRG, MXCFIP, MXFIP )  !  surrogate fractions
                
C...........   Mobile Sources Table
        
        INTEGER     IFIP  ( NMSRC )  !  source FIPS (county) ID
        INTEGER     IRCLAS( NMSRC )  !  EPA roadclass code
        INTEGER     ILINK ( NMSRC )  !  Link number
        REAL        XLOC1 ( NMSRC )  !   X coordinate for start of link
        REAL        XLOC2 ( NMSRC )  !   X coordinate for   end of link
        REAL        YLOC1 ( NMSRC )  !   Y coordinate for start of link
        REAL        YLOC2 ( NMSRC )  !   Y coordinate for   end of link
        
C...........   Scratch Gridding Matrix (subscripted by source-within-cell, cell)

        INTEGER     IS( MXSCEL , NGRID )
        REAL        CS( MXSCEL , NGRID )

C...........   Gridding Matrix

        INTEGER     NX( NGRID )
        INTEGER     IX( NMATX )
        REAL        CX( NMATX )

        COMMON  / GRIDMAT / NX, IX, CX

C...........   Scratch Ungridding Matrix (subscripted by cell, source)

        INTEGER     IT( MXCFIP , NMSRC )
        REAL        CT( MXCFIP , NMSRC )
        REAL        DN( NMSRC )           ! Source Normalization coefficient

C...........   Ungridding Matrix

        INTEGER     NU( NMSRC )
        INTEGER     IU( NMATX )
        REAL        CU( NMATX )

        COMMON  / UGRIDMAT / NU, IU, CU

C...........   Other local variables
        
        INTEGER         C, R, S, F, I, J, K  !  loop counters.

        CHARACTER*16    CNAME   !  coord sys name
        CHARACTER*16    ENAME   !  logical name for emission source  input file
        CHARACTER*16    GNAME   !  logical name for grid   matrix   output file
        CHARACTER*16    UNAME   !  logical name for UNgrid matrix   output file
        INTEGER         SDEV    !  for surrogate coeff file
        INTEGER         XDEV    !  for surrogate xref  file
        INTEGER         LDEV    !  for link defs file
        INTEGER         ILOG    !  for log-device

        INTEGER         NOSURG( NMSRC ) ! FIPs codes w/o surrogates
        
        INTEGER         IOS     !  I/O status
        INTEGER         IREC    !  input line (record) number
        INTEGER         I0      !  for counting FIP-independent srcs
        INTEGER         I1      !  for counting state-dependent   srcs
        INTEGER         I2      !  for counting FIP-dependent   srcs
        INTEGER         NOS     !  number of FIPs w/o surrogates

        LOGICAL         EFLAG   !  input error flat
        LOGICAL         EPSFLG  !  if true, the surrogates are in EPS format
        LOGICAL         GFLAG   !  generate output   gridding matrix if true
        LOGICAL         UFLAG   !  generate output UNgridding matrix if true
        LOGICAL         OFLAG   !  Did overflow occur
        CHARACTER*160   MESG



C...........   STATEMENT FUNCTION

        LOGICAL         FLTERR
        REAL            PP, QQ
        FLTERR( PP,QQ ) = 
     &      ( (PP - QQ)**2 .GT. 1.0E-10*( PP*PP + QQ*QQ + 1.0E-5 ) )


C***********************************************************************
C   begin body of program GRDMMAT
        
        ILOG = INIT3()
        
        WRITE( *,92000 ) 
     &  ' ',
     &  'Program GRDMMAT to take the sorted surrogate coefficient',
     &  'file, the EPS surrogate cross reference file, the MOBILE',
     &  'SOURCE file produced by RAWMOBIL, and produce the MOBILE',
     &  'SOURCE GRIDDING MATRIX file.', 
     &  ' '
        WRITE( *,92000 ) 
     &  'You will need to enter the logical names for the input and',
     &  'output files (and to have set them prior to program launch,',
     &  'using "setenv <logicalname> <pathname>").  Input files must',
     &  'have been sorted as indicated, prior to program execution.',
     &  ' '
        WRITE( *,92000 ) 
     &  'You may use END_OF-FILE (control-D) to quit the program',
     &  'during logical-name entry. Default responses are indicated',
     &  'in brackets [LIKE THIS].',
     &  ' '
        
        IF ( .NOT. GETYN( 'Continue with program?', .TRUE. ) ) THEN
            CALL M3EXIT( 'GRDMMAT', 0, 0, 'Ending program GRDMMAT', 2 )
        END IF


C.......   Get file name; open actual-scc file
        
        XDEV = PROMPTFFILE( 
     &  'Enter logical name for GRIDDING SURROGATE XREF file',
     &           .TRUE., .TRUE., 'MGREF', 'GRDMMAT' )

        EPSFLG = GETYN( 'EPS formatted surrogates?', .TRUE. )

        SDEV = PROMPTFFILE( 
     &  'Enter logical name for SURROGATE COEFF file',
     &           .TRUE., .TRUE., 'MGPRO', 'GRDMMAT' )

        LDEV = PROMPTFFILE( 
     &  'Enter logical name for LINK DEFS file (or "NONE")',
     &           .TRUE., .TRUE., 'LNKDEF', 'GRDMMAT' )

        ENAME = PROMPTMFILE( 
     &          'Enter logical name for MOBILE SOURCE input file',
     &          FSREAD3, 'MOBL', 'GRDMMAT' )

        IF( .NOT. DESC3( ENAME ) ) THEN
            CALL M3EXIT( 'GRDMMAT', 0, 0,
     &              'Could not get description of file ' // ENAME, 2 )
        END IF

        UZONE = NINT( P_ALP3D )

        IF( NROWS3D .NE. NMSRC ) THEN
            WRITE( *, 92000 ) ' '
            WRITE( *, 92010 )
     &             'Existing MOBL file has record count      ', NROWS3D,
     &             'But, GRDMMAT compiled using count (NMSRC)', NMSRC
            CALL M3EXIT( 'GRDMMAT', 0, 0,
     &        'Max record-count exceeded in MOBILE SOURCE file.', 2 )

        ENDIF

C...........   Get grid description from GRIDDESC file (Sets values
C...........   in FDESC3 commons for opening MGMAT, MUMAT files):
        
        IF ( .NOT. DSCGRID( GRDNM, CNAME, GDTYP3D, 
     &             P_ALP3D, P_BET3D, P_GAM3D, XCENT3D, YCENT3D,
     &             XORIG3D, YORIG3D, XCELL3D, YCELL3D, 
     &             NCOLS3D, NROWS3D, NTHIK3D ) ) THEN

            GNAME = GRDNM
            MESG  = 'Could not get description of grid "' //
     &              GNAME( 1:TRIMLEN( GNAME ) ) // 
     &              '" from GRIDDESC file.'
            CALL M3EXIT( 'GRDMMAT', 0, 0, MESG, 2 )

        END IF          !  if dscgrid() failed

        IF ( GDTYP3D .EQ. LAMGRD3 ) THEN

            IF ( .NOT. SETLAM( SNGL( P_ALP3D) , SNGL( P_BET3D ), 
     &                         SNGL( P_GAM3D), 
     &                         SNGL( XCENT3D), SNGL( YCENT3D ) ) ) THEN
                MESG  = 'Failure in SETLAM for grid "' //
     &                  GNAME( 1:TRIMLEN( GNAME ) ) // '"'
                CALL M3EXIT( 'GRDMMAT', 0, 0, MESG, 2 )
            END IF

        ELSE IF ( GDTYP3D .EQ. UTMGRD3 ) THEN

            GZONE = NINT( P_ALP3D )

        ELSE IF ( GDTYP3D .NE. LATGRD3 ) THEN

            WRITE( MESG,94010 ) 
     &      'Grid type', GDTYP3D, 'not yet supported'
            CALL M3EXIT( 'GRDMMAT', 0, 0, MESG, 2 )

        END IF


C.......   Get file name; open output gridding matrix file

        FTYPE3D = SMATRX3
        SDATE3D = 0
        STIME3D = 0
        TSTEP3D = 0
        NVARS3D = 1
        NCOLS3D = NMATX
        NROWS3D = NGRID
        NLAYS3D = 1
        NTHIK3D = NMSRC
        GDNAM3D = GRDNM
        VGTYP3D = IMISS3
        VGTOP3D = BADVAL3
        FDESC3D( 1 ) = 'Mobile Source gridding-coefficient matrix'
        DO  1  K = 2, MXDESC3
            FDESC3D( K ) = ' '
1       CONTINUE
        VNAME3D( 1 ) = 'MGRDMAT'
        UNITS3D( 1 ) = 'n/a'
        VDESC3D( 1 ) = 'Mobile Source gridding-coefficients'
        VTYPE3D( 1 ) = M3REAL

        GNAME = PROMPTMFILE( 
     &          'Enter logical name for GRIDDING MATRIX output file',
     &          FSUNKN3, 'MGMAT', 'GRDMMAT' )

        GFLAG = ( GNAME( 1:5 ) .NE. 'NONE ' )


C.......   Get file name; open output UNgridding matrix file

        FTYPE3D = SMATRX3
        SDATE3D = 0
        STIME3D = 0
        TSTEP3D = 0
        NVARS3D = 1
        NCOLS3D = NMATX
        NROWS3D = NMSRC
        NLAYS3D = 1
        NTHIK3D = NGRID
        GDNAM3D = GRDNM
        VGTYP3D = IMISS3
        VGTOP3D = BADVAL3
        FDESC3D( 1 ) = 'Mobile Source UNgridding-coefficient matrix'
        DO  2  K = 2, MXDESC3
            FDESC3D( K ) = ' '
2       CONTINUE
        VNAME3D( 1 ) = 'MUGRDMAT'
        UNITS3D( 1 ) = 'n/a'
        VDESC3D( 1 ) = 'Mobile Source UNgridding-coefficients'
        VTYPE3D( 1 ) = M3REAL

        UNAME = PROMPTMFILE( 
     &          'Enter logical name for UNGRIDDING MATRIX output file',
     &          FSUNKN3, 'MUMAT', 'GRDMMAT' )

        UFLAG = ( UNAME( 1:5 ) .NE. 'NONE ' )


C.......   If it exists, read the link definitions file:
C...............   ASSUMPTION:  sorted by FIP and road class:

        NUMLNK = 0

        IF ( LDEV .LT. 0 )  GO TO  44

        EFLAG = .FALSE.
        IREC  =  0
        S     =  0
        LFIP  = -1
        LRDT  = -1
        
11      CONTINUE        !  head of loop reading LINK DEFS file
            
            READ( LDEV, *, END = 44, IOSTAT = IOS )
     &          FIP,
     &          RDT,
     &          XBEG, YBEG, XEND, YEND,
     &          ZONE

            IREC = IREC + 1
        
            IF ( IOS .GT. 0 ) THEN

                WRITE( MESG,94010 ) 
     &              'Error', IOS, 
     &              'reading LINK DEFS file at line', IREC
                EFLAG = .TRUE.
                CALL M3MESG( MESG )
                GO TO  11               !  to head of loop:  read next line

            ELSE IF ( FIP .NE. LFIP  .OR.  RDT .NE. LRDT ) THEN

                S    = S + 1
                LFIP = FIP
                LRDT = RDT
                IF ( S .LE. MXLNKD ) THEN
                    LNKFIP( S ) = FIP
                    LNKRDT( S ) = RDT
                    LNKCNT( S ) = 0
                    LNKTOT( S ) = 0.0
                ELSE
                    EFLAG = .TRUE.
                    WRITE( MESG,94010 )
     &              'Link table overflow at line', IREC
                    CALL M3MESG( MESG )
                    GO TO  11               !  to head of loop:  read next line
                END IF

            END IF

C...........   Perform coordinate transformations if necessary:
            
            IF ( GDTYP3D .EQ. LATGRD3 ) THEN
                CALL UTM2LL( XBEG, YBEG, ZONE,  XDUM, YDUM )
                XBEG = XDUM
                YBEG = YDUM

                CALL UTM2LL( XEND, YEND, ZONE,  XDUM, YDUM )
                XEND = XDUM
                YEND = YDUM

            ELSEIF ( GDTYP3D .EQ. LAMGRD3 ) THEN

                IF ( .NOT.  UTM2LAM( XBEG, YBEG, ZONE, 
     &                               XDUM, YDUM ) ) THEN
                    EFLAG = .TRUE.
                    WRITE( MESG,94010 ) 
     &              'Error converting <XBEG,YBEG> to Lambert at line', 
     &              IREC
                    CALL M3MESG( MESG )
                    GO TO  11               !  to head of loop:  read next line
                END IF
                XBEG = XDUM
                YBEG = YDUM

                IF ( .NOT.  UTM2LAM( XEND, YEND, ZONE, 
     &                               XDUM, YDUM ) ) THEN
                    EFLAG = .TRUE.
                    WRITE( MESG,94010 ) 
     &              'Error converting <XEND,YEND> to Lambert at line', 
     &              IREC
                    CALL M3MESG( MESG )
                    GO TO  11               !  to head of loop:  read next line
                END IF
                XEND = XDUM
                YEND = YDUM

            ELSE IF ( GDTYP3D .EQ. UTMGRD3 ) THEN 

                IF ( ZONE .NE. GZONE ) THEN
                    CALL UTM2LL( XBEG, YBEG, ZONE,  XDUM, YDUM )
                    CALL LL2UTM( XDUM, YDUM, GZONE, XBEG, YBEG )

                    CALL UTM2LL( XEND, YEND, ZONE,  XDUM, YDUM )
                    CALL LL2UTM( XDUM, YDUM, GZONE, XEND, YEND )
                END IF

            END IF

C...........   Compute conversion fractions:

            CALL LNK2GRD( XBEG, YBEG, XEND, YEND, 
     &                    NCEL, ACEL, AFAC, ALEN )

            IF ( NCEL .GT. MXCFIP ) THEN
                EFLAG = .TRUE.
                WRITE( MESG,94010 )
     &          'LNK:CEL table overflow for source', S, 'at line', IREC
                CALL M3MESG( MESG )
                GO TO  11               !  to head of loop:  read next line
            END IF
            
            K = LNKCNT( S )

            DO  33  I = 1, NCEL

                LNKTOT( S ) = LNKTOT( S ) + ALEN
                
                DO  22  J = 1, K
                    IF ( ACEL( I ) .EQ. LNKCEL( J,S ) ) THEN
                        LNKFAC( J,S ) = LNKFAC( J,S ) + AFAC( I )
                        GO TO  33
                    END IF
22              CONTINUE

                K = K + 1       !  if you get to here:  new cell for this source
                IF ( K .LE. MXCFIP ) THEN
                    LNKCNT(   S ) = K
                    LNKFAC( K,S ) = AFAC( I )
                ELSE
                    EFLAG = .TRUE.
                    WRITE( MESG,94010 ) 
     &              'Link table overflow for source', S, 'at line', IREC
                    CALL M3MESG( MESG )
                    GO TO  11
                END IF

33          CONTINUE
            
            GO TO  11           !  end of loop:  "read LINK DEFS file"
            

44      CONTINUE        !  exit from "read LINK DEFS file"

        IF ( EFLAG ) THEN
             CALL M3EXIT( 'GRDMMAT', 0,0, 
     &                    'Error reading LINK DEFS file', 2 )
        END IF
        
        NUMLNK = S
        
C.......   Rescale (normalize) the LNKFAC's:
        
        DO  55  S = 1, NUMLNK
            IF ( LNKCNT( S ) .GT. 0 ) THEN
                DIV = 1.0 / LNKTOT( S )       !  OK by construction of LNK2GRD
                DO  54  C = 1, LNKCNT( S )
                    LNKFAC( C,S ) = DIV * LNKFAC( C,S )
54              CONTINUE
            END IF
55      CONTINUE        !  end loop on S normalizing LNKFACs

C.......   Get grid description for use in check
        IF ( .NOT. DSCGRID( GRDNM, CNAME, GDTYP3D, 
     &             P_ALP3D, P_BET3D, P_GAM3D, XCENT3D, YCENT3D,
     &             XORIG3D, YORIG3D, XCELL3D, YCELL3D, 
     &             NCOLS3D, NROWS3D, NTHIK3D ) ) THEN

            GNAME = GRDNM
            MESG  = 'Could not get description of grid "' //
     &              GNAME( 1:TRIMLEN( GNAME ) ) // 
     &              '" from GRIDDESC file.'
            CALL M3EXIT( 'GRDMMAT', 0, 0, MESG, 2 )

        END IF          !  if dscgrid() failed

C...........   Read the surrogate coefficient file
C.......   Initializations:


        DO  66  C = 1, NGRID
            NX( C ) = 0         !  for Gridding Matrix
66      CONTINUE

        DO  67  S = 1, MXFIP
            NCELLS( S ) = 0     !  for Surrogate  cell::FIPS table
67      CONTINUE


C.......   ASSUMPTION:  sorted by FIP

        WRITE( *,92000 ) 
     &  ' ', 'Reading GRIDDING SURROGATE COEFFICIENT file...', ' '

C.........  Read EPS format file header
        IF( EPSFLG ) THEN  
            READ( SDEV, 93010, IOSTAT=IOS ) !  read file header
     &              X00, Y00,
     &              XDUM, YDUM, ZONE,
     &              XCELL, YCELL
        
            IF ( IOS .NE. 0 ) THEN
                WRITE( MESG,94010 ) 
     &              'I/O error ', IOS, 
     &              'reading SURG COEFF file header at line', IREC
                CALL M3EXIT( 'GRDMMAT', 0, 0, MESG, 1 )
            END IF
                           
            XX = ( XORIG3D - X00 ) / XCELL
            YY = ( YORIG3D - Y00 ) / YCELL
            IF ( IOS .NE. 0 ) THEN
                WRITE( MESG,94010 ) 
     &          'Error', IOS, 'reading header of SURROGATE file'
                CALL M3EXIT( 'GRDMMAT', 0, 0, MESG, 2 )
            ELSE IF ( GDTYP3D .EQ. UTMGRD3  .AND.
     &                ZONE .NE. GZONE ) THEN
                CALL M3EXIT( 'GRDMMAT', 0, 0,
     &                       'Bad UTM ZONE in SURROGATE file', 2 )
            ELSE IF ( FLTERR( XCELL, SNGL( XCELL3D ) ) ) THEN
                CALL M3EXIT( 'GRDMMAT', 0, 0,
     &                       'Bad XCELL in SURROGATE file', 2 )
            ELSE IF ( FLTERR( YCELL, SNGL( YCELL3D ) ) ) THEN
                CALL M3EXIT( 'GRDMMAT', 0, 0,
     &                       'Bad YCELL in SURROGATE file', 2 )
            ELSE IF ( GDTYP3D .EQ. UTMGRD3  .AND.
     &                FLTERR( XX, FLOAT( NINT( XX ) ) ) ) THEN
                CALL M3EXIT( 'GRDMMAT', 0, 0,
     &                   'Bad X-origin alignment in SURROGATE file', 2 )
            ELSE IF ( GDTYP3D .EQ. UTMGRD3  .AND.
     &                FLTERR( YY, FLOAT( NINT( YY ) ) ) ) THEN
            CALL M3EXIT( 'GRDMMAT', 0, 0,
     &                   'Bad Y-origin alignment in SURROGATE file', 2 )
            ELSE IF ( XX .LT. -1.0E-4 * XCELL ) THEN
                CALL M3EXIT( 'GRDMMAT', 0, 0,
     &                   'Bad X-origin placement in SURROGATE file', 2 )
            ELSE IF ( YY .LT. -1.0E-4 * YCELL ) THEN
                CALL M3EXIT( 'GRDMMAT', 0, 0,
     &                   'Bad Y-origin placement in SURROGATE file', 2 )
            ELSE IF ( XDUM .GT. 
     &                XORIG3D + ( NCOLS3D + 1.0E-4 ) * XCELL ) THEN
                CALL M3EXIT( 'GRDMMAT', 0, 0,
     &                   'Bad grid X-placement in SURROGATE file', 2 )
            ELSE IF ( YDUM .GT.
     &                YORIG3D + ( NROWS3D + 1.0E-4 ) * YCELL ) THEN
                CALL M3EXIT( 'GRDMMAT', 0, 0,
     &                   'Bad grid Y-placement in SURROGATE file', 2 )
            END IF

        ENDIF ! EPS Surrogates file

        IREC     =  1
        NFIPLIST =  0
        S        =  0
        LFIP     = -1
        LRDT     = -1
        EFLAG    = .FALSE.
        DDX      = SNGL( 1.0D0 / XCELL3D )
        DDY      = SNGL( 1.0D0 / YCELL3D )
        NUMSCR   =  0      !  for Surrogates::ASC table
        NUMSCS   =  0      !  for Surrogates::ASC table
        NUMSCF   =  0      !  for Surrogates::ASC table
        
70      CONTINUE        !  head of the SDEV-read loop
      
C.............  Read EPS format 
            IF( EPSFLG ) THEN  

                READ( SDEV, 93020 , END=79, IOSTAT=IOS ) 
     &                FIP, XX, YY, ZONE, ( FRAC( C ), C = 1, NMSRG )

                IREC = IREC + 1

                IF ( IOS .NE. 0 ) THEN

                    WRITE( MESG,94010 ) 
     &                  'I/O error ', IOS, 
     &                  'reading SURG COEFF file at line', IREC
                    CALL M3MESG( MESG )
                    EFLAG = .TRUE.
                    GO TO  70

                END IF
            
                COL = 1 + NINT( DDX * ( XX - XORIG3D ) )
                ROW = 1 + NINT( DDY * ( YY - YORIG3D ) )
                IF ( COL .LT. 1  .OR. COL .GT. NCOLS ) THEN

                    WRITE( MESG,94010 )
     &                  'Column ', COL,
     &                  'out of range in SURG COEFF file at line', IREC
                    CALL M3MESG( MESG )
                    GO TO  70

                ELSE IF ( ROW .LT. 1  .OR. ROW .GT. NROWS ) THEN

                    WRITE( MESG,94010 )
     &                  'Row ', ROW,
     &                  'out of range in SURG COEFF file at line', IREC
                    CALL M3MESG( MESG )
                    GO TO  70

                END IF

C.............  Read default SMOKE format
            ELSE

                READ( SDEV, * , END=79, IOSTAT=IOS )
     &                FIP, COL, ROW, ( FRAC( C ), C = 1, NMSRG )
 
                IREC = IREC + 1
 
                IF( IOS .LT. 0 ) THEN
 
                    WRITE( MESG,94010 )
     &              'I/O error, ', IOS,
     &              'reading EMS-95 SURROGATE DATA at line', IREC
                    CALL M3EXIT( 'GRDMMAT', 0, 0, MESG, 2 )
 
                ENDIF

            ENDIF

            C = COL  +  NCOLS * ( ROW - 1 )
            
C...........   Find this FIP and record this cell for it
            IF ( LFIP .GT. FIP ) THEN
 
                EFLAG = .TRUE.
                WRITE( MESG,94010 )
     &              'SURG COEFF file out of order at line', IREC
                CALL M3MESG( MESG )
 
            ELSEIF ( FIP .NE. LFIP ) THEN
                
                IF ( S .GE. MXFIP ) THEN
                    WRITE( MESG,94010 ) 
     &                  'Max number FIPs', MXFIP, 
     &                  'exceeded in SURG COEF file at line', IREC
                    CALL M3MESG( MESG )
                    EFLAG = .TRUE.
                    GO TO  70
                END IF
            
                LFIP = FIP
                S    = S + 1
                FIPLIST( S ) = FIP

            END IF      !  if fip not lfip
            
            K = NCELLS( S ) + 1
       
            IF ( K .GT. MXFIP ) THEN
                WRITE( MESG,94010 )
     &             'Max number FIPs', MXFIP, 
     &             'exceeded in SURG COEF file for FIP', FIP,
     &             'Cell', C,
     &             ' at line', IREC
                CALL M3MESG( MESG )
                EFLAG = .TRUE.
                GO TO  70
            
            ELSE

                NCELLS (   S ) = K
                FIPCELL( K,S ) = C
                DO  78  F = 1, NMSRG
                    SRGFRAC( F,K,S ) = FRAC( F )
78              CONTINUE

            END IF

            GO TO  70  !  to head of loop
      
79      CONTINUE        !  end of the SDEV-read loop

        IF ( EFLAG ) THEN
            CALL M3EXIT( 'GRDMMAT', 0, 0, 
     &                   'Error reading SURG COEF file.', 2 )
        END IF

        NFIPLIST = S


C...........   Read the surrogate cross-reference file

        WRITE( *,92000 ) 
     &      ' ', 'Reading GRIDDING SURROGATE XREF file...', ' '

        IREC  = 0
        I     = 0
        EFLAG = .FALSE.

88     CONTINUE        !  head of the XDEV-read loop

            READ( XDEV, *, END=89, IOSTAT=IOS ) FIP, SID, RDT
            IREC = IREC + 1

            IF ( IOS .NE. 0 ) THEN

                WRITE( MESG,94010 ) 
     &              'I/O Status ', IOS, 
     &              'reading SURG XREF file at line', IREC
                CALL M3MESG( MESG )
                EFLAG = .TRUE.

            ELSE IF ( I .LT. NMMAX ) THEN  !  dimension-check OK

                I = I + 1
                INDEXA( I ) = I
                SRFIPA( I ) = FIP
                SRGIDA( I ) = SID
                SRRDTA( I ) = RDT

            ELSE    !  overflow

                EFLAG = .TRUE.
                WRITE( MESG,94010 ) 
     &              'Max rec-cnt exceeded in XREF file at line', IREC
                CALL M3MESG( MESG )

            END IF

            GO TO  88


89      CONTINUE        !  end of the SDEV-read loop

        IF ( EFLAG ) THEN
            CALL M3EXIT( 'GRDMMAT', 0, 0, 
     &                   'Error reading SURG XREF file.', 2 )
        END IF

        NUMSCA = I


C.......   Sort and process the surrogate cross-reference file

        WRITE( *,92000 ) 
     &      ' ', 'Processing GRIDDING SURROGATE XREF file...', ' '
        
        CALL SORTI3( NUMSCA, INDEXA, SRFIPA, SRRDTA, SRLNKA )

        IREC  = 0
        I0    = 0
        I1    = 0
        I2    = 0
        EFLAG = .FALSE.
        SRGID0= IMISS3

        DO  99  I = 1, NUMSCA

            J   = INDEXA( I )
            FIP = SRFIPA( J )
            SID = SRGIDA( J )
            RDT = SRRDTA( J )

            IF ( FIP .EQ. 0 ) THEN  !  FIP-independent

                IF ( RDT .EQ. 0 ) THEN
                    SRGID0 = SID
                ELSE IF ( I0 .LT. NRCLAS ) THEN
                    I0 = I0 + 1
                    SRGIDR( I0 ) = SID 
                    SRRDTR( I0 ) = RDT
                ELSE
                    EFLAG = .TRUE.
                    WRITE( MESG,94010 )
     &                  'Max FIP-ind XREF exceeded for RDT', RDT
                    CALL M3MESG( MESG )
                END IF

            ELSE IF ( MOD( FIP, 1000 ) .EQ. 0 ) THEN   ! state-dependent case

                IF ( I1 .LT. MXSREF ) THEN  !  FIP-specific
                    I1 = I1 + 1
                    SRFIPS( I1 ) = FIP / 1000
                    SRGIDS( I1 ) = SID 
                    SRRDTS( I1 ) = RDT
                ELSE
                    EFLAG = .TRUE.
                    WRITE( MESG,94010 )
     &                  'Max STATE-dep XREF exceeded for FIP', FIP*1000,
     &                  'and RDT', RDT
                    CALL M3MESG( MESG )
                END IF

            ELSE        !  FIP > 0:  fip-dependent case

                IF ( I2 .LT. MXFREF ) THEN  !  FIP-specific
                    I2 = I2 + 1
                    SRFIPF( I2 ) = FIP  
                    SRGIDF( I2 ) = SID 
                    SRRDTF( I2 ) = RDT
                ELSE
                    EFLAG = .TRUE.
                    WRITE( MESG,94010 )
     &                  'Max FIP-dep XREF exceeded for FIP', FIP,
     &                  ' and RDT', RDT
                    CALL M3MESG( MESG )
                END IF

            END IF      !  if fip zero, or not

99      CONTINUE        !  end of the xref-processing loop

        IF ( EFLAG ) THEN
            CALL M3EXIT( 'GRDMMAT', 0, 0, 
     &                   'Error processing SURG XREF file.', 2 )
        END IF

        NUMSCR = I0
        NUMSCS = I1
        NUMSCF = I2


C.......   Read in the mobile source emissions tables:

        WRITE( *,92000 ) ' ', 'Reading in MOBILE SOURCES file...', ' '

        IF ( .NOT. READ3( ENAME, 'IFIP', ALLAYS3, 0,0, IFIP ) ) THEN
            CALL M3EXIT( 'GRDMMAT', 0, 0, 
     &      'Error reading variable "FIP" from MOBILE file', 2 )
        END IF

        IF ( .NOT. READ3( ENAME, 'IRCLAS', ALLAYS3,0,0, IRCLAS ) ) THEN
            CALL M3EXIT( 'GRDMMAT', 0, 0, 
     &      'Error reading variable "IRCLAS" from MOBILE file', 2 )
        END IF

        IF ( .NOT. READ3( ENAME, 'ILINK', ALLAYS3, 0,0, ILINK ) ) THEN
            CALL M3EXIT( 'GRDMMAT', 0, 0, 
     &      'Error reading variable "ILINK" from MOBILE file', 2 )
        END IF

        IF ( .NOT. READ3( ENAME, 'XLOC1', ALLAYS3, 0,0, XLOC1 ) ) THEN
            CALL M3EXIT( 'GRDMMAT', 0, 0, 
     &      'Error reading variable "XLOC1" from MOBILE file', 2 )
        END IF

        IF ( .NOT. READ3( ENAME, 'XLOC2', ALLAYS3, 0,0, XLOC2 ) ) THEN
            CALL M3EXIT( 'GRDMMAT', 0, 0, 
     &      'Error reading variable "XLOC2" from MOBILE file', 2 )
        END IF

        IF ( .NOT. READ3( ENAME, 'YLOC1', ALLAYS3, 0,0, YLOC1 ) ) THEN
            CALL M3EXIT( 'GRDMMAT', 0, 0, 
     &      'Error reading variable "YLOC1" from MOBILE file', 2 )
        END IF

        IF ( .NOT. READ3( ENAME, 'YLOC2', ALLAYS3, 0,0, YLOC2 ) ) THEN
            CALL M3EXIT( 'GRDMMAT', 0, 0, 
     &      'Error reading variable "YLOC2" from MOBILE file', 2 )
        END IF


C.......   Compute gridding matrix:
C.......       First case:   explicit link (ILINK > 0
C.......       Second case:  some LNKDEF entry applies
C.......       Third case:   FIP+roadtype cross-reference match
C.......       Fourth case:  state+roadtype cross-reference match
C.......       fifth case:   roadtype cross-reference match
C.......       sixth case:   fallback default

        WRITE( *,92000 ) ' ', 'Computing gridding matrix...', ' '

        OFLAG = .FALSE.         !  flag to detect overflow
        EFLAG = .FALSE.         !  flag to detect other errors
        LFIP  = 0
        NOS   = 0

        DO  155  I = 1, NMSRC
            
            FIP = IFIP  ( I )
            RDT = IRCLAS( I )
            LNK = ILINK ( I )
            
            IF ( LNK .GT. 0 ) THEN	!  explicit calculation for this link

C.................  Convert coordinate system from UTM to required system

                XBEG = XLOC1 ( I )
                YBEG = YLOC1 ( I )
                XEND = XLOC2 ( I )
                YEND = YLOC2 ( I )

                IF ( GDTYP3D .EQ. LATGRD3 ) THEN
                    CALL UTM2LL( XBEG, YBEG, UZONE,  XDUM, YDUM )
                    XBEG = XDUM
                    YBEG = YDUM

                    CALL UTM2LL( XEND, YEND, UZONE,  XDUM, YDUM )
                    XEND = XDUM
                    YEND = YDUM

                ELSEIF ( GDTYP3D .EQ. LAMGRD3 ) THEN

                    IF ( .NOT.  UTM2LAM( XBEG, YBEG, UZONE, 
     &                                   XDUM, YDUM ) ) THEN
                        EFLAG = .TRUE.
                        WRITE( MESG,94010 ) 
     &                  'Error converting <XBEG,YBEG> to Lambert ' //
     &                  'at line', IREC
                        CALL M3MESG( MESG )
                        GO TO  11     !  to head of loop:  read next line
                    END IF
                    XBEG = XDUM
                    YBEG = YDUM

                    IF ( .NOT.  UTM2LAM( XEND, YEND, UZONE, 
     &                                   XDUM, YDUM ) ) THEN
                        EFLAG = .TRUE.
                        WRITE( MESG,94010 ) 
     &                  'Error converting <XEND,YEND> to Lambert ' //
     &                  'at line', IREC
                        CALL M3MESG( MESG )
                        GO TO  11     !  to head of loop:  read next line
                    END IF
                    XEND = XDUM
                    YEND = YDUM

                END IF

                CALL LNK2GRD( XBEG, YBEG, XEND, YEND,
     &                        NCEL, ACEL, AFAC, ALEN )

                IF ( NCEL .GT. MXCFIP ) THEN
                    EFLAG = .TRUE.
                    WRITE( MESG,94010 )
     &                  'LNK:CEL table overflow for source', S
                    CALL M3MESG( MESG )
                    GO TO  155
                END IF

                IF( ALEN .NE. 0 ) THEN
                    DD = 1.0 / ALEN

                ELSE
                    EFLAG = .TRUE.
                    WRITE( MESG,94010 )
     &                  'Bad coordinates for source', I
                    CALL M3MESG( MESG )
                    GO TO  155
                    
                ENDIF

                DO  111  K = 1, NCEL
                    C = ACEL( K )
                    J = NX( C ) + 1
                    IF ( J .GT. MXSCEL ) THEN
                        OFLAG = .TRUE.      !  overflow occurred.
                        WRITE( MESG,94010 )
     &                      'Overflow for cell', C
                        CALL M3MESG( MESG )
                    ELSE  
                        IS( J,C ) = I 
                        CS( J,C ) = DD * AFAC( K )
                    END IF
                    NX( C )   = J
111             CONTINUE

                GO TO 155	!  to end of loop

            END IF
            
            S = FIND2( FIP, RDT, NUMLNK, LNKFIP, LNKRDT )  !  use LNKDEF factor
            IF ( S .GT. 0 ) THEN
                
                DO  122  K = 1, LNKCNT( S )
                    C = LNKCEL( J,S )
                    J = NX( C ) + 1
                    IF ( J .GT. MXSCEL ) THEN
                        OFLAG = .TRUE.      !  overflow occurred.
                        WRITE( MESG,94010 )
     &                      'Overflow for cell', C
                        CALL M3MESG( MESG )
                    ELSE  
                        IS( J,C ) = I 
                        CS( J,C ) = LNKFAC( K,S )
                    END IF
                    NX( C )   = J
122             CONTINUE

                GO TO 155       !  to end of loop

            END IF

            S = FIND2( FIP, RDT, NUMSCF, SRFIPF, SRRDTF )  !  FIP+roadtype
            IF ( S .GT. 0 ) THEN
                
                SID = SRGIDF( S )
                GO TO 133       !  to XREF processing

            END IF

            S = FIND2( FIP, RDT, NUMSCS, SRFIPS, SRRDTS )  !  state+roadtype
            IF ( S .GT. 0 ) THEN
                
                SID = SRGIDS( S )
                GO TO 133       !  to XREF processing

            END IF

            S = FIND1( RDT, NUMSCR, SRRDTR )  !  roadtype only
            IF ( S .GT. 0 ) THEN
                
                SID = SRGIDR( S )
                GO TO 133       !  to XREF processing

            ELSE                !  ultimate fallback
                
                SID = SRGID0
                IF ( SID .LT. 0 ) THEN
                    EFLAG = .TRUE.
                    WRITE( MESG,94010 ) 
     &                  'No surrogate XREF found for FIP', FIP,
     &                  'road type', RDT
                    CALL M3MESG( MESG )
                    GO TO  155
                END IF
            END IF    

C...........   Find all cells intersecting this FIPS and fill in 
C...........   corresponding gridding-matrix coefficients.
C...........   ASSUMPTION:  list is sorted by FIP (and then by col-row),
C...........   with multiple entries per FIP.
            
133         CONTINUE    !  Process XREFs

C...........   Check SID against maximum to make sure no overflow
            IF( SID .GT. NMSRG ) THEN
                EFLAG = .TRUE.
                WRITE( MESG,94010 ) 
     &                 'Surrogate index', SID, 
     &                 'used in MGREF, but maximum (NMSRG) is', NMSRG
                CALL M3MESG( MESG )
                GO TO  155

            ENDIF

            IF( FIP .NE. LFIP ) THEN

                LFIP = FIP
                F = FIND1( FIP, NFIPLIST, FIPLIST )

C.................  Report no surrogates found            
                IF ( F .LT. 0 ) THEN

                    WRITE( MESG,94010 ) 
     &                     'No surrogates for FIP code ', FIP
                    CALL M3MESG( MESG )

C.....................  Store list of unfound FIPs
                    NOS = NOS + 1
                    NOSURG( NOS ) = FIP

                    GO TO 155

                END IF

            ELSEIF( F .LT. 0 ) THEN
                GO TO 155

            ENDIF

            DO  144  K = 1, NCELLS( F )
            
                C = FIPCELL( K,F )
                J = NX( C ) + 1

                IF ( J .GT. MXSCEL ) THEN
                    OFLAG = .TRUE.      !  overflow occurred.
                    WRITE( MESG,94010 )  
     &              'Overflow for FIP  code ', FIP,
     &              'RDT code', RDT,
     &              'cell code', C, ' '
                    CALL M3MESG( MESG )
                ELSE
                    IS( J,C ) = I
                    CS( J,C ) = SRGFRAC( SID,K,F )
                END IF
                NX( C )   = J

144         CONTINUE    !  end of loop on cells K for this FIP

155     CONTINUE        !  end loop on sources I, computing gridding matrix.
        
        CMAX = NX( 1 )
        CMIN = CMAX
        
        IF ( OFLAG ) THEN   !  overflow occurred:  just compute statistics.
            
            DO  211  R = 2, NGRID
            
                J = NX( R )
                IF      ( J .GT. CMAX ) THEN
                    CMAX = J
                ELSE IF ( J .LT. CMIN ) THEN
                    CMIN = J
                END IF

211         CONTINUE    !  end of loop on cells K for this FIP

            WRITE( *,92000 )
     &          ' ',
     &          'WARNING:  ' // GNAME( 1:TRIMLEN( GNAME ) ) // 
     &          ' not written.',
     &          'Internal arrays would have overflowed.'
            WRITE( *,92010 ) 
     &          'Sources per cell dimensioned', MXSCEL

        ELSE IF ( EFLAG ) THEN          ! error in the XREF lookup

            CALL M3EXIT( 'GRDMMAT', 0, 0, 'Missing FIP codes', 3 ) 

        ELSE IF ( GFLAG ) THEN          ! write gridding matrix to file

C.......   Compress matrix into I/O representation from scratch representation.
C.......   Compute statistics.

            K = 0
                
            DO 222  C = 1, NX( 1 )
                K       = K + 1
                IX( K ) = IS( C,1 )
                CX( K ) = CS( C,1 )
222         CONTINUE
                   
            DO  233  R = 2, NGRID
            
                J = NX( R )
                   
                IF (      J .GT. CMAX ) THEN
                    CMAX = J
                ELSE IF ( J .LT. CMIN ) THEN
                    CMIN = J
                END IF
                   
                DO 232  C = 1, J
                    K       = K + 1
                    IF ( K .LE. NMATX ) THEN
                        IX( K ) = IS( C,R )
                        CX( K ) = CS( C,R )
                    END IF
232             CONTINUE
                   
233         CONTINUE    !  end of loop on cells K for this FIP
            
            IF ( K .GT. NMATX ) THEN
                WRITE( MESG,94010 ) 
     &                  'Number of coeffs', K, 
     &                  'exceeds NMATX=', NMATX
                CALL M3EXIT( 'GRDMMAT', 0, 0, MESG, 3 ) 
            END IF
                   
C.......   Write out the gridding matrix:
            
            WRITE( *,92000 ) ' ', 'Writing out GRIDDING MATRIX file...'

            IF ( .NOT. WRITE3( GNAME, 'ALL', 0, 0, NX ) ) THEN
                CALL M3EXIT( 'GRDMMAT', 0, 0, 
     &              'Error writing GRIDDING MATRIX file.', 2 )
            END IF

        ELSE    !  not GFLAG:  report matrix stats only
        
C.......   Just compute statistics.
            
            K = 0
                
            DO  244  R = 2, NGRID
            
                J = NX( R )
                K = K + J
                   
                IF      ( J .GT. CMAX ) THEN
                    CMAX = J
                ELSE IF ( J .LT. CMIN ) THEN
                    CMIN = J
                END IF

244         CONTINUE    !  end of loop on cells K for this FIP

        END IF  !  if gflag, or not
                        
C.......   Report statistics:

        XDUM = FLOAT( K ) / FLOAT( NGRID )
        WRITE( *,92000 ) ' ', 'GRIDDING-MATRIX statistics', ' '
        WRITE( *,92010 ) 
     &      'Total number of coefficients   ', K   ,
     &      'Max  number of sources per cell', CMAX,
     &      'Min  number of sources per cell', CMIN
        WRITE( *,92020 ) 
     &      'Mean number of sources per cell', XDUM,
     &      ' '

C...........   IF UFLAG and no errors, generate the UNgridding matrix:
        
        IF (       OFLAG ) GO TO  999	!  to exit from program
        IF (       EFLAG ) GO TO  999	!  to exit from program
        IF ( .NOT. UFLAG ) GO TO  999	!  to exit from program
        
C.......   Initialize ungridding matrix related arrays
        DO  301  S = 1, NMSRC
            NU( S ) = 0
            DN( S ) = 0.
301     CONTINUE

C.......   Create renormalization factors for each source
        DO  311  C = 1, NGRID
            DO 310 K = 1, NX( C )

                S = IS( K,C )
                DN( S ) = DN( S ) + CS( K,C )

310         CONTINUE
311     CONTINUE

C.......   Pre-divide normalization factors
        DO  312 S = 1, NMSRC

C.............  Look for FIP in list of FIPs w/ no surrogates
            F = FIND1( IFIP( S ), NOS, NOSURG )

C.............  Check for no gridding factors and make sure source
C               hasn't already been flagged as having no surrogates
            IF( DN( S ) .EQ. 0. .AND. F .LE. 0 ) THEN

                WRITE( MESG,94010 ) 
     &                'No gridding factors. FIP:', IFIP( S ),
     &                'RDT:', IRCLAS( S ), 'LNK:', ILINK( S )
                CALL M3WARN( 'GRDMMAT', 0, 0, MESG )

            ELSEIF( F .LE. 0 ) THEN
                DN( S ) = 1.0 / DN( S )

            ENDIF

312     CONTINUE

C.......   Renormalize and Transpose ((IS,CS) into scratch arrays (IT,CT):
        DO  333  C = 1, NGRID
            
            DO  322  K = 1, NX( C )

                S = IS( K,C )
                J = NU( S ) + 1

                IF ( J .LE. MXCFIP ) THEN
                    IT( J,S ) = C
                    CT( J,S ) = DN( S ) * CS( K,C )
                ELSE
                    OFLAG = .TRUE.
                    WRITE( MESG, 94010 ) 'Overflow for source', S
                    CALL M3MESG( MESG )
                END IF

                NU( S ) = J

322         CONTINUE

333     CONTINUE


        CMAX = NU( 1 )
        CMIN = CMAX

        IF ( OFLAG ) THEN		!  just compute statistics

            K = CMIN
            DO  344  S = 2, NMSRC
                J = NU( S )
                K = K + J
                IF ( J .GT. CMAX ) THEN
                    CMAX = J
                ELSE IF ( J .LT. CMIN ) THEN
                    CMIN = J
                END IF
344         CONTINUE        

            WRITE( *,92000 )
     &          ' ',
     &          'WARNING:  ' // UNAME( 1:TRIMLEN( UNAME ) ) // 
     &          ' not written.',
     &          'Internal arrays would have overflowed.'
            WRITE( *,92010 ) 
     &          'Cells per source dimensioned', MXCFIP

        ELSE
        
C.......   Compress matrix into I/O representation from scratch representation.
C.......   Compute statistics.
          
            K  = 0
            
            DO  355  I = 1, NU( 1 )
                K = K + 1
                IU( K ) = IT( I,1 )
                CU( K ) = CT( I,1 )
355         CONTINUE
                   
            DO  377  S = 2, NMSRC

                J = NU( S )
                IF ( J .GT. CMAX ) THEN
                    CMAX = J
                ELSE IF ( J .LT. CMIN ) THEN
                    CMIN = J
                END IF
                
                DO  366 I = 1, J
                    K = K + 1
                    IU( K ) = IT( I,S )
                    CU( K ) = CT( I,S )
366             CONTINUE

377         CONTINUE

            WRITE( *,92000 ) ' ','Writing out UNGRIDDING MATRIX file...'

            IF ( .NOT. WRITE3( UNAME, 'ALL', 0, 0, NU ) ) THEN
                CALL M3EXIT( 'GRDMMAT', 0, 0, 
     &              'Error writing UNGRIDDING MATRIX file.', 2 )
            END IF

        END IF		!  if oflag:  overflow for ungridding matrix, or not


C.......   Report statistics for ungridding matrix:

        XDUM = FLOAT( K ) / FLOAT( NMSRC )
        WRITE( *,92000 ) ' ', 'UNGRIDDING-MATRIX statistics', ' '
        WRITE( *,92010 ) 
     &      'Total number of coefficients   ', K   ,
     &      'Max  number of cells per source', CMAX,
     &      'Min  number of cells per source', CMIN
        WRITE( *,92020 ) 
     &      'Mean number of cells per source', XDUM,
     &      ' '

                  
C...............   End of program

999     CONTINUE
      
        CALL M3EXIT( 'GRDMMAT', 0, 0, 
     &               'SUCCESSFUL COMPLETION of program GRDMMAT', 0 )


C******************  FORMAT  STATEMENTS   ******************************

C...........   Informational (LOG) message formats... 92xxx

92000   FORMAT( 5X, A )

92010   FORMAT( 5X, A, :, I12 )

92020   FORMAT( 5X, A, :, F17.4 )


C...........   Formatted file I/O formats............ 93xxx

93010   FORMAT( 4 F10.0, I10, 2F10.0 )

93020   FORMAT( I5, 2F10.0, I10, 30F10.8 )

93030   FORMAT( I5, I4, I8, I3 )        !  for surrogate XREF file


C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10 ( A, :, I8, :, 2X ) )


        END

