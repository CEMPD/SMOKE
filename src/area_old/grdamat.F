
C Version "@(#)$Id$ $Source$ $Date$ 

        PROGRAM GRDAMAT

C***********************************************************************
C  program body starts at line 171
C
C  DESCRIPTION:
C       Construct area source gridding matrix from data contained in
C       EPS-style surrogates file.
C
C  PRECONDITIONS REQUIRED:
C       Sorted, cut-down input data for surrogate coeffs.
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C       Models-3 I/O
C       FIND1, FIND2, FIND3, GETYN, PROMPTFFILE, PROMPTMFILE, TRIMLEN
C
C  REVISION  HISTORY:
C       Prototype  2/95 by CJC.
C
C       Version   11/95 by CJC sorts GREF on the fly; more sophisticated
C       error trapping.
C
C       Version   1/96 by CJC reads FIP-X-Y surrogate coeff files
C***********************************************************************

      IMPLICIT NONE

C...........   INCLUDES:
        
        INCLUDE 'ARDIMS3.EXT'   !  area-source dimensioning parameters
        INCLUDE 'CHDIMS3.EXT'   !  emis chem parms (inventory + model)
        INCLUDE 'GRDIMS3.EXT'   !  grid parameters
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
        INCLUDE 'FDESC3.EXT'    !  I/O API file description data structures.


C...........   EXTERNAL FUNCTIONS and their descriptions:
        
        INTEGER         FIND1, FIND2, FIND3
        LOGICAL         GETYN
        INTEGER         PROMPTFFILE
        CHARACTER*16    PROMPTMFILE
        LOGICAL         DSCGRID
        INTEGER         STR2INT
        INTEGER         TRIMLEN
        
        EXTERNAL        DSCGRID, FIND1, FIND2, FIND3, GETYN, 
     &                  PROMPTFFILE, PROMPTMFILE, STR2INT, TRIMLEN


C...........   PARAMETER:

        INTEGER         NAMAX
        PARAMETER     ( NAMAX = NASCC + NASRC )

C...........   LOCAL VARIABLES and their descriptions:
C...........   NOTE that ASC (area-source-category) ID's are 10-digit
C...........   unsigned integers which may be treated as a leading 7-digit
C...........   field, and a trailing 3-digit field.  *7 and *3 arrays below
C...........   follow this scheme with parallel arrays

        INTEGER     FIP
        INTEGER     SID
        INTEGER     ID7
        INTEGER     ID3
        INTEGER     LFIP, LID7, LID3
        INTEGER     COL
        INTEGER     ROW
        INTEGER     ZONE
        INTEGER     CMAX, CMIN
        REAL        XDUM, YDUM, SDUM, TDUM, UDUM, VDUM
        REAL        DDX, DDY
        REAL        FRAC( NASRG )         !  (scratch) surrogate fractions

        CHARACTER*7 ID7C
        CHARACTER*3 ID3C
 
C.......   Actually-occurring ASC table (used for filtering out 
C.......   irrelevant xrefs

        INTEGER     NASC
        INTEGER     ASCA7( NASCC )
        INTEGER     ASCA3( NASCC )

C...........   XREF before sorting:                                 

        INTEGER     NUMSCA           !  current surrogates-table count
        INTEGER     INDEXA( NAMAX )  !  FIP associated with this surrogate
        INTEGER     SRFIPA( NAMAX )  !  FIP associated with this surrogate
        INTEGER     SRGIDA( NAMAX )  !  surrogate munber (1...15)
        INTEGER     ASCI7A( NAMAX )  !  leading-7  digits of ASC     
        INTEGER     ASCI3A( NAMAX )  !  trailing-3 digits of ASC

C...........   XREF after sorting and processing: Surrogates::ASC tables 
C...........   (TWO tables, one for FIP==0; the other for FIP-specific data.

        INTEGER     NUMSC0          !  current surrogates-table count
        INTEGER     SRGID0( NASCC ) !  surrogate munber (1...15)
        INTEGER     ASCI70( NASCC ) !  leading-7  digits of ASC     
        INTEGER     ASCI30( NASCC ) !  trailing-3 digits of ASC

        INTEGER     NUMSCC          !  current surrogates-table count
        INTEGER     SRFIP( NASRC )  !  FIP associated with this surrogate
        INTEGER     SRGID( NASRC )  !  surrogate munber (1...15)
        INTEGER     ASCI7( NASRC )  !  leading-7  digits of ASC     
        INTEGER     ASCI3( NASRC )  !  trailing-3 digits of ASC

C...........   Surrogate-cell::FIPS table
        
        INTEGER     NFIPLIST                        !  # of FIPS
        INTEGER     NCELLS (               NAFIP )  !  # cells per FIP
        INTEGER     FIPLIST(               NAFIP )  !  FIP-list
        INTEGER     FIPCELL(        NCFIP, NAFIP )  !  cty codes::cell number
        REAL        SRGFRAC( NASRG, NCFIP, NAFIP )  !  surrogate fractions
                
C...........   Area Sources Table
        
        INTEGER     IFIPS( NASRC )  !  source FIPS (county) ID
        INTEGER     ISCC7( NASRC )  !  leading-7  digits of source ASC
        INTEGER     ISCC3( NASRC )  !  trailing-3 digits of source ASC
        
C...........   Scratch Gridding Matrix (subscripted by source-within-cell, cell)

        INTEGER     IS( NACEL , NGRID )
        REAL        CS( NACEL , NGRID )

C...........   Gridding Matrix

        INTEGER     NX( NGRID )
        INTEGER     IX( NMATX )
        REAL        CX( NMATX )

        COMMON  / GRIDMAT / NX, IX, CX

C...........   Other local variables
        
        INTEGER         C, R, S, F, I, J, K  !  loop counters.

        CHARACTER*16    ENAME   !  logical name for emission source  input file
        CHARACTER*16    GNAME   !  logical name for grid matrix     output file
        INTEGER         ADEV    !  for actual-ASC file
        INTEGER         SDEV    !  for surrogate coeff file
        INTEGER         XDEV    !  for surrogate xref  file
        INTEGER         LDEV    !  log-device
        
        INTEGER         IOS     !  I/O status
        INTEGER         IREC    !  input line (record) number
        INTEGER         I0      !  for counting FIP-independent srcs
        INTEGER         I1      !  for counting FIP-dependent   srcs

        LOGICAL         GFLAG   !  generate output gridding matrix if true
        LOGICAL         OFLAG   !  Did overflow occur
        CHARACTER*16    SCRBUF
        CHARACTER*256   MESG

        LOGICAL         EFLAG   !  input error flat


C...........   STATEMENT FUNCTIONS:  
C.......   floating point "unequal" -- true iff 
C.......   | P - Q | > 1e-5 * sqrt( p*p + q*q + 1e-5 )

        REAL            P, Q
        LOGICAL         FLTERR
        FLTERR( P, Q ) = 
     &  ( (P - Q)**2  .GT.  1.0E-10*( P*P + Q*Q + 1.0E-5 ) )

C***********************************************************************
C   begin body of program GRDAMAT
        
        LDEV = INIT3()
        
        WRITE( *,92000 ) 
     &  ' ',
     &'Program GRDAMAT to take the sorted EPS surrogate coefficient',
     &'file, the EPS surrogate cross reference file, the AREA SOURCE',
     &'file produced by RAWAREA, and the ACTUAL-ASC file, and produce',
     &'the AREA SOURCE GRIDDING MATRIX file.', 
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
            CALL M3EXIT( 'GRDAMAT', 0, 0, 'Ending program GRDAMAT', 2 )
        END IF


C.......   Get file name; open actual-scc file
        
        ADEV = PROMPTFFILE( 
     &           'Enter logical name for ACTUAL ASC file',
     &           .TRUE., .TRUE., 'ASCC', 'GRDAMAT' )

        XDEV = PROMPTFFILE( 
     &  'Enter logical name for GRIDDING SURROGATE XREF file',
     &           .TRUE., .TRUE., 'AGREF', 'GRDAMAT' )

        SDEV = PROMPTFFILE( 
     &  'Enter logical name for GRIDDING SURROGATE COEFF file',
     &           .TRUE., .TRUE., 'AGPRO', 'GRDAMAT' )

        ENAME = PROMPTMFILE( 
     &          'Enter logical name for AREA SOURCE input file',
     &          FSREAD3, 'AREA', 'GRDAMAT' )


C.......   Initializations:
        
        NUMSCC = 0      !  for Surrogates::ASC table

        DO  11  C = 1, NGRID
            NX( C ) = 0         !  for Gridding Matrix
11      CONTINUE

        DO  12  S = 1, NAFIP
            NCELLS( S ) = 0     !  for Surrogate  cell::FIPS table
12      CONTINUE

cmh add: Get grid description for use in check

        IF ( .NOT. DSCGRID( GRDNM, SCRBUF, GDTYP3D, 
     &              P_ALP3D, P_BET3D, P_GAM3D, XCENT3D, YCENT3D,
     &              XORIG3D, YORIG3D, XCELL3D, YCELL3D, 
     &              NCOLS3D, NROWS3D, NTHIK3D ) ) THEN
            SCRBUF = GRDNM
            MESG   = '"' // SCRBUF( 1:TRIMLEN( SCRBUF ) ) //
     &               '" not found in GRIDDESC file'
            CALL M3EXIT( 'GRDAMAT', 0, 0, MESG, 2 )

        END IF          ! if dscgrid() failed

C.......   Read the surrogate coefficient file
C...............   ASSUMPTION:  sorted by FIP
        
        WRITE( *,92000 ) 
     &  ' ', 'Reading GRIDDING SURROGATE COEFFICIENT file...', ' '

        READ( SDEV, *, IOSTAT=IOS ) !  read file header
     &          SDUM, TDUM, 
     &          UDUM, VDUM, ZONE,
     &          XDUM, YDUM
        
        IF ( IOS .NE. 0 ) THEN
            WRITE( MESG,94010 ) 
     &          'I/O error ', IOS, 
     &          'reading SURG COEFF file header at line', IREC
            CALL M3EXIT( 'GRDAMAT', 0, 0, MESG, 1 )
        END IF

C.......   Check grid compatibility:

        DDX  = 1.0 / XCELL3D
        DDY  = 1.0 / YCELL3D
        SDUM = DDX * ( SDUM - XORIG3D )
        TDUM = DDY * ( TDUM - YORIG3D )
        IF ( FLTERR( XDUM, XCELL3D ) ) THEN
            CALL M3EXIT( 'GRDAMAT', 0, 0, 
     &                   'Bad XCELL in SURROGATE file', 2 )
        ELSE IF ( FLTERR( YDUM, YCELL3D ) ) THEN
            CALL M3EXIT( 'GRDAMAT', 0, 0, 
     &                   'Bad YCELL in SURROGATE file', 2 )
        ELSE IF ( FLTERR( SDUM, FLOAT( NINT( SDUM ) ) ) ) THEN
            CALL M3EXIT( 'GRDAMAT', 0, 0, 
     &                   'Bad X alignment in SURROGATE file', 2 )
        ELSE IF ( FLTERR( TDUM, FLOAT( NINT( TDUM ) ) ) ) THEN
            CALL M3EXIT( 'GRDAMAT', 0, 0, 
     &                   'Bad Y alignment in SURROGATE file', 2 )
        END IF
                           
        IREC     =  1
        NFIPLIST =  0
        S        =  0
        LFIP     = -1
        EFLAG    = .FALSE.
        
22      CONTINUE        !  head of the SDEV-read loop
        
            IREC = IREC + 1
            READ( SDEV, *, END=24, IOSTAT=IOS ) 
     &          FIP, XDUM, YDUM, ZONE, ( FRAC( C ), C = 1, NASRG )

            IF ( IOS .NE. 0 ) THEN

                WRITE( MESG,94010 ) 
     &              'I/O error ', IOS, 
     &              'reading SURG COEFF file at line', IREC
                CALL M3MESG( MESG )
                EFLAG = .TRUE.
                GO TO  22

            END IF

            COL = NINT( DDX * ( XDUM - XORIG3D ) )  +  1
            ROW = NINT( DDY * ( YDUM - YORIG3D ) )  +  1
            IF ( COL .LT. 1  .OR. COL .GT. NCOLS ) THEN

                WRITE( MESG,94010 )
     &              'Column ', COL,
     &              'out of range in SURG COEFF file at line', IREC
                CALL M3MESG( MESG )
                EFLAG = .TRUE.
                GO TO  22

            ELSE IF ( ROW .LT. 1  .OR. ROW .GT. NROWS ) THEN

                WRITE( MESG,94010 )
     &              'Row ', ROW,
     &              'out of range in SURG COEFF file at line', IREC
                CALL M3MESG( MESG )
                EFLAG = .TRUE.
                GO TO  22

            END IF

            C = COL  +  NCOLS * ( ROW - 1 )
            

C...........   Find this FIP and record this cell for it
        
            IF ( LFIP .GT. FIP ) THEN

                EFLAG = .TRUE.
                WRITE( MESG,94010 )
     &              'SURG COEFF file out of order at line', IREC
                CALL M3MESG( MESG )

            ELSEIF( FIP .NE. LFIP ) THEN
                
                IF ( S .GE. NAFIP ) THEN
                    WRITE( MESG,94010 ) 
     &                  'FIP max', NAFIP, 
     &                  'exceeded in SURG COEF file at line', IREC
                    CALL M3MESG( MESG )
                    EFLAG = .TRUE.
                    GO TO  22
                END IF
            
                LFIP = FIP
                S    = S + 1
                FIPLIST( S ) = FIP

            END IF      !  if fip not lfip
            
            K = NCELLS( S ) + 1
       
            IF ( K .GT. NCFIP ) THEN
                WRITE( MESG,94010 )
     &             'Cell:FIP max', NCFIP, 
     &             'exceeded in SURG COEF file for FIP', FIP,
     &             'Cell', C,
     &             ' at line', IREC
                CALL M3MESG( MESG )
                EFLAG = .TRUE.
                GO TO  22
            
            ELSE

                NCELLS (   S ) = K
                FIPCELL( K,S ) = C
                DO  23  F = 1, NASRG
                    SRGFRAC( F,K,S ) = FRAC( F )
23              CONTINUE

            END IF

            GO TO  22  !  to head of loop
      
24      CONTINUE        !  end of the SDEV-read loop

        IF ( EFLAG ) THEN
            CALL M3EXIT( 'GRDAMAT', 0, 0, 
     &                   'Error reading SURG COEF file.', 2 )
        END IF

        NFIPLIST = S


C.......   Read the ACTUAL ASC FILE
        
        WRITE( *,92000 ) 
     &       ' ', 'Reading ACTUAL ASC file...', ' '

        IREC  =  0
        I     =  0
        NASC  =  0
        LID3  = -1
        LID7  = -1
        EFLAG = .FALSE.
        
33      CONTINUE        !  head of the SDEV-read loop
        
            READ( ADEV, 93020, END=34, IOSTAT=IOS ) ID7, ID3

            IREC = IREC + 1

            IF ( IOS .NE. 0 ) THEN

                EFLAG = .TRUE.
                WRITE( MESG,94010 ) 
     &              'I/O error', IOS, 
     &              'reading ACTUAL ASC file at line', IREC
                CALL M3MESG( MESG )

            ELSE IF ( I .GE. NASCC ) THEN 

                EFLAG = .TRUE.
                WRITE( MESG,94010 )
     &              'Max ASC TABLE size', NASCC, 
     &              'exceeded at line', IREC
                CALL M3MESG( MESG )

            ELSE IF ( ( LID7 .GT. ID7 ) .OR.
     &                ( LID7 .EQ. ID7 .AND. LID3 .GE. ID3 ) ) THEN 

                EFLAG = .TRUE.
                WRITE( MESG,94010 )
     &              'ASC table out of order at line', IREC
                CALL M3MESG( MESG )

            ELSE

                I          = I + 1
                LID7       = ID7
                LID3       = ID3
                ASCA7( I ) = ID7
                ASCA3( I ) = ID3

            END IF

            GO TO  33  !  to head of loop
      
34      CONTINUE        !  end of the SDEV-read loop

        IF ( EFLAG ) THEN
            CALL M3EXIT( 'GRDAMAT', 0, 0, 
     &                   'Error reading ACTUAL ASC file.', 2 )
        END IF

        NASC = I


C.......   Read the surrogate cross-reference file

        WRITE( *,92000 ) 
     &      ' ', 'Reading GRIDDING SURROGATE XREF file...', ' '

        IREC  = 0
        I     = 0
        EFLAG = .FALSE.

44     CONTINUE        !  head of the XDEV-read loop

            READ( XDEV, 93030, END=45, IOSTAT=IOS ) FIP, SID, ID7C, ID3C
            IREC = IREC + 1

            IF ( IOS .NE. 0 ) THEN

                WRITE( MESG,94010 ) 
     &              'I/O Status ', IOS, 
     &              'reading SURG XREF file at line', IREC
                CALL M3MESG( MESG )
                EFLAG = .TRUE.

            ELSE IF ( ID7C(1:1) .GT. '9' ) THEN  ! Motor vehicles x-ref entry

                GO TO 44

            ELSE

                ID7 = STR2INT( ID7C )
                ID3 = STR2INT( ID3C )

                IF( ID7 .GT. 0 .AND. ID3 .GE. 0 ) THEN
                    J = FIND2( ID7, ID3, NASC, ASCA7, ASCA3 )
                ELSE
                    J = 1
                ENDIF

                IF ( J .LE. 0 ) THEN  ! ASC code not found in ASCC LIST file

                    GO TO 44

                ELSEIF ( I .GE. NAMAX ) THEN  !  Overflow

                    WRITE( MESG,94010 ) 
     &                     'XREF table overflow at record', IREC
                    CALL M3MESG( MESG )
                    EFLAG = .TRUE.

                ELSE 

                    I = I + 1
                    INDEXA( I ) = I
                    SRFIPA( I ) = FIP
                    SRGIDA( I ) = SID
                    ASCI7A( I ) = ID7
                    ASCI3A( I ) = ID3

                END IF  !  No overflow, ASC found in list

            END IF      !  Not motor vehicles record

        GO TO  44       !  to head of XDEV-read loop

45      CONTINUE        !  end of the XDEV-read loop

        IF ( EFLAG ) THEN
            CALL M3EXIT( 'GRDAMAT', 0, 0, 
     &                   'Error reading SURG XREF file.', 2 )
        END IF

        NUMSCA = I


C.......   Sort and process the surrogate cross-reference file

        WRITE( *,92000 ) 
     &      ' ', 'Processing GRIDDING SURROGATE XREF file...', ' '
        
        CALL SORTI3( NUMSCA, INDEXA, SRFIPA, ASCI7A, ASCI3A )

        IREC  = 0
        I0    = 0
        I1    = 0
        EFLAG = .FALSE.

        DO  55  I = 1, NUMSCA

            J   = INDEXA( I )
            FIP = SRFIPA( J )
            SID = SRGIDA( J )
            ID7 = ASCI7A( J )
            ID3 = ASCI3A( J )

            IF ( FIP .EQ. 0 ) THEN  !  FIP-independent

                IF ( I0 .LT. NASCC ) THEN
                    I0 = I0 + 1
                    SRGID0( I0 ) = SID 
                    ASCI70( I0 ) = ID7 
                    ASCI30( I0 ) = ID3
                ELSE
                    EFLAG = .TRUE.
                    WRITE( MESG,94010 )
     &                  'Max FIP-ind XREF exceeded for ID7', ID7,
     &                  'and ID3', ID3
                    CALL M3MESG( MESG )
                END IF

            ELSE        !  FIP > 0:  fip-dependent case

                IF ( I1 .LT. NASRC ) THEN  !  FIP-specific
                    I1 = I1 + 1
                    SRFIP( I1 ) = FIP  
                    SRGID( I1 ) = SID 
                    ASCI7( I1 ) = ID7 
                    ASCI3( I1 ) = ID3
                ELSE
                    EFLAG = .TRUE.
                    WRITE( MESG,94010 )
     &                  'Max FIP-dep XREF exceeded for FIP', FIP,
     &                  'ID7', ID7,
     &                  'and ID3', ID3
                    CALL M3MESG( MESG )
                END IF

            END IF      !  if fip zero, or not

55      CONTINUE        !  end of the xref-processing loop

        IF ( EFLAG ) THEN
            CALL M3EXIT( 'GRDAMAT', 0, 0, 
     &                   'Error processing SURG XREF file.', 2 )
        END IF

        NUMSC0 = I0
        NUMSCC = I1


C.......   Get the grid name and parameters
C.......   to put into file header:

        IF ( .NOT. DSCGRID( GRDNM, SCRBUF, GDTYP3D, 
     &              P_ALP3D, P_BET3D, P_GAM3D, XCENT3D, YCENT3D,
     &              XORIG3D, YORIG3D, XCELL3D, YCELL3D, 
     &              NCOLS3D, NROWS3D, NTHIK3D ) ) THEN

            SCRBUF = GRDNM
            MESG   = '"' // SCRBUF( 1:TRIMLEN( SCRBUF ) ) //
     &               '" not found in GRIDDESC file'
            CALL M3EXIT( 'GRDAMAT', 0, 0, MESG, 2 )

        END IF          ! if dscgrid() failed


C.......   Get file name; open output gridding matrix file

        FTYPE3D = SMATRX3
        SDATE3D = 0
        STIME3D = 0
        TSTEP3D = 0
        NVARS3D = 1
        NCOLS3D = NMATX
        NROWS3D = NGRID
        NLAYS3D = 1
        NTHIK3D = NASRC
        GDNAM3D = GRDNM
        VGTYP3D = IMISS3
        VGTOP3D = BADVAL3
        FDESC3D( 1 ) = 'NC Area Source gridding-coefficient matrix'
        DO  167  K = 2, MXDESC3
            FDESC3D( K ) = ' '
167     CONTINUE
        VNAME3D( 1 ) = 'AGRDMAT'
        UNITS3D( 1 ) = 'n/a'
        VDESC3D( 1 ) = 'NC Area Source gridding-coefficient matrix'
        VTYPE3D( 1 ) = M3REAL
           
        GNAME = PROMPTMFILE( 
     &          'Enter logical name for GRIDDING MATRIX output file',
     &          FSUNKN3, 'AGMAT', 'GRDAMAT' )

        GFLAG = ( GNAME( 1:5 ) .NE. 'NONE ' )


C.......   Read in the area source emissions FIP and ASC tables:
        
        WRITE( *,92000 ) ' ', 'Reading in AREA SOURCES file...', ' '

        IF ( .NOT. READ3( ENAME, 'FIP', ALLAYS3, 0, 0,  IFIPS ) ) THEN
            CALL M3EXIT( 'GRDAMAT', 0, 0, 
     &      'Error reading variable "FIP" from AREA SOURCES file', 2 )
        END IF

        IF ( .NOT. READ3( ENAME, 'ASC7', ALLAYS3, 0, 0,  ISCC7 ) ) THEN
            CALL M3EXIT( 'GRDAMAT', 0, 0, 
     &      'Error reading variable "ASC7" from AREA SOURCES file', 2 )
        END IF

        IF ( .NOT. READ3( ENAME, 'ASC3', ALLAYS3, 0, 0,  ISCC3 ) ) THEN
            CALL M3EXIT( 'GRDAMAT', 0, 0, 
     &      'Error reading variable "ASC3" from AREA SOURCES file', 2 )
        END IF


C...........   Compute gridding matrix
        
        WRITE( *,92000 ) ' ', 'Computing gridding matrix...', ' '

        OFLAG = .FALSE.         !  flag to detect overflow
        EFLAG = .FALSE.         !  flag to detect other errors

        DO  266  I = 1, NASRC
            
            FIP = IFIPS( I )
            ID7 = ISCC7( I )
            ID3 = ISCC3( I )

C...........   Get surrogate ID for this source: 2 cases.
C...........   (Possibly empty) search:  FIP-specific-surrogate case
            
            S = FIND3( FIP, ID7, ID3, NUMSCC, SRFIP, ASCI7, ASCI3 )
                
            IF ( S .GT. 0 ) THEN
                
                SID = SRGID( S )
                
            ELSE        !  FIP-independent case

                S = FIND2( ID7, ID3, NUMSC0, ASCI70, ASCI30 )

                IF ( S .LT. 0 )  THEN !  keep searching

                    WRITE( *,92010 ) 
     &                  'FIP code ', FIP,
     &                  'ID7 code ', ID7,
     &                  'ID3 code ', ID3
                    CALL M3EXIT( 'GRDAMAT', 0, 0, 
     &                           'Surrogate not found.', 2 )

                END IF  !  if <ID7,ID3> not found

                SID = SRGID0( S )
                
            END IF


C...........   Find all cells intersecting this FIPS and fill in 
C...........   corresponding gridding-matrix coefficients.
C...........   ASSUMPTION:  list is sorted by FIP (and then by col-row),
C...........   with multiple entries per FIP.
            
            F = FIND1( FIP, NFIPLIST, FIPLIST )
            
            IF ( F .LT. 0 ) THEN
                WRITE( MESG,94010 )
     &          'FIP code ', FIP, 'not found in surrogate list'
                CALL M3MESG( MESG )
                GO TO  266
            END IF

            DO  255  K = 1, NCELLS( F )
            
                C = FIPCELL( K,F )
                J = NX( C ) + 1

                IF ( J .GT. NACEL ) THEN
                    OFLAG = .TRUE.      !  overflow occurred.
                    WRITE( *,92010 )  
     &              'FIP  code ', FIP,
     &              'ASC7 code ', ID7,
     &              'ASC3 code ', ID3,
     &              'cell code ', C, ' '
                ELSE
                    IS( J,C ) = I
                    CS( J,C ) = SRGFRAC( SID,K,F )
                END IF
                NX( C )   = J

255         CONTINUE    !  end of loop on cells K for this FIP

266     CONTINUE        !  end loop on sources I, computing gridding matrix.
        
        CMAX = NX( 1 )
        CMIN = CMAX
        
        IF ( OFLAG ) THEN   !  overflow occurred:  just compute statistics.
            
            K = 0
                
            DO  277  R = 2, NGRID
            
                J = NX( R )
                K = K + J
                   
                IF      ( J .GT. CMAX ) THEN
                    CMAX = J
                ELSE IF ( J .LT. CMIN ) THEN
                    CMIN = J
                END IF

277         CONTINUE    !  end of loop on cells K for this FIP

            WRITE( *,92000 )
     &          ' ',
     &          'WARNING:  ' // GNAME( 1:TRIMLEN( GNAME ) ) // 
     &          ' not written.',
     &          'Internal arrays would have overflowed.'
            WRITE( *,92010 ) 
     &          'Sources per cell dimensioned', NACEL

        ELSE IF ( EFLAG ) THEN          ! write gridding matrix to file

            CALL M3EXIT( 'GRDAMAT', 0, 0, 'Missing FIP codes', 3 ) 

        ELSE IF ( GFLAG ) THEN          ! write gridding matrix to file

C.......   Compress matrix into I/O representation from scratch representation.
C.......   Compute statistics.
            
            K = 0
                
            DO 288  C = 1, NX( 1 )
                K       = K + 1
                IX( K ) = IS( C,1 )
                CX( K ) = CS( C,1 )
288         CONTINUE
                   
            DO  299  R = 2, NGRID
            
                J = NX( R )
                   
                IF (      J .GT. CMAX ) THEN
                    CMAX = J
                ELSE IF ( J .LT. CMIN ) THEN
                    CMIN = J
                END IF
                   
                DO 298  C = 1, J
                    K       = K + 1
                    IF ( K .LE. NMATX ) THEN
                        IX( K ) = IS( C,R )
                        CX( K ) = CS( C,R )
                    END IF
298             CONTINUE
                   
299         CONTINUE    !  end of loop on cells K for this FIP
            
            IF ( K .GT. NMATX ) THEN
                WRITE( MESG,94010 ) 
     &                  'Number of coeffs', K, 
     &                  'exceeds NMATX=', NMATX
                CALL M3EXIT( 'GRDAMAT', 0, 0, MESG, 3 ) 
            END IF
                   
C.......   Write out the gridding matrix:
            
            WRITE( *,92000 ) ' ', 'Writing out GRIDDING MATRIX file...'

            IF ( .NOT. WRITE3( GNAME, 'ALL', 0, 0, NX ) ) THEN
                CALL M3EXIT( 'GRDAMAT', 0, 0, 
     &              'Error writing GRIDDING MATRIX file.', 2 )
            END IF

        ELSE    !  not GFLAG:  report matrix stats only
        
C.......   Just compute statistics.
            
            K = 0
                
            DO  303  R = 2, NGRID
            
                J = NX( R )
                K = K + J
                   
                IF      ( J .GT. CMAX ) THEN
                    CMAX = J
                ELSE IF ( J .LT. CMIN ) THEN
                    CMIN = J
                END IF

303         CONTINUE    !  end of loop on cells K for this FIP

        END IF  !  if gflag, or not
                        
C.......   Report statistics:

        XDUM = FLOAT( K ) / FLOAT( NGRID )
        WRITE( *,92010 ) 
     &      'Total number of coefficients   ', K   ,
     &      'Max  number of sources per cell', CMAX,
     &      'Min  number of sources per cell', CMIN
        WRITE( *,92020 ) 
     &      'Mean number of sources per cell', XDUM,
     &      ' '
           

C...............   End of program

999     CONTINUE
      
        CALL M3EXIT( 'GRDAMAT', 0, 0, 
     &               'SUCCESSFUL COMPLETION of program GRDAMAT', 0 )


C******************  FORMAT  STATEMENTS   ******************************

C...........   Informational (LOG) message formats... 92xxx

92000   FORMAT( 5X, A )

92010   FORMAT( 5X, A, :, I12 )

92020   FORMAT( 5X, A, :, F17.4 )


C...........   Formatted file I/O formats............ 93xxx

93020   FORMAT( I7, I3 )

93030   FORMAT( I5, I4, 1X, A7, A3 )        !  for surrogate XREF file


C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I10, :, 2X ) )


        END


