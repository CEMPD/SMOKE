
        PROGRAM SRGGRID

C***********************************************************************
C  program body starts at line 178
C
C  DESCRIPTION:
C       Construct by-FIP and sum-total gridded files containing data from
C       EPS-style surrogates file.
C
C  PRECONDITIONS REQUIRED:
C       Sorted, cut-down input data for surrogate coeffs.
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C       Models-3 I/O
C       FIND1, GETYN, PROMPTFFILE, PROMPTMFILE, TRIMLEN
C
C  REVISION  HISTORY:
C       Prototype  9/96 by CJC.
C       NOTE: County totals processing not QA'd
C
C***********************************************************************
C
C Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE)
C                Modeling System
C File Version  @(#)$Id$
C
C COPYRIGHT (C) 1996, MCNC--North Carolina Supercomputing Center
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
C Pathname:     $Source$
C Last updated: $Date$ 
C
C***************************************************************************/

      IMPLICIT NONE

C...........   INCLUDES:

        INCLUDE 'ARDIMS3.EXT'   !  area-source dimensioning parameters
        INCLUDE 'GRDIMS3.EXT'   !  grid parameters
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
        INCLUDE 'FDESC3.EXT'    !  I/O API file description data structures.


C...........   EXTERNAL FUNCTIONS and their descriptions:

        CHARACTER*2     CRLF
        LOGICAL         DSCGRID
        INTEGER         FIND1
        INTEGER         GETMENU
        INTEGER         GETNUM
        LOGICAL         GETYN
        INTEGER         LBLANK
        INTEGER         PROMPTFFILE
        CHARACTER*16    PROMPTMFILE
        INTEGER         STR2INT
        REAL            STR2REAL
        INTEGER         TRIMLEN
        
        EXTERNAL  CRLF, DSCGRID, FIND1, GETMENU, GETNUM, GETYN, LBLANK,
     &            PROMPTFFILE, PROMPTMFILE, STR2INT, STR2REAL, TRIMLEN


C...........   PARAMETERS:

        INTEGER          NRMAX                  !  max # of surrogate records
        CHARACTER*5      BLANK5
        PARAMETER     (  NRMAX = 275 000
     &                 , BLANK5= ' '     )


C...........   LOCAL VARIABLES and their descriptions:
C...........   NOTE that ASC (area-source-category) ID's are 10-digit
C...........   FIPS/areas tables

        INTEGER     NSTALIST              !  # of FIPS
        INTEGER     STALIST( NAFIP )      !  STATE-codes list
        INTEGER     STAINDX( NAFIP )      !  1st loc of state in FIPLIST

        INTEGER     NFIPLIST                     !  # of FIPS
        INTEGER     FIPLIST(        NAFIP )      !  FIP-codes list
        REAL        FIPINDX(        NAFIP )      !  1st loc of FIP in FIPS()
        REAL        FIPSIZE( NASRG, NAFIP )      !  1 / sum_cells( srg( fip ) )

C...........   Surrogate-cell::FIPS tables (unsorted; sorted and processed)

        INTEGER     NSREC                 !  number of (unsorted) entries
        INTEGER     INDXA( NRMAX )        !  subscripts for sorti1()
        INTEGER     CELLA( NRMAX )        !  cell number,        before sorting
        INTEGER     FIPSA( NRMAX )        !  FIP code,           before sorting
        INTEGER     SSCSA( NRMAX )        !  surrogate ID,       before sorting
        REAL        FRACA( NRMAX )        !  surrogate fraction, before sorting

        INTEGER     CELL( NRMAX )        !  cell number,        before sorting
        INTEGER     FIPS( NRMAX )        !  FIP code,           before sorting
        INTEGER     SSCS( NRMAX )        !  surrogate ID,       before sorting
        REAL        FRAC( NRMAX )        !  surrogate fraction, before sorting

C...........   Output grids:

        REAL            GRID( NGRID, NASRG )

C...........   Logical names and unit numbers

        INTEGER         LDEV    !  log-device
        INTEGER         SDEV    !  for surrogate coeff file
        INTEGER         TDEV    !  for surrogate total file

        CHARACTER*16    GNAME   !  logical name for grid matrix     output file

C...........   Other local variables

        REAL            DDX, DDY                !  inverse of cell widths
        REAL            SRG
        REAL            SRGT( NASRG )
        REAL            X1, X2, Y1, Y2, XC, YC  !  surrogate file header
        REAL            XOFF, YOFF              !  surrogate offsets from grid
        REAL            XX, YY                  !  cell number or l-l position

        INTEGER         C, S, F, J, K, N, V  !  counters, subscripts
        INTEGER         CHOICE
        INTEGER         COL,  ROW
        INTEGER         FIP, STA, SSC, SID, CID
        INTEGER         HCOL, HROW              !  surrogate hdr cols and rows
        INTEGER         JDATE, JTIME
        INTEGER         IOS     !  I/O status
        INTEGER         IREC    !  input line (record) number
        INTEGER         L1, L2
        INTEGER         LFIP, LSTA
        INTEGER         LO, HI
        INTEGER         LOLIST, HILIST
        INTEGER         NSRG    !  Count of srgts in EPS2 srgts
        INTEGER         NVAROUT !  Number of output surrogates
        INTEGER         SSCT ( NASRG )
        INTEGER         VINDX( NASRG )

        LOGICAL         EFLAG   !  input error flat
        LOGICAL         TFLAG   !  iff totals file provided

        CHARACTER*5     FFORMAT !  temporary indicator for input formats
        CHARACTER*16    SCRBUF
        CHARACTER*16    NAME, VNAME( NASRG )
        CHARACTER*80    DESC, VDESC( NASRG )
        CHARACTER*256   LINE    !  input buffer 
        CHARACTER*256   MESG

        CHARACTER*72    MENU( 3 )
     &  /
     &  'Compute gridded surrogates for ALL COUNTIES',	!  choice 1
     &  '... for all counties in a selected STATE',	!  choice 2
     &  '... for a selected COUNTY in a selected STATE'	!  choice 3
     &  /

C...........   STATEMENT FUNCTIONS:
C.......   floating point "unequal" -- true iff
C.......   | P - Q | > 1e-5 * sqrt( p*p + q*q + 1e-5 )
 
        REAL            P, Q
        LOGICAL         FLTERR
        FLTERR( P, Q ) =
     &  ( (P - Q)**2  .GT.  1.0E-10*( P*P + Q*Q + 1.0E-5 ) )


C***********************************************************************
C   begin body of program SRGGRID

        LDEV = INIT3()

        WRITE( *,92000 ) 
     &  ' ',
     &  'Program SRGGRID to take the EMS-95 SURROGATE COEFFICIENT,',
     &  'file and the SURROGATE TOTAL BY FIP file, and produce gridded',
     &  'SURROGATE TOTAL and SURROGATE BY-FIPS file.', 
     &  ' ',
     &  'NOTE: If the EPS2.0 surrogate coefficient file is used, it',
     &  '      must use the ZONE field of the header line (line 1,',
     &  '      columns 41-50) to indicate the number of surrogates',
     &  '      in the file, and each surrogates column must be',
     &  '      compatible with an f10.8 format.  In addition, SMOKE',
     &  '      permits the x-cell and y-cell numbers to be used in',
     &  '      place of the (x,y) coordinates.',
     &  ' ',
     &  'NOTE: The EMS95-formatted surrogates file refers to a',
     &  '      format created by writing the SAS srgratio.ssd EMS95',
     &  '      file to an ASCII file.  The file format is delimited',
     &  '      columns of surrogate ID code, state code, county code,',
     &  '      grid column number, grid row number, and surrogate',
     &  '      fraction, in 59 columns or less.',
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
            CALL M3EXIT( 'SRGGRID', 0, 0, 'Ending program SRGGRID', 2 )
        END IF


C.......   Get file name; open actual-scc file

        SDEV = PROMPTFFILE( 
     &  'Enter logical name for SURROGATE COEFF file',
     &           .TRUE., .TRUE., 'AGPRO', 'SRGGRID' )

        TDEV = PROMPTFFILE( 
     &  'Enter logical name for SURROGATE TOTAL file or "NONE"',
     &           .TRUE., .TRUE., 'AGTOT', 'SRGGRID' )

        TFLAG = ( TDEV .NE. -2 )

C.......   Get the grid name and parameters
C.......   to put into file header:

        IF ( .NOT. DSCGRID( GRDNM, SCRBUF, GDTYP3D, 
     &              P_ALP3D, P_BET3D, P_GAM3D, XCENT3D, YCENT3D,
     &              XORIG3D, YORIG3D, XCELL3D, YCELL3D, 
     &              NCOLS3D, NROWS3D, NTHIK3D ) ) THEN

            SCRBUF = GRDNM
            MESG   = '"' // SCRBUF( 1:TRIMLEN( SCRBUF ) ) //
     &               '" not found in GRIDDESC file'
            CALL M3EXIT( 'SRGGRID', 0, 0, MESG, 2 )

        END IF                  ! if dscgrid() failed


C........  Determine input format of gridding surrogates file
C........  For EPS2 format, treat 1st line as header line, and QA file
C........  For Lat-lon grids, header will be irrelevant, but can't check either
 
        FFORMAT = 'EPS2'
        READ( SDEV, 93000, END=27, IOSTAT=IOS ) LINE
 
        IF( TRIMLEN( LINE ) .LT. 60 ) THEN  ! EMS-95 format
 
            FFORMAT = 'EMS95'
            NSRG    = 1
            REWIND( SDEV )
 
        ELSE                                ! EPS2   format
 
            IREC = IREC + 1
 
            IF( GDTYP3D .NE. LATGRD3 ) THEN
 
                X1  = STR2REAL( LINE(  1:10 ) )
                Y1  = STR2REAL( LINE( 11:20 ) )
                X2  = STR2REAL( LINE( 21:30 ) )
                Y2  = STR2REAL( LINE( 31:40 ) )
                NSRG= STR2INT ( LINE( 41:50 ) )
                XC  = STR2REAL( LINE( 51:60 ) )
                YC  = STR2REAL( LINE( 61:70 ) )
 
                IF( XC .EQ. 0 ) THEN
                    CALL M3EXIT( 'GRDAMAT', 0, 0,
     &                   'XCELL is zero in SURROGATE file header', 2 )
 
                ELSEIF( YC .EQ. 0 ) THEN
                    CALL M3EXIT( 'GRDAMAT', 0, 0,
     &                   'YCELL is zero in SURROGATE file header', 2 )
 
                ELSE
                    DDX = 1.0 / XC
                    DDY = 1.0 / YC
                ENDIF
 
                XOFF= DDX * ( X1 - XORIG3D )
                YOFF= DDY * ( Y1 - YORIG3D )
                HCOL= NINT( DDX * ( X2 - X1 ) )
                HCOL= NINT( DDX * ( X2 - X1 ) )
                HROW= NINT( DDY * ( Y2 - Y1 ) )
 
                IF ( FLTERR( XC, SNGL( XCELL3D ) ) ) THEN
                    CALL M3EXIT( 'GRDAMAT', 0, 0,
     &                           'Bad XCELL in SURROGATE file', 2 )
                ELSEIF ( FLTERR( YC, SNGL( YCELL3D ) ) ) THEN
                    CALL M3EXIT( 'GRDAMAT', 0, 0,
     &                           'Bad YCELL in SURROGATE file', 2 )
                ELSEIF ( FLTERR( XOFF, FLOAT( NINT( XOFF ) ) ) ) THEN
                    CALL M3EXIT( 'GRDAMAT', 0, 0,
     &                     'Bad X alignment in SURROGATE file', 2 )
                ELSEIF ( FLTERR( YOFF, FLOAT( NINT( YOFF ) ) ) ) THEN
                    CALL M3EXIT( 'GRDAMAT', 0, 0,
     &                     'Bad Y alignment in SURROGATE file', 2 )
                ELSEIF ( FLTERR( REAL( HCOL ), REAL( NCOLS ) ) ) THEN
                    CALL M3EXIT( 'GRDAMAT', 0, 0,
     &                     'Bad ending X-coord in SURROGATE file', 2 )
                ELSEIF ( FLTERR( REAL( HROW ), REAL( NROWS ) ) ) THEN
                    CALL M3EXIT( 'GRDAMAT', 0, 0,
     &                     'Bad ending Y-coord in SURROGATE file', 2 )
                ELSEIF ( NSRG .LE. 0 ) THEN
                    WRITE( MESG, 94010 )
     &               'NOTE: Use ZONE portion of EPS2.0 surrogates '//
     &               'header' // CRLF() // BLANK5 //
     &               '      to specify number of surrogates ' //
     &               'in the file.'
                    CALL M3MSG2( MESG )
 
                    WRITE( MESG, 94010 )
     &               'Number of surrogates in header was ', NSRG
                    CALL M3EXIT( 'GRDAMAT', 0, 0, MESG, 2 )
 
                ELSEIF ( NSRG .GT. NASRG ) THEN
                    WRITE( MESG, 94010 )
     &               'Number of surrogates in file =', NSRG,
     &               'but maximum (NASRG) =', NASRG
                    CALL M3EXIT( 'GRDAMAT', 0, 0, MESG, 2 )
 
                END IF
 
            ELSE  ! Grid type is LAT-LON
 
                MESG= 'Using LAT-LON grid, with EPS2.0 surrogates file.'
                CALL M3WARN( 'GRDAMAT', 0, 0, MESG )
 
            END IF
 
        ENDIF

        WRITE( *,92000 ) 
     &      ' ', 'Reading GRIDDING SURROGATE COEFF file...', ' '

C...........   Start loop for actually reading surrogates data
 
22      CONTINUE                !  head of the SDEV-read loop
 
C...............  Read EMS-95 formatted file
            IF( FFORMAT .EQ. 'EMS95' ) THEN
 
                READ( SDEV, *, END=27, IOSTAT=IOS )
     &              SSCT( NSRG ),   !  surrogate ID
     &              SID, CID,       !  state, county FIP codes
     &              COL, ROW,       !  grid col and row
     &              SRGT( NSRG )    !  surrogate fraction
 
                IREC = IREC + 1
 
                FIP  = 1000 * SID + CID
 
C...............  Read EPS2 file.  This read is made somewhat more complex
C...............  because we're conforming to EMS95 file structure.  However,
C...............  since the EPS2.0 format does not specify exact column
C...............  numbers for each surrogate class (just must start on
C...............  column 36 and continue to column 186), and the number of
C...............  surrogates is not pre-defined, the EPS2.0 format is ill-
C...............  defined anyway.  Here, we add the constraint that the
C...............  zone field is used to define the number of surrogate
C...............  entries. Also, all file columns are as described in the
C...............  EPS2 documentation, AND the surrogates entries are 10
C...............  columns wide.  These constraints allow us to read in the
C...............  file more assuredly.  MRH.
            ELSEIF( FFORMAT .EQ. 'EPS2' ) THEN
 
                READ( SDEV, 93010, END=27, IOSTAT=IOS )
     &                FIP, XX, YY, ( SRGT( J ), J = 1, NSRG )
 
                IREC = IREC + 1
 
C.................  First try col/row format
                COL = INT( XX )
                ROW = INT( YY )
 
C.................  If no good, then try coordinates format
                IF( COL .LT. 1  .OR. COL .GT. NCOLS .OR.
     &              ROW .LT. 1  .OR. ROW .GT. NROWS      ) THEN
 
                    COL = 1 + NINT( DDX * ( XX - XORIG3D ) )
                    ROW = 1 + NINT( DDY * ( YY - YORIG3D ) )
 
                ENDIF
 
                DO 23 J = 1, NSRG
                    SSCT( J ) = J
23              CONTINUE
 
            END IF         ! end input format type
 
            IF ( IOS .NE. 0 ) THEN
 
                WRITE( MESG,94010 )
     &              'I/O error ', IOS,
     &              'reading SURG COEFF file at line', IREC
                CALL M3MESG( MESG )
                EFLAG = .TRUE.
                GO TO  22   ! to head of SDEV read loop
 
C.................  If columns and rows are still a problem, give message,
C.................  but don't error because may want to use bigger
C.................  surrogates file to do smaller domain
            ELSE IF ( COL .LT. 1  .OR. COL .GT. NCOLS ) THEN
 
                WRITE( MESG,94010 )
     &              'Column ', COL,
     &              'out of range in SURG COEFF file at line', IREC
                CALL M3MESG( MESG )
                GO TO  22   ! to head of SDEV read loop
 
            ELSE IF ( ROW .LT. 1  .OR. ROW .GT. NROWS ) THEN
 
                WRITE( MESG,94010 )
     &              'Row ', ROW,
     &              'out of range in SURG COEFF file at line', IREC
                CALL M3MESG( MESG )
                GO TO  22   ! to head of SDEV read loop
 
            END IF              !  if IOS bad, or col, or row out of range

C...........   Record surrogate fraction for this cell.
C...........   Loop through the number of surrogates in current record
            DO 25 J = 1, NSRG
 
                N = N + 1
 
                IF( SSCT( J ) .LT. 1     .OR.
     &              SSCT( J ) .GT. NASRG      ) THEN
 
                    WRITE( MESG,94010 )
     &                  'Surrogate code ', SSCT( J ),
     &                  'out of range (NASRG) in ASPRO at line', IREC
                    CALL M3MESG( MESG )
                    EFLAG = .TRUE.
                    GO TO  22   ! to head of SDEV read loop
 
                ELSEIF( N .LE. NRMAX ) THEN
 
                    INDXA( N ) = N
                    CELLA( N ) = COL  +  NCOLS * ( ROW - 1 )
                    SSCSA( N ) = SSCT( J )
                    FIPSA( N ) = FIP
                    FRACA( N ) = SRGT( J )
 
                END IF                      ! if n in bounds
 
25          CONTINUE
 
            GO TO  22  !  to head of SDEV read loop
 
27      CONTINUE       !  end of the SDEV-read loop

        CALL M3MSG2( 'NOTE: File read in as ' // FFORMAT // ' format.' )



C...........   Read and process the surrogate coefficient file:
C...............   Must be sorted by FIP (state ID, county ID), surrogate ID,
C...............   col, and row.  
C...............   All actually-occurring area source FIP codes must be present.

c        IREC    =  1     
c        N       =  0
c        EFLAG   = .FALSE.
c
c22      CONTINUE                !  head of the SDEV-read loop
c
c            IREC = IREC + 1
c            READ( SDEV, *, END=23, IOSTAT=IOS ) 
c     &          SSC,            !  surrogate ID
c     &          SID, CID,       !  state, county FIP codes
c     &          COL, ROW,       !  grid col and row
c     &          SRG             !  surrogate fraction
c
c            IF ( IOS .NE. 0 ) THEN
c
c                WRITE( MESG,94010 ) 
c     &              'I/O error ', IOS, 
c     &              'reading SURG COEFF file at line', IREC
c                CALL M3MESG( MESG )
c                EFLAG = .TRUE.
c                GO TO  22
c
c            ELSE IF ( COL .LT. 1  .OR. COL .GT. NCOLS ) THEN
c
c                WRITE( MESG,94010 )
c     &              'Column ', COL,
c     &              'out of range in SURG COEFF file at line', IREC
c                CALL M3MESG( MESG )
c                EFLAG = .TRUE.
c                GO TO  22
c
c            ELSE IF ( ROW .LT. 1  .OR. ROW .GT. NROWS ) THEN
c
c                WRITE( MESG,94010 )
c     &              'Row ', ROW,
c     &              'out of range in SURG COEFF file at line', IREC
c                CALL M3MESG( MESG )
c                EFLAG = .TRUE.
c                GO TO  22
c
c            ELSE IF ( SSC .LT. 1  .OR. SSC .GT. NASRG ) THEN
c
c                WRITE( MESG,94010 )
c     &              'Surrogate code ', SSC,
c     &              'out of range in SURG COEFF file at line', IREC
c                CALL M3MESG( MESG )
c                EFLAG = .TRUE.
c                GO TO  22
c
c            END IF              !  if IOS bad, or col,row, or ssc out of range
c
c
cC...........   Record surrogate fraction for this cell.
c
c            N = N + 1
c
c            IF ( N .LE. NRMAX ) THEN
c
c                INDXA( N ) = N
c                CELLA( N ) = COL  +  NCOLS * ( ROW - 1 )
c                SSCSA( N ) = SSC
c                FIPSA( N ) = 1000 * SID + CID
c                FRACA( N ) = SRG
c
c            END IF                      ! if n in bounds
c
c            GO TO  22                   !  to head of loop
c
c23      CONTINUE                        !  end of the SDEV-read loop


        NSREC = N
        WRITE( MESG,94010 ) 
     &      'Max SURG COEF records - actual:', N,
     &      'dimensioned:', NRMAX
        CALL M3MESG( MESG )

        IF ( N .GT. NRMAX ) THEN
            CALL M3EXIT( 'SRGGRID', 0, 0, 
     &                   'Table overflow reading SURG COEFF file', 2 )
        END IF

        IF ( EFLAG ) THEN
            MESG = 'Error reading SURG COEF file.'
            CALL M3EXIT( 'SRGGRID', 0, 0, MESG, 2 )
        END IF


C...........   Sort and process the surrogates:

        WRITE( *,92000 ) 
     &      ' ', 'Processing GRIDDING SURROGATE COEFF file...', ' '

        CALL SORTI3( NSREC, INDXA, FIPSA, CELLA, SSCSA )

        F    = 0
        S    = 0
        LFIP = 0
        LSTA = 0

        DO  33  J = 1, NSREC            !  shuffle SURROGATE COEFFs

            K   = INDXA( J )
            C   = CELLA( K )
            FIP = FIPSA( K )
            SID = SSCSA( K )
            SRG = FRACA( K )

            CELL( J ) = C
            FIPS( J ) = FIP
            SSCS( J ) = SID
            FRAC( J ) = SRG

            IF ( FIP .NE. LFIP ) THEN

                F = F + 1

                IF ( F .LT. NAFIP ) THEN
                    FIPLIST( F ) = FIP
                    FIPINDX( F ) = J
                    STA  = FIP / 1000
                    IF ( STA .NE. LSTA ) THEN
                        S = S + 1
                        STALIST( S ) = STA
                        STAINDX( S ) = F
                        LSTA = STA
                    END IF
                END IF          !  if fip found, or not

                LFIP = FIP
                
            END IF              !  if new FIP encountered

33      CONTINUE                !  end loop shuffling SURROGATE COEFFs

        FIPINDX( F+1 ) = NSREC + 1
        STAINDX( S+1 ) = F     + 1
        NFIPLIST = F
        NSTALIST = S

        WRITE( MESG,94010 ) 
     &      'SURG COEF table - actual:', NFIPLIST,
     &      'dimensioned:', NAFIP
        CALL M3MESG( MESG )

        IF ( NFIPLIST .GT. NAFIP ) THEN
            CALL M3EXIT( 'SRGGRID', 0, 0, 
     &      'Table overflow processing SURG COEFF file', 2 )
        END IF


C...........  If surrogates file present, process SURROGATE-TOTAL-BY-FIP file:

        IF( TFLAG ) THEN

            DO  44  F = 1, NFIPLIST
            DO  43  S = 1, NASRG        
                FIPSIZE( S,F ) = 0.0
43          CONTINUE                !  end loop on SID's
44          CONTINUE                !  end loop on FIPS

            IREC    = 0
            NVAROUT = 0
            SSC     = IMISS3
            EFLAG   = .FALSE.

55          CONTINUE                !  head of loop reading TDEV

                READ( TDEV, 93000, END=66, IOSTAT=IOS ) LINE 

                IREC = IREC + 1

                L1 = LBLANK ( LINE )
                L2 = TRIMLEN( LINE )

                IF ( IOS .NE. 0 ) THEN

                    WRITE( MESG,94010 ) 
     &                  'I/O error ', IOS, 
     &                  'reading SURG TOTAL file at line', IREC
                    CALL M3MESG( MESG )
                    EFLAG = .TRUE.
                    GO TO  55

                ELSEIF( LINE( L1+1:L1+7 ) .EQ. '/STYPE/' ) THEN

                    READ( LINE( L1+8:L2 ), * ) SSC, NAME, DESC

                    IF ( SSC .LT. 1  .OR. SSC .GT. NASRG ) THEN

                        WRITE( MESG,94010 )
     &                    'Surrogate code ', SSC, 'out of range' //
     &                    'in SURG TOTAL file at line', IREC
                        CALL M3MESG( MESG )
                        EFLAG = .TRUE.

                    ELSE

                        NVAROUT = NVAROUT + 1
                        VNAME( NVAROUT ) = NAME 
                        VDESC( NVAROUT ) = DESC
                        VINDX( SSC ) = NVAROUT

                    END IF  !  if IO error, or if SID out of range

                    GO TO  55

                ELSEIF( SSC .EQ. IMISS3 ) THEN

                    MESG = 'Total surrogates file must use ' //
     &                     '/STYPE/ packet in the form: ' // CRLF() //
     &                     BLANK5 // '/STYPE/ <SSC code> ' //
     &                     '"<OUTPUT VAR NAME>" "<OUTPUT VAR DESC>"'
                    CALL M3EXIT( 'SRGGRID', 0, 0, MESG, 2 )

                ENDIF

                FIP = STR2INT( LINE( 1:2 ) ) * 1000 +
     &                STR2INT( LINE( 3:5 ) )

                SRG = STR2REAL( LINE( 6:L2 ) )


                F = FIND1( FIP, NFIPLIST, FIPLIST )

                IF ( F .GE. 0 ) THEN
                    FIPSIZE( SSC,F ) = SRG
                END IF

                GO TO  55           !  to head of loop

66          CONTINUE                !  end of loop reading TDEV

            IF ( EFLAG ) THEN
                MESG = 'Error reading SURG TOTAL file.'
                CALL M3EXIT( 'SRGGRID', 0, 0, MESG, 2 )
            END IF

C...............   Process gridded surrogate-totals:

            DO  122  S = 1, NASRG
            DO  111  C = 1, NGRID
                GRID( C,S ) = 0.0
111         CONTINUE
122         CONTINUE

            LFIP = -1
            DO  133  K = 1, NSREC
                FIP = FIPS( K )
                C   = CELL( K )
                S   = SSCS( K )
                V   = VINDX( S )

                IF ( FIP .NE. LFIP ) THEN
                    F   = FIND1( FIP, NFIPLIST, FIPLIST )
                    LFIP = FIP
                END IF

                IF ( V .GT. 0 ) THEN

                    GRID( C,V )= GRID( C,V ) + FRAC( K )*FIPSIZE( S,F )

                END IF

133         CONTINUE


        ENDIF   ! End if surrogates file present or not

C...........   Select sublist to compute:

        CHOICE = GETMENU( 3, 2, 'Enter type of county-selection',
     &                    MENU )

        IF ( CHOICE .EQ. 1 ) THEN

            LOLIST = 1
            HILIST = NFIPLIST

        ELSE IF ( CHOICE .EQ. 2 ) THEN			!  select a state

            WRITE( *,92010 )
     &          'The available state codes are:',
     &          ( S, STALIST( S ), S = 1, NSTALIST )
            S = GETNUM( 1, NSTALIST, 1, 
     &                  'Enter number for state selected' )
            LOLIST = STAINDX( S )
            HILIST = STAINDX( S+1 ) - 1

        ELSE IF ( CHOICE .EQ. 3 ) THEN			!  select a fip

            WRITE( *,92010 )
     &          'The available state codes are:',
     &          ( S, STALIST( S ), S = 1, NSTALIST )
            S = GETNUM( 1, NSTALIST, 1, 
     &                  'Enter number for state selected' )
            LO = STAINDX( S )
            HI = STAINDX( S+1 ) - 1
            WRITE( *,92010 )
     &          'The available county codes for this state are:',
     &          ( F, FIPLIST( F ), F = LO, HI )
            F = STAINDX( S ) + GETNUM( 1, NSTALIST, 1, 
     &                                 'Enter number for state' )
            LOLIST = F
            HILIST = F

        END IF			!  if all counties, else if state, else if fip


        FTYPE3D = GRDDED3
        SDATE3D = 1
        STIME3D = 0
        TSTEP3D = 10000
        NVARS3D = NVAROUT
        NLAYS3D = 1
        GDNAM3D = GRDNM
        VGTYP3D = IMISS3
        VGTOP3D = BADVAL3
        FDESC3D( 1 ) = 'Gridding surrogates, by FIP'
        FDESC3D( 2 ) = 
     &  'First record is totals, SUM_fip( surg( fip )*tot( fip ) )'
        FDESC3D( 3 ) = 
     &  'Other records are FIP-specific gridded surrogate fractions'
        FDESC3D( 4 ) = 
     &  'with the FIPS in consecutive sorted order.'

        DO  166  K = 5, MXDESC3
            FDESC3D( K ) = ' '
166     CONTINUE

        DO  188  K = 1, NVAROUT

            VNAME3D( K ) = VNAME( K )
            UNITS3D( K ) = 'n/a'
            VDESC3D( K ) = VDESC( K )
            VTYPE3D( K ) = M3REAL

188     CONTINUE

        GNAME = PROMPTMFILE( 
     &      'Enter logical name for GRIDDING SURROGATES output file',
     &      FSUNKN3, 'ASURG', 'SRGGRID' )


C...........   Process gridded county-level surrogates

        WRITE( *,92000 ) 
     &  ' ', 
     &  'Constructing GRIDDED SURROGATE FRACTIONS output file...',
     &  ' '

        JDATE = 1
        JTIME = 0

        IF ( .NOT. WRITE3( GNAME, 'ALL', JDATE, JTIME, GRID ) ) THEN

            MESG = 'Could not write SURROGATE TOTALS to "' //
     &             GNAME( 1 : TRIMLEN( GNAME ) ) // '"'
            CALL M3EXIT( 'SRGGRID', 0, 0, MESG, 2 )

        END IF

        CALL NEXTIME( JDATE, JTIME, 10000 )

        DO  299  F = LOLIST, HILIST                !  loop on counties

            DO  222  S = 1, NASRG
            DO  211  C = 1, NGRID
                GRID( C,S ) = 0.0
211         CONTINUE
222         CONTINUE

            LO = FIPINDX( F   )
            HI = FIPINDX( F+1 ) - 1

            DO  288  K = LO, HI         !  loop on records K for this fip

                C = CELL( K )
                S = VINDX( SSCS( K ) )

                IF( S .GT. 0 ) THEN
                    GRID( C,S ) = GRID( C,S )  +  FRAC( K )
                ENDIF

288         CONTINUE                    !  end loop on records K for this fip

            IF ( WRITE3( GNAME, 'ALL', JDATE, JTIME, GRID ) ) THEN

                WRITE( MESG,94020 ) 
     &              'FIP', FIPLIST( F ),  
     &              'written to "' // GNAME( 1 : TRIMLEN( GNAME ) ) //
     &              '" as record', F - LOLIST + 2
                CALL M3MESG( MESG )

            ELSE                        !  else write3() failed

                WRITE( MESG,94020 ) 
     &              'Could not write', FIP,  
     &              'to "' // GNAME( 1 : TRIMLEN( GNAME ) ) // '"'
                CALL M3EXIT( 'SRGGRID', 0, 0, MESG, 2 )

            END IF                      !  if write3() succeeded, or not
            
            CALL NEXTIME( JDATE, JTIME, 10000 )

299     CONTINUE                        !  end loop on counties F


C...............   End of program

        CALL M3EXIT( 'SRGGRID', 0, 0, 
     &               'SUCCESSFUL COMPLETION of program SRGGRID', 0 )


C******************  FORMAT  STATEMENTS   ******************************

C...........   Informational (LOG) message formats... 92xxx

92000   FORMAT( 5X, A )

92010   FORMAT( /5X, A , /, 10 ( /, 5( 4X, I3, ':', I5, : ) ) )


C...........   Formatted file I/O formats............ 93xxx

93000   FORMAT( A )

93010   FORMAT( I5, 2F10.0, 10X, <NSRG>( F10.8 ) )

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I10, :, 2X ) )

94020   FORMAT( A, :, I7, :, 2X, A, :, I10 )


        END

