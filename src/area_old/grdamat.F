
        PROGRAM GRDAMAT

C***********************************************************************
C  program body starts at line 209
C
C  DESCRIPTION:
C       Construct area source gridding matrix from data contained in
C       EPS-style surrogates file.
C
C  PRECONDITIONS REQUIRED:
C       Sorted, cut-down input data for surrogate coeffs.
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C       RDASCC()
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
C
C       Version   9/96 by CJC reads EMS-95 formatted files
C***********************************************************************
C
C Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
C                System
C File: @(#)$Id$
C
C COPYRIGHT (C) 1998, MCNC--North Carolina Supercomputing Center
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
C***************************************************************************/

      IMPLICIT NONE

C...........   INCLUDES:

        INCLUDE 'ARDIMS3.EXT'   !  area-source dimensioning parameters
        INCLUDE 'CHDIMS3.EXT'   !  emis chem parms (inventory + model)
        INCLUDE 'GRDIMS3.EXT'   !  grid parameters
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
        INCLUDE 'FDESC3.EXT'    !  I/O API file description data structures.


C...........   EXTERNAL FUNCTIONS and their descriptions:

        CHARACTER*2     CRLF
        LOGICAL         DSCGRID
        INTEGER         FIND1, FIND2, FIND3
        LOGICAL         GETYN
        INTEGER         PROMPTFFILE
        CHARACTER*16    PROMPTMFILE
        INTEGER         STR2INT
        REAL            STR2REAL
        INTEGER         TRIMLEN
        
        EXTERNAL        CRLF, DSCGRID, FIND1, FIND2, FIND3, GETYN, 
     &                  PROMPTFFILE, PROMPTMFILE, STR2INT, STR2REAL, 
     &                  TRIMLEN

C...........   PARAMETERS:

        INTEGER         NRMAX                   !  max # of surrogate records
        INTEGER         NSMAX                   !  max # of surrogate coeffs
        INTEGER         NXMAX                   !  max # of XREF entries
        CHARACTER*5     BLANK5

        PARAMETER     ( NRMAX = 275 000 ,
     &                  NSMAX =  45 000 ,
     &                  NXMAX = 240 000 ,
     &                  BLANK5= ' '      )

C...........   LOCAL VARIABLES and their descriptions:
C...........   NOTE that ASC (area-source-category) ID's are 10-digit
C...........   unsigned integers which may be treated as a leading 7-digit
C...........   field, and a trailing 3-digit field.  *7 and *3 arrays below
C...........   follow this scheme with parallel arrays
C...........   Area Sources Table

        INTEGER     IFIPS( NASRC )  	  !  source FIPS (county) ID
        INTEGER     ISCC7( NASRC )  	  !  leading-7  digits of source ASC
        INTEGER     ISCC3( NASRC )  	  !  trailing-3 digits of source ASC
        INTEGER     SRGID( NASRC )  	  !  surrogate number (1...NASRG)

C...........   Actual-FIPS table

        INTEGER     NFIPLIST              !  # of actual FIPS
        INTEGER     FIPLIST( NAFIP     )  !  Actually-occurring FIP-codes list
        INTEGER     SRCLIST( NAFIP + 1 )  !  source index   for FIP, & sentinel

C.......   Actually-occurring ASC table
 
        INTEGER     NASC
        INTEGER     ASCA7( NASCC )
        INTEGER     ASCA3( NASCC )

C...........   Area GRIDDING XREF before sorting

        INTEGER     NAREF		  !  number of XREF entries
        INTEGER     INDEXA( NXMAX )       !  sorting index
        INTEGER     AFIPSA( NXMAX )  	  !  FIPs codes in XREF
        INTEGER     ASCX7A( NXMAX )  	  !  leading-7  digits of ASC in XREF
        INTEGER     ASCX3A( NXMAX )  	  !  trailing-3 digits of ASC in XREF
        INTEGER     SSCXDA( NXMAX )  	  !  surrogate number (1...NASRG)

C...........   Area GRIDDING XREF after sorting and processing
C...........   (TWO tables, one for FIP==0; the other for FIP-specific data.

        INTEGER     NPRNF                 !  no-FIP surrogates-table count
        INTEGER     ASCX70( NASCC )       !  leading-7  digits of ASC
        INTEGER     ASCX30( NASCC )       !  trailing-3 digits of ASC
        INTEGER     SSCXD0( NASCC )       !  surrogate munber (1...NASRG)

        INTEGER     NPRFC                 !  FIP:ASCT surrogates-table count
        INTEGER     AFIPS ( NXMAX )       !  leading-7  digits of ASC
        INTEGER     ASCX7 ( NXMAX )       !  leading-7  digits of ASC
        INTEGER     ASCX3 ( NXMAX )       !  trailing-3 digits of ASC
        INTEGER     SSCXD ( NXMAX )       !  surrogate number (1...NASRG) by ASC

C...........   Surrogate-cell::FIPS tables (unsorted; sorted and processed)

        INTEGER     NSREC		  !  number of (unsorted) entries
        INTEGER     INDXA( NRMAX )	  !  subscripts for sorti1()
        INTEGER     CELLA( NRMAX )        !  cell number,        before sorting
        INTEGER     FIPSA( NRMAX )        !  FIP code,           before sorting
        INTEGER     SSCSA( NRMAX )        !  surrogate ID,       before sorting
        REAL        FRACA( NRMAX )        !  surrogate fraction, before sorting

        INTEGER     NFIP( NGRID )         !  number of fips per cell
        INTEGER     CFIP( NSMAX )         !  cell number::index in FIPLIST
        REAL        SFRC( NASRG, NSMAX )  !  surrogate fractions

C...........   Gridding Matrix

        INTEGER     NX( NGRID )		  !  number of sources per cell
        INTEGER     IX( NMATX )		  !  list of sources per cell
        REAL        CX( NMATX )		  !  coefficients for sources

        COMMON  / GRIDMAT / NX, IX, CX

C...........   Logical file names and unit numbers

        INTEGER         ADEV    !  for actual ASCS file
        INTEGER         LDEV    !  log-device
        INTEGER         SDEV    !  for surrogate coeff file
        INTEGER         XDEV    !  for surrogate xref  file
        CHARACTER*16    ANAME   !  logical name for emission source  input file
        CHARACTER*16    GNAME   !  logical name for grid matrix     output file

C...........   Other local variables

        REAL            DDX, DDY                !  inverse of cell widths
        REAL            SRG
        REAL            SRGT( NASRG )
        REAL            X1, X2, Y1, Y2, XC, YC  !  surrogate file header
        REAL            XOFF, YOFF              !  surrogate offsets from grid
        REAL            XX, YY                  !  cell number or l-l position

        INTEGER         C, S, F, I, J, K, L, M, N, T  !  counters, subscripts
        INTEGER         CMAX, CMIN
        INTEGER         COL,  ROW
        INTEGER         FIP, SSC, SID, CID
        INTEGER         HCOL, HROW              !  surrogate hdr cols and rows
        INTEGER         ID7,  ID3
        INTEGER         IOS                     !  I/O status
        INTEGER         IREC                    !  input line (record) number
        INTEGER         LFIP, LCEL, LAST
        INTEGER         LO, HI
        INTEGER         NSRG                    ! Count of srgts in EPS2 srgts
        INTEGER         SSCT( NASRG )

        LOGICAL         EFLAG   !  input error flat

        CHARACTER*5     FFORMAT !  temporary indicator for input formats
        CHARACTER*16    SCRBUF
        CHARACTER*256   LINE
        CHARACTER*256   MESG

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

        CALL INITEM( LDEV )

        WRITE( *,92000 ) 
     &  ' ',
     &  'Program GRDAMAT to take the EPS2.0 or EMS-95 surrogate',
     &  'coefficient file, the EPS2.0 or EMS95 surrogate cross-',
     &  'reference file, the AREA file produced by EMSAREA, the',
     &  'list of ASCs produced by EMSAREA, and produce the AREA',
     &  'GRIDDING MATRIX file.', 
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
            CALL M3EXIT( 'GRDAMAT', 0, 0, 'Ending program GRDAMAT', 2 )
        END IF


C.......   Get file name; open actual-scc file

        ANAME = PROMPTMFILE( 
     &          'Enter logical name for AREA SOURCE input file',
     &          FSREAD3, 'AREA', 'GRDAMAT' )

        ADEV = PROMPTFFILE(
     &           'Enter logical name for ACTUAL ASCs file',
     &           .TRUE., .TRUE., 'ASCS', 'GRDAMAT' )

        XDEV = PROMPTFFILE( 
     &  'Enter logical name for GRIDDING SURROGATE XREF file',
     &           .TRUE., .TRUE., 'AGREF', 'GRDAMAT' )

        SDEV = PROMPTFFILE( 
     &  'Enter logical name for GRIDDING SURROGATE COEFF file',
     &           .TRUE., .TRUE., 'AGPRO', 'GRDAMAT' )


C.......   Read in the area source emissions FIP and ASC tables:

        CALL M3MSG2( 'Reading in AREA SOURCES file...' )

        IF ( .NOT. DESC3( ANAME ) ) THEN
            CALL M3EXIT( 'GRDAMAT', 0, 0, 
     &      'Error reading HEADER from AREA SOURCES file', 2 )
        ELSE IF ( NROWS3D .NE. NASRC ) THEN
            WRITE( MESG, 94010 )
     &      'Dimension mismatch.  AREA SOURCES file:', NROWS3D,
     &      'program:', NASRC
            CALL M3EXIT( 'GRDAMAT', 0, 0, MESG, 2 )
        END IF

        IF ( .NOT. READ3( ANAME, 'FIP', ALLAYS3, 0, 0,  IFIPS ) ) THEN
            CALL M3EXIT( 'GRDAMAT', 0, 0, 
     &      'Error reading variable "FIP" from AREA SOURCES file', 2 )
        END IF

        IF ( .NOT. READ3( ANAME, 'ASC7', ALLAYS3, 0, 0,  ISCC7 ) ) THEN
            CALL M3EXIT( 'GRDAMAT', 0, 0, 
     &      'Error reading variable "ASC7" from AREA SOURCES file', 2 )
        END IF

        IF ( .NOT. READ3( ANAME, 'ASC3', ALLAYS3, 0, 0,  ISCC3 ) ) THEN
            CALL M3EXIT( 'GRDAMAT', 0, 0, 
     &      'Error reading variable "ASC3" from AREA SOURCES file', 2 )
        END IF                                                              


C.......   Read the grid parameters from GRIDDESC for use in checking
C.......   EPS2.0-formatted surrogate files

        IF ( .NOT. DSCGRID( GRDNM, SCRBUF, GDTYP3D, 
     &              P_ALP3D, P_BET3D,P_GAM3D, XCENT3D, YCENT3D,
     &              XORIG3D, YORIG3D, XCELL3D, YCELL3D, 
     &              NCOLS3D, NROWS3D, NTHIK3D ) ) THEN

            SCRBUF = GRDNM
            MESG   = '"' // SCRBUF( 1:TRIMLEN( SCRBUF ) ) //
     &               '" not found in GRIDDESC file'
            CALL M3EXIT( 'GRDAMAT', 0, 0, MESG, 2 )

        END IF                  ! if dscgrid() failed

C.......   Read the ACTUAL ASCs file

        CALL M3MSG2( 'Reading ACTUAL ASCs file...' )

        CALL RDASCC( ADEV, NASCC, NASC, ASCA7, ASCA3 )

C.......   Read the GRIDDING XREF file

        CALL M3MSG2( 'Reading GRIDDING XREF file...' )

C........  Can read either EMS-95 _or_ EPS2.0 gridding cross-reference
C........  files.  The distinguishing factor is that EMS95 files do not
C........  have a blank space in column 6

        IREC    =  0
        I       =  0
        EFLAG   = .FALSE.
        FFORMAT = 'EPS2'   ! Default XREF input format

11      CONTINUE                        !  head of the XDEV-read loop

            READ( XDEV, 93000, END=22, IOSTAT=IOS ) LINE
            IREC = IREC + 1

            IF ( IOS .NE. 0 ) THEN

                WRITE( MESG,94010 ) 
     &              'I/O error', IOS, 
     &              'reading GRIDDING XREF file at line', IREC
                CALL M3MESG( MESG )
                EFLAG = .TRUE.
                GO TO  11

            END IF

C.............  Determine input format of cross-reference file
            IF( IREC .EQ. 1 .AND. LINE(6:6) .NE. ' ' ) FFORMAT = 'EMS95'

            IF( FFORMAT .EQ. 'EMS95' ) THEN

                SID = STR2INT( LINE(  1: 2 ) )
                CID = STR2INT( LINE(  3: 5 ) )
                ID7 = STR2INT( LINE(  6:12 ) )
                ID3 = STR2INT( LINE( 13:15 ) )
                SSC = STR2INT( LINE( 16:32 ) )
                FIP = 1000 * SID + CID

            ELSEIF( FFORMAT .EQ. 'EPS2' ) THEN

C.................  Skip point and mbile source ASCTs
                IF( LINE( 11:11 ) .GT. '9' .OR.
     &              LINE( 11:11 ) .EQ. ' ' .OR.
     &              LINE( 20:20 ) .EQ. ' '       ) GO TO 11  ! to head loop

                FIP = STR2INT( LINE(  1:5  ) )
                SSC = STR2INT( LINE(  7:9  ) )
                ID7 = STR2INT( LINE( 11:17 ) )
                ID3 = STR2INT( LINE( 18:20 ) )
                SID = FIP/1000
                CID = FIP - SID*1000

            ENDIF

C.............  Skip record if not in actual ASCS list.  This IF statement
C.............  cannot be merged with the one below because of the counter I
            IF ( ID7 .GT. 0 .AND. ID3 .GE. 0 .AND.
     &           FIND2( ID7, ID3, NASC, ASCA7, ASCA3 ) .LE. 0 ) GO TO 11

            I = I + 1

C.............  If input record doesn't have any blank/bad values, or
C.............  if there is not an overflow, then store record
            IF ( SID .LT. 0  .OR.
     &           CID .LT. 0  .OR.
     &           ID7 .LT. 0  .OR.
     &           ID3 .LT. 0  .OR.
     &           SSC .LT. 0  .OR.
     &           SSC .GT. NASRG   ) THEN

                WRITE( MESG,94010 )
     &              'Bad XREF record', IREC,
     &              'state',  SID,
     &              'county', CID,
     &              'ASC7',   ID7,
     &              'ASC3',   ID3,
     &              'code',   SSC
                CALL M3MESG( MESG )
                EFLAG = .TRUE.

            ELSE IF ( I .LE. NXMAX ) THEN 

                    INDEXA( I ) = I
                    AFIPSA( I ) = FIP
                    ASCX7A( I ) = ID7
                    ASCX3A( I ) = ID3
                    SSCXDA( I ) = SSC

            END IF

            GO TO  11           !  to head of loop

22      CONTINUE                !  end of the XDEV-read loop

        CALL M3MSG2( 'NOTE: File read in as ' // FFORMAT // ' format.' )

        NAREF = I
        WRITE( MESG,94010 ) 
     &      'SURG XREF records - actual:', NAREF,
     &      'dimensioned max:', NXMAX
        CALL M3MESG( MESG )

        IF ( EFLAG ) THEN
            CALL M3EXIT( 'GRDAMAT', 0, 0, 
     &                   'Error reading ACTUAL ASC file.', 2 )
        ELSE IF ( NAREF .GT. NXMAX ) THEN
            CALL M3EXIT( 'GRDAMAT', 0, 0, 
     &                   'Table overflow reading ACTUAL ASC file.', 2 )
        END IF


        CALL M3MSG2( 'Processing GRIDDING XREF file...' )

C.......   Sort input XREF by FIPs codes and ASC codes

        CALL SORTI3( NAREF, INDEXA, AFIPSA, ASCX7A, ASCX3A )

C.......   Store unsorted XREF as sorted lists

        NPRNF = 0
        NPRFC = 0
        DO 27 I = 1, NAREF

            J   = INDEXA( I )
            FIP = AFIPSA( J )
            ID7 = ASCX7A( J )
            ID3 = ASCX3A( J )
            SSC = SSCXDA( J )

            IF ( FIP .EQ. 0 ) THEN  !  FIP-independent

                IF ( NPRNF .LT. NASCC ) THEN
                    NPRNF = NPRNF + 1
                    ASCX70( NPRNF ) = ID7
                    ASCX30( NPRNF ) = ID3
                    SSCXD0( NPRNF ) = SSC

                ELSE
                    EFLAG = .TRUE.
                    WRITE( MESG,94010 )
     &                  'Max FIP-ind XREF exceeded for ID7', ID7,
     &                  'and ID3', ID3
                    CALL M3MESG( MESG )

                END IF

            ELSE                    !  FIP-dependent case
 
                IF ( NPRFC .LT. NXMAX ) THEN 
                    NPRFC = NPRFC + 1
                    AFIPS( NPRFC ) = FIP
                    ASCX7( NPRFC ) = ID7
                    ASCX3( NPRFC ) = ID3
                    SSCXD( NPRFC ) = SSC

                ELSE
                    EFLAG = .TRUE.
                    WRITE( MESG,94010 )
     &                  'Max FIP-dep XREF exceeded for FIP', FIP,
     &                  'ID7', ID7, 'and ID3', ID3
                    CALL M3MESG( MESG )

                END IF
 
            END IF      !  if fip zero, or not

27      CONTINUE

C.......   Initializations:

        DO 33  C = 1, NGRID
            NFIP( C ) = 0
33      CONTINUE

        DO 44  K = 1, NSMAX
        DO 43  J = 1, NASRG
            SFRC( J,K ) = BADVAL3
43      CONTINUE
44      CONTINUE

        LFIP  = -1
        F     =  0
        EFLAG = .FALSE.

        DO  55  S = 1, NASRC    !  initializations for Surrogates::ASC table

            FIP = IFIPS( S )

            IF ( FIP .NE. LFIP ) THEN

                LFIP = FIP
                F    = F + 1

                IF ( F .LE. NAFIP ) THEN

                    FIPLIST( F ) = FIP
                    SRCLIST( F ) = S

                END IF          !  if F overflows, else stays in bounds


            END IF              !  if encountered new fip

            ID7 = ISCC7( S )
            ID3 = ISCC3( S )

C...........   Get surrogate ID for this source; flag bad xrefs

            M = FIND3( FIP, ID7, ID3, NPRFC, AFIPS, ASCX7, ASCX3 )

            IF ( M .LT. 0 )  THEN           !  <FIP,ID7,ID3> not found

                M = FIND2( FIP, ID7, NPRFC, AFIPS, ASCX7 )

                IF ( M .LT. 0 )  THEN       !  <FIP,ID7> not found

                    M = FIND2( ID7, ID3, NPRNF, ASCX70, ASCX30 )

                    IF( M .LT. 0 ) THEN     !  <ID7,ID3> not found
                        EFLAG = .TRUE.
                        WRITE( MESG,94020 ) 
     &                      'SRG XREF not found for FIP', FIP,
     &                      'ASC', ID7, ID3
                        CALL M3MESG( MESG )
                        GO TO  55

                    ELSE
                        N = SSCXD0( M )
                    END IF                  !  <FIP,ID7> not found

                ELSE
                    N = SSCXD( M )
                END IF                      !  <FIP,ID7> not found

            ELSE
                N = SSCXD( M )
            END IF                          !  <FIP,ID7,ID3> not found

            IF ( N .LT. 0 )  THEN       !  xref never set.

                EFLAG = .TRUE.
                WRITE( MESG,94020 ) 
     &              'Bad SURG XREF for FIP', FIP,
     &              'ASC', ID7, ID3
                CALL M3MESG( MESG )

            ELSE

                SRGID( S ) = N

            END IF                      !  if <ID7,ID3> not found

55      CONTINUE        !  end loop constructing FIP list from actual sources

        IF ( F .GT. NAFIP ) THEN

            WRITE( MESG,94010 ) 
     &              'Overflow:  actual # of FIPS:', F,
     &              'dimensioned max:', NAFIP
            CALL M3EXIT( 'GRDAMAT', 0, 0, MESG, 2 )

        ELSE IF ( EFLAG ) THEN

            CALL M3EXIT( 'GRDAMAT', 0, 0, 
     &           'Error in XREF mapping sources to surrogates', 2 )

        END IF                  !  if overflow

        SRCLIST( F+1 ) = NASRC + 1
        NFIPLIST       = F


C...........   Read and process the surrogate coefficient file:
C...............   Must be sorted by FIP (state ID, county ID), surrogate ID,
C...............   col, and row.  
C...............   All actually-occurring area source FIP codes must be present.

        CALL M3MSG2( 'Reading GRIDDING SURROGATE COEFF file...' )

        FFORMAT = 'EPS2'

C........  Can read EMS-95 _or_ EPS2.0 formatted gridding surrogate
C........  files.  The distinguishing factor is that EMS95 files do not
C........  have a blank space in column 6

C...........   Initialize parameters for reading surrogates data
        IREC    =  0
        N       =  0
        EFLAG   = .FALSE.

C........  Determine input format of gridding surrogates file
C........  For EPS2 format, treat 1st line as header line, and QA file
C........  For Lat-lon grids, header will be irrelevant, but can't check either

        READ( SDEV, 93000, END=77, IOSTAT=IOS ) LINE

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
                    CALL M3WARN( 'GRDAMAT', 0, 0, MESG )

                    CALL M3MSG2( 'Resetting actual number of ' //
     &                           'surrogates to maximum.' )
 
                    NSRG = NASRG

                END IF

            ELSE  ! Grid type is LAT-LON

                MESG= 'Using LAT-LON grid, with EPS2.0 surrogates file.'
                CALL M3WARN( 'GRDAMAT', 0, 0, MESG )

            END IF

        ENDIF

C...........   Start loop for actually reading surrogates data

66      CONTINUE                !  head of the SDEV-read loop

C...............  Read EMS-95 formatted file
            IF( FFORMAT .EQ. 'EMS95' ) THEN

                READ( SDEV, *, END=77, IOSTAT=IOS ) 
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

                READ( SDEV, 93010, END=77, IOSTAT=IOS ) 
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

                DO 73 J = 1, NSRG
                    SSCT( J ) = J
   73           CONTINUE

            END IF         ! end input format type

            IF ( IOS .NE. 0 ) THEN

                WRITE( MESG,94010 ) 
     &              'I/O error ', IOS, 
     &              'reading SURG COEFF file at line', IREC
                CALL M3MESG( MESG )
                EFLAG = .TRUE.
                GO TO  66   ! to head of SDEV read loop

C.................  If columns and rows are still a problem, give message,
C.................  but don't error because may want to use bigger
C.................  surrogates file to do smaller domain
            ELSE IF ( COL .LT. 1  .OR. COL .GT. NCOLS ) THEN

                WRITE( MESG,94010 )
     &              'Column ', COL,
     &              'out of range in SURG COEFF file at line', IREC
                CALL M3MESG( MESG )
                GO TO  66   ! to head of SDEV read loop

            ELSE IF ( ROW .LT. 1  .OR. ROW .GT. NROWS ) THEN

                WRITE( MESG,94010 )
     &              'Row ', ROW,
     &              'out of range in SURG COEFF file at line', IREC
                CALL M3MESG( MESG )
                GO TO  66   ! to head of SDEV read loop

            END IF              !  if IOS bad, or col, or row out of range


C...........   Record surrogate fraction for this cell.
C...........   Loop through the number of surrogates in current record
            DO 64 J = 1, NSRG

                N = N + 1

                IF( SSCT( J ) .LT. 1     .OR. 
     &              SSCT( J ) .GT. NASRG      ) THEN

                    WRITE( MESG,94010 )
     &                  'Surrogate code ', SSCT( J ),
     &                  'out of range (NASRG) in ASPRO at line', IREC
                    CALL M3MESG( MESG )
                    EFLAG = .TRUE.
                    GO TO  66   ! to head of SDEV read loop

                ELSEIF( N .LE. NRMAX ) THEN

                    INDXA( N ) = N
                    CELLA( N ) = COL  +  NCOLS * ( ROW - 1 )
                    SSCSA( N ) = SSCT( J )
                    FIPSA( N ) = FIP
                    FRACA( N ) = SRGT( J )

                END IF                      ! if n in bounds

64          CONTINUE

            GO TO  66  !  to head of SDEV read loop

77      CONTINUE       !  end of the SDEV-read loop

        CALL M3MSG2( 'NOTE: File read in as ' // FFORMAT // ' format.' )

        NSREC = N
        WRITE( MESG,94010 ) 
     &      'Max SURG COEF records - actual:', N,
     &      'dimensioned:', NRMAX
        CALL M3MESG( MESG )

        IF ( N .GT. NRMAX ) THEN
            CALL M3EXIT( 'GRDAMAT', 0, 0, 
     &                   'Table overflow reading SURG COEFF file', 2 )
        END IF

        IF ( EFLAG ) THEN
            MESG = 'Error reading SURG COEF file.'
            CALL M3EXIT( 'GRDAMAT', 0, 0, MESG, 2 )
        END IF

C...........   Sort and process the surrogates:

        CALL M3MSG2 ( 'Processing GRIDDING SURROGATE COEFF file...' )

        CALL SORTI3( NSREC, INDXA, CELLA, FIPSA, SSCSA )

        F    = 0
        S    = 0
        T    = 0
        LCEL = IMISS3
        LFIP = IMISS3

        DO  99  J = 1, NSREC

            K   = INDXA( J )
            C   = CELLA( K )
            FIP = FIPSA( K )
            SSC = SSCSA( K )

            IF ( FIP .NE. LFIP  .OR.  C .NE. LCEL ) THEN

                S = FIND1( FIP, NFIPLIST, FIPLIST )

                IF ( S .GT. 0 ) THEN
                    LFIP = FIP
                    LCEL = C
                    F    = F + 1
                    NFIP( C ) = NFIP( C ) + 1
                    IF ( F .LE. NSMAX ) THEN
                        CFIP( F ) = S
                    END IF
                ELSE			!  else fip not found in inventory
                    WRITE( MESG,94010 ) 
     &                  'Surrogate FIP', FIP, 'not in inventory'
                    CALL M3MESG( MESG )
                    GO TO  99
                END IF		!  if fip found, or not

            END IF		!  if new FIP encountered

            IF ( S .GT. 0 ) THEN
                SFRC( SSC , F ) = FRACA( K )
            END IF

99      CONTINUE

        WRITE( MESG,94010 ) 
     &      'SURG COEF table - actual:', F,
     &      'dimensioned:', NSMAX
        CALL M3MESG( MESG )


        IF ( F .GT. NSMAX ) THEN
            CALL M3EXIT( 'GRDAMAT', 0, 0, 
     &      'Table overflow processing SURG COEFF file', 2 )
        END IF


C...........   Compute gridding matrix

        CALL M3MSG2( 'Computing gridding matrix...' )

        EFLAG = .FALSE.         !  flag to detect other errors
        LAST  = 0		!  coeff count for last cell
        N     = 0		!  CFIP counter
        T     = 0		!  coefficient counter
        CMAX  = 0		!  cell max:  # of coeffs
        CMIN  = NSREC		!  cell min:  # of coeffs

        DO  244  C = 1, NGRID

            L = NFIP( C )	!  # of FIPS for this cell

            DO  233  I = 1, L		!  process this fip-cell pair:

                N  = N + 1
                F  = CFIP( N )		!  index into FIPLIST, etc.
                LO = SRCLIST( F   )
                HI = SRCLIST( F+1 ) - 1

                DO  222  S = LO, HI	!  loop on AREA source subscript

                    SSC = SRGID( S )
                    SRG = SFRC( SSC, N )

                    IF ( SRG .GT. 0.0 ) THEN

                        T = T + 1		!  increment coeff counter

                        IF ( T .LT. NMATX ) THEN

                            IX( T ) = S
                            CX( T ) = SRG

                        END IF		!  if T in bounds

                    END IF		!  if SRG ok

222             CONTINUE		!  end loop on S:  sources for this FIP
                     
233         CONTINUE            !  end loop on I:  FIPs for this cell
            
            M       = T - LAST
            NX( C ) = M
            LAST    = T
            IF ( CMAX .LT. M )  THEN
                CMAX = M
            ELSE IF ( CMIN .GT. M )  THEN
                CMIN = M
            END IF

244     CONTINUE                !  end loop on cells C


C.......   Report statistics:

        WRITE( MESG,94010 ) 
     &      'Total number of coefficients   :', T   , CRLF()// BLANK5//
     &      'Max number of sources per cell :', CMAX, CRLF()// BLANK5//
     &      'Min number of sources per cell :', CMIN

        WRITE( MESG,94040 ) MESG( 1:TRIMLEN( MESG ) ) // 
     &      CRLF() // BLANK5 // 'Mean number of sources per cell:',
     &      FLOAT( T ) / FLOAT( NGRID )

        CALL M3MSG2( MESG )

C...........   Report errors, or if GFLAG, write out matrix to file:

        IF ( T .GT. NMATX ) THEN  ! overflow occurred:  just compute statistics.

            WRITE( MESG,94010 )
     &       'WARNING:  ' // GNAME( 1:TRIMLEN( GNAME ) ) //
     &       ' not written.' // CRLF() // BLANK5 //
     &       'Arrays would have overflowed.'// CRLF() // BLANK5 //
     &       'NMATX dimensioned:', NMATX, CRLF() // BLANK5 //
     &       'NMATX    required:', T

            CALL M3MSG2( MESG )
            CALL M3EXIT( 'GRDAMAT', 0,0, 'Bad executable program.', 3 )

        ELSE IF ( EFLAG ) THEN          ! other errors:  missing codes?

            CALL M3EXIT( 'GRDAMAT', 0,0, 'Missing surrogate codes', 3 )

        ELSE 			! open and write gridding matrix to file

C.......   Get the grid name and parameters to put into file header:
C.......   This is called above, but repeat here for safety from future mods.

            IF ( .NOT. DSCGRID( GRDNM, SCRBUF, GDTYP3D, 
     &                  P_ALP3D, P_BET3D,P_GAM3D, XCENT3D, YCENT3D,
     &                  XORIG3D, YORIG3D, XCELL3D, YCELL3D, 
     &                  NCOLS3D, NROWS3D, NTHIK3D ) ) THEN

                SCRBUF = GRDNM
                MESG   = '"' // SCRBUF( 1:TRIMLEN( SCRBUF ) ) //
     &                   '" not found in GRIDDESC file'
                CALL M3EXIT( 'GRDAMAT', 0, 0, MESG, 2 )

            END IF                  ! if dscgrid() failed

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
            DO  255  K = 2, MXDESC3
                FDESC3D( K ) = ' '
255         CONTINUE
            VNAME3D( 1 ) = 'AGRDMAT'
            UNITS3D( 1 ) = 'n/a'
            VDESC3D( 1 ) = 'NC Area Source gridding-coefficient matrix'
            VTYPE3D( 1 ) = M3REAL

            GNAME = PROMPTMFILE( 
     &        'Enter name for GRIDDING MATRIX output file, or "NONE"',
     &        FSUNKN3, 'AGMAT', 'GRDAMAT' )
     
            IF ( GNAME( 1:5 ) .NE. 'NONE ' ) THEN

                CALL M3MSG2( 'Writing out GRIDDING MATRIX file...' )

                IF ( .NOT. WRITE3( GNAME, 'ALL', 0, 0, NX ) ) THEN
                    CALL M3EXIT( 'GRDAMAT', 0, 0, 
     &                  'Error writing GRIDDING MATRIX file.', 2 )
                END IF
            END IF			!  if gname not "NONE"

        END IF                          !  if overflow, or if eflag, or if gflag

C...............   End of program

        CALL M3EXIT( 'GRDAMAT', 0, 0, 
     &               'SUCCESSFUL COMPLETION of program GRDAMAT', 0 )

C******************  FFORMAT  STATEMENTS   ******************************

C...........   Informational (LOG) message formats... 92xxx

92000   FORMAT( 5X, A )

92010   FORMAT( 5X, A, :, I12 )


C...........   Formatted file I/O formats............ 93xxx

93000   FORMAT( A )

93010   FORMAT( I5, 2F10.0, 10X, 100( F10.8 ) )


C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I10, :, 2X ) )

94020   FORMAT( 10( A, :, I6.5, :, 1X, A, :, I8.7, I3.3, :, 1X, A ) )

94030   FORMAT( A, :, I4, :, 2X, 
     &          A, I7, :, 2X, 
     &          A, :, I9.7, I3.3, :, 2X, A )

94040   FORMAT( A, :, F10.2 )

        END
