
        PROGRAM SPCAMAT

C***********************************************************************
C  program body starts at line 231
C
C  DESCRIPTION:
C       Construct area source chemical speciation matrix from data 
C       contained in EPS-style chem-split and chem-xref files.
C       Uses fixed NOtoNOX  factor of 0.0197  moles  NO/gm NOX (from EPS)
C       Uses fixed NO2toNOX factor of 0.00217 moles NO2/gm NOX  "
C
C  PRECONDITIONS REQUIRED:
C       Sorted, cut-down input data for chem-split and chem-xref files;
C       M3IO area-sources file from program GRDAMAT
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C       RDASCC(), RDSPROF(), RDACONV()
C       Models-3 I/O
C
C  REVISION  HISTORY:
C       Prototype  3/95 by CJC.
C
C       Version 12/96 by CJC with adaptations toward generalized chemistry
C       mechanism treatment for use with EMS-95-derived inventories.
C       Optional volatile conversion file.  Now uses RDSPROF, RDACONV, RDASCC.
C
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
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
        INCLUDE 'FDESC3.EXT'    !  I/O API file description data structures.

C...........   EXTERNAL FUNCTIONS and their descriptions:

        CHARACTER*2     CRLF
        INTEGER         FIND1, FIND2, FIND3, INDEX1
        INTEGER         LBLANK
        LOGICAL         ENVYN
        LOGICAL         GETYN
        INTEGER         PROMPTFFILE
        CHARACTER*16    PROMPTMFILE
        INTEGER         RDSPROF
        INTEGER         STR2INT
        INTEGER         TRIMLEN

        EXTERNAL    CRLF, FIND1, FIND2, FIND3, ENVYN, GETYN, INDEX1, 
     &              LBLANK, PROMPTFFILE, PROMPTMFILE, RDSPROF, STR2INT,
     &              TRIMLEN

C...........   PARAMETERS and their descriptions:

        CHARACTER*5     BLANK5
        REAL            TON2GM
        PARAMETER (     BLANK5 = ' ',
     &                  TON2GM = 907184.74 )

C...........   LOCAL VARIABLES and their descriptions:
C...........   NOTE that ASC (area-source-category) ID's are 10-digit
C...........   unsigned integers which may be treated as a leading 7-digit
C...........   field, and a trailing 3-digit field.  *7 and *3 arrays below
C...........   follow this scheme with parallel arrays

        INTEGER     FIP
        INTEGER     PLANT
        INTEGER     SID
        INTEGER     ID7
        INTEGER     ID3

C...........   Speciation profiles table

        REAL        VREACT
        REAL        ETHADJ
        REAL        FRMADJ
        REAL        ALLADJ

        INTEGER     NPROFS                        !  actual number of profiles
        INTEGER     SXREF( MXSPRO )               !  profile indexes (sorted)
        REAL        SCOEF( MXSPRO, MXMPOL )       !  profile coefficients

        CHARACTER*7 ID7C
        CHARACTER*3 ID3C

C.......   Actually-occurring ASC table

        INTEGER     NASC
        INTEGER     ASCA7( NASCC )
        INTEGER     ASCA3( NASCC )

C...........   Sorted Volatile Pollutant Type ("ROG-to-TOG" Conversion Table:
C........ Two parts.  First part depends on scc, state fip code

        INTEGER     NCST             !  actual table-size
        INTEGER     FIPST( MXSCNV )  !  conversion factor STATE ID
        INTEGER     ASC7S( MXSCNV )  !  conversion factor ASC leading-7 
        INTEGER     ASC3S( MXSCNV )  !  conversion factor ASC trailing-3
        REAL        CONVS( MXSCNV )  !  conversion factor

C........ Second part depends on scc only (independent of FIP)

        INTEGER     NCNN             !  actual table-size
        INTEGER     ASC7N( MXSCNV )  !  conversion factor ASC leading-7 
        INTEGER     ASC3N( MXSCNV )  !  conversion factor ASC trailing-3
        REAL        CONVN( MXSCNV )  !  conversion factor

C...........   Speciation cross-reference table:  
C.......   Before sorting:

        INTEGER     NREFA
        INTEGER     INDEXA( NASRC )
        INTEGER     FIPSXA( NASRC )
        INTEGER     ASC7XA( NASRC )
        INTEGER     ASC3XA( NASRC )
        INTEGER     SPPRFA( NASRC )
        REAL        VREACA( NASRC )
        REAL        ETHADA( NASRC )
        REAL        FRMADA( NASRC )
        REAL        ALLADA( NASRC )

C.......   After sorting -- three parts:   first part is FIP-independent part

        INTEGER     NREF0               !  number of actual entries
        INTEGER     ASC7X0( NASCC )     !  leading  7 digits of ASC
        INTEGER     ASC3X0( NASCC )     !  trailing 3 digits of ASC
        INTEGER     SPPRF0( NASCC )     !  XREF for lookup in ASPRO
        REAL        VREAC0( NASCC )     !  
        REAL        ETHAD0( NASCC )     !  
        REAL        FRMAD0( NASCC )     !  
        REAL        ALLAD0( NASCC )     !  

C.......   Second part is FIP-state-code-dependent part

        INTEGER     NREFS
        INTEGER     FIPSXS( NASCC )
        INTEGER     ASC7XS( NASCC )
        INTEGER     ASC3XS( NASCC )
        INTEGER     SPPRFS( NASCC )
        REAL        VREACS( NASCC )
        REAL        ETHADS( NASCC )
        REAL        FRMADS( NASCC )
        REAL        ALLADS( NASCC )

C.......   Third part is FIP-dependent part

        INTEGER     NREFF
        INTEGER     FIPSXF( NASCC )
        INTEGER     ASC7XF( NASCC )
        INTEGER     ASC3XF( NASCC )
        INTEGER     SPPRFF( NASCC )
        REAL        VREACF( NASCC )
        REAL        ETHADF( NASCC )
        REAL        FRMADF( NASCC )
        REAL        ALLADF( NASCC )

C.......   Speciation matrix variable names:

        REAL        SPCFACS( NASRC, NMPOL )    !  speciation coefficients

C...........   Area Sources Table

        INTEGER      FIPS( NASRC )  !  source FIPS (county) ID
        INTEGER      ASC7( NASRC )  !  leading-7  digits of source ASC
        INTEGER      ASC3( NASRC )  !  trailing-3 digits of source ASC

C...........   Other local variables
        
        INTEGER         S, F, I, J, V, K!  loop counters.
        INTEGER         ADEV            !  for actual-ASC file
        INTEGER         KDEV            !  volatile form conversion file
        INTEGER         SDEV            !  for surrogate coeff file
        INTEGER         XDEV            !  for surrogate xref  file
        INTEGER         LDEV            !  log-device
        INTEGER         SPCPRF          !  scratch xref value
        INTEGER         PRF0            !  subscript for profile 0
        INTEGER         ECNT, XCNT      !  error count, XREF-missing count
        INTEGER         IOS             !  I/O status
        INTEGER         IREC            !  input line (record) number
        INTEGER         I0              !  for counting FIP-independent srcs
        INTEGER         I1              !  for counting STATE-dependent srcs
        INTEGER         I2              !  for counting FIP-dependent   srcs
        INTEGER         ICO             !  sub. for CO   in EMNAM
        INTEGER         INO             !  sub. for NO   in EMNAM
        INTEGER         INO2            !  sub. for NO2  in EMNAM
        INTEGER         IETH            !  sub. for ETH  in EMNAM
        INTEGER         IFRM            !  sub. for FORM in EMNAM
        INTEGER         IORG            !  smallest sub. for VOC's
        INTEGER         VID             !  id index of volatile pollutant
 
        REAL            FAC, T          !  common factor subexpressions

        CHARACTER*8     INVPRO  !  volatile inventory pollutant in spec profile
        CHARACTER*8     POLV    !  volatile inventory pollutant in inventory
        CHARACTER*16    ANAME   !  logical name for area-source  input file
        CHARACTER*16    MNAME   !  logical name for spec matrix output file
        CHARACTER*256   MESG    !  message buffer for M3MESG() and M3EXIT()
        
        LOGICAL         EFLAG   !  input error flag
        LOGICAL         KFLAG   !  if volatile format conversion factor or not
        LOGICAL         VFLAG   !  verbose conversion messages

C***********************************************************************
C   begin body of program SPCAMAT

        LDEV = INIT3()

        CALL INITEM( LDEV )

        WRITE( *,92000 ) 
     &  ' ',
     &  'Program SPCAMAT to take the netCDF AREA SOURCES file, the',
     &  'sorted SPECIATION PROFILES file, and the sorted SPECIATION',
     &  'CROSS-REFERENCE file, and produce the SPECIATION MATRIX ',
     &  'file.',
     &  ' ',
     &  'You will need to enter the logical names for the input and',
     &  'output files (and to have set them prior to program launch,',
     &  'using "setenv <logicalname> <pathname>").  Input files must',
     &  'have been sorted as indicated, prior to program execution.',
     &  ' ',
     &  'You may use END_OF-FILE (control-D) to quit the program',
     &  'during logical-name entry. Default responses are indicated',
     &  'in brackets [LIKE THIS].',
     &  ' '

        IF ( .NOT. GETYN( 'Continue with program?', .TRUE. ) ) THEN
            CALL M3EXIT( 'SPCAMAT', 0, 0, 'Ending program SPCAMAT', 2 )
        END IF

C.......   Get file name; open input area sources file

        ANAME = PROMPTMFILE( 
     &          'Enter logical name for AREA SOURCE input file',
     &          FSREAD3, 'AREA', 'SPCAMAT' )

        ADEV = PROMPTFFILE( 
     &           'Enter logical name for ACTUAL ASC file',
     &           .TRUE., .TRUE., 'ASCS', 'SPCAMAT' )

        XDEV = PROMPTFFILE( 
     &           'Enter logical name for SPECIATION XREF file',
     &           .TRUE., .TRUE., 'ASREF', 'SPCAMAT' )

        SDEV = PROMPTFFILE( 
     &           'Enter logical name for SPECIATION PROFILES file',
     &           .TRUE., .TRUE., 'ASPRO', 'SPCAMAT' )

        VFLAG = ENVYN( 'SPCAMAT_VERBOSE', 
     &                 'Verbose reporting of XXX-to-VOC conversion',
     &                 .FALSE., IOS )

C...........   Read in emissions indexing data IFIP, ASC7, ASC3.

        IF ( .NOT. DESC3( ANAME ) ) THEN
            CALL M3EXIT( 'SPCAMAT', 0, 0, 
     &      'Error reading HEADER from AREA SOURCES file', 2 )
        ELSE IF ( NROWS3D .NE. NASRC ) THEN
            WRITE( MESG, 94010 )
     &      'Dimension mismatch.  AREA SOURCES file:', NROWS3D,
     &      'program:', NASRC
            CALL M3EXIT( 'SPCAMAT', 0, 0, MESG, 2 )
        END IF

        IF ( .NOT. READ3( ANAME, 'FIP', ALLAYS3, 0, 0, FIPS ) ) THEN
            CALL M3EXIT( 'SPCAMAT', 0, 0, 
     &                   'Error reading FIP from file' // ANAME, 2 )
        END IF

        IF ( .NOT. READ3( ANAME, 'ASC7', ALLAYS3, 0, 0, ASC7 ) ) THEN
            CALL M3EXIT( 'SPCAMAT', 0, 0, 
     &                   'Error reading ASC7 from file' // ANAME, 2 )
        END IF

        IF ( .NOT. READ3( ANAME, 'ASC3', ALLAYS3, 0, 0, ASC3 ) ) THEN
            CALL M3EXIT( 'SPCAMAT', 0, 0, 
     &                   'Error reading ASC3 from file' // ANAME, 2 )
        END IF

C.......   Read the ACTUAL ASC FILE
        
        CALL M3MSG2( 'Reading ACTUAL ASC file...' )
        CALL RDASCC( ADEV, NASCC, NASC, ASCA7, ASCA3 )

C.......   Read the SPECIATION PROFILES FILE

        CALL M3MSG2( 'Reading SPECIATION PROFILES file...' )

        NPROFS = RDSPROF( SDEV , NMPOL, EMNAM, MXSPRO, 
     &                    SXREF, SCOEF, INVPRO )

        IF ( NPROFS .LE. 0 ) THEN
            CALL M3EXIT( 'SPCAMAT', 0, 0,
     &                   'Error reading SPECIATION PROFILE file.', 2 )
        END IF

C.......   Compare inventory pollutant in the speciation profile with
C.......   volatile pollutant in inventory.

        VID  = 0
        I    = 0
        POLV = ' '
77      CONTINUE  ! Loop to find volatile inventory pollutant

            I = I + 1
            VID = INDEX1( EINAM( I ), NUMVOL, VOLTYP )
            IF( VID .GT. 0 ) THEN
 
                POLV = EINAM( I )
                GO TO 88
 
            ELSEIF( I .EQ. NIPOL ) THEN
 
                I = I + 1
                CALL M3WARN( 'SPCMMAT', 0, 0,
     &                       'No valid volatile inventory pollutant' )
                GO TO 88

            ENDIF
 
        GO TO 77  ! To head of finding volatile pollutant type
 
88      CONTINUE

        KFLAG = ( INVPRO .NE. POLV .AND. I .LE. NIPOL )

        IF( KFLAG ) THEN

            MESG = 'Enter logical name for ' //
     &              POLV ( 1:TRIMLEN( POLV ) ) // ' TO ' //
     &              INVPRO ( 1:TRIMLEN( INVPRO ) ) // 
     &             ' CONVERSION file'

            KDEV = PROMPTFFILE( 
     &             MESG, .TRUE., .TRUE., 'ASCNV', 'SPCAMAT' )

            CALL RDACONV( KDEV, FIPST, ASC7S, ASC3S, CONVS, NCST,
     &                                 ASC7N, ASC3N, CONVN, NCNN )

            IF ( NCST .LT. 0  .OR.  NCNN .LT. 0 ) THEN
                CALL M3EXIT( 'SPCAMAT', 0, 0,
     &                       'Error in RDACONV()', 2 )
            ELSE IF ( NCST .EQ. 0  .AND. NCNN .EQ. 0 ) THEN
                CALL M3MSG2( 'ASCNV file empty' )
                KFLAG = .FALSE.
            END IF

        END IF                  !  if invpro does not match polv

C.......   Read and process the speciation cross-reference file

        CALL M3MSG2( 'Reading SPECIATION XREF file...' )

        IREC  = 0
        I     = 0
        EFLAG = .FALSE.

101     CONTINUE                !  head of the XDEV-read loop

            IREC = IREC + 1
            READ( XDEV, 93030, END=111, IOSTAT=IOS ) 
     &          FIP, PLANT, ID7C, ID3C, SID, VREACT, ETHADJ, 
     &          FRMADJ, ALLADJ

            IF( IOS .NE. 0 ) THEN

                EFLAG = .TRUE.
                WRITE( MESG,94010 ) 
     &              'Error ', IOS, 
     &              'reading SPECIATION XREF file at line', IREC
                CALL M3MESG( MESG )
                GO TO 101

            ELSE IF( PLANT .GT. 0 ) THEN        ! Point source x-ref entry

                GO TO 101

            ELSE IF( ID7C(1:1) .GT. '9' ) THEN  ! Motor vehicle x-ref entry

                GO TO 101

            ELSE

                ID7 = STR2INT( ID7C )
                ID3 = STR2INT( ID3C )

                IF( ID7 .GT. 0  .AND.  ID3 .GE. 0 ) THEN
                    J = FIND2( ID7, ID3, NASC, ASCA7, ASCA3 )
                    IF( J .LE. 0 ) THEN   ! ASC code not found in ASCS LIST file
                        GO TO 101
                    END IF
                ELSE IF ( ID7 .LT. 0  .OR.  ID3 .LT. 0 ) THEN
                    EFLAG = .TRUE.
                    WRITE( MESG,94010 ) 
     &                  'Bad ASC:  ID7=', ID7, 'ID3=', ID3,
     &                  'in SPECIATION XREF file at line', IREC
                    CALL M3MESG( MESG )
                    GO TO 101
                END IF

                I = I + 1
                IF( I .LE. NASRC ) THEN

                    INDEXA( I ) = I
                    FIPSXA( I ) = FIP
                    ASC7XA( I ) = ID7
                    ASC3XA( I ) = ID3
                    SPPRFA( I ) = SID
                    VREACA( I ) = VREACT
                    ETHADA( I ) = ETHADJ
                    FRMADA( I ) = FRMADJ
                    ALLADA( I ) = ALLADJ

                END IF                  !  No overflow, ASC found in list

            END IF                      !  Area source record

            GO TO  101                  !  to head of the XDEV-read loop

111     CONTINUE                        !  end of the XDEV-read loop

        IF ( I .GT. NASRC ) THEN
            WRITE( MESG,94010 )
     &          'XREF table overflow; actual', I,
     &          'dimensioned NASRC=', NASRC
            CALL M3EXIT( 'SPCAMAT', 0, 0, MESG, 2 )
        ELSE IF ( EFLAG ) THEN
            CALL M3EXIT( 'SPCAMAT', 0,0, 'Error reading XREF file', 2 )
        END IF

        NREFA = I


C...........   Now sort the XREF inputs and place into appropriate tables:

        CALL SORTI3( NREFA, INDEXA, FIPSXA, ASC7XA, ASC3XA )

        CALL M3MSG2( 'Processing SPECIATION XREF file...' )

        PRF0 = FIND1( 0, NPROFS, SXREF )		!  fallback profile

        I0    = 0
        I1    = 0
        I2    = 0
        ECNT  = 0
        EFLAG = .FALSE.

        DO  122  I = 1, NREFA

            J   = INDEXA( I )
            FIP = FIPSXA( J ) 
            ID7 = ASC7XA( J )
            ID3 = ASC3XA( J )
            SID = SPPRFA( J )

C...........   Find the correct SPECIATION PROFILE for this cross-reference:
          
            K = FIND1( SID, NPROFS, SXREF )

            IF ( K .LE. 0 ) THEN        !  error:  XREF index not found

                ECNT   = ECNT + 1
                WRITE( MESG,94020 ) 
     &              'FIP', FIP,
     &              'ASC', ID7, ID3, 
     &              'PROF ID', SID,
     &              'not found; using profile 0000'
                CALL M3MESG( MESG )

                K = PRF0

            END IF              !  if k nonpositive (cross-reference not found)

            IF ( FIP .EQ. 0  ) THEN             !  FIP independent refs

                I0 = I0 + 1
                IF ( I0 .LE. NASCC ) THEN

                    SPPRF0( I0 ) = K
                    ASC7X0( I0 ) = ID7 
                    ASC3X0( I0 ) = ID3
                    VREAC0( I0 ) = VREACA( J )
                    ETHAD0( I0 ) = ETHADA( J )
                    FRMAD0( I0 ) = FRMADA( J )
                    ALLAD0( I0 ) = ALLADA( J )

                END IF

            ELSE IF ( MOD( FIP, 1000 ) .EQ. 0 )  THEN ! state-dependent refs

                I1 = I1 + 1
                IF ( I1 .LE. NASCC ) THEN

                    FIPSXS( I1 ) = FIP
                    SPPRFS( I1 ) = K
                    ASC7XS( I1 ) = ID7 
                    ASC3XS( I1 ) = ID3
                    VREACS( I1 ) = VREACA( J )
                    ETHADS( I1 ) = ETHADA( J )
                    FRMADS( I1 ) = FRMADA( J )
                    ALLADS( I1 ) = ALLADA( J )

                END IF

            ELSE                        ! fip > 0:  FIP-dep refs

                I2 = I2 + 1
                IF ( I2 .LE. NASCC ) THEN

                    FIPSXF( I2 ) = FIP
                    SPPRFF( I2 ) = K
                    ASC7XF( I2 ) = ID7 
                    ASC3XF( I2 ) = ID3
                    VREACF( I2 ) = VREACA( J )
                    ETHADF( I2 ) = ETHADA( J )
                    FRMADF( I2 ) = FRMADA( J )
                    ALLADF( I2 ) = ALLADA( J )

                END IF

            END IF		!  if fip-independent, or state-dep, or fip-dep

122     CONTINUE                !  end of the X-REF allocation loop

        IF ( ECNT .GT. 0 ) THEN

            WRITE( MESG,94010 ) 
     &      'Number of missing profiles from XREF:  ', ECNT,
     &      CRLF() // BLANK5 // '   Default profile 0000 used.'
            CALL M3MSG2( MESG )

            IF( PRF0 .LE. 0 ) THEN
                EFLAG = .TRUE.
                CALL M3MSG2( 'Default profile 0000 not found' )
            END IF                   !  if i<0

        END IF		!  if ecnt > 0 (defaults substituted)

        IF ( I0 .GT. NASCC ) THEN
            EFLAG = .TRUE.
            WRITE( MESG,94010 )
     &          'Max FIP-independent table-size (NASCC)=', NASCC, 
     &          CRLF()// BLANK5// 'exceeded in XREF file.  Needed:', I0
            CALL M3MSG2( MESG )
        END IF

        IF ( I1 .GT. NASCC ) THEN
            EFLAG = .TRUE.
            WRITE( MESG,94010 )
     &          'Max state-dependent table-size ', NASCC,  
     &          CRLF()// BLANK5// 'exceeded in XREF file.  Needed:', I1
            CALL M3MSG2( MESG )
        END IF

        IF ( I2 .GT. NASCC ) THEN
            EFLAG = .TRUE.
            WRITE( MESG,94010 )
     &          'Max FIP-dependent table-size ', NASCC,  
     &          CRLF()// BLANK5// 'exceeded in XREF file.  Needed', I2
            CALL M3MSG2( MESG )
        END IF

        IF ( EFLAG ) THEN
            CALL M3EXIT( 'SPCAMAT', 0, 0, 
     &                   'Error processing XREF file', 2 )
        END IF

        IETH  = INDEX1( 'ETH ', NMPOL, EMNAM )
        IFRM  = INDEX1( 'FORM', NMPOL, EMNAM )
        ICO   = INDEX1( 'CO  ', NMPOL, EMNAM )
        INO   = INDEX1( 'NO  ', NMPOL, EMNAM )
        INO2  = INDEX1( 'NO2 ', NMPOL, EMNAM )
        IORG  = 1 + MAX( ICO, INO, INO2 )
        NREF0 = I0              !  number of location-independent refs
        NREFS = I1              !  " state-only         dependent refs
        NREFF = I2              !  " state-and-county   dependent refs
        EFLAG = .FALSE.
        ECNT  = 0
        XCNT  = 0

C.......   Compute the speciation coefficient matrix
C...........   ASSUMPTION:  everything sorted by FIP and ASC, or by XREF ID
C...........   Two stage process:  first try for cross-referece by FIP and ASC
C...........   Then if that fails, try by ASC only.

        DO  155  S = 1, NASRC

            FIP = FIPS( S )
            ID7 = ASC7( S )
            ID3 = ASC3( S )

            F = FIND3( FIP, ID7, ID3, NREFF, FIPSXF, ASC7XF, ASC3XF )

            IF ( F .GT. 0 ) THEN                !  county-scc match:

                SPCPRF = SPPRFF( F )
                ALLADJ = ALLADF( F )
                FRMADJ = FRMADF( F )
                ETHADJ = ETHADF( F )
                VREACT = VREACF( F )

            ELSE                !  no county-scc match; try state-scc match

                SID = 1000 * ( FIP / 1000 )             !  state code
                F   = FIND3( SID, ID7, ID3, 
     &                       NREFS, FIPSXS, ASC7XS, ASC3XS )

                IF ( F .GT. 0 ) THEN

                    SPCPRF = SPPRFS( F )
                    ALLADJ = ALLADS( F )
                    FRMADJ = FRMADS( F )
                    ETHADJ = ETHADS( F )
                    VREACT = VREACS( F )

                ELSE                    !  no state-scc match; try just scc

                    F = FIND2( ID7, ID3, NREF0, ASC7X0, ASC3X0 )

                    IF ( F .GT. 0 ) THEN

                        SPCPRF = SPPRF0( F )
                        ALLADJ = ALLAD0( F )
                        FRMADJ = FRMAD0( F )
                        ETHADJ = ETHAD0( F )
                        VREACT = VREAC0( F )

                    ELSE                !  not found:  error

                        WRITE( MESG,94020 ) 
     &                      'Source:  FIP', FIP,
     &                      'ASC', ID7, ID3, 'not found in XREF'
                        CALL M3MESG( MESG )

                        XCNT   = XCNT + 1
                        SPCPRF = PRF0
                        ALLADJ = 1.0
                        FRMADJ = 0.0
                        ETHADJ = 1.0
                        VREACT = 1.0

                    END IF              !  if find2(scc) succeeded, or not

                END IF          !  if find3(state,scc) succeeded, or not

            END IF              !  if find3(county,scc) succeeded, or not

C...........   Processing for ORGANICS:
C...........   If KFLAG, find the the correct VOLATILE CONVERSION FACTOR 
C...........   and apply it via FAC -- a two stage process:
C...........       either by  FIP(state),SCC
C...........       or else by SCC only

            FAC = TON2GM * ALLADJ
            
            IF ( ICO  .GT. 0 ) SPCFACS( S,ICO  ) = FAC * COFAC
            IF ( INO  .GT. 0 ) SPCFACS( S,INO  ) = FAC * NOFAC
            IF ( INO2 .GT. 0 ) SPCFACS( S,INO2 ) = FAC * NO2FAC

            IF( KFLAG ) THEN

                SID = ( FIP / 1000 ) * 1000
                F   = FIND3( SID, ID7, ID3, NCST, FIPST, ASC7S, ASC3S )

                IF ( F .GT. 0 ) THEN

                    T   = CONVS( F )
                    FAC = T * FAC

                ELSE                            ! no state,scc match

                    F = FIND2( ID7, ID3, NCNN, ASC7N, ASC3N )

                    IF ( F .GT. 0 ) THEN

                        T   = CONVN( F )
                        FAC = T * FAC

                    ELSE 			!  no match at all

                        IF ( VFLAG ) THEN !  no match: WARNING; leave FAC alone

                            WRITE( MESG,94020 )
     &                       'Default ' // POLV( 1:TRIMLEN( POLV ) ) // 
     &                       ' to ' // INVPRO ( 1:TRIMLEN( INVPRO ) )// 
     &                       ' of 1.0 used for FIP= ', FIP, 
     &                       'ASC=',  ID7, ID3
                            CALL M3MESG( MESG )

                        END IF          ! if vflag:  verbose

                    END IF              ! scc only match or not

                END IF                  ! state and scc match or not

            END IF                      ! KFLAG is true, or not

            DO  144  V = IORG, NMPOL		!  loop on organics
                SPCFACS( S,V ) = FAC * SCOEF( SPCPRF,V )
144         CONTINUE

            IF( FRMADJ .NE. 0. ) THEN
                CALL M3MESG(
     &               'WARNING: Formaldehyde adjustment is non-zero' )
            ENDIF

            IF( IFRM .GT. 0 ) 
     &        SPCFACS( S, IFRM ) = ( 1. + FRMADJ ) * SPCFACS( S, IFRM )
            IF( IETH .GT. 0 ) 
     &        SPCFACS( S, IETH ) = ETHADJ * SPCFACS( S, IETH )

155     CONTINUE        ! end loop on area sources, computing speciation matrix

        IF ( XCNT .GT. 0 ) THEN
            WRITE( MESG,94010 ) 
     &          'Number of sources not found in XREF:  ', XCNT
            CALL M3MSG2( MESG ) 
        END IF

        IF ( EFLAG ) THEN
            CALL M3EXIT( 'SPCAMAT', 0, 0, 
     &                   'Error computing SPECIATION MATRIX', 2 )
        END IF

C.......   Open speciation matrix file:  reuse much of AREA file description

        IF ( .NOT. DESC3( ANAME ) ) THEN
            MESG = 'Could not get description of file ' // ANAME
            CALL M3EXIT( 'SPCAMAT', 0, 0, MESG, 2 )
        END IF

        NVARS3D = NMPOL
        FTYPE3D = GRDDED3

        DO  202  V = 1, NMPOL
            I = INDEX( SFNAM( V ), '_to_' )
            VNAME3D( V ) = SFNAM( V )
            UNITS3D( V ) = 'moles/ton'
            VDESC3D( V ) = 'Conversion factor, ' // 
     &              SFNAM( V )( 1:I-1 ) // ' to ' // 
     &              SFNAM( V )( I+4:16 )
            VTYPE3D( V ) = M3REAL
202     CONTINUE

        FDESC3D( 1 ) = 
     &  'Matrix of speciation factors, by source, for converting '
        FDESC3D( 2 ) = 
     &  'inventory species into model species for the CB-IV chemical'
        FDESC3D( 3 ) = 
     &  'mechanism.'
        DO  211  V = 4, MXDESC3
            FDESC3D( V ) = ' '
211     CONTINUE

        MNAME = PROMPTMFILE( 
     &  'Enter logical name for SPECIATION MATRIX output file',
     &          FSUNKN3, 'ASMAT', 'SPCAMAT' )


C.......   Write out the speciation matrix:

        CALL M3MSG2( 'Writing out SPECIATION MATRIX...' )

        IF ( .NOT. WRITE3( MNAME, 'ALL', 0, 0, SPCFACS ) ) THEN
            CALL M3EXIT( 'SPCAMAT', 0, 0, 
     &                   'Error writing SPECIATION MATRIX file.', 2 )
        END IF

      CALL M3EXIT( 'SPCAMAT', 0, 0, 'Normal completion.', 0 )

C******************  FORMAT  STATEMENTS   ******************************

C...........   Informational (LOG) message formats... 92xxx

92000   FORMAT( 5X, A )


C...........   Formatted file I/O formats............ 93xxx

93020   FORMAT( I7, I3 )

93030   FORMAT( I5, I6, 5X, A7, A3, 1X, I5, F6.3, 15X, 3 F6.3 ) ! for ASREF


C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10 ( A, :, I10, :, 2X ) )

94020   FORMAT( A, 2X, I5.5, 2X, A, 2X, I7.7, I3.3 , 
     &          10 ( :, 2X, A, :, 2X, I5 ) )

        END

