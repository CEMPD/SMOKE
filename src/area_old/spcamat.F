
C Version "@(#)$Id$ $Source$ $Date$ 

        PROGRAM SPCAMAT

C***********************************************************************
C  program body starts at line 155
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
C       Models-3 I/O
C       FIND1, FIND2, FIND3, GETEFILE, GETYN, TRIMLEN
C
C  REVISION  HISTORY:
C       Prototype  3/95 by CJC.
C
C***********************************************************************

      IMPLICIT NONE

C...........   INCLUDES:
        
        INCLUDE 'ARDIMS3.EXT'   !  area-source dimensioning parameters
        INCLUDE 'CHDIMS3.EXT'   !  emis chem parms (inventory + model)
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
        INCLUDE 'FDESC3.EXT'    !  I/O API file description data structures.


C...........   EXTERNAL FUNCTIONS and their descriptions:
        
        INTEGER         FIND1, FIND2, FIND3
        INTEGER         GETEFILE
        LOGICAL         GETYN
        INTEGER         PROMPTFFILE
        CHARACTER*16    PROMPTMFILE
        INTEGER         TRIMLEN
        
        EXTERNAL        FIND1, FIND2, FIND3, GETEFILE, GETYN, 
     &                  PROMPTFFILE, PROMPTMFILE, TRIMLEN


C...........   PARAMETERS and their descriptions:

        REAL        TON2GM
        PARAMETER ( TON2GM = 907184.74 )
        
C...........   LOCAL VARIABLES and their descriptions:
C...........   NOTE that ASC (area-source-category) ID's are 10-digit
C...........   unsigned integers which may be treated as a leading 7-digit
C...........   field, and a trailing 3-digit field.  *7 and *3 arrays below
C...........   follow this scheme with parallel arrays

        INTEGER     FIP
        INTEGER     SID
        INTEGER     ID7
        INTEGER     ID3
                                 
C...........   Speciation profiles table
        
        REAL        FACS( NMPOL )    !  input record:  speciation coefficients
        REAL        VREACT
        REAL        ETHADJ
        REAL        FRMADJ
        REAL        ALLADJ
        
        INTEGER     NPROFS                      !  actual number of profiles
        INTEGER     SXREF( NSPRA )              !  profile indexes (sorted)
        REAL        SCOEF( NSPRA, NMPOL )       !  profile coefficients

C...........   Speciation cross-reference table:  
C.......   three parts:   first part is FIP-independent part
            
        INTEGER     NREF0               !  number of actual entries
        INTEGER     ASC7X0( NASCC )	!  leading  7 digits of ASC
        INTEGER     ASC3X0( NASCC )	!  trailing 3 digits of ASC
        INTEGER     SPXRF0( NASCC )	!  XREF for lookup in ASPRO
        REAL        VREAC0( NASCC )	!  
        REAL        ETHAD0( NASCC )	!  
        REAL        FRMAD0( NASCC )	!  
        REAL        ALLAD0( NASCC )	!  

C.......   Second part is FIP-state-code-dependent part

        INTEGER     NREFS
        INTEGER     FIPSXS( NASCC )
        INTEGER     ASC7XS( NASCC )
        INTEGER     ASC3XS( NASCC )
        INTEGER     SPXRFS( NASCC )
        REAL        VREACS( NASCC )
        REAL        ETHADS( NASCC )
        REAL        FRMADS( NASCC )
        REAL        ALLADS( NASCC )

C.......   Third part is FIP-dependent part

        INTEGER     NREFF
        INTEGER     FIPSXF( NASCC )
        INTEGER     ASC7XF( NASCC )
        INTEGER     ASC3XF( NASCC )
        INTEGER     SPXRFF( NASCC )
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
        
        INTEGER         S, F, I, V   !  loop counters.
        INTEGER         SDEV            !  for surrogate coeff file
        INTEGER         XDEV            !  for surrogate xref  file
        INTEGER         LDEV            !  log-device
        INTEGER         SPCREF          !  scratch xref value
        INTEGER         IOS             !  I/O status
        INTEGER         IREC            !  input line (record) number
        INTEGER         I0              !  for counting FIP-independent srcs
        INTEGER         I1              !  for counting STATE-dependent srcs
        INTEGER         I2              !  for counting FIP-dependent   srcs
        REAL            FAC             !  common factor subexpression
        CHARACTER*16    ANAME   !  logical name for area-source  input file
        CHARACTER*16    MNAME   !  logical name for spec matrix output file
        CHARACTER*80    MSGBUF  !  message buffer for M3ERR()


C***********************************************************************
C   begin body of program SPCAMAT
        
        LDEV = INIT3()
        
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

        XDEV = PROMPTFFILE( 
     &           'Enter logical name for SPECIATION XREF file',
     &           .TRUE., .TRUE., 'ASREF', 'SPCAMAT' )

        SDEV = PROMPTFFILE( 
     &           'Enter logical name for SPECIATION PROFILES file',
     &           .TRUE., .TRUE., 'ASPRO', 'SPCAMAT' )

C.......   Read in description of area-source file
        
        IF ( .NOT. DESC3( ANAME ) ) THEN
            CALL M3ERR( 'SPCAMAT', 0, 0, 
     &                  'Could not get description of file ' // ANAME,
     &                  .TRUE. )
        END IF
        
        
C.......   Open speciation matrix file:  reuse much of AREA file description
        
        NVARS3D = NMPOL
        FTYPE3D = GRDDED3
        DO  44  V = 1, NMPOL

            I = INDEX( SFNAM( V ), '_to_' )
            VNAME3D( V ) = SFNAM( V )
            UNITS3D( V ) = 'moles/gm'
            VDESC3D( V ) = 'Conversion factor, ' // 
     &              SFNAM( V )( 1:I-1 ) // ' to ' // 
     &              SFNAM( V )( I+4:16 )
            VTYPE3D( V ) = M3REAL
44      CONTINUE

        FDESC3D( 1 ) = 
     &  'Matrix of speciation factors, by source, for converting '
        FDESC3D( 2 ) = 
     &  'inventory species into model species for the CB-IV chemical'
        FDESC3D( 3 ) = 
     &  'mechanism.'
        DO  45  V = 4, MXDESC3
            FDESC3D( V ) = ' '
45      CONTINUE

        MNAME = PROMPTMFILE( 
     &  'Enter logical name for SPECIATION MATRIX output file',
     &          FSUNKN3, 'ASMAT', 'SPCAMAT' )


C...........   Read in emissions indexing data IFIP, ASC7, ASC3.
        
        IF ( .NOT. READ3( ANAME, 'FIP', ALLAYS3, 0, 0, FIPS ) ) THEN
            CALL M3ERR( 'SPCAMAT', 0, 0, 
     &                  'Error reading FIP from file' // ANAME,
     &                  .TRUE. )
        END IF

        IF ( .NOT. READ3( ANAME, 'ASC7', ALLAYS3, 0, 0, ASC7 ) ) THEN
            CALL M3ERR( 'SPCAMAT', 0, 0, 
     &                  'Error reading ASC7 from file' // ANAME,
     &                  .TRUE. )
        END IF

        IF ( .NOT. READ3( ANAME, 'ASC3', ALLAYS3, 0, 0, ASC3 ) ) THEN
            CALL M3ERR( 'SPCAMAT', 0, 0, 
     &                  'Error reading ASC3 from file' // ANAME,
     &                  .TRUE. )
        END IF


C.......   Read the SPECIATION PROFILES FILE
C...............   ASSUMPTION:  sorted by FIP, ASC
        
        WRITE( *,92000 ) 
     &       ' ', 'Reading SPECIATION PROFILES file...', ' '

        IREC  =  0
        
101     CONTINUE        !  head of the SDEV-read loop
        
            READ( SDEV, *, END=122, IOSTAT=IOS ) 
     &          SID, ( FACS( V ), V = 4, NMPOL )

            IREC = IREC + 1
            IF ( IOS .NE. 0 ) THEN
                WRITE( *,92010 ) 
     &              'Unit # ', SDEV, 
     &              'IOSTAT ', IOS, 
     &              'Line # ', IREC, ' '
                CALL M3ERR( 'SPCAMAT', 0, 0, 
     &              'I/O error reading SPECIATION PROFILE file.',
     &              .TRUE. )
            ELSE IF ( IREC .GT. NSPRA ) THEN 
                WRITE( *,92010 ) 
     &              'Record number  ', IREC,
     &              'Dimensioned max', NSPRA, ' '
                CALL M3ERR( 'SPCAMAT', 0, 0, 
     &              'Max speciation profile table size NSPRA exceeded',
     &              .TRUE. )
            END IF
            
            SXREF( IREC ) = SID
            DO  111  V = 4, NMPOL
                SCOEF( IREC,V ) = FACS( V )
111         CONTINUE

            GO TO  101  !  to head of loop
      
122     CONTINUE        !  end of the SDEV-read loop

        NPROFS = IREC



C.......   Read the speciation cross-reference file
C...............   ASSUMPTION:  sorted by FIP and ASC (i.e., <ID7,ID3> pair)

        WRITE( *,92000 ) 
     &      ' ', 'Reading SPECIATION XREF file...', ' '

        IREC = 0
        I0   = 0
        I1   = 0
        I2   = 0

133     CONTINUE        !  head of the XDEV-read loop
        
            READ( XDEV, 93030, END=135, IOSTAT=IOS ) 
     &          FIP, ID7, ID3, SID, VREACT, ETHADJ, FRMADJ, ALLADJ
            IREC = IREC + 1

            IF ( IOS .NE. 0 ) THEN
            
                WRITE( MSGBUF,94050 ) 
     &              'Error ', IOS, 
     &              'reading SPECIATION XREF file at line', IREC
                CALL M3EXIT( 'SPCAMAT', 0, 0, MSGBUF, 2 )
            
            END IF
            
            IF ( FIP .EQ. 0  .AND.  I0 .LT. NASCC ) THEN !  FIP independent refs

                I0 = I0 + 1
                SPXRF0( I0 ) = SID 
                ASC7X0( I0 ) = ID7 
                ASC3X0( I0 ) = ID3
                VREAC0( I0 ) = VREACT
                ETHAD0( I0 ) = ETHADJ
                FRMAD0( I0 ) = FRMADJ
                ALLAD0( I0 ) = ALLADJ

                GO TO  133

            ELSE IF ( ( MOD( FIP, 1000 ) .EQ. 0 )
     &                .AND.     ( I1 .LT. NASCC ) ) THEN ! state-dependent refs

                I1 = I1 + 1
                FIPSXS( I1 ) = FIP
                SPXRFS( I1 ) = SID 
                ASC7XS( I1 ) = ID7 
                ASC3XS( I1 ) = ID3
                VREACS( I1 ) = VREACT
                ETHADS( I1 ) = ETHADJ
                FRMADS( I1 ) = FRMADJ
                ALLADS( I1 ) = ALLADJ

                GO TO  133

            ELSE IF ( FIP .GT. 0 .AND. I2 .LT. NASCC ) THEN  ! FIP-dep refs

                I2 = I2 + 1
                FIPSXF( I2 ) = FIP
                SPXRFF( I2 ) = SID 
                ASC7XF( I2 ) = ID7 
                ASC3XF( I2 ) = ID3
                VREACF( I2 ) = VREACT
                ETHADF( I2 ) = ETHADJ
                FRMADF( I2 ) = FRMADJ
                ALLADF( I2 ) = ALLADJ

                GO TO  133

            ELSE
            
                WRITE( MSGBUF,94050 ) 
     &          'Max table-size ', NASCC,  
     &          'exceeded in XREF file at line', IREC
                CALL M3EXIT( 'SPCAMAT', 0, 0, MSGBUF, 2 )
            
            END IF

135     CONTINUE        !  end of the SDEV-read loop
        
        NREF0 = I0              !  number of location-independent refs
        NREFS = I1              !  " state-only         dependent refs
        NREFF = I2              !  " state-and-county   dependent refs
        
        
C.......   Compute the speciation coefficient matrix
C...........   ASSUMPTION:  everything sorted by FIP and ASC, or by XREF ID
C...........   Two stage process:  first try for cross-referece by FIP and ASC
C...........   Then if that fails, try by ASC only.

        DO  144  S = 1, NASRC

            FIP = FIPS( S )
            ID7 = ASC7( S )
            ID3 = ASC3( S )

            F = FIND3( FIP, ID7, ID3, NREFF, FIPSXF, ASC7XF, ASC3XF )

            IF ( F .GT. 0 ) THEN        !  county-scc match:

                SPCREF = SPXRFF( F )
                ALLADJ = ALLADF( F ) * TON2GM
                FRMADJ = FRMADF( F )
                ETHADJ = ETHADF( F )
                VREACT = VREACF( F )

            ELSE        !  no county-scc match; try state-scc match

                SID = 1000 * ( FIP / 1000 )     !  state code
                F   = FIND3( SID, ID7, ID3, 
     &                       NREFS, FIPSXS, ASC7XS, ASC3XS )

                IF ( F .GT. 0 ) THEN

                    SPCREF = SPXRFS( F )
                    ALLADJ = ALLADS( F ) * TON2GM
                    FRMADJ = FRMADS( F )
                    ETHADJ = ETHADS( F )
                    VREACT = VREACS( F )

                ELSE            !  no state-scc match; try just scc

                    F = FIND2( ID7, ID3, NREF0, ASC7X0, ASC3X0 )

                    IF ( F .GT. 0 ) THEN

                        SPCREF = SPXRF0( F )
                        ALLADJ = ALLAD0( F ) * TON2GM
                        FRMADJ = FRMAD0( F )
                        ETHADJ = ETHAD0( F )
                        VREACT = VREAC0( F )

                    ELSE        !  not found:  error

                        WRITE( MSGBUF,93000 ) 
     &                      'Error computing speciation matrix:  '//
     &                      'source not found'
                        WRITE( *,92010 ) 
     &                      'Source #', S, 
     &                      'FIP code', FIP,
     &                      'ASC7    ', ID7,
     &                      'ASC3    ', ID3
                        CALL M3EXIT( 'SPCAMAT', 0, 0, MSGBUF, 2 )

                    END IF      !  if find2(scc) succeeded, or not

                END IF  !  if find3(state,scc) succeeded, or not

            END IF      !  if find3(county,scc) succeeded, or not

            I = FIND1( SPCREF, NPROFS, SXREF )

            IF ( I .LE. 0 ) THEN        !  error:  XREF index not found

                WRITE( MSGBUF,93000 ) 
     &              'Error computing spc matrix:   XREF not found.'
                WRITE( *,92010 ) 
     &              'FIPS code', FIP,
     &              'ASC7 code', ID7,
     &              'ASC3 code', ID3, 
     &              'XREF code', SPXRF0( F )
                CALL M3EXIT( 'SPCAMAT', 0,0, MSGBUF, 2 )

            END IF      !  if i nonpositive (cross-reference not found)

            SPCFACS( S, 1 ) = ALLADJ * COFAC      ! 1/molwt for CO
            SPCFACS( S, 2 ) = ALLADJ * NOFAC      ! 90% of 1/molwt for NO
            SPCFACS( S, 3 ) = ALLADJ * NO2FAC     ! 10% of 1/molwt for NO2
            FAC = ALLADJ	!  if THC used, need VREACT * ALLADJ
            DO  143  V = 4, NMPOL
                SPCFACS( S,V ) = FAC * SCOEF( I,V )
143         CONTINUE
            SPCFACS( S, 8 ) = FRMADJ * SPCFACS( S, 8 )
            SPCFACS( S,10 ) = ETHADJ * SPCFACS( S,10 )

144     CONTINUE        ! end loop on area sources, computing speciation matrix

        
C.......   Write out the speciation matrix:
        
        WRITE( *,92000 ) ' ', 'Writing out SPECIATION MATRIX...', ' '

        IF ( .NOT. WRITE3( MNAME, 'ALL', 0, 0, SPCFACS ) ) THEN
            CALL M3EXIT( 'SPCAMAT', 0, 0, 
     &                   'Error writing SPECIATION MATRIX file.', 2 )
        END IF



999   CONTINUE          !  exit program
      
      CALL M3EXIT( 'SPCAMAT', 0, 0, 'Normal completion.', 0 )


C******************  FORMAT  STATEMENTS   ******************************

C...........   Informational (LOG) message formats... 92xxx

92000   FORMAT( 5X, A )

92010   FORMAT ( 5X , A, :, I10 )


C...........   Formatted file I/O formats............ 93xxx

93000   FORMAT( A )

93010   FORMAT( A16 )

93020   FORMAT( I7, I3 )

93030   FORMAT( I5, I8, I3, I6, F6.3, 15X, 3 F6.3 ) !  for SPC XREF file


C...........   Internal buffering formats............ 94xxx

94010   FORMAT( I5 )

94020   FORMAT( I7 )

94030   FORMAT( I3 )

94040   FORMAT( G10.3 )

94050   FORMAT( 3( A, :, I10, :, 2X ) )


        END

