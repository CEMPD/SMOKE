
C Version "@(#)$Id$ $Source$ $Date$ 

        PROGRAM SPCMMAT

C***********************************************************************
C  program body starts at line 180
C
C  DESCRIPTION:
C       Construct mobile source chemical speciation matrix from data 
C       contained in EPS-style chem-split and chem-xref files.
C       Uses fixed NOtoNOX  factor of 0.0197  moles  NO/gm NOX (from EPS)
C       Uses fixed NO2toNOX factor of 0.00217 moles NO2/gm NOX  "
C
C  PRECONDITIONS REQUIRED:
C       Sorted, cut-down input data for chem-split and chem-xref files;
C       M3IO mobile-sources file from program GRDAMAT
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C       Models-3 I/O
C       FIND1, FIND2, FIND3, GETEFILE, GETYN, TRIMLEN
C
C  REVISION  HISTORY:
C       Prototype  1/96 by CJC.
C
C***********************************************************************

      IMPLICIT NONE

C...........   INCLUDES:
        
        INCLUDE 'CHDIMS3.EXT'   !  emis chem parms (inventory + model)
        INCLUDE 'MBDIMS3.EXT'   !  mobile-source dimensioning parameters
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
        INCLUDE 'FDESC3.EXT'    !  I/O API file description data structures.


C...........   EXTERNAL FUNCTIONS and their descriptions:
        
        INTEGER         FIND1, FIND2, FIND3, FIND4
        INTEGER         GETEFILE
        LOGICAL         GETYN
        INTEGER         PROMPTFFILE
        CHARACTER*16    PROMPTMFILE
        INTEGER         RDSPROF
        INTEGER         TRIMLEN
        
        EXTERNAL        FIND1, FIND2, FIND3, FIND4, GETEFILE, GETYN, 
     &                  PROMPTFFILE, PROMPTMFILE, RDSPROF, TRIMLEN


C...........   PARAMETERS and their descriptions:

        REAL        TON2GM
        PARAMETER ( TON2GM = 907184.74 )
        
C...........   LOCAL VARIABLES and their descriptions:
C...........   NOTE that ASC (mobile-source-category) ID's are 10-digit
C...........   unsigned integers which may be treated as a leading 7-digit
C...........   field, and a trailing 3-digit field.  *7 and *3 arrays below
C...........   follow this scheme with parallel arrays

        INTEGER     FIP
        INTEGER     SID
        INTEGER     RDT
        INTEGER     EMT
        INTEGER     VHT

C...........   Mobile Source attributes from ANAME                                 

        INTEGER         IFIP  ( NMSRC )        ! FIP codes for sources
        INTEGER         IRCLAS( NMSRC )        ! EPA road classification code
        INTEGER         ILINK ( NMSRC )        ! link number or 0 for non-link

C...........   Speciation profiles table
        
        REAL        VREACT
        REAL        ETHADJ
        REAL        FRMADJ
        REAL        ALLADJ
        
        INTEGER     NPROFS                      !  actual number of profiles
        INTEGER     SPROF( NSPRM )              !  profile indexes (sorted)
        REAL        SCOEF( NSPRM, NMPOL )       !  profile coefficients

C...........   Speciation cross-reference table:  
C.......   Before sorting:

        INTEGER     NREFA
        INTEGER     INDEXA( MXREF )
        INTEGER     FIPSXA( MXREF )
        INTEGER     RDTSXA( MXREF )
        INTEGER     EMTSXA( MXREF )
        INTEGER     VHTSXA( MXREF )
        INTEGER     SPPROA( MXREF )
        REAL        VREACA( MXREF )
        REAL        ETHADA( MXREF )
        REAL        FRMADA( MXREF )
        REAL        ALLADA( MXREF )

C.......   After sorting -- four parts:   first part independent of road class
            
        INTEGER     NREFV               !  number of actual entries
        INTEGER     EMTSXV( MXVREF )     !
        INTEGER     VHTSXV( MXVREF )     !
        INTEGER     SPPROV( MXVREF )    !  XREF for lookup in ASPRO
        REAL        VREACV( MXVREF )    !  
        REAL        ETHADV( MXVREF )    !  
        REAL        FRMADV( MXVREF )    !  
        REAL        ALLADV( MXVREF )    !  

C.......   Second part depends on road class (and vehicle type) only

        INTEGER     NREFR                !  number of actual entries
        INTEGER     RDTSXR( MXRREF )     !
        INTEGER     EMTSXR( MXRREF )     !
        INTEGER     VHTSXR( MXRREF )     !
        INTEGER     SPPROR( MXRREF )     !
        REAL        VREACR( MXRREF )     !  
        REAL        ETHADR( MXRREF )     !  
        REAL        FRMADR( MXRREF )     !  
        REAL        ALLADR( MXRREF )     !  

C.......   Third part depends on state, road class, and vehicle type

        INTEGER     NREFS
        INTEGER     FIPSXS( MXSREF )
        INTEGER     RDTSXS( MXSREF )
        INTEGER     EMTSXS( MXSREF )
        INTEGER     VHTSXS( MXSREF )
        INTEGER     SPPROS( MXSREF )
        REAL        VREACS( MXSREF )
        REAL        ETHADS( MXSREF )
        REAL        FRMADS( MXSREF )
        REAL        ALLADS( MXSREF )

C.......   Fourth part depends on FIP, road class, and vehicle type

        INTEGER     NREFF
        INTEGER     FIPSXF( MXFREF )
        INTEGER     RDTSXF( MXFREF )
        INTEGER     EMTSXF( MXFREF )
        INTEGER     VHTSXF( MXFREF )
        INTEGER     SPPROF( MXFREF )
        REAL        VREACF( MXFREF )
        REAL        ETHADF( MXFREF )
        REAL        FRMADF( MXFREF )
        REAL        ALLADF( MXFREF )

C.......   Speciation matrix variables and their names:
        
        REAL          SPCFACS( NMSRC, NVTYPE, NSMATX ) !  speciation coefficients


C...........   Other local variables
        
        INTEGER         S, E, F, I, J, K, V   !  loop counters.
        INTEGER         SDEV            !  for surrogate coeff file
        INTEGER         XDEV            !  for surrogate xref  file
        INTEGER         LDEV            !  log-device
        INTEGER         SPCPRO          !  scratch xref value
        INTEGER         IOS             !  I/O status
        INTEGER         IREC            !  input line (record) number
        INTEGER         I0              !  for counting RDT-independent refs
        INTEGER         I1              !  for counting RDT-dependent   refs
        INTEGER         I2              !  for counting state-dependent refs
        INTEGER         I3              !  for counting FIP-dependent   refs
        REAL            FAC             !  common factor subexpression
        CHARACTER*16    ANAME   !  logical name for mobile-source  input file
        CHARACTER*16    MNAME   !  logical name for spec matrix output file
        CHARACTER*256   MESG    !  message buffer for M3MESG() and M3EXIT()
        
        LOGICAL         EFLAG   !  input error flat


C***********************************************************************
C   begin body of program SPCMMAT
        
        LDEV = INIT3()
        
        WRITE( *,92000 ) 
     &  ' ',
     &  'Program SPCMMAT to take the netCDF mobile SOURCES file, the',
     &  'SPECIATION PROFILES file, and the SPECIATION CROSS-REFERENCE',
     &  'file, and produce the SPECIATION MATRIX file.',
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
            CALL M3EXIT( 'SPCMMAT', 0, 0, 'Ending program SPCMMAT', 2 )
        END IF


C.......   Get file name; open input mobile sources file

        ANAME = PROMPTMFILE( 
     &          'Enter logical name for MOBILE SOURCE input file',
     &          FSREAD3, 'MOBL', 'SPCMMAT' )

        XDEV = PROMPTFFILE( 
     &           'Enter logical name for SPECIATION XREF file',
     &           .TRUE., .TRUE., 'MSREF', 'SPCMMAT' )

        SDEV = PROMPTFFILE( 
     &           'Enter logical name for SPECIATION PROFILES file',
     &           .TRUE., .TRUE., 'MSPRO', 'SPCMMAT' )

C.......   Read in description of mobile-source file
        
        IF ( .NOT. DESC3( ANAME ) ) THEN
            CALL M3EXIT( 'SPCMMAT', 0, 0, 
     &                   'Could not get description of file ' // ANAME, 
     &                   2 )
        END IF
        
        
C.......   Open speciation matrix file:  reuse much of MOBILE file description
        
        NLAYS3D = NVTYPE
        NVARS3D = NSMATX
        FTYPE3D = GRDDED3
            
C.......   Create table of variable names.  First set of names
C.......   are CO_to_CO, ..., VOC_to_ISOP (from CHDIMS3.EXT's SFNAM)

        J = 0
        DO  44  V = 1, NMPOL

            J = J + 1
            VNAME3D( J ) = SFNAM( V )
            UNITS3D( J ) = 'moles/gm'
            VDESC3D( J ) = 'Conversion factor, ' // 
     &              EMTNAM( V )( 1:TRIMLEN( EMTNAM( V ) ) )// 
     &              ' to ' // EMNAM( V )
            VTYPE3D( J ) = M3REAL
44      CONTINUE

C.......   Remaining names are of the form [EVP|RNL|RST]_to_[CO|...|ISOP]
C.......   (non-running emissions-type _to_ VOC-derived model species)

        DO  55  K = NIPOL + 1, NETYPE
        DO  54  V = 4, NMPOL
            J = J + 1
            VNAME3D( J ) = EMTNAM( K )( 1:3 ) // '_to_' // EMNAM( V )
            UNITS3D( J ) = 'moles/gm'
            VDESC3D( J ) = 'Conversion factor, ' //
     &              EMTNAM( K )( 1:TRIMLEN( EMTNAM( K ) ) ) // 
     &              ' to ' // EMNAM( V )
            VTYPE3D( J ) = M3REAL
54      CONTINUE
55      CONTINUE

        FDESC3D( 1 ) = 
     &  'Matrix of speciation factors, by source, for converting '
        FDESC3D( 2 ) = 
     &  'mobile source emissions type species into model species for'
        FDESC3D( 3 ) = 
     &  'the CB-IV chemical mechanism.'
        DO  66  V = 4, MXDESC3
            FDESC3D( V ) = ' '
66      CONTINUE

        MNAME = PROMPTMFILE( 
     &  'Enter logical name for SPECIATION MATRIX output file',
     &          FSUNKN3, 'MSMAT', 'SPCMMAT' )


C...........   Read in emissions indexing data IFIP, IRCLAS, ILINK
        
        IF ( .NOT. READ3( ANAME, 'FIP', ALLAYS3, 0, 0, IFIP ) ) THEN
            CALL M3EXIT( 'SPCMMAT', 0, 0, 
     &                   'Error reading IFIP from file' // ANAME, 2 )
        END IF

        IF ( .NOT. READ3( ANAME, 'IRCLAS', ALLAYS3,0,0, IRCLAS ) ) THEN
            CALL M3EXIT( 'SPCMMAT', 0, 0, 
     &                   'Error reading IRCLAS from file' // ANAME, 2 )
        END IF

        IF ( .NOT. READ3( ANAME, 'ILINK', ALLAYS3,0,0, ILINK ) ) THEN
            CALL M3EXIT( 'SPCMMAT', 0, 0, 
     &                   'Error reading ILINK from file' // ANAME, 2 )
        END IF


C.......   Read/sort the SPECIATION PROFILES FILE
C...............   ASSUMPTION:  sorted by profile ID

        WRITE( *,92000 ) 
     &       ' ', 'Reading SPECIATION PROFILES file...', ' '

        
        NPROFS = RDSPROF( SDEV, NSPRM, SPROF, SCOEF )

        IF ( NPROFS .LE. 0 ) THEN
            CALL M3EXIT( 'SPCMMAT', 0, 0, 
     &                   'Error reading SPECIATION PROFILE file.', 2 )
        END IF



C.......   Read the speciation cross-reference file

        WRITE( *,92000 ) 
     &      ' ', 'Processing SPECIATION XREF file...', ' '

        IREC = 0
        I     =  0
        EFLAG = .FALSE.

101     CONTINUE        !  head of the XDEV-read loop

            IREC = IREC + 1
            READ( XDEV, *, END=111, IOSTAT=IOS ) 
     &          FIP, RDT, EMT, VHT, SID, VREACT, ETHADJ, FRMADJ, ALLADJ

            IF ( IOS .NE. 0 ) THEN

                EFLAG = .TRUE.
                WRITE( MESG,94050 ) 
     &              'Error ', IOS, 
     &              'reading SPECIATION XREF file at line', IREC
                CALL M3MESG( MESG )

            ELSE IF ( I .GE. MXREF ) THEN

                EFLAG = .TRUE.
                WRITE( MESG,94050 ) 
     &              'Max table-size ', MXREF,  
     &              'exceeded reading XREF file at line', IREC
                CALL M3MESG( MESG )

            ELSE

                I = I + 1
                INDEXA( I ) = I
                FIPSXA( I ) = FIP
                RDTSXA( I ) = RDT
                EMTSXA( I ) = EMT
                VHTSXA( I ) = VHT
                SPPROA( I ) = SID
                VREACA( I ) = VREACT
                ETHADA( I ) = ETHADJ
                FRMADA( I ) = FRMADJ
                ALLADA( I ) = ALLADJ

            END IF
            
            GO TO  101  !  to head of the XDEV-read loop

111     CONTINUE        !  end of the XDEV-read loop

        IF ( EFLAG ) THEN
            CALL M3EXIT( 'SPCMMAT', 0, 0, 
     &                   'Error reading ACTUAL SCC file.', 2 )
        END IF

        NREFA = I


C...........   Now sort the XREF inputs and place into appropriate tables:

        CALL SORTI4( NREFA, INDEXA, FIPSXA, RDTSXA, EMTSXA, VHTSXA )

        I0    = 0
        I1    = 0
        I2    = 0
        I3    = 0
        EFLAG = .FALSE.

        DO  122  I = 1, NREFA
        
            J   = INDEXA( I )
            FIP = FIPSXA( J ) 
            RDT = RDTSXA( J )
            EMT = EMTSXA( J )
            VHT = VHTSXA( J )
            
            IF ( RDT .EQ. 0  ) THEN !  FIP independent refs

                IF ( I0 .LT. MXVREF ) THEN

                    I0 = I0 + 1
                    SPPROV( I0 ) = SPPROA( J )
                    EMTSXV( I0 ) = EMT
                    VHTSXV( I0 ) = VHT
                    VREACV( I0 ) = VREACA( J )
                    ETHADV( I0 ) = ETHADA( J )
                    FRMADV( I0 ) = FRMADA( J )
                    ALLADV( I0 ) = ALLADA( J )

                ELSE

                    EFLAG = .TRUE.
                    WRITE( MESG,94050 )
     &                  'Max RDT-independent table-size ', NVTYPE,  
     &                  'exceeded in XREF file at line', IREC
                    CALL M3MESG( MESG )

                END IF

            ELSE IF ( FIP .EQ. 0  ) THEN !  FIP independent refs

                IF ( I1 .LT. MXRREF ) THEN

                    I1 = I1 + 1
                    SPPROR( I1 ) = SPPROA( J )
                    RDTSXR( I1 ) = RDT
                    EMTSXR( I1 ) = EMT
                    VHTSXR( I1 ) = VHT
                    VREACR( I1 ) = VREACA( J )
                    ETHADR( I1 ) = ETHADA( J )
                    FRMADR( I1 ) = FRMADA( J )
                    ALLADR( I1 ) = ALLADA( J )

                ELSE

                    EFLAG = .TRUE.
                    WRITE( MESG,94050 )
     &                  'Max FIP-independent table-size ', MXRREF,  
     &                  'exceeded in XREF file at line', IREC
                    CALL M3MESG( MESG )

                END IF

            ELSE IF ( MOD( FIP, 1000 ) .EQ. 0 )  THEN ! state-dependent refs

                IF ( I2 .LT. MXSREF ) THEN

                    I2 = I2 + 1
                    FIPSXS( I2 ) = FIP/1000
                    SPPROS( I2 ) = SPPROA( J )
                    RDTSXS( I2 ) = RDT
                    EMTSXS( I2 ) = EMT
                    VHTSXS( I2 ) = VHT
                    VREACS( I2 ) = VREACA( J )
                    ETHADS( I2 ) = ETHADA( J )
                    FRMADS( I2 ) = FRMADA( J )
                    ALLADS( I2 ) = ALLADA( J )

                ELSE

                    EFLAG = .TRUE.
                    WRITE( MESG,94050 )
     &                  'Max state-dependent table-size ', MXSREF,  
     &                  'exceeded in XREF file at line', IREC
                    CALL M3MESG( MESG )

                END IF

            ELSE                ! fip > 0:  FIP-dep refs

                IF ( I3 .LT. MXFREF ) THEN

                    I3 = I3 + 1
                    FIPSXF( I3 ) = FIP
                    SPPROF( I3 ) = SPPROA( J )
                    RDTSXF( I3 ) = RDT
                    EMTSXF( I3 ) = EMT
                    VHTSXF( I3 ) = VHT
                    VREACF( I3 ) = VREACA( J )
                    ETHADF( I3 ) = ETHADA( J )
                    FRMADF( I3 ) = FRMADA( J )
                    ALLADF( I3 ) = ALLADA( J )

                ELSE

                    EFLAG = .TRUE.
                    WRITE( MESG,94050 )
     &                  'Max FIP-dependent table-size ', MXFREF,  
     &                  'exceeded in XREF file at line', IREC
                    CALL M3MESG( MESG )

                END IF

            END IF

122     CONTINUE        !  end of the SDEV-read loop
        
        IF ( EFLAG ) THEN
                CALL M3EXIT( 'SPCMMAT', 0, 0, 
     &                       'Max table-size exceeded in XREF file', 2 )
        END IF
        
        NREFV = I0              !  number of location-independent refs
        NREFR = I1              !  " road class-only    dependent refs
        NREFS = I2              !  " state-only         dependent refs
        NREFF = I3              !  " state-and-county   dependent refs
        EFLAG = .FALSE.
        
        
C.......   Compute the speciation coefficient matrix
C...........   ASSUMPTION:  everything sorted by FIP and ASC, or by XREF ID
C...........   Two stage process:  first try for cross-referece by FIP and ASC
C...........   Then if that fails, try by ASC only.

        K = 0                   !  "vars" subscript into SPCFACS

        DO  146  EMT = 1, NEMODE
        DO  145  VHT = 1, NVTYPE
        DO  144  S   = 1, NMSRC

            FIP = IFIP  ( S )
            RDT = IRCLAS( S )

            F = FIND4( FIP, RDT, EMT, VHT, 
     &                 NREFF, FIPSXF, RDTSXF, EMTSXF, VHTSXF )

            IF ( F .GT. 0 ) THEN        !  full match:

                SPCPRO = SPPROF( F )
                ALLADJ = ALLADF( F ) * TON2GM
                FRMADJ = FRMADF( F )
                ETHADJ = ETHADF( F )
                VREACT = VREACF( F )
                GO TO  131

            END IF              !  if full match

            SID = 1000 * ( FIP / 1000 )     !  state code
            F   = FIND4( SID, RDT, EMT, VHT, 
     &                   NREFS, FIPSXS, RDTSXS, EMTSXS, VHTSXS )
           
            IF ( F .GT. 0 ) THEN        !  state match
           
                SPCPRO = SPPROS( F )
                ALLADJ = ALLADS( F ) * TON2GM
                FRMADJ = FRMADS( F )
                ETHADJ = ETHADS( F )
                VREACT = VREACS( F )
                GO TO  131

            END IF              !  if state match

     
            F = FIND3( RDT, EMT, VHT, NREFR, RDTSXR, EMTSXR, VHTSXR )

            IF ( F .GT. 0 ) THEN        !  roadclass match

                SPCPRO = SPPROR( F )
                ALLADJ = ALLADR( F ) * TON2GM
                FRMADJ = FRMADR( F )
                ETHADJ = ETHADR( F )
                VREACT = VREACR( F )
                GO TO  131

            END IF              !  if roadclass matcjh

            F = FIND2( EMT, VHT, NREFV, EMTSXV, VHTSXV )

            IF ( F .GT. 0 ) THEN        !  vehicletype/emissionmode-only match

                SPCPRO = SPPROV( F )
                ALLADJ = ALLADV( F ) * TON2GM
                FRMADJ = FRMADV( F )
                ETHADJ = ETHADV( F )
                VREACT = VREACV( F )
                GO TO  131

            END IF              !  if roadclass matcjh

            EFLAG = .TRUE.      !  if you get to here:  no match found
            WRITE( MESG,94010 ) 
     &          'Source', S, 
     &          '(FIP', FIP,
     &          'RDT    ', RDT,
     &          'VHT    ', V, ') not found'
            CALL M3MESG( MESG )
      
            GO TO 144      !  to end of loop
      

131         CONTINUE            !  Cross-reference lookups succeeded

            I = FIND1( SPCPRO, NPROFS, SPROF )

            IF ( I .LE. 0 ) THEN        !  error:  XREF index not found

                EFLAG = .TRUE.
                WRITE( MESG,94010 ) 
     &              'XREF for FIPS', FIP,
     &              'RDT', RDT,
     &              'VHT', V, 
     &              'profile code', SPCPRO, 'not found'
                CALL M3MESG( MESG )

                GO TO  144      !  to end of loop

            END IF      !  if i nonpositive (cross-reference not found)

C...........   If E=1 (EXH), all spcs; otherwise, only VOC-related spcs

            IF ( E .EQ. 1 ) THEN                        !  EXH
                SPCFACS( S,V,1 ) = ALLADJ * COFAC       ! 1/molwt for CO
                SPCFACS( S,V,2 ) = ALLADJ * NOFAC       ! 90% of 1/molwt for NO
                SPCFACS( S,V,3 ) = ALLADJ * NO2FAC      ! 10% of 1/molwt for NO2
            END IF
            
            FAC = ALLADJ            !  if THC used, need VREACT * ALLADJ

            DO  141  I = 4, NMPOL
                K = K + 1
                IF ( I .EQ. 8 ) THEN
                    SPCFACS( S,V,K ) = FRMADJ * FAC * SCOEF( I,V )
                ELSE IF ( I .EQ. 10 ) THEN
                    SPCFACS( S,V,K ) = ETHADJ * FAC * SCOEF( I,V )
                ELSE
                    SPCFACS( S,V,K ) = FAC * SCOEF( I,V )
                END IF
141         CONTINUE      ! end loop on VOC-related species

144     CONTINUE      ! end loop on mobile sources
145     CONTINUE      ! end loop on vehicle types
146     CONTINUE      ! end loop on emissions classes, computing spc matrix

        IF ( EFLAG ) THEN
            CALL M3EXIT( 'SPCMMAT', 0, 0, 
     &                   'Error computing SPECIATION MATRIX', 2 )
        END IF

        
C.......   Write out the speciation matrix:
        
        WRITE( *,92000 ) ' ', 'Writing out SPECIATION MATRIX...', ' '

        IF ( .NOT. WRITE3( MNAME, 'ALL', 0, 0, SPCFACS ) ) THEN
            CALL M3EXIT( 'SPCMMAT', 0, 0, 
     &                   'Error writing SPECIATION MATRIX file.', 2 )
        END IF



999   CONTINUE          !  exit program
      
      CALL M3EXIT( 'SPCMMAT', 0, 0, 'Normal completion.', 0 )


C******************  FORMAT  STATEMENTS   ******************************

C...........   Informational (LOG) message formats... 92xxx

92000   FORMAT( 5X, A )

92010   FORMAT ( 5X , A, :, I10 )


C...........   Formatted file I/O formats............ 93xxx

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( I5 )

94020   FORMAT( I7 )

94030   FORMAT( I3 )

94040   FORMAT( G10.3 )

94050   FORMAT( 10 ( A, :, I10, :, 2X ) )


        END

