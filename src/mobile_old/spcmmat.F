
C Version "@(#)$Id$ $Source$ $Date$ 

        PROGRAM SPCMMAT

C***********************************************************************
C  program body starts at line 237
C
C  DESCRIPTION:
C       Construct mobile source chemical speciation matrix from data 
C       contained in EPS-style chem-split and chem-xref files.
C       Uses fixed NO to NOX  factor of 0.0197  moles  NO/gm NOX (from EPS)
C       Uses fixed NO2 to NOX factor of 0.00217 moles NO2/gm NOX  "
C
C  PRECONDITIONS REQUIRED:
C       M3IO mobile-sources file from program GRDAMAT
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C       Models-3 I/O
C       CSVNAME, FIND1, FIND2, FIND3, GETEFILE, GETYN, INDEX1, 
C       RDSPROF, TRIMLEN
C
C  REVISION  HISTORY:
C       Prototype  1/96 by CJC.
C       Revisions 5/96 by MRH.
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

        INTEGER         CSVNAME
        INTEGER         FIND1, FIND2, FIND3, FIND4
        INTEGER         GETEFILE
        INTEGER         GETMENU
        LOGICAL         GETYN
        INTEGER         INDEX1
        INTEGER         PROMPTFFILE
        CHARACTER*16    PROMPTMFILE
        INTEGER         RDSPROF
        INTEGER         STR2INT
        REAL            STR2REAL
        INTEGER         TRIMLEN
        
        EXTERNAL        FIND1, FIND2, FIND3, FIND4, GETEFILE, GETMENU,
     &                  GETYN, INDEX1, PROMPTFFILE, PROMPTMFILE, 
     &                  RDSPROF, STR2INT, STR2REAL, TRIMLEN


C...........   PARAMETERS and their descriptions:

        
C...........   LOCAL VARIABLES and their descriptions:
C...........   NOTE that ASC (mobile-source-category) ID's are 10-digit
C...........   unsigned integers which may be treated as a leading 7-digit
C...........   field, and a trailing 3-digit field.  *7 and *3 arrays below
C...........   follow this scheme with parallel arrays

        INTEGER     FIP    ! State and County code
        INTEGER     SID    ! State code only
        INTEGER     RDT    ! Road class number
        INTEGER     EMT    ! Emission type number
        INTEGER     VHT    ! Vehicle type number

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
        INTEGER     SXREF( NSPRM )              !  profile indexes (sorted)
        REAL        SCOEF( NSPRM, MXMPOL )      !  profile coefficients

C...........   Speciation naming and indexing lists

        INTEGER     IDXSPC  ( MXEPRC, MXMPOL )!  Pointer to actual model species
        INTEGER     IDXVAR  ( MXEPRC, MXMPOL )!  Pointer to SMAT variables
        INTEGER     NEMODE                    !  No. MV processes (emis modes)
        INTEGER     NIPOL                     !  No. inventory pollutants
        INTEGER     NMPOL                     !  No. model species
        INTEGER     NSNAME                    !  No. of variables in SMAT
        INTEGER     NSP     ( MXEPRC )        !  No. model spcs per MV process
        CHARACTER*8 EINAM   ( MXIPOL )        !  Inventory species names
        CHARACTER*8 EMNAM   ( MXMPOL )        !  Model species names
        CHARACTER*8 ISPC_PRC( MXEPRC, MXMPOL )!  Inv species per MV process
        CHARACTER*8 MSPC_PRC( MXEPRC, MXMPOL )!  Model species per MV process
        CHARACTER*8 MEMODE  ( MXEPRC )        !  Emission modes
 
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
        INTEGER     EMTSXV( MXVREF )    !
        INTEGER     VHTSXV( MXVREF )    !
        INTEGER     SPPROV( MXVREF )    !  XREF for lookup in MSPRO
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

C...........   Volatile Pollutant Type Conversion Table
C........ After sorting, two parts
 
C........ Default
 
        REAL        CONVV0( NVTYPE, MXEPRC )  !  conversion factor
 
C........ First part depends on FIP,road class
        INTEGER     NCNV1
        INTEGER     FIPSV1( MXSCNV )  !  conversion factor FIP
        INTEGER     RDCLV1( MXSCNV )  !  conversion factor road class
        REAL        CONVV1( MXSCNV, NVTYPE, MXEPRC )  !  conversion factor

C........ Second part depends on State,road class
        INTEGER     NCNV2
        INTEGER     FIPSV2( MXSCNV )  !  conversion factor FIP
        INTEGER     RDCLV2( MXSCNV )  !  conversion factor road class
        REAL        CONVV2( MXSCNV, NVTYPE, MXEPRC )  !  conversion factor

C........ Third part depends on road class only
        INTEGER     NCNV3
        INTEGER     RDCLV3( NRCLAS )  !  conversion factor road class
        REAL        CONVV3( NRCLAS, NVTYPE, MXEPRC )  !  conversion factor

C.......   Speciation matrix variables and their names:
        
        REAL          SPCFACS( NMSRC, NVTYPE, NSMATX ) ! speciation coefficients

C...........   Other local variables
        
        INTEGER         S, F, I, J, K, L, P, V !  index counters
        INTEGER         CHOICE          !  selection from menu 
        INTEGER         SDEV            !  unit for surrogate coeff file
        INTEGER         XDEV            !  unit for surrogate xref  file
        INTEGER         KDEV            !  unit for conversion file
        INTEGER         LDEV            !  log-device
        INTEGER         TDEV            !  for list file
        INTEGER         SPCPRO          !  scratch xref value
        INTEGER         IOS             !  I/O status
        INTEGER         IREC            !  input line (record) number
        INTEGER         I0              !  for counting RDT-independent refs
        INTEGER         I1              !  for counting RDT-dependent   refs
        INTEGER         I2              !  for counting state-dependent refs
        INTEGER         I3              !  for counting FIP-dependent   refs
        INTEGER         NRDT            !  index for RDT in MROADS3

        REAL            FAC             !  common factor subexpression
        REAL            VOLADJ          !  common volatile adjustment 

        CHARACTER*3     INV     !  inventory species name
        CHARACTER*3     PROC    !  MV process name
        CHARACTER*5     INVPRO  !  Inventory pollutant from speciation profile
        CHARACTER*5     POLV    !  volatile inventory pollutant in inventory
        CHARACTER*8     MODEL   !  model species name
        CHARACTER*16    ANAME   !  logical name for mobile-source  input file
        CHARACTER*16    MNAME   !  logical name for spec matrix output file
        CHARACTER*256   MESG    !  message buffer for M3MESG() and M3EXIT()
        CHARACTER*80    LINE    !  input line from XREF file
        
        LOGICAL         EFLAG   !  input error flat
        LOGICAL         KFLAG   !  if volatile format conversion factor or not 


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

        TDEV = PROMPTFFILE(
     &           'Enter logical name for SPECIATION LIST file',
     &           .TRUE., .TRUE., 'MSLST', 'SPCMMAT' )

C.......   Read in description of mobile-source file
        
        IF ( .NOT. DESC3( ANAME ) ) THEN
            CALL M3EXIT( 'SPCMMAT', 0, 0, 
     &                   'Could not get description of file ' // ANAME, 
     &                   2 )
        END IF
     
C.......   Read speciation list file to find out what species and emission
C.......   processes are actually being used in this scenario

        NIPOL  = CSVNAME( TDEV, .TRUE. , .FALSE., .FALSE., EINAM  )

C.......   For first time CSVNAME is called, see if there were any errors
        IF( NIPOL .LE. 0 ) THEN
            WRITE( MESG,94050 ) 
     &             'Error reading SPECIATION LIST file'
            CALL M3EXIT( 'SPCMMAT', 0, 0, MESG, 2 )

        ENDIF

        NMPOL  = CSVNAME( TDEV, .FALSE., .TRUE. , .FALSE., EMNAM  )
        NEMODE = CSVNAME( TDEV, .FALSE., .FALSE., .TRUE. , MEMODE )

C.......   Get speciation factor names and the number of names 
C.......   from the SPECIATION LIST file

        NSNAME = CSVNAME( TDEV, .TRUE., .TRUE., .TRUE., VNAME3D )

C.......   Intialize species per MV process counter
        DO 44 I = 1, NEMODE
            NSP( I ) = 0
44      CONTINUE

C.......   Process the SMAT variable names
        DO 55 I = 1, NSNAME

            L   = TRIMLEN( VNAME3D( I ) )
            J   = INDEX( VNAME3D( I )(   1:L ), '_' )
            K   = INDEX( VNAME3D( I )( J+1:L ), '_' ) + J

            PROC = VNAME3D( I )(  1 : J-1)
            INV  = VNAME3D( I )(J+1 : K-1)
            MODEL= VNAME3D( I )(K+1 : L  )

C.............  Set other variable properties
            UNITS3D( I ) = 'moles/gm'
            VDESC3D( I ) = 'Conversion factor for ' // 
     &         INV  ( 1:TRIMLEN( INV  ) ) // ' from process ' // 
     &         PROC ( 1:TRIMLEN( PROC ) ) // ' to species '   //
     &         MODEL( 1:TRIMLEN( MODEL) )
            VTYPE3D( I ) = M3REAL

C.............  Create MV process to model species indicies
            P = INDEX1( PROC, NEMODE, MEMODE )
            NSP( P ) = NSP( P ) + 1
            MSPC_PRC( P, NSP( P ) ) = MODEL
            ISPC_PRC( P, NSP( P ) ) = INV
            IDXVAR  ( P, NSP( P ) ) = I
            IDXSPC  ( P, NSP( P ) ) = INDEX1( MODEL, NMPOL, EMNAM )
            
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

C.......   Open speciation matrix file:  reuse much of MOBILE file description
        
        NLAYS3D = NVTYPE
        NVARS3D = NSNAME
        FTYPE3D = GRDDED3
            
        MNAME = PROMPTMFILE( 
     &  'Enter logical name for SPECIATION MATRIX output file',
     &          FSUNKN3, 'MSMAT', 'SPCMMAT' )

C...........   Read in emissions indexing data IFIP, IRCLAS, ILINK
        
        IF ( .NOT. READ3( ANAME, 'IFIP', ALLAYS3, 0, 0, IFIP ) ) THEN
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
C...............   ASSUMPTIONS: CSVNAME has provided model species list (EMNAM)

        WRITE( *,92000 ) 
     &       ' ', 'Reading SPECIATION PROFILES file...', ' '

        
        NPROFS = RDSPROF( SDEV , NMPOL, EMNAM, NSPRM, 
     &                    SXREF, SCOEF, INVPRO )

        IF ( NPROFS .LE. 0 ) THEN
            CALL M3EXIT( 'SPCMMAT', 0, 0, 
     &                   'Error reading SPECIATION PROFILE file.', 2 )
        END IF

C.......   Compare inventory pollutant in the speciation profile with
C.......   volatile pollutant in inventory.  If not the same, prompt for
C.......   a conversion file.
 
C.......   Prompt for which inventory pollutant is volatile
 
        CHOICE = GETMENU( NIPOL, 1,
     &                    'Which inventory pollutant is volatile?',
     &                    EINAM )
 
        POLV = EINAM( CHOICE )
 
        IF( INVPRO .NE. POLV ) THEN
 
            MESG = 'Enter logical name for ' //
     &              POLV ( 1:TRIMLEN( POLV ) ) // ' TO ' //
     &              INVPRO ( 1:TRIMLEN( INVPRO ) ) //
     &             ' CONVERSION file >> '
            KDEV = PROMPTFFILE(
     &             MESG, .TRUE., .TRUE., 'MSCNV', 'SPCPMAT' )
 
            KFLAG = .TRUE.
 
        ENDIF

C.......   If needed, read the volatile pollutant conversion file
 
        IF( KFLAG ) THEN

            CALL RDMCONV( KDEV, TDEV, 
     &                    CONVV0, FIPSV1, RDCLV1, CONVV1, NCNV1,
     &                            FIPSV2, RDCLV2, CONVV2, NCNV2,
     &                                    RDCLV3, CONVV3, NCNV3  )
 
        ENDIF

C.......   Read and process the speciation cross-reference file

        WRITE( *,92000 ) 
     &      ' ', 'Processing SPECIATION XREF file...', ' '

        IREC = 0
        I     =  0
        EFLAG = .FALSE.

101     CONTINUE        !  head of the XDEV-read loop

            IREC = IREC + 1
            READ( XDEV, 93010, END=111, IOSTAT=IOS ) LINE

            IF ( IOS .NE. 0 ) THEN

                EFLAG = .TRUE.
                WRITE( MESG,94050 ) 
     &              'Error ', IOS, 
     &              'reading SPECIATION XREF file at line', IREC
                CALL M3MESG( MESG )

            ENDIF

C.............  Skip x-ref line if not for MV
            IF ( LINE(17:17) .LE. '9' ) THEN
                GO TO 101

C.............  Otherwise convert character variables to integer
            ELSE
                VHT = STR2INT( LINE( 19:20 ) )
                RDT = STR2INT( LINE( 21:23 ) )
                EMT = INDEX1 ( LINE( 24:26 ), NEMODE, MEMODE )

                IF( RDT .NE. 0 ) THEN  ! Make sure we keep the defaults
                    NRDT = FIND1( RDT, NRCLAS, MROADS3 )
                ELSE
                    NRDT = 1
                ENDIF

            ENDIF

C.............  Warn and skip profile if MV-ASCT is blank
            IF( VHT .LT. 0 .AND. RDT .LT. 0 .AND. EMT .EQ. 0 ) THEN

                WRITE( MESG,94050 ) 'Speciation x-ref record', IREC,
     &                              'has blank motor-vehicles code.'
                CALL M3MESG( MESG )

                GO TO 101  !  to head of the XDEV-read loop

C.............  Skip profile if road class not used
            ELSE IF( NRDT .LT. 0 ) THEN
                GO TO 101  !  to head of the XDEV-read loop
                 
C.............  Skip profile if emission mode not used
            ELSE IF( EMT .EQ. 0 ) THEN
                GO TO 101  !  to head of the XDEV-read loop

            ELSE IF ( I .GE. MXREF ) THEN
                EFLAG = .TRUE.
                WRITE( MESG,94050 ) 
     &              'Max XREF table size ', MXREF,  
     &              'exceeded at line', IREC
                CALL M3MESG( MESG )

            ELSE
                I = I + 1
                INDEXA( I ) = I
                FIPSXA( I ) = STR2INT ( LINE( 1:5 ) )
                RDTSXA( I ) = RDT
                EMTSXA( I ) = EMT
                VHTSXA( I ) = VHT
                SPPROA( I ) = STR2INT ( LINE( 28:32 ) )
                VREACA( I ) = STR2REAL( LINE( 34:38 ) )
                ETHADA( I ) = STR2REAL( LINE( 55:59 ) )
                FRMADA( I ) = STR2REAL( LINE( 61:65 ) )
                ALLADA( I ) = STR2REAL( LINE( 67:71 ) )

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

                I0 = I0 + 1
                IF ( I0 .LE. MXVREF ) THEN
                    SPPROV( I0 ) = SPPROA( J )
                    EMTSXV( I0 ) = EMT
                    VHTSXV( I0 ) = VHT
                    VREACV( I0 ) = VREACA( J )
                    ETHADV( I0 ) = ETHADA( J )
                    FRMADV( I0 ) = FRMADA( J )
                    ALLADV( I0 ) = ALLADA( J )
                END IF

            ELSE IF ( FIP .EQ. 0  ) THEN !  FIP independent refs

                I1 = I1 + 1
                IF ( I1 .LE. MXRREF ) THEN
                    SPPROR( I1 ) = SPPROA( J )
                    RDTSXR( I1 ) = RDT
                    EMTSXR( I1 ) = EMT
                    VHTSXR( I1 ) = VHT
                    VREACR( I1 ) = VREACA( J )
                    ETHADR( I1 ) = ETHADA( J )
                    FRMADR( I1 ) = FRMADA( J )
                    ALLADR( I1 ) = ALLADA( J )
                END IF

            ELSE IF ( MOD( FIP, 1000 ) .EQ. 0 )  THEN ! state-dependent refs

                I2 = I2 + 1
                IF ( I2 .LE. MXSREF ) THEN
                    FIPSXS( I2 ) = FIP/1000
                    SPPROS( I2 ) = SPPROA( J )
                    RDTSXS( I2 ) = RDT
                    EMTSXS( I2 ) = EMT
                    VHTSXS( I2 ) = VHT
                    VREACS( I2 ) = VREACA( J )
                    ETHADS( I2 ) = ETHADA( J )
                    FRMADS( I2 ) = FRMADA( J )
                    ALLADS( I2 ) = ALLADA( J )
                END IF

            ELSE                ! fip > 0:  FIP-dep refs

                I3 = I3 + 1
                IF ( I3 .LE. MXFREF ) THEN
                    FIPSXF( I3 ) = FIP
                    SPPROF( I3 ) = SPPROA( J )
                    RDTSXF( I3 ) = RDT
                    EMTSXF( I3 ) = EMT
                    VHTSXF( I3 ) = VHT
                    VREACF( I3 ) = VREACA( J )
                    ETHADF( I3 ) = ETHADA( J )
                    FRMADF( I3 ) = FRMADA( J )
                    ALLADF( I3 ) = ALLADA( J )
                END IF

            END IF

122     CONTINUE        !  X-REF factor tables loop
        
        IF ( I0 .GT. MXVREF ) THEN

             WRITE( MESG,94050 )
     &           'Max RDT-independent table-size ', NVTYPE,  
     &           'exceeded by XREF file: count=', I0
             CALL M3EXIT( 'SPCMMAT', 0, 0, MESG, 2 )

        ELSE IF ( I1 .GT. MXRREF ) THEN

             WRITE( MESG,94050 )
     &           'Max FIP-independent table-size ', MXRREF,  
     &           'exceeded by XREF file: count=', I1
             CALL M3EXIT( 'SPCMMAT', 0, 0, MESG, 2 )

        ELSE IF ( I2 .GT. MXSREF ) THEN

             WRITE( MESG,94050 )
     &           'Max state-dependent table-size ', MXSREF,  
     &           'exceeded by XREF file: count=', I2
             CALL M3EXIT( 'SPCMMAT', 0, 0, MESG, 2 )

        ELSE IF ( I3 .GT. MXFREF ) THEN

             WRITE( MESG,94050 )
     &           'Max FIP-dependent table-size ', MXFREF,  
     &           'exceeded by XREF file: count=', I3
             CALL M3EXIT( 'SPCMMAT', 0, 0, MESG, 2 )

        ENDIF

        NREFV = I0              !  number of location-independent refs
        NREFR = I1              !  " road class-only    dependent refs
        NREFS = I2              !  " state-only         dependent refs
        NREFF = I3              !  " state-and-county   dependent refs
        EFLAG = .FALSE.
        
C.......   Compute the speciation coefficient matrix
C...........   ASSUMPTION:  everything sorted by FIP and ASC, or by XREF ID
C...........   Two stage process:  first try for cross-referece by FIP and ASC
C...........   Then if that fails, try by ASC only.

        DO  146  EMT = 1, NEMODE  ! By MV process
        DO  145  VHT = 1, NVTYPE  ! By VTYPE (X-REF handles bad EMT/VHT combo)
        DO  144  S   = 1, NMSRC   ! By MV source

            FIP = IFIP  ( S )
            RDT = IRCLAS( S )

            F = FIND4( FIP, RDT, EMT, VHT, 
     &                 NREFF, FIPSXF, RDTSXF, EMTSXF, VHTSXF )

            IF ( F .GT. 0 ) THEN        !  full match:

                SPCPRO = SPPROF( F )
                ALLADJ = ALLADF( F )
                FRMADJ = FRMADF( F )
                ETHADJ = ETHADF( F )
                VREACT = VREACF( F )
                GO TO  131

            END IF                      !  if full match

            SID = 1000 * ( FIP / 1000 ) !  state code
            F   = FIND4( SID, RDT, EMT, VHT, 
     &                   NREFS, FIPSXS, RDTSXS, EMTSXS, VHTSXS )
           
            IF ( F .GT. 0 ) THEN        !  state match
           
                SPCPRO = SPPROS( F )
                ALLADJ = ALLADS( F )
                FRMADJ = FRMADS( F )
                ETHADJ = ETHADS( F )
                VREACT = VREACS( F )
                GO TO  131

            END IF                      !  if state match
     
            F = FIND3( RDT, EMT, VHT, NREFR, RDTSXR, EMTSXR, VHTSXR )

            IF ( F .GT. 0 ) THEN        !  roadclass match

                SPCPRO = SPPROR( F )
                ALLADJ = ALLADR( F )
                FRMADJ = FRMADR( F )
                ETHADJ = ETHADR( F )
                VREACT = VREACR( F )
                GO TO  131

            END IF                      !  if roadclass match

            F = FIND2( EMT, VHT, NREFV, EMTSXV, VHTSXV )

            IF ( F .GT. 0 ) THEN        !  vehicletype/emissionmode-only match

                SPCPRO = SPPROV( F )
                ALLADJ = ALLADV( F )
                FRMADJ = FRMADV( F )
                ETHADJ = ETHADV( F )
                VREACT = VREACV( F )
                GO TO  131

            END IF              !  if roadclass match

            EFLAG = .TRUE.      !  if you get to here:  no match found
            WRITE( MESG,94050 ) 
     &          'Source', S, 
     &          '(FIP', FIP,
     &          'RDT', RDT,
     &          'VHT', VHT, ') not found'
            CALL M3MESG( MESG )
      
            GO TO 144      !  to end of loop

  131       CONTINUE            !  Cross-reference lookups succeeded

C...........   Finding the the correct volatile conversion factor is a
C...........       Four stage process:
C...........       by      FIP,Roadclass
C...........       by      State,Roadclass
C...........       by      Roadclass
C...........       by      Default

C...........   Initialize volatile conversion factor 
            VOLADJ = 1.0
 
C...........   Process volatile conversion factors if necessary
            IF( KFLAG ) THEN

                F = FIND2( FIP, RDT, NCNV1, FIPSV1, RDCLV1 )

                IF ( F .GT. 0 ) THEN      ! FIP, road class match
 
                    VOLADJ = CONVV1( F,VHT,EMT )
 
                ELSE                      ! no FIP, road class match

                    F = FIND2( (FIP/1000)*1000, RDT, NCNV2, 
     &                         FIPSV2, RDCLV2 )

                    IF ( F .GT. 0 ) THEN
 
                        VOLADJ = CONVV2( F,VHT,EMT )
 
                    ELSE                  !  no state, road class match

                        F = FIND2( RDT, NCNV3, RDCLV3 )

                        IF ( F .GT. 0 ) THEN
 
                            VOLADJ = CONVV3( F,VHT,EMT )
 
                        ELSE              !  no road class match

                            VOLADJ = CONVV0( VHT,EMT )

                        ENDIF  ! road class match
                    ENDIF      ! state-road class match
                ENDIF          ! FIP-road class match

            ENDIF          ! KFLAG is true

            I = FIND1( SPCPRO, NPROFS, SXREF )

            IF ( I .LE. 0 ) THEN        !  error:  XREF index not found

                EFLAG = .TRUE.
                WRITE( MESG,94050 ) 
     &              'XREF for FIPS', FIP,
     &              'RDT', RDT,
     &              'VHT', VHT, 
     &              'profile code', SPCPRO, 'not found'
                CALL M3MESG( MESG )

                GO TO  144      !  to end of loop

            END IF      !  if i nonpositive (cross-reference not found)

C...........   Loop through species that exist for this MV process
            DO 141 P = 1, NSP( EMT )

                K = IDXVAR( EMT, P )
                V = IDXSPC( EMT, P )

C.................  Set inventory-species specific variables
C.................  Use VREACT only if the volatile inventory pollutant is THC
                IF( ISPC_PRC( EMT,P ) .EQ. 'THC' ) THEN
                    FAC = VREACT * ALLADJ * VOLADJ

C.................  Use the volatile correct factor only if the current
C.................  inventory species is volatile
                ELSEIF( ISPC_PRC( EMT,P ) .EQ. POLV ) THEN
                    FAC = ALLADJ * VOLADJ

                ELSE
                    FAC = ALLADJ

                ENDIF

C.................  Set factors with consideration to model species

                IF( MSPC_PRC( EMT,P ) .EQ. 'FORM' ) THEN
                    SPCFACS( S,VHT,K ) = 
     &                          FRMADJ * FAC * SCOEF( I,V )
            
                ELSEIF( MSPC_PRC( EMT,P ) .EQ. 'ETH' ) THEN
                    SPCFACS( S,VHT,K ) = 
     &                          ETHADJ * FAC * SCOEF( I,V )

                ELSE
                    SPCFACS( S,VHT,K ) = FAC * SCOEF( I,V )

                END IF

141         CONTINUE      ! end loop on VOC-related species

144     CONTINUE      ! end loop on mobile sources (S)
145     CONTINUE      ! end loop on vehicle types (VHT)
146     CONTINUE      ! end loop on emissions modes, computing spc matrix (EMT)

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

93010   FORMAT( A )

C...........   Internal buffering formats............ 94xxx

94020   FORMAT( I7 )

94030   FORMAT( I3 )

94040   FORMAT( G10.3 )

94050   FORMAT( 10 ( A, :, I10, :, 2X ) )

        END

