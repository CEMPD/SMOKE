C ****************************************************************************
C *
C * Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
C *                System
C * File: @(#)$Id$
C *
C * COPYRIGHT (C) 1996, MCNC--North Carolina Supercomputing Center
C * All Rights Reserved
C *
C * See file COPYRIGHT for conditions of use.
C *
C * Environmental Programs Group
C * MCNC--North Carolina Supercomputing Center
C * P.O. Box 12889
C * Research Triangle Park, NC  27709-2889
C *
C * env_progs@mcnc.org
C *
C * Pathname: $Source$
C * Last updated: $Date$ 
C *
C ****************************************************************************/

        PROGRAM CTLPMAT

C***********************************************************************
C  program body starts at line 171
C
C  DESCRIPTION:
C       Construct point source control matrix from data 
C       contained in prototype POINT and EPS-style control files.
C
C  PRECONDITIONS REQUIRED:
C       Sorted, cut-down control file.
C       M3IO point-sources file from program RAWPOINT
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C       Models-3 I/O
C       FIND1, FIND2, PROMPTFFILE, PROMPTMFILE, GETYN, TRIMLEN
C
C  REVISION  HISTORY:
C       Prototype  6/95 by CJC.
C
C***********************************************************************

      IMPLICIT NONE

C...........   INCLUDES:
        
        INCLUDE 'PTDIMS3.EXT'   !  point-source dimensioning parameters
        INCLUDE 'CHDIMS3.EXT'   !  emis chem parms (inventory + model)
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
        INCLUDE 'FDESC3.EXT'    !  I/O API file description data structures.


C...........   EXTERNAL FUNCTIONS and their descriptions:
        
        INTEGER         FIND1, FIND2, FIND3, FIND4
        LOGICAL         GETYN
        LOGICAL         GETPKT
        INTEGER         INDEX1
        INTEGER         PROMPTFFILE
        CHARACTER*16    PROMPTMFILE
        LOGICAL         SKIPCTL
        INTEGER         STR2INT
        REAL            STR2REAL
        INTEGER         TRIMLEN
        
        EXTERNAL        FIND1, FIND2, FIND3, FIND4, GETYN, GETPKT,
     &                  INDEX1, PROMPTFFILE, PROMPTMFILE, 
     &                  SKIPCTL, STR2INT, STR2REAL, TRIMLEN


C...........   PARAMETERS and their descriptions:

        INTEGER		MAXPKT		!  max # of packets in control file
        PARAMETER     ( MAXPKT = 20 )


C...........   LOCAL VARIABLES and their descriptions:
C...........   NOTE that ASC (point-source-category) ID's are 10-digit
C...........   unsigned integers which may be treated as a leading 7-digit
C...........   field, and a trailing 3-digit field.  *7 and *3 arrays below
C...........   follow this scheme with parallel arrays

        INTEGER     FIP         !  FIP (state+county) code
        INTEGER     SID         !  state ID (leading-2 digits of FIP code)
        INTEGER     SCC         !  source category code
        INTEGER     SIC         !  standard industrial category
        INTEGER     PLT         !  plant ID
        INTEGER     STK         !  stack ID
        REAL        CNOX        !  NOX control factor
        REAL        CVOC        !  VOC control factor
        REAL        CCO         !  OC control factor
        REAL        FNOX        !  NOX control factor
        REAL        FVOC        !  VOC control factor
        REAL        FCO         !  OC control factor
        REAL        PFAC        !  projection factor
        REAL        TNOX        !  NOX control threshold
        REAL        TVOC        !  VOC control threshold
        INTEGER     LFIP        !  last FIP (state+county) code
        INTEGER     LSID        !  last state ID
        INTEGER     NFIP        !  subscript for JFIP
        INTEGER     NSID        !  subscript for JSID
        INTEGER     DYR         !  desination year for projection factor
        INTEGER     YEAR        !  inventory year for projection factor
                                 
C.......   Control matrix:
        
        REAL        CMAT( NPSRC, NIPOL )    !  control coefficients


C...........   Point Sources Table
        
        INTEGER     IFIP  ( NPSRC )  !  source FIPS (county) ID
        INTEGER     ISCC  ( NPSRC )  !  source SCC
        INTEGER     ISIC  ( NPSRC )  !  source SIC
        INTEGER     ISI3  ( NPSRC )  !  source SIC (leading 3 digits)
        INTEGER     ISI2  ( NPSRC )  !  source SIC (leading 2 digits)
        INTEGER     IPLANT( NPSRC )  !  Plant ID
        INTEGER     ISTACK( NPSRC )  !  Stack ID
        INTEGER     INVYR ( NPSRC )  !  Inventory year
        REAL        RULPEN( NPSRC, NIPOL )  !  rule penetration   fraction
        REAL        RULEFF( NPSRC, NIPOL )  !  rule effectiveness fraction
        REAL        CTLEFF( NPSRC, NIPOL )  !  rule effectiveness fraction
        REAL        PRJFAC( NPSRC )  !  projection factor.
        REAL        NOX   ( NPSRC )  !  NOx emissions values (tons/yr)
        REAL        VOC   ( NPSRC )  !  VOC emissions values (tons/yr)
        REAL        CO    ( NPSRC )  !  VOC emissions values (tons/yr)

C.......   Auxilliary tables:

        INTEGER     SIDS  ( NPSID )  !  source FIPS (state-only) ID
        INTEGER     JSID  ( NPSID )  !  Subscript at which state ID occurs
        INTEGER     FFIP  ( NPFIP )  !  source FIPS list
        INTEGER     JFIP  ( NPFIP )  !  Subscript at which FIP occurs


C...........   Other local variables
        
        INTEGER         S, I, J      !  loop counters.
        INTEGER         CDEV         !  for control file
        INTEGER         LDEV         !  log-device
        INTEGER         RDEV         !  report-file device
        INTEGER         TDEV         !  device for report-these-sources file
        LOGICAL         PFLAG        !  true if projection entries have years
        LOGICAL         RFLAG        !  flag:  make report-file
        LOGICAL         TFLAG( NPSRC )  !  flags:  track these sources
        
        INTEGER         IOS          !  I/O status
        INTEGER         IREC         !  input line (record) number

        REAL            T, U, V      !  scratch-subexpression variables

        CHARACTER*16    ENAME   !  logical name for point-source input file
        CHARACTER*16    CNAME   !  logical name for control     input file
        CHARACTER*16    MNAME   !  logical name for ctl matrix output file
        CHARACTER*16    VNAME   !  buffer for variable-names
        CHARACTER*240   LINE    !  input line buffer
        CHARACTER*256   MESG    !  error message buffer
        
        INTEGER         NPACKETS
        CHARACTER*20    PNAMES( MAXPKT )
        CHARACTER*20    PNAME

        
C***********************************************************************
C   begin body of program CTLPMAT
        
        LDEV = INIT3()
        
        WRITE( *,92000 ) 
     &  ' ',
     &  'Program CTLPMAT to take the netCDF POINT SOURCES file, an',
     &  '(ASCII) CONTROL FILE, an optional TRACK LIST of sources to',
     &  'be tracked in the report, and produce the CONTROL MATRIX ',
     &  'file and an optional REPORT file.',
     &  ' '
        WRITE( *,92000 ) 
     &  'You will need to have set the logical names for the input and',
     &  'output files (using "setenv <logicalname> <pathname>") prior',
     &  'to program launch, and will need to enter these logical names',
     &  'at the appropriate prompts.  Input files must have been ',
     &  'sorted as indicated, prior to program execution.',
     &  ' '
        WRITE( *,92000 ) 
     &  'You may use END_OF-FILE (control-D) to quit the program',
     &  'during logical-name entry.  Default responses to prompts',
     &  'are indicated in brackets [LIKE THIS].',
     &  ' '

        IF ( .NOT. GETYN( 'Continue with program?', .TRUE. ) ) THEN
            WRITE( *,92000 ) ' ', 'Ending program CTLPMAT', ' '
            CALL EXIT( 0 )
        END IF


C.......   Get file name; open input point sources, control definitions files

        ENAME = PROMPTMFILE( 
     &          'Enter logical name for POINT SOURCE input file',
     &          FSREAD3, 'PNTS', 'CTLPMAT' )

        CDEV = PROMPTFFILE( 
     &           'Enter logical name for CONTROL file',
     &           .TRUE., .TRUE., 'PCTRL', 'CTLPMAT' )

        RDEV = PROMPTFFILE( 
     &           '...for REPORT file, or "NONE"',
     &           .FALSE., .TRUE., 'REPORT', 'CTLPMAT' )

        RFLAG = ( RDEV .GT. 0 )

        IF ( RFLAG ) THEN

            TDEV = PROMPTFFILE( 
     &           '...for TRACKING LIST, or "ALL"',
     &           .TRUE., .TRUE., 'TRACK', 'CTLPMAT' )

        END IF

C.......   Read in description of point-source file
        
        IF ( .NOT. DESC3( ENAME ) ) THEN
            CALL M3EXIT( 'CTLPMAT', 0, 0, 
     &          'Could not get description of file ' // ENAME, 3 )
        END IF
        
        
C.......   Open control matrix file:  reuse much of POINT file description
        
        NVARS3D = 3
        FTYPE3D = GRDDED3

        DO 11 V = 1, NIPOL
            VNAME3D( V ) = 'CONTROL_' // EINAM( V )
            UNITS3D( V ) = 'n/a'
            VDESC3D( V ) = 'Control factor for this source, variable ' 
     &                     // EINAM( V )( 1:TRIMLEN( EINAM( V ) ) )
            VTYPE3D( V ) = M3REAL
11      CONTINUE

        FDESC3D( 1 ) = 'Matrix of control factors, by source'
        DO  22  I = 2, MXDESC3
            FDESC3D( I ) = ' '
22      CONTINUE

        MNAME = PROMPTMFILE( 
     &          '...for CONTROL MATRIX output file',
     &          FSUNKN3, 'PCMAT', 'CTLPMAT' )


C............................................................................
C...........   Read in emissions indexing data IFIP, ISCC, ISIC;
C...........   existing-control data CTLEFF, RULEFF, RULPEN:
        
        IF ( .NOT. READ3( ENAME, 'IFIP', ALLAYS3, 0, 0, IFIP ) ) THEN
            CALL M3EXIT( 'CTLPMAT', 0, 0, 
     &                  'Error reading FIP from file' // ENAME, 2 )
        END IF

        IF ( .NOT. READ3( ENAME, 'ISCC', ALLAYS3, 0, 0, ISCC ) ) THEN
            CALL M3EXIT( 'CTLPMAT', 0, 0, 
     &                  'Error reading ISCC from file' // ENAME, 2 )
        END IF

        IF ( .NOT. READ3( ENAME, 'ISIC', ALLAYS3, 0, 0, ISIC ) ) THEN
            CALL M3EXIT( 'CTLPMAT', 0, 0, 
     &                  'Error reading ISIC from file' // ENAME, 2 )
        END IF

        IF ( .NOT. READ3( ENAME, 'IPLANT', ALLAYS3,0,0, IPLANT ) ) THEN
            CALL M3EXIT( 'CTLPMAT', 0, 0, 
     &                  'Error reading IPLANT from file' // ENAME, 2 )
        END IF

        IF ( .NOT. READ3( ENAME, 'ISTACK', ALLAYS3,0,0, ISTACK ) ) THEN
            CALL M3EXIT( 'CTLPMAT', 0, 0, 
     &                  'Error reading ISTACK from file' // ENAME, 2 )
        END IF

        IF ( .NOT. READ3( ENAME, 'INVYR', ALLAYS3,0,0, INVYR ) ) THEN
            CALL M3EXIT( 'CTLPMAT', 0, 0, 
     &                  'Error reading INVYR from file' // ENAME, 2 )
        END IF

        DO  33  I = 1, NIPOL

            VNAME = 'CTLEFF_' // EINAM( I )
            IF ( .NOT. READ3( ENAME, VNAME, ALLAYS3, 
     &                        0, 0, CTLEFF( 1,I ) ) ) THEN
                CALL M3EXIT( 'CTLPMAT', 0, 0, 
     &               'Error reading CTLEFF from file' // ENAME, 2 )
            END IF

            VNAME = 'RULEFF_' // EINAM( I )
            IF ( .NOT. READ3( ENAME, VNAME, ALLAYS3, 
     &                        0, 0, RULEFF( 1,I ) ) ) THEN
                CALL M3EXIT( 'CTLPMAT', 0, 0, 
     &               'Error reading RULEFF from file' // ENAME, 2 )
            END IF

            VNAME = 'RULPEN_' // EINAM( I )
            IF ( .NOT. READ3( ENAME, VNAME, ALLAYS3, 
     &                        0, 0, RULPEN( 1,I ) ) ) THEN
                CALL M3EXIT( 'CTLPMAT', 0, 0, 
     &               'Error reading RULPEN from file' // ENAME, 2 )
            END IF

33      CONTINUE                !  end loop reading ctleff, ruleff, rulpen of V

        IF ( .NOT. READ3( ENAME, 'CO', ALLAYS3, 0, 0, CO ) ) THEN
            CALL M3EXIT( 'CTLPMAT', 0, 0, 
     &                  'Error reading CO from file' // ENAME, 2 )
        END IF

        IF ( .NOT. READ3( ENAME, 'NOX', ALLAYS3, 0, 0, NOX ) ) THEN
            CALL M3EXIT( 'CTLPMAT', 0, 0, 
     &                   'Error reading NOX from file' // ENAME, 2 )
        END IF

        IF ( .NOT. READ3( ENAME, 'VOC', ALLAYS3, 0, 0, VOC ) ) THEN
            CALL M3EXIT( 'CTLPMAT', 0, 0, 
     &                  'Error reading VOC from file' // ENAME, 2 )
        END IF

        
C............................................................................
C.......   Initialization; precompute states for each source
C.......   Note also:  convert RULEFF, RULPEN, CTLEFF from
C.......   percent to purely-numeric fraction
            
        DO  44  I = 1, NIPOL
        DO  43  S = 1, NPSRC

            U = RULEFF( S,I )
            IF ( U .LT. AMISS3 ) THEN
                RULEFF( S,I ) = 1.0
            ELSE
                RULEFF( S,I ) = 0.01 * U
            END IF
            V = RULPEN( S,I )
            IF ( V .LT. AMISS3 ) THEN
                RULPEN( S,I ) = 1.0
            ELSE
                RULPEN( S,I ) = 0.01 * V
            END IF
            T = CTLEFF( S,I )
            IF ( T .LT. AMISS3 ) THEN
                CTLEFF( S,I ) = 0.0
            ELSE
                CTLEFF( S,I ) = 0.01 * T
            END IF
            
            CMAT( S,I ) = 1.0 / 
     &          ( 1.0 - CTLEFF( S,I ) * RULEFF( S,I ) * RULPEN( S,I ) )
 
43      CONTINUE        !  end loop on sources S
44      CONTINUE        !  end loop on sources S

        LFIP = -1
        LSID = -1
        NFIP =  0
        NSID =  0

        DO  55  S = 1, NPSRC

            PRJFAC( S ) = 1.0

            FIP       = IFIP( S )
            SID       = FIP / 1000
            
            IF ( FIP .NE. LFIP ) THEN
                NFIP = NFIP + 1
                IF ( NFIP .GE. NPFIP ) THEN
                    CALL M3EXIT( 'CTLPMAT', 0, 0, 
     &                           'FIP-table maximum size exceeded', 2 )
                ELSE
                    JFIP( NFIP ) = S
                    FFIP( NFIP ) = FIP
                    LFIP         = FIP
                END IF
            END IF

            IF ( SID .NE. LSID ) THEN
                NSID = NSID + 1
                IF ( NSID .GE. NPSID ) THEN
                    CALL M3EXIT( 'CTLPMAT', 0, 0, 
     &                           'SID-table maximum size exceeded', 2 )
                ELSE
                    JSID( NSID ) = S
                    SIDS( NSID ) = SID
                    LSID         = SID
                END IF
            END IF
            
C...........   Compute 3-digit, 2-digit SIC hierarchy codes:
            
            ISI3( S ) = ISIC( S ) /10
            ISI2( S ) = ISIC( S ) /100

55      CONTINUE        !  end loop on sources S


        JSID( NSID+1 ) = NPSRC + 1      !  sentinel
        JFIP( NSID+1 ) = NPSRC + 1      !  sentinel


C............................................................................
C.......   Process list of FIP-SCC pairs to be tracked in the 
C.......   REPORT file (if any).
            
        DO  66  S = 1, NPSRC
            TFLAG( S ) = .FALSE.
66      CONTINUE        !  end loop on sources S

        IF ( RFLAG ) THEN

            IF ( TDEV .GT. 0 ) THEN

                IREC = 0

75              CONTINUE

                    IREC = IREC + 1
                    READ( TDEV,*, END = 77, IOSTAT=IOS ) FIP, SCC

                    IF ( IOS .NE. 0 ) THEN
                        WRITE( MESG,94910 )
     &                      'Error reading TRACKING LIST', IREC, IOS
                        CALL M3EXIT( 'CTLPMAT', 0, 0, MESG, 1 )
                    END IF
                    
                    I = FIND1( FIP, NFIP, FFIP )

                    IF ( I .EQ. 0 )  GO TO 75

                    DO  76  S = JFIP( I ), JFIP( I+1 ) - 1
                        IF( SCC .EQ. ISCC ( S ) )  TFLAG( S ) = .TRUE.
76                  CONTINUE        !  end loop on sources S

                    GO TO 75

77              CONTINUE        !  end loop reading file:  target of "END=..."

            ELSE IF ( TDEV .EQ. -3 ) THEN       ! "ALL"

                DO  88  S = 1, NPSRC
                    TFLAG( S ) = .TRUE.
88              CONTINUE        !  end loop on sources S

            ELSE

                CALL M3EXIT( 'CTLPMAT', 0, 0, 
     &                       'Error opening TRACKING LIST file', 3 )

            END IF              ! if tdev > 0, or tdev = -3, or not.

        END IF                  !  if rflag


C............................................................................
C...............   CONTROL PACKET PROCESSING
C.......   Process packets using "painters algorithm" -- most specific
C.......   controls LAST (overwriting less specific controls as appropriate)
C.......   in the order
C.......        SUBGRID packet       -- IGNORED
C.......        SOURCE CATEGORY packet
C.......        FIPS CODE packet
C.......        PROCESS CODE packet  -- IGNORED
C.......        PROFILE CODE packet  -- IGNORED
C.......        CONTROL CODE packet  -- IGNORED
C.......        POD CODE packet      -- IGNORED
C.......        ACTIVITY CODE packet -- IGNORED
C.......        ALLOWABLE packet
C.......        PROJECT AMS packet   -- IGNORED
C.......        PROJECT PTS packet
C.......        RULE PENETRATION packet
C.......        RULE EFFECTIVENESS packet
C.......        CONTROL EFFICIENCY packet
C.......        RACT packet
C.......        MACT packet
C.......        CTG packet
C.......   
C.......   First, get list of packet names in the control file
        
        IREC     = 1
        NPACKETS = 0

90      CONTINUE

            READ( CDEV,93000, IOSTAT=IOS, END=99 ) LINE
            IF( IOS .NE. 0 ) THEN
                WRITE( MESG,94920 ) 
     &                  'Error reading', CNAME( 1:TRIMLEN( CNAME ) ),
     &                  IREC, IOS
                CALL M3EXIT( 'CTLPMAT', 0, 0, MESG, 1 )
            END IF
            CALL UPCASE( LINE )
            IF ( 0 .NE. INDEX( LINE, '/END/' ) ) GO TO 90

            I = INDEX( LINE, '/' )
            IF ( I .EQ. 0 ) GO TO 90

            J = INDEX( LINE( I+1 : 240 ), '/' )

            IF ( J .EQ. 0 ) THEN

                WRITE( MESG,94920 )
     &              'Malformed packet in ',
     &              CNAME( 1:TRIMLEN( CNAME ) ),
     &              IREC
                CALL M3EXIT( 'CTLPMAT', 0, 0, MESG, 1 )

            ELSE IF ( J .GT. 20 ) THEN

                WRITE( MESG,94920 )
     &              'Packet name too long in ',
     &              CNAME( 1:TRIMLEN( CNAME ) ),
     &              IREC
                CALL M3EXIT( 'CTLPMAT', 0, 0, MESG, 1 )

            END IF

            IF ( NPACKETS .LT. MAXPKT ) THEN
                NPACKETS = NPACKETS + 1
                PNAMES( NPACKETS ) = LINE( I:I+J )
                WRITE( LDEV,92000 )
     &              'Control packet ', PNAMES( NPACKETS ), ' present'
                GO TO 90
            ELSE
                WRITE( MESG,94920 )
     &              'Max packet count exceeded in ',
     &              CNAME( 1:TRIMLEN( CNAME ) ), IREC
                CALL M3EXIT( 'CTLPMAT', 0, 0, MESG, 1 )
            END IF
            
99      CONTINUE        !  end loop reading packet names.
            

C............................................................................
C.......   Process SOURCE CATEGORY packet of control file:

        PNAME = '/SOURCE CATEGORY/'
        IF ( 0 .EQ. INDEX1( PNAME, NPACKETS, PNAMES ) )  GO TO 109

        IF ( .NOT.  GETPKT( PNAME, CDEV ) ) THEN
            CALL M3EXIT( 'CTLPMAT', 0, 0,
     &                   'Could not find ' // PNAME // ' in ' // CNAME, 
     &                   2 )
        END IF

101     CONTINUE        !  head of loop reading lines in the SOURCE CATEGORY packet

            READ( CDEV,93000,IOSTAT=IOS,END=900 ) LINE
            IREC = IREC + 1
            IF( IOS .NE. 0 ) THEN
                WRITE( MESG,94920 ) 
     &              'Error reading', CNAME( 1:TRIMLEN( CNAME ) ), 
     &              IREC, IOS
                CALL M3EXIT( 'CTLPMAT', 0, 0, MESG, 1 )
            END IF

            CALL UPCASE( LINE )
            IF( 0 .NE. INDEX( LINE, 'END' ) )  GO TO 109

            CALL LUSTR( LINE( 21:30 ) )
            FIP  = STR2INT ( LINE(  1: 5 ) )
            SCC  = STR2INT ( LINE( 21:30 ) )
            CNOX = 1.0 - STR2REAL( LINE( 32:41 ) )
            CVOC = 1.0 - STR2REAL( LINE( 42:51 ) )
            CCO  = 1.0 - STR2REAL( LINE( 52:61 ) )

            IF ( FIP .LE. 0 ) THEN      !  all states.

                DO  102 S = 1, NPSRC

                    IF( SCC .NE. ISCC( S ) ) GO TO 102  ! n/a control

                    IF ( TFLAG( S ) ) 
     &                  WRITE( RDEV,93010 ) 
     &                      PNAME, IFIP( S ), ISCC( S ), 
     &                      IPLANT( S ), ISTACK( S ),
     &                      CTLEFF( S,1 ), CCO,
     &                      CTLEFF( S,2 ), CNOX,
     &                      CTLEFF( S,3 ), CVOC

                    CTLEFF( S,1 ) = CCO
                    CTLEFF( S,2 ) = CNOX
                    CTLEFF( S,3 ) = CVOC

102             CONTINUE

            ELSE IF ( MOD( FIP, 1000 ) .EQ. 0 ) THEN

                SID = FIP / 1000        !  state code only
                I   = FIND1( SID, NSID, SIDS )

                IF ( I .EQ. 0 ) GO TO 101

                DO  103 S = JSID( I ), JSID( I+1 ) - 1

                    IF( SCC .NE. ISCC( S ) ) GO TO 103  ! n/a control

                    IF ( TFLAG( S ) ) 
     &                  WRITE( RDEV,93010 ) 
     &                      PNAME, IFIP( S ), ISCC( S ), 
     &                      IPLANT( S ), ISTACK( S ),
     &                      CTLEFF( S,1 ), CCO,
     &                      CTLEFF( S,2 ), CNOX,
     &                      CTLEFF( S,3 ), CVOC

                    CTLEFF( S,1 ) = CCO
                    CTLEFF( S,2 ) = CNOX
                    CTLEFF( S,3 ) = CVOC

103             CONTINUE

            ELSE        !  fip-SCC specific control

                I = FIND1( FIP, NFIP, FFIP )

                IF ( I .LE. 0 ) GO TO 101

                DO  104 S = JFIP( I ), JFIP( I+1 ) - 1

                    IF( SCC .NE. ISCC( S ) ) GO TO 104  ! n/a control

                    IF ( TFLAG( S ) ) 
     &                  WRITE( RDEV,93010 ) 
     &                      PNAME, IFIP( S ), ISCC( S ), 
     &                      IPLANT( S ), ISTACK( S ),
     &                      CTLEFF( S,1 ), CCO,
     &                      CTLEFF( S,2 ), CNOX,
     &                      CTLEFF( S,3 ), CVOC

                    CTLEFF( S,1 ) = CCO
                    CTLEFF( S,2 ) = CNOX
                    CTLEFF( S,3 ) = CVOC

104             CONTINUE

            END IF      ! if fip is 0, or nn000, or not.

            GO TO  101  !  end loop reading control-line in SOURCE CATEGORY packet

109     CONTINUE        !  exit from / skip around loop 101


C............................................................................
C.......   Process FIPS CODE packet of control file:

        PNAME = '/FIPS CODE/'
        IF( 0 .EQ. INDEX1( PNAME, NPACKETS, PNAMES ) ) GO TO 119

        IF ( .NOT. GETPKT( PNAME, CDEV ) ) THEN
            CALL M3EXIT( 'CTLPMAT', 0, 0,
     &                   'Could not find ' // PNAME // ' in ' // CNAME, 
     &                   2 )
        END IF

111     CONTINUE        !  head of loop reading lines in the IFIP CODE packet

            READ( CDEV,93000,IOSTAT=IOS,END=900 ) LINE
            IREC = IREC + 1
            IF( IOS .NE. 0 ) THEN
                WRITE( MESG,94920 ) 
     &              'Error reading', CNAME( 1:TRIMLEN( CNAME ) ), 
     &              IREC, IOS
                CALL M3EXIT( 'CTLPMAT', 0, 0, MESG, 1 )
            END IF

            CALL UPCASE( LINE )
            IF( 0 .NE. INDEX( LINE, 'END' ) )  GO TO 119

            FIP  = STR2INT ( LINE(  1: 5 ) )
            CNOX = 1.0 - STR2REAL( LINE( 16:25 ) )
            CVOC = 1.0 - STR2REAL( LINE( 26:35 ) )
            CCO  = 1.0 - STR2REAL( LINE( 36:45 ) )

            IF ( FIP .LE. 0 ) THEN      !  all states.

                DO  112 S = 1, NPSRC

                    IF ( TFLAG( S ) ) 
     &                  WRITE( RDEV,93010 ) 
     &                      PNAME, IFIP( S ), ISCC( S ), 
     &                      IPLANT( S ), ISTACK( S ),
     &                      CTLEFF( S,1 ), CCO,
     &                      CTLEFF( S,2 ), CNOX,
     &                      CTLEFF( S,3 ), CVOC

                    CTLEFF( S,1 ) = CCO
                    CTLEFF( S,2 ) = CNOX
                    CTLEFF( S,3 ) = CVOC

112             CONTINUE

            ELSE IF ( MOD( FIP, 1000 ) .EQ. 0 ) THEN

                SID = FIP / 1000        !  state code only
                I   = FIND1( SID, NSID, SIDS )

                IF ( I .EQ. 0 ) GO TO 111

                DO  113 S = JSID( I ), JSID( I+1 ) - 1

                    IF ( TFLAG( S ) ) 
     &                  WRITE( RDEV,93010 ) 
     &                      PNAME, IFIP( S ), ISCC( S ), 
     &                      IPLANT( S ), ISTACK( S ),
     &                      CTLEFF( S,1 ), CCO,
     &                      CTLEFF( S,2 ), CNOX,
     &                      CTLEFF( S,3 ), CVOC

                    CTLEFF( S,1 ) = CCO
                    CTLEFF( S,2 ) = CNOX
                    CTLEFF( S,3 ) = CVOC

113             CONTINUE

            ELSE        !  fip-ASC specific control

                I = FIND1( FIP, NFIP, FFIP )

                IF ( I .EQ. 0 ) GO TO 111

                DO  114 S = JFIP( I ), JFIP( I+1 ) - 1

                    IF ( TFLAG( S ) ) 
     &                  WRITE( RDEV,93010 ) 
     &                      PNAME, IFIP( S ), ISCC( S ), 
     &                      IPLANT( S ), ISTACK( S ),
     &                      CTLEFF( S,1 ), CCO,
     &                      CTLEFF( S,2 ), CNOX,
     &                      CTLEFF( S,3 ), CVOC

                    CTLEFF( S,1 ) = CCO
                    CTLEFF( S,2 ) = CNOX
                    CTLEFF( S,3 ) = CVOC

114             CONTINUE

            END IF      ! if fip is 0, or nn000, or not.

            GO TO  111  !  end loop reading control-lines in FIPS CODE packet

119     CONTINUE        !  exit from / skip around loop 261


C............................................................................
C.......   Process ALLOWABLE packet of control file:

        PNAME = '/ALLOWABLE/'
        IF( 0 .EQ. INDEX1( PNAME, NPACKETS, PNAMES ) ) GO TO 140

        IF ( .NOT. GETPKT( PNAME, CDEV ) ) THEN
            CALL M3EXIT( 'CTLPMAT', 0, 0,
     &                   'Could not find ' // PNAME // ' in ' // CNAME, 
     &                   2 )
        END IF

121     CONTINUE        !  head of loop reading lines in the RULPEN packet

            READ( CDEV,93000,IOSTAT=IOS,END=900 ) LINE
            IREC = IREC + 1
            IF( IOS .NE. 0 ) THEN
                WRITE( MESG,94920 ) 
     &              'Error reading', CNAME( 1:TRIMLEN( CNAME ) ), 
     &              IREC, IOS
                CALL M3EXIT( 'CTLPMAT', 0, 0, MESG, 1 )
            END IF

            CALL UPCASE( LINE )
            IF( 0 .NE. INDEX( LINE, 'END' ) )  GO TO 140

            IF( LINE( 43:44 ) .NE. 'AE' ) THEN
                WRITE( MESG,94910 ) 
     &              'ALLOWABLE control type "'//LINE( 43:44 ) //
     &              '" not implemented', IREC
                CALL M3MESG( MESG )
                IF ( RFLAG ) WRITE( RDEV,92000 ) MESG( 1:80 )
                GO TO 121
            END IF

            FIP  = STR2INT ( LINE(  1: 5 ) )
            PLT  = STR2INT ( LINE( 13:17 ) )
            STK  = STR2INT ( LINE( 19:21 ) )
            SIC  = STR2INT ( LINE( 27:30 ) )
            SCC  = STR2INT ( LINE( 32:41 ) )
            CNOX = 1.0 - STR2REAL( LINE( 46:55 ) )
            CVOC = 1.0 - STR2REAL( LINE( 56:65 ) )
            CCO  = 1.0 - STR2REAL( LINE( 66:75 ) )

            IF ( SCC .GT. 0 ) THEN

                IF ( SIC .GT. 0 ) THEN

                    IF ( FIP .LE. 0 ) THEN      !  all states.

                        DO  122 S = 1, NPSRC

                            IF( SCC .NE. ISCC( S ) ) GO TO 122  ! n/a control
                            IF( SIC .NE. ISIC( S ) .AND.
     &                          SIC .NE. ISI3( S ) .AND.
     &                          SIC .NE. ISI2( S ) ) GO TO 122

                            FNOX = ( CNOX - NOX( S ) ) /
     &                             ( RULPEN( S,1 ) 
     &                             * RULEFF( S,1 ) 
     &                             * CNOX )

                            FVOC = ( CVOC - VOC( S ) ) /
     &                             ( RULPEN( S,2 ) 
     &                             * RULEFF( S,2 ) 
     &                             * CVOC )

                            FCO  = ( CCO  - CO ( S ) ) /
     &                             ( RULPEN( S,3 ) 
     &                             * RULEFF( S,3 ) 
     &                             * CCO  )

                            IF ( TFLAG( S ) ) 
     &                          WRITE( RDEV,93010 ) 
     &                              PNAME, IFIP( S ), ISCC( S ), 
     &                              IPLANT( S ), ISTACK( S ),
     &                              CTLEFF( S,1 ), FCO,
     &                              CTLEFF( S,2 ), FNOX,
     &                              CTLEFF( S,3 ), FVOC

                            CTLEFF( S,1 ) = FNOX
                            CTLEFF( S,2 ) = FVOC
                            CTLEFF( S,3 ) = FCO

122                     CONTINUE

                    ELSE IF ( MOD( FIP, 1000 ) .EQ. 0 ) THEN	!  this state

                        SID = FIP / 1000        !  state code only
                        I   = FIND1( SID, NSID, SIDS )
                        IF ( I .EQ. 0 ) GO TO 121

                        DO  123 S = JSID( I ), JSID( I+1 ) -1

                            IF( SCC .NE. ISCC ( S ) ) GO TO 123  ! n/a control
                            IF( SIC .NE. ISIC( S ) .AND.
     &                          SIC .NE. ISI3( S ) .AND.
     &                          SIC .NE. ISI2( S ) ) GO TO 123

                            FNOX = ( CNOX - NOX( S ) ) /
     &                             ( RULPEN( S,1 ) 
     &                             * RULEFF( S,1 ) 
     &                             * CNOX )

                            FVOC = ( CVOC - VOC( S ) ) /
     &                             ( RULPEN( S,2 ) 
     &                             * RULEFF( S,2 ) 
     &                             * CVOC )

                            FCO  = ( CCO  - CO ( S ) ) /
     &                             ( RULPEN( S,3 ) 
     &                             * RULEFF( S,3 ) 
     &                             * CCO  )

                            IF ( TFLAG( S ) ) 
     &                          WRITE( RDEV,93010 ) 
     &                              PNAME, IFIP( S ), ISCC( S ), 
     &                              IPLANT( S ), ISTACK( S ),
     &                              CTLEFF( S,1 ), FCO,
     &                              CTLEFF( S,2 ), FNOX,
     &                              CTLEFF( S,3 ), FVOC

                            CTLEFF( S,1 ) = FNOX
                            CTLEFF( S,2 ) = FVOC
                            CTLEFF( S,3 ) = FCO

123                     CONTINUE

                    ELSE IF ( STK .GT. 0 ) THEN   !  match plant,stack,fip,scc

                        S = FIND4( FIP, SCC, PLT, STK, 
     &                             NPSRC, IFIP, ISCC, IPLANT, ISTACK )

                        IF ( S .GT. 0 ) THEN

                            FNOX = ( CNOX - NOX( S ) ) /
     &                             ( RULPEN( S,1 ) 
     &                             * RULEFF( S,1 ) 
     &                             * CNOX )

                            FVOC = ( CVOC - VOC( S ) ) /
     &                             ( RULPEN( S,2 ) 
     &                             * RULEFF( S,2 ) 
     &                             * CVOC )

                            FCO  = ( CCO  - CO ( S ) ) /
     &                             ( RULPEN( S,3 ) 
     &                             * RULEFF( S,3 ) 
     &                             * CCO  )

                            IF ( TFLAG( S ) ) 
     &                          WRITE( RDEV,93010 ) 
     &                              PNAME, IFIP( S ), ISCC( S ), 
     &                              IPLANT( S ), ISTACK( S ),
     &                              CTLEFF( S,1 ), FCO,
     &                              CTLEFF( S,2 ), FNOX,
     &                              CTLEFF( S,3 ), FVOC

                            CTLEFF( S,1 ) = FNOX
                            CTLEFF( S,2 ) = FVOC
                            CTLEFF( S,3 ) = FCO

                        END IF  ! if s > 0

                    ELSE IF ( PLT .GT. 0 ) THEN !  match plant,fip,scc,sic

                        I = FIND1( FIP, NFIP, FFIP )
                        IF ( I .EQ. 0 ) GO TO 121

                        DO  124 S = JFIP( I ), JFIP( I+1 ) -1

                            IF( PLT .NE. IPLANT( S ) ) GO TO 124  ! n/a control
                            IF( SCC .NE. ISCC  ( S ) ) GO TO 124  ! n/a control
                            IF( SIC .NE. ISIC( S ) .AND.
     &                          SIC .NE. ISI3( S ) .AND.
     &                          SIC .NE. ISI2( S ) ) GO TO 124

                            FNOX = ( CNOX - NOX( S ) ) /
     &                             ( RULPEN( S,1 ) 
     &                             * RULEFF( S,1 ) 
     &                             * CNOX )

                            FVOC = ( CVOC - VOC( S ) ) /
     &                             ( RULPEN( S,2 ) 
     &                             * RULEFF( S,2 ) 
     &                             * CVOC )

                            FCO  = ( CCO  - CO ( S ) ) /
     &                             ( RULPEN( S,3 ) 
     &                             * RULEFF( S,3 ) 
     &                             * CCO  )

                            IF ( TFLAG( S ) ) 
     &                          WRITE( RDEV,93010 ) 
     &                              PNAME, IFIP( S ), ISCC( S ), 
     &                              IPLANT( S ), ISTACK( S ),
     &                              CTLEFF( S,1 ), FCO,
     &                              CTLEFF( S,2 ), FNOX,
     &                              CTLEFF( S,3 ), FVOC

                            CTLEFF( S,1 ) = FNOX
                            CTLEFF( S,2 ) = FVOC
                            CTLEFF( S,3 ) = FCO

124                     CONTINUE

                    ELSE        !  plant=o:  all matches for this fip,scc,sic

                        I = FIND1( FIP, NFIP, FFIP )
                        IF ( I .EQ. 0 ) GO TO 121

                        DO  125 S = JFIP( I ), JFIP( I+1 ) -1

                            IF( SCC .NE. ISCC  ( S ) ) GO TO 125  ! n/a control
                            IF( SIC .NE. ISIC( S ) .AND.
     &                          SIC .NE. ISI3( S ) .AND.
     &                          SIC .NE. ISI2( S ) ) GO TO 125

                            FNOX = ( CNOX - NOX( S ) ) /
     &                             ( RULPEN( S,1 ) 
     &                             * RULEFF( S,1 ) 
     &                             * CNOX )

                            FVOC = ( CVOC - VOC( S ) ) /
     &                             ( RULPEN( S,2 ) 
     &                             * RULEFF( S,2 ) 
     &                             * CVOC )

                            FCO  = ( CCO  - CO ( S ) ) /
     &                             ( RULPEN( S,3 ) 
     &                             * RULEFF( S,3 ) 
     &                             * CCO  )

                            IF ( TFLAG( S ) ) 
     &                          WRITE( RDEV,93010 ) 
     &                              PNAME, IFIP( S ), ISCC( S ), 
     &                              IPLANT( S ), ISTACK( S ),
     &                              CTLEFF( S,1 ), FCO,
     &                              CTLEFF( S,2 ), FNOX,
     &                              CTLEFF( S,3 ), FVOC

                            CTLEFF( S,1 ) = FNOX
                            CTLEFF( S,2 ) = FVOC
                            CTLEFF( S,3 ) = FCO

125                     CONTINUE

                    END IF      ! if fip is 0, or nn000, or not.

                ELSE		!  sic missing (where scc present)

                    IF ( FIP .LE. 0 ) THEN      !  all states.

                        DO  126 S = 1, NPSRC

                            IF( SCC .NE. ISCC( S ) ) GO TO 126  ! n/a control

                            FNOX = ( CNOX - NOX( S ) ) /
     &                             ( RULPEN( S,1 ) 
     &                             * RULEFF( S,1 ) 
     &                             * CNOX )

                            FVOC = ( CVOC - VOC( S ) ) /
     &                             ( RULPEN( S,2 ) 
     &                             * RULEFF( S,2 ) 
     &                             * CVOC )

                            FCO  = ( CCO  - CO ( S ) ) /
     &                             ( RULPEN( S,3 ) 
     &                             * RULEFF( S,3 ) 
     &                             * CCO  )

                            IF ( TFLAG( S ) ) 
     &                          WRITE( RDEV,93010 ) 
     &                              PNAME, IFIP( S ), ISCC( S ), 
     &                              IPLANT( S ), ISTACK( S ),
     &                              CTLEFF( S,1 ), FCO,
     &                              CTLEFF( S,2 ), FNOX,
     &                              CTLEFF( S,3 ), FVOC

                            CTLEFF( S,1 ) = FNOX
                            CTLEFF( S,2 ) = FVOC
                            CTLEFF( S,3 ) = FCO

126                     CONTINUE

                    ELSE IF ( MOD( FIP, 1000 ) .EQ. 0 ) THEN

                        SID = FIP / 1000        !  state code only
                        I   = FIND1( SID, NSID, SIDS )
                        IF ( I .EQ. 0 ) GO TO 121

                        DO  127 S = JSID( I ), JSID( I+1 ) -1

                            IF( SCC .NE. ISCC ( S ) ) GO TO 127  ! n/a control

                            FNOX = ( CNOX - NOX( S ) ) /
     &                             ( RULPEN( S,1 ) 
     &                             * RULEFF( S,1 ) 
     &                             * CNOX )

                            FVOC = ( CVOC - VOC( S ) ) /
     &                             ( RULPEN( S,2 ) 
     &                             * RULEFF( S,2 ) 
     &                             * CVOC )

                            FCO  = ( CCO  - CO ( S ) ) /
     &                             ( RULPEN( S,3 ) 
     &                             * RULEFF( S,3 ) 
     &                             * CCO  )

                            IF ( TFLAG( S ) ) 
     &                          WRITE( RDEV,93010 ) 
     &                              PNAME, IFIP( S ), ISCC( S ), 
     &                              IPLANT( S ), ISTACK( S ),
     &                              CTLEFF( S,1 ), FCO,
     &                              CTLEFF( S,2 ), FNOX,
     &                              CTLEFF( S,3 ), FVOC

                            CTLEFF( S,1 ) = FNOX
                            CTLEFF( S,2 ) = FVOC
                            CTLEFF( S,3 ) = FCO

127                     CONTINUE

                    ELSE IF ( STK .GT. 0 ) THEN         !  match plant,stack,fip,scc

                        S = FIND4( FIP, SCC, PLT, STK, 
     &                             NPSRC, IFIP, ISCC, IPLANT, ISTACK )

                        IF ( S .GT. 0 ) THEN

                            FNOX = ( CNOX - NOX( S ) ) /
     &                             ( RULPEN( S,1 ) 
     &                             * RULEFF( S,1 ) 
     &                             * CNOX )

                            FVOC = ( CVOC - VOC( S ) ) /
     &                             ( RULPEN( S,2 ) 
     &                             * RULEFF( S,2 ) 
     &                             * CVOC )

                            FCO  = ( CCO  - CO ( S ) ) /
     &                             ( RULPEN( S,3 ) 
     &                             * RULEFF( S,3 ) 
     &                             * CCO  )

                            IF ( TFLAG( S ) ) 
     &                          WRITE( RDEV,93010 ) 
     &                              PNAME, IFIP( S ), ISCC( S ), 
     &                              IPLANT( S ), ISTACK( S ),
     &                              CTLEFF( S,1 ), FCO,
     &                              CTLEFF( S,2 ), FNOX,
     &                              CTLEFF( S,3 ), FVOC

                            CTLEFF( S,1 ) = FNOX
                            CTLEFF( S,2 ) = FVOC
                            CTLEFF( S,3 ) = FCO

                        END IF  ! if s > 0

                    ELSE IF ( PLT .GT. 0 ) THEN !  match plant,fip,scc,sic

                        I = FIND1( FIP, NFIP, FFIP )
                        IF ( I .EQ. 0 ) GO TO 121

                        DO  128 S = JFIP( I ), JFIP( I+1 ) -1

                            IF( PLT .NE. IPLANT( S ) ) GO TO 128  ! n/a control
                            IF( SCC .NE. ISCC  ( S ) ) GO TO 128  ! n/a control

                            FNOX = ( CNOX - NOX( S ) ) /
     &                             ( RULPEN( S,1 ) 
     &                             * RULEFF( S,1 ) 
     &                             * CNOX )

                            FVOC = ( CVOC - VOC( S ) ) /
     &                             ( RULPEN( S,2 ) 
     &                             * RULEFF( S,2 ) 
     &                             * CVOC )

                            FCO  = ( CCO  - CO ( S ) ) /
     &                             ( RULPEN( S,3 ) 
     &                             * RULEFF( S,3 ) 
     &                             * CCO  )

                            IF ( TFLAG( S ) ) 
     &                          WRITE( RDEV,93010 ) 
     &                              PNAME, IFIP( S ), ISCC( S ), 
     &                              IPLANT( S ), ISTACK( S ),
     &                              CTLEFF( S,1 ), FCO,
     &                              CTLEFF( S,2 ), FNOX,
     &                              CTLEFF( S,3 ), FVOC

                            CTLEFF( S,1 ) = FNOX
                            CTLEFF( S,2 ) = FVOC
                            CTLEFF( S,3 ) = FCO

128                     CONTINUE

                    ELSE        !  plant=o:  all matches for this fip,scc,sic

                        I = FIND1( FIP, NFIP, FFIP )
                        IF ( I .EQ. 0 ) GO TO 121

                        DO  129 S = JFIP( I ), JFIP( I+1 ) -1

                            IF( SCC .NE. ISCC  ( S ) ) GO TO 129  ! n/a control
   
                            FNOX = ( CNOX - NOX( S ) ) /
     &                             ( RULPEN( S,1 ) 
     &                             * RULEFF( S,1 ) 
     &                             * CNOX )

                            FVOC = ( CVOC - VOC( S ) ) /
     &                             ( RULPEN( S,2 ) 
     &                             * RULEFF( S,2 ) 
     &                             * CVOC )

                            FCO  = ( CCO  - CO ( S ) ) /
     &                             ( RULPEN( S,3 ) 
     &                             * RULEFF( S,3 ) 
     &                             * CCO  )

                            IF ( TFLAG( S ) ) 
     &                          WRITE( RDEV,93010 ) 
     &                              PNAME, IFIP( S ), ISCC( S ), 
     &                              IPLANT( S ), ISTACK( S ),
     &                              CTLEFF( S,1 ), FCO,
     &                              CTLEFF( S,2 ), FNOX,
     &                              CTLEFF( S,3 ), FVOC

                            CTLEFF( S,1 ) = FNOX
                            CTLEFF( S,2 ) = FVOC
                            CTLEFF( S,3 ) = FCO

129                     CONTINUE

                    END IF      ! if fip is 0, or nn000, or not.

                END IF		!  sic > 0, or not (where scc present)

            ELSE		!  scc missing (where scc present)

                IF ( SIC .GT. 0 ) THEN

                    IF ( FIP .LE. 0 ) THEN      !  all states.

                        DO  130 S = 1, NPSRC

                            IF( SIC .NE. ISIC( S ) .AND.
     &                          SIC .NE. ISI3( S ) .AND.
     &                          SIC .NE. ISI2( S ) ) GO TO 130

                            FNOX = ( CNOX - NOX( S ) ) /
     &                             ( RULPEN( S,1 ) 
     &                             * RULEFF( S,1 ) 
     &                             * CNOX )

                            FVOC = ( CVOC - VOC( S ) ) /
     &                             ( RULPEN( S,2 ) 
     &                             * RULEFF( S,2 ) 
     &                             * CVOC )

                            FCO  = ( CCO  - CO ( S ) ) /
     &                             ( RULPEN( S,3 ) 
     &                             * RULEFF( S,3 ) 
     &                             * CCO  )

                            IF ( TFLAG( S ) ) 
     &                          WRITE( RDEV,93010 ) 
     &                              PNAME, IFIP( S ), ISCC( S ), 
     &                              IPLANT( S ), ISTACK( S ),
     &                              CTLEFF( S,1 ), FCO,
     &                              CTLEFF( S,2 ), FNOX,
     &                              CTLEFF( S,3 ), FVOC

                            CTLEFF( S,1 ) = FNOX
                            CTLEFF( S,2 ) = FVOC
                            CTLEFF( S,3 ) = FCO

130                     CONTINUE

                    ELSE IF ( MOD( FIP, 1000 ) .EQ. 0 ) THEN

                        SID = FIP / 1000        !  state code only
                        I   = FIND1( SID, NSID, SIDS )
                        IF ( I .EQ. 0 ) GO TO 121

                        DO  131 S = JSID( I ), JSID( I+1 ) -1

                            IF( SIC .NE. ISIC( S ) .AND.
     &                          SIC .NE. ISI3( S ) .AND.
     &                          SIC .NE. ISI2( S ) ) GO TO 131

                            FNOX = ( CNOX - NOX( S ) ) /
     &                             ( RULPEN( S,1 ) 
     &                             * RULEFF( S,1 ) 
     &                             * CNOX )

                            FVOC = ( CVOC - VOC( S ) ) /
     &                             ( RULPEN( S,2 ) 
     &                             * RULEFF( S,2 ) 
     &                             * CVOC )

                            FCO  = ( CCO  - CO ( S ) ) /
     &                             ( RULPEN( S,3 ) 
     &                             * RULEFF( S,3 ) 
     &                             * CCO  )

                            IF ( TFLAG( S ) ) 
     &                          WRITE( RDEV,93010 ) 
     &                              PNAME, IFIP( S ), ISCC( S ), 
     &                              IPLANT( S ), ISTACK( S ),
     &                              CTLEFF( S,1 ), FCO,
     &                              CTLEFF( S,2 ), FNOX,
     &                              CTLEFF( S,3 ), FVOC

                            CTLEFF( S,1 ) = FNOX
                            CTLEFF( S,2 ) = FVOC
                            CTLEFF( S,3 ) = FCO

131                     CONTINUE

                    ELSE IF ( STK .GT. 0 ) THEN !  match plant,stack,fip

                        I = FIND1( FIP, NFIP, FFIP )
                        IF ( I .EQ. 0 ) GO TO 121

                        DO  132 S = JFIP( I ), JFIP( I+1 ) - 1

                            IF( SIC .NE. ISIC  ( S ) .AND.
     &                          SIC .NE. ISI3  ( S ) .AND.
     &                          SIC .NE. ISI2  ( S ) ) GO TO 132
                            IF( PLT .NE. IPLANT( S ) ) GO TO 132  ! n/a control
                            IF( STK .NE. ISTACK( S ) ) GO TO 132  ! n/a control

                            FNOX = ( CNOX - NOX( S ) ) /
     &                             ( RULPEN( S,1 ) 
     &                             * RULEFF( S,1 ) 
     &                             * CNOX )

                            FVOC = ( CVOC - VOC( S ) ) /
     &                             ( RULPEN( S,2 ) 
     &                             * RULEFF( S,2 ) 
     &                             * CVOC )

                            FCO  = ( CCO  - CO ( S ) ) /
     &                             ( RULPEN( S,3 ) 
     &                             * RULEFF( S,3 ) 
     &                             * CCO  )

                            IF ( TFLAG( S ) ) 
     &                          WRITE( RDEV,93010 ) 
     &                              PNAME, IFIP( S ), ISCC( S ), 
     &                              IPLANT( S ), ISTACK( S ),
     &                              CTLEFF( S,1 ), FCO,
     &                              CTLEFF( S,2 ), FNOX,
     &                              CTLEFF( S,3 ), FVOC

                            CTLEFF( S,1 ) = FNOX
                            CTLEFF( S,2 ) = FVOC
                            CTLEFF( S,3 ) = FCO

132                     CONTINUE

                    ELSE IF ( PLT .GT. 0 ) THEN !  match plant,fip,scc,sic

                        I = FIND1( FIP, NFIP, FFIP )
                        IF ( I .EQ. 0 ) GO TO 121

                        DO  133 S = JFIP( I ), JFIP( I+1 ) -1

                            IF( PLT .NE. IPLANT( S ) ) GO TO 133  ! n/a control
                            IF( SIC .NE. ISIC  ( S ) .AND.
     &                          SIC .NE. ISI3  ( S ) .AND.
     &                          SIC .NE. ISI2  ( S ) ) GO TO 133

                            FNOX = ( CNOX - NOX( S ) ) /
     &                             ( RULPEN( S,1 ) 
     &                             * RULEFF( S,1 ) 
     &                             * CNOX )

                            FVOC = ( CVOC - VOC( S ) ) /
     &                             ( RULPEN( S,2 ) 
     &                             * RULEFF( S,2 ) 
     &                             * CVOC )

                            FCO  = ( CCO  - CO ( S ) ) /
     &                             ( RULPEN( S,3 ) 
     &                             * RULEFF( S,3 ) 
     &                             * CCO  )

                            IF ( TFLAG( S ) ) 
     &                          WRITE( RDEV,93010 ) 
     &                              PNAME, IFIP( S ), ISCC( S ), 
     &                              IPLANT( S ), ISTACK( S ),
     &                              CTLEFF( S,1 ), FCO,
     &                              CTLEFF( S,2 ), FNOX,
     &                              CTLEFF( S,3 ), FVOC

                            CTLEFF( S,1 ) = FNOX
                            CTLEFF( S,2 ) = FVOC
                            CTLEFF( S,3 ) = FCO

133                     CONTINUE

                    ELSE        !  plant=o:  all matches for this fip,sic

                        I = FIND1( FIP, NFIP, FFIP )
                        IF ( I .EQ. 0 ) GO TO 121

                        DO  134 S = JFIP( I ), JFIP( I+1 ) -1

                            IF( SIC .NE. ISIC( S ) .AND.
     &                          SIC .NE. ISI3( S ) .AND.
     &                          SIC .NE. ISI2( S ) ) GO TO 134

                            FNOX = ( CNOX - NOX( S ) ) /
     &                             ( RULPEN( S,1 ) 
     &                             * RULEFF( S,1 ) 
     &                             * CNOX )

                            FVOC = ( CVOC - VOC( S ) ) /
     &                             ( RULPEN( S,2 ) 
     &                             * RULEFF( S,2 ) 
     &                             * CVOC )

                            FCO  = ( CCO  - CO ( S ) ) /
     &                             ( RULPEN( S,3 ) 
     &                             * RULEFF( S,3 ) 
     &                             * CCO  )

                            IF ( TFLAG( S ) ) 
     &                          WRITE( RDEV,93010 ) 
     &                              PNAME, IFIP( S ), ISCC( S ), 
     &                              IPLANT( S ), ISTACK( S ),
     &                              CTLEFF( S,1 ), FCO,
     &                              CTLEFF( S,2 ), FNOX,
     &                              CTLEFF( S,3 ), FVOC

                            CTLEFF( S,1 ) = FNOX
                            CTLEFF( S,2 ) = FVOC
                            CTLEFF( S,3 ) = FCO

134                     CONTINUE

                    END IF      ! if fip is 0, or nn000, or not.

                ELSE		!  sic missing (where scc missing)

                    IF ( FIP .LE. 0 ) THEN      !  all states.

                        DO  135 S = 1, NPSRC

                            FNOX = ( CNOX - NOX( S ) ) /
     &                             ( RULPEN( S,1 ) 
     &                             * RULEFF( S,1 ) 
     &                             * CNOX )

                            FVOC = ( CVOC - VOC( S ) ) /
     &                             ( RULPEN( S,2 ) 
     &                             * RULEFF( S,2 ) 
     &                             * CVOC )

                            FCO  = ( CCO  - CO ( S ) ) /
     &                             ( RULPEN( S,3 ) 
     &                             * RULEFF( S,3 ) 
     &                             * CCO  )

                            IF ( TFLAG( S ) ) 
     &                          WRITE( RDEV,93010 ) 
     &                              PNAME, IFIP( S ), ISCC( S ), 
     &                              IPLANT( S ), ISTACK( S ),
     &                              CTLEFF( S,1 ), FCO,
     &                              CTLEFF( S,2 ), FNOX,
     &                              CTLEFF( S,3 ), FVOC

                            CTLEFF( S,1 ) = FNOX
                            CTLEFF( S,2 ) = FVOC
                            CTLEFF( S,3 ) = FCO

135                     CONTINUE

                    ELSE IF ( MOD( FIP, 1000 ) .EQ. 0 ) THEN

                        SID = FIP / 1000        !  state code only
                        I   = FIND1( SID, NSID, SIDS )
                        IF ( I .EQ. 0 ) GO TO 121

                        DO  136 S = JSID( I ), JSID( I+1 ) -1

                            FNOX = ( CNOX - NOX( S ) ) /
     &                             ( RULPEN( S,1 ) 
     &                             * RULEFF( S,1 ) 
     &                             * CNOX )

                            FVOC = ( CVOC - VOC( S ) ) /
     &                             ( RULPEN( S,2 ) 
     &                             * RULEFF( S,2 ) 
     &                             * CVOC )

                            FCO  = ( CCO  - CO ( S ) ) /
     &                             ( RULPEN( S,3 ) 
     &                             * RULEFF( S,3 ) 
     &                             * CCO  )

                            IF ( TFLAG( S ) ) 
     &                          WRITE( RDEV,93010 ) 
     &                              PNAME, IFIP( S ), ISCC( S ), 
     &                              IPLANT( S ), ISTACK( S ),
     &                              CTLEFF( S,1 ), FCO,
     &                              CTLEFF( S,2 ), FNOX,
     &                              CTLEFF( S,3 ), FVOC

                            CTLEFF( S,1 ) = FNOX
                            CTLEFF( S,2 ) = FVOC
                            CTLEFF( S,3 ) = FCO

136                     CONTINUE

                    ELSE IF ( STK .GT. 0 ) THEN !  match plant,stack,fip

                        I   = FIND1( FIP, NFIP, FFIP )
                        IF ( I .EQ. 0 ) GO TO 121

                        DO  137 S = JFIP( I ), JFIP( I+1 ) - 1

                            IF( PLT .NE. IPLANT( S ) ) GO TO 137  ! n/a control
                            IF( STK .NE. ISTACK( S ) ) GO TO 137  ! n/a control

                            FNOX = ( CNOX - NOX( S ) ) /
     &                             ( RULPEN( S,1 ) 
     &                             * RULEFF( S,1 ) 
     &                             * CNOX )

                            FVOC = ( CVOC - VOC( S ) ) /
     &                             ( RULPEN( S,2 ) 
     &                             * RULEFF( S,2 ) 
     &                             * CVOC )

                            FCO  = ( CCO  - CO ( S ) ) /
     &                             ( RULPEN( S,3 ) 
     &                             * RULEFF( S,3 ) 
     &                             * CCO  )

                            IF ( TFLAG( S ) ) 
     &                          WRITE( RDEV,93010 ) 
     &                              PNAME, IFIP( S ), ISCC( S ), 
     &                              IPLANT( S ), ISTACK( S ),
     &                              CTLEFF( S,1 ), FCO,
     &                              CTLEFF( S,2 ), FNOX,
     &                              CTLEFF( S,3 ), FVOC

                            CTLEFF( S,1 ) = FNOX
                            CTLEFF( S,2 ) = FVOC
                            CTLEFF( S,3 ) = FCO

137                     CONTINUE

                    ELSE IF ( PLT .GT. 0 ) THEN !  match plant,fip,scc,sic

                        I = FIND1( FIP, NFIP, FFIP )
                        IF ( I .EQ. 0 ) GO TO 121

                        DO  138 S = JFIP( I ), JFIP( I+1 ) -1

                            IF( PLT .NE. IPLANT( S ) ) GO TO 138  ! n/a control

                            FNOX = ( CNOX - NOX( S ) ) /
     &                             ( RULPEN( S,1 ) 
     &                             * RULEFF( S,1 ) 
     &                             * CNOX )

                            FVOC = ( CVOC - VOC( S ) ) /
     &                             ( RULPEN( S,2 ) 
     &                             * RULEFF( S,2 ) 
     &                             * CVOC )

                            FCO  = ( CCO  - CO ( S ) ) /
     &                             ( RULPEN( S,3 ) 
     &                             * RULEFF( S,3 ) 
     &                             * CCO  )

                            IF ( TFLAG( S ) ) 
     &                          WRITE( RDEV,93010 ) 
     &                              PNAME, IFIP( S ), ISCC( S ), 
     &                              IPLANT( S ), ISTACK( S ),
     &                              CTLEFF( S,1 ), FCO,
     &                              CTLEFF( S,2 ), FNOX,
     &                              CTLEFF( S,3 ), FVOC

                            CTLEFF( S,1 ) = FNOX
                            CTLEFF( S,2 ) = FVOC
                            CTLEFF( S,3 ) = FCO

138                     CONTINUE

                    ELSE        !  plant=o:  all matches for this fip,scc,sic

                        I = FIND1( FIP, NFIP, FFIP )
                        IF ( I .EQ. 0 ) GO TO 121

                        DO  139 S = JFIP( I ), JFIP( I+1 ) -1

                            FNOX = ( CNOX - NOX( S ) ) /
     &                             ( RULPEN( S,1 ) 
     &                             * RULEFF( S,1 ) 
     &                             * CNOX )

                            FVOC = ( CVOC - VOC( S ) ) /
     &                             ( RULPEN( S,2 ) 
     &                             * RULEFF( S,2 ) 
     &                             * CVOC )

                            FCO  = ( CCO  - CO ( S ) ) /
     &                             ( RULPEN( S,3 ) 
     &                             * RULEFF( S,3 ) 
     &                             * CCO  )

                            IF ( TFLAG( S ) ) 
     &                          WRITE( RDEV,93010 ) 
     &                              PNAME, IFIP( S ), ISCC( S ), 
     &                              IPLANT( S ), ISTACK( S ),
     &                              CTLEFF( S,1 ), FCO,
     &                              CTLEFF( S,2 ), FNOX,
     &                              CTLEFF( S,3 ), FVOC

                            CTLEFF( S,1 ) = FNOX
                            CTLEFF( S,2 ) = FVOC
                            CTLEFF( S,3 ) = FCO

139                     CONTINUE

                    END IF      ! if fip is 0, or nn000, or not.

                END IF		!  sic > 0, or not (where scc missing)

            END IF		!  scc > 0, or not where scc missing)

            GO TO  121  !  end loop on control-lines in this RULE PENETRATION packet

140     CONTINUE        !  exit from / skip around loop 121


C............................................................................
C.......   Process PROJECT PTS packet of control file:

        PNAME = '/PROJECT PTS/'
        IF( 0 .EQ. INDEX1( PNAME, NPACKETS, PNAMES ) ) GO TO 149

        IF ( .NOT. GETPKT( PNAME, CDEV ) ) THEN
            CALL M3EXIT( 'CTLPMAT', 0, 0,
     &                   'Could not find ' // PNAME // ' in ' // CNAME, 
     &                   2 )
        END IF

C.........  Prompt for type of projection entries: default (EPS) or 
C           modified EPS (with year)
        PFLAG = GETYN(
     &        'Are projection entries in year-specific format?', .TRUE.)

141     CONTINUE        !  head of loop reading lines in the PROJECT PTS packet

            READ( CDEV,93000,IOSTAT=IOS,END=900 ) LINE
            IREC = IREC + 1
            IF( IOS .NE. 0 ) THEN
                WRITE( MESG,94920 ) 
     &              'Error reading', CNAME( 1:TRIMLEN( CNAME ) ),
     &              IREC, IOS
                CALL M3EXIT( 'CTLPMAT', 0, 0, MESG, 1 )
            END IF

            CALL UPCASE( LINE )
            IF( 0 .NE. INDEX( LINE, 'END' ) ) GO TO 149

            IF( PFLAG ) THEN
                YEAR = STR2INT ( LINE(  1: 4 ) )
                DYR  = STR2INT ( LINE(  6: 9 ) )
                FIP  = STR2INT ( LINE( 11:15 ) )
                SIC  = STR2INT ( LINE( 26:35 ) )
                PFAC = STR2REAL( LINE( 37:46 ) )

            ELSE
                FIP  = STR2INT ( LINE(  1: 5 ) )
                SIC  = STR2INT ( LINE( 16:25 ) )
                PFAC = STR2REAL( LINE( 27:36 ) )

            ENDIF

            IF ( FIP .LE. 0 ) THEN      !  all states.

                DO  142 S = 1, NPSRC

                    IF( SIC .NE. ISIC( S ) .AND.
     &                  SIC .NE. ISI3( S ) .AND.
     &                  SIC .NE. ISI2( S ) ) GO TO 142

                    IF( PFLAG .AND. INVYR( S ) .NE. YEAR ) GO TO 142

                    IF ( TFLAG( S ) ) WRITE( RDEV,93020 ) 
     &                     PNAME, YEAR, DYR, FIP, SIC, PRJFAC( S ), PFAC

                    PRJFAC( S ) = PFAC

142             CONTINUE       !  end loop on S

            ELSE IF ( MOD( FIP, 1000 ) .EQ. 0 ) THEN

                SID = FIP / 1000        !  state code only
                I   = FIND1( SID, NSID, SIDS )

                IF ( I .EQ. 0 ) GO TO 141

                DO  143 S = JSID( I ), JSID( I+1 ) - 1

                    IF( SIC .NE. ISIC( S ) .AND.
     &                  SIC .NE. ISI3( S ) .AND.
     &                  SIC .NE. ISI2( S ) ) GO TO 143  ! n/a control

                    IF( PFLAG .AND. INVYR( S ) .NE. YEAR ) GO TO 143

                    IF ( TFLAG( S ) ) WRITE( RDEV,93020 ) 
     &                PNAME, YEAR, DYR, SID*1000, SIC, PRJFAC( S ), PFAC

                    PRJFAC( S ) = PFAC

143             CONTINUE       !  end loop on S

            ELSE        !  fip-specific control

                I = FIND1( FIP, NFIP, FFIP )

                DO  144 S = JFIP( I ), JFIP( I+1 ) - 1

                    IF( SIC .NE. ISIC( S ) .AND.
     &                  SIC .NE. ISI3( S ) .AND.
     &                  SIC .NE. ISI2( S ) ) GO TO 144  ! n/a control

                    IF( PFLAG .AND. INVYR( S ) .NE. YEAR ) GO TO 144

                    IF ( TFLAG( S ) ) WRITE( RDEV,93020 )
     &                     PNAME, YEAR, DYR, FIP, SIC, PRJFAC( S ), PFAC

                    PRJFAC( S ) = PFAC

144             CONTINUE       !  end loop on S

            END IF      ! if fip is 0, or nn000, or not.

            GO TO  141  !  end loop reading control-line in PROJECT PTS packet

149     CONTINUE        !  exit from / skip around loop 141


C............................................................................
C.......   Process RULE PENETRATION packet of control file:

        PNAME = '/RULE PENETRATION/'
        IF( 0 .EQ. INDEX1( PNAME, NPACKETS, PNAMES ) ) GO TO 159

        IF ( .NOT. GETPKT( PNAME, CDEV ) ) THEN
            CALL M3EXIT( 'CTLPMAT', 0, 0,
     &                   'Could not find ' // PNAME // ' in ' // CNAME, 
     &                   2 )
        END IF

151     CONTINUE        !  head of loop reading lines in the RULPEN packet

            READ( CDEV,93000,IOSTAT=IOS,END=900 ) LINE
            IREC = IREC + 1
            IF( IOS .NE. 0 ) THEN
                WRITE( MESG,94920 ) 
     &              'Error reading', CNAME( 1:TRIMLEN( CNAME ) ), 
     &              IREC, IOS
                CALL M3EXIT( 'CTLPMAT', 0, 0, MESG, 1 )
            END IF

            CALL UPCASE( LINE )
            IF( 0 .NE. INDEX( LINE, 'END' ) )  GO TO 159

            FIP  = STR2INT ( LINE(  1: 5 ) )
            PLT  = STR2INT ( LINE( 16:20 ) )
            STK  = STR2INT ( LINE( 22:24 ) )
            SIC  = STR2INT ( LINE( 26:29 ) )
            SCC  = STR2INT ( LINE( 31:40 ) )
            CNOX = 0.01 * STR2REAL( LINE( 46:50 ) )
            CVOC = 0.01 * STR2REAL( LINE( 51:55 ) )
            CCO  = 0.01 * STR2REAL( LINE( 56:60 ) )

            IF ( FIP .LE. 0 ) THEN      !  all states.
                
                DO  152 S = 1, NPSRC

                    IF( SCC .NE. ISCC( S ) ) GO TO 152  ! n/a control
                    IF( SIC .NE. ISIC( S ) .AND.
     &                  SIC .NE. ISI3( S ) .AND.
     &                  SIC .NE. ISI2( S ) ) GO TO 152

                    IF ( TFLAG( S ) ) 
     &                  WRITE( RDEV,93010 ) 
     &                      PNAME, IFIP( S ), ISCC( S ), 
     &                      IPLANT( S ), ISTACK( S ),
     &                      RULPEN( S,1 ), CCO,
     &                      RULPEN( S,2 ), CNOX,
     &                      RULPEN( S,3 ), CVOC

                    RULPEN( S,1 ) = CNOX
                    RULPEN( S,2 ) = CVOC
                    RULPEN( S,3 ) = CCO

152             CONTINUE

            ELSE IF ( MOD( FIP, 1000 ) .EQ. 0 ) THEN
                
                SID = FIP / 1000        !  state code only
                I   = FIND1( SID, NSID, SIDS )
                IF ( I .EQ. 0 ) GO TO 151

                DO  153 S = JSID( I ), JSID( I+1 ) -1

                    IF( SCC .NE. ISCC ( S ) ) GO TO 153  ! n/a control
                    IF( SIC .NE. ISIC( S ) .AND.
     &                  SIC .NE. ISI3( S ) .AND.
     &                  SIC .NE. ISI2( S ) ) GO TO 153

                    IF ( TFLAG( S ) ) 
     &                  WRITE( RDEV,93010 ) 
     &                      PNAME, IFIP( S ), ISCC( S ), 
     &                      IPLANT( S ), ISTACK( S ),
     &                      RULPEN( S,1 ), CCO,
     &                      RULPEN( S,2 ), CNOX,
     &                      RULPEN( S,3 ), CVOC

                    RULPEN( S,1 ) = CNOX
                    RULPEN( S,2 ) = CVOC
                    RULPEN( S,3 ) = CCO

153             CONTINUE

            ELSE IF ( STK .GT. 0 ) THEN         !  match plant,stack,fip,scc

                S = FIND4( FIP, SCC, PLT, STK, 
     &                     NPSRC, IFIP, ISCC, IPLANT, ISTACK )

                IF ( S .GT. 0 ) THEN

                    IF ( TFLAG( S ) ) 
     &                  WRITE( RDEV,93010 ) 
     &                      PNAME, IFIP( S ), ISCC( S ), 
     &                      IPLANT( S ), ISTACK( S ),
     &                      RULPEN( S,1 ), CCO,
     &                      RULPEN( S,2 ), CNOX,
     &                      RULPEN( S,3 ), CVOC

                    RULPEN( S,1 ) = CNOX
                    RULPEN( S,2 ) = CVOC
                    RULPEN( S,3 ) = CCO

                END IF  ! if s > 0

            ELSE IF ( PLT .GT. 0 ) THEN !  match plant,fip,scc,sic

                I = FIND1( FIP, NFIP, FFIP )
                IF ( I .EQ. 0 ) GO TO 151

                DO  154 S = JFIP( I ), JFIP( I+1 ) -1

                    IF( PLT .NE. IPLANT( S ) ) GO TO 154  ! n/a control
                    IF( SCC .NE. ISCC  ( S ) ) GO TO 154  ! n/a control
                    IF( SIC .NE. ISIC  ( S ) .AND.
     &                  SIC .NE. ISI3  ( S ) .AND.
     &                  SIC .NE. ISI2  ( S ) ) GO TO 154

                    IF ( TFLAG( S ) ) 
     &                  WRITE( RDEV,93010 ) 
     &                      PNAME, IFIP( S ), ISCC( S ), 
     &                      IPLANT( S ), ISTACK( S ),
     &                      RULPEN( S,1 ), CCO,
     &                      RULPEN( S,2 ), CNOX,
     &                      RULPEN( S,3 ), CVOC

                    RULPEN( S,1 ) = CNOX
                    RULPEN( S,2 ) = CVOC
                    RULPEN( S,3 ) = CCO

154             CONTINUE

            ELSE        !  plant=o:  all matches for this fip,scc,sic

                I = FIND1( FIP, NFIP, FFIP )
                IF ( I .EQ. 0 ) GO TO 151

                DO  155 S = JFIP( I ), JFIP( I+1 ) -1

                    IF( SCC .NE. ISCC  ( S ) ) GO TO 155  ! n/a control
                    IF( SIC .NE. ISIC( S ) .AND.
     &                  SIC .NE. ISI3( S ) .AND.
     &                  SIC .NE. ISI2( S ) ) GO TO 155

                    IF ( TFLAG( S ) ) 
     &                  WRITE( RDEV,93010 ) 
     &                      PNAME, IFIP( S ), ISCC( S ), 
     &                      IPLANT( S ), ISTACK( S ),
     &                      RULPEN( S,1 ), CCO,
     &                      RULPEN( S,2 ), CNOX,
     &                      RULPEN( S,3 ), CVOC

                    RULPEN( S,1 ) = CNOX
                    RULPEN( S,2 ) = CVOC
                    RULPEN( S,3 ) = CCO

155             CONTINUE

            END IF      ! if fip is 0, or nn000, or not.

            GO TO  151  !  end loop on control-lines in this RULE PENETRATION packet

159     CONTINUE        !  exit from / skip around loop 151


C............................................................................
C.......   Process RULE EFFECTIVENESS packet of control file:

        PNAME = '/RULE EFFECTIVENESS/'
        IF( 0 .EQ. INDEX1( PNAME, NPACKETS, PNAMES ) ) GO TO 169

        IF ( .NOT. GETPKT( PNAME, CDEV ) ) THEN
            CALL M3EXIT( 'CTLPMAT', 0, 0,
     &                   'Could not find ' // PNAME // ' in ' // CNAME, 
     &                   2 )
        END IF

161     CONTINUE        !  head of loop reading lines in the RULEFF packet

            READ( CDEV,93000,IOSTAT=IOS,END=900 ) LINE
            IREC = IREC + 1
            IF( IOS .NE. 0 ) THEN
                WRITE( MESG,94920 ) 
     &              'Error reading', CNAME( 1:TRIMLEN( CNAME ) ), 
     &              IREC, IOS
                CALL M3EXIT( 'CTLPMAT', 0, 0, MESG, 1 )
            END IF

            CALL UPCASE( LINE )
            IF( 0 .NE. INDEX( LINE, 'END' ) )  GO TO 169

            FIP  = STR2INT ( LINE(  1: 5 ) )
            PLT  = STR2INT ( LINE( 16:20 ) )
            STK  = STR2INT ( LINE( 22:24 ) )
            SIC  = STR2INT ( LINE( 26:29 ) )
            SCC  = STR2INT ( LINE( 31:40 ) )
            CNOX = 0.01 * STR2REAL( LINE( 46:50 ) )
            CVOC = 0.01 * STR2REAL( LINE( 51:55 ) )
            CCO  = 0.01 * STR2REAL( LINE( 56:60 ) )

            IF ( FIP .LE. 0 ) THEN      !  all states.
                
                DO  162 S = 1, NPSRC

                    IF( SCC .NE. ISCC( S ) ) GO TO 162  ! n/a control
                    IF( SIC .NE. ISIC( S ) .AND.
     &                  SIC .NE. ISI3( S ) .AND.
     &                  SIC .NE. ISI2( S ) ) GO TO 162

                    IF ( TFLAG( S ) ) 
     &                  WRITE( RDEV,93010 ) 
     &                      PNAME, IFIP( S ), ISCC( S ), 
     &                      IPLANT( S ), ISTACK( S ),
     &                      RULEFF( S,1 ), CCO,
     &                      RULEFF( S,2 ), CNOX,
     &                      RULEFF( S,3 ), CVOC

                    RULEFF( S,1 ) = CNOX
                    RULEFF( S,2 ) = CVOC
                    RULEFF( S,3 ) = CCO

162             CONTINUE

            ELSE IF ( MOD( FIP, 1000 ) .EQ. 0 ) THEN
                
                SID = FIP / 1000        !  state code only
                I   = FIND1( SID, NSID, SIDS )
                IF ( I .EQ. 0 ) GO TO 161

                DO  163 S = JSID( I ), JSID( I+1 ) -1

                    IF( SCC .NE. ISCC( S ) ) GO TO 163  ! n/a control
                    IF( SIC .NE. ISIC( S ) .AND.
     &                  SIC .NE. ISI3( S ) .AND.
     &                  SIC .NE. ISI2( S ) ) GO TO 163

                    IF ( TFLAG( S ) ) 
     &                  WRITE( RDEV,93010 ) 
     &                      PNAME, IFIP( S ), ISCC( S ), 
     &                      IPLANT( S ), ISTACK( S ),
     &                      RULEFF( S,1 ), CCO,
     &                      RULEFF( S,2 ), CNOX,
     &                      RULEFF( S,3 ), CVOC

                    RULEFF( S,1 ) = CNOX
                    RULEFF( S,2 ) = CVOC
                    RULEFF( S,3 ) = CCO

163             CONTINUE

            ELSE IF ( STK .GT. 0 ) THEN         !  match plant,stack,fip,scc

                S = FIND4( FIP, SCC, PLT, STK, 
     &                     NPSRC, IFIP, ISCC, IPLANT, ISTACK )

                IF ( S .GT. 0 ) THEN

                    IF ( TFLAG( S ) ) 
     &                  WRITE( RDEV,93010 ) 
     &                      PNAME, IFIP( S ), ISCC( S ), 
     &                      IPLANT( S ), ISTACK( S ),
     &                      RULEFF( S,1 ), CCO,
     &                      RULEFF( S,2 ), CNOX,
     &                      RULEFF( S,3 ), CVOC

                    RULEFF( S,1 ) = CNOX
                    RULEFF( S,2 ) = CVOC
                    RULEFF( S,3 ) = CCO

                END IF  ! if s > 0

            ELSE IF ( PLT .GT. 0 ) THEN !  match plant,fip,scc,sic

                I = FIND1( FIP, NFIP, FFIP )
                IF ( I .EQ. 0 ) GO TO 161

                DO  164 S = JFIP( I ), JFIP( I+1 ) -1

                    IF( PLT .NE. IPLANT( S ) ) GO TO 164  ! n/a control
                    IF( SCC .NE. ISCC  ( S ) ) GO TO 164  ! n/a control
                    IF( SIC .NE. ISIC  ( S ) .AND.
     &                  SIC .NE. ISI3  ( S ) .AND.
     &                  SIC .NE. ISI2  ( S ) ) GO TO 164

                    IF ( TFLAG( S ) ) 
     &                  WRITE( RDEV,93010 ) 
     &                      PNAME, IFIP( S ), ISCC( S ), 
     &                      IPLANT( S ), ISTACK( S ),
     &                      RULEFF( S,1 ), CCO,
     &                      RULEFF( S,2 ), CNOX,
     &                      RULEFF( S,3 ), CVOC

                    RULEFF( S,1 ) = CNOX
                    RULEFF( S,2 ) = CVOC
                    RULEFF( S,3 ) = CCO

164             CONTINUE

            ELSE        !  plant=o:  all matches for this fip,scc,sic

                I = FIND1( FIP, NFIP, FFIP )
                IF ( I .EQ. 0 ) GO TO 161

                DO  165 S = JFIP( I ), JFIP( I+1 ) -1

                    IF( SCC .NE. ISCC( S ) ) GO TO 165  ! n/a control
                    IF( SIC .NE. ISIC( S ) .AND.
     &                  SIC .NE. ISI3( S ) .AND.
     &                  SIC .NE. ISI2( S ) ) GO TO 165

                    IF ( TFLAG( S ) ) 
     &                  WRITE( RDEV,93010 ) 
     &                      PNAME, IFIP( S ), ISCC( S ), 
     &                      IPLANT( S ), ISTACK( S ),
     &                      RULEFF( S,1 ), CCO,
     &                      RULEFF( S,2 ), CNOX,
     &                      RULEFF( S,3 ), CVOC

                    RULEFF( S,1 ) = CNOX
                    RULEFF( S,2 ) = CVOC
                    RULEFF( S,3 ) = CCO

165             CONTINUE

            END IF      ! if fip is 0, or nn000, or not.

            GO TO  161  !  end loop on control-lines in this RULEFF packet

169     CONTINUE        !  exit from / skip around loop 161


C............................................................................
C.......   Process CONTROL EFFICIENCY packet of control file:

        PNAME = '/CONTROL EFFICIENCY/'
        IF( 0 .EQ. INDEX1( PNAME, NPACKETS, PNAMES ) ) GO TO 179

        IF ( .NOT. GETPKT( PNAME, CDEV ) ) THEN
            CALL M3EXIT( 'CTLPMAT', 0, 0,
     &                   'Could not find ' // PNAME // ' in ' // CNAME, 
     &                   2 )
        END IF

171     CONTINUE        !  head of loop reading lines in the CTLEFF packet

            READ( CDEV,93000,IOSTAT=IOS,END=900 ) LINE
            IREC = IREC + 1
            IF( IOS .NE. 0 ) THEN
                WRITE( MESG,94920 ) 
     &              'Error reading', CNAME( 1:TRIMLEN( CNAME ) ), 
     &              IREC, IOS
                CALL M3EXIT( 'CTLPMAT', 0, 0, MESG, 1 )
            END IF

            CALL UPCASE( LINE )
            IF( 0 .NE. INDEX( LINE, 'END' ) )  GO TO 179

            FIP  = STR2INT ( LINE(  1: 5 ) )
            PLT  = STR2INT ( LINE( 16:20 ) )
            STK  = STR2INT ( LINE( 22:24 ) )
            SIC  = STR2INT ( LINE( 26:29 ) )
            SCC  = STR2INT ( LINE( 31:40 ) )
            CNOX = 0.01 * STR2REAL( LINE( 46:50 ) )
            CVOC = 0.01 * STR2REAL( LINE( 51:55 ) )
            CCO  = 0.01 * STR2REAL( LINE( 56:60 ) )

            IF ( FIP .LE. 0 ) THEN      !  all states.
                
                DO  172 S = 1, NPSRC

                    IF( SCC .NE. ISCC( S ) ) GO TO 172  ! n/a control
                    IF( SIC .NE. ISIC( S ) .AND.
     &                  SIC .NE. ISI3( S ) .AND.
     &                  SIC .NE. ISI2( S ) ) GO TO 172

                    IF ( TFLAG( S ) ) 
     &                  WRITE( RDEV,93010 ) 
     &                      PNAME, IFIP( S ), ISCC( S ), 
     &                      IPLANT( S ), ISTACK( S ),
     &                      CTLEFF( S,1 ), CCO,
     &                      CTLEFF( S,2 ), CNOX,
     &                      CTLEFF( S,3 ), CVOC

                    CTLEFF( S,1 ) = CCO
                    CTLEFF( S,2 ) = CNOX
                    CTLEFF( S,3 ) = CVOC

172             CONTINUE

            ELSE IF ( MOD( FIP, 1000 ) .EQ. 0 ) THEN
                
                SID = FIP / 1000        !  state code only
                I   = FIND1( SID, NSID, SIDS )
                IF ( I .EQ. 0 ) GO TO 171

                DO  173 S = JSID( I ), JSID( I+1 ) -1

                    IF( SCC .NE. ISCC( S ) ) GO TO 173  ! n/a control
                    IF( SIC .NE. ISIC( S ) .AND.
     &                  SIC .NE. ISI3( S ) .AND.
     &                  SIC .NE. ISI2( S ) ) GO TO 173

                    IF ( TFLAG( S ) ) 
     &                  WRITE( RDEV,93010 ) 
     &                      PNAME, IFIP( S ), ISCC( S ), 
     &                      IPLANT( S ), ISTACK( S ),
     &                      CTLEFF( S,1 ), CCO,
     &                      CTLEFF( S,2 ), CNOX,
     &                      CTLEFF( S,3 ), CVOC

                    CTLEFF( S,1 ) = CCO
                    CTLEFF( S,2 ) = CNOX
                    CTLEFF( S,3 ) = CVOC

173             CONTINUE

            ELSE IF ( STK .GT. 0 ) THEN         !  match plant,stack,fip,scc

                S = FIND4( FIP, SCC, PLT, STK, 
     &                     NPSRC, IFIP, ISCC, IPLANT, ISTACK )

                IF ( S .GT. 0 ) THEN

                    IF ( TFLAG( S ) ) 
     &                  WRITE( RDEV,93010 ) 
     &                      PNAME, IFIP( S ), ISCC( S ), 
     &                      IPLANT( S ), ISTACK( S ),
     &                      CTLEFF( S,1 ), CCO,
     &                      CTLEFF( S,2 ), CNOX,
     &                      CTLEFF( S,3 ), CVOC

                    CTLEFF( S,1 ) = CCO
                    CTLEFF( S,2 ) = CNOX
                    CTLEFF( S,3 ) = CVOC

                END IF  ! if s > 0

            ELSE IF ( PLT .GT. 0 ) THEN !  match plant,fip,scc,sic

                I = FIND1( FIP, NFIP, FFIP )
                IF ( I .EQ. 0 ) GO TO 171

                DO  174 S = JFIP( I ), JFIP( I+1 ) -1

                    IF( PLT .NE. IPLANT( S ) ) GO TO 174  ! n/a control
                    IF( SCC .NE. ISCC  ( S ) ) GO TO 174  ! n/a control
                    IF( SIC .NE. ISIC  ( S ) .AND.
     &                  SIC .NE. ISI3  ( S ) .AND.
     &                  SIC .NE. ISI2  ( S ) ) GO TO 174

                    IF ( TFLAG( S ) ) 
     &                  WRITE( RDEV,93010 ) 
     &                      PNAME, IFIP( S ), ISCC( S ), 
     &                      IPLANT( S ), ISTACK( S ),
     &                      CTLEFF( S,1 ), CCO,
     &                      CTLEFF( S,2 ), CNOX,
     &                      CTLEFF( S,3 ), CVOC

                    CTLEFF( S,1 ) = CCO
                    CTLEFF( S,2 ) = CNOX
                    CTLEFF( S,3 ) = CVOC

174             CONTINUE

            ELSE        !  plant=o:  all matches for this fip,scc,sic

                I = FIND1( FIP, NFIP, FFIP )
                IF ( I .EQ. 0 ) GO TO 171

                DO  175 S = JFIP( I ), JFIP( I+1 ) -1

                    IF( SCC .NE. ISCC( S ) ) GO TO 175  ! n/a control
                    IF( SIC .NE. ISIC( S ) .AND.
     &                  SIC .NE. ISI3( S ) .AND.
     &                  SIC .NE. ISI2( S ) ) GO TO 175

                    IF ( TFLAG( S ) ) 
     &                  WRITE( RDEV,93010 ) 
     &                      PNAME, IFIP( S ), ISCC( S ), 
     &                      IPLANT( S ), ISTACK( S ),
     &                      CTLEFF( S,1 ), CCO,
     &                      CTLEFF( S,2 ), CNOX,
     &                      CTLEFF( S,3 ), CVOC

                    CTLEFF( S,1 ) = CCO
                    CTLEFF( S,2 ) = CNOX
                    CTLEFF( S,3 ) = CVOC

175             CONTINUE

            END IF      ! if fip is 0, or nn000, or not.

            GO TO  171  !  end loop on control-lines in this CTLEFF packet

179     CONTINUE        !  exit from / skip around loop 171


C............................................................................
C.......   Process RACT packet of control file:

        PNAME = '/RACT/'
        IF( 0 .EQ. INDEX1( PNAME, NPACKETS, PNAMES ) ) GO TO 189

        IF ( .NOT. GETPKT( PNAME, CDEV ) ) THEN
            CALL M3EXIT( 'CTLPMAT', 0, 0,
     &                   'Could not find ' // PNAME // ' in ' // CNAME, 
     &                   2 )
        END IF

181     CONTINUE        !  head of loop reading lines in the RACT packet

            READ( CDEV,93000,IOSTAT=IOS,END=900 ) LINE
            IREC = IREC + 1
            IF( IOS .NE. 0 ) THEN
                WRITE( MESG,94920 ) 
     &              'Error reading', CNAME( 1:TRIMLEN( CNAME ) ), 
     &              IREC, IOS
                CALL M3EXIT( 'CTLPMAT', 0, 0, MESG, 1 )
            END IF

            CALL UPCASE( LINE )
            IF( 0 .NE. INDEX( LINE, 'END' ) )  GO TO 189

            FIP  = STR2INT ( LINE(  1: 5 ) )
            SCC  = STR2INT ( LINE( 16:25 ) )
            CNOX = 1.0 - STR2REAL( LINE( 26:30 ) )
            CVOC = 1.0 - STR2REAL( LINE( 32:36 ) )

            IF ( FIP .LE. 0 ) THEN      !  all states.
                
                DO  182 S = 1, NPSRC

                    IF( SCC .NE. ISCC( S ) ) GO TO 182  ! n/a control

                    IF ( TFLAG( S ) ) 
     &                  WRITE( RDEV,93010 ) 
     &                      PNAME, IFIP( S ), ISCC( S ), 
     &                      IPLANT( S ), ISTACK( S ),
     &                      CTLEFF( S,2 ), CNOX,
     &                      CTLEFF( S,3 ), CVOC

                    CTLEFF( S,2 ) = CNOX
                    CTLEFF( S,3 ) = CVOC

182             CONTINUE

            ELSE IF ( MOD( FIP, 1000 ) .EQ. 0 ) THEN
                
                SID = FIP / 1000        !  state code only
                I = FIND1( SID, NSID, SIDS )
                IF ( I .EQ. 0 ) GO TO 181

                DO  183 S = JSID( I ), JSID( I+1 ) - 1

                    IF( SCC .NE. ISCC ( S ) ) GO TO 183  ! n/a control

                    IF ( TFLAG( S ) ) 
     &                  WRITE( RDEV,93010 ) 
     &                      PNAME, IFIP( S ), ISCC( S ), 
     &                      IPLANT( S ), ISTACK( S ),
     &                      CTLEFF( S,2 ), CNOX,
     &                      CTLEFF( S,3 ), CVOC

                    CTLEFF( S,2 ) = CNOX
                    CTLEFF( S,3 ) = CVOC

183             CONTINUE

            ELSE        !  fip-ASC specific control
            
                I = FIND1( FIP, NFIP, FFIP )
                IF ( I .EQ. 0 ) GO TO 181

                DO  184 S = JFIP( I ), JFIP( I+1 ) - 1

                    IF( SCC .NE. ISCC ( S ) ) GO TO 184  ! n/a control

                    IF ( TFLAG( S ) ) 
     &                  WRITE( RDEV,93010 ) 
     &                      PNAME, IFIP( S ), ISCC( S ), 
     &                      IPLANT( S ), ISTACK( S ),
     &                      CTLEFF( S,2 ), CNOX,
     &                      CTLEFF( S,3 ), CVOC

                    CTLEFF( S,2 ) = CNOX
                    CTLEFF( S,3 ) = CVOC

184             CONTINUE

            END IF      ! if fip is 0, or nn000, or not.

            GO TO  181  !  end loop on control-lines in this RACT packet

189     CONTINUE        !  exit from / skip around loop 181


C............................................................................
C.......   Process MACT packet of control file:

        PNAME = '/MACT/'
        IF( 0 .EQ. INDEX1( PNAME, NPACKETS, PNAMES ) ) GO TO 199

        IF ( .NOT. GETPKT( PNAME, CDEV ) ) THEN
            CALL M3EXIT( 'CTLPMAT', 0, 0,
     &                   'Could not find ' // PNAME // ' in ' // CNAME, 
     &                   2 )
        END IF

191     CONTINUE        !  head of loop reading lines in the MACT packet

            READ( CDEV,93000,IOSTAT=IOS,END=900 ) LINE
            IREC = IREC + 1
            IF( IOS .NE. 0 ) THEN
                WRITE( MESG,94920 ) 
     &              'Error reading', CNAME( 1:TRIMLEN( CNAME ) ), 
     &              IREC, IOS
                CALL M3EXIT( 'CTLPMAT', 0, 0, MESG, 1 )
            END IF

            CALL UPCASE( LINE )
            IF( 0 .NE. INDEX( LINE, 'END' ) )  GO TO 199

            FIP  = STR2INT ( LINE(  1: 5 ) )
            SCC  = STR2INT ( LINE( 16:25 ) )
            CVOC = 1.0 - STR2REAL( LINE( 26:30 ) )

            IF ( FIP .LE. 0 ) THEN      !  all states.

                DO  193 S = 1, NPSRC

                    IF( SCC .NE. ISCC( S ) ) GO TO 193  ! n/a control

                    IF ( TFLAG( S ) ) 
     &                  WRITE( RDEV,93010 ) 
     &                      PNAME, IFIP( S ), ISCC( S ), 
     &                      IPLANT( S ), ISTACK( S ),
     &                      CTLEFF( S,3 ), CVOC

                    CTLEFF( S,3 ) = CVOC

193             CONTINUE

            ELSE IF ( MOD( FIP, 1000 ) .EQ. 0 ) THEN

                SID = FIP / 1000        !  state code only
                I   = FIND1( SID, NSID, SIDS )

                IF ( I .EQ. 0 ) GO TO 191

                DO  195 S = JSID( I ), JSID( I+1 ) -1

                    IF( SCC .NE. ISCC( S ) ) GO TO 195  ! n/a control

                    IF ( TFLAG( S ) ) 
     &                  WRITE( RDEV,93010 ) 
     &                      PNAME, IFIP( S ), ISCC( S ), 
     &                      IPLANT( S ), ISTACK( S ),
     &                      CTLEFF( S,3 ), CVOC

                    CTLEFF( S,3 ) = CVOC

195             CONTINUE

            ELSE        !  fip-ASC specific control

                I = FIND1( FIP, NFIP, FFIP )

                IF ( I .EQ. 0 ) GO TO 191

                DO  197 S = JFIP( I ), JFIP( I+1 ) - 1

                    IF( SCC .NE. ISCC( S ) ) GO TO 197  ! n/a control

                    IF ( TFLAG( S ) ) 
     &                  WRITE( RDEV,93010 ) 
     &                      PNAME, IFIP( S ), ISCC( S ), 
     &                      IPLANT( S ), ISTACK( S ),
     &                      CTLEFF( S,3 ), CVOC

                    CTLEFF( S,3 ) = CVOC

197             CONTINUE

            END IF      ! if fip is 0, or nn000, or not.

            GO TO  191  !  end loop reading control-line in MACT packet

199     CONTINUE        !  exit from / skip around loop 191


C............................................................................
C.......   Process CTG packet of control file:
        
        PNAME = '/CTG/'
        IF( 0 .EQ. INDEX1( PNAME, NPACKETS, PNAMES ) ) GO TO 209

        IF ( .NOT. GETPKT( PNAME, CDEV ) ) THEN
            CALL M3EXIT( 'CTLPMAT', 0, 0,
     &                   'Could not find ' // PNAME // ' in ' // CNAME, 
     &                   2 )
        END IF

201     CONTINUE        !  head of loop reading this category of packet

            READ( CDEV,93000,IOSTAT=IOS,END=900 ) LINE
            IREC = IREC + 1
            IF( IOS .NE. 0 ) THEN
                WRITE( MESG,94920 ) 
     &              'Error reading', CNAME( 1:TRIMLEN( CNAME ) ), 
     &              IREC, IOS
                CALL M3EXIT( 'CTLPMAT', 0, 0, MESG, 1 )
            END IF
                   
            CALL UPCASE( LINE )
            IF( 0 .NE. INDEX( LINE, 'END' ) )  GO TO 209
            
            FIP  = STR2INT ( LINE(  1: 5 ) )
            SCC  = STR2INT ( LINE( 16:25 ) )
            CNOX = 1.0 - STR2REAL( LINE( 27:31 ) )
            CVOC = 1.0 - STR2REAL( LINE( 33:37 ) )
            TNOX = STR2REAL( LINE( 39:48 ) )
            TVOC = STR2REAL( LINE( 50:59 ) )

            IF ( FIP .LE. 0 ) THEN      !  all states.
                
                DO  202 S = 1, NPSRC

                    IF( SCC .NE. ISCC( S ) ) GO TO 202  ! n/a control

                    IF ( NOX( S ) .GT. TNOX ) THEN

                        IF( VOC( S ) .GT. TVOC ) THEN   ! NOX,VOC both apply

                            IF ( TFLAG( S ) ) 
     &                          WRITE( RDEV,93010 ) 
     &                              PNAME, IFIP( S ), ISCC( S ), 
     &                              IPLANT( S ), ISTACK( S ),
     &                              CTLEFF( S,2 ), CNOX,
     &                              CTLEFF( S,3 ), CVOC

                            CTLEFF( S,2 ) = CNOX
                            CTLEFF( S,3 ) = CVOC

                        ELSE                            ! NOX control applies

                            IF ( TFLAG( S ) ) 
     &                          WRITE( RDEV,93010 ) 
     &                              PNAME//'NOX', IFIP( S ), ISCC( S ), 
     &                              IPLANT( S ), ISTACK( S ),
     &                              CTLEFF( S,2 ), CNOX

                            CTLEFF( S,2 ) = CNOX

                        END IF

                    ELSE IF( VOC( S ) .GT. TVOC ) THEN  ! VOC control applies                            T = RULEFF( S ) * RULPEN( S )

                        IF ( TFLAG( S ) ) 
     &                      WRITE( RDEV,93010 ) 
     &                          PNAME//'VOC', IFIP( S ), ISCC( S ), 
     &                          IPLANT( S ), ISTACK( S ),
     &                          CTLEFF( S,3 ), CVOC

                        CTLEFF( S,3 ) = CVOC

                    END IF      !  if strongest

202             CONTINUE        !  end loop on sources S, state-ind ctls

            ELSE IF ( MOD( FIP, 1000 ) .EQ. 0 ) THEN
                
                SID = FIP / 1000        !  state code only
                I   = FIND1( SID, NSID, SIDS )

                IF ( I .EQ. 0 )  GO TO  201

                DO  203 S = JSID( I ), JSID( I+1 ) -1

                    IF( SCC .NE. ISCC ( S ) ) GO TO 203  ! n/a control

                    IF ( NOX( S ) .GT. TNOX ) THEN
                        IF( VOC( S ) .GT. TVOC ) THEN   ! NOX,VOC both apply
                             
                            IF ( TFLAG( S ) ) 
     &                          WRITE( RDEV,93010 ) 
     &                              PNAME, IFIP( S ), ISCC( S ), 
     &                              IPLANT( S ), ISTACK( S ),
     &                              CTLEFF( S,2 ), CNOX,
     &                              CTLEFF( S,3 ), CVOC

                            CTLEFF( S,2 ) = CNOX
                            CTLEFF( S,3 ) = CVOC

                        ELSE                            ! NOX control applies

                            IF ( TFLAG( S ) ) 
     &                          WRITE( RDEV,93010 ) 
     &                              PNAME//'NOX', IFIP( S ), ISCC( S ), 
     &                              IPLANT( S ), ISTACK( S ),
     &                              CTLEFF( S,2 ), CNOX

                            CTLEFF( S,2 ) = CNOX

                        END IF

                    ELSE IF( VOC( S ) .GT. TVOC ) THEN  ! VOC control applies                            T = RULEFF( S ) * RULPEN( S )

                        IF ( TFLAG( S ) ) 
     &                      WRITE( RDEV,93010 ) 
     &                          PNAME//'VOC', IFIP( S ), ISCC( S ), 
     &                          IPLANT( S ), ISTACK( S ),
     &                          CTLEFF( S,3 ), CVOC

                        CTLEFF( S,3 ) = CVOC

                    END IF      !  if strongest

203             CONTINUE        !  end loop on sources S, state code CTG ctls

            ELSE        !  fip-ASC specific control
            
                S = FIND1( FIP, NPSRC, IFIP )

                I   = FIND1( FIP, NFIP, FFIP )

                IF ( I .EQ. 0 )  GO TO  201

                DO  204 S = JFIP( I ), JFIP( I+1 ) - 1

                    IF( SCC .NE. ISCC ( S ) ) GO TO 204  ! n/a control

                    IF ( NOX( S ) .GT. TNOX ) THEN

                        IF( VOC( S ) .GT. TVOC ) THEN   ! NOX,VOC both apply

                            IF ( TFLAG( S ) ) 
     &                          WRITE( RDEV,93010 ) 
     &                              PNAME, IFIP( S ), ISCC( S ), 
     &                              IPLANT( S ), ISTACK( S ),
     &                              CTLEFF( S,2 ), CNOX,
     &                              CTLEFF( S,3 ), CVOC

                            CTLEFF( S,2 ) = CNOX
                            CTLEFF( S,3 ) = CVOC

                        ELSE                            ! NOX control applies

                            IF ( TFLAG( S ) ) 
     &                          WRITE( RDEV,93010 ) 
     &                              PNAME//'NOX', IFIP( S ), ISCC( S ), 
     &                              IPLANT( S ), ISTACK( S ),
     &                              CTLEFF( S,2 ), CNOX

                            CTLEFF( S,2 ) = CNOX

                        END IF

                    ELSE IF( VOC( S ) .GT. TVOC ) THEN  ! VOC control applies                            T = RULEFF( S ) * RULPEN( S )

                        IF ( TFLAG( S ) ) 
     &                      WRITE( RDEV,93010 ) 
     &                          PNAME//'VOC', IFIP( S ), ISCC( S ), 
     &                          IPLANT( S ), ISTACK( S ),
     &                          CTLEFF( S,3 ), CVOC

                        CTLEFF( S,3 ) = CVOC

                    END IF      !  if strongest
         
204             CONTINUE       !  end loop on S

            END IF      ! if fip is 0, or nn000, or not.

            GO TO  201  !  end loop on control-lines in this CTG packet

209     CONTINUE        !  exit from / skip around loop 201


C............................................................................
C.......   Compute the control matrix:

        DO  222  I = 1, NIPOL
        DO  211  S = 1, NPSRC

            CMAT( S,I ) = CMAT( S,I ) * PRJFAC( S ) *
     &          ( 1.0 - CTLEFF( S,I ) * RULEFF( S,I ) * RULPEN( S,I ) )

211     CONTINUE
222     CONTINUE

C............................................................................
C.......   Write out the control matrix:

900   CONTINUE          !  Done processing; write output file:

        WRITE( *,92000 ) ' ', 'Writing out CONTROL MATRIX...', ' '

        IF ( .NOT. WRITE3( MNAME, 'ALL', 0, 0,CMAT ) ) THEN
            CALL M3EXIT( 'CTLPMAT', 0, 0, 
     &                   'Error writing CONTROL MATRIX file.', 2 )
        END IF



999   CONTINUE          !  exit program

      CALL M3EXIT( 'CTLPMAT', 0, 0,
     &             '====  END PROGRAM CTLPMAT  ====', 0 )

C******************  FORMAT  STATEMENTS   ******************************

C...........   Informational (LOG) message formats... 92xxx

92000   FORMAT( 5X, A )


C...........   Formatted file I/O formats............ 93xxx

93000   FORMAT( A )

93010   FORMAT( A, ':', I7, ',', I8, ',', I4.4, ',', I3.3, ' -- ', 
     &          3( F8.5, '~>', F8.5, : '; ' )  )

93020   FORMAT( A, ':', I4, ' to ', I4, 1X, I5.5, 1X, I4,
     &          3( F8.5, '~>', F8.5, : '; ' )  )


C...........   Internal buffering formats............ 94xxx

94910   FORMAT( A, :, ' at line', I7, :,  '; IOSTAT=', I7 )

94920   FORMAT( A, ' "', A, :,  '" at line', I7, :, '; IOSTAT=', I7 )

94930   FORMAT( A )


        END

