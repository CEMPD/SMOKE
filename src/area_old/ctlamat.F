
        PROGRAM CTLAMAT

C***********************************************************************
C  program body starts at line 160
C
C  DESCRIPTION:
C       Construct area source control matrix from data 
C       contained in prototype AREA and EPS-style control files.
C
C  PRECONDITIONS REQUIRED:
C       Sorted, cut-down control file.
C       M3IO area-sources file from program RAWAREA
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C       Models-3 I/O
C       FIND1, FIND2, FIND3, GETYN, TRIMLEN
C
C  REVISION  HISTORY:
C       Prototype  6/95 by CJC.
C
C***********************************************************************
C
C Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE)
C                Modeling System
C File Version @(#)$Id$
C Pathname:    $Source$
C Last updated: $Date$ 
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
C***************************************************************************

      IMPLICIT NONE

C...........   INCLUDES:
        
        INCLUDE 'ARDIMS3.EXT'   !  area-source dimensioning parameters
        INCLUDE 'CHDIMS3.EXT'   !  emis chem parms (inventory + model)
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
        INCLUDE 'FDESC3.EXT'    !  I/O API file description data structures.


C...........   EXTERNAL FUNCTIONS and their descriptions:
        
        INTEGER         FIND1, FIND2, FIND3
        LOGICAL         GETPKT
        LOGICAL         GETYN
        INTEGER         INDEX1
        INTEGER         LBLANK
        INTEGER         PROMPTFFILE
        CHARACTER*16    PROMPTMFILE
        INTEGER         STR2INT
        REAL            STR2REAL
        INTEGER         TRIMLEN
        
        EXTERNAL        FIND1, FIND2, FIND3, GETPKT, GETYN, INDEX1, 
     &                  LBLANK, PROMPTFFILE, PROMPTMFILE, 
     &                  STR2INT, STR2REAL, TRIMLEN


C...........   PARAMETERS and their descriptions:
 
        INTEGER         MAXPKT          !  max # of packets in control file
        PARAMETER     ( MAXPKT = 20 )
 
 
C...........   LOCAL VARIABLES and their descriptions:
C...........   NOTE that ASC (area-source-category) ID's are 10-digit
C...........   unsigned integers which may be treated as a leading 7-digit
C...........   field, and a trailing 3-digit field.  *7 and *3 arrays below
C...........   follow this scheme with parallel arrays

        INTEGER     FIP         !  FIP (state+county) code
        INTEGER     SID         !  state ID (leading-2 digits of FIP code)
        INTEGER     ID4         !  leading-4 digits of ASC
        INTEGER     ID7         !  leading-7 digits of ASC
        INTEGER     ID3         !  trailing-3 digits of ASC
        REAL        CNOX        !  NOX control factor
        REAL        CVOC        !  VOC control factor
        REAL        CCO         !  CO control factor
        REAL        PFAC        !  projection factor
        REAL        TNOX        !  NOX control threshold
        REAL        TVOC        !  VOC control threshold
        INTEGER     LFIP        !  last FIP (state+county) code
        INTEGER     LSID        !  last state ID
        INTEGER     NFIP        !  subscript and count for CFIP
        INTEGER     NSID        !  subscript and count for ISID
        INTEGER     DYR         !  desination year for projection factor
        INTEGER     YEAR        !  inventory year for projection factor
                                 
C.......   Control matrix:
        
        REAL        CMAT( NASRC, NIPOL )    !  control coefficients


C...........   Area Sources Table
        
        INTEGER     IFIP  ( NASRC )  !  source FIPS (state+county) ID
        INTEGER     ASC7  ( NASRC )  !  leading-7  digits of source ASC
        INTEGER     ASC3  ( NASRC )  !  trailing-3 digits of source ASC
        INTEGER     INVYR ( NASRC )  !  Inventory year
        REAL        CTLEFF( NASRC,NIPOL )  !  control-factor (efficiency)
        REAL        RULEFF( NASRC,NIPOL )  !  control-rule effectiveness
        REAL        RULPEN( NASRC,NIPOL )  !  control-rule pentration
        REAL        PRJFAC( NASRC )  !  projection factors
        REAL        NOX   ( NASRC )  !  NOX emissions
        REAL        VOC   ( NASRC )  !  VOC emissions
        REAL        CO    ( NASRC )  !  VOC emissions

C.......   Auxilliary tables:

        INTEGER     ASC4  ( NASRC )  !  leading-4  digits of source ASC
        INTEGER     SIDS  ( NASID )  !  source FIPS (state-only) ID
        INTEGER     ISID  ( NASID+1) !  Subscript at which state ID occurs
        INTEGER     FFIP  ( NAFIP )  !  source FIPS list
        INTEGER     CFIP  ( NAFIP+1) !  Subscript at which FIP occurs


C...........   Other local variables
        
        INTEGER         S, I, J, L, V!  loop counters.
        INTEGER         CDEV         !  for control file
        INTEGER         LDEV         !  log-device
        INTEGER         RDEV         !  report-file device
        INTEGER         TDEV         !  device for report-these-sources file
        LOGICAL         PFLAG        !  true if projection entries have years
        LOGICAL         RFLAG        !  flag:  make report-file
        LOGICAL         TFLAG( NASRC )  !  flags:  track these sources
        
        INTEGER         IOS          !  I/O status
        INTEGER         IREC         !  input line (record) number

        CHARACTER*16    ENAME   !  logical name for area-source input file
        CHARACTER*16    CNAME   !  logical name for control     input file
        CHARACTER*16    MNAME   !  logical name for ctl matrix output file
        CHARACTER*16    VNAME   !  scratch buffer for variable-names
        CHARACTER*240   LINE    !  input line buffer
        CHARACTER*256   MESG    !  error message buffer

        INTEGER         NPACKETS
        CHARACTER*20    PNAMES( MAXPKT )
        CHARACTER*20    PNAME


C***********************************************************************
C   begin body of program CTLAMAT
        
        LDEV = INIT3()
        
        WRITE( *,92000 ) 
     &  ' ',
     &  'Program CTLAMAT to take the netCDF AREA SOURCES file, an',
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
            CALL M3EXIT( 'CTLAMAT', 0, 0, 'Ending program.', 2 )
        END IF


C.......   Get file name; open input area sources, control definitions files

        ENAME = PROMPTMFILE( 
     &          'Enter logical name for AREA SOURCE input file',
     &          FSREAD3, 'AREA', 'CTLAMAT' )

        CDEV = PROMPTFFILE( 
     &           'Enter logical name for CONTROL file',
     &           .TRUE., .TRUE., 'ACTRL', 'CTLAMAT' )

        RDEV = PROMPTFFILE( 
     &           '...for REPORT file, or "NONE"',
     &           .FALSE., .TRUE., 'ACREPT', 'CTLAMAT' )

        RFLAG = ( RDEV .GT. 0 )

        IF ( RFLAG ) THEN

            TDEV = PROMPTFFILE( 
     &           '...for TRACKING LIST, or "ALL"',
     &           .TRUE., .TRUE., 'ACTRACK', 'CTLAMAT' )

        END IF


C.......   Read in description of area-source file
        
        IF ( .NOT. DESC3( ENAME ) ) THEN
            CALL M3ERR( 'CTLAMAT', 0, 0, 
     &                  'Could not get description of file ' // ENAME,
     &                  .TRUE. )
        END IF
        
        
C.......   Open control matrix file:  reuse much of AREA file description
        
        NVARS3D = 3
        FTYPE3D = GRDDED3

        DO 11 V = 1, NIPOL
            VNAME3D( V ) = 'CONTROL_' // EINAM( V )
            UNITS3D( V ) = 'n/a'
            VDESC3D( V ) = 'Control factor for this source, variable ' 
     &                     // EINAM( V )
            VTYPE3D( V ) = M3REAL
11      CONTINUE

        FDESC3D( 1 ) = 'Matrix of control factors, by source'
        DO  22  I = 2, MXDESC3
            FDESC3D( I ) = ' '
22      CONTINUE

        MNAME = PROMPTMFILE( 
     &          'Enter logical name for CONTROL MATRIX output file',
     &          FSUNKN3, 'ACMAT', 'CTLAMAT' )


C...........   Read in emissions indexing data IFIP, ASC7, ASC3, INVYR;
C...........   existing-control data CTLEFF, RULEFF, RULPEN:
        
        IF ( .NOT. READ3( ENAME, 'FIP', ALLAYS3, 0, 0, IFIP ) ) THEN
            CALL M3EXIT( 'CTLAMAT', 0, 0, 
     &                   'Error reading FIP from file' // ENAME, 2 )
        END IF

        IF ( .NOT. READ3( ENAME, 'ASC7', ALLAYS3, 0, 0, ASC7 ) ) THEN
            CALL M3EXIT( 'CTLAMAT', 0, 0, 
     &                   'Error reading ASC7 from file' // ENAME, 2 )
        END IF

        IF ( .NOT. READ3( ENAME, 'ASC3', ALLAYS3, 0, 0, ASC3 ) ) THEN
            CALL M3EXIT( 'CTLAMAT', 0, 0, 
     &                    'Error reading ASC3 from file' // ENAME, 2 )
        END IF

        IF ( .NOT. READ3( ENAME, 'INVYR', ALLAYS3, 0, 0, INVYR ) ) THEN
            CALL M3EXIT( 'CTLAMAT', 0, 0, 
     &                   'Error reading INVYR from file' // ENAME, 2 )
        END IF

        DO  33  I = 1, NIPOL
 
            VNAME = 'CTLEFF_' // EINAM( I )
            IF ( .NOT. READ3( ENAME, VNAME, ALLAYS3, 
     &                        0, 0, CTLEFF( 1,I ) ) ) THEN
                CALL M3EXIT( 'CTLAMAT', 0, 0, 
     &            'Error reading '// VNAME //' from file' // ENAME, 2 )
            END IF
 
            VNAME = 'RULEFF_' // EINAM( I )
            IF ( .NOT. READ3( ENAME, VNAME, ALLAYS3, 
     &                        0, 0, RULEFF( 1,I ) ) ) THEN
                CALL M3EXIT( 'CTLAMAT', 0, 0, 
     &            'Error reading '// VNAME //' from file' // ENAME, 2 )
            END IF
 
            VNAME = 'RULPEN_' // EINAM( I )
            IF ( .NOT. READ3( ENAME, VNAME, ALLAYS3, 
     &                        0, 0, RULPEN( 1,I ) ) ) THEN
                CALL M3EXIT( 'CTLAMAT', 0, 0, 
     &            'Error reading '// VNAME //' from file' // ENAME, 2 )
            END IF
 
33      CONTINUE                !  end loop reading ctleff, ruleff, rulpen of V
 
        IF ( .NOT. READ3( ENAME, EINAM( 1 ), ALLAYS3, 
     &                    0, 0, CO ) ) THEN
            CALL M3ERR( 'CTLAMAT', 0, 0, 
     &                  'Error reading CO from file' // ENAME,
     &                  .TRUE. )
        END IF

        IF ( .NOT. READ3( ENAME, EINAM( 2 ), ALLAYS3, 
     &                    0, 0, NOX ) ) THEN
            CALL M3ERR( 'CTLAMAT', 0, 0, 
     &                  'Error reading NOX from file' // ENAME,
     &                  .TRUE. )
        END IF

        IF ( .NOT. READ3( ENAME, EINAM( 3 ), ALLAYS3, 
     &                    0, 0, VOC ) ) THEN
            CALL M3ERR( 'CTLAMAT', 0, 0, 
     &                  'Error reading VOC from file' // ENAME,
     &                  .TRUE. )
        END IF

        
C...........   Compute the control-coefficient matrix.
C...........   Note conversion from percent to fraction for
C...........   non-missing values.
C...........   Initialization; precompute states, ASC4's for each source
            
        DO  44  I = 1, NIPOL
        DO  43  S = 1, NASRC

            IF( RULEFF( S,I ) .LT. AMISS3 ) THEN
                RULEFF( S,I ) = 1.0
            ELSE
                RULEFF( S,I ) = 0.01 * RULEFF( S,I )
            END IF
            IF( RULPEN( S,I ) .LT. AMISS3 ) THEN
                RULPEN( S,I ) = 1.0
            ELSE
                RULPEN( S,I ) = 0.01 * RULPEN( S,I )
            END IF
            IF( CTLEFF( S,I ) .LT. AMISS3 ) THEN
                CTLEFF( S,I ) = 0.0
            ELSE
                CTLEFF( S,I ) = 0.01 * CTLEFF( S,I )
            END IF
            
            CMAT( S,I ) = 1.0 / 
     &          ( 1.0 - CTLEFF( S,I ) * RULEFF( S,I ) * RULPEN( S,I ) )
 
43      CONTINUE	!  end loop on sources S
44      CONTINUE	!  end loop on sources S

        LFIP = -1
        LSID = -1
        NFIP =  0
        NSID =  0

        DO  55  S = 1, NASRC

            PRJFAC( S ) = 1.0

            FIP        = IFIP( S )
            SID        = FIP / 1000
            
            IF ( FIP .NE. LFIP ) THEN
                NFIP = NFIP + 1
                IF ( NFIP .GE. NAFIP ) THEN
                    CALL M3EXIT( 'CTLAMAT', 0, 0, 
     &                           'FIP-table maximum size exceeded', 2 )
                ELSE
                    CFIP( NFIP ) = S
                    FFIP( NFIP ) = FIP
                    LFIP         = FIP
                END IF
            END IF

            IF ( SID .NE. LSID ) THEN
                NSID = NSID + 1
                IF ( NSID .GE. NASID ) THEN
                    CALL M3EXIT( 'CTLAMAT', 0, 0, 
     &                           'SID-table maximum size exceeded', 2 )
                ELSE
                    ISID( NSID ) = S
                    SIDS( NSID ) = SID
                    LSID         = SID
                END IF
            END IF

            ASC4 ( S ) = ASC7( S ) / 1000

55      CONTINUE	!  end loop on sources S


        CFIP( NFIP+1 ) = NASRC + 1  	!  sentinel
        ISID( NSID+1 ) = NASRC + 1  	!  sentinel

C............................................................................
C.......   Process list of FIP-SCC pairs to be tracked in the 
C.......   REPORT file (if any).
            
        DO  66  S = 1, NASRC
            TFLAG( S ) = .FALSE.
66      CONTINUE        !  end loop on sources S
 
        IF ( RFLAG ) THEN
 
            IF ( TDEV .GT. 0 ) THEN
 
                IREC = 0
 
75              CONTINUE
 
                    IREC = IREC + 1
                    READ( TDEV, 93020, END = 77, IOSTAT=IOS ) 
     &                  FIP, ID7, ID3
 
                    IF ( IOS .NE. 0 ) THEN
                        WRITE( MESG,94910 )
     &                      'Error reading TRACKING LIST', IREC, IOS
                        CALL M3EXIT( 'CTLAMAT', 0, 0, MESG, 1 )
                    END IF
                    
                    I = FIND3( FIP, ID7, ID3, NASRC, IFIP, ASC7, ASC3 )
 
                    IF ( I .GT. 0 )   TFLAG( I ) = .TRUE.
 
                    GO TO 75
 
77              CONTINUE        !  end loop reading file:  target of "END=..."
 
            ELSE IF ( TDEV .EQ. -3 ) THEN       ! "ALL"
 
                DO  88  S = 1, NASRC
                    TFLAG( S ) = .TRUE.
88              CONTINUE        !  end loop on sources S
 
            ELSE
 
                CALL M3EXIT( 'CTLAMAT', 0, 0, 
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
C.......        ALLOWABLE packet     -- IGNORED
C.......        PROJECT AMS packet
C.......        PROJECT PTS packet   -- IGNORED
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
                WRITE( MESG,94910 ) 
     &                  'Error reading', CNAME( 1:TRIMLEN( CNAME ) ),
     &                  IREC, IOS
                CALL M3EXIT( 'CTLAMAT', 0, 0, MESG, 1 )
            END IF
            CALL UPCASE( LINE )
            IF ( 0 .NE. INDEX( LINE, '/END/' ) ) GO TO 90
 
            I = INDEX( LINE, '/' )
            IF ( I .EQ. 0 ) GO TO 90
 
            J = INDEX( LINE( I+1 : 240 ), '/' )
 
            IF ( J .EQ. 0 ) THEN
 
                WRITE( MESG,94910 )
     &              'Malformed packet in ',
     &              CNAME( 1:TRIMLEN( CNAME ) ),
     &              IREC
                CALL M3EXIT( 'CTLAMAT', 0, 0, MESG, 1 )
 
            ELSE IF ( J .GT. 20 ) THEN
 
                WRITE( MESG,94910 )
     &              'Packet name too long in ',
     &              CNAME( 1:TRIMLEN( CNAME ) ),
     &              IREC
                CALL M3EXIT( 'CTLAMAT', 0, 0, MESG, 1 )
 
            END IF
 
            IF ( NPACKETS .LT. MAXPKT ) THEN
                NPACKETS = NPACKETS + 1
                PNAMES( NPACKETS ) = LINE( I:I+J )
                WRITE( LDEV,92000 )
     &              'Control packet ', PNAMES( NPACKETS ), ' present'
                GO TO 90
            ELSE
                WRITE( MESG,94910 )
     &              'Max packet count exceeded in ',
     &              CNAME( 1:TRIMLEN( CNAME ) ), IREC
                CALL M3EXIT( 'CTLAMAT', 0, 0, MESG, 1 )
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

            READ( CDEV,93000,IOSTAT=IOS,END=109 ) LINE
            IREC = IREC + 1
            IF( IOS .NE. 0 ) THEN
                WRITE( MESG,94910 ) 
     &              'Error reading',
     &              CNAME( 1:TRIMLEN( CNAME ) ), 
     &              IREC, IOS
                CALL M3EXIT( 'CTLAMAT', 0, 0, MESG, 1 )
            END IF

            CALL UPCASE( LINE )
            IF( 0 .NE. INDEX( LINE, 'END' ) )  GO TO 109

            CALL LUSTR( LINE( 21:30 ) )
            FIP  = STR2INT ( LINE(  1: 5 ) )
            CNOX = 1.0 - STR2REAL( LINE( 32:41 ) )
            CVOC = 1.0 - STR2REAL( LINE( 42:51 ) )
            CCO  = 1.0 - STR2REAL( LINE( 52:61 ) )
            L    = TRIMLEN ( LINE( 21:30 ) ) - LBLANK( LINE( 21:30 ) )
            IF ( L .EQ. 4 ) THEN
                ID4  = STR2INT ( LINE( 21:30 ) )
                ID7  = ID4 * 1000
                ID3  = 0
            ELSE
                ID7  = STR2INT ( LINE( 21:27 ) )
                ID3  = STR2INT ( LINE( 28:30 ) )
                ID4  = ID7 / 1000
            END IF


            IF ( FIP .EQ. 0 ) THEN      !  all states.

                DO  102 S = 1, NASRC

                    IF( ID7 .NE. ASC7( S ) ) GO TO 102  ! n/a control
                    IF( ID3 .NE. ASC3( S ) ) GO TO 102  ! n/a control

                    IF ( TFLAG( S ) ) 
     &                  WRITE( RDEV,93030 ) 
     &                      PNAME, IFIP( S ), ASC7( S ), ASC3( S ),
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

                DO  103 S = ISID( I ), ISID( I+1 ) - 1

                    IF( ID7 .NE. ASC7( S ) ) GO TO 103  ! n/a control
                    IF( ID3 .NE. ASC3( S ) ) GO TO 103  ! n/a control

                    IF ( TFLAG( S ) ) 
     &                  WRITE( RDEV,93030 ) 
     &                      PNAME, IFIP( S ), ASC7( S ), ASC3( S ),
     &                      CTLEFF( S,1 ), CCO,
     &                      CTLEFF( S,2 ), CNOX,
     &                      CTLEFF( S,3 ), CVOC

                    CTLEFF( S,1 ) = CCO
                    CTLEFF( S,2 ) = CNOX
                    CTLEFF( S,3 ) = CVOC

103             CONTINUE

            ELSE	!  fip-ASC specific control

                S = FIND3( FIP, ID7, ID3, NASRC, IFIP, ASC7, ASC3 )

                IF ( S .LE. 0 ) GO TO 101

                IF( ID7 .NE. ASC7( S ) ) GO TO 101  ! n/a control
                IF( ID3 .NE. ASC3( S ) ) GO TO 101  ! n/a control

                IF ( TFLAG( S ) ) 
     &              WRITE( RDEV,93030 ) 
     &                      PNAME, IFIP( S ), ASC7( S ), ASC3( S ),
     &                      CTLEFF( S,1 ), CCO,
     &                      CTLEFF( S,2 ), CNOX,
     &                      CTLEFF( S,3 ), CVOC

                CTLEFF( S,1 ) = CCO
                CTLEFF( S,2 ) = CNOX
                CTLEFF( S,3 ) = CVOC

            END IF	! if fip is 0, or nn000, or not.

            GO TO  101  !  end loop reading lines in SOURCE CATEGORY packet

109     CONTINUE        !  exit from / skip around loop 101


C............................................................................
C.......   Process FIPS CODE packet of control file:

        PNAME = '/FIPS CODE/'
        IF( 0 .EQ. INDEX1( PNAME, NPACKETS, PNAMES ) ) GO TO 119

        IF ( .NOT.  GETPKT( PNAME, CDEV ) ) THEN
            CALL M3EXIT( 'CTLAMAT', 0, 0,
     &                   'Could not find ' // PNAME // ' in ' // CNAME, 
     &                   2 )
        END IF

111     CONTINUE        !  head of loop reading lines in the FIPS CODE packet

            READ( CDEV,93000,IOSTAT=IOS,END=119 ) LINE
            IREC = IREC + 1
            IF( IOS .NE. 0 ) THEN
                WRITE( MESG,94910 ) 
     &              'Error reading',
     &              CNAME( 1:TRIMLEN( CNAME ) ), 
     &              IREC, IOS
                CALL M3EXIT( 'CTLAMAT', 0, 0, MESG, 1 )
            END IF

            CALL UPCASE( LINE )
            IF( 0 .NE. INDEX( LINE, 'END' ) )  GO TO 119

            FIP  = STR2INT ( LINE(  1: 5 ) )
            CNOX = 1.0 - STR2REAL( LINE( 16:25 ) )
            CVOC = 1.0 - STR2REAL( LINE( 26:35 ) )
            CCO  = 1.0 - STR2REAL( LINE( 36:45 ) )

            IF ( FIP .EQ. 0 ) THEN      !  all states.

                DO  112 S = 1, NASRC

                    IF ( TFLAG( S ) ) 
     &                  WRITE( RDEV,93030 ) 
     &                      PNAME, IFIP( S ), ASC7( S ), ASC3( S ),
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

                DO  113 S = ISID( I ), ISID( I+1 ) - 1

                    IF ( TFLAG( S ) ) 
     &                  WRITE( RDEV,93030 ) 
     &                      PNAME, IFIP( S ), ASC7( S ), ASC3( S ),
     &                      CTLEFF( S,1 ), CCO,
     &                      CTLEFF( S,2 ), CNOX,
     &                      CTLEFF( S,3 ), CVOC

                    CTLEFF( S,1 ) = CCO
                    CTLEFF( S,2 ) = CNOX
                    CTLEFF( S,3 ) = CVOC

113             CONTINUE

            ELSE	!  fip-ASC specific control

                I = FIND1( FIP, NFIP, FFIP )

                IF ( I .EQ. 0 ) GO TO 111

                DO  114 S = CFIP( I ), CFIP( I+1 ) - 1

                    IF ( TFLAG( S ) ) 
     &                  WRITE( RDEV,93030 ) 
     &                      PNAME, IFIP( S ), ASC7( S ), ASC3( S ),
     &                      CTLEFF( S,1 ), CCO,
     &                      CTLEFF( S,2 ), CNOX,
     &                      CTLEFF( S,3 ), CVOC

                    CTLEFF( S,1 ) = CCO
                    CTLEFF( S,2 ) = CNOX
                    CTLEFF( S,3 ) = CVOC

114             CONTINUE

            END IF	! if fip is 0, or nn000, or not.

            GO TO  111  !  end loop reading lines in FIPS CODE packet

119     CONTINUE        !  exit from / skip around loop 111


C............................................................................
C.......   Process PROJECT AMS packet of control file:

        PNAME = '/PROJECT AMS/'
        IF( 0 .EQ. INDEX1( PNAME, NPACKETS, PNAMES ) ) GO TO 129

        IF ( .NOT.  GETPKT( PNAME, CDEV ) ) THEN
            CALL M3EXIT( 'CTLAMAT', 0, 0,
     &                   'Could not find ' // PNAME // ' in ' // CNAME, 
     &                   2 )
        END IF

C.........  Prompt for type of projection entries: default (EPS) or
C           modified EPS (with year)
        PFLAG = GETYN(
     &        'Are projection entries in year-specific format?', .TRUE.)

121     CONTINUE        !  head of loop reading lines in the PROJECT AMS packet

            READ( CDEV,93000,IOSTAT=IOS,END=129 ) LINE
            IREC = IREC + 1
            IF( IOS .NE. 0 ) THEN
                WRITE( MESG,94910 ) 
     &              'Error reading',
     &              CNAME( 1:TRIMLEN( CNAME ) ), 
     &              IREC, IOS
                CALL M3EXIT( 'CTLAMAT', 0, 0, MESG, 1 )
            END IF

            CALL UPCASE( LINE )
            IF( 0 .NE. INDEX( LINE, 'END' ) )  GO TO 129

            IF( PFLAG ) THEN
                YEAR = STR2INT ( LINE(  1: 4 ) )
                DYR  = STR2INT ( LINE(  6: 9 ) )
                FIP  = STR2INT ( LINE( 11:15 ) )
                PFAC = STR2REAL( LINE( 37:46 ) )

                L    = TRIMLEN( LINE( 23:32 ) ) - 
     &                 LBLANK ( LINE( 23:32 ) )

                IF ( L .EQ. 4 ) THEN
                    ID4  = STR2INT ( LINE( 23:32 ) )
                    ID7  = ID4 * 1000
                    ID3  = 0
                ELSE
                    ID7  = STR2INT ( LINE( 23:29 ) )
                    ID3  = STR2INT ( LINE( 30:32 ) )
                    ID4  = ID7 / 1000
                END IF

            ELSE
                FIP  = STR2INT ( LINE(  1: 5 ) )
                PFAC = STR2REAL( LINE( 27:36 ) )

                L    = TRIMLEN( LINE( 16:25 ) ) - 
     &                 LBLANK ( LINE( 16:25 ) )

                IF ( L .EQ. 4 ) THEN
                    ID4  = STR2INT ( LINE( 16:25 ) )
                    ID7  = ID4 * 1000
                    ID3  = 0
                ELSE
                    ID7  = STR2INT ( LINE( 16:22 ) )
                    ID3  = STR2INT ( LINE( 23:25 ) )
                    ID4  = ID7 / 1000
                END IF
            ENDIF

            IF ( L .EQ. 4 ) THEN	!  using ID4
            
                IF ( FIP .EQ. 0 ) THEN      !  all states.

                    DO  122 S = 1, NASRC

                        IF( ID4 .NE. ASC4( S ) ) THEN
                            GO TO 122  ! n/a control

                        ELSEIF( PFLAG .AND. INVYR( S ) .NE. YEAR ) THEN
                            GO TO 122
 
                        ELSEIF( PFLAG .AND. TFLAG( S ) ) THEN
                            WRITE( RDEV,93040 )
     &                          PNAME, YEAR, DYR, IFIP( S ), ASC7( S ),
     &                          ASC3( S ), PRJFAC( S ), PFAC
                         
                        ELSEIF ( TFLAG( S ) ) THEN
                            WRITE( RDEV,93030 ) 
     &                          PNAME, IFIP( S ), ASC7( S ), ASC3( S ),
     &                          PRJFAC( S ), PFAC
                        ENDIF

                        PRJFAC( S ) = PFAC

122                 CONTINUE

                ELSE IF ( MOD( FIP, 1000 ) .EQ. 0 ) THEN

                    SID = FIP / 1000        !  state code only
                    I = FIND1( SID, NSID, SIDS )
                    IF ( I .EQ. 0 ) GO TO 121

                    DO  123 S = ISID( I ), ISID( I+1 ) -1

                        IF( ID4 .NE. ASC4( S ) ) THEN
                            GO TO 123  ! n/a control
 
                        ELSEIF ( PFLAG .AND. INVYR( S ) .NE. YEAR ) THEN
                            GO TO 123

                        ELSEIF ( PFLAG .AND. TFLAG( S ) ) THEN
                            WRITE( RDEV,93040 ) 
     &                          PNAME, YEAR, DYR, IFIP( S ), ASC7( S ),
     &                          ASC3( S ), PRJFAC( S ), PFAC

                        ELSEIF ( TFLAG( S ) ) THEN
                            WRITE( RDEV,93030 ) 
     &                          PNAME, IFIP( S ), ASC7( S ), ASC3( S ),
     &                          PRJFAC( S ), PFAC
                        ENDIF

                        PRJFAC( S ) = PFAC

123                 CONTINUE

                ELSE	!  fip-ASC specific control

                    I = FIND1( FIP, NFIP, FFIP )

                    IF ( I .EQ. 0 ) GO TO 121

                    DO 124 S = CFIP( I ), CFIP( I+1 ) - 1

                        IF( ID4 .NE. ASC4( S ) ) THEN
                            GO TO 124 ! n/a control
                         
                        ELSEIF ( PFLAG .AND. INVYR( S ) .NE. YEAR ) THEN
                            GO TO 124

                        ELSEIF ( PFLAG .AND. TFLAG( S ) ) THEN
                            WRITE( RDEV,93040 ) 
     &                          PNAME, YEAR, DYR, IFIP( S ), ASC7( S ),
     &                          ASC3( S ), PRJFAC( S ), PFAC

                        ELSEIF ( TFLAG( S ) ) THEN
                            WRITE( RDEV,93030 ) 
     &                          PNAME, IFIP( S ), ASC7( S ), ASC3( S ),
     &                          PRJFAC( S ), PFAC
                        ENDIF

                        PRJFAC( S ) = PFAC

124                 CONTINUE

                END IF	! if fip is 0, or nn000, or not.

            ELSE	! using all of ID7,ID3 for SCC

                IF ( FIP .EQ. 0 ) THEN      !  all states.

                    DO  125 S = 1, NASRC

                        IF( ID7 .NE. ASC7( S ) .OR. 
     &                      ID3 .NE. ASC3( S )      ) THEN
                            GO TO 125  ! n/a control

                        ELSEIF ( PFLAG .AND. INVYR( S ) .NE. YEAR ) THEN
                            GO TO 125

                        ELSEIF ( PFLAG .AND. TFLAG( S ) ) THEN
                            WRITE( RDEV,93040 ) 
     &                          PNAME, YEAR, DYR, IFIP( S ), ASC7( S ),
     &                          ASC3( S ), PRJFAC( S ), PFAC

                        ELSEIF ( TFLAG( S ) ) THEN
                            WRITE( RDEV,93030 ) 
     &                          PNAME, IFIP( S ), ASC7( S ), ASC3( S ),
     &                          PRJFAC( S ), PFAC
                        ENDIF

                        PRJFAC( S ) = PFAC

125                 CONTINUE

                ELSE IF ( MOD( FIP, 1000 ) .EQ. 0 ) THEN

                    SID = FIP / 1000        !  state code only
                    I   = FIND1( SID, NSID, SIDS )

                    IF ( I .EQ. 0 ) GO TO 121

                    DO  126 S = ISID( I ), ISID( I+1 ) -1

                        IF( ID7 .NE. ASC7( S ) .OR. 
     &                      ID3 .NE. ASC3( S )      ) THEN
                            GO TO 126  ! n/a control

                        ELSEIF ( PFLAG .AND. INVYR( S ) .NE. YEAR ) THEN
                            GO TO 126

                        ELSEIF ( PFLAG .AND. TFLAG( S ) ) THEN
                            WRITE( RDEV,93040 ) 
     &                          PNAME, YEAR, DYR, IFIP( S ), ASC7( S ),
     &                          ASC3( S ), PRJFAC( S ), PFAC

                        ELSEIF ( TFLAG( S ) ) THEN
                            WRITE( RDEV,93030 ) 
     &                          PNAME, IFIP( S ), ASC7( S ), ASC3( S ),
     &                          PRJFAC( S ), PFAC
                        ENDIF

                        PRJFAC( S ) = PFAC

126                 CONTINUE

                ELSE	!  fip-ASC specific control

                    S = FIND3( FIP, ID7, ID3, NASRC, IFIP, ASC7, ASC3 )

                    IF ( S .GT. 0 ) THEN
                         
                        IF ( PFLAG .AND. INVYR( S ) .NE. YEAR ) THEN
                            GO TO 121

                        ELSEIF ( PFLAG .AND. TFLAG( S ) ) THEN
                            WRITE( RDEV,93040 ) 
     &                          PNAME, YEAR, DYR, IFIP( S ), ASC7( S ),
     &                          ASC3( S ), PRJFAC( S ), PFAC

                        ELSEIF ( TFLAG( S ) ) THEN
                            WRITE( RDEV,93030 ) 
     &                          PNAME, IFIP( S ), ASC7( S ), ASC3( S ),
     &                          PRJFAC( S ), PFAC
                        ENDIF

                        PRJFAC( S ) = PFAC

                    END IF		! if s > 0

                END IF	! if fip is 0, or nn000, or not.

            END IF	!  if using ID4 or ID7
                   
        GO TO  121      !  end loop reading lines in PROJECT AMS packet

129     CONTINUE        !  exit from / skip around loop 121


C............................................................................
C.......   Process the RULE PENETRATION packet of control file

        PNAME = '/RULE PENETRATION/'
        IF( 0 .EQ. INDEX1( PNAME, NPACKETS, PNAMES ) ) GO TO 139

        IF ( .NOT.  GETPKT( PNAME, CDEV ) ) THEN
            CALL M3EXIT( 'CTLAMAT', 0, 0,
     &                   'Could not find ' // PNAME // ' in ' // CNAME, 
     &                   2 )
        END IF

131     CONTINUE

            READ( CDEV,93000,IOSTAT=IOS,END=139 ) LINE
            IREC = IREC + 1
            IF( IOS .NE. 0 ) THEN
                WRITE( MESG,94910 ) 
     &              'Error reading',
     &              CNAME( 1:TRIMLEN( CNAME ) ), 
     &              IREC, IOS
                CALL M3EXIT( 'CTLAMAT', 0, 0, MESG, 1 )
            END IF

            CALL UPCASE( LINE )
            IF( 0 .NE. INDEX( LINE, 'END' ) )  GO TO 139

            FIP  = STR2INT ( LINE(  1: 5 ) )
            CNOX = 0.01 * STR2REAL( LINE( 46:50 ) )
            CVOC = 0.01 * STR2REAL( LINE( 51:55 ) )
            CCO  = 0.01 * STR2REAL( LINE( 56:60 ) )
            L    = TRIMLEN( LINE( 31:40 ) ) - LBLANK( LINE( 31:40 ) )
            IF ( L .EQ. 4 ) THEN
                ID4  = STR2INT ( LINE( 31:40 ) )
                ID7  = ID4 * 1000
                ID3  = 0
            ELSE
                ID7  = STR2INT ( LINE( 31:37 ) )
                ID3  = STR2INT ( LINE( 38:40 ) )
                ID4  = ID7 / 1000
            END IF


            IF ( FIP .EQ. 0 ) THEN      !  all states.

                DO  132 S = 1, NASRC

                    IF( ID7 .NE. ASC7( S ) ) GO TO 132  ! n/a control
                    IF( ID3 .NE. ASC3( S ) ) GO TO 132  ! n/a control
                     
                    IF ( TFLAG( S ) ) 
     &                  WRITE( RDEV,93030 ) 
     &                      PNAME, IFIP( S ), ASC7( S ), ASC3( S ),
     &                      RULPEN( S,1 ), CCO,
     &                      RULPEN( S,2 ), CNOX,
     &                      RULPEN( S,3 ), CVOC

                    RULPEN( S,1 ) = CCO
                    RULPEN( S,2 ) = CNOX
                    RULPEN( S,3 ) = CVOC

132             CONTINUE

            ELSE IF ( MOD( FIP, 1000 ) .EQ. 0 ) THEN

                SID = FIP / 1000        !  state code only
                I   = FIND1( SID, NSID, SIDS )

                IF ( I .EQ. 0 ) GO TO 131

                DO  133 S = ISID( I ), ISID( I+1 ) -1

                    IF( ID7 .NE. ASC7( S ) ) GO TO 133  ! n/a control
                    IF( ID3 .NE. ASC3( S ) ) GO TO 133  ! n/a control

                    IF ( TFLAG( S ) ) 
     &                  WRITE( RDEV,93030 ) 
     &                      PNAME, IFIP( S ), ASC7( S ), ASC3( S ),
     &                      RULPEN( S,1 ), CCO,
     &                      RULPEN( S,2 ), CNOX,
     &                      RULPEN( S,3 ), CVOC

                    RULPEN( S,1 ) = CCO
                    RULPEN( S,2 ) = CNOX
                    RULPEN( S,3 ) = CVOC

133             CONTINUE

            ELSE IF ( MOD( FIP, 1000 ) .EQ. 0 ) THEN

                I   = FIND1( FIP, NFIP, FFIP )

                IF ( I .EQ. 0 ) GO TO 131

                DO  134 S = CFIP( I ), CFIP( I+1 ) - 1

                    IF ( TFLAG( S ) ) 
     &                  WRITE( RDEV,93030 ) 
     &                      PNAME, IFIP( S ), ASC7( S ), ASC3( S ),
     &                      RULPEN( S,1 ), CCO,
     &                      RULPEN( S,2 ), CNOX,
     &                      RULPEN( S,3 ), CVOC

                    RULPEN( S,1 ) = CCO
                    RULPEN( S,2 ) = CNOX
                    RULPEN( S,3 ) = CVOC

134             CONTINUE

            ELSE	!  fip-ASC specific control

                S = FIND3( FIP, ID7, ID3, NASRC, IFIP, ASC7, ASC3 )

                IF ( S .GT. 0 ) THEN

                    IF ( TFLAG( S ) ) 
     &                  WRITE( RDEV,93030 ) 
     &                      PNAME, IFIP( S ), ASC7( S ), ASC3( S ),
     &                      RULPEN( S,1 ), CCO,
     &                      RULPEN( S,2 ), CNOX,
     &                      RULPEN( S,3 ), CVOC

                    RULPEN( S,1 ) = CCO
                    RULPEN( S,2 ) = CNOX
                    RULPEN( S,3 ) = CVOC

                END IF		! if s > 0

            END IF	! if fip is 0, or nn000, or not.

            GO TO  131  !  end loop reading lines in MACT packet

139     CONTINUE  !  end loop reading lines in RULE PENETRATION packet


C............................................................................
C.......   Process the RULE EFFECTIVENESS packet of control file

        PNAME = '/RULE EFFECTIVENESS/'
        IF( 0 .EQ. INDEX1( PNAME, NPACKETS, PNAMES ) ) GO TO 149

        IF ( .NOT.  GETPKT( PNAME, CDEV ) ) THEN
            CALL M3EXIT( 'CTLAMAT', 0, 0,
     &                   'Could not find ' // PNAME // ' in ' // CNAME, 
     &                   2 )
        END IF

141     CONTINUE

            READ( CDEV,93000,IOSTAT=IOS,END=149 ) LINE
            IREC = IREC + 1
            IF( IOS .NE. 0 ) THEN
                WRITE( MESG,94910 ) 
     &              'Error reading',
     &              CNAME( 1:TRIMLEN( CNAME ) ), 
     &              IREC, IOS
                CALL M3EXIT( 'CTLAMAT', 0, 0, MESG, 1 )
            END IF

            CALL UPCASE( LINE )
            IF( 0 .NE. INDEX( LINE, 'END' ) )  GO TO 149

            FIP  = STR2INT ( LINE(  1: 5 ) )
            CNOX = 0.01 * STR2REAL( LINE( 46:50 ) )
            CVOC = 0.01 * STR2REAL( LINE( 51:55 ) )
            CCO  = 0.01 * STR2REAL( LINE( 56:60 ) )
            L    = TRIMLEN( LINE( 31:40 ) ) - LBLANK( LINE( 31:40 ) )
            IF ( L .EQ. 4 ) THEN
                ID4  = STR2INT ( LINE( 31:40 ) )
                ID7  = ID4 * 1000
                ID3  = 0
            ELSE
                ID7  = STR2INT ( LINE( 31:37 ) )
                ID3  = STR2INT ( LINE( 38:40 ) )
                ID4  = ID7 / 1000
            END IF

            IF ( FIP .EQ. 0 ) THEN      !  all states.

                DO  142 S = 1, NASRC

                    IF( ID7 .NE. ASC7( S ) ) GO TO 142  ! n/a control
                    IF( ID3 .NE. ASC3( S ) ) GO TO 142  ! n/a control
                     
                    IF ( TFLAG( S ) ) 
     &                  WRITE( RDEV,93030 ) 
     &                      PNAME, IFIP( S ), ASC7( S ), ASC3( S ),
     &                      RULEFF( S,1 ), CCO,
     &                      RULEFF( S,2 ), CNOX,
     &                      RULEFF( S,3 ), CVOC

                    RULEFF( S,1 ) = CCO
                    RULEFF( S,2 ) = CNOX
                    RULEFF( S,3 ) = CVOC

142             CONTINUE

            ELSE IF ( MOD( FIP, 1000 ) .EQ. 0 ) THEN

                SID = FIP / 1000        !  state code only
                I   = FIND1( SID, NSID, SIDS )

                IF ( I .EQ. 0 ) GO TO 141

                DO  143 S = ISID( I ), ISID( I+1 ) -1

                    IF( ID7 .NE. ASC7( S ) ) GO TO 143  ! n/a control
                    IF( ID3 .NE. ASC3( S ) ) GO TO 143  ! n/a control

                    IF ( TFLAG( S ) ) 
     &                  WRITE( RDEV,93030 ) 
     &                      PNAME, IFIP( S ), ASC7( S ), ASC3( S ),
     &                      RULEFF( S,1 ), CCO,
     &                      RULEFF( S,2 ), CNOX,
     &                      RULEFF( S,3 ), CVOC

                    RULEFF( S,1 ) = CCO
                    RULEFF( S,2 ) = CNOX
                    RULEFF( S,3 ) = CVOC

143             CONTINUE

            ELSE IF ( MOD( FIP, 1000 ) .EQ. 0 ) THEN

                I   = FIND1( FIP, NFIP, FFIP )

                IF ( I .EQ. 0 ) GO TO 141

                DO  144 S = CFIP( I ), CFIP( I+1 ) - 1

                    IF ( TFLAG( S ) ) 
     &                  WRITE( RDEV,93030 ) 
     &                      PNAME, IFIP( S ), ASC7( S ), ASC3( S ),
     &                      RULEFF( S,1 ), CCO,
     &                      RULEFF( S,2 ), CNOX,
     &                      RULEFF( S,3 ), CVOC

                    RULEFF( S,1 ) = CCO
                    RULEFF( S,2 ) = CNOX
                    RULEFF( S,3 ) = CVOC

144             CONTINUE

            ELSE	!  fip-ASC specific control

                S = FIND3( FIP, ID7, ID3, NASRC, IFIP, ASC7, ASC3 )

                IF ( S .GT. 0 ) THEN

                    IF ( TFLAG( S ) ) 
     &                  WRITE( RDEV,93030 ) 
     &                      PNAME, IFIP( S ), ASC7( S ), ASC3( S ),
     &                      RULEFF( S,1 ), CCO,
     &                      RULEFF( S,2 ), CNOX,
     &                      RULEFF( S,3 ), CVOC

                    RULEFF( S,1 ) = CCO
                    RULEFF( S,2 ) = CNOX
                    RULEFF( S,3 ) = CVOC

                END IF		! if s > 0

            END IF	! if fip is 0, or nn000, or not.

            GO TO  141  !  end loop reading lines

149     CONTINUE  !  end loop reading lines in RULE EFFECTIVENESS packet


C............................................................................
C.......   Process the CONTROL EFFICIENCY packet of control file

        PNAME = '/CONTROL EFFICIENCY/'
        IF( 0 .EQ. INDEX1( PNAME, NPACKETS, PNAMES ) ) GO TO 159

        IF ( .NOT.  GETPKT( PNAME, CDEV ) ) THEN
            CALL M3EXIT( 'CTLAMAT', 0, 0,
     &                   'Could not find ' // PNAME // ' in ' // CNAME, 
     &                   2 )
        END IF

151     CONTINUE

            READ( CDEV,93000,IOSTAT=IOS,END=159 ) LINE
            IREC = IREC + 1
            IF( IOS .NE. 0 ) THEN
                WRITE( MESG,94910 ) 
     &              'Error reading',
     &              CNAME( 1:TRIMLEN( CNAME ) ), 
     &              IREC, IOS
                CALL M3EXIT( 'CTLAMAT', 0, 0, MESG, 1 )
            END IF

            CALL UPCASE( LINE )
            IF( 0 .NE. INDEX( LINE, 'END' ) )  GO TO 159

            FIP  = STR2INT ( LINE(  1: 5 ) )
            CNOX = 0.01 * STR2REAL( LINE( 46:50 ) )
            CVOC = 0.01 * STR2REAL( LINE( 51:55 ) )
            CCO  = 0.01 * STR2REAL( LINE( 56:60 ) )
            L    = TRIMLEN( LINE( 31:40 ) ) - LBLANK( LINE( 31:40 ) )
            IF ( L .EQ. 4 ) THEN
                ID4  = STR2INT ( LINE( 31:40 ) )
                ID7  = ID4 * 1000
                ID3  = 0
            ELSE
                ID7  = STR2INT ( LINE( 31:37 ) )
                ID3  = STR2INT ( LINE( 38:40 ) )
                ID4  = ID7 / 1000
            END IF

            IF ( FIP .EQ. 0 ) THEN      !  all states.

                DO  152 S = 1, NASRC

                    IF( ID7 .NE. ASC7( S ) ) GO TO 152  ! n/a control
                    IF( ID3 .NE. ASC3( S ) ) GO TO 152  ! n/a control
                     
                    IF ( TFLAG( S ) ) 
     &                  WRITE( RDEV,93030 ) 
     &                      PNAME, IFIP( S ), ASC7( S ), ASC3( S ),
     &                      CTLEFF( S,1 ), CCO,
     &                      CTLEFF( S,2 ), CNOX,
     &                      CTLEFF( S,3 ), CVOC

                    CTLEFF( S,1 ) = CCO
                    CTLEFF( S,2 ) = CNOX
                    CTLEFF( S,3 ) = CVOC

152             CONTINUE

            ELSE IF ( MOD( FIP, 1000 ) .EQ. 0 ) THEN

                SID = FIP / 1000        !  state code only
                I   = FIND1( SID, NSID, SIDS )

                IF ( I .EQ. 0 ) GO TO 151

                DO  153 S = ISID( I ), ISID( I+1 ) -1

                    IF( ID7 .NE. ASC7( S ) ) GO TO 153  ! n/a control
                    IF( ID3 .NE. ASC3( S ) ) GO TO 153  ! n/a control

                    IF ( TFLAG( S ) ) 
     &                  WRITE( RDEV,93030 ) 
     &                      PNAME, IFIP( S ), ASC7( S ), ASC3( S ),
     &                      CTLEFF( S,1 ), CCO,
     &                      CTLEFF( S,2 ), CNOX,
     &                      CTLEFF( S,3 ), CVOC

                    CTLEFF( S,1 ) = CCO
                    CTLEFF( S,2 ) = CNOX
                    CTLEFF( S,3 ) = CVOC

153             CONTINUE

            ELSE IF ( MOD( FIP, 1000 ) .EQ. 0 ) THEN

                I   = FIND1( FIP, NFIP, FFIP )

                IF ( I .EQ. 0 ) GO TO 151

                DO  154 S = CFIP( I ), CFIP( I+1 ) - 1

                    IF ( TFLAG( S ) ) 
     &                  WRITE( RDEV,93030 ) 
     &                      PNAME, IFIP( S ), ASC7( S ), ASC3( S ),
     &                      CTLEFF( S,1 ), CCO,
     &                      CTLEFF( S,2 ), CNOX,
     &                      CTLEFF( S,3 ), CVOC

                    CTLEFF( S,1 ) = CCO
                    CTLEFF( S,2 ) = CNOX
                    CTLEFF( S,3 ) = CVOC

154             CONTINUE

            ELSE	!  fip-ASC specific control

                S = FIND3( FIP, ID7, ID3, NASRC, IFIP, ASC7, ASC3 )

                IF ( S .GT. 0 ) THEN

                    IF ( TFLAG( S ) ) 
     &                  WRITE( RDEV,93030 ) 
     &                      PNAME, IFIP( S ), ASC7( S ), ASC3( S ),
     &                      CTLEFF( S,1 ), CCO,
     &                      CTLEFF( S,2 ), CNOX,
     &                      CTLEFF( S,3 ), CVOC

                    CTLEFF( S,1 ) = CCO
                    CTLEFF( S,2 ) = CNOX
                    CTLEFF( S,3 ) = CVOC

                END IF		! if s > 0

            END IF	! if fip is 0, or nn000, or not.

            GO TO  151  !  end loop reading lines

159     CONTINUE  !  end loop reading lines in CONTROL EFFICIENCY packet


C............................................................................
C.......   Process RACT packet of control file:

        PNAME = '/RACT/'
        IF( 0 .EQ. INDEX1( PNAME, NPACKETS, PNAMES ) ) GO TO 169

        IF ( .NOT.  GETPKT( PNAME, CDEV ) ) THEN
            CALL M3EXIT( 'CTLAMAT', 0, 0,
     &                   'Could not find ' // PNAME // ' in ' // CNAME, 
     &                   2 )
        END IF

161     CONTINUE        !  head of loop reading lines in the RACT packet

            READ( CDEV,93000,IOSTAT=IOS,END=169 ) LINE
            IREC = IREC + 1
            IF( IOS .NE. 0 ) THEN
                WRITE( MESG,94910 ) 
     &              'Error reading',
     &              CNAME( 1:TRIMLEN( CNAME ) ), 
     &              IREC, IOS
                CALL M3EXIT( 'CTLAMAT', 0, 0, MESG, 1 )
            END IF

            CALL UPCASE( LINE )
            IF( 0 .NE. INDEX( LINE, 'END' ) )  GO TO 169

            FIP  = STR2INT ( LINE(  1: 5 ) )
            CNOX = 1.0 - STR2REAL( LINE( 26:30 ) )
            CVOC = 1.0 - STR2REAL( LINE( 32:36 ) )
            L    = TRIMLEN ( LINE( 16:25 ) ) - LBLANK( LINE( 16: 25 ) )
            IF ( L .EQ. 4 ) THEN
                ID4  = STR2INT ( LINE( 16:25 ) )
                ID7  = ID4 * 1000
                ID3  = 0
            ELSE
                ID7  = STR2INT ( LINE( 16:22 ) )
                ID3  = STR2INT ( LINE( 23:25 ) )
                ID4  = ID7 / 1000
            END IF


            IF ( FIP .EQ. 0 ) THEN      !  all states.
                 
            DO  162 S = 1, NASRC

                    IF( ID7 .NE. ASC7( S ) ) GO TO 162  ! n/a control
                    IF( ID3 .NE. ASC3( S ) ) GO TO 162  ! n/a control

                    IF ( TFLAG( S ) ) 
     &                  WRITE( RDEV,93030 ) 
     &                      PNAME, IFIP( S ), ASC7( S ), ASC3( S ),
     &                      CTLEFF( S,2 ), CNOX,
     &                      CTLEFF( S,3 ), CVOC

                    CTLEFF( S,2 ) = CNOX
                    CTLEFF( S,3 ) = CVOC

162             CONTINUE

            ELSE IF ( MOD( FIP, 1000 ) .EQ. 0 ) THEN
                
                SID = FIP / 1000        !  state code only
                I = FIND1( SID, NSID, SIDS )
                IF ( I .EQ. 0 ) GO TO 161
                
                DO  163 S = ISID( I ), ISID( I+1 ) -1

                    IF( ID7 .NE. ASC7 ( S ) ) GO TO 163  ! n/a control
                    IF( ID3 .NE. ASC3 ( S ) ) GO TO 163  ! n/a control

                    IF ( TFLAG( S ) ) 
     &                  WRITE( RDEV,93030 ) 
     &                      PNAME, IFIP( S ), ASC7( S ), ASC3( S ),
     &                      CTLEFF( S,2 ), CNOX,
     &                      CTLEFF( S,3 ), CVOC

                    CTLEFF( S,2 ) = CNOX
                    CTLEFF( S,3 ) = CVOC

163             CONTINUE

            ELSE	!  fip-ASC specific control
            
                S = FIND3( FIP, ID7, ID3, NASRC, IFIP, ASC7, ASC3 )

                IF ( S .GT. 0 ) THEN

                    IF ( TFLAG( S ) ) 
     &                  WRITE( RDEV,93030 ) 
     &                      PNAME, IFIP( S ), ASC7( S ), ASC3( S ),
     &                      CTLEFF( S,2 ), CNOX,
     &                      CTLEFF( S,3 ), CVOC

                    CTLEFF( S,2 ) = CNOX
                    CTLEFF( S,3 ) = CVOC

                END IF		! if s > 0

            END IF	! if fip is 0, or nn000, or not.

            GO TO  161  !  end loop on lines in this RACT packet

169     CONTINUE        !  exit from / skip around loop 161


C............................................................................
C.......   Process MACT packet of control file:

        PNAME = '/MACT/'
        IF( 0 .EQ. INDEX1( PNAME, NPACKETS, PNAMES ) ) GO TO 179

        IF ( .NOT.  GETPKT( PNAME, CDEV ) ) THEN
            CALL M3EXIT( 'CTLAMAT', 0, 0,
     &                   'Could not find ' // PNAME // ' in ' // CNAME, 
     &                   2 )
        END IF

171     CONTINUE        !  head of loop reading lines in the MACT packet

            READ( CDEV,93000,IOSTAT=IOS,END=179 ) LINE
            IREC = IREC + 1
            IF( IOS .NE. 0 ) THEN
                WRITE( MESG,94910 ) 
     &              'Error reading',
     &              CNAME( 1:TRIMLEN( CNAME ) ), 
     &              IREC, IOS
                CALL M3EXIT( 'CTLAMAT', 0, 0, MESG, 1 )
            END IF

            CALL UPCASE( LINE )
            IF( 0 .NE. INDEX( LINE, 'END' ) )  GO TO 179

            FIP  = STR2INT ( LINE(  1: 5 ) )
            CVOC = 1.0 - STR2REAL( LINE( 26:30 ) )
            L    = TRIMLEN ( LINE( 16:25 ) ) - LBLANK( LINE( 16:25 ) )
            IF ( L .EQ. 4 ) THEN
                ID4  = STR2INT ( LINE( 16:25 ) )
                ID7  = ID4 * 1000
                ID3  = 0
            ELSE
                ID7  = STR2INT ( LINE( 16:22 ) )
                ID3  = STR2INT ( LINE( 23:25 ) )
                ID4  = ID7 / 1000
            END IF

            IF ( FIP .EQ. 0 ) THEN      !  all states.

                DO  172 S = 1, NASRC

                    IF( ID7 .NE. ASC7( S ) ) GO TO 172  ! n/a control
                    IF( ID3 .NE. ASC3( S ) ) GO TO 172  ! n/a control
                     
                    IF ( TFLAG( S ) ) 
     &                  WRITE( RDEV,93030 ) 
     &                      PNAME, IFIP( S ), ASC7( S ), ASC3( S ),
     &                      CTLEFF( S,3 ), CVOC

                    CTLEFF( S,3 ) = CVOC

172             CONTINUE

            ELSE IF ( MOD( FIP, 1000 ) .EQ. 0 ) THEN

                SID = FIP / 1000        !  state code only
                I   = FIND1( SID, NSID, SIDS )

                IF ( I .EQ. 0 ) GO TO 171

                DO  173 S = ISID( I ), ISID( I+1 ) -1

                    IF( ID7 .NE. ASC7( S ) ) GO TO 173  ! n/a control
                    IF( ID3 .NE. ASC3( S ) ) GO TO 173  ! n/a control

                    IF ( TFLAG( S ) ) 
     &                  WRITE( RDEV,93030 ) 
     &                      PNAME, IFIP( S ), ASC7( S ), ASC3( S ),
     &                      CTLEFF( S,3 ), CVOC

                    CTLEFF( S,3 ) = CVOC

173             CONTINUE

            ELSE	!  fip-ASC specific control

                S = FIND3( FIP, ID7, ID3, NASRC, IFIP, ASC7, ASC3 )

                IF ( S .GT. 0 ) THEN

                    IF ( TFLAG( S ) ) 
     &                  WRITE( RDEV,93030 ) 
     &                      PNAME, IFIP( S ), ASC7( S ), ASC3( S ),
     &                      CTLEFF( S,3 ), CVOC

                    CTLEFF( S,3 ) = CVOC

                END IF		! if s > 0

            END IF	! if fip is 0, or nn000, or not.

            GO TO  171  !  end loop reading lines in MACT packet

179     CONTINUE        !  exit from / skip around loop 171


C............................................................................
C.......   Process CTG packet of control file:
        
        PNAME = '/CTG/'
        IF( 0 .EQ. INDEX1( PNAME, NPACKETS, PNAMES ) ) GO TO 189

        IF ( .NOT.  GETPKT( PNAME, CDEV ) ) THEN
            CALL M3EXIT( 'CTLAMAT', 0, 0,
     &                   'Could not find ' // PNAME // ' in ' // CNAME, 
     &                   2 )
        END IF

        CALL M3MSG2( 'Reading CONTROL PACKETS file...' )

181     CONTINUE        !  head of loop reading this category of packet

            READ( CDEV,93000,IOSTAT=IOS,END=189 ) LINE
            IREC = IREC + 1
            IF( IOS .NE. 0 ) THEN
                WRITE( MESG,94910 ) 
     &              'Error reading',
     &              CNAME( 1:TRIMLEN( CNAME ) ), 
     &              IREC, IOS
                CALL M3EXIT( 'CTLAMAT', 0, 0, MESG, 1 )
            END IF
                   
            CALL UPCASE( LINE )
            IF( 0 .NE. INDEX( LINE, 'END' ) )  GO TO 189
            
            FIP  = STR2INT ( LINE(  1: 5 ) )
            CNOX = 1.0 - STR2REAL( LINE( 27:31 ) )
            CVOC = 1.0 - STR2REAL( LINE( 33:37 ) )
            TNOX = STR2REAL( LINE( 39:48 ) )
            TVOC = STR2REAL( LINE( 50:59 ) )
            L    = TRIMLEN ( LINE( 16:25 ) ) - LBLANK( LINE( 16:25 ) )
            IF ( L .EQ. 4 ) THEN
                ID4  = STR2INT ( LINE( 16:25 ) )
                ID7  = ID4 * 1000
                ID3  = 0
            ELSE
                ID7  = STR2INT ( LINE( 16:22 ) )
                ID3  = STR2INT ( LINE( 23:25 ) )
                ID4  = ID7 / 1000
            END IF

            IF ( FIP .EQ. 0 ) THEN      !  all states.
                
                DO  182 S = 1, NASRC

                    IF( ID7 .NE. ASC7( S ) ) GO TO 182  ! n/a control
                    IF( ID3 .NE. ASC3( S ) ) GO TO 182  ! n/a control

                    IF ( NOX( S ) .GT. TNOX ) THEN

                        IF( VOC( S ) .GT. TVOC ) THEN   ! NOX,VOC both apply

                            IF ( TFLAG( S ) ) 
     &                          WRITE( RDEV,93030 ) 
     &                              PNAME, 
     &                              IFIP( S ), ASC7( S ), ASC3( S ),
     &                              CTLEFF( S,2 ), CNOX,
     &                              CTLEFF( S,3 ), CVOC

                            CTLEFF( S,2 ) = CNOX
                            CTLEFF( S,3 ) = CVOC

                        ELSE                            ! NOX control applies

                            IF ( TFLAG( S ) ) 
     &                          WRITE( RDEV,93030 ) 
     &                              PNAME, 
     &                              IFIP( S ), ASC7( S ), ASC3( S ),
     &                              CTLEFF( S,2 ), CNOX

                            CTLEFF( S,2 ) = CNOX
                            
                        END IF

                    ELSE IF( VOC( S ) .GT. TVOC ) THEN! VOC control applies
                                                      ! T = RULEFF(S)*RULPEN(S)

                        IF ( TFLAG( S ) ) 
     &                      WRITE( RDEV,93030 ) 
     &                          PNAME, IFIP( S ), ASC7( S ), ASC3( S ),
     &                          CTLEFF( S,3 ), CVOC
           
                        CTLEFF( S,3 ) = CVOC
     
                    END IF      !  if strongest

182             CONTINUE	!  end loop on sources S, state-ind ctls

            ELSE IF ( MOD( FIP, 1000 ) .EQ. 0 ) THEN
                
                SID = FIP / 1000        !  state code only
                I   = FIND1( SID, NSID, SIDS )

                IF ( I .EQ. 0 )  GO TO  181
                
                DO  183 S = ISID( I ), ISID( I+1 ) -1

                    IF( ID7 .NE. ASC7 ( S ) ) GO TO 183  ! n/a control
                    IF( ID3 .NE. ASC3 ( S ) ) GO TO 183  ! n/a control

                    IF ( NOX( S ) .GT. TNOX ) THEN
                        IF( VOC( S ) .GT. TVOC ) THEN   ! NOX,VOC both apply
                             
                            IF ( TFLAG( S ) ) 
     &                          WRITE( RDEV,93030 ) 
     &                              PNAME, 
     &                              IFIP( S ), ASC7( S ), ASC3( S ),
     &                              CTLEFF( S,2 ), CNOX,
     &                              CTLEFF( S,3 ), CVOC

                            CTLEFF( S,2 ) = CNOX
                            CTLEFF( S,3 ) = CVOC

                        ELSE                            ! NOX control applies

                            IF ( TFLAG( S ) ) 
     &                          WRITE( RDEV,93030 ) 
     &                              PNAME, 
     &                              IFIP( S ), ASC7( S ), ASC3( S ),
     &                              CTLEFF( S,2 ), CNOX

                            CTLEFF( S,2 ) = CNOX

                        END IF

                    ELSE IF( VOC( S ) .GT. TVOC ) THEN! VOC control applies
                                                      ! T = RULEFF(S)*RULPEN(S)

                        IF ( TFLAG( S ) ) 
     &                      WRITE( RDEV,93030 ) 
     &                          PNAME, IFIP( S ), ASC7( S ), ASC3( S ),
     &                          CTLEFF( S,3 ), CVOC

                        CTLEFF( S,3 ) = CVOC

                    END IF      !  if strongest

183             CONTINUE	!  end loop on sources S, state code CTG ctls

            ELSE	!  fip-ASC specific control
            
                S = FIND3( FIP, ID7, ID3, NASRC, IFIP, ASC7, ASC3 )

                IF ( S .GT. 0 ) THEN
                
                    IF ( NOX( S ) .GT. TNOX ) THEN

                        IF( VOC( S ) .GT. TVOC ) THEN   ! NOX,VOC both apply

                            IF ( TFLAG( S ) ) 
     &                          WRITE( RDEV,93030 ) 
     &                              PNAME, 
     &                              IFIP( S ), ASC7( S ), ASC3( S ),
     &                              CTLEFF( S,2 ), CNOX,
     &                              CTLEFF( S,3 ), CVOC

                            CTLEFF( S,2 ) = CNOX
                            CTLEFF( S,3 ) = CVOC

                        ELSE                            ! NOX control applies

                            IF ( TFLAG( S ) ) 
     &                          WRITE( RDEV,93030 ) 
     &                              PNAME, IFIP( S ), 
     &                              ASC7( S ), ASC3( S ),
     &                              CTLEFF( S,2 ), CNOX

                            CTLEFF( S,2 ) = CNOX

                        END IF

                    ELSE IF( VOC( S ) .GT. TVOC ) THEN! VOC control applies
                                                      ! T = RULEFF(S)*RULPEN(S)

                        IF ( TFLAG( S ) ) 
     &                      WRITE( RDEV,93030 ) 
     &                          PNAME, IFIP( S ), ASC7( S ), ASC3( S ),
     &                          CTLEFF( S,3 ), CVOC

                        CTLEFF( S,3 ) = CVOC

                    END IF      !  if strongest
         
                END IF		! if s > 0:  source-specific controls

            END IF	! if fip is 0, or nn000, or not.

            GO TO  181  !  end loop on lines in this CTG packet

189     CONTINUE        !  exit from / skip around loop 181


C............................................................................
C.......   Compute the control matrix:
 
        DO  222  I = 1, NIPOL
        DO  211  S = 1, NASRC
 
            CMAT( S,I ) = PRJFAC( S ) * CMAT( S,I ) * 
     &          ( 1.0 - CTLEFF( S,I ) * RULEFF( S,I ) * RULPEN( S,I ) )
 
211     CONTINUE
222     CONTINUE
 
C............................................................................
C.......   Write out the control matrix:

        CALL M3MSG2( 'Writing out CONTROL MATRIX...' )

        IF ( .NOT. WRITE3( MNAME, 'ALL', 0, 0,CMAT ) ) THEN
            CALL M3EXIT( 'CTLAMAT', 0, 0, 
     &                   'Error writing CONTROL MATRIX file.', 2 )
        END IF



999   CONTINUE          !  exit program

      CALL M3EXIT( 'CTLAMAT', 0, 0, 
     &             'Normal completion of program CTLAMAT', 2 )

C******************  FORMAT  STATEMENTS   ******************************

C...........   Informational (LOG) message formats... 92xxx

92000   FORMAT( 5X, A )


C...........   Formatted file I/O formats............ 93xxx

93000   FORMAT( A )

93010   FORMAT( A16 )

93020   FORMAT( I5, 1X, I7, I3 )

93030   FORMAT( A, ':', 1X, I5.5, 1X, I7.7, I3.3,
     &          3( F8.5, '~>', F8.5, : '; ' )  )

93040   FORMAT( A, ':', I4, ' to ', I4, 1X, I5.5, 1X, I7.7, I3.3,
     &          3( F8.5, '~>', F8.5, : '; ' )  )

C...........   Internal buffering formats............ 94xxx

94910   FORMAT( A, ' "', A, :,  '" at line', I7, :, '; IOSTAT=', I7 )


        END

