
         PROGRAM EMSAREA

C***********************************************************************
C  program body starts at line 186
C
C  DESCRIPTION:
C       Construct Models3/EDSS area source file from data contained in
C       EPS/AMS-style area source file.
C
C  PRECONDITIONS REQUIRED:
C       Sorted, cut-down input data for area sources.
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C       Models-3 I/O
C       PROMPTFFILE, PROMPTMFILE, STR2INT, STR2REAL
C
C  REVISION  HISTORY:
C       Prototype adapted from rawarea.F 8/96 by CJC 
C	Version 9/96 by CJC now produces ASCS actual-SCC file.
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
C***************************************************************************/

        IMPLICIT NONE

        INCLUDE 'ARDIMS3.EXT'   !  area-source dimensioning parameters
        INCLUDE 'CHDIMS3.EXT'   !  emis chem parms (inventory + model)
        INCLUDE 'GRDIMS3.EXT'   !  grid parameters
        INCLUDE 'TMDIMS3.EXT'   !  time parameters
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
        INCLUDE 'FDESC3.EXT'    !  I/O API file description data structures.

C...........   EXTERNAL FUNCTIONS and their descriptions:

        CHARACTER*2     CRLF
        LOGICAL         DSCGRID
        LOGICAL         ENVYN
        INTEGER         FIND1
        INTEGER         GETNUM
        LOGICAL         GETYN
        INTEGER         INDEX1
        INTEGER         JULIAN
        INTEGER         JUNIT
        INTEGER         LBLANK
        INTEGER         PROMPTFFILE
        CHARACTER*16    PROMPTMFILE
        INTEGER         SECSDIFF
        INTEGER         STR2INT
        REAL            STR2REAL
        INTEGER         TRIMLEN

        EXTERNAL CRLF, DSCGRID, ENVYN, FIND1, GETNUM, GETYN, INDEX1,
     &           JULIAN, JUNIT, LBLANK, PROMPTFFILE, PROMPTMFILE, 
     &           SECSDIFF, STR2INT, STR2REAL, TRIMLEN


C...........   PARAMETER:

        CHARACTER*5     BLANK5
        INTEGER         NSMAX
        PARAMETER     ( BLANK5 = ' ',
     &                  NSMAX = NASRC * NIPOL )


C...........   LOCAL VARIABLES and their descriptions:
C...........   NOTE that ASC (area-source-category) ID's are 10-digit
C...........   unsigned integers which may be treated as a leading 7-digit
C...........   field, and a trailing 3-digit field.  *7 and *3 arrays below
C...........   follow this scheme with parallel arrays

        INTEGER     TZONE       !  time zone
        INTEGER     FIP
        INTEGER     ID7
        INTEGER     ID3
        INTEGER     IPOL
        INTEGER     TPF
        INTEGER     INY
        REAL        CEFF
        REAL        REFF
        REAL        RPEN
        REAL        EMIS
        REAL        VMISS
        INTEGER     LFIP, LID7, LID3

C...............   Time zone tables:  FIP-independent; state-only; state-county

        INTEGER         TZONE0

        INTEGER         NZS
        INTEGER         INDXSA( NASID )
        INTEGER         TZONSA( NASID )
        INTEGER         TFIPSA( NASID )
        INTEGER         TZONST( NASID )
        INTEGER         TFIPST( NASID )

        INTEGER         NZF
        INTEGER         INDXFA( NAFIP )
        INTEGER         TZONFA( NAFIP )
        INTEGER         TFIPFA( NAFIP )
        INTEGER         TZONEF( NAFIP )
        INTEGER         TFIPEF( NAFIP )

C...........   Area Sources Table input unsorted copy (*A); 

        INTEGER     NSRC             !  record-count
        INTEGER     INDEXA( NSMAX )  !  subscript table
        INTEGER     IFIPSA( NSMAX )  !  source FIPS (county) ID
        INTEGER     ISCC7A( NSMAX )  !  leading-7  digits of source ASC
        INTEGER     ISCC3A( NSMAX )  !  trailing-3 digits of source ASC
        INTEGER     ICODEA( NSMAX )  !  inventory species
        INTEGER     TPFLGA( NSMAX )  !  applicability of temporal profile types
        INTEGER     INVYRA( NSMAX )  !  year inventory taken
        REAL        CTLEFA( NSMAX )  !  control-factor (efficiency)
        REAL        RULEFA( NSMAX )  !  control-rule effectiveness
        REAL        RULPEA( NSMAX )  !  control-rule pentration
        REAL        EMISVA( NSMAX )  !  emissions values.

C...........   sorted final version (image EMISREC of data record)
C.......   Common EMISREC holds an entire output record.  
C.......   Order of arrays in EMISREC _must_ match order of
C.......   variables in the output file.

        INTEGER     NAREA            !  current source-count
        INTEGER     IFIPS ( NASRC )  !  source FIPS (county) ID
        INTEGER     ISCC7 ( NASRC )  !  leading-7  digits of source ASC
        INTEGER     ISCC3 ( NASRC )  !  trailing-3 digits of source ASC
        INTEGER     ZONES ( NASRC )  !  time zones
        INTEGER     TPFLAG( NASRC )  !  applicability of temporal profile types
        INTEGER     INVYR ( NASRC )  !  year inventory taken
        REAL        CTLEFF( NASRC, NIPOL )  !  control-factor (efficiency)
        REAL        RULEFF( NASRC, NIPOL )  !  control-rule effectiveness
        REAL        RULPEN( NASRC, NIPOL )  !  control-rule pentration
        REAL        EMISV ( NASRC, NIPOL )  !  emissions values.

        COMMON  / EMISREC / IFIPS,  ISCC7,  ISCC3,  ZONES, INVYR, 
     &                      TPFLAG, CTLEFF, RULEFF, RULPEN, EMISV

C...........   File unit numbers and logical names

        INTEGER         ADEV             !  for area-source file
        INTEGER         IDEV             !  Raw input file list file
        INTEGER         LDEV             !  log-device
        INTEGER         SDEV             !  for actual ASC file
        INTEGER         ZDEV             !  for time-zone file
        CHARACTER*16    ENAME   !  logical name for emission source output file
        CHARACTER*256   RNAME   !  phyical name for emission source input file

C...........   Other local variables

        INTEGER         C, S, I, J, K, V !  loop counters.
        INTEGER         IOS              !  I/O status
        INTEGER         IREC             !  input line (record) number
        INTEGER         L1, L2           !  start/end position of input LINE
        LOGICAL         VFLAG   !  input verification:  TRUE iff ERROR
        LOGICAL         MFLAG   !  input verification:  report missing species
        LOGICAL         DFLAG   !  input verification:  report duplicate records
        LOGICAL         AFLAG   !  accumulator: report missing species
        REAL            TMPFAC  !  factor for annualizing emissions
        CHARACTER*2     TMPAA   !  temporal basis 'AA' or 'AD'
        CHARACTER*16    SCRBUF  !  scratch name buffer
        CHARACTER*16    COORDN  !  coordinate system name
        CHARACTER*8     PNAME   !  pollutant name
        CHARACTER*240   LINE    !  input line from AREA file
        CHARACTER*256   MESG    !  for M3EXIT() output

C***********************************************************************
C   begin body of program EMSAREA

        LDEV = INIT3()

        WRITE( *,92000 ) 
     &  ' ',
     &  'Program EMSAREA to take multiple EMS-95 area source files',
     &  'and produce the AREA SOURCE EMISSIONS VECTOR file. ',
     &  ' ',
     &  'You will need to enter the logical names for the input and',
     &  'output files (and to have set them prior to program launch,',
     &  'using "setenv <logicalname> <pathname>").',
     &  'Optional checking that all species are reported for each ',
     &  'source may be turned on via "setenv RAW_SRC_CHECK Y".',
     &  ' ',
     &  'You may use END_OF-FILE (control-D) to quit the program',
     &  'during logical-name entry. Default responses are indicated',
     &  'in brackets [LIKE THIS].',
     &  ' '

        IF ( .NOT. GETYN( 'Continue with program?', .TRUE. ) ) THEN
            CALL M3EXIT( 'EMSAREA', 0, 0, 'Ending program', 2 )
        END IF

        MFLAG = ENVYN( 'RAW_MISS_CHECK', 
     &                 'EMSAREA check for missing species-records',
     &                 .FALSE.,
     &                 IOS )

        DFLAG = ENVYN( 'RAW_DUP_CHECK', 
     &                 'EMSAREA check for duplicate species-records',
     &                 .FALSE.,
     &                 IOS )


C.......   Get the coordinate system name and parameters
C.......   to put into file header:

        IF ( .NOT. DSCGRID( GRDNM, COORDN, GDTYP3D, 
     &              P_ALP3D, P_BET3D, P_GAM3D, XCENT3D, YCENT3D,
     &              XORIG3D, YORIG3D, XCELL3D, YCELL3D, 
     &              NCOLS3D, NROWS3D, NTHIK3D ) ) THEN

            SCRBUF = GRDNM
            MESG   = '"' // SCRBUF( 1:TRIMLEN( SCRBUF ) ) //
     &               '" not found in GRIDDESC file'
            CALL M3EXIT( 'EMSAREA', 0, 0, MESG, 2 )

        END IF


C.......   Get file name; open input time zone file

        ZDEV = PROMPTFFILE( 
     &           'Enter logical name for TIME ZONE file',
     &           .TRUE., .TRUE., 'ZONES', 'EMSAREA' )


        VFLAG  = .FALSE.
        TZONE0 = 5      !  default:  EST
        NZS    = 0
        NZF    = 0
        IREC   = 0
11      CONTINUE

            READ( ZDEV, *, END=12, IOSTAT=IOS ) FIP, TZONE
            IREC = IREC + 1

            IF ( IOS .NE. 0 ) THEN
                WRITE( MESG,94010 ) 
     &              'Unit number', ZDEV, 
     &              'I/O Status ', IOS, 
     &              'Line number', IREC
                CALL M3MSG2( MESG )
                CALL M3EXIT( 'EMSPOINT', 0, 0, 
     &              'Error reading TIME ZONE file.', 2 )
            END IF

            IF ( FIP .EQ. 0 ) THEN              !  fallback -- all sources

                TZONE0 = TZONE

            ELSE IF ( MOD( FIP, 1000 ) .EQ. 0 ) THEN     !  state-specific zone

                NZS = NZS + 1
                IF ( NZS .LE. NASID ) THEN
                    INDXSA( NZS ) = NZS
                    TFIPSA( NZS ) = FIP / 1000
                    TZONSA( NZS ) = TZONE
                END IF

            ELSE                                        !  county-specific zone

                NZF = NZF + 1
                IF ( NZF .LE. NAFIP ) THEN
                    INDXFA( NZF ) = NZF
                    TFIPFA( NZF ) = FIP
                    TZONFA( NZF ) = TZONE
                END IF

            END IF      !  if fip zero, or nn000, or not.

            GO TO  11

12      CONTINUE        !  exit from loop reading ZDEV

        CLOSE( ZDEV )
        
        IF ( NZS .GT. NASID .OR. NZF .GT. NAFIP ) THEN
            WRITE( MESG,94010 )
     &          'State-specific time zone table dim:', NASID,
     &          'actual:', NZS
            CALL M3MSG2( MESG )
            WRITE( MESG,94010 )
     &          'County-specific time zone table dim:', NAFIP,
     &          'actual:', NZF
            CALL M3MSG2( MESG )
            CALL M3EXIT( 'EMSPOINT', 0, 0, 
     &                   'OVERFLOW reading TIME ZONE file.', 2 )
        END IF
        

        CALL SORTI1( NZS, INDXSA, TFIPSA )
        DO  22  I = 1, NZS
            J = INDXSA( I )
            TZONST( I ) = TZONSA( J )
            TFIPST( I ) = TFIPSA( J )
22      CONTINUE

        CALL SORTI1( NZF, INDXFA, TFIPFA )
        DO  33  I = 1, NZF
            J = INDXFA( I )
            TZONEF( I ) = TZONFA( J )
            TFIPEF( I ) = TFIPFA( J )
33      CONTINUE


C.......   OTAG initializations (doesn't use REFF, RPEN):

        REFF = 1.0
        RPEN = 1.0

C.......   Get file name for opening input raw area source file
        IDEV = PROMPTFFILE(
     &         'Enter the name of the RAW DATA FILENAME LISTING',
     &          .TRUE., .TRUE., 'ANLST', 'EMSAREA' )

        CALL M3MSG2( 'Reading RAW AREA SOURCE files...' )

        S   = 0
        INY = IMISS3
101     CONTINUE

C.............  Read file names
C.............  Exit to line 141 if read is EOF
            READ( IDEV, 93000, END=141  ) LINE

            L1 = INDEX( LINE, 'INVYEAR' )
            L2 = TRIMLEN( LINE )

            IF( L1 .GT. 0 ) THEN

                INY = STR2INT( LINE( L1+7:L2 ) )
                GO TO 101                       ! To head of ANLST read loop

            ELSEIF( INY .LE. 0 ) THEN

                CALL M3EXIT( 'EMSAREA', 0, 0, 
     &               'Must set inventory year using INVYEAR packet ' //
     &               'in ANLST file', 2 )

            ELSE

                RNAME = LINE( 1:L2 )   ! Set inventory file name

            ENDIF

            ADEV   = JUNIT()
            OPEN( ADEV, ERR=1006, FILE=RNAME, STATUS='OLD' )

            WRITE( MESG,94010 ) 
     &             'Successful OPEN using year', INY, 
     &             'for inventory file:' // CRLF() // BLANK5 //
     &             RNAME( 1:TRIMLEN( RNAME ) )
            CALL M3MSG2( MESG )

C.......   Factors for annualizing the emissions from daily totals
C.......   (corrected for leap-year)

            IF (      MOD( INY,   4 ) .NE. 0 ) THEN
                TMPFAC = 365.0
            ELSE IF ( MOD( INY, 100 ) .NE. 0 ) THEN
                TMPFAC = 366.0
            ELSE IF ( MOD( INY, 400 ) .NE. 0 ) THEN
                TMPFAC = 365.0
            ELSE
                TMPFAC = 366.0
            END IF

C.......   Process this area source file

            IREC = 0

111         CONTINUE        !  head of the ADEV-read loop

                READ( ADEV, 93000, END=133, IOSTAT=IOS ) LINE

                IF ( IOS .NE. 0 ) THEN
                    WRITE( MESG,94010 ) 
     &                  'Unit number', ADEV, 
     &                  'I/O Status ', IOS, 
     &                  'Line number', IREC
                    CALL M3MSG2( MESG )
                    CALL M3EXIT( 'EMSAREA', 0, 0, 
     &                  'Error reading "' //
     &                  RNAME( 1:TRIMLEN( RNAME ) ) // '"', 2 )
                END IF

                IREC = IREC + 1
                
                C     = 21 + LBLANK( LINE( 21:25 ) )
                PNAME = LINE( C : 25 )

                IPOL = INDEX1( PNAME, NIPOL, EINAM )

                IF ( IPOL .LE. 0 ) THEN
                    WRITE( MESG,94010 ) 
     &                  'Bad line', IREC, 
     &                  'Pollutant code "' // 
     &                  PNAME( 1:TRIMLEN( PNAME ) ) // 
     &                  '" in "' //
     &                  RNAME( 1:TRIMLEN( RNAME ) ) // '"'
                    CALL M3MESG( MESG )
                    GO TO  111
                END IF

                EMIS = STR2REAL( LINE( 52:65 ) )
                IF ( EMIS .LT. 0.0 )  THEN
                    WRITE( MESG,94010 ) 
     &                  'Bad line', IREC, 
     &                  'Emis value "' // LINE( 162:171 ) // 
     &                  '" in "' //
     &                  RNAME( 1:TRIMLEN( RNAME ) ) // '"'
                    CALL M3MESG( MESG )
                    EMIS  =  0.0
                    VFLAG = .TRUE.
                    GO TO  111
                END IF

                FIP  = 1000 * STR2INT( LINE( 1:2 ) ) + 
     &                        STR2INT( LINE( 3:5 ) )

                C    = LBLANK( LINE( 6:20 ) )
                ID7  = STR2INT( LINE(  6+C : 12+C ) )
                ID3  = STR2INT( LINE( 13+C : 15+C ) )
                
                CEFF = STR2REAL( LINE( 88 : 94 ) )
                IF ( CEFF .LT. 0 ) CEFF = 0.0

C-otag          REFF = STR2REAL( LINE( 86:91 ) )
C-otag          IF ( REFF .LT. 0 ) REFF = 100.0
C
C-otag          RPEN = STR2REAL( LINE( 93:97 ) )
C-otag          IF ( RPEN .LT. 0 ) RPEN = 100.0

                TMPAA = LINE( 95 : 96 )
                CALL UPCASE( TMPAA ) 
                
                IF      ( TMPAA .EQ. 'AD' ) THEN
                    EMIS = TMPFAC * EMIS
                    TPF  = WDTPFAC
                ELSE IF ( TMPAA .NE. 'AA' ) THEN
                    MESG = 'Bad temporal basis "' // TMPAA // '"'
                    CALL M3MESG( MESG )
                    VFLAG = .TRUE.
                    GO TO  111
                END IF

                S = S + 1
                IF ( S .LE. NSMAX ) THEN

                    INDEXA( S ) = S !  index table for later SORTI3()
                    IFIPSA( S ) = FIP  
                    ISCC7A( S ) = ID7 
                    ISCC3A( S ) = ID3
                    ICODEA( S ) = IPOL
                    TPFLGA( S ) = TPF
                    INVYRA( S ) = INY

                    EMISVA( S ) = EMIS
                    CTLEFA( S ) = CEFF
                    RULEFA( S ) = REFF
                    RULPEA( S ) = RPEN

                END IF              !  if S in bounds

                GO TO  111          !  to head of ADEV-read loop

133         CONTINUE        !  end of the ADEV-read loop
            
            CLOSE( ADEV )

            IF ( VFLAG ) THEN
                CALL M3EXIT( 'EMSAREA', 0, 0, 
     &                       'Error reading AREA SOURCE file.', 2 )
            END IF      !  if VFLAG

            GO TO 101   !  to process next input file

141     CONTINUE        !  all input files now processed

        NSRC = S
        IF ( S .GT. NSMAX ) THEN
            WRITE( MESG,94010 ) 
     &       'Actual record-count        :', NSRC, CRLF() // BLANK5 //
     &       'Max    record-count (NSMAX):', NSMAX
            CALL M3MSG2( MESG )
            CALL M3EXIT( 'EMSAREA', 0, 0, 
     &        'Max record-count exceeded in AREA SOURCE file.', 2 )

        END IF


C.......   Use SORTI3() to perform an indirect sort by FIPS,SCC7,SCC3;
C.......   then permute the records according to the result:

        CALL M3MSG2( 'Sorting RAW AREA SOURCE data...' )

        CALL SORTI3( NSRC, INDEXA, IFIPSA, ISCC7A, ISCC3A )

        IF ( MFLAG ) THEN
            VMISS = BADVAL3
        ELSE
            VMISS = 0.0
        END IF
        VFLAG = .FALSE.
        AFLAG = .FALSE.
        LFIP  = IMISS3
        LID7  = IMISS3
        LID3  = IMISS3
        J     = 0

        DO  144  S = 1, NSRC

            I   = INDEXA( S )
            FIP = IFIPSA( I )
            ID7 = ISCC7A( I )
            ID3 = ISCC3A( I )

            IF ( FIP .NE. LFIP  .OR.
     &           ID7 .NE. LID7  .OR.
     &           ID3 .NE. LID3 ) THEN           !  if new source encountered

                J = J + 1
                LFIP = FIP
                LID7 = ID7
                LID3 = ID3
                
                IF ( J .LE. NASRC ) THEN	!  if J in bounds

                    IFIPS ( J ) = FIP
                    ISCC7 ( J ) = ID7
                    ISCC3 ( J ) = ID3
                    TPFLAG( J ) = TPFLGA( I )
                    INVYR ( J ) = INVYRA( I )

                    K = FIND1( FIP, NZF, TFIPEF )
                    IF ( K .GT. 0 ) THEN
                        ZONES( J ) = TZONEF( K )
                    ELSE
                        K = FIND1( FIP/1000, NZS, TFIPST )
                        IF ( K .GT. 0 ) THEN
                            ZONES( J ) = TZONST( K )
                        ELSE
                            ZONES( J ) = TZONE0
                        END IF
                    END IF

                    IF ( MFLAG .AND. J .GT. 1 ) THEN !  verify last source and
						     !  initialize this source
                        DO  142  V = 1, NIPOL

                            IF ( EMISV ( J-1,V ) .LE. AMISS3 ) THEN
                                AFLAG = .TRUE.
                                WRITE( MESG,94020 )
     &                          'Missing emissions record:  FIP:', FIP,
     &                          'ASC:', ID7, ID3, 'Species:', EINAM( V )
                                CALL M3MESG( MESG )
                            END IF

                            EMISV ( J,V ) = BADVAL3
                            CTLEFF( J,V ) = 0.0
                            RULPEN( J,V ) = 0.0
                            RULEFF( J,V ) = 0.0

142                     CONTINUE

                    ELSE			!  initialize this source only

                        DO  143  V = 1, NIPOL

                            EMISV ( J,V ) = VMISS
                            CTLEFF( J,V ) = 0.0
                            RULPEN( J,V ) = 0.0
                            RULEFF( J,V ) = 0.0

143                     CONTINUE

                    END IF

                    V = ICODEA( I )

                    EMISV ( J,V ) = EMISVA( I )
                    CTLEFF( J,V ) = CTLEFA( I )
                    RULEFF( J,V ) = RULEFA( I )
                    RULPEN( J,V ) = RULPEA( I )

                END IF      !  if source number J in bounds

            ELSE IF ( J .LE. NASRC ) THEN

                V = ICODEA( I )

                IF( EMISV ( J,V ) .LE. 0.0 ) THEN

                    EMISV ( J,V ) = EMISVA( I )
                    CTLEFF( J,V ) = CTLEFA( I )
                    RULEFF( J,V ) = RULEFA( I )
                    RULPEN( J,V ) = RULPEA( I )

                ELSE IF ( EMISVA( I ) .GT. 0.0 ) THEN

                    VFLAG = .TRUE.
                    EMISV ( J,V ) = EMISV( J,V ) + EMISVA( I )
                    WRITE( MESG,94020 )
     &                  'Duplicate emissions record:  FIP:', FIP,
     &                  'ASC:', ID7, ID3, 
     &                   'Species:', EINAM( V )
                    CALL M3MESG( MESG )

                END IF

            END IF      !  if new source encountered; else if J in bounds

144     CONTINUE

        NAREA = J
        IF ( NAREA .NE. NASRC ) THEN
            WRITE( MESG,94010 ) 
     &          'Actual=', NAREA, ', Dimensioned=', NASRC,       
     &          'numbers of sources do not match!'
            CALL M3EXIT( 'EMSAREA', 0, 0, MESG, 2 )            
        END IF

        IF ( VFLAG ) THEN
            MESG = 'Duplicate records found in input file'
            IF ( DFLAG ) THEN
                CALL M3EXIT( 'EMSAREA', 0, 0, MESG, 2 )
            ELSE
                CALL M3WARN( 'EMSAREA', 0, 0, MESG )
            END IF
        END IF

        IF ( AFLAG ) THEN
            CALL M3EXIT( 'EMSAREA', 0, 0,
     &                   'Missing species recs in input file', 2 )
        END IF          !  if sflag and missing-record error


C........   Get file name; open and construct output actual-ASC file

        SDEV = PROMPTFFILE( 
     &      'Enter logical name for ACTUAL-ASC output file',
     &      .FALSE., .TRUE., 'ASCS', 'EMSAREA' )

        CALL SORTI2( NSRC, INDEXA, ISCC7A, ISCC3A )

        CALL M3MSG2( 'Writing ACTUAL-ASC output file...' )

        LID7 = -1
        LID3 = -1
        K    =  0

        DO  155  S = 1, NSRC
            J = INDEXA( S )
            ID7 = ISCC7A( J )
            ID3 = ISCC3A( J )
            IF ( ID7 .NE. LID7  .OR.  ID3 .NE. LID3 ) THEN
                LID7 = ID7
                LID3 = ID3
                K    = K + 1
                WRITE( SDEV,93010 ) ID7, ID3
            END IF		!  if new ASC pair encountered
155     CONTINUE

        WRITE( MESG,94010 )
     &      'Number of source ASC codes:  dim max', NASCC,
     &      'actual', K
        CALL M3MSG2( MESG )


C.......   Get file name; open output area sources file
C.......   Note that coordinate system definition was put into 
C.......   FDESC3.EXT data structures by DSCGRID().

        FTYPE3D = GRDDED3
        SDATE3D = 0 !  n/a
        STIME3D = 0 !  n/a
        TSTEP3D = 0             !  time independent
        NVARS3D = 4 * NIPOL + 6
        NCOLS3D = 1
        NROWS3D = NASRC     !  number of rows = # of area sources.
        NLAYS3D = 1
        NTHIK3D = 1
        XCELL3D = DBLE( NCOLS - 1 ) * XCELL3D	!  size of acceptance window
        YCELL3D = DBLE( NROWS - 1 ) * YCELL3D	!  size of acceptance window
        VGTYP3D = IMISS3
        VGTOP3D = AMISS3
        GDNAM3D = COORDN
        FDESC3D( 1 ) = 'Annual AREA SOURCE emissions values.'
        FDESC3D( 2 ) = 'Inventory  type "' // LINE(  1: 1 ) // '"'
        FDESC3D( 3 ) = 'Projection year "' // LINE(  2: 3 ) // '"'
        FDESC3D( 4 ) = 'Base year       "' // LINE(  5: 6 ) // '"'
        FDESC3D( 5 ) = 'Emission   type "' // LINE(  7: 8 ) // '"'
        FDESC3D( 6 ) = 'Period     type "' // LINE( 37:38 ) // '"'
        FDESC3D( 7 ) = 'Beginning  Date "' // LINE( 40:47 ) // '"'
        FDESC3D( 8 ) = 'Ending     Date "' // LINE( 49:56 ) // '"'

        J = 1
        VNAME3D( J ) = 'FIP'
        UNITS3D( J ) = 'n/a'
        VDESC3D( J ) = 'FIP code for counties'
        VTYPE3D( J ) = M3INT
        J = J + 1

        VNAME3D( J ) = 'ASC7'
        UNITS3D( J ) = 'n/a'
        VDESC3D( J ) = 'Area Source Category code digits 1-7'
        VTYPE3D( J ) = M3INT
        J = J + 1

        VNAME3D( J ) = 'ASC3'
        UNITS3D( J ) = 'n/a'
        VDESC3D( J ) = 'Area Source Category code digits 8-10'
        VTYPE3D( J ) = M3INT
        J = J + 1

        VNAME3D( J ) = 'ZONES'
        UNITS3D( J ) = 'hours from GMT'
        VDESC3D( J ) = 'Std. time zones (0 for GMT, 5 for Eastern)'
        VTYPE3D( J ) = M3INT
        J = J + 1

        VNAME3D( J ) = 'INVYR'
        UNITS3D( J ) = 'year AD'
        VDESC3D( J ) = 'Year of inventory for this record'
        VTYPE3D( J ) = M3INT
        J = J + 1

        VNAME3D( J ) = 'TPFLAG'
        UNITS3D( J ) = 'T|2? T|3?'
        VDESC3D( J ) = 'Use week(2), month(3) temporal profiles or not'
        VTYPE3D( J ) = M3INT
        J = J + 1

        DO  211  I = 1, NIPOL
            VNAME3D( J ) = 'CTLEFF_' // EINAM( I )
            UNITS3D( J ) = 'n/a'
            VDESC3D( J ) = 
     &      'control efficiency (in [0,100], or "MISSING": < -9.0E36)'
            VTYPE3D( J ) = M3REAL
            J = J + 1
211     CONTINUE

        DO  212  I = 1, NIPOL
            VNAME3D( J ) = 'RULEFF_' // EINAM( I )
            UNITS3D( J ) = 'n/a'
            VDESC3D( J ) = 
     &      'Rule effectiveness (in [0,100], or "MISSING": < -9.0E36)'
            VTYPE3D( J ) = M3REAL
            J = J + 1
212     CONTINUE

        DO  213  I = 1, NIPOL
            VNAME3D( J ) = 'RULPEN_' // EINAM( I )
            UNITS3D( J ) = 'n/a'
            VDESC3D( J ) = 
     &      'Rule penetration (in [0,100], or "MISSING": < -9.0E3J)'
            VTYPE3D( J ) = M3REAL
            J = J + 1
213     CONTINUE

        DO  214  I = 1, NIPOL
            VNAME3D( J ) = EINAM( I )
            UNITS3D( J ) = 'tons/year'
            VDESC3D( J ) = LINE( 1:1 ) // LINE( 37:38 ) 
     &                                 // ' emissions totals'
            VTYPE3D( J ) = M3REAL
            J = J + 1
214     CONTINUE

        ENAME = PROMPTMFILE( 
     &          'Enter logical name for AREA SOURCE output file',
     &          FSUNKN3, 'AREA', 'EMSAREA' )


C.......   Write out the area source emissions values:

        CALL M3MSG2( 'Writing out AREA SOURCES file...' )

        IF ( .NOT. WRITE3( ENAME, 'ALL', 0, 0,  IFIPS ) ) THEN
            CALL M3EXIT( 'EMSAREA', 0, 0,  
     &                   'Error writing AREA OUTPUT file' , 2 )
        END IF


999     CONTINUE          !  exit program

        CALL M3EXIT( 'EMSAREA', 0, 0,
     &               'Normal completion of Program EMSAREA', 0 )

1006    MESG = 'Error opening file ' // RNAME( 1:TRIMLEN( RNAME ) )
        CALL M3EXIT( 'EMSAREA', 0, 0, MESG, 2 )

C******************  FORMAT  STATEMENTS   ******************************

C...........   Informational (LOG) message formats... 92xxx

92000   FORMAT( 5X, A )

92010   FORMAT( 5X, A, :, I10 )


C...........   Formatted file I/O formats............ 93xxx

93000   FORMAT( A )

93010   FORMAT( I7.7, I3.3 )


C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I7, :, 1X ) )

94020   FORMAT( A, 1X, I5, 1X, A, I7.7, I3.3, 1X, A, A )

94040   FORMAT( A, I2.2 )

94041   FORMAT( A, I3.3 )

94042   FORMAT( A, I4.4 )

        END

