
         PROGRAM EMSPOINT

C***********************************************************************
C  program body starts at line 332
C
C  DESCRIPTION:
C       Sort data contained in EMS-95 raw point source files and 
C       construct Models3/EDSS point source file
C
C  PRECONDITIONS REQUIRED:
C       EMS-95 Input data for point sources.  Accepts multiple 
C       emission.pt, facility.pt, stack.pt, process.pt, and device.pt 
C       files.
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C       Models-3 I/O
C       FIND1, FIND2, FINDC, GETYN, PROMPTFFILE, PROMPTMFILE, SORTIC, TRIMLEN
C       FIXSTK
C
C  REVISION  HISTORY:
C       Prototype  5/96 by MR Houyoux (based on V1.3 rawpoint.F)
C       Optimizations 8/96 by CJC
C
C****************************************************************************/
C
C Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
C                System
C File: @(#)$Id$
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
C Pathname: $Source$
C Last updated: $Date$ 
C
C***************************************************************************

      IMPLICIT NONE

C...........   INCLUDES:

        INCLUDE 'PTDIMS3.EXT'   !  point-source dimensioning parameters
        INCLUDE 'GRDIMS3.EXT'   !  grid dimensioning parameters
        INCLUDE 'CHDIMS3.EXT'   !  emis chem info (both inventory and model)
        INCLUDE 'TMDIMS3.EXT'   !  time related parameters
        INCLUDE 'CONST3.EXT'    !  physical constants
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
        INCLUDE 'FDESC3.EXT'    !  I/O API file description data structures.


C...........   PARAMETER

        CHARACTER*5     BLANK5
        REAL            FTOC    ! Farenheit to Celcius ratio
        REAL            FT2M    ! Feet to meters conversion
        REAL            FLWE2M  ! Flow in cubic ft per min to cubic m per sec

        PARAMETER     ( BLANK5  = ' ',
     &                  FTOC   = 5./9., 
     &                  FT2M   = 0.3048,
     &                  FLWE2M = FT2M * FT2M * FT2M / 60.0 )

C...........   EXTERNAL FUNCTIONS and their descriptions:
        
        CHARACTER*2     CRLF
        LOGICAL         DSCGRID
        LOGICAL         ENVYN
        INTEGER         FIND1, FIND2, FINDC   !  returns -1 for failure
        LOGICAL         GETYN
        INTEGER         INDEX1
        INTEGER         JULIAN
        INTEGER         JUNIT
        LOGICAL         LAMBERT
        INTEGER         LBLANK
        LOGICAL         LL2LAM
        LOGICAL         UTM2LAM
        INTEGER         PROMPTFFILE
        CHARACTER*16    PROMPTMFILE
        INTEGER         SECSDIFF
        INTEGER         STR2INT
        REAL            STR2REAL
        INTEGER         TRIMLEN
        REAL            YR2DAY  

        EXTERNAL CRLF, DSCGRID, ENVYN, FIND1, FIND2, FINDC, 
     &           GETYN, INDEX1, JULIAN, JUNIT, LAMBERT, LBLANK, 
     &           LL2LAM, UTM2LAM, PROMPTFFILE, PROMPTMFILE, SECSDIFF, 
     &           STR2INT, STR2REAL, TRIMLEN, YR2DAY


C...........   LOCAL VARIABLES and their descriptions:
C...............   Time zone tables:  FIP-independent; state-only; state-county

        INTEGER         TZONE0

        INTEGER         NZS
        INTEGER         INDXSA( NPSID )
        INTEGER         TZONSA( NPSID )
        INTEGER         TFIPSA( NPSID )
        INTEGER         TZONST( NPSID )
        INTEGER         TFIPST( NPSID )

        INTEGER         NZF
        INTEGER         INDXFA( NPFIP )
        INTEGER         TZONFA( NPFIP )
        INTEGER         TFIPFA( NPFIP )
        INTEGER         TZONEF( NPFIP )
        INTEGER         TFIPEF( NPFIP )

C...........   Point sources table for facility.pt files (unsorted; sorted)

        INTEGER        FS                ! Counter & number of facility records
        INTEGER        INDXFCA( NPFMAX ) ! subscript array for sort
        INTEGER        FCZONEA( NPFMAX ) ! UTM zone from facility files
        REAL           FCCRDXA( NPFMAX ) ! Facility X coordinates
        REAL           FCCRDYA( NPFMAX ) ! Facility Y coordinates
        CHARACTER*20   FCKEYA ( NPFMAX ) ! FIP // facility ID
        CHARACTER*40   FCDESCA( NPFMAX ) ! facility name from facility files

        INTEGER        FCZONE( NPFMAX )  ! UTM zone from facility files
        REAL           FCCRDX( NPFMAX )  ! Facility X coordinates
        REAL           FCCRDY( NPFMAX )  ! Facility Y coordinates
        CHARACTER*20   FCKEY ( NPFMAX )  ! FIP // facility ID
        CHARACTER*40   FCDESC( NPFMAX )  ! facility name from facility files

C...........   Point sources table for stack.pt files  (unsorted; sorted)

        INTEGER        SS                ! Counter & number of stack records
        INTEGER        INDXSKA( NPSMAX ) ! subscript table for sort
        INTEGER        IFCKEYA( NPSMAX ) ! subscript table for sort
        REAL           SKDIAMA( NPSMAX ) ! stack diameter from stack files [m]
        REAL           SKHEITA( NPSMAX ) ! stack height from stack files [m]
        REAL           SKTEMPA( NPSMAX ) ! stack exit temp from stack files [C]
        REAL           SKCRDXA( NPSMAX ) ! stack UTM easting coord [m]
        REAL           SKCRDYA( NPSMAX ) ! stack UTM northing coord [m]
        REAL           SKVELOA( NPSMAX ) ! stack exit velocity [m/sec]
        REAL           SKFLOWA( NPSMAX ) ! stack flow when not match v [m^3/sec]
        CHARACTER*32   SKKEYA ( NPSMAX ) ! FIP // facility ID // stack ID

        INTEGER        IFCKEY( NPSMAX )  ! subscript table for sort
        REAL           SKDIAM( NPSMAX )  ! stack diameter from stack files [m]
        REAL           SKHEIT( NPSMAX )  ! stack height from stack files [m]
        REAL           SKTEMP( NPSMAX )  ! stack exit temp from stack files [C]
        REAL           SKCRDX( NPSMAX )  ! stack UTM easting coord [m]
        REAL           SKCRDY( NPSMAX )  ! stack UTM northing coord [m]
        REAL           SKVELO( NPSMAX )  ! stack exit velocity [m^3/sec]
        REAL           SKFLOW( NPSMAX )  ! new stack exit velocity [m^3/sec]
        CHARACTER*32   SKKEY ( NPSMAX )  ! FIP // facility ID // stack ID

C...........   Point sources table for process.pt files  (unsorted; sorted)

        INTEGER        PS                ! Counter & number of process records
        INTEGER        INDXPA( NPPMAX )  ! subscript table
        INTEGER        PSSCCA( NPPMAX )  ! SCC code from process files
        CHARACTER*56   PSKEYA( NPPMAX )  ! FIP // FCID // SKID // DVID // PRID

        INTEGER        PSSCC ( NPPMAX )  ! SCC code from process files
        CHARACTER*56   PSKEY ( NPPMAX )  ! FIP // FCID // SKID // DVID // PRID

C...........   Point sources table for device.pt files   (unsorted; sorted)

        INTEGER        DS                ! Counter & number of device records
        INTEGER        INDXDVA( NPDMAX ) ! subscript table for sort
        INTEGER        DVSICA ( NPDMAX ) ! SIC code from device files
        INTEGER        DVIWEKA( NPDMAX ) ! Weekly profile code from device files
        INTEGER        DVIDIUA( NPDMAX ) ! Hourly profile code...
        CHARACTER*44   DVKEYA ( NPDMAX ) ! FIP // FCID // SKID // DVID

        INTEGER        DVSIC ( NPDMAX )  ! SIC code from device files
        INTEGER        DVIWEK( NPDMAX )  ! Weekly profile code from device files
        INTEGER        DVIDIU( NPDMAX )  ! Hourly profile code...
        CHARACTER*44   DVKEY ( NPDMAX )  ! FIP // FCID // SKID // DVID


C...........   Point Sources Table input unsorted copy (*A); 
C...........   sorted final version (image EMISREC of data record)

        INTEGER        ES               !  current source-count
        INTEGER        NSRC             !  total source-count
        INTEGER        INDEXA( NPEMAX ) !  subscript table for SORTI4()
        INTEGER        IFIPA ( NPEMAX ) !  source FIPS (county) ID
        INTEGER        ISCCA ( NPEMAX ) !  source SCC
        INTEGER        ISICA ( NPEMAX ) !  source SIC
        INTEGER        IPLANA( NPEMAX ) !  Plant ID
        INTEGER        ISTACA( NPEMAX ) !  Stack ID
        INTEGER        ICODEA( NPEMAX ) !  species-subscript
        INTEGER        TPFLGA( NPEMAX ) !  applicability of temporal profile types
        INTEGER        INVYRA( NPEMAX ) !  inventory year
        INTEGER        IDIUA ( NPEMAX ) !  Hourly profile code for each source
        INTEGER        IWEKA ( NPEMAX ) !  Weekly profile code for each source
        REAL           XLOCAA( NPEMAX ) !  UTM X-location (m)
        REAL           YLOCAA( NPEMAX ) !  UTM Y-location (m)
        REAL           STKHTA( NPEMAX ) !  stack height   (m)
        REAL           STKDMA( NPEMAX ) !  stack diameter (m)
        REAL           STKTKA( NPEMAX ) !  exhaust temperature (deg K)
        REAL           STKVEA( NPEMAX ) !  exhaust velocity    (m/s)
        REAL           STKFLA( NPEMAX ) !  new exhaust velocity    (m/s)
        REAL           RULPEA( NPEMAX ) !  rule penetration   fraction
        REAL           RULEFA( NPEMAX ) !  rule effectiveness fraction
        REAL           CTLEFA( NPEMAX ) !  control efficiency fraction
        REAL           EMISVA( NPEMAX ) !  emissions values (tons/yr)
        CHARACTER*15   EMFCID( NPEMAX ) !  facility IDs
        CHARACTER*12   EMDVID( NPEMAX ) !  device IDs
        CHARACTER*12   EMSKID( NPEMAX ) !  stack IDs
        CHARACTER*40   EMDESC( NPEMAX ) !  facility name from facility files

C.......   Common EMISREC holds an entire output record.  
C.......   Order of arrays in EMISREC _must_ match order of
C.......   variables in the output file.

        INTEGER        NPOINT           !  current source-count
        INTEGER        IFIP  ( NPSRC )  !  source FIPS (county) ID
        INTEGER        ISCC  ( NPSRC )  !  source SCC
        INTEGER        ISIC  ( NPSRC )  !  source SIC
        INTEGER        IPLANT( NPSRC )  !  Plant ID
        INTEGER        ISTACK( NPSRC )  !  Stack ID
        INTEGER        TZONES( NPSRC )  !  time zones
        INTEGER        TPFLAG( NPSRC )  !  temporal profile types
        INTEGER        INVYR ( NPSRC )  !  inventory year for this record
        INTEGER        IDXSCC( NPSRC )  !  subscript table for SORTI1() for SCC
        REAL           XLOCA ( NPSRC )  !  UTM X-location (m)
        REAL           YLOCA ( NPSRC )  !  UTM Y-location (m)
        REAL           STKHT ( NPSRC )  !  stack height   (m)
        REAL           STKDM ( NPSRC )  !  stack diameter (m)
        REAL           STKTK ( NPSRC )  !  exhaust temperature (deg K)
        REAL           STKVE ( NPSRC )  !  exhaust velocity    (m/s)
        REAL           RULPEN( NPSRC, NIPOL )  !  rule penetration   fraction
        REAL           RULEFF( NPSRC, NIPOL )  !  rule effectiveness fraction
        REAL           CTLEFF( NPSRC, NIPOL )  !  control efficiency fraction
        REAL           EMISV ( NPSRC, NIPOL )  !  emissions values (tons/yr)

        COMMON /EMISREC / IFIP, ISIC, ISCC, IPLANT, ISTACK, TZONES,
     &                    TPFLAG, INVYR, XLOCA, YLOCA, 
     &                    STKHT, STKDM, STKTK, STKVE, 
     &                    CTLEFF, RULEFF, RULPEN, EMISV

C...........   File units and logical/physical names

        INTEGER         ADEV    !  Unit number for output actual SCCs
        INTEGER         CDEV    !  Unit number for EMS95/SMOKE facility/stack ID
        INTEGER         EDEV    !  Unit number for emission.pt file
        INTEGER         DDEV    !  Unit number for device.pt
        INTEGER         FDEV    !  Unit number for facility.pt file
        INTEGER         IDEV    !  List of physical raw input file names
        INTEGER         LDEV    !  log-device
        INTEGER         PDEV    !  Unit number for process.pt
        INTEGER         RDEV    !  Unit number for default stack parameters
        INTEGER         SDEV    !  Unit number for stack.pt file
        INTEGER         TDEV    !  Unit number for output temporal profile #s
        INTEGER         ZDEV    !  for time zone file

        CHARACTER*16    ENAME   !  emissions output file logical name
        CHARACTER*256   DVNAME  !  device.pt   file physical name buffer 
        CHARACTER*256   EMNAME  !  emission.pt file physical name buffer 
        CHARACTER*256   FCNAME  !  facility.pt file physical name buffer 
        CHARACTER*256   PSNAME  !  process.pt  file physical name buffer 
        CHARACTER*256   SKNAME  !  stack.pt    file physical name buffer 


C...........   Other local variables

        REAL            STKF, STKH, STKD, STKT, STKV  ! Temporary stack parms

        REAL            CEFF    !  Temporary control effectiveness
        REAL            DAY2YR  !  Local, leap-year-able, DAY to YEAR factor
        REAL            EDROP( NIPOL ) ! Sum of dropped emissions
        REAL            EMIS    !  Temporary emission value
        REAL            FSAV    !  Flow saved
        REAL            REFF    !  Temporary rule effectiveness
        REAL            RPEN    !  Temporary rule penetration
        REAL            VMISS   !  emissons initialization value
        REAL            VNEW    !  New velocity
        REAL            XLOC    !  Temporary X-coordinate
        REAL            YLOC    !  Temporary Y-coordinate
        REAL            XX, YY  !  scratch location variables
                                
        INTEGER         S, I, J, K, L, V !  counters and indices
        INTEGER         K1, K2, K3, K4   !  counters and indices
        INTEGER         L1, L2           !  counters and indices
        INTEGER         FIP,  SCC,  SIC, PLANT, STACK  ! Temporary fip, scc
        INTEGER         LFIP, LSCC, LLSCC, LPLT, LSTK  ! Last fip, scc, etc.

        INTEGER         COD     !  Temporary pollutant code number
        INTEGER         CSS     !  Start of non-blank character string
        INTEGER         IDIU    !  Temporary hourly profile code
        INTEGER         INY     !  Temporary inventory year
        INTEGER         IOS     !  I/O status
        INTEGER         IREC    !  input line (record) number
        INTEGER         IWEK    !  Temporary weekly profile code
        INTEGER         TZONE   !  time zone
        INTEGER         UZONE   !  UTM zone for the file as a whole
        INTEGER         NDROP   !  Count of number of dropped emis records
        INTEGER         NSCC    !  Number of total SCC codes 
        INTEGER         TPF     !  Temporary temporal ID
        INTEGER         ZONE    !  Temporary UTM zone

        LOGICAL         CFLAG   !  velocity recalc: TRUE iff VELOC_RECALC = Y
        LOGICAL         EFLAG   !  input verification:  TRUE iff ERROR
        LOGICAL         RULFLAG !  rule effective file(s): TRUE if exist
        LOGICAL         SFLAG   !  input verification:  report missing species

        CHARACTER*5     CFIP    !  Character FIP code
        CHARACTER*5     CPOL    !  Temporary pollutant code
        CHARACTER*12    DVID    !  Temporary device ID
        CHARACTER*12    PRID    !  Temporary process ID
        CHARACTER*12    SKID    !  Temporary stack ID
        CHARACTER*15    FCID    !  Temporary facility ID
        CHARACTER*16    SCRBUF  !  scratch buffer
        CHARACTER*16    COORDN  !  coordinate system name
        CHARACTER*20    FKEY
        CHARACTER*32    SKEY
        CHARACTER*44    DKEY
        CHARACTER*56    PKEY
        CHARACTER*256   LINE    !  input line from POINT file
        CHARACTER*256   MESG    !  text for M3EXIT()
        CHARACTER*256   NAMTMP  !  Temporary buffer for input file names

C***********************************************************************
C   begin body of program EMSPOINT

        LDEV = INIT3()

        WRITE( *,92000 ) 
     &  ' ',
     &  'Program EMSPOINT to take sequence of EMS-95 point source',
     &  'foundation files and produce the POINT SOURCE EMISSIONS',
     &  'VECTOR file. ',
     &  ' ',
     &  'You will need to enter the logical names for the input and',
     &  'output files (and to have set them prior to program launch,',
     &  'using "setenv <logicalname> <pathname>").',
     &  'Optional checking that all species are reported for each ',
     &  'source may be turned on via "setenv RAW_SRC_CHECK Y".',
     &  ' ',
     &  'You may use END_OF-FILE (control-D) to quit the program',
     &  'during logical-name entry. Default responses are given in',
     &  'brackets [LIKE THIS] and can be accepted by hitting the',
     &  '<RETURN> key.',
     &  ' '

        IF ( .NOT. GETYN( 'Continue with program?', .TRUE. ) ) THEN
            CALL M3EXIT( 'EMSPOINT', 0, 0, 'Ending program.', 2 )
        END IF

        SFLAG = ENVYN( 'RAW_SRC_CHECK', 
     &                 'EMSPOINT check for missing species-records',
     &                 .FALSE.,
     &                 IOS )

        CFLAG = ENVYN( 'VELOC_RECALC', 
     &                 'EMSPOINT setting for recalculating velocity',
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
            CALL M3EXIT( 'EMSPOINT', 0, 0, MESG, 2 )

        END IF          !  if dscgrid() failed

        IF      ( GDTYP3D .EQ. UTMGRD3 ) THEN
            UZONE = NINT( P_ALP3D )
        ELSE IF ( GDTYP3D .EQ. LAMGRD3 ) THEN
            IF ( .NOT. LAMBERT( COORDN, P_ALP3D, P_BET3D, P_GAM3D,
     &                          XCENT3D, YCENT3D ) ) THEN
                MESG = 'Bad grid: "' // GRDNM // '"'
                CALL M3EXIT( 'EMSPOINT', 0, 0, MESG, 2 )
            END IF      !  if lambert() failed
        ELSE IF ( GDTYP3D .EQ. LATGRD3 ) THEN
            !  do nothing -- no setup necessary
        ELSE
            WRITE( MESG, 94010 ) 
     &          'Unsupported coordinate system type', GDTYP3D
            CALL M3EXIT( 'EMSPOINT', 0, 0, MESG, 2 )
        END IF          !  if grid type utm, lambert, lat, or not

C.......   Get file name for opening input raw point source file

        IDEV = PROMPTFFILE( 
     &         'Enter the name of the RAW DATA FILENAME LISTING',
     &          .TRUE., .TRUE., 'PNLST', 'EMSPOINT' )

C.......   Get file name for input replacement stack parameters file

        RDEV = PROMPTFFILE(
     &          'Enter REPLACEMENT STACK PARAMETERS file',
     &          .TRUE., .TRUE., 'PSTK', 'EMSPOINT' )

C.......   Get file name for time zones files

        ZDEV = PROMPTFFILE( 
     &         'Enter logical name for TIME ZONE file',
     &         .TRUE., .TRUE., 'ZONES', 'EMSPOINT' )

C.......   Process input time zone file
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

            ELSE IF ( MOD( FIP, 100 ) .EQ. 0 ) THEN     !  state-specific zone

                NZS = NZS + 1
                IF ( NZS .LE. NPSID ) THEN
                    INDXSA( NZS ) = NZS
                    TFIPSA( NZS ) = FIP / 1000
                    TZONSA( NZS ) = TZONE
                END IF

            ELSE                                        !  county-specific zone

                NZF = NZF + 1
                IF ( NZF .LE. NPFIP ) THEN
                    INDXFA( NZF ) = NZF
                    TFIPFA( NZF ) = FIP
                    TZONFA( NZF ) = TZONE
                END IF

            END IF      !  if fip zero, or nn000, or not.

            GO TO  11

12      CONTINUE        !  exit from loop reading ZDEV

        CLOSE( ZDEV )
        IF ( NZS .GT. NPSID  .OR.  NZF .GT. NPFIP ) THEN
            WRITE( MESG,94010 ) 
     &        'Number of state-only records  :', NZS, CRLF()// BLANK5//
     &        'Number of state&county records:', NZF, CRLF()// BLANK5//
     &        'Dims-- NPSID:', NPSID, 'NPFIP:', NPFIP

            CALL M3MSG2( MESG )
            CALL M3EXIT( 'EMSPOINT', 0, 0,
     &                   'Overflow reading time zone file', 2 )
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


C.......   Read the point source files:
C.......   Initialize variables for keeping track of dropped emissions

        NDROP = 0
        DO 99  COD = 1, NIPOL
            EDROP( COD ) = 0.0
99      CONTINUE

C.........  Loop through all files in the raw point source list
        NSRC =    0
        ES   =    0     ! For sources (emission.pt)
        EFLAG = .FALSE.

101     CONTINUE     !  Head of IDEV read loop

C.............  Read file names, but OK to end only on first from group
C.............  Exit to line 202 if first read is EOF, or error to line 1005
            READ( IDEV, 93000, END=202  ) LINE    ! device.pt file

            L1 = INDEX( LINE, 'INVYEAR' )
            L2 = TRIMLEN( LINE )
 
            IF( L1 .GT. 0 ) THEN
 
                INY = STR2INT( LINE( L1+7:L2 ) )
 
                IF( INY .LE. 0 ) THEN
 
                    CALL M3EXIT( 'EMSPOINT', 0, 0,
     &               'Must set inventory year using INVYEAR packet ' //
     &               'in PNLST file', 2 ) 

                ELSEIF( INY .LT. 1970 ) THEN

                    CALL M3EXIT( 'EMSPOINT', 0, 0,
     &               'INVYEAR packet has set bad 4-digit year' //
     &               'in PNLST file', 2 ) 

                ENDIF

                GO TO 101                       ! To head of PNLST read loop

            ELSE
                DVNAME = LINE( 1:L2 )             ! device.pt file

            ENDIF

            READ( IDEV, 93000, END=1005 ) EMNAME  ! emission.pt file
            IF( INDEX( EMNAME, 'INVYEAR' ) .GT. 0 ) GO TO 1007

            READ( IDEV, 93000, END=1005 ) FCNAME  ! facility.pt file
            IF( INDEX( FCNAME, 'INVYEAR' ) .GT. 0 ) GO TO 1007

            READ( IDEV, 93000, END=1005 ) PSNAME  ! process.pt file
            IF( INDEX( PSNAME, 'INVYEAR' ) .GT. 0 ) GO TO 1007

            READ( IDEV, 93000, END=1005 ) SKNAME  ! stack.pt file
            IF( INDEX( SKNAME, 'INVYEAR' ) .GT. 0 ) GO TO 1007

C.............  Open files, and report status. Use NAMTMP for ease of
C.............  making changes
            DDEV   = JUNIT()
            NAMTMP = DVNAME
            OPEN( DDEV, ERR=1006, FILE=NAMTMP, STATUS='OLD' )

            WRITE( MESG,94010 ) 'Successful OPEN using year', INY,
     &             'for inventory files:' // CRLF() // BLANK5 //
     &             NAMTMP( 1:TRIMLEN( NAMTMP ) )
            CALL M3MSG2( MESG )

            EDEV   = JUNIT()
            NAMTMP = EMNAME
            OPEN( EDEV, ERR=1006, FILE=NAMTMP, STATUS='OLD' )
            CALL M3MSG2( NAMTMP( 1:TRIMLEN( NAMTMP ) ) )

            FDEV   = JUNIT()
            NAMTMP = FCNAME
            OPEN( FDEV, ERR=1006, FILE=NAMTMP, STATUS='OLD' )
            CALL M3MSG2( NAMTMP( 1:TRIMLEN( NAMTMP ) ) )

            PDEV   = JUNIT()
            NAMTMP = PSNAME
            OPEN( PDEV, ERR=1006, FILE=NAMTMP, STATUS='OLD' )
            CALL M3MSG2( NAMTMP( 1:TRIMLEN( NAMTMP ) ) )

            SDEV   = JUNIT()
            NAMTMP = SKNAME
            OPEN( SDEV, ERR=1006, FILE=NAMTMP, STATUS='OLD' )
            CALL M3MSG2( NAMTMP( 1:TRIMLEN( NAMTMP ) ) )

C.............   Calculate DAY2YR factor

            DAY2YR = 1. / YR2DAY( INY )

C........................................................................
C.............  Head of the FDEV-read loop  .............................
C........................................................................

            FS   = 0     ! For facility.pt's
            IREC = 0
111         CONTINUE

C.................  Read a line of facility.pt file; Check input status

                LINE = BLANK5
                IREC = IREC + 1
                READ( FDEV, 93000, END=112, IOSTAT=IOS ) LINE
                IF ( IOS .NE. 0 ) THEN
                    WRITE( MESG, 94010 ) 
     &                  'Error ', IOS,  'reading "' // 
     &                  FCNAME( 1:TRIMLEN( FCNAME ) ) //
     &                  '" at line', IREC
                    CALL M3EXIT( 'EMSPOINT', 0, 0, MESG, 2 )
                END IF

                ZONE = STR2INT( LINE( 43:44 ) )

C.................  Check zone conversion to integer

                IF( ZONE .EQ. IMISS3 ) THEN

                    WRITE( MESG,94010 )
     &                  'Zone blank at line', IREC,
     &                  'in "' // FCNAME( 1:TRIMLEN( FCNAME ) ) // 
     &                  '", so will assume LAT/LON coordinates.'
                    CALL M3MESG( MESG )
                    GO TO 111

                END IF

                FS = FS + 1

                IF ( FS .LE. NPFMAX ) THEN

                    WRITE ( CFIP,93020 ) 
     &                   1000 * STR2INT( LINE( 1:2 ) ) +
     &                          STR2INT( LINE( 3:5 ) )

                    CSS  = LBLANK( LINE( 6:20 ) )
                    FCID = LINE( MIN(CSS+6,20):20 )
                    IF( FCID .EQ. ' ' ) FCID = CMISS3

                    INDXFCA( FS ) = FS
                    FCKEYA ( FS ) = CFIP // FCID
                    FCDESCA( FS ) = LINE( 45:84 )
                    FCCRDXA( FS ) = STR2REAL( LINE( 25:33 ) )
                    FCCRDYA( FS ) = STR2REAL( LINE( 34:42 ) )
                    FCZONEA( FS ) = ZONE

                END IF          !  if fs in bounds

            GO TO  111      !  to head of FDEV-read loop
112         CONTINUE        !  end of the FDEV-read loop

            CLOSE( FDEV )

            WRITE( MESG,94010 ) 
     &       'FACILITY.PT FILE processed:' // CRLF() // BLANK5 //
     &       '   Actual FACILITY record count', FS, CRLF() // BLANK5 //
     &       '   Max dimensioned record-count (NPFMAX)=', NPFMAX

            CALL M3MSG2( MESG )

            IF ( FS .GT. NPFMAX ) THEN

                EFLAG = .TRUE.
                CALL M3MESG(
     &              'Max record-count exceeded in FACILITY file.' )

            ELSE	!  else sort input facility table:

                CALL SORTIC( FS, INDXFCA, FCKEYA )
                DO  113  I = 1, FS
                    J = INDXFCA( I )
                    FCKEY ( I ) = FCKEYA ( J )
                    FCZONE( I ) = FCZONEA( J )
                    FCCRDX( I ) = FCCRDXA( J )
                    FCCRDY( I ) = FCCRDYA( J )
                    FCDESC( I ) = FCDESCA( J )
113             CONTINUE

            END IF		!  if facility table overflow or not


C........................................................................
C.............  Head of the SDEV-read loop  .............................
C........................................................................

            SS   = 0     ! For stack.pt's
            IREC = 0
122         CONTINUE

C.................  Read a line of stack.pt file and check input status

                IREC = IREC + 1
                READ( SDEV, 93000, END=123, IOSTAT=IOS ) LINE
                IF ( IOS .NE. 0 ) THEN
                    WRITE( MESG, 94010 ) 
     &                  'Error ', IOS,  'reading "' // 
     &                  SKNAME( 1:TRIMLEN( SKNAME ) ) //
     &                  '" at line', IREC
                    CALL M3EXIT( 'EMSPOINT', 0, 0, MESG, 2 )
                END IF

C.................  Get lookup into facilities table

                WRITE( CFIP,93020 )  1000*STR2INT( LINE( 1:2 ) ) +
     &                               STR2INT( LINE( 3:5 ) )

C.................  Find source match in facilities table

                CSS  = LBLANK( LINE( 6:20 ) )
                FCID = LINE( MIN(CSS+6,20):20 )
                IF( FCID .EQ. ' ' ) FCID = CMISS3
                FKEY = CFIP // FCID
                K1   = FINDC( FKEY, FS, FCKEY )
                    
                CSS  = LBLANK( LINE( 21:32 ) )
                SKID = LINE( MIN(CSS+21,32):32 )
                IF( SKID .EQ. ' ' ) SKID = CMISS3

                IF ( K1 .LE. 0 ) THEN
                    WRITE( MESG,94010 )
     &                'Stack not in plant recs: FIP=' // CFIP // 
     &                '; Facility=' // FCID( 1:TRIMLEN( FCID ) ) //
     &                '; Stack='    // SKID //
     &                '; line= ', IREC
                    CALL M3MESG( MESG )
                    GO TO 122
                END IF

C.................  Convert units for stack parameters

                STKD = STR2REAL( LINE( 33:40 ) )
                STKH = STR2REAL( LINE( 41:47 ) )
                STKT = STR2REAL( LINE( 48:54 ) ) 
                STKV = STR2REAL( LINE( 55:61 ) )
                STKF = STR2REAL( LINE( 62:71 ) )

                IF ( STKD .GT. 0.0 ) STKD = STKD * FT2M          ! Diam ft to m
                IF ( STKH .GT. 0.0 ) STKH = STKH * FT2M          ! Ht,  ft to m
                IF ( STKT .GT. 0.0 ) STKT = STKT * FTOC  +  CTOK ! temp, degrees K
                IF ( STKV .GT. 0.0 ) STKV = STKV * FT2M          ! Veloc ft/s to m/s
                IF ( STKF .GT. 0.0 ) STKF = STKF * FLWE2M        ! Flow ft^3/min m^3/s
                FSAV = 0.

C.................  Calculate velocity from flow and diameter                
                IF ( STKV .LE. 0.0  .AND.
     &               STKD .GT. 0.0  .AND.
     &               STKF .GT. 0.0 ) THEN
                    STKV = STKF / ( 0.25 * PI * STKD * STKD )

C.................  Compare flow to velocity and diameter
                ELSEIF( ( STKF - STKV*0.25*PI*STKD*STKD ) / STKF
     &                  .GT. 0.001                               ) THEN

                    FSAV = STKF

                END IF

C.............   Resolve coordinates for this stack:

                XLOC = STR2REAL( LINE( 72:80 ) )! supposed to be UTM
                YLOC = STR2REAL( LINE( 81:89 ) )! supposed to be UTM
                
                IF ( XLOC .LE. 0.0  .OR. YLOC .LE. 0.0 ) THEN
                    
                    XLOC = FCCRDX( K1 )
                    YLOC = FCCRDY( K1 )

                    IF( XLOC .LE. 0. .OR. YLOC .LE. 0. ) THEN  ! Warning

                        WRITE( MESG,94010 )
     &                    'Source dropped: ' //
     &                    'bad coords for stack rec', IREC,
     &                     ': FIP='// CFIP // '; Facility=' // 
     &                     FCID( 1:TRIMLEN( FCID ) ) //
     &                    '; and Stack=' // SKID
                        CALL M3MESG( MESG )

                        GO TO 122            ! to end of do-loop

                    END IF              !  if xloc, yloc missing

                END IF

C.............   Create XLOC and YLOC for use at end of loop
C.............   For UTM Grid...

                ZONE = FCZONE( K1 )

                IF ( GDTYP3D .EQ. UTMGRD3 ) THEN   !  convert to UTM output

                    IF ( ZONE .NE. UZONE ) THEN !  UTM input w/ diff zone

                        CALL UTM2LL( XLOC, YLOC, ZONE, XX,   YY )
                        CALL LL2UTM( XX,   YY,  UZONE, XLOC, YLOC )

                    END IF  !  if zone missing (lat-lon), or different

                ELSE IF ( GDTYP3D .EQ. LAMGRD3 ) THEN   !  convert to lambert

                    XX = XLOC
                    YY = YLOC
           
                    IF ( .NOT. UTM2LAM( XX, YY, ZONE,
     &                                  XLOC, YLOC ) ) THEN
                        MESG = 
     &                    'Bad coordinates: UTM2LAM() failed for '//
     &                    'FIP='        // CFIP // 
     &                    '; facility=' // FCID //
     &                    '; Stack='    // SKID
                        CALL M3MESG( MESG )
           
                        WRITE( MESG,94060 )
     &                    '    Setting coordinates to ( ', 
     &                         XORIG3D, ',', YORIG3D, ')'
                        CALL M3MESG( MESG )
           
                        XLOC = XORIG3D
                        YLOC = YORIG3D
           
                    END IF      !  if utm2lam() failed

                ELSE IF ( GDTYP3D .EQ. LATGRD3 ) THEN   !  convert to lat-lon

                    XX = XLOC
                    YY = YLOC
                    CALL UTM2LL( XX, YY, ZONE, XLOC, YLOC )

                END IF       !  if output coord type UTM or Lambert or Lat-Lon


C.................  Store stack file information

                SS = SS + 1

                IF ( SS .LE. NPSMAX ) THEN

                    SKKEYA( SS ) = CFIP // 
     &                             FCID( 1:TRIMLEN( FCID ) ) //
     &                             SKID

                    INDXSKA( SS ) = SS
                    IFCKEYA( SS ) = K1
                    SKDIAMA( SS ) = STKD
                    SKHEITA( SS ) = STKH
                    SKTEMPA( SS ) = STKT
                    SKVELOA( SS ) = STKV
                    SKFLOWA( SS ) = FSAV
                    SKCRDXA( SS ) = XLOC
                    SKCRDYA( SS ) = YLOC

                END IF          !  if ss in bounds

            GO TO  122      !  to head of SDEV-read loop
123         CONTINUE        !  end of the SDEV-read loop

            CLOSE( SDEV )

            WRITE( MESG,94010 ) 
     &       'STACK.PT FILE processed:' // CRLF() // BLANK5 //
     &       '   Actual  STACK   record-count', SS, CRLF() // BLANK5 //
     &       '   Max dimensioned record-count (NPSMAX)=', NPSMAX

            CALL M3MSG2( MESG )

            IF ( SS .GT. NPSMAX ) THEN	!  overflow:

                EFLAG = .TRUE.
                CALL M3MESG(
     &              'Max record-count exceeded in STACK file.' )

            ELSE	!  else sort input stack parameter table:

                CALL SORTIC( SS, INDXSKA, SKKEYA )
                DO  124  I = 1, SS
                    J = INDXSKA( I )
                    SKDIAM( I ) = SKDIAMA( J )
                    IFCKEY( I ) = IFCKEYA( J )
                    SKHEIT( I ) = SKHEITA( J )
                    SKTEMP( I ) = SKTEMPA( J )
                    SKCRDX( I ) = SKCRDXA( J )
                    SKCRDY( I ) = SKCRDYA( J )
                    SKVELO( I ) = SKVELOA( J )
                    SKFLOW( I ) = SKFLOWA( J )
                    SKKEY ( I ) = SKKEYA ( J )
124             CONTINUE

            END IF		!  if stack parameter table overflow or not


C........................................................................
C.............  Head of the PDEV-read loop  .............................
C........................................................................

            PS   = 0     ! For process.pt's
            IREC = 0
133         CONTINUE

C.................  Read a line of process.pt file, and check input status

                IREC = IREC + 1
                READ( PDEV, 93000, END=134, IOSTAT=IOS ) LINE
                IF ( IOS .NE. 0 ) THEN
                    WRITE( MESG, 94010 ) 
     &                  'Error ', IOS,  'reading "' // 
     &                  PSNAME( 1:TRIMLEN( PSNAME ) ) //
     &                  '" at line', IREC
                    CALL M3EXIT( 'EMSPOINT', 0, 0, MESG, 2 )
                END IF

C.................  Convert and check SCC value

                SCC = STR2INT( LINE( 57:64 ) )
                IF ( SCC .LE. 0 ) THEN
                    WRITE( MESG,94010 )
     &                  'Bad SCC:', SCC, 'at line', IREC,
     &                  'in "' // PSNAME( 1:TRIMLEN( PSNAME ) ) // '"'
                    CALL M3MESG( MESG )
                    GO TO 133
                END IF

                PS = PS + 1

                IF ( PS .LE. NPPMAX ) THEN

                    WRITE ( CFIP,93020 ) 
     &                   1000 * STR2INT( LINE( 1:2 ) ) +
     &                          STR2INT( LINE( 3:5 ) )

                    CSS  = LBLANK( LINE( 6:20 ) )
                    FCID = LINE( MIN(CSS+6,20):20 )
                    IF( FCID .EQ. ' ' ) FCID = CMISS3

                    CSS  = LBLANK( LINE( 21:32 ) )
                    SKID = LINE( MIN(CSS+21,32):32 )
                    IF( SKID .EQ. ' ' ) SKID = CMISS3

                    CSS  = LBLANK( LINE( 33:44 ) )
                    DVID = LINE( MIN(CSS+33,44):44 )
                    IF( DVID .EQ. ' ' ) DVID = CMISS3

                    CSS  = LBLANK( LINE( 45:56 ) )
                    PRID = LINE( MIN(CSS+45,56):56 )
                    IF( PRID .EQ. ' ' ) PRID = CMISS3

                    PSKEYA( PS ) = CFIP // 
     &                             FCID( 1:TRIMLEN( FCID ) ) //
     &                             SKID( 1:TRIMLEN( SKID ) ) //
     &                             DVID( 1:TRIMLEN( DVID ) ) //
     &                             PRID( 1:TRIMLEN( PRID ) )

                    INDXPA( PS ) = PS
                    PSSCCA( PS ) = SCC

                END IF          !  if ps in bounds

            GO TO  133      !  to head of PDEV-read loop
134         CONTINUE        !  end of the PDEV-read loop

            CLOSE( PDEV )

            WRITE( MESG,94010 ) 
     &       'PROCESS.PT FILE processed' // CRLF() // BLANK5 //
     &       '   Actual PROCESS  record-count', PS, CRLF() // BLANK5 //
     &       '   Max dimensioned record-count (NPPMAX)=', NPPMAX

            CALL M3MSG2( MESG )

            IF ( PS .GT. NPPMAX ) THEN		!  overflow

                EFLAG = .TRUE.
                CALL M3MESG(
     &              'Max record-count exceeded in PROCESS file.' )

            ELSE	!  sort input process parameter table:

                CALL SORTIC( PS, INDXPA,  PSKEYA )
                DO  135  I = 1, PS
                    J = INDXPA( I )
                    PSSCC( I ) = PSSCCA( J )
                    PSKEY( I ) = PSKEYA( J )
135             CONTINUE

            END IF	!  if input process parameter table overflow or not


C........................................................................
C.............  Head of the DDEV-read loop  .............................
C........................................................................

            DS   = 0     ! For device.pt's
            IREC = 0
144         CONTINUE

                IREC = IREC + 1
                READ( DDEV, 93000, END=145, IOSTAT=IOS ) LINE
                IF ( IOS .NE. 0 ) THEN
                    WRITE( MESG, 94010 ) 
     &                  'Error ', IOS,  'reading "' // 
     &                  DVNAME( 1:TRIMLEN( DVNAME ) ) //
     &                  '" at line', IREC
                    CALL M3EXIT( 'EMSPOINT', 0, 0, MESG, 2 )
                END IF

C.................  Convert and check SIC value

                SIC = STR2INT( LINE( 45:48 ) )
                IF ( SIC .EQ. 0 ) THEN

                    WRITE( MESG,94010 )
     &                  'Default SIC:', SIC, 'at line', IREC,
     &                  'in "' // DVNAME( 1:TRIMLEN( DVNAME ) ) // '"'
                    CALL M3MESG( MESG )

                ELSE IF ( SIC .LT. 0 ) THEN

                    WRITE( MESG,94010 )
     &                  'Missing SIC:', SIC, 'at line', IREC,
     &                  'in "' // DVNAME( 1:TRIMLEN( DVNAME ) ) // '"'
                    CALL M3MESG( MESG )

                    WRITE( MESG,94010 )
     &                  '     Setting to default of 0'
                    CALL M3MESG( MESG )
                    SIC = 0

                END IF          !  if sic zero or not

C.................  Convert and check temporal profile numbers

C temp          IMON = STR2INT( LINE( ) )
                IWEK = STR2INT( LINE( 123:124 ) )
                IDIU = STR2INT( LINE( 121:122 ) )

                IF( IWEK .LE. 0 ) THEN
                    WRITE( MESG,94010 ) 
     &                  'Default weekly profile', IWEK, 'at line', IREC,
     &                  'in "' // DVNAME( 1:TRIMLEN( DVNAME ) ) // '"'
                    CALL M3MESG( MESG )
                END IF

                IF( IDIU .LE. 0 ) THEN
                    WRITE( MESG,94010 ) 
     &              'Default hourly profile', IDIU, 'at line', IREC,
     &              'in "' // DVNAME( 1:TRIMLEN( DVNAME ) ) // '"'
                    CALL M3MESG( MESG )
                END IF

                DS = DS + 1

                IF ( DS .LE. NPDMAX ) THEN

                    WRITE ( CFIP,93020 ) 
     &                   1000 * STR2INT( LINE( 1:2 ) ) + 
     &                          STR2INT( LINE( 3:5 ) )

                    CSS  = LBLANK( LINE( 6:20 ) )
                    FCID = LINE( MIN(CSS+6,20):20 )
                    IF( FCID .EQ. ' ' ) FCID = CMISS3

                    CSS  = LBLANK( LINE( 21:32 ) )
                    SKID = LINE( MIN(CSS+21,32):32 )
                    IF( SKID .EQ. ' ' ) SKID = CMISS3

                    CSS  = LBLANK( LINE( 33:44 ) )
                    DVID = LINE( MIN(CSS+33,44):44 )
                    IF( DVID .EQ. ' ' ) DVID = CMISS3

                    INDXDVA( DS ) = DS
                    DVKEYA ( DS ) = CFIP // 
     &                              FCID( 1:TRIMLEN( FCID ) ) //
     &                              SKID( 1:TRIMLEN( SKID ) ) //
     &                              DVID( 1:TRIMLEN( DVID ) )

                    DVSICA ( DS ) = SIC
C temp              DVIMONA( DS ) = IMON
                    DVIWEKA( DS ) = IWEK
                    DVIDIUA( DS ) = IDIU

                END IF          !  if ds in bounds

            GO TO  144      !  to head of DDEV-read loop
145         CONTINUE        !  end of the DDEV-read loop

            CLOSE( DDEV )

            WRITE( MESG,94010 ) 
     &       'DEVICE.PT FILE processed' // CRLF() // BLANK5 //
     &       '   Actual  DEVICE  record-count', DS, CRLF() // BLANK5 //
     &       '   Max dimensioned record-count (NPDMAX)=', NPDMAX

            CALL M3MSG2( MESG )

            IF ( DS .GT. NPDMAX ) THEN	!  overflow

                EFLAG = .TRUE.
                CALL M3MESG(
     &              'Max record-count exceeded in DEVICE file.' )

            ELSE	!  sort input device parameter table:

                CALL SORTIC( DS, INDXDVA, DVKEYA )
                DO  146  I = 1, DS
                    J = INDXDVA( I )
                    DVSIC ( I ) = DVSICA ( J )
                    DVIWEK( I ) = DVIWEKA( J )
                    DVIDIU( I ) = DVIDIUA( J )
                    DVKEY ( I ) = DVKEYA ( J )
146             CONTINUE

            END IF	!  if evice parameter table overflow or not


C............  Read rule effectiveness file
C................. NOTE: Will be added later.  Not needed for SESARM

            RULFLAG = .FALSE.   ! No rule effectiveness files


C........................................................................
C.............  Head of the EDEV-read loop  .............................
C........................................................................

            IREC =   0
            RPEN = 100.0
            REFF = 100.0
188         CONTINUE

C.................  Read a line of emission.pt file and check input status

                IREC = IREC + 1
                READ( EDEV, 93000, END=199, IOSTAT=IOS ) LINE
                IF ( IOS .NE. 0 ) THEN
                    WRITE( MESG, 94010 ) 
     &                  'Error ', IOS,  'reading "' // 
     &                  EMNAME( 1:TRIMLEN( EMNAME ) ) //
     &                  '" at line', IREC
                    CALL M3EXIT( 'EMSPOINT', 0, 0, MESG, 2 )
                END IF

C.................  Check pollutant code and set index I

                CSS  = LBLANK( LINE( 57:61 ) )
                CPOL = LINE( MIN(CSS+57,61):61 )
                COD  = INDEX1( CPOL( 1:TRIMLEN( CPOL ) ), NIPOL, EINAM )
                IF ( I  .LE. 0 ) THEN
                    WRITE( MESG,94010 ) 
     &                  'Bad line', IREC, 
     &                  'Pollutant code "' // LINE(57:61) // 
     &                  '" in "' // EMNAME( 1:TRIMLEN( EMNAME ) ) // '"'
                    CALL M3MESG( MESG )
                    GO TO  188      !  to head of loop
                END IF

C.................  Check and set emissions value

                EMIS  = STR2REAL( LINE( 88:100 ) )
                IF ( EMIS .LT. 0.0 )  THEN
                    WRITE( MESG,94010 ) 
     &                  'Bad line', IREC, '" in "' //
     &                  EMNAME( 1:TRIMLEN( EMNAME ) ) // 
     &                  '"--emis value "' // LINE( 88:100 ) // '"'
                    CALL M3MESG( MESG )
                    GO TO  188
                END IF

C.................  Set emission.pt file source arrays and temporary arrays
C.................  Assume that all sources will appear in emission.pt files
C.................  and give warning later if this is not the case.

                FIP  = 1000 * STR2INT( LINE( 1:2 ) ) +
     &                        STR2INT( LINE( 3:5 ) )
                WRITE( CFIP,93020 ) FIP
                CSS  = LBLANK( LINE( 6:20 ) )

                FCID = LINE( MIN(CSS+6,20):20 )
                FKEY = CFIP // FCID

                CSS  = LBLANK( LINE( 21:32 ) )
                SKID = LINE( MIN(CSS+21,32):32 )
                IF( SKID .EQ. ' ' ) SKID = CMISS3 
                SKEY = FKEY( 1:TRIMLEN( FKEY ) ) // SKID
                K2   = FINDC( SKEY, SS, SKKEY )

                IF( K2 .LE. 0 ) THEN   ! Key not found
                    WRITE( MESG,94010 )
     &                'Source dropped: not in stack recs: FIP=' 
     &                              // CFIP // 
     &                '; Facility=' // FCID( 1:TRIMLEN( FCID ) ) //
     &                '; Stack='    // SKID
                    CALL M3MESG( MESG )

                    NDROP = NDROP + 1
                    EDROP( COD ) = EDROP( COD ) + EMIS

                    GO TO 188  ! To next emissions source line

                END IF          !  if stack key not found

                CSS  = LBLANK( LINE( 33:44 ) )
                DVID = LINE( MIN(CSS+33,44):44 )
                IF( DVID .EQ. ' ' ) DVID = CMISS3
                DKEY = SKEY( 1:TRIMLEN( SKEY ) ) // DVID
                K3   = FINDC( DKEY, DS, DVKEY )

                IF( K3 .LE. 0 ) THEN    ! Key not found
                    WRITE( MESG,94010 ) 
     &                    'Source dropped: not in device recs: FIP=' 
     &                                  // CFIP // 
     &                    '; Facility=' // FCID( 1:TRIMLEN( FCID ) ) //
     &                    '; Stack='    // SKID( 1:TRIMLEN( SKID ) ) //
     &                    '; Device='   // DVID
                    CALL M3MESG( MESG )

                    NDROP = NDROP + 1
                    EDROP( COD ) = EDROP( COD ) + EMIS

                    GO TO 188  ! To next source line

                END IF          !  if device key not found

                CSS  = LBLANK( LINE( 45:56 ) )
                PRID = LINE( MIN(CSS+45,56):56 )
                IF( PRID .EQ. ' ' ) PRID = CMISS3
                PKEY = DKEY( 1:TRIMLEN( DKEY ) ) // PRID
                K4   = FINDC( PKEY, PS, PSKEY )

                IF( K4 .LE. 0 ) THEN    ! Key not found
                    WRITE( MESG,94010 )
     &                    'Source dropped: not in process recs: FIP=' //
     &                     CFIP // 
     &                    '; Facility=' // FCID( 1:TRIMLEN( FCID ) ) //
     &                    '; Stack='    // SKID( 1:TRIMLEN( SKID ) ) //
     &                    '; Device='   // DVID( 1:TRIMLEN( DVID ) ) //
     &                    '; Process='  // PRID
                    CALL M3MESG( MESG )

                    NDROP = NDROP + 1
                    EDROP( COD ) = EDROP( COD ) + EMIS

                    GO TO 188  ! To next emissions source line

                END IF          !   if process key not found

                SCC = PSSCC ( K4 )

C.................  Check and set time period type (Year/day/hourly)

                IF ( LINE( 114:115 ) .EQ. 'AA' ) THEN 

                    TPF = MTPRFAC * WTPRFAC     !  use month, week profiles

                ELSE IF ( LINE( 114:115 ) .EQ. 'AD' ) THEN 

                    TPF  = WTPRFAC              !  use week profiles
                    EMIS = DAY2YR * EMIS

                ELSE IF ( LINE( 114:115 ) .EQ. 'DS' ) THEN

                    TPF = 1                     !  use only hourly profiles

                ELSE                            !  unrecognized type

                    NDROP = NDROP + 1
                    EDROP( COD ) = EDROP( COD ) + EMIS
                    WRITE( MESG,94010 ) 
     &                  'Bad line', IREC, '" in "' //
     &                  EMNAME( 1:TRIMLEN( EMNAME ) ) // 
     &                  '"--Unsupported time period type "' // 
     &                  LINE( 57:58 ) // '"'
                    CALL M3MESG( MESG )
                    GO TO  188          !  to head of EDEV-read loop

                END IF          !  tests on record type line( 57:58 )

C.................  Set source control parameters (will get rule effectiveness
C.................  from another file.

                CEFF  = STR2REAL( LINE( 126:132 ) )
                IF ( CEFF .LT. 0.0 )  THEN 
                    CEFF = 0.0
                ELSE
                    CEFF = CEFF * 100.0
                END IF


C.............  (later:  Get rule pen, eff from prulefac records)

                IF( RULFLAG ) THEN   ! If rule effectiveness files exists
                    REFF = 100.0     ! (will do match on ASCII keys here)
                END IF

                ES = ES + 1

                IF ( ES .LE. NPEMAX ) THEN

                    K1 = IFCKEY( K2 )
                    INDEXA( ES ) = ES         !  index-table for use in SORTI4()
                    IFIPA ( ES ) = FIP
                    IPLANA( ES ) = K1
                    ISTACA( ES ) = K2
                    ISCCA ( ES ) = SCC
                    ICODEA( ES ) = COD       !  subscript into EINAM(*)
                    TPFLGA( ES ) = TPF
                    INVYRA( ES ) = INY
                    STKHTA( ES ) = SKHEIT( K2 )
                    STKDMA( ES ) = SKDIAM( K2 )
                    STKTKA( ES ) = SKTEMP( K2 )
                    STKVEA( ES ) = SKVELO( K2 )
                    STKFLA( ES ) = SKFLOW( K2 )
                    XLOCAA( ES ) = SKCRDX( K2 )
                    YLOCAA( ES ) = SKCRDY( K2 )
                    ISICA ( ES ) = DVSIC ( K3 )
                    IWEKA ( ES ) = DVIWEK( K3 )
                    IDIUA ( ES ) = DVIDIU( K3 )
                    CTLEFA( ES ) = CEFF
                    RULPEA( ES ) = RPEN
                    EMISVA( ES ) = EMIS
                    EMFCID( ES ) = FCID
                    EMDVID( ES ) = DVID
                    EMSKID( ES ) = SKID
                    EMDESC( ES ) = FCDESC( K1 )

                END IF          !  if S in range

            GO TO  188      !  to head of EDEV-read loop

199         CONTINUE        !  end of the EDEV-read loop

            CLOSE( EDEV )

            WRITE( MESG,94010 ) 
     &       'EMISSION.PT FILE processed:'  // CRLF() // BLANK5 //
     &       '   Actual     EMS-95 SOURCE  record-count', ES - NSRC,
     &        CRLF() // BLANK5 //
     &       '   Cumulative EMS-95 SOURCE  record-count', ES

            CALL M3MSG2( MESG )

            NSRC = ES        !  point source emissions lines

        GO TO 101   !  to head of loop on input files


C........................................................................

202     CONTINUE        !  exit from loop on input files

        WRITE( MESG,94010 ) 
     &      'Actual  EMS-95 SOURCE  record-count', NSRC, 
     &      CRLF() // BLANK5 //
     &      'Max dimensioned record-count (NPEMAX)=', NPEMAX

        CALL M3MSG2( MESG )

        IF ( NSRC .GT. NPEMAX ) THEN

            EFLAG = .TRUE.
            CALL M3EXIT( 'EMSPOINT', 0, 0,
     &          'Max record-count exceeded in EMISSION files.', 2 )

        ELSE IF ( EFLAG ) THEN

            CALL M3EXIT( 'EMSPOINT', 0, 0, 
     &           'Error reading POINT SOURCE FOUNDATION files.', 2 )

        END IF		!  if overflow or if errors

C.........  Report how many sources were dropped

        IF( NDROP .GT. 0 ) THEN

            WRITE( MESG,94010 ) 'NOTE:', NDROP, 
     &           'sources dropped. New source count is', NSRC
            CALL M3MSG2( MESG )

            DO 211 I = 1, NIPOL

                WRITE( MESG,94060 ) 
     &             '     Dropped ' 
     &             // EINAM( I )( 1:TRIMLEN(EINAM( I )) ) //
     &             ' emissions: ', EDROP( I ), 'tons/year'
                CALL M3MSG2( MESG )

211         CONTINUE

        END IF          !  if ndrop > 0


C.......   Open ASCII output files:

        ADEV = PROMPTFFILE( 
     &          'Enter the name of the ACTUAL SCC output file',
     &          .FALSE., .TRUE., 'PSCC', 'EMSPOINT' )

        TDEV = PROMPTFFILE( 
     &          'Enter the name of the TEMPORAL XREF output file',
     &          .FALSE., .TRUE., 'PTREF', 'EMSPOINT' )

        CDEV = PROMPTFFILE( 
     &          'Enter the name of the FACILITY/STACK output file',
     &          .FALSE., .TRUE., 'PSRC', 'EMSPOINT' )


C.......   Use SORTI4() to perform an indirect sort by FIPS,SCC,PLANT,STACK
C.......   then permute the records into output order according to the result:

        CALL SORTI4( NSRC, INDEXA, IFIPA, ISCCA, IPLANA, ISTACA )

        IF ( SFLAG ) THEN       !  emissions initialization value
            VMISS = BADVAL3     !  (if sflag, need to flag uninitialized values)
        ELSE
            VMISS = 0.0
        END IF

        LFIP  = -1
        LSCC  = -1
        LLSCC = -1
        LPLT  = -1
        LSTK  = -1
        J     =  0
        L     =  0

        DO  322  S = 1, NSRC

            I = INDEXA( S )

            FIP   = IFIPA ( I )
            SCC   = ISCCA ( I )
            PLANT = IPLANA( I )
            STACK = ISTACA( I )

C.............  If first time got to this source...

            IF ( FIP   .NE. LFIP  .OR.
     &           SCC   .NE. LSCC  .OR.
     &           PLANT .NE. LPLT  .OR.
     &           STACK .NE. LSTK       ) THEN

                LFIP = FIP
                LSCC = SCC
                LPLT = PLANT
                LSTK = STACK
                J    = J + 1		!  actual-sources counter

C.................  Write out ASCII EMS-95 / SMOKE facility-stack table,

                WRITE( CDEV, 93050 ) 
     &              FIP, SCC, PLANT, STACK,
     &              EMFCID( I ), EMSKID( I ), EMDVID( I ), EMDESC( I )

                IF( J .LE. NPSRC ) THEN

                    IF ( SCC .NE. LLSCC ) THEN
                        LLSCC       = SCC
                        L           = L + 1
                        IDXSCC( L ) = J
                    END IF

                    WRITE( TDEV, 93040 ) 1, IWEKA( I ), IDIUA( I )

                    DO  311  V = 1, NIPOL       !  initializations
                        EMISV ( J,V ) = VMISS
                        RULPEN( J,V ) = VMISS
                        RULEFF( J,V ) = VMISS
                        CTLEFF( J,V ) = VMISS
311                 CONTINUE


C.....................  Set final arrays

                    IFIP  ( J ) = FIP
                    ISCC  ( J ) = SCC
                    ISIC  ( J ) = ISICA ( I )
                    IPLANT( J ) = PLANT
                    ISTACK( J ) = STACK
                    TPFLAG( J ) = TPFLGA( I )
                    INVYR ( J ) = INVYRA( I )
                    XLOCA ( J ) = XLOCAA( I )
                    YLOCA ( J ) = YLOCAA( I )
                    STKHT ( J ) = STKHTA( I )
                    STKDM ( J ) = STKDMA( I )
                    STKTK ( J ) = STKTKA( I )

C.....................  Where flow is defined, and option is set to on, need
C.....................    to reset velocity based on flow, because they are 
C.....................    inconsistent.
C.....................  We've waited to do this so that the change can
C.....................    be reported with FIP, SCC, PLT, STK
                    IF( CFLAG .AND. STKFLA( I ) .GT. 0. ) THEN

                        VNEW = STKFLA( I ) / 
     &                         ( 0.25 * PI * STKDMA( I ) * STKDMA( I ) )

                        WRITE( MESG,94020 ) 
     &                    'Stack V reset using flow. FIP:', FIP, 
     &                    'SCC:', SCC, 'PLT:', PLANT, 'STK:', STACK
                        CALL M3MESG( MESG ) 

                        WRITE( LDEV,93060 ) 
     &                     BLANK5 // '  F:', STKFLA( I ), 'D:', 
     &                     STKDMA( I ), 'OLD V:', STKVEA( I ),
     &                     'NEW V:', VNEW

                        STKVE ( J ) = VNEW

                    ELSE
                        STKVE ( J ) = STKVEA( I )

                    ENDIF
                    
                    K = FIND1( FIP, NZF, TFIPEF )
                    IF ( K .GT. 0 ) THEN
                        TZONES( J ) = TZONEF( K )
                    ELSE
                        K = FIND1( FIP/1000, NZS, TFIPST )
                        IF ( K .GT. 0 ) THEN
                            TZONES( J ) = TZONST( K )
                        ELSE
                            TZONES( J ) = TZONE0
                        END IF
                    END IF

                ELSE            !  would have overflowed arrays:

                    WRITE( MESG,94010 ) 
     &                'Skipping FIP:', FIP, 
     &                'SCC:',          SCC, 
     &                'PLT:',          PLANT,
     &                'STK:',          STACK,
     &                'to avoid overflow.'

                    CALL M3MESG( MESG )

                    GO TO 322   ! To end of loop

                END IF          !  if j <= npsrc or not

C.............  Keep updating previous info so sources count
C.............  will be correct
            ELSE
                LFIP = FIP
                LSCC = SCC
                LPLT = PLANT
                LSTK = STACK

            END IF              !  if first encounter with this source

C.............  Now set emissions values, but only if we haven't exceeded
C.............  the correct number of sources

            IF ( J .LE. NPSRC ) THEN

                V    = ICODEA( I )
                EMIS = EMISVA( I )

                IF ( EMISV( J,V ) .LE. 0.0 ) THEN  !  initialize if necessary

                    EMISV ( J,V ) = EMIS
                    CTLEFF( J,V ) = CTLEFA( I )
                    RULEFF( J,V ) = RULEFA( I )
                    RULPEN( J,V ) = RULPEA( I )

                ELSE		!  else accumulate emissions

                    EMISV ( J,V ) = EMISV ( J,V ) + EMIS

                END IF

            END IF      !  if j <= npsrc

322     CONTINUE

        NSCC   = L      !  scc codes found
        NPOINT = J

        IF ( NPOINT .NE. NPSRC ) THEN

            WRITE( MESG,94010 ) 
     &       'Actual      number of sources', NPOINT, CRLF()// BLANK5//
     &       'Dimensioned number of sources (NPSRC)=', NPSRC,  
     &       CRLF()// BLANK5// 'Do not match!'

            CALL M3MSG2( MESG )

            CALL M3EXIT( 'EMSPOINT', 0, 0, 
     &                   'Error because of number of sources.', 2 )

        END IF

        IF ( SFLAG ) THEN       !  if missing records are fatal:

            EFLAG = .FALSE.

            DO  334  V = 1, NIPOL
            DO  333  S = 1, NPSRC

               IF ( EMISV( S,V ) .LT. AMISS3 ) THEN

                    EFLAG = .TRUE.
                    WRITE( MESG,94020 )
     &                  'Missing record:  FIP:', IFIP( S ),
     &                  'SCC:',     ISCC( S ), 
     &                  'Plant:',   IPLANT( S ),
     &                  'Stack:',   ISTACK( S ),
     &                  'Species: ' // EINAM( V )
                    CALL M3MESG( MESG )

                END IF

333         CONTINUE
334         CONTINUE

            IF ( EFLAG ) THEN
                CALL M3EXIT( 'EMSPOINT', 0, 0,
     &                       'Missing species recs in input file', 2 )
            END IF

        END IF          !  sflag or not:  check for non-set EMIS values


C...........   Pad rest of record with zeros, if appropriate:

        DO  342  I = 1, NIPOL   !  empty=record initializations
        DO  341  S = NPOINT+1, NPSRC
            EMISV ( S,I ) = 0.0
            RULPEN( S,I ) = 0.0
            RULEFF( S,I ) = 0.0
            CTLEFF( S,I ) = 0.0
341     CONTINUE
342     CONTINUE                !  end initializations

C.......   Fix stack parameters

        CALL FIXSTK( RDEV, NPOINT, IFIP, ISCC, IPLANT, ISTACK,
     &               STKHT, STKDM, STKTK, STKVE, NIPOL, EMISV )

C.......   Get file name; open output point sources file

        FTYPE3D = GRDDED3
        P_ALP3D = DBLE( UZONE )
        P_BET3D = DBLE( AMISS3 )
        P_GAM3D = DBLE( AMISS3 )
        XCENT3D = 0.0D0
        YCENT3D = 0.0D0
        XORIG3D = DBLE( AMISS3 )
        YORIG3D = DBLE( AMISS3 )
        SDATE3D = 0 !  n/a
        STIME3D = 0 !  n/a
        TSTEP3D = 0             !  time independent
        NVARS3D = 4 * NIPOL + 14
        NCOLS3D = 1
        NROWS3D = NPSRC    !  number of rows = # of point sources.
        NLAYS3D = 1
        NTHIK3D = 1
        GDTYP3D = GDTYP3D
        VGTYP3D = IMISS3
        VGTOP3D = AMISS3
        GDNAM3D = COORDN

        FDESC3D( 1 ) = 'Point Source emissions values.'

        DO  411  I = 2, MXDESC3
            FDESC3D( I ) = ' '
411     CONTINUE

        J = 1
        VNAME3D( J ) = 'IFIP'
        VTYPE3D( J ) = M3INT
        UNITS3D( J ) = 'n/a'
        VDESC3D( J ) = 'FIP code for counties'
        J = J + 1

        VNAME3D( J ) = 'ISIC'
        VTYPE3D( J ) = M3INT
        UNITS3D( J ) = 'n/a'
        VDESC3D( J ) = 'Source Industrial Code'
        J = J + 1

        VNAME3D( J ) = 'ISCC'
        VTYPE3D( J ) = M3INT
        UNITS3D( J ) = 'n/a'
        VDESC3D( J ) = 'Point Source Category Code'
        J = J + 1

        VNAME3D( J ) = 'IPLANT'
        VTYPE3D( J ) = M3INT
        UNITS3D( J ) = 'n/a'
        VDESC3D( J ) = 'Point Source plant-ID'
        J = J + 1

        VNAME3D( J ) = 'ISTACK'
        VTYPE3D( J ) = M3INT
        UNITS3D( J ) = 'n/a'
        VDESC3D( J ) = 'Point Source stack-ID'
        J = J + 1

        VNAME3D( J ) = 'TZONES'
        VTYPE3D( J ) = M3INT
        UNITS3D( J ) = 'n/a'
        VDESC3D( J ) = 'Time zone for site'
        J = J + 1

        VNAME3D( J ) = 'TPFLAG'
        VTYPE3D( J ) = M3INT
        UNITS3D( J ) = 'T|2? T|3?'
        VDESC3D( J ) = 'Use week(2), month(3) temporal profiles or not'
        J = J + 1

        VNAME3D( J ) = 'INVYR'
        VTYPE3D( J ) = M3INT
        UNITS3D( J ) = 'year AD'
        VDESC3D( J ) = 'Year of inventory for this record'
        J = J + 1

        VNAME3D( J ) = 'XLOCA'
        VTYPE3D( J ) = M3REAL
        IF( GDTYP3D .EQ. LATGRD3 ) THEN
            UNITS3D( J ) = 'deg LON'
        ELSE
            UNITS3D( J ) = 'meters'
        END IF
        VDESC3D( J ) = 'X coordinate for site'
        J = J + 1

        VNAME3D( J ) = 'YLOCA'
        VTYPE3D( J ) = M3REAL
        IF( GDTYP3D .EQ. LATGRD3 ) THEN
            UNITS3D( J ) = 'deg LAT'
        ELSE
            UNITS3D( J ) = 'meters'
        END IF
        VDESC3D( J ) = 'Y coordinate for site'
        J = J + 1

        VNAME3D( J ) = 'STKHT'
        VTYPE3D( J ) = M3REAL
        UNITS3D( J ) = 'm'
        VDESC3D( J ) = 'Stack height'
        J = J + 1

        VNAME3D( J ) = 'STKDM'
        VTYPE3D( J ) = M3REAL
        UNITS3D( J ) = 'm'
        VDESC3D( J ) = 'Stack diameter'
        J = J + 1

        VNAME3D( J ) = 'STKTK'
        VTYPE3D( J ) = M3REAL
        UNITS3D( J ) = 'deg K'
        VDESC3D( J ) = 'Stack exhaust temperature'
        J = J + 1

        VNAME3D( J ) = 'STKVE'
        VTYPE3D( J ) = M3REAL
        UNITS3D( J ) = 'm/s'
        VDESC3D( J ) = 'Stack exhaust velocity'
        J = J + 1

        DO  422  V = 1 , NIPOL

            VNAME3D( J ) = 'CTLEFF_' // EINAM( V )
            VTYPE3D( J ) = M3REAL
            UNITS3D( J ) = 'n/a'
            VDESC3D( J ) = 
     &      'Control efficiency (in [0,100], or "MISSING": < -9.0E36)'
            J = J + 1

422     CONTINUE        !  end loop on inventory pollutants I

        DO  433  V = 1, NIPOL

            VNAME3D( J ) = 'RULEFF_' // EINAM( V )
            VTYPE3D( J ) = M3REAL
            UNITS3D( J ) = 'n/a'
            VDESC3D( J ) = 
     &      'Rule Effectiveness  (in [0,100], or "MISSING": < -9.0E36)'
            J = J + 1

433     CONTINUE        !  end loop on inventory pollutants I

        DO  444  V = 1 , NIPOL

            VNAME3D( J ) = 'RULPEN_' // EINAM( V )
            VTYPE3D( J ) = M3REAL
            UNITS3D( J ) = 'n/a'
            VDESC3D( J ) = 
     &      'Rule penetration (in [0,100], or "MISSING": < -9.0E36)'
            J = J + 1

444     CONTINUE        !  end loop on inventory pollutants I

        DO  455  V = 1 , NIPOL

            VNAME3D( J ) = EINAM( V )
            VTYPE3D( J ) = M3REAL
            UNITS3D( J ) = 'tons/year'
            VDESC3D( J ) = 
     &        LINE( 1:1 ) // LINE( 57:58 )  // ' emissions totals'
            J = J + 1

455     CONTINUE        !  end loop on inventory pollutants I

        ENAME = PROMPTMFILE( 
     &          'Enter logical name for POINTS output file',
     &          FSUNKN3, 'PNTS', 'EMSPOINT' )

C.......   Write out the point source emissions values:

        CALL M3MSG2( 'Writing out POINT output file...' )

        IF ( .NOT. WRITE3( ENAME, ALLVAR3, 0, 0,  IFIP ) ) THEN
            CALL M3EXIT( 'EMSPOINT', 0, 0, 
     &                   'Error writing output file "' //
     &                   ENAME( 1:TRIMLEN( ENAME ) ) // '"',  2 )
        END IF


C.......   Sort and write out the point source scc values:

        CALL SORTI1( NSCC, IDXSCC, ISCC )

        LSCC = IMISS3
        DO 466 I = 1, NSCC

            SCC = ISCC( IDXSCC( I ) )

            IF( SCC .NE. LSCC ) THEN
                WRITE( ADEV, 93030 ) SCC
            END IF

            LSCC = SCC

466     CONTINUE


999     CONTINUE          !  exit program:  normal completion

        MESG = 'Successful completion of Program EMSPOINT'
        CALL M3EXIT( 'EMSPOINT', 0, 0, MESG, 0 )

C.........  Error because improper grouping of raw input files
1005    MESG = 'EMS95 input file list ended without complete set'
        CALL M3EXIT( 'EMSPOINT', 0, 0, MESG, 2 )

C.........  Error opening raw input file
1006    MESG = 'Error opening file ' // NAMTMP( 1:TRIMLEN( NAMTMP ) )
        CALL M3EXIT( 'EMSPOINT', 0, 0, MESG, 2 )

C.........  Error with INVYEAR packet read
1007    CALL M3EXIT( 'EMSPOINT', 0, 0, 'ERROR: Can only use INVYEAR ' //
     &               'packet in PNLST on increments of five files', 2 )

C******************  FORMAT  STATEMENTS   ******************************

C...........   Informational (LOG) message formats... 92xxx

92000   FORMAT( 5X, A )

92010   FORMAT( 5X, A, :, I10 )


C...........   Formatted file I/O formats............ 93xxx

93000   FORMAT( A )

93020   FORMAT( I5.5 )

93030   FORMAT( I8.8 )

93040   FORMAT( I5, ',' ,I5, ',' ,I5 )

C93041   FORMAT( I5, X, I5, X, I3, X, I8, X, I5.5, 3( X, I3 ) )

93050   FORMAT( I5.5, ',', I8.8, ',', 2 (I6, ','), 
     &          A15 , ',', A12 , ',', A12 , ',', A40 )

93060   FORMAT( 10( A, :, E10.3, :, 1X ) )

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I7, :, 1X ) )

94020   FORMAT( A, 1X, I5.5, 1X, A, 1X, I8.8, 1X,
     &          A, I6, X, A, I6, X, A, :, I6 )

94040   FORMAT( A, I2.2 )

94060   FORMAT( 10( A, :, E10.3, :, 1X ) )

94080   FORMAT( '************  ', A, I7, ' ,  ' , A, I12 )
 
        END

